(function () {
  'use strict';

  const _cs = getComputedStyle(document.documentElement);
  const _c   = n => _cs.getPropertyValue(n).trim();
  const _rgb  = n => { const h = _c(n).replace('#',''); const v = parseInt(h,16); return [(v>>16)&0xFF,(v>>8)&0xFF,v&0xFF]; };
  const _rgba = (n, a) => { const [r,g,b] = _rgb(n); return `rgba(${r},${g},${b},${a})`; };
  const _col  = n => { const [r,g,b] = _rgb(n); return new THREE.Color(r/255, g/255, b/255); };

  /* ── Inject CSS ── */
  (function () {
    if (document.getElementById('sho-styles')) return;
    const s = document.createElement('style');
    s.id = 'sho-styles';
    s.textContent = `
      #sho-ctrl-panel { display:flex; flex-direction:column; overflow:hidden; }
      #sho-qplot-section { flex:1; min-height:0; display:flex; flex-direction:column; padding:6px 12px 10px; }
      #sho-qplot { flex:1; min-height:0; width:100%; display:block; }
    `;
    document.head.appendChild(s);
  })();

  /* ── Parameters ── */
  let omega0 = 1.0;    // natural frequency √(k/m)
  let gammaD = 0.0;    // damping coefficient
  let driveW = 0.0;    // driving frequency (0 = wall at rest)

  /* ── State ── */
  let pos = 0, vel = 0, simTime = 0;
  let running = false, frameId = null;

  /* ── World geometry: wall at x = -HALF, mass equilibrium at x = 0 ── */
  const HALF     = 4;      // spring rest length (wall face → mass center)
  const X_MAX    = 3;      // max |displacement| in world units
  const DRIVE_A  = 0.1;    // wall vibration amplitude (small: 1/40 of rest length)
  const BALL_R   = 0.45;

  function wallXNow() {
    return -HALF + (driveW > 0 ? DRIVE_A * Math.sin(driveW * simTime) : 0);
  }
  function springRad() { return Math.max(0.02, Math.min(0.10, 0.05 * Math.pow(omega0, 0.9))); }

  /* ── Physics: base-excited damped oscillator
     ẍ = -ω₀²(x - x_w) - γ ẋ,   x_w = A sin(ω_d t) ── */
  const DT = 0.05, SUBSTEPS = 8;
  const dt = DT / SUBSTEPS;

  function step() {
    for (let s = 0; s < SUBSTEPS; s++) {
      const xw = driveW > 0 ? DRIVE_A * Math.sin(driveW * simTime) : 0;
      const a  = -omega0 * omega0 * (pos - xw) - gammaD * vel;
      vel += a * dt;
      pos += vel * dt;
      simTime += dt;
    }
    if (pos > X_MAX) { pos = X_MAX; if (vel > 0) vel = 0; }
    const xwMin = wallXNow() + BALL_R + 0.15;
    if (pos < xwMin) { pos = xwMin; if (vel < 0) vel = 0; }
  }

  function init() {
    pos = 0; vel = 0; simTime = 0;
    resetMeasurement(true);
  }

  /* ── three.js scene ── */
  let simCanvas, renderer, scene, camera;
  let massMesh = null, massGlow = null, springMesh = null, wallMesh = null;
  const LOOK = new (function(){ this.x = -0.5; })();
  const orbit = { dragging: false, lastX: 0, lastY: 0, theta: Math.PI / 2, phi: 1.30, radius: 8 };

  let _glowTexture = null;
  function glowTexture() {
    if (_glowTexture) return _glowTexture;
    const c = document.createElement('canvas');
    c.width = c.height = 64;
    const g = c.getContext('2d');
    const grad = g.createRadialGradient(32, 32, 0, 32, 32, 32);
    grad.addColorStop(0,   'rgba(255,255,255,1)');
    grad.addColorStop(0.4, 'rgba(255,255,255,0.45)');
    grad.addColorStop(1,   'rgba(255,255,255,0)');
    g.fillStyle = grad;
    g.fillRect(0, 0, 64, 64);
    _glowTexture = new THREE.CanvasTexture(c);
    return _glowTexture;
  }

  function clearScene() {
    if (!scene) return;
    scene.traverse(obj => {
      if (obj.geometry) obj.geometry.dispose();
      if (obj.material) obj.material.dispose();
    });
    while (scene.children.length) scene.remove(scene.children[0]);
    massMesh = massGlow = springMesh = wallMesh = null;
  }

  function buildScene() {
    clearScene();
    scene.add(new THREE.AmbientLight(0xffffff, 0.5));
    const dl = new THREE.DirectionalLight(0xffffff, 0.7);
    dl.position.set(4, 7, 5);
    scene.add(dl);

    /* Wall: glowing teal bar (vibrates when driven) */
    const tealL = _col('--teal-light');
    const wallMat = new THREE.MeshPhongMaterial({ color: tealL, emissive: tealL, emissiveIntensity: 0.85 });
    wallMesh = new THREE.Mesh(new THREE.BoxGeometry(0.25, 2.6, 2.6), wallMat);
    wallMesh.position.set(-HALF - 0.125, 0, 0);
    scene.add(wallMesh);

    /* Mass: white-hot sphere + additive teal halo */
    const mat = new THREE.MeshPhongMaterial({ color: 0xffffff, emissive: 0xffffff, emissiveIntensity: 1.0, shininess: 80 });
    massMesh = new THREE.Mesh(new THREE.SphereGeometry(BALL_R, 24, 18), mat);
    massMesh.position.set(0, 0, 0);
    scene.add(massMesh);
    const glowMat = new THREE.SpriteMaterial({ map: glowTexture(), color: _col('--teal-light'), transparent: true, depthWrite: false, blending: THREE.AdditiveBlending });
    massGlow = new THREE.Sprite(glowMat);
    massGlow.scale.set(BALL_R * 5.5, BALL_R * 5.5, 1);
    massGlow.position.copy(massMesh.position);
    scene.add(massGlow);

    /* Spring: neon cylinder (width ∝ ω₀), recolored by extension each frame */
    const sr = springRad();
    const sprMat = new THREE.MeshPhongMaterial({ color: _col('--teal-dark'), emissive: _col('--teal-dark'), emissiveIntensity: 0.8 });
    springMesh = new THREE.Mesh(new THREE.CylinderGeometry(sr, sr, 1, 8, 1), sprMat);
    springMesh.rotation.z = Math.PI / 2;
    scene.add(springMesh);

    updateScene();
  }

  let _cTealD, _cTealL, _cPinkD, _cPinkL, _cMid, _spring = null;
  function springColor(ext, scale) {
    if (!_cTealD) {
      _cTealD = _col('--teal-dark');  _cTealL = _col('--teal-light');
      _cPinkD = _col('--pink-dark');  _cPinkL = _col('--pink-light');
      _cMid   = _cTealD.clone().lerp(_cPinkD, 0.5);
      _spring = new THREE.Color();
    }
    const t = Math.max(0, Math.min(1, Math.abs(ext) / scale));
    const dark  = ext >= 0 ? _cTealD : _cPinkD;
    const light = ext >= 0 ? _cTealL : _cPinkL;
    if (t < 0.5) _spring.copy(_cMid).lerp(dark, t / 0.5);
    else         _spring.copy(dark).lerp(light, (t - 0.5) / 0.5);
    return _spring;
  }

  function updateScene() {
    const xw = wallXNow();
    wallMesh.position.x = xw - 0.125;
    massMesh.position.x = pos;
    massGlow.position.x = pos;
    /* spring spans wall face → mass center */
    const xa = xw, xb = pos;
    const len = Math.max(Math.abs(xb - xa), 1e-4);
    springMesh.scale.y = len;
    springMesh.position.x = (xa + xb) / 2;
    /* extension relative to rest length HALF */
    const col = springColor((xb - xa) - HALF, 1.5);
    springMesh.material.color.copy(col);
    springMesh.material.emissive.copy(col);
  }

  function initThree() {
    if (!renderer) {
      renderer = new THREE.WebGLRenderer({ canvas: simCanvas, antialias: true });
      renderer.setPixelRatio(window.devicePixelRatio);
      renderer.setClearColor(_col('--bg-void'), 1);
      camera = new THREE.PerspectiveCamera(50, 1, 0.1, 100);
      scene = new THREE.Scene();
      setupPointerControls();
      updateCamera();
    }
    buildScene();
    resizeRenderer();
  }

  function updateCamera() {
    const { theta, phi, radius } = orbit;
    camera.position.set(LOOK.x + radius*Math.sin(phi)*Math.cos(theta), radius*Math.cos(phi), radius*Math.sin(phi)*Math.sin(theta));
    camera.lookAt(LOOK.x, 0, 0);
  }

  function resizeRenderer() {
    if (!renderer) return;
    const w = simCanvas.clientWidth || simCanvas.width, h = simCanvas.clientHeight || simCanvas.height;
    renderer.setSize(w, h, false);
    camera.aspect = w / h;
    camera.updateProjectionMatrix();
  }

  /* ── Pointer: drag the mass if hit, otherwise orbit; wheel zooms (down = in) ── */
  let draggingMass = false;
  const _ray = { caster: null, ndc: null, plane: null, hit: null };

  function rayFromEvent(e) {
    if (!_ray.caster) {
      _ray.caster = new THREE.Raycaster();
      _ray.ndc    = new THREE.Vector2();
      _ray.plane  = new THREE.Plane();
      _ray.hit    = new THREE.Vector3();
    }
    const rect = simCanvas.getBoundingClientRect();
    _ray.ndc.x = ((e.clientX - rect.left) / rect.width) * 2 - 1;
    _ray.ndc.y = -((e.clientY - rect.top) / rect.height) * 2 + 1;
    _ray.caster.setFromCamera(_ray.ndc, camera);
    return _ray.caster;
  }

  function setupPointerControls() {
    simCanvas.addEventListener('pointerdown', e => {
      const caster = rayFromEvent(e);
      const hits = massMesh ? caster.intersectObject(massMesh) : [];
      if (hits.length > 0) {
        draggingMass = true;
        vel = 0;
        resetMeasurement(false);
      } else {
        orbit.dragging = true;
        orbit.lastX = e.clientX; orbit.lastY = e.clientY;
      }
      simCanvas.setPointerCapture(e.pointerId);
    });
    simCanvas.addEventListener('pointermove', e => {
      if (draggingMass) {
        const caster = rayFromEvent(e);
        const n = camera.getWorldDirection(new THREE.Vector3());
        _ray.plane.setFromNormalAndCoplanarPoint(n, new THREE.Vector3(LOOK.x, 0, 0));
        if (caster.ray.intersectPlane(_ray.plane, _ray.hit)) {
          pos = Math.max(wallXNow() + BALL_R + 0.3, Math.min(X_MAX, _ray.hit.x));
          vel = 0;
          if (!running) { updateScene(); renderer.render(scene, camera); }
        }
      } else if (orbit.dragging) {
        orbit.theta -= (e.clientX - orbit.lastX) * 0.008;
        orbit.phi = Math.max(0.15, Math.min(Math.PI - 0.15, orbit.phi + (e.clientY - orbit.lastY) * 0.008));
        orbit.lastX = e.clientX; orbit.lastY = e.clientY;
        updateCamera();
        if (!running) renderer.render(scene, camera);
      }
    });
    const release = () => { draggingMass = false; orbit.dragging = false; };
    simCanvas.addEventListener('pointerup', release);
    simCanvas.addEventListener('pointercancel', release);
    simCanvas.addEventListener('wheel', e => {
      e.preventDefault();
      /* scroll down = zoom in */
      orbit.radius = Math.max(3, Math.min(22, orbit.radius - e.deltaY * 0.01));
      updateCamera();
      if (!running) renderer.render(scene, camera);
    }, { passive: false });
  }

  /* ── Resonance-response plot: gain G(ω_d) = X_steady / A ──
     Theory curve for the current ω₀, γ (its peak ≈ Q = ω₀/γ) plus live
     measured dots: after each change the sim settles, then max|x| is
     measured over a few drive periods and plotted at the current ω_d. */
  let qCanvas, qCtx;
  const W_PLOT_MAX = 4;                  // matches the drive slider range
  const measured = new Map();            // round(ω·20) → measured gain
  const T_SETTLE = 12, T_MEAS_MIN = 8;   // sim-time seconds
  let measClock = 0, measuring = false, measAmp = 0;

  function resetMeasurement(clearAll) {
    if (clearAll) measured.clear();
    measClock = 0; measAmp = 0; measuring = false;
  }

  function updateMeasurement() {
    if (driveW <= 0 || draggingMass) return;
    measClock += DT;
    if (!measuring) {
      if (measClock >= T_SETTLE) { measuring = true; measAmp = 0; measClock = 0; }
    } else {
      if (Math.abs(pos) > measAmp) measAmp = Math.abs(pos);
      const T_MEAS = Math.max(T_MEAS_MIN, 3 * 2 * Math.PI / driveW);
      if (measClock >= T_MEAS) {
        measured.set(Math.round(driveW * 20), measAmp / DRIVE_A);
        measAmp = 0; measClock = 0;   // keep re-measuring, updating the dot
      }
    }
  }

  function gainTheory(w) {
    const d = omega0 * omega0 - w * w;
    const den = Math.sqrt(d * d + gammaD * gammaD * w * w);
    return den > 1e-9 ? (omega0 * omega0) / den : 1e9;
  }

  function resizeQCanvas() {
    if (!qCanvas) return;
    const sect = qCanvas.closest('#sho-qplot-section');
    if (!sect) return;
    const w = sect.clientWidth - 24, h = sect.clientHeight - 30;
    if (w > 0 && h > 0 && (qCanvas.width !== w || qCanvas.height !== h)) {
      qCanvas.width = w; qCanvas.height = h;
    }
  }

  function renderQPlot() {
    if (!qCtx) return;
    resizeQCanvas();
    const W = qCanvas.width, H = qCanvas.height;
    if (W === 0 || H === 0) return;
    qCtx.fillStyle = _c('--bg-dark');
    qCtx.fillRect(0, 0, W, H);
    /* readable font + real margins around the plot area */
    const fs = Math.max(14, Math.round(H * 0.075));
    const PL = Math.round(fs * 2.2), PR = Math.round(fs * 1.2);
    const PT = Math.round(fs * 1.6), PB = Math.round(fs * 2.0);
    const pw = W - PL - PR, ph = H - PT - PB;

    /* y-scale from theory peak (capped) and any measured overshoot */
    let gMax = 1;
    for (let s = 0; s <= 100; s++) {
      const g = gainTheory((s / 100) * W_PLOT_MAX);
      if (g > gMax) gMax = g;
    }
    for (const g of measured.values()) if (g > gMax) gMax = g;
    gMax = Math.min(40, gMax) * 1.15;

    const toX = w => PL + (w / W_PLOT_MAX) * pw;
    const toY = g => PT + ph * (1 - Math.min(g, gMax) / gMax);

    /* axes */
    qCtx.strokeStyle = _rgba('--text-dim', 0.35);
    qCtx.lineWidth = 1;
    qCtx.beginPath();
    qCtx.moveTo(PL, PT); qCtx.lineTo(PL, PT + ph); qCtx.lineTo(PL + pw, PT + ph);
    qCtx.stroke();
    /* G = 1 reference */
    qCtx.setLineDash([3, 5]);
    qCtx.beginPath();
    qCtx.moveTo(PL, toY(1)); qCtx.lineTo(PL + pw, toY(1));
    qCtx.stroke();
    qCtx.setLineDash([]);
    qCtx.fillStyle = _rgba('--text-mid', 0.9);
    qCtx.font = `${fs}px 'EB Garamond', Georgia, serif`;
    qCtx.textAlign = 'center'; qCtx.textBaseline = 'top';
    qCtx.fillText('ω_d', PL + pw / 2, PT + ph + Math.round(fs * 0.35));
    qCtx.save();
    qCtx.translate(Math.round(fs * 0.9), PT + ph / 2); qCtx.rotate(-Math.PI / 2);
    qCtx.textAlign = 'center';
    qCtx.textBaseline = 'middle';
    qCtx.fillText('gain', 0, 0);
    qCtx.restore();
    if (gammaD > 0.005) {
      qCtx.textAlign = 'right';
      qCtx.textBaseline = 'top';
      qCtx.fillStyle = _rgba('--teal-light', 0.9);
      qCtx.fillText('Q ≈ ' + (omega0 / gammaD).toFixed(1), W - PR, Math.round(fs * 0.25));
    }

    /* theory curve — neon */
    qCtx.save();
    qCtx.shadowColor = _c('--teal-dark');
    qCtx.shadowBlur  = 10;
    qCtx.strokeStyle = _c('--teal-light');
    qCtx.lineWidth   = 1.6;
    qCtx.beginPath();
    for (let s = 0; s <= 140; s++) {
      const w = (s / 140) * W_PLOT_MAX;
      const x = toX(w), y = toY(gainTheory(w));
      s === 0 ? qCtx.moveTo(x, y) : qCtx.lineTo(x, y);
    }
    qCtx.stroke();
    qCtx.restore();

    /* current drive frequency marker */
    if (driveW > 0) {
      qCtx.save();
      qCtx.setLineDash([4, 5]);
      qCtx.strokeStyle = _rgba('--pink-light', 0.5);
      qCtx.beginPath();
      qCtx.moveTo(toX(driveW), PT); qCtx.lineTo(toX(driveW), PT + ph);
      qCtx.stroke();
      qCtx.restore();
    }

    /* measured dots — glowing pink */
    qCtx.save();
    qCtx.shadowColor = _c('--pink-dark');
    qCtx.shadowBlur  = 8;
    qCtx.fillStyle   = _c('--pink-light');
    for (const [key, g] of measured) {
      const w = key / 20;
      qCtx.beginPath();
      qCtx.arc(toX(w), toY(g), 2.6, 0, Math.PI * 2);
      qCtx.fill();
    }
    qCtx.restore();
  }

  /* ── Loop ── */
  let plotSkip = 0;
  function loop() {
    if (running) {
      if (!draggingMass) step();
      else simTime += DT;              // wall keeps vibrating while mass is held
      updateMeasurement();
      updateScene();
      renderer.render(scene, camera);
      if (++plotSkip >= 3) { plotSkip = 0; renderQPlot(); }
    }
    frameId = requestAnimationFrame(loop);
  }

  /* ── Shell wiring ── */
  const shell = new AppletShell({
    id:    'sho',
    title: 'SHO &mdash; Damped &amp; Driven',
    gap:   0,

    headerBtns: `<button class="applet-shell-header-btn" onclick="shoReset()">Reset</button><button class="applet-shell-header-btn" id="sho-pause-btn" onclick="shoTogglePause()">Pause</button>`,

    ctrlHTML: `
      <div class="applet-shell-ctrl-section">
        <div class="applet-shell-ctrl-title">Natural frequency &omega;&#8320;</div>
        <div class="applet-shell-slider-row">
          <span class="applet-shell-side">Low</span>
          <input type="range" id="sho-omega0" min="0.2" max="4" step="0.05" value="1.0">
          <span class="applet-shell-side">High</span>
          <span class="applet-shell-val" id="sho-omega0-val">1.00</span>
        </div>
      </div>
      <div class="applet-shell-ctrl-section">
        <div class="applet-shell-ctrl-title">Damping &gamma;</div>
        <div class="applet-shell-slider-row">
          <span class="applet-shell-side">None</span>
          <input type="range" id="sho-damping" min="0" max="2" step="0.01" value="0">
          <span class="applet-shell-side">Heavy</span>
          <span class="applet-shell-val" id="sho-damping-val">0.00</span>
        </div>
      </div>
      <div class="applet-shell-ctrl-section">
        <div class="applet-shell-ctrl-title">Driving frequency &omega;<sub>d</sub></div>
        <div class="applet-shell-slider-row">
          <span class="applet-shell-side">Off</span>
          <input type="range" id="sho-drive" min="0" max="4" step="0.05" value="0">
          <span class="applet-shell-side">Fast</span>
          <span class="applet-shell-val" id="sho-drive-val">Off</span>
        </div>
      </div>
      <div id="sho-qplot-section">
        <div class="applet-shell-ctrl-title">Resonance response</div>
        <canvas id="sho-qplot"></canvas>
      </div>
    `,

    onOpen: function ({ canvas: c }) {
      simCanvas = c;
      init();
      const pb = document.getElementById('sho-pause-btn');
      if (pb) { pb.textContent = 'Pause'; pb.classList.remove('active'); }
      function start() {
        setTimeout(() => {
          qCanvas = document.getElementById('sho-qplot');
          qCtx    = qCanvas.getContext('2d');
          initThree();
          renderQPlot();
          running = true;
          if (!frameId) frameId = requestAnimationFrame(loop);
        }, 80);
      }
      if (window.THREE) {
        start();
      } else {
        const s = document.createElement('script');
        s.src = 'https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js';
        s.onload = start;
        document.head.appendChild(s);
      }
    },

    onClose: function () {
      running = false;
      if (frameId) { cancelAnimationFrame(frameId); frameId = null; }
      const pb = document.getElementById('sho-pause-btn');
      if (pb) { pb.textContent = 'Pause'; pb.classList.remove('active'); }
    },

    onResize: function () {
      resizeRenderer();
      if (renderer && !running) renderer.render(scene, camera);
    },
  });

  window.shoOpen  = () => shell.open();
  window.shoClose = () => shell.close();

  window.shoReset = function () {
    init();
    if (scene) { updateScene(); if (!running) renderer.render(scene, camera); }
  };

  window.shoTogglePause = function () {
    running = !running;
    const pb = document.getElementById('sho-pause-btn');
    if (pb) {
      pb.textContent = running ? 'Pause' : 'Resume';
      pb.classList.toggle('active', !running);
    }
  };

  /* ── Sliders ── */
  document.getElementById('sho-omega0').addEventListener('input', function () {
    omega0 = parseFloat(this.value);
    document.getElementById('sho-omega0-val').textContent = omega0.toFixed(2);
    resetMeasurement(true);       // theory curve changed — old dots invalid
    if (scene) buildScene();      // spring width follows ω₀
    renderQPlot();
  });
  document.getElementById('sho-damping').addEventListener('input', function () {
    gammaD = parseFloat(this.value);
    document.getElementById('sho-damping-val').textContent = gammaD.toFixed(2);
    resetMeasurement(true);       // theory curve changed — old dots invalid
    renderQPlot();
  });
  document.getElementById('sho-drive').addEventListener('input', function () {
    driveW = parseFloat(this.value);
    document.getElementById('sho-drive-val').textContent = driveW > 0 ? driveW.toFixed(2) : 'Off';
    resetMeasurement(false);      // keep dots from other frequencies
    renderQPlot();
  });

})();
