(function () {
  'use strict';

  const _cs = getComputedStyle(document.documentElement);
  const _c   = n => _cs.getPropertyValue(n).trim();
  const _rgb  = n => { const h = _c(n).replace('#',''); const v = parseInt(h,16); return [(v>>16)&0xFF,(v>>8)&0xFF,v&0xFF]; };
  const _rgba = (n, a) => { const [r,g,b] = _rgb(n); return `rgba(${r},${g},${b},${a})`; };
  const _col  = n => { const [r,g,b] = _rgb(n); return new THREE.Color(r/255, g/255, b/255); };

  const [_TLR, _TLG, _TLB] = _rgb('--teal-light');
  const [_PLR, _PLG, _PLB] = _rgb('--pink-light');

  /* ── Inject CSS (phase-plot section fills the control panel) ── */
  (function () {
    if (document.getElementById('co-styles')) return;
    const s = document.createElement('style');
    s.id = 'co-styles';
    s.textContent = `
      #co-ctrl-panel { display:flex; flex-direction:column; overflow:hidden; }
      #co-phase-section { flex:1; min-height:0; display:flex; flex-direction:column; padding:6px 12px 10px; }
      #co-phase { flex:1; min-height:0; width:100%; display:block; }
    `;
    document.head.appendChild(s);
  })();

  /* ── Parameters ── */
  let numMasses = 5;
  let kspring   = 1.0;
  let mass      = 1.0;
  let damping   = 0.0;

  /* ── State ── */
  let lPos = null, lVel = null;
  let running = false, frameId = null;
  let wasRunning = false;   // sim state stashed while the docs panel is open

  /* ── World geometry: chain along x between two walls ── */
  const SPAN = 8;                                  // wall inner face to wall inner face
  const HALF = SPAN / 2;
  function pitch()     { return SPAN / (numMasses + 1); }
  function ballR()     {
    /* radius grows with mass like a solid ball, capped so neighbors can't touch at rest */
    const base = Math.min(0.42, pitch() * 0.30);
    return Math.max(0.08, Math.min(pitch() * 0.45, base * Math.cbrt(mass)));
  }
  function springRad() { return Math.max(0.015, Math.min(0.10, 0.04 * Math.pow(kspring, 0.45))); }
  function dispScale() { return pitch() * 0.5; }   // world units per displacement unit
  function eqX(i)      { return -HALF + (i + 1) * pitch(); }

  /* ── Physics (symplectic Euler + hard-core collisions, clamped walls) ── */
  const DT = 0.075, SUBSTEPS = 8;
  const dt = DT / SUBSTEPS;

  function step() {
    const N = numMasses;
    const dscl   = dispScale();
    const minSep = (2 * ballR()) / dscl;
    const eqSep  = pitch() / dscl;
    for (let s = 0; s < SUBSTEPS; s++) {
      for (let i = 0; i < N; i++) {
        const left  = i > 0     ? lPos[i-1] : 0;
        const right = i < N - 1 ? lPos[i+1] : 0;
        const F = kspring * (right - lPos[i]) - kspring * (lPos[i] - left);
        const a = F / mass - (damping / mass) * lVel[i];
        lVel[i] += a * dt;
      }
      for (let i = 0; i < N; i++) lPos[i] += lVel[i] * dt;
      /* hard-core collisions between adjacent masses */
      for (let pass = 0; pass < 2; pass++) {
        const i0 = pass === 0 ? 0 : N - 2;
        const di = pass === 0 ? 1 : -1;
        for (let ii = 0; ii < N - 1; ii++) {
          const i = i0 + ii * di, j = i + 1;
          const sep = eqSep + (lPos[j] - lPos[i]);
          if (sep < minSep) {
            const ov = minSep - sep;
            lPos[i] -= ov * 0.5;
            lPos[j] += ov * 0.5;
            const vi = lVel[i], vj = lVel[j];
            if (vi > vj) { lVel[i] = vj; lVel[j] = vi; }
          }
        }
      }
    }
  }

  function init() {
    lPos = new Float64Array(numMasses);
    lVel = new Float64Array(numMasses);
    phaseHist = [];
    resetPhaseLimits();
  }

  /* ── three.js scene ── */
  let simCanvas, renderer, scene, camera;
  let massMeshes = [], massGlows = [], springMeshes = [], wallMeshes = [];
  const orbit = { dragging: false, lastX: 0, lastY: 0, theta: Math.PI / 2, phi: 1.30, radius: 9 };

  /* Soft radial sprite texture for halos (as in the Kuramoto applet) */
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

  function chainColor(i) {
    const t = numMasses > 1 ? i / (numMasses - 1) : 0;
    return new THREE.Color(
      (_TLR + (_PLR - _TLR) * t) / 255,
      (_TLG + (_PLG - _TLG) * t) / 255,
      (_TLB + (_PLB - _TLB) * t) / 255
    );
  }

  function clearScene() {
    if (!scene) return;
    scene.traverse(obj => {
      if (obj.geometry) obj.geometry.dispose();
      if (obj.material) obj.material.dispose();
    });
    while (scene.children.length) scene.remove(scene.children[0]);
    massMeshes = []; massGlows = []; springMeshes = []; wallMeshes = [];
  }

  function buildScene() {
    clearScene();
    scene.add(new THREE.AmbientLight(0xffffff, 0.5));
    const dl = new THREE.DirectionalLight(0xffffff, 0.7);
    dl.position.set(4, 7, 5);
    scene.add(dl);

    /* Walls: glowing teal bars */
    const wallGeo = new THREE.BoxGeometry(0.25, 2.6, 2.6);
    const tealL = _col('--teal-light');
    for (const sx of [-1, 1]) {
      const mat = new THREE.MeshPhongMaterial({ color: tealL, emissive: tealL, emissiveIntensity: 0.85 });
      const wall = new THREE.Mesh(wallGeo.clone(), mat);
      wall.position.set(sx * (HALF + 0.125), 0, 0);
      scene.add(wall); wallMeshes.push(wall);
    }

    /* Masses: white-hot spheres + additive chain-gradient halos */
    const r = ballR();
    const sphGeo = new THREE.SphereGeometry(r, 20, 16);
    for (let i = 0; i < numMasses; i++) {
      const mat = new THREE.MeshPhongMaterial({ color: 0xffffff, emissive: 0xffffff, emissiveIntensity: 1.0, shininess: 80 });
      const mesh = new THREE.Mesh(sphGeo.clone(), mat);
      mesh.position.set(eqX(i), 0, 0);
      mesh.userData.massIndex = i;
      scene.add(mesh); massMeshes.push(mesh);
      const glowMat = new THREE.SpriteMaterial({ map: glowTexture(), color: chainColor(i), transparent: true, depthWrite: false, blending: THREE.AdditiveBlending });
      const glow = new THREE.Sprite(glowMat);
      glow.scale.set(r * 5.5, r * 5.5, 1);
      glow.position.copy(mesh.position);
      scene.add(glow); massGlows.push(glow);
    }

    /* Springs: N+1 neon cylinders (width ∝ k^0.45), recolored by extension each frame */
    const sr = springRad();
    const sprGeo = new THREE.CylinderGeometry(sr, sr, 1, 8, 1);
    for (let s = 0; s <= numMasses; s++) {
      const mat = new THREE.MeshPhongMaterial({ color: _col('--teal-dark'), emissive: _col('--teal-dark'), emissiveIntensity: 0.8 });
      const cyl = new THREE.Mesh(sprGeo.clone(), mat);
      scene.add(cyl); springMeshes.push(cyl);
    }
    updateScene();
  }

  function positionCylinderX(cyl, xa, xb) {
    const len = Math.max(Math.abs(xb - xa), 1e-4);
    cyl.scale.y = len;
    cyl.position.set((xa + xb) / 2, 0, 0);
    cyl.rotation.z = Math.PI / 2;
  }

  let _cTealD, _cTealL, _cPinkD, _cPinkL, _cMid, _spring = null;
  function springColor(ext, scale) {
    if (!_cTealD) {
      _cTealD = _col('--teal-dark');  _cTealL = _col('--teal-light');
      _cPinkD = _col('--pink-dark');  _cPinkL = _col('--pink-light');
      _cMid   = _cTealD.clone().lerp(_cPinkD, 0.5);   // neutral dark at zero extension
      _spring = new THREE.Color();
    }
    /* Continuous ramp: pink-light ← pink-dark ← neutral → teal-dark → teal-light,
       so the color no longer jumps as a spring crosses zero extension. */
    const t = Math.max(0, Math.min(1, Math.abs(ext) / scale));
    const dark  = ext >= 0 ? _cTealD : _cPinkD;
    const light = ext >= 0 ? _cTealL : _cPinkL;
    if (t < 0.5) _spring.copy(_cMid).lerp(dark, t / 0.5);
    else         _spring.copy(dark).lerp(light, (t - 0.5) / 0.5);
    return _spring;
  }

  function updateScene() {
    const dscl = dispScale();
    const N = numMasses;
    for (let i = 0; i < N; i++) {
      const x = eqX(i) + lPos[i] * dscl;
      massMeshes[i].position.x = x;
      massGlows[i].position.x  = x;
    }
    /* springs: wall→m0, m_i→m_{i+1}, m_{N-1}→wall */
    const extScale = 1.0;   // displacement units at which spring color saturates
    for (let s = 0; s <= N; s++) {
      const xa = s === 0 ? -HALF : massMeshes[s-1].position.x;
      const xb = s === N ?  HALF : massMeshes[s].position.x;
      positionCylinderX(springMeshes[s], xa, xb);
      const left  = s === 0 ? 0 : lPos[s-1];
      const right = s === N ? 0 : lPos[s];
      const col = springColor(right - left, extScale);
      springMeshes[s].material.color.copy(col);
      springMeshes[s].material.emissive.copy(col);
    }
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
    camera.position.set(radius*Math.sin(phi)*Math.cos(theta), radius*Math.cos(phi), radius*Math.sin(phi)*Math.sin(theta));
    camera.lookAt(0, 0, 0);
  }

  function resizeRenderer() {
    if (!renderer) return;
    const w = simCanvas.clientWidth || simCanvas.width, h = simCanvas.clientHeight || simCanvas.height;
    renderer.setSize(w, h, false);
    camera.aspect = w / h;
    camera.updateProjectionMatrix();
  }

  /* ── Pointer: drag a mass if hit, otherwise orbit; wheel zooms (down = in) ── */
  let dragIdx = -1;
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
      const hits = caster.intersectObjects(massMeshes);
      if (hits.length > 0) {
        dragIdx = hits[0].object.userData.massIndex;
        lVel[dragIdx] = 0;
        resetPhaseLimits();
      } else {
        orbit.dragging = true;
        orbit.lastX = e.clientX; orbit.lastY = e.clientY;
      }
      simCanvas.setPointerCapture(e.pointerId);
    });
    simCanvas.addEventListener('pointermove', e => {
      if (dragIdx >= 0) {
        const caster = rayFromEvent(e);
        /* drag in the camera-facing plane containing the chain axis */
        const n = camera.getWorldDirection(new THREE.Vector3());
        _ray.plane.setFromNormalAndCoplanarPoint(n, new THREE.Vector3(0, 0, 0));
        if (caster.ray.intersectPlane(_ray.plane, _ray.hit)) {
          const dscl = dispScale();
          let d = (_ray.hit.x - eqX(dragIdx)) / dscl;
          /* keep dragged mass between its neighbors' current positions (minus contact) */
          const minSep = (2 * ballR()) / dscl, eqSep = pitch() / dscl;
          if (dragIdx > 0)             d = Math.max(d, lPos[dragIdx-1] - eqSep + minSep);
          else                         d = Math.max(d, -eqSep + minSep);
          if (dragIdx < numMasses - 1) d = Math.min(d, lPos[dragIdx+1] + eqSep - minSep);
          else                         d = Math.min(d, eqSep - minSep);
          lPos[dragIdx] = d;
          lVel[dragIdx] = 0;
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
    const release = () => { dragIdx = -1; orbit.dragging = false; };
    simCanvas.addEventListener('pointerup', release);
    simCanvas.addEventListener('pointercancel', release);
    simCanvas.addEventListener('wheel', e => {
      e.preventDefault();
      /* scroll down = zoom in */
      orbit.radius = Math.max(3, Math.min(25, orbit.radius - e.deltaY * 0.01));
      updateCamera();
      if (!running) renderer.render(scene, camera);
    }, { passive: false });
  }

  /* ── Phase-space plot (control panel) ── */
  let phaseCanvas, phaseCtx;
  let phaseHist = [];
  const PHASE_HIST = 90;

  /* Axis limits ratchet up to the max encountered; they reset to defaults
     only when the user changes something (sliders, reset, grabbing a mass). */
  const PHASE_LIM_DEFAULT = 0.2;
  let phaseMaxX = PHASE_LIM_DEFAULT, phaseMaxP = PHASE_LIM_DEFAULT;

  function resetPhaseLimits() {
    phaseMaxX = PHASE_LIM_DEFAULT;
    phaseMaxP = PHASE_LIM_DEFAULT;
  }

  function recordPhase() {
    const frame = new Float64Array(2 * numMasses);
    for (let i = 0; i < numMasses; i++) {
      frame[2*i]   = lPos[i];
      frame[2*i+1] = mass * lVel[i];
      const ax = Math.abs(frame[2*i]), ap = Math.abs(frame[2*i+1]);
      if (ax > phaseMaxX) phaseMaxX = ax;
      if (ap > phaseMaxP) phaseMaxP = ap;
    }
    phaseHist.push(frame);
    if (phaseHist.length > PHASE_HIST) phaseHist.shift();
  }

  /* Cached phase-blob sprites: chain gradient, rebuilt when radius changes */
  const N_BLOB_C = 16;
  let blobSprites = null, blobSpriteR = -1;

  function buildBlobSprites(blobR) {
    blobSprites = new Array(N_BLOB_C);
    blobSpriteR = blobR;
    const R3 = blobR * 3;
    const size = Math.ceil(2 * R3) + 2;
    for (let ci = 0; ci < N_BLOB_C; ci++) {
      const t = ci / (N_BLOB_C - 1);
      const r = Math.round(_TLR + (_PLR - _TLR) * t);
      const g = Math.round(_TLG + (_PLG - _TLG) * t);
      const b = Math.round(_TLB + (_PLB - _TLB) * t);
      const cnv = document.createElement('canvas');
      cnv.width = cnv.height = size;
      const c2 = cnv.getContext('2d');
      const cc = size / 2;
      const grd = c2.createRadialGradient(cc, cc, 0, cc, cc, R3);
      grd.addColorStop(0,    `rgba(${r},${g},${b},1)`);
      grd.addColorStop(0.35, `rgba(${r},${g},${b},0.4)`);
      grd.addColorStop(1,    `rgba(${r},${g},${b},0)`);
      c2.beginPath(); c2.arc(cc, cc, R3, 0, Math.PI * 2);
      c2.fillStyle = grd; c2.fill();
      blobSprites[ci] = cnv;
    }
  }

  function resizePhaseCanvas() {
    if (!phaseCanvas) return;
    const sect = phaseCanvas.closest('#co-phase-section');
    if (!sect) return;
    const w = sect.clientWidth - 24, h = sect.clientHeight - 30;
    if (w > 0 && h > 0 && (phaseCanvas.width !== w || phaseCanvas.height !== h)) {
      phaseCanvas.width = w; phaseCanvas.height = h;
    }
  }

  function renderPhase() {
    if (!phaseCtx) return;
    resizePhaseCanvas();
    const W = phaseCanvas.width, H = phaseCanvas.height;
    if (W === 0 || H === 0) return;
    phaseCtx.fillStyle = _c('--bg-dark');
    phaseCtx.fillRect(0, 0, W, H);
    const cx = W / 2, cy = H / 2;

    /* axes */
    phaseCtx.strokeStyle = _rgba('--text-dim', 0.35);
    phaseCtx.lineWidth = 1;
    phaseCtx.beginPath();
    phaseCtx.moveTo(6, cy); phaseCtx.lineTo(W - 6, cy);
    phaseCtx.moveTo(cx, 6); phaseCtx.lineTo(cx, H - 6);
    phaseCtx.stroke();
    phaseCtx.fillStyle = _rgba('--text-dim', 0.6);
    phaseCtx.font = `11px 'EB Garamond', Georgia, serif`;
    phaseCtx.textAlign = 'right'; phaseCtx.textBaseline = 'top';
    phaseCtx.fillText('x', W - 8, cy + 4);
    phaseCtx.textAlign = 'left';
    phaseCtx.fillText('p', cx + 5, 8);

    if (phaseHist.length === 0) return;

    const maxX = phaseMaxX, maxP = phaseMaxP;

    const nM = numMasses;
    const baseAlpha = Math.max(0.10, Math.min(0.55, 0.55 / Math.sqrt(nM)));
    const blobR = Math.max(1.4, Math.min(2.2, W * 0.011));
    if (blobSpriteR !== blobR) buildBlobSprites(blobR);
    const half = blobSprites[0].width / 2;

    const nFrames = phaseHist.length;
    for (let fi = 0; fi < nFrames; fi++) {
      const fr = phaseHist[fi];
      const ageFrac = (fi + 1) / nFrames;
      phaseCtx.globalAlpha = Math.min(1, baseAlpha * ageFrac * ageFrac * 2.5);
      for (let mi = 0; mi < nM; mi++) {
        const bx = cx + (fr[2*mi]   / maxX) * (W / 2 - 12);
        const by = cy - (fr[2*mi+1] / maxP) * (H / 2 - 12);
        const ci = nM === 1 ? 0 : ((mi / (nM - 1)) * (N_BLOB_C - 1) + 0.5) | 0;
        phaseCtx.drawImage(blobSprites[ci], bx - half, by - half);
      }
    }
    phaseCtx.globalAlpha = 1;
  }

  /* ── Loop ── */
  let phaseSkip = 0;
  function loop() {
    if (running) {
      step();
      if (dragIdx >= 0) lVel[dragIdx] = 0;   // pinned while held
      updateScene();
      renderer.render(scene, camera);
      recordPhase();
      if (++phaseSkip >= 2) { phaseSkip = 0; renderPhase(); }
    }
    frameId = requestAnimationFrame(loop);
  }

  /* ── Shell wiring ── */
  const shell = new AppletShell({
    id:    'co',
    title: 'Coupled Oscillators',
    gap:   0,

    headerBtns: `<button class="applet-shell-header-btn" onclick="coReset()">Reset</button><button class="applet-shell-header-btn" id="co-pause-btn" onclick="coTogglePause()">Pause</button>`,

    docs: {
      whatis: `Pull one mass aside and release it, and its motion does not stay put: it leaks along the chain, sloshing from mass to mass until, some time later, the disturbance returns. A chain of masses joined by springs between two walls is the simplest system where oscillators must negotiate, and it is the bridge between one particle and a continuum.¶The negotiation was settled by Daniel Bernoulli in 1753: however complicated the motion looks, it is a superposition of normal modes, collective patterns in which every mass oscillates at one shared frequency [bernoulli1753]. For $N$ equal masses $m$ and springs $k$, mode $n$ is a standing sine wave across the chain with frequency
$$\\omega_n = 2\\sqrt{\\frac{k}{m}} \\, \\sin\\!\\left( \\frac{n \\pi}{2(N+1)} \\right), \\qquad n = 1, \\dots, N.$$
Two masses give the classic beat phenomenon — energy sloshing wholly from one to the other and back — while fifteen begin to resemble a vibrating string, the mode frequencies crowding toward a maximum cutoff. This is how sound lives in a crystal.¶Two departures from the textbook keep things interesting here. The masses are hard spheres: displace them violently enough and neighbors collide, exchanging velocities in a way no linear theory contains: the same door to nonlinearity that Fermi, Pasta, and Ulam opened in 1955 when they asked how such a chain shares energy among its modes and got an answer that founded nonlinear science [fermi1955]. And the phase-space panel plots every mass as a point in the $(x, p)$ plane; a single undamped mass traces an ellipse, while the full chain weaves through its projections.`,

      howto: `Each mass glows with its own hue from teal to pink along the chain; springs tint teal when stretched, pink when compressed, with thickness following $k$ and ball size following $m$. Drag a mass to displace it (it pins while held), drag empty space to orbit, and scroll to zoom.¶Masses sets the chain length from a single oscillator to fifteen. Spring constant and Mass set $\\omega_0 = \\sqrt{k/m}$ for the whole chain, and Damping drains energy from every mass. Displace one end mass gently to launch a traveling pulse, or displace the middle one to watch symmetric modes only.¶The Phase space panel traces recent $(x, p)$ history for every mass, colored by chain position and fading with age: ellipses for clean modes, spirals under damping, and scribbles once collisions set in. Reset zeroes the chain; Pause freezes it.`,

      references: ['bernoulli1753', 'fermi1955'],
    },

    onDocsOpen:  function () { wasRunning = running; running = false; },
    onDocsClose: function () { running = wasRunning; },

    ctrlHTML: `
      <div class="applet-shell-ctrl-section">
        <div class="applet-shell-ctrl-title">Masses</div>
        <div class="applet-shell-slider-row">
          <span class="applet-shell-side">1</span>
          <input type="range" id="co-nmass" min="1" max="15" step="1" value="5">
          <span class="applet-shell-side">15</span>
          <span class="applet-shell-val" id="co-nmass-val">5</span>
        </div>
      </div>
      <div class="applet-shell-ctrl-section">
        <div class="applet-shell-ctrl-title">Spring constant k</div>
        <div class="applet-shell-slider-row">
          <span class="applet-shell-side">Soft</span>
          <input type="range" id="co-kspring" min="0.05" max="10" step="0.05" value="1.0">
          <span class="applet-shell-side">Stiff</span>
        </div>
      </div>
      <div class="applet-shell-ctrl-section">
        <div class="applet-shell-ctrl-title">Mass m</div>
        <div class="applet-shell-slider-row">
          <span class="applet-shell-side">Light</span>
          <input type="range" id="co-mass" min="0.05" max="10" step="0.05" value="1.0">
          <span class="applet-shell-side">Heavy</span>
        </div>
      </div>
      <div class="applet-shell-ctrl-section">
        <div class="applet-shell-ctrl-title">Damping</div>
        <div class="applet-shell-slider-row">
          <span class="applet-shell-side">None</span>
          <input type="range" id="co-damping" min="0" max="5" step="0.01" value="0">
          <span class="applet-shell-side">Heavy</span>
        </div>
      </div>
      <div id="co-phase-section">
        <div class="applet-shell-ctrl-title">Phase space</div>
        <canvas id="co-phase"></canvas>
      </div>
    `,

    onOpen: function ({ canvas: c }) {
      simCanvas = c;
      init();
      const pb = document.getElementById('co-pause-btn');
      if (pb) { pb.textContent = 'Pause'; pb.classList.remove('active'); }
      function start() {
        setTimeout(() => {
          phaseCanvas = document.getElementById('co-phase');
          phaseCtx    = phaseCanvas.getContext('2d');
          initThree();
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
      const pb = document.getElementById('co-pause-btn');
      if (pb) { pb.textContent = 'Pause'; pb.classList.remove('active'); }
    },

    onResize: function () {
      resizeRenderer();
      if (renderer && !running) renderer.render(scene, camera);
    },
  });

  window.coOpen  = () => shell.open();
  window.coClose = () => shell.close();

  window.coReset = function () {
    init();
    if (scene) { updateScene(); if (!running) renderer.render(scene, camera); }
  };

  window.coTogglePause = function () {
    running = !running;
    const pb = document.getElementById('co-pause-btn');
    if (pb) {
      pb.textContent = running ? 'Pause' : 'Resume';
      pb.classList.toggle('active', !running);
    }
  };

  /* ── Sliders ── */
  document.getElementById('co-nmass').addEventListener('input', function () {
    numMasses = parseInt(this.value);
    document.getElementById('co-nmass-val').textContent = numMasses;
    init();
    if (scene) buildScene();
  });
  document.getElementById('co-kspring').addEventListener('input', function () {
    kspring = parseFloat(this.value);
    resetPhaseLimits();
    if (scene) buildScene();      // spring width follows k
  });
  document.getElementById('co-mass').addEventListener('input', function () {
    mass = parseFloat(this.value);
    resetPhaseLimits();
    if (scene) buildScene();      // ball size follows m
  });
  document.getElementById('co-damping').addEventListener('input', function () {
    damping = parseFloat(this.value);
    resetPhaseLimits();
  });

})();
