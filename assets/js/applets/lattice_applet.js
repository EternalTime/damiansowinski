(function () {
  'use strict';

  const _cs = getComputedStyle(document.documentElement);
  const _c   = n => _cs.getPropertyValue(n).trim();
  const _rgb  = n => { const h = _c(n).replace('#',''); const v = parseInt(h,16); return [(v>>16)&0xFF,(v>>8)&0xFF,v&0xFF]; };

  const [_TDR, _TDG, _TDB] = _rgb('--teal-dark');
  const [_PDR, _PDG, _PDB] = _rgb('--pink-dark');

  /*
   * COORDINATE CONVENTION
   * ---------------------
   * The observer moves in the -Z world direction (Three.js camera default look).
   * Rest-frame lattice: integer multiples of SPACING in x, y, z.
   * Camera fixed at world origin; lattice streams past in +Z as observer moves forward.
   *
   * For a sphere at rest-frame offset (rx, ry, rz) relative to camera:
   *   rz < 0  →  ahead   (observer moving toward it)
   *   rz > 0  →  behind
   *
   * Lorentz contraction: oz = rz / gamma  (z-separations shrink)
   * Aberration: apparent direction bunches toward forward (-Z) at high beta
   * Doppler:    forward spheres blueshift, rear redshift
   */

  /* ── Grid parameters ── */
  const NXY     = 9;    // transverse half-extent (x and y)
  const NZ_F    = 40;   // planes ahead  (negative z)
  const NZ_B    = 6;    // planes behind (positive z)
  const SPACING = 3.0;  // rest-frame lattice spacing

  const PERIOD_XY = (2*NXY + 1) * SPACING;
  const PERIOD_Z  = (NZ_F + NZ_B + 1) * SPACING;

  /* Sphere half-offset so camera sits at body-centre of a cell (no sphere on any axis) */
  const OFF = SPACING * 0.5;

  let beta = 0.0;
  let offsetX = 0.0, offsetY = 0.0;  // camera translation within one cell, in [-SPACING/2, SPACING/2]
  let running = false, frameId = null;
  let renderer, scene, camera;
  let meshes = [];

  /* free-look */
  let yaw = 0, pitch = 0;
  let dragging = false, mx = 0, my = 0;
  const SENS = 0.003;

  /* travelZ: how far observer has moved in -Z (world). Increases with time. */
  let travelZ = 0;

  /* ── Relativistic helpers ── */
  function gamma(b) { return 1 / Math.sqrt(1 - b*b); }

  /* Aberration: given unit vector toward source in rest frame (with motion along -Z),
   * return apparent unit vector in observer frame.
   * Standard formula with cos θ measured from motion direction.
   * Motion = -Z, so "forward" component of direction = -dz.
   * cosEm  = -dz   (angle between source direction and motion axis -Z)
   * cosObs = (cosEm + beta) / (1 + beta*cosEm)    [source approaching → cosEm>0 → blueshift]
   */
  function aberrate(dx, dy, dz, b) {
    if (b === 0) return [dx, dy, dz];
    const g      = gamma(b);
    const cosEm  = -dz;                            // projection onto motion axis (-Z)
    const denom  = 1 + b * cosEm;
    const cosObs = (cosEm + b) / denom;            // observed angle from motion axis
    // transverse components scale by 1/(g*denom)
    const sc  = 1 / (g * denom);
    const ox  = dx * sc;
    const oy  = dy * sc;
    // cosObs is along motion axis (-Z), so observed z-component = -cosObs
    const oz_obs = -cosObs;
    const len = Math.sqrt(ox*ox + oy*oy + oz_obs*oz_obs);
    return [ox/len, oy/len, oz_obs/len];
  }

  /* Doppler ratio f_obs/f_em.
   * cosEm = -dz = projection of source direction onto motion axis (-Z).
   * Head-on approach (dz<0, cosEm>0): blueshift.
   */
  function doppler(dz, b) {
    if (b === 0) return 1;
    const cosEm = -dz;
    return Math.sqrt((1 + b*cosEm) / (1 - b*cosEm));
  }

  /* Map Doppler ratio to RGB hex.
   * ratio > 1 (blueshift): white → teal-dark  #23bbad = (35, 187, 173)
   * ratio < 1 (redshift):  white → pink-dark  #f92672 = (249, 38, 114)
   */
  function dopplerColor(ratio) {
    const t  = Math.log2(ratio);
    const tc = Math.max(-2.5, Math.min(2.5, t));
    if (tc >= 0) {
      const s = tc / 2.5;
      const r = Math.round(255 + s * (_TDR - 255));
      const g = Math.round(255 + s * (_TDG - 255));
      const b = Math.round(255 + s * (_TDB - 255));
      return (r << 16) | (g << 8) | b;
    } else {
      const s = -tc / 2.5;
      const r = Math.round(255 + s * (_PDR - 255));
      const g = Math.round(255 + s * (_PDG - 255));
      const b = Math.round(255 + s * (_PDB - 255));
      return (r << 16) | (g << 8) | b;
    }
  }

  /* ── Three.js init ── */
  function initThree(canvas, S) {
    renderer = new THREE.WebGLRenderer({ canvas, antialias: true });
    renderer.setSize(S, S, false);
    renderer.setClearColor(new THREE.Color(_c('--bg-void')), 1);

    scene  = new THREE.Scene();
    camera = new THREE.PerspectiveCamera(70, 1, 0.1, 800);
    camera.position.set(0, 0, 0);

    scene.add(new THREE.AmbientLight(0xffffff, 0.45));
    const sun = new THREE.DirectionalLight(0xffffff, 0.8);
    sun.position.set(0.5, 1, -1);   // from slightly forward
    scene.add(sun);

    const geo = new THREE.SphereGeometry(SPACING * 0.015, 7, 5);
    meshes = [];

    for (let ix = -NXY; ix <= NXY; ix++) {
      for (let iy = -NXY; iy <= NXY; iy++) {
        for (let iz = -NZ_F; iz <= NZ_B; iz++) {
          const mat  = new THREE.MeshPhongMaterial({ color: 0xffffff });
          const mesh = new THREE.Mesh(geo, mat);
          mesh.userData.ix = ix;
          mesh.userData.iy = iy;
          mesh.userData.iz = iz;
          scene.add(mesh);
          meshes.push(mesh);
        }
      }
    }

    travelZ = 0;
    yaw = 0; pitch = 0;
    placeAll();
  }

  /* ── Place all spheres for current state ── */
  function placeAll() {
    const b = beta;
    const g = gamma(b);

    const cosY = Math.cos(yaw), sinY = Math.sin(yaw);
    const cosP = Math.cos(pitch), sinP = Math.sin(pitch);

    for (let i = 0; i < meshes.length; i++) {
      const m = meshes[i];

      /* Rest-frame position of this lattice slot relative to camera.
       * travelZ increases as camera moves in -Z, so lattice shifts in +Z relative to camera. */
      let rx = m.userData.ix * SPACING + OFF - offsetX;
      let ry = m.userData.iy * SPACING + OFF - offsetY;
      let rz = m.userData.iz * SPACING + OFF + travelZ;

      /* Periodic tiling: wrap x,y symmetrically, wrap z into forward-biased window */
      rx -= PERIOD_XY * Math.round(rx / PERIOD_XY);
      ry -= PERIOD_XY * Math.round(ry / PERIOD_XY);
      /* Z window: [-NZ_F*SPACING, NZ_B*SPACING) */
      const Z_MIN = -NZ_F * SPACING;
      rz -= PERIOD_Z * Math.floor((rz - Z_MIN) / PERIOD_Z);

      /* Lorentz contraction along Z */
      const ox = rx;
      const oy = ry;
      const oz = rz / g;   // rz<0 ahead, contraction pulls forward spheres closer

      const dist = Math.sqrt(ox*ox + oy*oy + oz*oz);
      if (dist < 0.01) { m.visible = false; continue; }
      m.visible = true;

      /* Unit direction in world frame */
      const dx = ox/dist, dy = oy/dist, dz = oz/dist;

      /* Doppler (uses world-frame dz; dz<0 = ahead = blueshift) */
      const col = dopplerColor(doppler(dz, b));
      m.material.color.setHex(col);
      m.material.emissive.setHex(col);
      m.material.emissiveIntensity = 0.15;

      /* Aberration: shifts apparent position toward forward (-Z) at high beta */
      const [adx, ady, adz] = aberrate(dx, dy, dz, b);

      /* Place mesh along aberrated direction at contracted distance */
      m.position.set(adx*dist, ady*dist, adz*dist);
    }

    /* Camera orientation: yaw around Y, pitch around local X */
    camera.rotation.order = 'YXZ';
    camera.rotation.y = yaw;
    camera.rotation.x = pitch;
  }

  /* ── Animation loop ── */
  let lastTime = null;

  /* Visual traversal speed: travelZ grows by beta*gamma*SPACING per "lattice crossing time".
   * We pick a base rate so at low beta the lattice drifts gently, at high beta it streams. */
  const BASE_RATE = SPACING * 0.8;   // world-units per second at beta=1,gamma=1

  function loop(ts) {
    frameId = requestAnimationFrame(loop);
    if (!running) return;

    if (lastTime !== null) {
      const dt = Math.min((ts - lastTime) / 1000, 0.05);
      if (beta > 0) {
        /* In observer frame, lattice planes pass at rate beta*gamma*c/SPACING.
         * Visually we scale by BASE_RATE so it looks good on screen. */
        travelZ += beta * gamma(beta) * BASE_RATE * dt;
        /* Prevent float drift */
        travelZ  = travelZ % PERIOD_Z;
      }
    }
    lastTime = ts;

    placeAll();
    renderer.render(scene, camera);
  }

  /* ── Input handlers ── */
  function onDown(e)  { dragging = true;  mx = e.clientX; my = e.clientY; e.preventDefault(); }
  function onMove(e)  {
    if (!dragging) return;
    yaw   -= (e.clientX - mx) * SENS;
    pitch -= (e.clientY - my) * SENS;
    pitch  = Math.max(-Math.PI/2+0.01, Math.min(Math.PI/2-0.01, pitch));
    mx = e.clientX; my = e.clientY;
  }
  function onUp()     { dragging = false; }

  let tx = 0, ty = 0;
  function onTouchStart(e) { if (e.touches.length===1) { tx=e.touches[0].clientX; ty=e.touches[0].clientY; } }
  function onTouchMove(e)  {
    if (e.touches.length===1) {
      yaw   -= (e.touches[0].clientX - tx) * SENS;
      pitch -= (e.touches[0].clientY - ty) * SENS;
      pitch  = Math.max(-Math.PI/2+0.01, Math.min(Math.PI/2-0.01, pitch));
      tx=e.touches[0].clientX; ty=e.touches[0].clientY;
    }
    e.preventDefault();
  }

  let simPanel = null;

  /* ── Shell wiring ── */
  const shell = new AppletShell({
    id:    'lt',
    title: 'Relativistic Lattice',
    gap:   0,

    ctrlHTML: `
      <div class="applet-shell-ctrl-section">
        <div class="applet-shell-ctrl-title">Actions</div>
        <div class="applet-shell-btn-row">
          <button class="applet-shell-btn" onclick="ltReset()">Reset</button>
          <button class="applet-shell-btn" id="lt-pause-btn" onclick="ltTogglePause()">Pause</button>
        </div>
      </div>
      <div class="applet-shell-ctrl-section">
        <div class="applet-shell-ctrl-title">Velocity</div>
        <div class="applet-shell-slider-row">
          <span class="applet-shell-side">0</span>
          <input type="range" id="lt-vel" min="0" max="1000" step="1" value="0">
          <span class="applet-shell-side">0.999<i>c</i></span>
        </div>
        <div style="text-align:center;margin-top:4px;">
          <span class="applet-shell-val" id="lt-vel-val">&beta; = 0.000000</span>
        </div>
      </div>
      <div class="applet-shell-ctrl-section">
        <div class="applet-shell-ctrl-title">Position within cell</div>
        <div class="applet-shell-slider-row">
          <span class="applet-shell-side">L</span>
          <input type="range" id="lt-offx" min="-500" max="500" step="1" value="0">
          <span class="applet-shell-side">R</span>
        </div>
        <div class="applet-shell-slider-row">
          <span class="applet-shell-side">D</span>
          <input type="range" id="lt-offy" min="-500" max="500" step="1" value="0">
          <span class="applet-shell-side">U</span>
        </div>
      </div>
      <div class="applet-shell-ctrl-section">
        <div class="applet-shell-ctrl-title">Relativistic effects</div>
        <div style="font-size:0.78em;color:var(--text-dim);line-height:1.9;">
          <div id="lt-gamma-val">&gamma; = 1.000</div>
          <div id="lt-contract-val">Contraction: 1.000&times;</div>
          <div style="margin-top:6px;color:var(--teal-dark);">&#9632; Blueshift (ahead)</div>
          <div style="color:var(--pink-dark);">&#9632; Redshift (behind)</div>
        </div>
      </div>
    `,

    onOpen: function ({ canvas: c, S }) {
      setTimeout(() => {
        function doInit() {
          initThree(c, S);
          running  = true;
          lastTime = null;
          const pb = document.getElementById('lt-pause-btn');
          if (pb) { pb.textContent = 'Pause'; pb.classList.remove('active'); }
          if (!frameId) frameId = requestAnimationFrame(loop);

          simPanel = document.getElementById('lt-sim-panel') || c.parentElement;
          simPanel.addEventListener('mousedown',  onDown);
          window.addEventListener('mousemove',    onMove);
          window.addEventListener('mouseup',      onUp);
          simPanel.addEventListener('touchstart', onTouchStart, { passive: true });
          simPanel.addEventListener('touchmove',  onTouchMove,  { passive: false });
          simPanel.addEventListener('touchend',   onUp,         { passive: true });
        }
        if (typeof THREE !== 'undefined') { doInit(); }
        else {
          const s = document.createElement('script');
          s.src = 'https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js';
          s.onload = doInit;
          document.head.appendChild(s);
        }
      }, 80);
    },

    onClose: function () {
      running = false;
      if (frameId) { cancelAnimationFrame(frameId); frameId = null; }
      const pb = document.getElementById('lt-pause-btn');
      if (pb) { pb.textContent = 'Pause'; pb.classList.remove('active'); }
      if (simPanel) {
        simPanel.removeEventListener('mousedown',  onDown);
        simPanel.removeEventListener('touchstart', onTouchStart);
        simPanel.removeEventListener('touchmove',  onTouchMove);
        simPanel.removeEventListener('touchend',   onUp);
      }
      window.removeEventListener('mousemove', onMove);
      window.removeEventListener('mouseup',   onUp);
    },

    onResize: function ({ canvas: c, S }) {
      if (renderer) renderer.setSize(S, S, false);
      if (camera)   { camera.aspect = 1; camera.updateProjectionMatrix(); }
    },
  });

  window.ltOpen  = () => shell.open();
  window.ltClose = () => shell.close();

  window.ltReset = function () {
    travelZ = 0; yaw = 0; pitch = 0; beta = 0; offsetX = 0; offsetY = 0;
    const sl = document.getElementById('lt-vel');
    if (sl) sl.value = 0;
    const sx = document.getElementById('lt-offx');
    if (sx) sx.value = 0;
    const sy = document.getElementById('lt-offy');
    if (sy) sy.value = 0;
    updateReadout();
    if (renderer) placeAll();
  };

  window.ltTogglePause = function () {
    running = !running;
    if (running) lastTime = null;
    const pb = document.getElementById('lt-pause-btn');
    if (pb) { pb.textContent = running ? 'Pause' : 'Resume'; pb.classList.toggle('active', !running); }
  };

  /* Linear slider: x in [0,1000] → beta in [0, 0.999] */
  function sliderToBeta(x) {
    return x * 0.999 / 1000;
  }

  function updateReadout() {
    const b = beta, g = gamma(b);
    const eb = document.getElementById('lt-vel-val');
    const eg = document.getElementById('lt-gamma-val');
    const ec = document.getElementById('lt-contract-val');
    if (eb) eb.textContent = '\u03b2 = ' + b.toFixed(6);
    if (eg) eg.textContent = '\u03b3 = ' + g.toFixed(4);
    if (ec) ec.textContent = 'Contraction: ' + g.toFixed(4) + '\u00d7';
  }

  document.getElementById('lt-vel').addEventListener('input', function () {
    beta = sliderToBeta(parseFloat(this.value));
    updateReadout();
  });

  // Position sliders: range [-500,500] maps to [-SPACING/2, SPACING/2]
  document.getElementById('lt-offx').addEventListener('input', function () {
    offsetX = parseFloat(this.value) / 1000 * SPACING;
  });
  document.getElementById('lt-offy').addEventListener('input', function () {
    offsetY = parseFloat(this.value) / 1000 * SPACING;
  });

})();
