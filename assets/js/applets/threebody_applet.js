(function () {
'use strict';

const _cs  = getComputedStyle(document.documentElement);
const _c   = n => _cs.getPropertyValue(n).trim();
const _rgb  = n => { const h = _c(n).replace('#',''); const v = parseInt(h,16); return [(v>>16)&0xFF,(v>>8)&0xFF,v&0xFF]; };

// ── Simulation parameters ─────────────────────────────────────────────────────
const G = 1.0;       // gravitational constant (normalised)
const DT_BASE = 5e-4; // base timestep
let   speed = 20;    // steps per frame
let   trailLen = 400; // max trail points per body
let   r12 = 1.0;     // m1/m2
let   r23 = 1.0;     // m2/m3

// ── Bodies state ──────────────────────────────────────────────────────────────
// Each body: { m, r:[x,y,z], v:[vx,vy,vz] }
let bodies = [];
let trails = [[], [], []]; // ring-buffer of THREE.Vector3

// ── Three.js state ────────────────────────────────────────────────────────────
let renderer, scene, camera;
let simCanvas;
let bodyMeshes = [];
let trailLines = [];
let trailGeos  = [];
const BASE_RADIUS = 7.0;
let orbit = { dragging: false, lastX: 0, lastY: 0, theta: 0.4, phi: 1.1, radius: BASE_RADIUS, tx: 0, ty: 0, tz: 0 };

let running = false, frameId = null;

// ── Palette ───────────────────────────────────────────────────────────────────
// Three body colors: teal, pink, amber
const BODY_COLORS = ['--teal-light', '--pink-light', '--amber'];

// ── Presets ───────────────────────────────────────────────────────────────────
// Each preset has ratios: canonical m1/m2 and m2/m3
const PRESETS = {
  'Figure-8':    { r12: 1.00, r23: 1.00, init: function() {
    const m = 1.0;
    const x1 = 0.97000436, y1 = -0.24308753;
    const vx3 = -0.93240737, vy3 = -0.86473146;
    return [
      { m, r: [ x1,  y1, 0], v: [-vx3/2, -vy3/2, 0] },
      { m, r: [-x1, -y1, 0], v: [-vx3/2, -vy3/2, 0] },
      { m, r: [0, 0, 0],     v: [vx3, vy3, 0] },
    ];
  }},
  'Lagrange':    { r12: 1.00, r23: 1.00, init: function() {
    const m = 1.0;
    const a = 2.0, rc = a / Math.sqrt(3);
    const w = Math.sqrt(G * m / (a * a * a) * Math.sqrt(3));
    const p = [[rc,0,0],[-rc/2,rc*Math.sqrt(3)/2,0],[-rc/2,-rc*Math.sqrt(3)/2,0]];
    return p.map(([x,y,z]) => ({ m, r: [x,y,z], v: [-w*y, w*x, 0] }));
  }},
  'Broucke-Henon': { r12: 1.00, r23: 1.00, init: function() {
    const m = 1.0;
    return [
      { m, r: [-1.0, 0, 0], v: [0,  0.347113, 0] },
      { m, r: [ 1.0, 0, 0], v: [0, -0.347113, 0] },
      { m, r: [ 0.0, 0, 0], v: [0,  0, 0] },
    ];
  }},
  'Chaotic':     { r12: 0.67, r23: 0.53, init: function() {
    const bs = [
      { m: 1.0, r: [-1.2,  0.5, 0.3] },
      { m: 1.0, r: [ 1.0, -0.8, 0.0] },
      { m: 1.0, r: [ 0.3,  1.2,-0.4] },
    ];
    let M=0, cx=0, cy=0, cz=0;
    for (const b of bs) { M+=b.m; cx+=b.m*b.r[0]; cy+=b.m*b.r[1]; cz+=b.m*b.r[2]; }
    cx/=M; cy/=M; cz/=M;
    const omega = 0.6;
    return bs.map(b => {
      const dx = b.r[0]-cx, dy = b.r[1]-cy;
      return { m: b.m, r: [...b.r], v: [-omega*dy, omega*dx, 0] };
    });
  }},
};

// ── RK4 integrator ────────────────────────────────────────────────────────────
function accel(bs) {
  // Returns array of [ax, ay, az] for each body
  const n = bs.length;
  const a = bs.map(() => [0,0,0]);
  for (let i=0; i<n; i++) {
    for (let j=i+1; j<n; j++) {
      const dx = bs[j].r[0] - bs[i].r[0];
      const dy = bs[j].r[1] - bs[i].r[1];
      const dz = bs[j].r[2] - bs[i].r[2];
      const r2 = dx*dx + dy*dy + dz*dz + 1e-6; // softening
      const r3 = r2 * Math.sqrt(r2);
      const fij = G / r3;
      a[i][0] += fij * bs[j].m * dx;
      a[i][1] += fij * bs[j].m * dy;
      a[i][2] += fij * bs[j].m * dz;
      a[j][0] -= fij * bs[i].m * dx;
      a[j][1] -= fij * bs[i].m * dy;
      a[j][2] -= fij * bs[i].m * dz;
    }
  }
  return a;
}

function rk4Step(bs, dt) {
  const n = bs.length;

  // State: flat array of [r, v] per body
  function deriv(state) {
    // state: [{r,v,m}]
    const a = accel(state);
    return state.map((b,i) => ({ dr: [...b.v], dv: a[i] }));
  }

  function addScaled(state, d, h) {
    return state.map((b, i) => ({
      m: b.m,
      r: [b.r[0]+h*d[i].dr[0], b.r[1]+h*d[i].dr[1], b.r[2]+h*d[i].dr[2]],
      v: [b.v[0]+h*d[i].dv[0], b.v[1]+h*d[i].dv[1], b.v[2]+h*d[i].dv[2]],
    }));
  }

  const k1 = deriv(bs);
  const k2 = deriv(addScaled(bs, k1, dt/2));
  const k3 = deriv(addScaled(bs, k2, dt/2));
  const k4 = deriv(addScaled(bs, k3, dt));

  return bs.map((b, i) => ({
    m: b.m,
    r: [
      b.r[0] + dt/6*(k1[i].dr[0] + 2*k2[i].dr[0] + 2*k3[i].dr[0] + k4[i].dr[0]),
      b.r[1] + dt/6*(k1[i].dr[1] + 2*k2[i].dr[1] + 2*k3[i].dr[1] + k4[i].dr[1]),
      b.r[2] + dt/6*(k1[i].dr[2] + 2*k2[i].dr[2] + 2*k3[i].dr[2] + k4[i].dr[2]),
    ],
    v: [
      b.v[0] + dt/6*(k1[i].dv[0] + 2*k2[i].dv[0] + 2*k3[i].dv[0] + k4[i].dv[0]),
      b.v[1] + dt/6*(k1[i].dv[1] + 2*k2[i].dv[1] + 2*k3[i].dv[1] + k4[i].dv[1]),
      b.v[2] + dt/6*(k1[i].dv[2] + 2*k2[i].dv[2] + 2*k3[i].dv[2] + k4[i].dv[2]),
    ],
  }));
}

// ── Initialise from preset ────────────────────────────────────────────────────
function applyMassRatios(bs) {
  // Set m2 = m1/r12, m3 = m2/r23, normalise so total mass = original total
  const M0 = bs.reduce((s,b) => s+b.m, 0);
  const m1 = 1.0, m2 = m1/r12, m3 = m2/r23;
  const Mnew = m1 + m2 + m3;
  const scale = M0 / Mnew;
  bs[0].m = m1 * scale;
  bs[1].m = m2 * scale;
  bs[2].m = m3 * scale;
}

function initPreset(name, syncSliders = true) {
  const preset = PRESETS[name];
  if (syncSliders) {
    r12 = preset.r12; r23 = preset.r23;
    const s12 = document.getElementById('tb-r12'), s23 = document.getElementById('tb-r23');
    if (s12) s12.value = r12; if (s23) s23.value = r23;
  }
  bodies = preset.init();
  applyMassRatios(bodies);
  trails = bodies.map(() => []);
  // Move to center of mass frame
  let mx=0, my=0, mz=0, pvx=0, pvy=0, pvz=0, M=0;
  for (const b of bodies) {
    M += b.m;
    mx += b.m*b.r[0]; my += b.m*b.r[1]; mz += b.m*b.r[2];
    pvx += b.m*b.v[0]; pvy += b.m*b.v[1]; pvz += b.m*b.v[2];
  }
  for (const b of bodies) {
    b.r[0] -= mx/M; b.r[1] -= my/M; b.r[2] -= mz/M;
    b.v[0] -= pvx/M; b.v[1] -= pvy/M; b.v[2] -= pvz/M;
  }
}

// ── Circular point sprite texture ─────────────────────────────────────────────
function makeCircleTexture(size) {
  const c = document.createElement('canvas');
  c.width = c.height = size;
  const ctx = c.getContext('2d');
  const r = size / 2;
  const grad = ctx.createRadialGradient(r, r, 0, r, r, r);
  grad.addColorStop(0,   'rgba(255,255,255,1)');
  grad.addColorStop(0.4, 'rgba(255,255,255,1)');
  grad.addColorStop(1,   'rgba(255,255,255,0)');
  ctx.fillStyle = grad;
  ctx.fillRect(0, 0, size, size);
  return new THREE.CanvasTexture(c);
}
let spriteTexture = null;

// ── Three.js scene ────────────────────────────────────────────────────────────
function initThree() {
  simCanvas = document.getElementById('tb-canvas');
  renderer = new THREE.WebGLRenderer({ canvas: simCanvas, antialias: true });
  renderer.setPixelRatio(window.devicePixelRatio);
  renderer.setClearColor(new THREE.Color(_c('--bg-void')), 1);

  camera = new THREE.PerspectiveCamera(45, 1, 0.01, 200);
  updateCamera();
  setupOrbit();

  scene = new THREE.Scene();
  scene.add(new THREE.AmbientLight(0xffffff, 0.6));
  const dl = new THREE.DirectionalLight(0xffffff, 0.8);
  dl.position.set(5, 8, 6); scene.add(dl);

  spriteTexture = makeCircleTexture(64);

  resizeRenderer();
  buildSceneObjects();
}

function buildSceneObjects() {
  // Remove old
  bodyMeshes.forEach(m => { scene.remove(m); m.geometry.dispose(); m.material.dispose(); });
  trailLines.forEach(l => { scene.remove(l); l.geometry.dispose(); l.material.dispose(); });
  bodyMeshes = []; trailLines = []; trailGeos = [];

  bodies.forEach((b, i) => {
    const [r,g,bv] = _rgb(BODY_COLORS[i]);
    const col = new THREE.Color(r/255, g/255, bv/255);

    // Sphere — radius proportional to mass
    const rad = 0.08 * Math.cbrt(b.m);
    const geo = new THREE.SphereGeometry(rad, 16, 16);
    const mat = new THREE.MeshPhongMaterial({ color: col, emissive: col, emissiveIntensity: 0.4, shininess: 80 });
    const mesh = new THREE.Mesh(geo, mat);
    scene.add(mesh);
    bodyMeshes.push(mesh);

    // Trail as Points (thick dots, works cross-platform unlike linewidth)
    const trailGeo = new THREE.BufferGeometry();
    const positions = new Float32Array(trailLen * 3);
    trailGeo.setAttribute('position', new THREE.BufferAttribute(positions, 3));
    trailGeo.setDrawRange(0, 0);
    const trailMat = new THREE.PointsMaterial({
      color: col,
      size: 0.12,
      map: spriteTexture,
      transparent: true,
      opacity: 0.5,
      depthWrite: false,
      sizeAttenuation: true,
    });
    const line = new THREE.Points(trailGeo, trailMat);
    scene.add(line);
    trailLines.push(line);
    trailGeos.push(trailGeo);
  });
}

function updateSceneObjects() {
  bodies.forEach((b, i) => {
    // Update body position
    bodyMeshes[i].position.set(b.r[0], b.r[1], b.r[2]);

    // Append to trail
    trails[i].push(new THREE.Vector3(b.r[0], b.r[1], b.r[2]));
    if (trails[i].length > trailLen) trails[i].shift();

    // Update trail geometry
    const pts = trails[i];
    const attr = trailGeos[i].attributes.position;
    for (let k=0; k<pts.length; k++) {
      attr.array[k*3]   = pts[k].x;
      attr.array[k*3+1] = pts[k].y;
      attr.array[k*3+2] = pts[k].z;
    }
    attr.needsUpdate = true;
    trailGeos[i].setDrawRange(0, pts.length);
  });
}

function updateCamera() {
  const { theta, phi, radius, tx, ty, tz } = orbit;
  camera.position.set(
    tx + radius * Math.sin(phi) * Math.cos(theta),
    ty + radius * Math.cos(phi),
    tz + radius * Math.sin(phi) * Math.sin(theta)
  );
  camera.lookAt(tx, ty, tz);
}

function resizeRenderer() {
  const w = simCanvas.clientWidth, h = simCanvas.clientHeight;
  renderer.setSize(w, h, false);
  camera.aspect = w / h;
  camera.updateProjectionMatrix();
}

function setupOrbit() {
  simCanvas.addEventListener('pointerdown', e => {
    orbit.dragging = true; orbit.lastX = e.clientX; orbit.lastY = e.clientY;
    simCanvas.setPointerCapture(e.pointerId);
  });
  simCanvas.addEventListener('pointermove', e => {
    if (!orbit.dragging) return;
    orbit.theta -= (e.clientX - orbit.lastX) * 0.008;
    orbit.phi = Math.max(0.05, Math.min(Math.PI-0.05, orbit.phi + (e.clientY - orbit.lastY) * 0.008));
    orbit.lastX = e.clientX; orbit.lastY = e.clientY;
    updateCamera();
  });
  simCanvas.addEventListener('pointerup', () => { orbit.dragging = false; });
  simCanvas.addEventListener('wheel', e => {
    e.preventDefault();
    const oldRadius = orbit.radius;
    const newRadius = Math.max(0.5, Math.min(40, oldRadius + e.deltaY * 0.02));
    const delta = newRadius - oldRadius;
    orbit.radius = newRadius;

    // Compute mouse NDC [-1,1]
    const rect = simCanvas.getBoundingClientRect();
    const nx = ((e.clientX - rect.left) / rect.width)  * 2 - 1;
    const ny = -((e.clientY - rect.top)  / rect.height) * 2 + 1;

    // Ray direction in world space
    const fov = camera.fov * Math.PI / 180;
    const aspect = camera.aspect;
    // Camera basis vectors
    const { theta, phi, tx, ty, tz } = orbit;
    const cx = Math.sin(phi)*Math.cos(theta), cy = Math.cos(phi), cz = Math.sin(phi)*Math.sin(theta);
    // Forward = -camDir (toward target)
    const fx = -cx, fy = -cy, fz = -cz;
    // Right = forward × world-up, then normalise
    const ux = 0, uy = 1, uz = 0;
    let rx = fy*uz - fz*uy, ry = fz*ux - fx*uz, rz = fx*uy - fy*ux;
    const rl = Math.sqrt(rx*rx+ry*ry+rz*rz); rx/=rl; ry/=rl; rz/=rl;
    // Up = right × forward
    const upx = ry*fz - rz*fy, upy = rz*fx - rx*fz, upz = rx*fy - ry*fx;

    const hh = Math.tan(fov/2);
    const rdx = nx*aspect*hh*rx + ny*hh*upx + fx;
    const rdy = nx*aspect*hh*ry + ny*hh*upy + fy;
    const rdz = nx*aspect*hh*rz + ny*hh*upz + fz;
    const rl2 = Math.sqrt(rdx*rdx+rdy*rdy+rdz*rdz);

    // Shift target by fraction of zoom delta along ray
    const shift = -delta * 0.4;
    orbit.tx += shift * rdx/rl2;
    orbit.ty += shift * rdy/rl2;
    orbit.tz += shift * rdz/rl2;

    updateCamera();
  }, { passive: false });
}

// ── Animation loop ────────────────────────────────────────────────────────────
const DT = DT_BASE;
function loop() {
  if (running) {
    for (let s = 0; s < speed; s++) {
      bodies = rk4Step(bodies, DT);
    }
    updateSceneObjects();
    renderer.render(scene, camera);
  }
  frameId = requestAnimationFrame(loop);
}

// ── Shell wiring ──────────────────────────────────────────────────────────────
let currentPreset = 'Figure-8';

const shell = new AppletShell({
  id:    'tb',
  title: 'Three-Body Problem',
  gap:   0,

  ctrlHTML: `
    <div class="applet-shell-ctrl-section">
      <div class="applet-shell-ctrl-title">Actions</div>
      <div class="applet-shell-btn-row">
        <button class="applet-shell-btn" onclick="tbReset()">Reset</button>
        <button class="applet-shell-btn" id="tb-pause-btn" onclick="tbTogglePause()">Pause</button>
        <button class="applet-shell-btn" onclick="tbRandomize()">Randomize</button>
      </div>
    </div>
    <div class="applet-shell-ctrl-section">
      <div class="applet-shell-ctrl-title">Preset</div>
      <div class="applet-shell-btn-row" id="tb-preset-btns"></div>
    </div>
    <div class="applet-shell-ctrl-section">
      <div class="applet-shell-ctrl-title">Speed</div>
      <div class="applet-shell-slider-row">
        <span class="applet-shell-side">Slow</span>
        <input type="range" id="tb-speed" min="1" max="60" step="1" value="20">
        <span class="applet-shell-side">Fast</span>
      </div>
    </div>
    <div class="applet-shell-ctrl-section">
      <div class="applet-shell-ctrl-title">Mass Ratios</div>
      <div class="applet-shell-slider-row">
        <span class="applet-shell-side">m&#8321;/m&#8322;</span>
        <input type="range" id="tb-r12" min="0.01" max="1" step="0.01" value="1">
        <span class="applet-shell-side">Equal</span>
      </div>
      <div class="applet-shell-slider-row">
        <span class="applet-shell-side">m&#8322;/m&#8323;</span>
        <input type="range" id="tb-r23" min="0.01" max="1" step="0.01" value="1">
        <span class="applet-shell-side">Equal</span>
      </div>
    </div>
    <div class="applet-shell-ctrl-section">
      <div class="applet-shell-ctrl-title">Trail Length</div>
      <div class="applet-shell-slider-row">
        <span class="applet-shell-side">Short</span>
        <input type="range" id="tb-trail" min="20" max="1000" step="20" value="400">
        <span class="applet-shell-side">Long</span>
      </div>
    </div>
  `,

  onOpen: function ({ canvas: c, S }) {
    currentPreset = 'Figure-8';
    initPreset(currentPreset);

    // Build preset buttons
    const row = document.getElementById('tb-preset-btns');
    if (row) {
      row.innerHTML = '';
      Object.keys(PRESETS).forEach(name => {
        const btn = document.createElement('button');
        btn.className = 'applet-shell-btn' + (name === currentPreset ? ' active' : '');
        btn.textContent = name;
        btn.addEventListener('click', () => {
          currentPreset = name;
          initPreset(name);
          trails = bodies.map(() => []);
          buildSceneObjects();
          document.querySelectorAll('#tb-preset-btns .applet-shell-btn')
            .forEach(b => b.classList.toggle('active', b.textContent === name));
        });
        row.appendChild(btn);
      });
    }

    const pb = document.getElementById('tb-pause-btn');
    if (pb) { pb.textContent = 'Pause'; pb.classList.remove('active'); }

    function startThree() {
      setTimeout(() => {
        initThree();
        running = true;
        if (!frameId) frameId = requestAnimationFrame(loop);
      }, 80);
    }

    if (window.THREE) { startThree(); }
    else {
      const s = document.createElement('script');
      s.src = 'https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js';
      s.onload = startThree;
      document.head.appendChild(s);
    }
  },

  onClose: function () {
    running = false;
    if (frameId) { cancelAnimationFrame(frameId); frameId = null; }
  },

  onResize: function () {
    if (renderer) resizeRenderer();
  },
});

window.tbOpen  = () => shell.open();
window.tbClose = () => shell.close();

window.tbReset = function () {
  orbit.tx = 0; orbit.ty = 0; orbit.tz = 0; orbit.radius = BASE_RADIUS;
  initPreset(currentPreset, true);
  trails = bodies.map(() => []);
  if (scene) buildSceneObjects();
};

window.tbRandomize = function () {
  // Keep current positions, assign random 3D velocities
  // Scale: v ~ sqrt(G * M_total / R_rms) for rough orbital speed
  let M = 0, R2 = 0;
  for (const b of bodies) {
    M += b.m;
    R2 += b.r[0]**2 + b.r[1]**2 + b.r[2]**2;
  }
  const vscale = Math.sqrt(G * M / Math.sqrt(R2 / bodies.length));
  for (const b of bodies) {
    // Random point on unit sphere (Marsaglia method)
    let x, y, s;
    do { x = Math.random()*2-1; y = Math.random()*2-1; s = x*x+y*y; } while (s >= 1);
    const f = 2*Math.sqrt(1-s);
    const nx = x*f, ny = y*f, nz = 1-2*s;
    const mag = (0.3 + Math.random() * 0.7) * vscale;
    b.v = [nx*mag, ny*mag, nz*mag];
  }
  // Boost to CoM frame
  let pvx=0, pvy=0, pvz=0;
  for (const b of bodies) { pvx+=b.m*b.v[0]; pvy+=b.m*b.v[1]; pvz+=b.m*b.v[2]; }
  for (const b of bodies) { b.v[0]-=pvx/M; b.v[1]-=pvy/M; b.v[2]-=pvz/M; }
  trails = bodies.map(() => []);
  if (scene) buildSceneObjects();
};

window.tbTogglePause = function () {
  running = !running;
  const pb = document.getElementById('tb-pause-btn');
  if (pb) { pb.textContent = running ? 'Pause' : 'Resume'; pb.classList.toggle('active', !running); }
};

document.getElementById('tb-speed').addEventListener('input', function () {
  speed = parseInt(this.value);
});
document.getElementById('tb-trail').addEventListener('input', function () {
  trailLen = parseInt(this.value);
  trails = bodies.map(() => []);
  if (scene) buildSceneObjects();
});
document.getElementById('tb-r12').addEventListener('input', function () {
  r12 = parseFloat(this.value);
  initPreset(currentPreset, false);
  if (scene) buildSceneObjects();
});
document.getElementById('tb-r23').addEventListener('input', function () {
  r23 = parseFloat(this.value);
  initPreset(currentPreset, false);
  if (scene) buildSceneObjects();
});

})();
