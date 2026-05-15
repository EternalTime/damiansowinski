(function () {
'use strict';

/* ── Palette from CSS variables ── */
const _cs  = getComputedStyle(document.documentElement);
const _c   = n => _cs.getPropertyValue(n).trim();
const _hex = n => _c(n).replace('#', '');
function hexToRGB01(hex) {
  const v = parseInt(hex, 16);
  return [(v >> 16 & 0xFF) / 255, (v >> 8 & 0xFF) / 255, (v & 0xFF) / 255];
}
const C = {
  bgDeep:     hexToRGB01(_hex('--blue-mid')),
  tealDark:   hexToRGB01(_hex('--teal-dark')),
  tealLight:  hexToRGB01(_hex('--teal-light')),
  cyan:       hexToRGB01(_hex('--cyan')),
  pinkLight:  hexToRGB01(_hex('--pink-light')),
  greenLight: hexToRGB01(_hex('--green-light')),
  greenDark:  hexToRGB01(_hex('--green-dark')),
  blueMid:    hexToRGB01(_hex('--blue-mid')),
  blueLight:  hexToRGB01(_hex('--blue-light')),
  amber:      hexToRGB01(_hex('--amber')),
};

/* ══════════════════════════════════════════════════════════════════
   MODEL PARAMETERS — "Behavioral Speciation in the Neuroevolution
   of Recurrent Neural Net Foragers" (Sowinski 2026), adapted to 3D.
   ──────────────────────────────────────────────────────────────────
   ENVIRONMENT
   ───────────
   Resources are a Poisson point process in R³.  Each resource
   carries energy ε and decays at rate τE⁻¹.  Energy influx Γ
   maintains a mean resource energy density ρeq = Γ·τE.

   Mean resource spacing (3D generalisation of paper eq. 1):
     ℓE³ = ε / ρeq  →  ℓE = (ε / (Γ·τE))^(1/3)

   Agent length scale = diameter:
     ℓA = 2·R_COL

   Homogeneity parameter (paper eq. 1, generalised to 3D):
     Ξ = ℓA / ℓE = ℓA · (Γ·τE / ε)^(1/3)

   Inverting for ε given Ξ:
     ε = Γ · τE · ℓA³ / Ξ³  =  Γ · τE · (2·R_COL)³ / Ξ³

   Larger Ξ → more homogeneous → smaller ε → more resources.

   No-agent resource equilibrium (paper eq. 6, 3D):
     NRmax = Γ · τE · DOMAIN³ / ε  =  Ξ³ · DOMAIN³ / ℓA³

   Variable map:
     GAMMA      ↔  Γ   (energy influx density, energy·vol⁻¹·time⁻¹)
     GAMMA_DECAY↔  γ = τE⁻¹  (resource decay rate)
     TAU_E      ↔  τE = 1/γ
     epsilon    ↔  ε   (energy per resource; derived from Ξ)
     Xi         ↔  Ξ   (homogeneity; slider range 0.1–2)
     R_COL      ↔  r   (agent collection radius; ℓA = 2r)

   AGENT
   ─────
   ds/dt = -µ₀ + h  (paper eq. 3)
   Reproduces at s ≥ S_REP (daughters share fuel equally).
   Dies at s ≤ 0.

   Variable map:
     MU0     ↔  µ₀  (basal metabolic rate)
     S_MAX   ↔  S   (max fuel supply)
     S_REP       (reproduction threshold)
     vel     ↔  v   (speed)
     sigma   ↔  σ   (tumble parameter)
     R_sense     (sensor range)

   PHASE PLANE
   ───────────
   Axes are NR/NRmax and NA/NAmax, both on a log scale.
     NRmax = Γ·τE·DOMAIN³/ε   (varies with Ξ)
     NAmax = Γ·DOMAIN³/µ₀     (constant; paper eq. 7)
   At equilibrium: NR/NRmax + NA/NAmax = 1  (paper eq. 9).
══════════════════════════════════════════════════════════════════ */

/* ── Simulation constants ── */
const DOMAIN      = 100;
const GAMMA       = 1e-4;
const MU0         = 0.1;
const S_MAX       = 1.0;
const S_REP       = 0.8;
const R_COL       = 1.0;
const DT          = 0.002;
const GAMMA_DECAY = 0.1;
const TAU_E       = 1.0 / GAMMA_DECAY;
// N_A_MAX = Γ·DOMAIN³/µ₀ — constant, independent of Ξ (eq. 7)
const N_A_MAX = Math.round(GAMMA * DOMAIN ** 3 / MU0);

let R_sense = 6.0, vel = 20.0, sigma = 0.1, Xi = 0.5;
let epsilon, stepsPerFrame = 5;

const L_A = 2 * R_COL;  // agent length scale = diameter

// ε = Γ·τE·ℓA³/Ξ³  (3D generalisation of paper eq. 1)
function epsilonFromXi(xi) { return GAMMA * TAU_E * L_A ** 3 / xi ** 3; }
// NRmax = Γ·τE·DOMAIN³/ε  (paper eq. 6, no-agent resource equilibrium)
function nRMax() { return GAMMA * TAU_E * DOMAIN ** 3 / epsilonFromXi(Xi); }
// NReq = NRmax at current Xi (used to seed initial resources)
function nEqResources() { return nRMax(); }

/* ── State ── */
let resources = [], agents = [], simTime = 0;
let phaseTrail = [];
let running = false, frameId = null, initialized = false, extinct = false;

/* ── Spatial grid (3D) ── */
let grid = null, gridNC = 0, gridCS = 0;

function buildGrid() {
  gridNC = Math.max(1, Math.ceil(DOMAIN / R_sense));
  gridCS = DOMAIN / gridNC;
  const total = gridNC * gridNC * gridNC;
  grid = new Array(total);
  for (let i = 0; i < total; i++) grid[i] = [];
  for (let i = 0; i < resources.length; i++) {
    const r = resources[i];
    const ix = Math.min(gridNC - 1, Math.floor(r.x / gridCS));
    const iy = Math.min(gridNC - 1, Math.floor(r.y / gridCS));
    const iz = Math.min(gridNC - 1, Math.floor(r.z / gridCS));
    grid[iz * gridNC * gridNC + iy * gridNC + ix].push(i);
  }
}

function gridNeighbors3(ax, ay, az) {
  const cx = Math.min(gridNC - 1, Math.floor(ax / gridCS));
  const cy = Math.min(gridNC - 1, Math.floor(ay / gridCS));
  const cz = Math.min(gridNC - 1, Math.floor(az / gridCS));
  const out = [];
  for (let dz = -1; dz <= 1; dz++) {
    const iz = ((cz + dz) % gridNC + gridNC) % gridNC;
    for (let dy = -1; dy <= 1; dy++) {
      const iy = ((cy + dy) % gridNC + gridNC) % gridNC;
      for (let dx = -1; dx <= 1; dx++) {
        const ix = ((cx + dx) % gridNC + gridNC) % gridNC;
        const cell = grid[iz * gridNC * gridNC + iy * gridNC + ix];
        for (let k = 0; k < cell.length; k++) out.push(cell[k]);
      }
    }
  }
  return out;
}

/* ── Math helpers ── */
function randn() {
  return Math.sqrt(-2 * Math.log(Math.random() || 1e-20)) * Math.cos(2 * Math.PI * Math.random());
}
function wrap(v) { return ((v % DOMAIN) + DOMAIN) % DOMAIN; }
function toroidalDisp3(ax, ay, az, bx, by, bz) {
  let dx = bx - ax, dy = by - ay, dz = bz - az;
  const h = DOMAIN / 2;
  if (dx >  h) dx -= DOMAIN; if (dx < -h) dx += DOMAIN;
  if (dy >  h) dy -= DOMAIN; if (dy < -h) dy += DOMAIN;
  if (dz >  h) dz -= DOMAIN; if (dz < -h) dz += DOMAIN;
  return [dx, dy, dz];
}
function poissonSample(lam) {
  if (lam <= 0) return 0;
  if (lam > 30) return Math.max(0, Math.round(lam + Math.sqrt(lam) * randn()));
  const L = Math.exp(-lam); let k = 0, p = 1;
  do { k++; p *= Math.random(); } while (p > L);
  return k - 1;
}

/* ── Simulation ── */
function initResources() {
  resources = [];
  const n = Math.round(nEqResources());
  for (let i = 0; i < n; i++)
    resources.push({ x: Math.random()*DOMAIN, y: Math.random()*DOMAIN, z: Math.random()*DOMAIN });
}

function makeAgent(x, y, z, theta, phi, s) {
  const t = theta !== undefined ? theta : Math.random() * 2 * Math.PI;
  const p = phi   !== undefined ? phi   : Math.acos(2 * Math.random() - 1);
  return {
    x: x !== undefined ? x : Math.random() * DOMAIN,
    y: y !== undefined ? y : Math.random() * DOMAIN,
    z: z !== undefined ? z : Math.random() * DOMAIN,
    theta: t, phi: p,
    s: s !== undefined ? s : S_MAX * (0.3 + Math.random() * 0.5),
  };
}

function agentStep(a) {
  const cands = gridNeighbors3(a.x, a.y, a.z);
  let best2 = R_sense * R_sense, bestIdx = -1;
  for (let k = 0; k < cands.length; k++) {
    const i = cands[k];
    if (i >= resources.length) continue;
    const [dx, dy, dz] = toroidalDisp3(a.x, a.y, a.z, resources[i].x, resources[i].y, resources[i].z);
    const d2 = dx*dx + dy*dy + dz*dz;
    if (d2 < best2) { best2 = d2; bestIdx = i; }
  }
  if (bestIdx >= 0) {
    const [dx, dy, dz] = toroidalDisp3(a.x, a.y, a.z, resources[bestIdx].x, resources[bestIdx].y, resources[bestIdx].z);
    a.theta = Math.atan2(dy, dx);
    a.phi   = Math.acos(Math.max(-1, Math.min(1, dz / Math.sqrt(dx*dx + dy*dy + dz*dz))));
  } else {
    // random walk on sphere via small angular perturbations
    a.theta += sigma * randn() * Math.sqrt(DT);
    a.phi   = Math.max(0.01, Math.min(Math.PI - 0.01, a.phi + sigma * randn() * Math.sqrt(DT)));
  }
  const sp = Math.sin(a.phi);
  a.x = wrap(a.x + vel * sp * Math.cos(a.theta) * DT);
  a.y = wrap(a.y + vel * sp * Math.sin(a.theta) * DT);
  a.z = wrap(a.z + vel * Math.cos(a.phi) * DT);
  let harvested = 0;
  const colCands = gridNeighbors3(a.x, a.y, a.z);
  for (let k = colCands.length - 1; k >= 0; k--) {
    const i = colCands[k];
    if (i >= resources.length) continue;
    const [dx, dy, dz] = toroidalDisp3(a.x, a.y, a.z, resources[i].x, resources[i].y, resources[i].z);
    if (dx*dx + dy*dy + dz*dz <= R_COL * R_COL) {
      resources[i] = resources[resources.length - 1]; resources.pop(); harvested++;
    }
  }
  a.s = Math.min(S_MAX, a.s - MU0 * DT + harvested * epsilon);
}

function updateResources() {
  const birthMean = GAMMA * DOMAIN ** 3 * DT / epsilon;
  const nBirth = poissonSample(birthMean);
  for (let i = 0; i < nBirth; i++)
    resources.push({ x: Math.random()*DOMAIN, y: Math.random()*DOMAIN, z: Math.random()*DOMAIN });
  const pDecay = 1 - Math.exp(-GAMMA_DECAY * DT);
  for (let i = resources.length - 1; i >= 0; i--)
    if (Math.random() < pDecay) { resources[i] = resources[resources.length - 1]; resources.pop(); }
}

function simStep() {
  updateResources();
  buildGrid();
  const newborns = [];
  for (let i = agents.length - 1; i >= 0; i--) {
    const a = agents[i];
    agentStep(a);
    if (a.s >= S_REP) {
      const t = Math.random() * 2 * Math.PI;
      const p = Math.acos(2 * Math.random() - 1);
      a.s /= 2;
      const sp = Math.sin(p);
      newborns.push(makeAgent(
        wrap(a.x + Math.cos(t)*sp*R_COL),
        wrap(a.y + Math.sin(t)*sp*R_COL),
        wrap(a.z + Math.cos(p)*R_COL),
        Math.random()*2*Math.PI, Math.acos(2*Math.random()-1), a.s
      ));
    }
    if (a.s <= 0) agents.splice(i, 1);
  }
  agents.push(...newborns);
  simTime++;

  const nrNorm = resources.length / Math.max(1, nRMax());
  const naNorm = agents.length / Math.max(1, N_A_MAX);
  phaseTrail.push({ nr: nrNorm, na: naNorm });
  if (phaseTrail.length > 200) phaseTrail.shift();

  if (agents.length === 0 && !extinct) {
    extinct = true; running = false;
    const btn = document.getElementById('fa3-btn-run');
    if (btn) { btn.textContent = 'Run'; btn.classList.remove('active'); }
    const stamp = document.getElementById('fa3-extinction-stamp');
    if (stamp) stamp.style.display = 'block';
  }
}

function fa3InitSim() {
  simTime = 0; phaseTrail = []; extinct = false;
  epsilon = epsilonFromXi(Xi);
  const stamp = document.getElementById('fa3-extinction-stamp');
  if (stamp) stamp.style.display = 'none';
  initResources();
  const nA = Math.max(1, Math.round(N_A_MAX * 0.1)); // start at 10% of N_A_MAX
  agents = Array.from({ length: nA }, () => makeAgent());
  phaseTrail.push({ nr: 1, na: 1 });
}

/* ══════════════════════════════════════════════════════════════════
   WebGL Renderer
══════════════════════════════════════════════════════════════════ */
let gl = null, simCanvas3 = null;
let progSphere = null, progArrow = null, progCube = null, progSense = null, progAgent = null;
let sphereUnitBuf = null, sphereUnitCount = 0;  // unit icosphere triangles
let senseCenterBuf = null;                       // per-agent centers+radii
let agentInstBuf = null;                         // per-agent: x,y,z,radius (4 floats)

/* Orbit state */
let orbit = { theta: 0.6, phi: 1.1, dist: 90, dragging: false, lastX: 0, lastY: 0 };

/* ── Shader sources ── */
const VS_SPHERE = `
attribute vec3 aPos;
attribute vec3 aOffset;      // instance world position
attribute float aSize;       // point radius in world units
attribute vec3 aColor;       // outer color
attribute vec3 aColorInner;  // inner color
uniform mat4 uMVP;
uniform vec2 uViewport;
varying vec3 vColor;
varying vec3 vColorInner;
varying float vEyeDist;
void main() {
  vec4 clip = uMVP * vec4(aOffset, 1.0);
  gl_Position = clip;
  float eyeDist = clip.w;
  vEyeDist = eyeDist;
  gl_PointSize = clamp(aSize * uViewport.y / eyeDist, 2.0, 48.0);
  vColor = aColor;
  vColorInner = aColorInner;
}`;

const FS_SPHERE = `
precision mediump float;
varying vec3 vColor;
varying vec3 vColorInner;
varying float vEyeDist;
uniform float uNear;
uniform float uFar;
void main() {
  vec2 uv = gl_PointCoord * 2.0 - 1.0;
  float r2 = dot(uv, uv);
  if (r2 > 1.0) discard;
  float r = sqrt(r2);
  vec3 col = mix(vColorInner, vColor, r * r);
  float alpha = smoothstep(1.0, 0.5, r);
  float fog = 1.0 - clamp((vEyeDist - uNear) / (uFar - uNear), 0.0, 1.0);
  gl_FragColor = vec4(col, alpha * fog);
}`;

const VS_LINE = `
attribute vec3 aPos;
uniform mat4 uMVP;
uniform vec3 uColor;
varying vec3 vColor;
void main() {
  gl_Position = uMVP * vec4(aPos, 1.0);
  vColor = uColor;
}`;

const FS_LINE = `
precision mediump float;
varying vec3 vColor;
void main() { gl_FragColor = vec4(vColor, 0.6); }`;

/* Sense-sphere mesh shader — translucent wireframe sphere per agent */
const VS_SENSE = `
attribute vec3 aPos;       // unit sphere vertex
attribute vec3 aCenter;    // agent world position (instance)
attribute float aRadius;   // R_sense
uniform mat4 uMVP;
void main() {
  gl_Position = uMVP * vec4(aCenter + aPos * aRadius, 1.0);
}`;

const FS_SENSE = `
precision mediump float;
uniform vec3 uColor;
void main() { gl_FragColor = vec4(uColor, 0.01); }`;

/* Agent mesh shader — instanced icosphere with rim lighting */
const VS_AGENT = `
attribute vec3 aPos;       // unit sphere normal/position
attribute vec3 aCenter;    // agent world position
attribute float aRadius;   // world-space radius
uniform mat4 uMVP;
uniform vec3 uEye;
varying float vRim;
varying float vEyeDist;
void main() {
  vec3 worldPos = aCenter + aPos * aRadius;
  vec4 clip = uMVP * vec4(worldPos, 1.0);
  gl_Position = clip;
  vEyeDist = clip.w;
  vec3 toEye = normalize(uEye - aCenter);
  vRim = 1.0 - abs(dot(aPos, toEye));
}`;

const FS_AGENT = `
precision mediump float;
uniform vec3 uColorFace;
uniform vec3 uColorRim;
uniform float uNear;
uniform float uFar;
varying float vRim;
varying float vEyeDist;
void main() {
  vec3 col = mix(uColorFace, uColorRim, pow(vRim, 2.0));
  float fog = 1.0 - clamp((vEyeDist - uNear) / (uFar - uNear), 0.0, 1.0);
  gl_FragColor = vec4(col, fog);
}`;

/* ── Icosphere builder (1 subdivision) ── */
function buildUnitSphere() {
  const t = (1 + Math.sqrt(5)) / 2;
  const vraw = [
    [-1,t,0],[1,t,0],[-1,-t,0],[1,-t,0],
    [0,-1,t],[0,1,t],[0,-1,-t],[0,1,-t],
    [t,0,-1],[t,0,1],[-t,0,-1],[-t,0,1]
  ].map(v => { const l=Math.sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]); return [v[0]/l,v[1]/l,v[2]/l]; });
  const faces = [
    [0,11,5],[0,5,1],[0,1,7],[0,7,10],[0,10,11],
    [1,5,9],[5,11,4],[11,10,2],[10,7,6],[7,1,8],
    [3,9,4],[3,4,2],[3,2,6],[3,6,8],[3,8,9],
    [4,9,5],[2,4,11],[6,2,10],[8,6,7],[9,8,1]
  ];
  function mid(a,b) {
    const v=[a[0]+b[0],a[1]+b[1],a[2]+b[2]];
    const l=Math.sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
    return [v[0]/l,v[1]/l,v[2]/l];
  }
  const verts = [];
  for (const [a,b,c] of faces) {
    const va=vraw[a], vb=vraw[b], vc=vraw[c];
    const ab=mid(va,vb), bc=mid(vb,vc), ca=mid(vc,va);
    for (const tri of [[va,ab,ca],[vb,bc,ab],[vc,ca,bc],[ab,bc,ca]])
      for (const p of tri) verts.push(...p);
  }
  return new Float32Array(verts);
}

/* ── GL helpers ── */
function compileShader(type, src) {
  const s = gl.createShader(type);
  gl.shaderSource(s, src); gl.compileShader(s);
  if (!gl.getShaderParameter(s, gl.COMPILE_STATUS))
    console.error('Shader error:', gl.getShaderInfoLog(s));
  return s;
}
function makeProgram(vs, fs) {
  const p = gl.createProgram();
  gl.attachShader(p, compileShader(gl.VERTEX_SHADER,   vs));
  gl.attachShader(p, compileShader(gl.FRAGMENT_SHADER, fs));
  gl.linkProgram(p);
  return p;
}

/* ── Matrix math (column-major) ── */
function mat4Identity() { return new Float32Array([1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1]); }
function mat4Mul(a, b) {
  const o = new Float32Array(16);
  for (let i = 0; i < 4; i++)
    for (let j = 0; j < 4; j++)
      for (let k = 0; k < 4; k++)
        o[j*4+i] += a[k*4+i] * b[j*4+k];
  return o;
}
function mat4Perspective(fovy, aspect, near, far) {
  const f = 1 / Math.tan(fovy / 2);
  const nf = 1 / (near - far);
  return new Float32Array([
    f/aspect, 0, 0, 0,
    0, f, 0, 0,
    0, 0, (far+near)*nf, -1,
    0, 0, 2*far*near*nf, 0
  ]);
}
function mat4LookAt(eye, center, up) {
  let fx = center[0]-eye[0], fy = center[1]-eye[1], fz = center[2]-eye[2];
  const fl = Math.sqrt(fx*fx+fy*fy+fz*fz);
  fx/=fl; fy/=fl; fz/=fl;
  let rx = fy*up[2]-fz*up[1], ry = fz*up[0]-fx*up[2], rz = fx*up[1]-fy*up[0];
  const rl = Math.sqrt(rx*rx+ry*ry+rz*rz);
  rx/=rl; ry/=rl; rz/=rl;
  const ux = ry*fz-rz*fy, uy = rz*fx-rx*fz, uz = rx*fy-ry*fx;
  return new Float32Array([
    rx, ux, -fx, 0,
    ry, uy, -fy, 0,
    rz, uz, -fz, 0,
    -(rx*eye[0]+ry*eye[1]+rz*eye[2]),
    -(ux*eye[0]+uy*eye[1]+uz*eye[2]),
    fx*eye[0]+fy*eye[1]+fz*eye[2], 1
  ]);
}

/* ── Buffers ── */
let sphereBuf = null;   // interleaved: x,y,z,size,r,g,b  per point
let arrowBuf  = null;   // line segments: x0,y0,z0, x1,y1,z1 per arrow
let cubeBuf   = null;   // 12 edges × 2 verts × 3 floats

function initGL(canvas) {
  gl = canvas.getContext('webgl', { antialias: true, alpha: false });
  if (!gl) { console.error('WebGL not supported'); return; }
  progSphere = makeProgram(VS_SPHERE, FS_SPHERE);
  progArrow  = makeProgram(VS_LINE,   FS_LINE);
  progCube   = makeProgram(VS_LINE,   FS_LINE);
  progSense  = makeProgram(VS_SENSE,  FS_SENSE);
  progAgent  = makeProgram(VS_AGENT,  FS_AGENT);
  sphereBuf     = gl.createBuffer();
  arrowBuf      = gl.createBuffer();
  cubeBuf       = gl.createBuffer();
  sphereUnitBuf = gl.createBuffer();
  senseCenterBuf= gl.createBuffer();
  agentInstBuf  = gl.createBuffer();
  buildCube();
  const unitVerts = buildUnitSphere();
  sphereUnitCount = unitVerts.length / 3;
  gl.bindBuffer(gl.ARRAY_BUFFER, sphereUnitBuf);
  gl.bufferData(gl.ARRAY_BUFFER, unitVerts, gl.STATIC_DRAW);
  gl.enable(gl.BLEND);
  gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA);
  gl.enable(gl.DEPTH_TEST);
  gl.depthFunc(gl.LEQUAL);
  gl.clearColor(...C.bgDeep, 1.0);
}

function buildCube() {
  const d = DOMAIN;
  const verts = [
    0,0,0, d,0,0,  d,0,0, d,d,0,  d,d,0, 0,d,0,  0,d,0, 0,0,0,
    0,0,d, d,0,d,  d,0,d, d,d,d,  d,d,d, 0,d,d,  0,d,d, 0,0,d,
    0,0,0, 0,0,d,  d,0,0, d,0,d,  d,d,0, d,d,d,  0,d,0, 0,d,d,
  ];
  gl.bindBuffer(gl.ARRAY_BUFFER, cubeBuf);
  gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(verts), gl.STATIC_DRAW);
}

function getMVP() {
  const W = simCanvas3.width, H = simCanvas3.height;
  const cx = DOMAIN/2, cy = DOMAIN/2, cz = DOMAIN/2;
  const eye = [
    cx + orbit.dist * Math.sin(orbit.phi) * Math.cos(orbit.theta),
    cy + orbit.dist * Math.sin(orbit.phi) * Math.sin(orbit.theta),
    cz + orbit.dist * Math.cos(orbit.phi),
  ];
  const view = mat4LookAt(eye, [cx, cy, cz], [0, 0, 1]);
  const proj = mat4Perspective(Math.PI / 4, W / H, 1, 2000);
  return { mvp: mat4Mul(proj, view), eye };
}

function drawScene() {
  if (!gl) return;
  const W = simCanvas3.width, H = simCanvas3.height;
  gl.viewport(0, 0, W, H);
  gl.depthMask(true);
  gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);

  const { mvp, eye } = getMVP();

  /* ── Cube wireframe ── */
  gl.useProgram(progCube);
  gl.uniformMatrix4fv(gl.getUniformLocation(progCube, 'uMVP'), false, mvp);
  gl.uniform3f(gl.getUniformLocation(progCube, 'uColor'), ...C.blueMid);
  gl.bindBuffer(gl.ARRAY_BUFFER, cubeBuf);
  const aPosCube = gl.getAttribLocation(progCube, 'aPos');
  gl.enableVertexAttribArray(aPosCube);
  gl.vertexAttribPointer(aPosCube, 3, gl.FLOAT, false, 0, 0);
  gl.drawArrays(gl.LINES, 0, 24);

  /* ── Point sprites: 10 floats per point: x,y,z, size, r,g,b (outer), r,g,b (inner) ── */
  const stride = 10 * 4;
  gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA);
  gl.useProgram(progSphere);
  gl.uniformMatrix4fv(gl.getUniformLocation(progSphere, 'uMVP'), false, mvp);
  gl.uniform2f(gl.getUniformLocation(progSphere, 'uViewport'), W, H);
  gl.uniform1f(gl.getUniformLocation(progSphere, 'uNear'), orbit.dist * 0.3);
  gl.uniform1f(gl.getUniformLocation(progSphere, 'uFar'),  orbit.dist * 1.5);
  const aPosS      = gl.getAttribLocation(progSphere, 'aOffset');
  const aSizeS     = gl.getAttribLocation(progSphere, 'aSize');
  const aColS      = gl.getAttribLocation(progSphere, 'aColor');
  const aColInnerS = gl.getAttribLocation(progSphere, 'aColorInner');

  /* Resources */
  const rData = new Float32Array(resources.length * 10);
  for (let i = 0; i < resources.length; i++) {
    const r = resources[i];
    rData[i*10+0] = r.x; rData[i*10+1] = r.y; rData[i*10+2] = r.z;
    rData[i*10+3] = 1.25 * R_COL / Xi;  // sprite radius ∝ ε^(1/3) ∝ 1/Ξ
    rData[i*10+4] = C.greenLight[0]; rData[i*10+5] = C.greenLight[1]; rData[i*10+6] = C.greenLight[2];
    rData[i*10+7] = C.greenDark[0];  rData[i*10+8] = C.greenDark[1];  rData[i*10+9] = C.greenDark[2];
  }
  gl.bindBuffer(gl.ARRAY_BUFFER, sphereBuf);
  gl.bufferData(gl.ARRAY_BUFFER, rData, gl.DYNAMIC_DRAW);
  gl.enableVertexAttribArray(aPosS);      gl.vertexAttribPointer(aPosS,      3, gl.FLOAT, false, stride, 0);
  gl.enableVertexAttribArray(aSizeS);     gl.vertexAttribPointer(aSizeS,     1, gl.FLOAT, false, stride, 12);
  gl.enableVertexAttribArray(aColS);      gl.vertexAttribPointer(aColS,      3, gl.FLOAT, false, stride, 16);
  gl.enableVertexAttribArray(aColInnerS); gl.vertexAttribPointer(aColInnerS, 3, gl.FLOAT, false, stride, 28);
  if (resources.length > 0) gl.drawArrays(gl.POINTS, 0, resources.length);

  /* ── Agents: instanced icosphere meshes + arrow lines ── */
  if (agents.length > 0) {
    const extA = gl.getExtension('ANGLE_instanced_arrays');
    const AGENT_BASE_R = R_COL * 1.2;
    const agInstData = new Float32Array(agents.length * 4);
    const arrowVerts  = new Float32Array(agents.length * 6);
    const ARROW_LEN   = AGENT_BASE_R * 2.5;
    for (let i = 0; i < agents.length; i++) {
      const a = agents[i];
      const hf = Math.max(0, a.s / S_MAX);
      agInstData[i*4+0] = a.x;
      agInstData[i*4+1] = a.y;
      agInstData[i*4+2] = a.z;
      agInstData[i*4+3] = AGENT_BASE_R * (0.35 + 0.65 * hf);
      const sp = Math.sin(a.phi);
      arrowVerts[i*6+0] = a.x; arrowVerts[i*6+1] = a.y; arrowVerts[i*6+2] = a.z;
      arrowVerts[i*6+3] = a.x + Math.cos(a.theta)*sp*ARROW_LEN;
      arrowVerts[i*6+4] = a.y + Math.sin(a.theta)*sp*ARROW_LEN;
      arrowVerts[i*6+5] = a.z + Math.cos(a.phi)*ARROW_LEN;
    }

    gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA);
    gl.useProgram(progAgent);
    gl.uniformMatrix4fv(gl.getUniformLocation(progAgent, 'uMVP'), false, mvp);
    gl.uniform3fv(gl.getUniformLocation(progAgent, 'uEye'), eye);
    gl.uniform3f(gl.getUniformLocation(progAgent, 'uColorFace'), ...C.pinkLight);
    gl.uniform3f(gl.getUniformLocation(progAgent, 'uColorRim'),  ...C.blueMid);
    gl.uniform1f(gl.getUniformLocation(progAgent, 'uNear'), orbit.dist * 0.3);
    gl.uniform1f(gl.getUniformLocation(progAgent, 'uFar'),  orbit.dist * 1.5);

    const aPosAg  = gl.getAttribLocation(progAgent, 'aPos');
    const aCenAg  = gl.getAttribLocation(progAgent, 'aCenter');
    const aRadAg  = gl.getAttribLocation(progAgent, 'aRadius');

    gl.bindBuffer(gl.ARRAY_BUFFER, sphereUnitBuf);
    gl.enableVertexAttribArray(aPosAg);
    gl.vertexAttribPointer(aPosAg, 3, gl.FLOAT, false, 0, 0);
    if (extA) extA.vertexAttribDivisorANGLE(aPosAg, 0);

    gl.bindBuffer(gl.ARRAY_BUFFER, agentInstBuf);
    gl.bufferData(gl.ARRAY_BUFFER, agInstData, gl.DYNAMIC_DRAW);
    gl.enableVertexAttribArray(aCenAg);
    gl.vertexAttribPointer(aCenAg, 3, gl.FLOAT, false, 16, 0);
    gl.enableVertexAttribArray(aRadAg);
    gl.vertexAttribPointer(aRadAg, 1, gl.FLOAT, false, 16, 12);

    if (extA) {
      extA.vertexAttribDivisorANGLE(aCenAg, 1);
      extA.vertexAttribDivisorANGLE(aRadAg, 1);
      extA.drawArraysInstancedANGLE(gl.TRIANGLES, 0, sphereUnitCount, agents.length);
      extA.vertexAttribDivisorANGLE(aCenAg, 0);
      extA.vertexAttribDivisorANGLE(aRadAg, 0);
    }

    /* Arrow lines */
    gl.blendFunc(gl.SRC_ALPHA, gl.ONE);
    gl.useProgram(progArrow);
    gl.uniformMatrix4fv(gl.getUniformLocation(progArrow, 'uMVP'), false, mvp);
    gl.uniform3f(gl.getUniformLocation(progArrow, 'uColor'), ...C.tealLight);
    gl.bindBuffer(gl.ARRAY_BUFFER, arrowBuf);
    gl.bufferData(gl.ARRAY_BUFFER, arrowVerts, gl.DYNAMIC_DRAW);
    const aPosA = gl.getAttribLocation(progArrow, 'aPos');
    gl.enableVertexAttribArray(aPosA);
    gl.vertexAttribPointer(aPosA, 3, gl.FLOAT, false, 0, 0);
    gl.drawArrays(gl.LINES, 0, agents.length * 2);
  }

  /* Sense spheres — standard alpha blend, no depth writes */
  gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA);
  gl.depthMask(false);
  /* ── Sense spheres (translucent mesh, one per agent) ── */
  if (agents.length > 0) {
    // Pack center+radius per agent: 4 floats each
    const senseData = new Float32Array(agents.length * 4);
    for (let i = 0; i < agents.length; i++) {
      senseData[i*4+0] = agents[i].x;
      senseData[i*4+1] = agents[i].y;
      senseData[i*4+2] = agents[i].z;
      senseData[i*4+3] = R_sense;
    }
    gl.useProgram(progSense);
    gl.uniformMatrix4fv(gl.getUniformLocation(progSense, 'uMVP'), false, mvp);
    gl.uniform3f(gl.getUniformLocation(progSense, 'uColor'), ...C.blueLight);

    const aPosU = gl.getAttribLocation(progSense, 'aPos');
    const aCen  = gl.getAttribLocation(progSense, 'aCenter');
    const aRad  = gl.getAttribLocation(progSense, 'aRadius');

    const ext = gl.getExtension('ANGLE_instanced_arrays');
    gl.bindBuffer(gl.ARRAY_BUFFER, sphereUnitBuf);
    gl.enableVertexAttribArray(aPosU);
    gl.vertexAttribPointer(aPosU, 3, gl.FLOAT, false, 0, 0);
    if (ext) ext.vertexAttribDivisorANGLE(aPosU, 0);

    gl.bindBuffer(gl.ARRAY_BUFFER, senseCenterBuf);
    gl.bufferData(gl.ARRAY_BUFFER, senseData, gl.DYNAMIC_DRAW);
    gl.enableVertexAttribArray(aCen);
    gl.vertexAttribPointer(aCen, 3, gl.FLOAT, false, 16, 0);
    gl.enableVertexAttribArray(aRad);
    gl.vertexAttribPointer(aRad, 1, gl.FLOAT, false, 16, 12);

    if (ext) {
      ext.vertexAttribDivisorANGLE(aCen, 1);
      ext.vertexAttribDivisorANGLE(aRad, 1);
      ext.drawArraysInstancedANGLE(gl.TRIANGLES, 0, sphereUnitCount, agents.length);
      ext.vertexAttribDivisorANGLE(aCen, 0);
      ext.vertexAttribDivisorANGLE(aRad, 0);
    }
  }
  gl.depthMask(true);
}

/* ══════════════════════════════════════════════════════════════════
   Phase plane (2D canvas)
══════════════════════════════════════════════════════════════════ */
const _rgba = (n, a) => {
  const h = _c(n).replace('#', '');
  const v = parseInt(h, 16);
  return `rgba(${(v>>16)&0xFF},${(v>>8)&0xFF},${v&0xFF},${a})`;
};
const [_GLR, _GLG, _GLB] = [
  Math.round(C.greenLight[0]*255),
  Math.round(C.greenLight[1]*255),
  Math.round(C.greenLight[2]*255),
];

let phaseCanvas3, phaseCtx3, accumCanvas3, accumCtx3;
const BLOB_A = 0.12, PHASE_LEN = 200;
const MARGIN = 0.15, MARGIN_INNER = 0.075;
const LOG_MIN = 0.01, LOG_MAX = 10.0;
const LOG_SPAN = Math.log10(LOG_MAX / LOG_MIN);

function logScale(v) {
  return (Math.log10(Math.max(LOG_MIN, Math.min(LOG_MAX, v))) - Math.log10(LOG_MIN)) / LOG_SPAN;
}
function phaseToPixel(nr, na, W, H) {
  const x0=W*MARGIN, x1=W*(1-MARGIN_INNER), y0=H*MARGIN_INNER, y1=H*(1-MARGIN);
  return [x0 + logScale(nr)*(x1-x0), y1 - logScale(na)*(y1-y0)];
}

function initAccum3() {
  accumCanvas3 = document.createElement('canvas');
  accumCanvas3.width  = phaseCanvas3.width;
  accumCanvas3.height = phaseCanvas3.height;
  accumCtx3 = accumCanvas3.getContext('2d');
}
function resizeAccum3() {
  const tmp = document.createElement('canvas');
  tmp.width = phaseCanvas3.width; tmp.height = phaseCanvas3.height;
  const tCtx = tmp.getContext('2d');
  if (accumCanvas3) tCtx.drawImage(accumCanvas3, 0, 0, tmp.width, tmp.height);
  accumCanvas3 = tmp; accumCtx3 = tCtx;
}
function clearAccum3() {
  if (accumCtx3) accumCtx3.clearRect(0, 0, accumCanvas3.width, accumCanvas3.height);
}
function addBlob3(nr, na, alphaMul) {
  const W=accumCanvas3.width, H=accumCanvas3.height;
  const [cx, cy] = phaseToPixel(nr, na, W, H);
  const r = W * 0.03, a = BLOB_A * alphaMul;
  const grad = accumCtx3.createRadialGradient(cx, cy, 0, cx, cy, r);
  grad.addColorStop(0,   `rgba(${_GLR},${_GLG},${_GLB},${a})`);
  grad.addColorStop(0.5, `rgba(${_GLR},${_GLG},${_GLB},${a * 0.4})`);
  grad.addColorStop(1,   `rgba(${_GLR},${_GLG},${_GLB},0)`);
  accumCtx3.save();
  accumCtx3.globalCompositeOperation = 'lighter';
  accumCtx3.fillStyle = grad;
  accumCtx3.beginPath(); accumCtx3.arc(cx, cy, r, 0, 2*Math.PI); accumCtx3.fill();
  accumCtx3.restore();
}

function drawPhase3() {
  const W=phaseCanvas3.width, H=phaseCanvas3.height;
  if (W===0||H===0) return;
  const x0=W*MARGIN, x1=W*(1-MARGIN_INNER), y0=H*MARGIN_INNER, y1=H*(1-MARGIN);
  clearAccum3();
  const decay = Math.log(200) / 100;
  for (let i = 0; i < phaseTrail.length; i++)
    addBlob3(phaseTrail[i].nr, phaseTrail[i].na, Math.exp(-decay*(phaseTrail.length-1-i)));
  phaseCtx3.fillStyle = _c('--bg-dark'); phaseCtx3.fillRect(0, 0, W, H);
  phaseCtx3.drawImage(accumCanvas3, 0, 0);
  phaseCtx3.save();
  phaseCtx3.setLineDash([5, 4]);
  phaseCtx3.strokeStyle = _rgba('--pink-light', 0.7);
  phaseCtx3.lineWidth = 1;
  phaseCtx3.beginPath();
  let started = false;
  for (let k = 0; k <= 200; k++) {
    const nr = LOG_MIN + (1 - 2*LOG_MIN)*k/200, na = 1 - nr;
    if (na < LOG_MIN || nr < LOG_MIN) continue;
    const [px, py] = phaseToPixel(nr, na, W, H);
    if (!started) { phaseCtx3.moveTo(px, py); started = true; } else phaseCtx3.lineTo(px, py);
  }
  phaseCtx3.stroke(); phaseCtx3.setLineDash([]); phaseCtx3.restore();
  phaseCtx3.save(); phaseCtx3.lineWidth = 1;
  const gx = phaseCtx3.createLinearGradient(x0,y1,x1,y1);
  gx.addColorStop(0,_rgba('--cyan',0.85)); gx.addColorStop(1,_rgba('--cyan',0));
  phaseCtx3.strokeStyle = gx;
  phaseCtx3.beginPath(); phaseCtx3.moveTo(x0,y1); phaseCtx3.lineTo(x1,y1); phaseCtx3.stroke();
  const gy = phaseCtx3.createLinearGradient(x0,y1,x0,y0);
  gy.addColorStop(0,_rgba('--cyan',0.85)); gy.addColorStop(1,_rgba('--cyan',0));
  phaseCtx3.strokeStyle = gy;
  phaseCtx3.beginPath(); phaseCtx3.moveTo(x0,y1); phaseCtx3.lineTo(x0,y0); phaseCtx3.stroke();
  phaseCtx3.restore();
  const fontSize = Math.round(W * 0.05);
  phaseCtx3.save();
  phaseCtx3.font = `${fontSize}px 'EB Garamond', Georgia, serif`;
  phaseCtx3.textAlign='center'; phaseCtx3.textBaseline='top'; phaseCtx3.fillStyle=_c('--green-light');
  phaseCtx3.fillText('Resources', (x0+x1)/2, y1+fontSize*0.3);
  phaseCtx3.translate(x0-fontSize*0.3, (y0+y1)/2);
  phaseCtx3.rotate(-Math.PI/2); phaseCtx3.textAlign='center'; phaseCtx3.textBaseline='bottom';
  phaseCtx3.fillStyle=_c('--pink-light'); phaseCtx3.fillText('Foragers', 0, 0);
  phaseCtx3.restore();
}

function resizePhaseCanvas3() {
  if (!phaseCanvas3) return;
  const ctrl = document.getElementById('fa3-ctrl-panel');
  const side = ctrl ? Math.floor(ctrl.offsetWidth * 0.85) : 180;
  phaseCanvas3.width  = side;
  phaseCanvas3.height = side;
  if (accumCanvas3) resizeAccum3();
}

/* ── Main loop ── */
function loop3() {
  if (running) {
    for (let i = 0; i < stepsPerFrame; i++) simStep();
    drawPhase3();
  }
  drawScene();
  frameId = requestAnimationFrame(loop3);
}

/* ── Orbit controls ── */
function attachOrbitControls(canvas) {
  canvas.addEventListener('mousedown', e => {
    orbit.dragging = true;
    orbit.lastX = e.clientX;
    orbit.lastY = e.clientY;
    e.preventDefault();
  });
  window.addEventListener('mouseup', () => { orbit.dragging = false; });
  window.addEventListener('mousemove', e => {
    if (!orbit.dragging) return;
    const dx = e.clientX - orbit.lastX;
    const dy = e.clientY - orbit.lastY;
    orbit.theta -= dx * 0.01;
    orbit.phi    = Math.max(0.05, Math.min(Math.PI - 0.05, orbit.phi + dy * 0.01));
    orbit.lastX = e.clientX;
    orbit.lastY = e.clientY;
  });
  canvas.addEventListener('wheel', e => {
    orbit.dist = Math.max(50, Math.min(120, orbit.dist + e.deltaY * 0.3));
    e.preventDefault();
  }, { passive: false });
  // Touch
  let lastTouchDist = null;
  canvas.addEventListener('touchstart', e => {
    if (e.touches.length === 1) {
      orbit.dragging = true;
      orbit.lastX = e.touches[0].clientX;
      orbit.lastY = e.touches[0].clientY;
    } else if (e.touches.length === 2) {
      const dx = e.touches[0].clientX - e.touches[1].clientX;
      const dy = e.touches[0].clientY - e.touches[1].clientY;
      lastTouchDist = Math.sqrt(dx*dx+dy*dy);
    }
    e.preventDefault();
  }, { passive: false });
  canvas.addEventListener('touchend', () => { orbit.dragging = false; lastTouchDist = null; });
  canvas.addEventListener('touchmove', e => {
    if (e.touches.length === 1 && orbit.dragging) {
      const dx = e.touches[0].clientX - orbit.lastX;
      const dy = e.touches[0].clientY - orbit.lastY;
      orbit.theta -= dx * 0.01;
      orbit.phi    = Math.max(0.05, Math.min(Math.PI-0.05, orbit.phi + dy * 0.01));
      orbit.lastX = e.touches[0].clientX;
      orbit.lastY = e.touches[0].clientY;
    } else if (e.touches.length === 2 && lastTouchDist !== null) {
      const dx = e.touches[0].clientX - e.touches[1].clientX;
      const dy = e.touches[0].clientY - e.touches[1].clientY;
      const d = Math.sqrt(dx*dx+dy*dy);
      orbit.dist = Math.max(50, Math.min(120, orbit.dist - (d - lastTouchDist) * 0.5));
      lastTouchDist = d;
    }
    e.preventDefault();
  }, { passive: false });
}

/* ── CSS injection ── */
(function () {
  if (document.getElementById('fa3-styles')) return;
  const s = document.createElement('style');
  s.id = 'fa3-styles';
  s.textContent = `
    #fa3-ctrl-panel {
      overflow: hidden;
    }
    #fa3-ctrl-inner {
      display: flex;
      flex-direction: column;
      flex: 1;
      min-height: 0;
      overflow: hidden;
    }
    #fa3-ctrl-fixed { flex-shrink: 0; }
    #fa3-sliders-scroll { flex: 1; overflow-y: auto; min-height: 0; }
    #fa3-phase-section {
      flex-shrink: 0;
      display: flex;
      flex-direction: column;
      align-items: center;
      padding: 6px 10px 8px;
      border-top: 1px solid var(--border-dark);
    }
    #fa3-phase-title {
      font-size: calc(14px * var(--shell-fs, 1));
      text-transform: uppercase;
      letter-spacing: 1.5px;
      color: var(--text-dim);
      margin-bottom: 4px;
      flex-shrink: 0;
      align-self: flex-start;
    }
    #fa3-extinction-stamp {
      display: none;
      position: absolute;
      top: 50%; left: 50%;
      transform: translate(-50%, -50%) rotate(-20deg);
      font-size: 3em;
      color: var(--red);
      opacity: 0.7;
      pointer-events: none;
      letter-spacing: 4px;
      text-transform: uppercase;
      font-family: 'EB Garamond', Georgia, serif;
    }
  `;
  document.head.appendChild(s);
})();

/* ── AppletShell wiring ── */
const shell = new AppletShell({
  id:    'fa3',
  title: 'Forager 3D &mdash; Active Energy Harvesting',
  gap:   0,

  headerBtns: `<button class="applet-shell-header-btn" id="fa3-btn-reset">Reset</button><button class="applet-shell-header-btn" id="fa3-btn-run">Run</button>`,

  ctrlHTML: `
    <div id="fa3-ctrl-inner">
      <div id="fa3-ctrl-fixed">
      </div>
      <div id="fa3-sliders-scroll">
        <div class="applet-shell-ctrl-section">
          <div class="applet-shell-ctrl-title">Speed</div>
          <div class="applet-shell-slider-row">
            <span class="applet-shell-side">Slow</span>
            <input type="range" id="fa3-sl-speed" min="1" max="20" step="1" value="5">
            <span class="applet-shell-side">Fast</span>
          </div>
        </div>
        <div class="applet-shell-ctrl-section">
          <div class="applet-shell-ctrl-title">Sense Radius</div>
          <div class="applet-shell-slider-row">
            <span class="applet-shell-side">Small</span>
            <input type="range" id="fa3-sl-sense" min="1" max="20" step="0.5" value="6">
            <span class="applet-shell-side">Large</span>
          </div>
        </div>
        <div class="applet-shell-ctrl-section">
          <div class="applet-shell-ctrl-title">Velocity</div>
          <div class="applet-shell-slider-row">
            <span class="applet-shell-side">Slow</span>
            <input type="range" id="fa3-sl-vel" min="1" max="40" step="1" value="20">
            <span class="applet-shell-side">Fast</span>
          </div>
        </div>
        <div class="applet-shell-ctrl-section">
          <div class="applet-shell-ctrl-title">Tumble</div>
          <div class="applet-shell-slider-row">
            <span class="applet-shell-side">Ballistic</span>
            <input type="range" id="fa3-sl-sigma" min="0" max="2" step="0.05" value="0.1">
            <span class="applet-shell-side">Diffusive</span>
          </div>
        </div>
        <div class="applet-shell-ctrl-section">
          <div class="applet-shell-ctrl-title">Homogeneity</div>
          <div class="applet-shell-slider-row">
            <span class="applet-shell-side">Sparse</span>
            <input type="range" id="fa3-sl-xi" min="0.05" max="2.0" step="0.05" value="0.5">
            <span class="applet-shell-side">Dense</span>
          </div>
        </div>
      </div>
      <div id="fa3-phase-section">
        <div id="fa3-phase-title">Phase Plane</div>
        <canvas id="fa3-phase-canvas"></canvas>
      </div>
    </div>
  `,

  onOpen: function ({ canvas: c, W, H, S }) {
    const simW = W || S, simH = H || S;

    // Extinction stamp
    const simPanel = document.getElementById('fa3-sim-panel');
    if (simPanel && !document.getElementById('fa3-extinction-stamp')) {
      const stamp = document.createElement('div');
      stamp.id = 'fa3-extinction-stamp';
      stamp.textContent = 'EXTINCTION';
      simPanel.appendChild(stamp);
    }

    simCanvas3 = c;
    simCanvas3.width  = simW;
    simCanvas3.height = simH;

    initGL(simCanvas3);
    attachOrbitControls(simCanvas3);

    phaseCanvas3 = document.getElementById('fa3-phase-canvas');
    phaseCtx3    = phaseCanvas3.getContext('2d');

    setTimeout(() => {
      resizePhaseCanvas3();
      if (!initialized) { fa3InitSim(); initialized = true; }
      initAccum3();
      running = true;
      const btn = document.getElementById('fa3-btn-run');
      if (btn) { btn.textContent = 'Pause'; btn.classList.add('active'); }
      if (!frameId) frameId = requestAnimationFrame(loop3);
    }, 80);
  },

  onClose: function () {
    running = false;
    if (frameId) { cancelAnimationFrame(frameId); frameId = null; }
    const btn = document.getElementById('fa3-btn-run');
    if (btn) { btn.textContent = 'Run'; btn.classList.remove('active'); }
  },

  onResize: function ({ canvas: c, W, H, S }) {
    const simW = W || S, simH = H || S;
    simCanvas3 = c;
    simCanvas3.width  = simW;
    simCanvas3.height = simH;
    resizePhaseCanvas3();
  },
});

window.fa3Open  = () => shell.open();
window.fa3Close = () => shell.close();

document.getElementById('fa3-btn-run').addEventListener('click', () => {
  running = !running;
  const btn = document.getElementById('fa3-btn-run');
  btn.textContent = running ? 'Pause' : 'Run';
  btn.classList.toggle('active', running);
  if (running && !frameId) frameId = requestAnimationFrame(loop3);
});

document.getElementById('fa3-btn-reset').addEventListener('click', () => {
  fa3InitSim(); clearAccum3(); drawPhase3();
});

document.getElementById('fa3-sl-speed').addEventListener('input', function () {
  stepsPerFrame = parseInt(this.value);
});
document.getElementById('fa3-sl-sense').addEventListener('input', function () {
  R_sense = parseFloat(this.value);
});
document.getElementById('fa3-sl-vel').addEventListener('input', function () {
  vel = parseFloat(this.value);
});
document.getElementById('fa3-sl-sigma').addEventListener('input', function () {
  sigma = parseFloat(this.value);
});
document.getElementById('fa3-sl-xi').addEventListener('input', function () {
  Xi = parseFloat(this.value); epsilon = epsilonFromXi(Xi);
  initResources();
});

})();
