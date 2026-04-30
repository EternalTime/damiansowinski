(function () {

// ── Inject kuramoto-specific styles ──────────────────────────────────────────
(function() {
  const s = document.createElement('style');
  s.textContent = `
    #km-hist-section {
      flex: 1; display: flex; flex-direction: column;
      padding: 14px 16px 10px; min-height: 0;
    }
    #km-hist-canvas { flex: 1; display: block; width: 100%; min-height: 0; }
  `;
  document.head.appendChild(s);
})();

const _cs = getComputedStyle(document.documentElement);
const _c   = n => _cs.getPropertyValue(n).trim();
const _rgb  = n => { const h = _c(n).replace('#',''); const v = parseInt(h,16); return [(v>>16)&0xFF,(v>>8)&0xFF,v&0xFF]; };
const _rgba = (n, a) => { const [r,g,b] = _rgb(n); return `rgba(${r},${g},${b},${a})`; };

const GRAPHS = [
  {
    name: "Full",
    N: 25,
    Kmin: 0.0,
    Kmax: 4.0,
    camera: { radius: 4.05, theta: 0, phi: 1.08 },
    layout: { seed: 1, idealEdge: 1.4, relaxFrames: 300, edgeOpacity: 0.40, edgeWidth: 1.0, edgeStyle: 'line', nodeRadiusMin: 0.13, nodeRadiusMax: 0.13 },
    edges: [[0,1],[0,2],[0,3],[0,4],[0,5],[0,6],[0,7],[0,8],[0,9],[0,10],[0,11],[0,12],[0,13],[0,14],[0,15],[0,16],[0,17],[0,18],[0,19],[0,20],[0,21],[0,22],[0,23],[0,24],[1,2],[1,3],[1,4],[1,5],[1,6],[1,7],[1,8],[1,9],[1,10],[1,11],[1,12],[1,13],[1,14],[1,15],[1,16],[1,17],[1,18],[1,19],[1,20],[1,21],[1,22],[1,23],[1,24],[2,3],[2,4],[2,5],[2,6],[2,7],[2,8],[2,9],[2,10],[2,11],[2,12],[2,13],[2,14],[2,15],[2,16],[2,17],[2,18],[2,19],[2,20],[2,21],[2,22],[2,23],[2,24],[3,4],[3,5],[3,6],[3,7],[3,8],[3,9],[3,10],[3,11],[3,12],[3,13],[3,14],[3,15],[3,16],[3,17],[3,18],[3,19],[3,20],[3,21],[3,22],[3,23],[3,24],[4,5],[4,6],[4,7],[4,8],[4,9],[4,10],[4,11],[4,12],[4,13],[4,14],[4,15],[4,16],[4,17],[4,18],[4,19],[4,20],[4,21],[4,22],[4,23],[4,24],[5,6],[5,7],[5,8],[5,9],[5,10],[5,11],[5,12],[5,13],[5,14],[5,15],[5,16],[5,17],[5,18],[5,19],[5,20],[5,21],[5,22],[5,23],[5,24],[6,7],[6,8],[6,9],[6,10],[6,11],[6,12],[6,13],[6,14],[6,15],[6,16],[6,17],[6,18],[6,19],[6,20],[6,21],[6,22],[6,23],[6,24],[7,8],[7,9],[7,10],[7,11],[7,12],[7,13],[7,14],[7,15],[7,16],[7,17],[7,18],[7,19],[7,20],[7,21],[7,22],[7,23],[7,24],[8,9],[8,10],[8,11],[8,12],[8,13],[8,14],[8,15],[8,16],[8,17],[8,18],[8,19],[8,20],[8,21],[8,22],[8,23],[8,24],[9,10],[9,11],[9,12],[9,13],[9,14],[9,15],[9,16],[9,17],[9,18],[9,19],[9,20],[9,21],[9,22],[9,23],[9,24],[10,11],[10,12],[10,13],[10,14],[10,15],[10,16],[10,17],[10,18],[10,19],[10,20],[10,21],[10,22],[10,23],[10,24],[11,12],[11,13],[11,14],[11,15],[11,16],[11,17],[11,18],[11,19],[11,20],[11,21],[11,22],[11,23],[11,24],[12,13],[12,14],[12,15],[12,16],[12,17],[12,18],[12,19],[12,20],[12,21],[12,22],[12,23],[12,24],[13,14],[13,15],[13,16],[13,17],[13,18],[13,19],[13,20],[13,21],[13,22],[13,23],[13,24],[14,15],[14,16],[14,17],[14,18],[14,19],[14,20],[14,21],[14,22],[14,23],[14,24],[15,16],[15,17],[15,18],[15,19],[15,20],[15,21],[15,22],[15,23],[15,24],[16,17],[16,18],[16,19],[16,20],[16,21],[16,22],[16,23],[16,24],[17,18],[17,19],[17,20],[17,21],[17,22],[17,23],[17,24],[18,19],[18,20],[18,21],[18,22],[18,23],[18,24],[19,20],[19,21],[19,22],[19,23],[19,24],[20,21],[20,22],[20,23],[20,24],[21,22],[21,23],[21,24],[22,23],[22,24],[23,24]]
  },
  {
    name: "Ring",
    N: 60,
    Kmin: 0.0,
    Kmax: 8.0,
    camera: { radius: 4.54, theta: -3.32, phi: 0.835 },
    layout: { seed: 2, idealEdge: 0.1, relaxFrames: 800, edgeOpacity: 0.75, edgeWidth: 3.0, edgeStyle: 'cylinder', edgeRadius: 0.025, nodeRadiusMin: 0.13, nodeRadiusMax: 0.13 },
    edges: [[0,1],[0,19],[0,20],[1,2],[1,21],[2,3],[2,22],[3,4],[3,23],[4,5],[4,24],[5,6],[5,25],[6,7],[6,26],[7,8],[7,27],[8,9],[8,28],[9,10],[9,29],[10,11],[10,30],[11,12],[11,31],[12,13],[12,32],[13,14],[13,33],[14,15],[14,34],[15,16],[15,35],[16,17],[16,36],[17,18],[17,37],[18,19],[18,38],[19,39],[20,21],[20,39],[20,40],[21,22],[21,41],[22,23],[22,42],[23,24],[23,43],[24,25],[24,44],[25,26],[25,45],[26,27],[26,46],[27,28],[27,47],[28,29],[28,48],[29,30],[29,49],[30,31],[30,50],[31,32],[31,51],[32,33],[32,52],[33,34],[33,53],[34,35],[34,54],[35,36],[35,55],[36,37],[36,56],[37,38],[37,57],[38,39],[38,58],[39,59],[40,41],[40,59],[41,42],[42,43],[43,44],[44,45],[45,46],[46,47],[47,48],[48,49],[49,50],[50,51],[51,52],[52,53],[53,54],[54,55],[55,56],[56,57],[57,58],[58,59]]
  },
  {
    name: "Scale-free",
    N: 80,
    Kmin: 0.0,
    Kmax: 30.0,
    camera: { radius: 6.72, theta: -3.74, phi: 0.7 },
    layout: { seed: 3, idealEdge: 0.15, relaxFrames: 300, edgeOpacity: 0.60, edgeWidth: 2.0, edgeStyle: 'cylinder', edgeRadius: 0.025, nodeRadiusMin: 0.15, nodeRadiusMax: 0.6 },
    edges: [[0,1],[0,2],[0,3],[0,4],[0,6],[0,7],[0,8],[0,9],[0,11],[0,13],[0,14],[0,15],[0,22],[0,23],[0,24],[0,34],[0,44],[0,45],[0,47],[0,50],[0,62],[0,63],[0,65],[0,67],[0,68],[0,79],[1,2],[1,3],[1,29],[1,45],[1,53],[2,4],[2,5],[2,6],[2,7],[2,12],[2,14],[2,22],[2,23],[2,24],[2,25],[2,27],[2,28],[2,33],[2,35],[2,56],[2,59],[2,69],[2,77],[3,8],[3,9],[3,10],[3,18],[3,26],[3,35],[3,39],[3,40],[3,41],[3,43],[3,72],[3,74],[4,5],[4,26],[4,30],[4,34],[4,37],[4,52],[4,55],[4,77],[4,79],[6,10],[6,12],[6,16],[6,17],[6,19],[6,38],[6,58],[6,71],[7,11],[7,42],[7,48],[7,62],[7,65],[7,74],[8,20],[8,21],[9,39],[9,47],[10,17],[10,31],[10,78],[11,15],[11,16],[11,21],[11,25],[11,28],[11,44],[11,54],[12,13],[12,31],[12,40],[13,19],[13,49],[13,72],[14,41],[14,56],[15,20],[15,29],[15,32],[15,52],[15,76],[16,37],[16,67],[16,70],[17,18],[17,46],[17,60],[17,61],[20,33],[22,76],[23,27],[23,32],[23,53],[23,59],[23,70],[23,73],[24,51],[24,58],[27,36],[27,51],[28,78],[29,30],[29,46],[32,36],[32,38],[33,42],[33,75],[36,60],[37,50],[38,61],[41,73],[42,43],[42,68],[42,69],[43,55],[44,57],[45,64],[46,49],[47,48],[47,57],[49,54],[49,64],[50,66],[54,63],[54,71],[59,66],[74,75]]
  },
  {
    name: "Small-world",
    N: 80,
    Kmin: 0.0,
    Kmax: 12.0,
    camera: { radius: 6.89, theta: 1.221, phi: 0.595 },
    layout: { seed: 4, idealEdge: 0.2, relaxFrames: 400, edgeOpacity: 0.50, edgeWidth: 1.5, edgeStyle: 'cylinder', edgeRadius: 0.025, nodeRadiusMin: 0.1, nodeRadiusMax: 0.20 },
    edges: [[0,1],[0,5],[0,55],[0,78],[0,79],[1,2],[1,3],[1,79],[2,3],[2,4],[3,4],[3,5],[4,5],[4,6],[5,6],[5,7],[6,7],[6,8],[6,32],[7,8],[7,24],[8,9],[8,10],[9,10],[9,11],[10,11],[10,12],[11,12],[11,13],[12,13],[12,14],[13,14],[13,15],[14,15],[14,16],[15,16],[15,17],[16,17],[16,18],[17,18],[17,19],[18,19],[18,20],[18,72],[19,20],[19,21],[20,21],[20,22],[21,22],[21,23],[22,23],[22,24],[23,24],[23,25],[24,25],[24,26],[25,26],[25,27],[26,27],[26,28],[26,35],[27,29],[27,46],[27,53],[28,29],[28,30],[29,31],[29,42],[30,31],[30,32],[31,32],[31,33],[32,33],[33,34],[33,35],[34,35],[34,36],[35,36],[36,37],[36,38],[37,39],[37,72],[38,39],[38,40],[39,40],[39,41],[40,41],[40,42],[41,42],[41,43],[42,43],[42,44],[43,44],[43,45],[44,45],[44,46],[45,46],[45,47],[46,48],[47,48],[47,49],[48,50],[48,71],[49,50],[49,51],[50,52],[50,56],[51,52],[51,53],[52,53],[52,54],[53,54],[53,55],[54,55],[54,56],[55,57],[56,57],[56,58],[56,76],[57,58],[57,59],[58,59],[58,60],[59,60],[59,61],[60,61],[60,62],[60,66],[61,62],[61,63],[62,63],[62,64],[63,64],[63,65],[64,65],[64,66],[65,66],[65,67],[66,67],[67,68],[67,69],[68,69],[68,70],[69,70],[69,71],[70,71],[70,72],[71,72],[71,73],[72,73],[73,74],[73,75],[74,75],[74,76],[75,76],[75,77],[76,77],[77,78],[77,79],[78,79]]
  },
  {
    name: "Random",
    N: 60,
    Kmin: 0.0,
    Kmax: 18.0,
    camera: { radius: 6.5, theta: 0.66, phi: 0.113 },
    layout: { seed: 5, idealEdge: 0.2, relaxFrames: 350, edgeOpacity: 0.55, edgeWidth: 1.5, edgeStyle: 'cylinder', edgeRadius: 0.02, nodeRadiusMin: 0.1, nodeRadiusMax: 0.3 },
    edges: [[0,4],[0,7],[0,9],[0,11],[0,22],[0,34],[0,35],[0,57],[1,13],[1,21],[1,25],[1,31],[1,49],[1,58],[2,10],[2,11],[2,15],[2,16],[2,20],[2,32],[2,38],[2,43],[3,9],[3,10],[3,30],[3,57],[4,28],[4,31],[5,15],[5,17],[5,40],[5,43],[5,47],[5,53],[6,12],[6,13],[6,16],[6,20],[6,22],[6,34],[6,43],[6,45],[7,14],[7,16],[7,20],[7,38],[7,44],[7,46],[7,48],[7,49],[8,13],[8,30],[8,40],[8,46],[8,49],[9,12],[9,18],[9,19],[9,33],[9,54],[9,55],[10,14],[10,29],[10,47],[10,53],[11,42],[11,57],[11,58],[12,16],[12,33],[12,44],[12,50],[13,16],[13,24],[13,52],[13,53],[14,19],[14,25],[14,36],[14,43],[14,47],[14,48],[14,49],[15,18],[15,23],[15,43],[15,47],[15,48],[15,55],[16,49],[17,31],[17,39],[18,20],[18,28],[18,32],[18,49],[19,30],[19,33],[19,46],[20,22],[20,24],[20,28],[20,45],[20,59],[21,50],[21,58],[21,59],[22,31],[22,32],[22,52],[23,26],[23,30],[23,32],[23,44],[23,48],[23,57],[24,38],[25,29],[25,32],[25,37],[25,38],[25,43],[25,49],[25,56],[26,36],[26,46],[26,58],[27,30],[27,46],[27,55],[27,57],[29,50],[31,38],[32,34],[32,37],[33,49],[34,36],[34,43],[34,44],[35,51],[35,58],[36,44],[36,47],[36,56],[38,42],[38,46],[38,48],[39,47],[40,43],[40,45],[40,51],[41,44],[41,45],[41,49],[41,52],[41,57],[42,45],[42,53],[42,54],[42,57],[43,47],[44,49],[45,53],[45,59],[46,57],[48,49],[49,58],[50,59],[51,52],[52,58],[53,55],[53,56],[54,57]]
  }
];

let gIdx   = 0;
let phases = null, omegas = null, nbrs = null, N = 0;
let K      = 1.0;
let sigma  = 2 * Math.PI * 0.125;
const STEPS = 4;
let lastTime = null;

const N_BINS = 30;
let smoothBins = new Float64Array(N_BINS);
const HIST_ALPHA = 0.12;
let hctx;

let renderer, scene, camera;
let nodeMeshes = [], edgeLines = null, edgeMeshes = [];
let simCanvas;

let orbit = { dragging: false, lastX: 0, lastY: 0, theta: 0.4, phi: 1.1, radius: 9.0 };
let nodePos = null, nodeDeg = null;
let relaxIter = 0, RELAX_FRAMES = 120;

const [_TLR, _TLG, _TLB] = _rgb('--teal-light');
const [_TDR, _TDG, _TDB] = _rgb('--teal-dark');
const [_CYR, _CYG, _CYB] = _rgb('--cyan');
const [_PDR, _PDG, _PDB] = _rgb('--pink-dark');
const [_PLR, _PLG, _PLB] = _rgb('--pink-light');

const PALETTE_STOPS = [
  { u: 0.00, c: [_TLR/255, _TLG/255, _TLB/255] },
  { u: 0.12, c: [_TDR/255, _TDG/255, _TDB/255] },
  { u: 0.50, c: [_CYR/255, _CYG/255, _CYB/255] },
  { u: 0.88, c: [_PDR/255, _PDG/255, _PDB/255] },
  { u: 1.00, c: [_PLR/255, _PLG/255, _PLB/255] },
];
function angleToColor(th) {
  const TWO_PI = 2 * Math.PI;
  const u = ((th % TWO_PI) + TWO_PI) % TWO_PI / TWO_PI;
  let i = 0;
  while (i < PALETTE_STOPS.length - 2 && PALETTE_STOPS[i+1].u <= u) i++;
  const s0 = PALETTE_STOPS[i], s1 = PALETTE_STOPS[Math.min(i+1, PALETTE_STOPS.length-1)];
  const t = (u - s0.u) / (s1.u - s0.u);
  const c0 = s0.c, c1 = s1.c;
  return new THREE.Color(c0[0]+t*(c1[0]-c0[0]), c0[1]+t*(c1[1]-c0[1]), c0[2]+t*(c1[2]-c0[2]));
}

function buildGraph(idx) {
  gIdx = idx;
  const g = GRAPHS[idx];
  N = g.N;
  nbrs = Array.from({length: N}, () => []);
  for (const [a, b] of g.edges) { nbrs[a].push(b); nbrs[b].push(a); }
  nbrs = nbrs.map(arr => new Int32Array(arr));
  let seed = g.layout.seed || 42;
  function seededRand() {
    seed |= 0; seed = seed + 0x6D2B79F5 | 0;
    let t = Math.imul(seed ^ seed >>> 15, 1 | seed);
    t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t;
    return ((t ^ t >>> 14) >>> 0) / 4294967296;
  }
  const pos = Array.from({length: N}, () => {
    const th = Math.acos(2*seededRand()-1);
    const ph = seededRand()*2*Math.PI;
    return [Math.sin(th)*Math.cos(ph), Math.sin(th)*Math.sin(ph), Math.cos(th)];
  });
  nodeDeg = new Float32Array(N);
  for (const [a, b] of g.edges) { nodeDeg[a]++; nodeDeg[b]++; }
  const minDeg2 = Math.max(1, Math.min(...nodeDeg));
  for (let i = 0; i < N; i++) nodeDeg[i] = (nodeDeg[i] / minDeg2) ** 2;
  nodePos = pos; relaxIter = 0; RELAX_FRAMES = g.layout.relaxFrames;
  const cam = g.camera;
  orbit.radius = cam.radius; orbit.theta = cam.theta; orbit.phi = cam.phi;
  if (camera) updateCamera();
  return pos;
}

function relaxStep(pos, edges, t, idealEdge, nodeDeg) {
  const n = pos.length, stepSz = 0.12 * (1 - t);
  const force = Array.from({length: n}, () => [0,0,0]);
  const baseRepulse = idealEdge * idealEdge * 0.5;
  for (let i = 0; i < n; i++) {
    for (let j = i+1; j < n; j++) {
      const dx = pos[i][0]-pos[j][0], dy = pos[i][1]-pos[j][1], dz = pos[i][2]-pos[j][2];
      const d = Math.sqrt(dx*dx+dy*dy+dz*dz) + 1e-6;
      const charge = nodeDeg ? (nodeDeg[i]+nodeDeg[j])*0.5 : 1.0;
      const f = baseRepulse * charge / (d*d);
      force[i][0]+=f*dx/d; force[i][1]+=f*dy/d; force[i][2]+=f*dz/d;
      force[j][0]-=f*dx/d; force[j][1]-=f*dy/d; force[j][2]-=f*dz/d;
    }
  }
  for (const [a, b] of edges) {
    const dx = pos[b][0]-pos[a][0], dy = pos[b][1]-pos[a][1], dz = pos[b][2]-pos[a][2];
    const d = Math.sqrt(dx*dx+dy*dy+dz*dz) + 1e-6;
    const f = (d - idealEdge) * 0.3;
    force[a][0]+=f*dx/d; force[a][1]+=f*dy/d; force[a][2]+=f*dz/d;
    force[b][0]-=f*dx/d; force[b][1]-=f*dy/d; force[b][2]-=f*dz/d;
  }
  for (let i = 0; i < n; i++) {
    pos[i][0]+=stepSz*force[i][0]; pos[i][1]+=stepSz*force[i][1]; pos[i][2]+=stepSz*force[i][2];
  }
}

function initPhases() {
  phases = new Float64Array(N);
  for (let i = 0; i < N; i++) phases[i] = Math.random() * 2 * Math.PI;
}

function initOmegas() {
  omegas = new Float64Array(N);
  const TWO_PI = 2 * Math.PI, OMEGA_MEAN = Math.PI;
  for (let i = 0; i < N; i += 2) {
    const u1 = Math.random(), u2 = Math.random();
    const mag = sigma * Math.sqrt(-2 * Math.log(u1 + 1e-15));
    omegas[i]   = OMEGA_MEAN + mag * Math.cos(TWO_PI*u2);
    if (i+1 < N) omegas[i+1] = OMEGA_MEAN + mag * Math.sin(TWO_PI*u2);
  }
}

function step(dt) {
  const TWO_PI = 2 * Math.PI, KoverN = K / N;
  const newPhases = new Float64Array(N);
  for (let i = 0; i < N; i++) {
    let coupling = 0;
    const nb = nbrs[i];
    for (let k = 0; k < nb.length; k++) coupling += Math.sin(phases[nb[k]] - phases[i]);
    newPhases[i] = phases[i] + dt * (omegas[i] + KoverN * coupling);
  }
  for (let i = 0; i < N; i++) phases[i] = ((newPhases[i] % (2*Math.PI)) + 2*Math.PI) % (2*Math.PI);
}

function initThree(pos) {
  if (!renderer) {
    simCanvas = document.getElementById('km-canvas');
    renderer = new THREE.WebGLRenderer({ canvas: simCanvas, antialias: true });
    renderer.setPixelRatio(window.devicePixelRatio);
    renderer.setClearColor(new THREE.Color(_c('--bg-void')), 1);
    camera = new THREE.PerspectiveCamera(55, 1, 0.1, 100);
    updateCamera();
    setupOrbitControls();
  }
  if (scene) {
    scene.traverse(obj => {
      if (obj.geometry) obj.geometry.dispose();
      if (obj.material) obj.material.dispose();
    });
  }
  scene = new THREE.Scene();
  nodeMeshes = []; edgeLines = null; edgeMeshes = [];
  scene.add(new THREE.AmbientLight(0xffffff, 0.4));
  const dLight = new THREE.DirectionalLight(0xffffff, 0.8);
  dLight.position.set(5, 8, 6); scene.add(dLight);
  const layout = GRAPHS[gIdx].layout;
  const deg = new Int32Array(N);
  for (const [a, b] of GRAPHS[gIdx].edges) { deg[a]++; deg[b]++; }
  const minDeg = Math.min(...deg), maxDeg = Math.max(...deg);
  const rMin = layout.nodeRadiusMin || 0.10, rMax = layout.nodeRadiusMax || 0.22;
  for (let i = 0; i < N; i++) {
    const t = maxDeg > minDeg ? (deg[i]-minDeg)/(maxDeg-minDeg) : 0.5;
    const r = rMin + t * (rMax - rMin);
    const geo = new THREE.SphereGeometry(r, 16, 12);
    const mat = new THREE.MeshPhongMaterial({ color: angleToColor(phases[i]), emissive: angleToColor(phases[i]), emissiveIntensity: 0.45, shininess: 60 });
    const mesh = new THREE.Mesh(geo, mat);
    mesh.position.set(pos[i][0], pos[i][1], pos[i][2]);
    scene.add(mesh); nodeMeshes.push(mesh);
  }
  if (layout.edgeStyle === 'cylinder') {
    const cylGeo = new THREE.CylinderGeometry(layout.edgeRadius, layout.edgeRadius, 1, 6, 1);
    for (const [a, b] of GRAPHS[gIdx].edges) {
      const ca = angleToColor(phases[a]), cb = angleToColor(phases[b]);
      const col = new THREE.Color((ca.r+cb.r)/2, (ca.g+cb.g)/2, (ca.b+cb.b)/2);
      const mat = new THREE.MeshPhongMaterial({ color: col, emissive: col, emissiveIntensity: 0.3, transparent: true, opacity: layout.edgeOpacity });
      const cyl = new THREE.Mesh(cylGeo, mat);
      positionCylinder(cyl, pos[a], pos[b]);
      scene.add(cyl); edgeMeshes.push({ mesh: cyl, a, b });
    }
  } else {
    const edgeGeo = new THREE.BufferGeometry();
    const verts = [], edgeColors = [];
    for (const [a, b] of GRAPHS[gIdx].edges) {
      verts.push(...pos[a], ...pos[b]);
      const ca = angleToColor(phases[a]), cb = angleToColor(phases[b]);
      edgeColors.push(ca.r, ca.g, ca.b, cb.r, cb.g, cb.b);
    }
    edgeGeo.setAttribute('position', new THREE.Float32BufferAttribute(verts, 3));
    edgeGeo.setAttribute('color',    new THREE.Float32BufferAttribute(edgeColors, 3));
    const edgeMat = new THREE.LineBasicMaterial({ vertexColors: true, transparent: true, opacity: layout.edgeOpacity, linewidth: layout.edgeWidth });
    edgeLines = new THREE.LineSegments(edgeGeo, edgeMat);
    scene.add(edgeLines);
  }
  resizeRenderer();
}

function positionCylinder(cyl, pa, pb) {
  const dx = pb[0]-pa[0], dy = pb[1]-pa[1], dz = pb[2]-pa[2];
  const len = Math.sqrt(dx*dx+dy*dy+dz*dz);
  cyl.scale.y = len;
  cyl.position.set((pa[0]+pb[0])/2, (pa[1]+pb[1])/2, (pa[2]+pb[2])/2);
  const dir = new THREE.Vector3(dx/len, dy/len, dz/len);
  cyl.setRotationFromQuaternion(new THREE.Quaternion().setFromUnitVectors(new THREE.Vector3(0,1,0), dir));
}

function updateCamera() {
  const { theta, phi, radius } = orbit;
  camera.position.set(radius*Math.sin(phi)*Math.cos(theta), radius*Math.cos(phi), radius*Math.sin(phi)*Math.sin(theta));
  camera.lookAt(0, 0, 0);
}

function resizeRenderer() {
  const w = simCanvas.clientWidth, h = simCanvas.clientHeight;
  renderer.setSize(w, h, false);
  camera.aspect = w / h;
  camera.updateProjectionMatrix();
}

function setupOrbitControls() {
  simCanvas.addEventListener('pointerdown', e => { orbit.dragging=true; orbit.lastX=e.clientX; orbit.lastY=e.clientY; simCanvas.setPointerCapture(e.pointerId); });
  simCanvas.addEventListener('pointermove', e => {
    if (!orbit.dragging) return;
    orbit.theta -= (e.clientX-orbit.lastX)*0.008;
    orbit.phi = Math.max(0.1, Math.min(Math.PI-0.1, orbit.phi+(e.clientY-orbit.lastY)*0.008));
    orbit.lastX=e.clientX; orbit.lastY=e.clientY; updateCamera();
  });
  simCanvas.addEventListener('pointerup', () => { orbit.dragging=false; });
  simCanvas.addEventListener('wheel', e => { e.preventDefault(); orbit.radius=Math.max(2,Math.min(18,orbit.radius+e.deltaY*0.01)); updateCamera(); }, { passive: false });
}

function renderHistogram() {
  const hc = document.getElementById('km-hist-canvas');
  const W = hc.clientWidth||200, H = hc.clientHeight||120;
  if (hc.width!==W||hc.height!==H) { hc.width=W; hc.height=H; }
  const TWO_PI = 2*Math.PI;
  const raw = new Float64Array(N_BINS);
  for (let i = 0; i < N; i++) {
    const b = Math.min((phases[i]/TWO_PI*N_BINS)|0, N_BINS-1);
    raw[b]++;
  }
  for (let b = 0; b < N_BINS; b++) raw[b] /= N;
  for (let b = 0; b < N_BINS; b++) smoothBins[b] += HIST_ALPHA*(raw[b]-smoothBins[b]);
  const ymax = (1/N_BINS)*8;
  const PL=4,PR=4,PT=6,PB=4, pw=W-PL-PR, ph=H-PT-PB, bw=pw/N_BINS;
  hctx.clearRect(0,0,W,H);
  for (let b = 0; b < N_BINS; b++) {
    const bh = Math.min(smoothBins[b]/ymax,1)*ph;
    const col = angleToColor((b+0.5)/N_BINS*TWO_PI);
    hctx.fillStyle = `rgba(${Math.round(col.r*255)},${Math.round(col.g*255)},${Math.round(col.b*255)},0.8)`;
    hctx.fillRect(PL+b*bw, PT+ph-bh, bw-1, bh);
  }
}

let running = false, frameId = null;

function loop(now) {
  if (running) {
    if (lastTime===null) lastTime=now;
    const dtReal = Math.min((now-lastTime)/1000, 0.05);
    lastTime = now;
    const dtSub = dtReal/STEPS;
    if (relaxIter < RELAX_FRAMES) {
      relaxStep(nodePos, GRAPHS[gIdx].edges, relaxIter/RELAX_FRAMES, GRAPHS[gIdx].layout.idealEdge, nodeDeg);
      for (let i = 0; i < N; i++) nodeMeshes[i].position.set(nodePos[i][0], nodePos[i][1], nodePos[i][2]);
      if (edgeLines) {
        const posAttr = edgeLines.geometry.attributes.position;
        let vi = 0;
        for (const [a, b] of GRAPHS[gIdx].edges) {
          posAttr.array[vi++]=nodePos[a][0]; posAttr.array[vi++]=nodePos[a][1]; posAttr.array[vi++]=nodePos[a][2];
          posAttr.array[vi++]=nodePos[b][0]; posAttr.array[vi++]=nodePos[b][1]; posAttr.array[vi++]=nodePos[b][2];
        }
        posAttr.needsUpdate = true;
      } else {
        for (const { mesh, a, b } of edgeMeshes) positionCylinder(mesh, nodePos[a], nodePos[b]);
      }
      relaxIter++;
    }
    for (let s = 0; s < STEPS; s++) step(dtSub);
    for (let i = 0; i < N; i++) {
      const col = angleToColor(phases[i]);
      nodeMeshes[i].material.color.copy(col);
      nodeMeshes[i].material.emissive.copy(col);
    }
    if (edgeLines) {
      const colAttr = edgeLines.geometry.attributes.color;
      let ci = 0;
      for (const [a, b] of GRAPHS[gIdx].edges) {
        const ca = angleToColor(phases[a]), cb = angleToColor(phases[b]);
        colAttr.array[ci++]=ca.r; colAttr.array[ci++]=ca.g; colAttr.array[ci++]=ca.b;
        colAttr.array[ci++]=cb.r; colAttr.array[ci++]=cb.g; colAttr.array[ci++]=cb.b;
      }
      colAttr.needsUpdate = true;
    } else {
      for (const { mesh, a, b } of edgeMeshes) {
        const ca = angleToColor(phases[a]), cb = angleToColor(phases[b]);
        const col = new THREE.Color((ca.r+cb.r)/2,(ca.g+cb.g)/2,(ca.b+cb.b)/2);
        mesh.material.color.copy(col); mesh.material.emissive.copy(col);
      }
    }
    renderer.render(scene, camera);
    renderHistogram();
    // const ro = document.getElementById('km-cam-readout');
    // if (ro) ro.textContent = `θ: ${orbit.theta.toFixed(3)}  φ: ${orbit.phi.toFixed(3)}  r: ${orbit.radius.toFixed(2)}`;
  }
  frameId = requestAnimationFrame(loop);
}

function buildButtons() {
  const row = document.getElementById('km-graph-btns');
  row.innerHTML = '';
  GRAPHS.forEach((g, i) => {
    const btn = document.createElement('button');
    btn.className = 'applet-shell-btn km-graph-btn' + (i===gIdx ? ' active' : '');
    btn.textContent = g.name;
    btn.addEventListener('click', () => switchGraph(i));
    row.appendChild(btn);
  });
}

function switchGraph(idx) {
  if (idx===gIdx && running) return;
  document.querySelectorAll('.km-graph-btn').forEach((b,i) => b.classList.toggle('active', i===idx));
  const oldG = GRAPHS[gIdx], newG = GRAPHS[idx];
  const t = (K-oldG.Kmin)/(oldG.Kmax-oldG.Kmin);
  K = newG.Kmin + t*(newG.Kmax-newG.Kmin);
  document.getElementById('km-K').value = (K-newG.Kmin)/(newG.Kmax-newG.Kmin);
  const pos = buildGraph(idx);
  initPhases(); initOmegas();
  smoothBins = new Float64Array(N_BINS);
  lastTime = null;
  initThree(pos);
}

/* ── Shell wiring ── */
const shell = new AppletShell({
  id:    'km',
  title: 'Kuramoto Model &mdash; Synchronisation',
  gap:   0,

  ctrlHTML: `
    <div class="applet-shell-ctrl-section">
      <div class="applet-shell-ctrl-title">Actions</div>
      <div class="applet-shell-btn-row">
        <button class="applet-shell-btn" onclick="kmRandomiseFreqs()">Reset</button>
        <button class="applet-shell-btn" id="km-pause-btn" onclick="kmTogglePause()">Pause</button>
      </div>
    </div>
    <div class="applet-shell-ctrl-section">
      <div class="applet-shell-ctrl-title">Graph</div>
      <div class="applet-shell-btn-row" id="km-graph-btns"></div>
    </div>
    <div class="applet-shell-ctrl-section">
      <div class="applet-shell-ctrl-title">Coupling</div>
      <div class="applet-shell-slider-row">
        <span class="applet-shell-side">Weak</span>
        <input type="range" id="km-K" min="0" max="1" step="0.001" value="0.5">
        <span class="applet-shell-side">Strong</span>
      </div>
    </div>
    <div class="applet-shell-ctrl-section">
      <div class="applet-shell-ctrl-title">Frequency Spread</div>
      <div class="applet-shell-slider-row">
        <span class="applet-shell-side">Narrow</span>
        <input type="range" id="km-sigma" min="0.01" max="1" step="0.01" value="0.5">
        <span class="applet-shell-side">Wide</span>
      </div>
    </div>
    <div id="km-hist-section">
      <div class="applet-shell-ctrl-title">Phase Distribution</div>
      <canvas id="km-hist-canvas"></canvas>
    </div>
  `,

  onOpen: function ({ canvas: c, S }) {
    // // Inject camera readout div into the sim panel
    // const simPanel = document.getElementById('km-sim-panel') || document.getElementById('km-asm-sim');
    // if (simPanel && !document.getElementById('km-cam-readout')) {
    //   const ro = document.createElement('div');
    //   ro.id = 'km-cam-readout';
    //   ro.style.cssText = 'position:absolute;bottom:8px;left:10px;font-family:monospace;font-size:11px;color:' + _rgba('--white', 0.55) + ';pointer-events:none;';
    //   simPanel.appendChild(ro);
    // }

    buildButtons();
    const pos = buildGraph(0);
    initPhases(); initOmegas();
    smoothBins = new Float64Array(N_BINS);
    const g = GRAPHS[0];
    const t = parseFloat(document.getElementById('km-K').value);
    K = g.Kmin + t*(g.Kmax-g.Kmin);

    const pb = document.getElementById('km-pause-btn');
    if (pb) { pb.textContent = 'Pause'; pb.classList.remove('active'); }
    function startThree() {
      setTimeout(() => {
        hctx = document.getElementById('km-hist-canvas').getContext('2d');
        initThree(pos);
        lastTime = null;
        running = true;
        if (!frameId) frameId = requestAnimationFrame(loop);
      }, 80);
    }
    if (window.THREE) {
      startThree();
    } else {
      const s = document.createElement('script');
      s.src = 'https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js';
      s.onload = startThree;
      document.head.appendChild(s);
    }
  },

  onClose: function () {
    running = false;
    if (frameId) { cancelAnimationFrame(frameId); frameId = null; }
    const pb = document.getElementById('km-pause-btn');
    if (pb) { pb.textContent = 'Pause'; pb.classList.remove('active'); }
  },

  onResize: function () {
    if (renderer) resizeRenderer();
  },
});

window.kmOpen  = () => shell.open();
window.kmClose = () => shell.close();
window.kmRandomiseFreqs = function () {
  initPhases(); initOmegas();
  smoothBins = new Float64Array(N_BINS);
};
window.kmTogglePause = function () {
  running = !running;
  if (running) lastTime = null;
  const pb = document.getElementById('km-pause-btn');
  if (pb) {
    pb.textContent = running ? 'Pause' : 'Resume';
    pb.classList.toggle('active', !running);
  }
};

document.getElementById('km-K').addEventListener('input', function () {
  const g = GRAPHS[gIdx];
  K = g.Kmin + parseFloat(this.value)*(g.Kmax-g.Kmin);
});
document.getElementById('km-sigma').addEventListener('input', function () {
  sigma = 2*Math.PI*(0.025+parseFloat(this.value)*0.225);
  initOmegas(); smoothBins = new Float64Array(N_BINS);
});

})();
