(function () {
  'use strict';

  const _cs = getComputedStyle(document.documentElement);
  const _c   = n => _cs.getPropertyValue(n).trim();
  const _rgb  = n => { const h = _c(n).replace('#',''); const v = parseInt(h,16); return [(v>>16)&0xFF,(v>>8)&0xFF,v&0xFF]; };
  const _rgba = (n, a) => { const [r,g,b] = _rgb(n); return `rgba(${r},${g},${b},${a})`; };

  const [_BVR, _BVG, _BVB] = _rgb('--bg-void');
  const [_TDR, _TDG, _TDB] = _rgb('--teal-dark');
  const [_TLR, _TLG, _TLB] = _rgb('--teal-light');
  const [_PDR, _PDG, _PDB] = _rgb('--pink-dark');
  const [_PLR, _PLG, _PLB] = _rgb('--pink-light');

  /* ── Parameters ── */
  let springK     = 200.0;
  let massM       = 1.0;
  let damping     = 0.5;
  let viewMode    = 'tile';
  let latticeType = 'square';
  let Nbase       = 22;   // interior+2 in x; Ny adjusted for triangular aspect ratio

  let Nx, Ny;             // actual grid dimensions
  let spacing = 1.0;
  let px, py, vx, vy, rx, ry;
  let running = false, frameId = null;
  let wasRunning = false;   // sim state stashed while the docs panel is open
  let canvas, ctx, S;
  let dragIdx = -1, simPanel = null;

  function gridDims() {
    if (latticeType === 'square') {
      return [Nbase, Nbase];
    } else {
      // Row pitch is spacing*√3/2, so to cover the same height need Ny = Nx * 2/√3
      return [Nbase, Math.round(Nbase * 2 / Math.sqrt(3))];
    }
  }

  /* ── Init ── */
  function init(canvasSize) {
    S = canvasSize;
    [Nx, Ny] = gridDims();
    // Interior spans full canvas: spacing = S / (Nx - 1)
    spacing = S / (Nx - 1);

    const n2 = Nx * Ny;
    px = new Float64Array(n2); py = new Float64Array(n2);
    vx = new Float64Array(n2); vy = new Float64Array(n2);
    rx = new Float64Array(n2); ry = new Float64Array(n2);

    const rowPitch = (latticeType === 'triangular') ? spacing * Math.sqrt(3) / 2 : spacing;

    for (let j = 0; j < Ny; j++) {
      for (let i = 0; i < Nx; i++) {
        const idx = j * Nx + i;
        const xOff = (latticeType === 'triangular') ? (j % 2) * spacing * 0.5 : 0;
        rx[idx] = i * spacing + xOff;
        ry[idx] = j * rowPitch;
        px[idx] = rx[idx];
        py[idx] = ry[idx];
        vx[idx] = 0;
        vy[idx] = 0;
      }
    }
    dragIdx = -1;
    buildNeighbors();
  }

  /* ── Precomputed neighbor tables + persistent accel buffers ── */
  let ax, ay, nbrStart, nbrIdx, nbrRest;

  function buildNeighbors() {
    const n2 = Nx * Ny;
    const squareDirs = [[1,0],[-1,0],[0,1],[0,-1]];
    const triEven    = [[1,0],[-1,0],[0,1],[-1,1],[0,-1],[-1,-1]];
    const triOdd     = [[1,0],[-1,0],[0,1],[1,1],[0,-1],[1,-1]];
    nbrStart = new Int32Array(n2 + 1);
    const idxTmp = [], restTmp = [];
    for (let j = 0; j < Ny; j++) {
      for (let i = 0; i < Nx; i++) {
        const a = j * Nx + i;
        const dirs = (latticeType === 'triangular')
          ? (j % 2 === 0 ? triEven : triOdd)
          : squareDirs;
        for (let d = 0; d < dirs.length; d++) {
          const ni = i + dirs[d][0], nj = j + dirs[d][1];
          if (ni < 0 || ni >= Nx || nj < 0 || nj >= Ny) continue;
          const b = nj * Nx + ni;
          const rdx = rx[b] - rx[a], rdy = ry[b] - ry[a];
          idxTmp.push(b);
          restTmp.push(Math.sqrt(rdx*rdx + rdy*rdy));
        }
        nbrStart[a + 1] = idxTmp.length;
      }
    }
    nbrIdx  = Int32Array.from(idxTmp);
    nbrRest = Float64Array.from(restTmp);
    ax = new Float64Array(n2);
    ay = new Float64Array(n2);
  }

  /* ── Physics ── */
  const DT = 0.002;
  const SUBSTEPS = 8;

  function stepOnce(dt) {
    const n2 = Nx * Ny;
    ax.fill(0);
    ay.fill(0);
    const kOverM = springK / massM;

    for (let a = 0; a < n2; a++) {
      const s0 = nbrStart[a], s1 = nbrStart[a + 1];
      const pax = px[a], pay = py[a];
      for (let s = s0; s < s1; s++) {
        const b = nbrIdx[s];
        const cdx = px[b] - pax, cdy = py[b] - pay;
        const curLen = Math.sqrt(cdx*cdx + cdy*cdy);
        if (curLen < 1e-9) continue;
        const f = kOverM * (curLen - nbrRest[s]) / curLen;
        ax[a] += f * cdx;
        ay[a] += f * cdy;
      }
    }

    for (let idx = 0; idx < n2; idx++) {
      const i = idx % Nx, j = (idx / Nx) | 0;
      if (i === 0 || i === Nx-1 || j === 0 || j === Ny-1) continue;
      if (idx === dragIdx) continue;
      ax[idx] -= (damping / massM) * vx[idx];
      ay[idx] -= (damping / massM) * vy[idx];
      vx[idx] += ax[idx] * dt;
      vy[idx] += ay[idx] * dt;
      px[idx] += vx[idx] * dt;
      py[idx] += vy[idx] * dt;
    }
  }

  function step() { for (let s = 0; s < SUBSTEPS; s++) stepOnce(DT); }

  /* ── Color ramp ── */
  function lerpCh(a, b, t) { return Math.round(a + t * (b - a)); }

  /* Neutral midpoint between the two dark hues — used by the network view so
     tension and compression meet in a shared color at zero strain. The tile
     view keeps the original bg-void midpoint. */
  const _MDR = Math.round((_TDR + _PDR) / 2);
  const _MDG = Math.round((_TDG + _PDG) / 2);
  const _MDB = Math.round((_TDB + _PDB) / 2);

  function paletteColor(signed, scale, smoothMid) {
    const mR = smoothMid ? _MDR : _BVR;
    const mG = smoothMid ? _MDG : _BVG;
    const mB = smoothMid ? _MDB : _BVB;
    const t = Math.max(0, Math.min(1, Math.abs(signed) / scale));
    let r, g, b;
    if (signed >= 0) {
      if (t < 0.5) { const s=t/0.5;      r=lerpCh(mR,_TDR,s); g=lerpCh(mG,_TDG,s); b=lerpCh(mB,_TDB,s); }
      else         { const s=(t-0.5)/0.5; r=lerpCh(_TDR,_TLR,s); g=lerpCh(_TDG,_TLG,s); b=lerpCh(_TDB,_TLB,s); }
    } else {
      if (t < 0.5) { const s=t/0.5;      r=lerpCh(mR,_PDR,s); g=lerpCh(mG,_PDG,s); b=lerpCh(mB,_PDB,s); }
      else         { const s=(t-0.5)/0.5; r=lerpCh(_PDR,_PLR,s); g=lerpCh(_PDG,_PLG,s); b=lerpCh(_PDB,_PLB,s); }
    }
    return `rgb(${r},${g},${b})`;
  }

  /* ── Divergence ── */
  function divergence(i, j) {
    if (i < 1 || i >= Nx-1 || j < 1 || j >= Ny-1) return 0;
    const ip=j*Nx+(i+1), im=j*Nx+(i-1), jp=(j+1)*Nx+i, jm=(j-1)*Nx+i;
    const dux = ((px[ip]-rx[ip]) - (px[im]-rx[im])) / (2*spacing);
    const duy = ((py[jp]-ry[jp]) - (py[jm]-ry[jm])) / (2*spacing);
    return dux + duy;
  }

  /* ── Spring extension ── */
  function springExt(a, b) {
    const rdx=rx[b]-rx[a], rdy=ry[b]-ry[a];
    const restLen = Math.sqrt(rdx*rdx + rdy*rdy);
    const cdx=px[b]-px[a], cdy=py[b]-py[a];
    const curLen  = Math.sqrt(cdx*cdx + cdy*cdy);
    return (curLen - restLen) / restLen;
  }

  /* ── Render ── */
  function render() {
    ctx.fillStyle = _c('--bg-void');
    ctx.fillRect(0, 0, S, S);
    if (viewMode === 'tile') renderTile(); else renderNetwork();
    if (dragIdx >= 0) {
      const dotR = Math.max(2, spacing * 0.1);
      ctx.beginPath();
      ctx.arc(px[dragIdx], py[dragIdx], dotR*2.5, 0, Math.PI*2);
      ctx.strokeStyle=_c('--teal-dark'); ctx.lineWidth=1.5; ctx.stroke();
    }
  }

  /* ── Neon mass dots: white-hot core sprite + soft additive halo,
     pre-baked so ~400 dots stay cheap ── */
  let _dotSprite = null, _dotHalo = null, _dotSpriteR = -1;

  function buildDotSprites(dotR) {
    _dotSpriteR = dotR;
    const [wr, wg, wb] = _rgb('--white');
    _dotSprite = document.createElement('canvas');
    const bs = Math.ceil(2 * dotR) + 2;
    _dotSprite.width = _dotSprite.height = bs;
    const bctx = _dotSprite.getContext('2d');
    const bc = bs / 2;
    const bg = bctx.createRadialGradient(bc, bc, 0, bc, bc, dotR);
    bg.addColorStop(0,   `rgba(${wr},${wg},${wb},1)`);
    bg.addColorStop(0.6, `rgba(${wr},${wg},${wb},0.9)`);
    bg.addColorStop(1,   `rgba(${wr},${wg},${wb},0)`);
    bctx.beginPath(); bctx.arc(bc, bc, dotR, 0, Math.PI * 2);
    bctx.fillStyle = bg; bctx.fill();
    const hR = dotR * 3;
    _dotHalo = document.createElement('canvas');
    const hs = Math.ceil(2 * hR) + 2;
    _dotHalo.width = _dotHalo.height = hs;
    const hctx2 = _dotHalo.getContext('2d');
    const hc = hs / 2;
    const hg = hctx2.createRadialGradient(hc, hc, 0, hc, hc, hR);
    hg.addColorStop(0, `rgba(${wr},${wg},${wb},0.35)`);
    hg.addColorStop(1, `rgba(${wr},${wg},${wb},0)`);
    hctx2.beginPath(); hctx2.arc(hc, hc, hR, 0, Math.PI * 2);
    hctx2.fillStyle = hg; hctx2.fill();
  }

  function drawNeonDots(dotR) {
    if (_dotSpriteR !== dotR) buildDotSprites(dotR);
    const hHalf = _dotHalo.width / 2, bHalf = _dotSprite.width / 2;
    ctx.save();
    ctx.globalCompositeOperation = 'lighter';
    for (let j = 1; j < Ny - 1; j++)
      for (let i = 1; i < Nx - 1; i++) {
        const idx = j * Nx + i;
        ctx.drawImage(_dotHalo, px[idx] - hHalf, py[idx] - hHalf);
      }
    ctx.restore();
    for (let j = 1; j < Ny - 1; j++)
      for (let i = 1; i < Nx - 1; i++) {
        const idx = j * Nx + i;
        ctx.drawImage(_dotSprite, px[idx] - bHalf, py[idx] - bHalf);
      }
  }

  function renderTile() {
    let maxDiv = 0.001;
    for (let j=1; j<Ny-1; j++)
      for (let i=1; i<Nx-1; i++) { const d=Math.abs(divergence(i,j)); if(d>maxDiv) maxDiv=d; }
    const scale = Math.max(maxDiv*0.5, 0.02);

    if (latticeType === 'square') {
      const half = spacing * 0.5;
      for (let j=1; j<Ny-1; j++)
        for (let i=1; i<Nx-1; i++) {
          const idx = j*Nx+i;
          ctx.fillStyle = paletteColor(divergence(i,j), scale);
          ctx.fillRect(px[idx]-half, py[idx]-half, spacing, spacing);
        }
    } else {
      // Triangular lattice: Voronoi cells are regular hexagons.
      // Circumradius = spacing / sqrt(3), flat-top orientation (angle offset π/6).
      const R = spacing / Math.sqrt(3);
      for (let j=1; j<Ny-1; j++)
        for (let i=1; i<Nx-1; i++) {
          const idx = j*Nx+i;
          ctx.fillStyle = paletteColor(divergence(i,j), scale);
          ctx.beginPath();
          for (let k=0; k<6; k++) {
            const a = Math.PI/6 + k * Math.PI/3;
            const hx = px[idx] + R * Math.cos(a);
            const hy = py[idx] + R * Math.sin(a);
            k === 0 ? ctx.moveTo(hx, hy) : ctx.lineTo(hx, hy);
          }
          ctx.closePath();
          ctx.fill();
        }
    }

    drawNeonDots(Math.max(1.5, spacing*0.08));
  }

  function renderNetwork() {
    let maxExt = 0.001;
    for (let j=1; j<Ny-1; j++)
      for (let i=1; i<Nx-1; i++) {
        const a=j*Nx+i;
        if (i<Nx-2) maxExt=Math.max(maxExt,Math.abs(springExt(a,j*Nx+(i+1))));
        if (j<Ny-2) {
          if (latticeType==='square') {
            maxExt=Math.max(maxExt,Math.abs(springExt(a,(j+1)*Nx+i)));
          } else {
            const di0=(j%2===0)?0:1, di1=(j%2===0)?-1:0;
            const ni0=i+di0, ni1=i+di1;
            if (ni0>=1&&ni0<=Nx-2) maxExt=Math.max(maxExt,Math.abs(springExt(a,(j+1)*Nx+ni0)));
            if (ni1>=1&&ni1<=Nx-2) maxExt=Math.max(maxExt,Math.abs(springExt(a,(j+1)*Nx+ni1)));
          }
        }
      }
    const scale=Math.max(maxExt*0.5,0.001);

    ctx.lineWidth=1.5;
    for (let j=1; j<Ny-1; j++) {
      for (let i=1; i<Nx-1; i++) {
        const a=j*Nx+i;
        const ax=px[a], ay=py[a];

        if (i<Nx-2) {
          const b=j*Nx+(i+1);
          ctx.strokeStyle=paletteColor(springExt(a,b),scale,true);
          ctx.beginPath(); ctx.moveTo(ax,ay); ctx.lineTo(px[b],py[b]); ctx.stroke();
        }
        if (j<Ny-2) {
          if (latticeType==='square') {
            const b=(j+1)*Nx+i;
            ctx.strokeStyle=paletteColor(springExt(a,b),scale,true);
            ctx.beginPath(); ctx.moveTo(ax,ay); ctx.lineTo(px[b],py[b]); ctx.stroke();
          } else {
            const di0=(j%2===0)?0:1, di1=(j%2===0)?-1:0;
            const ni0=i+di0, ni1=i+di1;
            if (ni0>=1&&ni0<=Nx-2) {
              const b=(j+1)*Nx+ni0;
              ctx.strokeStyle=paletteColor(springExt(a,b),scale,true);
              ctx.beginPath(); ctx.moveTo(ax,ay); ctx.lineTo(px[b],py[b]); ctx.stroke();
            }
            if (ni1>=1&&ni1<=Nx-2) {
              const b=(j+1)*Nx+ni1;
              ctx.strokeStyle=paletteColor(springExt(a,b),scale,true);
              ctx.beginPath(); ctx.moveTo(ax,ay); ctx.lineTo(px[b],py[b]); ctx.stroke();
            }
          }
        }
      }
    }

    drawNeonDots(Math.max(2, spacing*0.1));
  }

  /* ── Loop ── */
  function loop() {
    if (running) { step(); render(); }
    frameId = requestAnimationFrame(loop);
  }

  /* ── Input ── */
  function canvasXY(e) {
    const rect=canvas.getBoundingClientRect();
    const sx=canvas.width/rect.width, sy=canvas.height/rect.height;
    if (e.touches) return [(e.touches[0].clientX-rect.left)*sx,(e.touches[0].clientY-rect.top)*sy];
    return [(e.clientX-rect.left)*sx,(e.clientY-rect.top)*sy];
  }

  function nearestInteriorMass(cx,cy) {
    let best=-1, bestD2=Infinity;
    const thresh=spacing*1.5;
    for (let j=1; j<Ny-1; j++)
      for (let i=1; i<Nx-1; i++) {
        const idx=j*Nx+i, dx=px[idx]-cx, dy=py[idx]-cy, d2=dx*dx+dy*dy;
        if (d2<thresh*thresh && d2<bestD2) { bestD2=d2; best=idx; }
      }
    return best;
  }

  function onDown(e) { const [cx,cy]=canvasXY(e); dragIdx=nearestInteriorMass(cx,cy); e.preventDefault(); }
  function onMove(e) {
    if (dragIdx<0) return;
    const [cx,cy]=canvasXY(e);
    px[dragIdx]=cx; py[dragIdx]=cy; vx[dragIdx]=0; vy[dragIdx]=0;
    e.preventDefault();
  }
  function onUp() { dragIdx=-1; }

  /* ── Shell ── */
  const shell = new AppletShell({
    id: 'sl', title: 'Spring Lattice', gap: 0,

    headerBtns: `<button class="applet-shell-header-btn" onclick="slReset()">Reset</button><button class="applet-shell-header-btn" id="sl-pause-btn" onclick="slTogglePause()">Pause</button>`,

    docs: {
      whatis: `Look closely enough and every solid is a lattice of springs. The atoms of a crystal sit in a regular array, each tethered to its neighbors by interatomic forces that, for small displacements, are as good as Hookean [hooke1678]; sound, heat capacity, and elasticity all live in the vibrations of that network. Navier and Cauchy built the first theories of elasticity in the 1820s from exactly this picture, and Born and von Kármán made lattice vibrations quantitative in 1912, founding what became phonon physics [born1912].¶Each interior mass $m$ here is connected to its neighbors by springs of stiffness $k$ and natural length $a$, with the boundary pinned. Displace a mass and the restoring forces launch waves that cross the lattice at speed
$$c = a \\sqrt{\\frac{k}{m}},$$
reflecting from the pinned edges and interfering on the way back. Damping bleeds the motion away; without it the lattice rings indefinitely.¶The geometry of the network matters. The square lattice is anisotropic: waves and stiffness differ along the axes and the diagonals, and with no diagonal braces it shears far too easily. The triangular lattice, with six bonds per site, is the minimal two-dimensional network that behaves like an isotropic elastic solid at long wavelengths: the same rigidity in every direction, and a well-defined resistance to both compression and shear. Switching between them is a lesson in how continuum elasticity remembers its microscopic scaffolding.`,

      howto: `Grab any interior mass and drag; release and the disturbance propagates. In the Tile view each site's cell is colored by the local divergence of the displacement field — teal for local expansion, pink for compression — so sound waves appear as traveling color fronts. The Network view draws the springs themselves, each colored by its own stretch or compression.¶Lattice switches square versus triangular connectivity, and Resolution sets the grid from $20^2$ up to $80^2$ (both re-initialize the lattice). Spring constant and Mass set the wave speed $c = a\\sqrt{k/m}$: stiffen the springs and watch the fronts quicken. Damping controls how long the ringing persists; set it to zero and drag hard for standing-wave chaos, or high to watch deformations relax quasi-statically.¶Reset returns every mass to its rest position; Pause freezes the lattice mid-wave.`,

      references: ['hooke1678', 'born1912'],
    },

    onDocsOpen:  function () { wasRunning = running; running = false; },
    onDocsClose: function () { running = wasRunning; },

    ctrlHTML: `
      <div id="sl-ctrl-fixed">
        <div class="applet-shell-ctrl-section">
          <div class="applet-shell-ctrl-title">Lattice</div>
          <div class="applet-shell-btn-row">
            <button class="applet-shell-btn active" id="sl-sq-btn"  onclick="slSetLattice('square')">Square</button>
            <button class="applet-shell-btn"         id="sl-tri-btn" onclick="slSetLattice('triangular')">Triangular</button>
          </div>
        </div>
        <div class="applet-shell-ctrl-section">
          <div class="applet-shell-ctrl-title">Resolution</div>
          <div class="applet-shell-btn-row">
            <button class="applet-shell-btn active" id="sl-res-20" onclick="slSetRes(22)">20</button>
            <button class="applet-shell-btn"         id="sl-res-40" onclick="slSetRes(42)">40</button>
            <button class="applet-shell-btn"         id="sl-res-80" onclick="slSetRes(82)">80</button>
          </div>
        </div>
        <div class="applet-shell-ctrl-section">
          <div class="applet-shell-ctrl-title">View</div>
          <div class="applet-shell-btn-row">
            <button class="applet-shell-btn active" id="sl-tile-btn"    onclick="slSetView('tile')">Tile</button>
            <button class="applet-shell-btn"         id="sl-network-btn" onclick="slSetView('network')">Network</button>
          </div>
        </div>
      </div>
      <div id="sl-ctrl-scroll">
        <div class="applet-shell-ctrl-section">
          <div class="applet-shell-ctrl-title">Spring constant</div>
          <div class="applet-shell-slider-row">
            <span class="applet-shell-side">Soft</span>
            <input type="range" id="sl-k" min="10" max="800" step="5" value="200">
            <span class="applet-shell-side">Stiff</span>
            <span class="applet-shell-val" id="sl-k-val">200</span>
          </div>
        </div>
        <div class="applet-shell-ctrl-section">
          <div class="applet-shell-ctrl-title">Mass</div>
          <div class="applet-shell-slider-row">
            <span class="applet-shell-side">Light</span>
            <input type="range" id="sl-m" min="0.1" max="5" step="0.1" value="1.0">
            <span class="applet-shell-side">Heavy</span>
            <span class="applet-shell-val" id="sl-m-val">1.0</span>
          </div>
        </div>
        <div class="applet-shell-ctrl-section">
          <div class="applet-shell-ctrl-title">Damping</div>
          <div class="applet-shell-slider-row">
            <span class="applet-shell-side">None</span>
            <input type="range" id="sl-damp" min="0" max="5" step="0.05" value="0.5">
            <span class="applet-shell-side">High</span>
            <span class="applet-shell-val" id="sl-damp-val">0.50</span>
          </div>
        </div>
        <div class="applet-shell-ctrl-section">
          <div style="font-size:0.78em;color:var(--text-dim);line-height:2.0;">
            <div style="color:var(--teal-dark);">&#9632; Expansion / stretched</div>
            <div style="color:var(--pink-dark);">&#9632; Compression / compressed</div>
          </div>
        </div>
      </div>
    `,

    onOpen: function ({ canvas: c, S: canvasS }) {
      canvas=c; ctx=canvas.getContext('2d');
      init(canvasS);
      running=true;
      const pb=document.getElementById('sl-pause-btn');
      if (pb) { pb.textContent='Pause'; pb.classList.remove('active'); }
      if (!frameId) loop();
      simPanel=document.getElementById('sl-sim-panel')||canvas.parentElement;
      simPanel.addEventListener('mousedown',  onDown);
      simPanel.addEventListener('touchstart', onDown, {passive:false});
      window.addEventListener('mousemove',   onMove);
      window.addEventListener('touchmove',   onMove, {passive:false});
      window.addEventListener('mouseup',     onUp);
      window.addEventListener('touchend',    onUp);
    },

    onClose: function () {
      running=false;
      if (frameId) { cancelAnimationFrame(frameId); frameId=null; }
      const pb=document.getElementById('sl-pause-btn');
      if (pb) { pb.textContent='Pause'; pb.classList.remove('active'); }
      if (simPanel) {
        simPanel.removeEventListener('mousedown',  onDown);
        simPanel.removeEventListener('touchstart', onDown);
      }
      window.removeEventListener('mousemove', onMove);
      window.removeEventListener('touchmove', onMove);
      window.removeEventListener('mouseup',   onUp);
      window.removeEventListener('touchend',  onUp);
    },

    onResize: function ({ canvas: c, S: canvasS }) {
      canvas=c; ctx=canvas.getContext('2d'); init(canvasS);
    },
  });

  /* ── Globals ── */
  window.slOpen  = () => shell.open();
  window.slClose = () => shell.close();
  window.slReset = () => init(S);

  window.slTogglePause = function () {
    running=!running;
    const pb=document.getElementById('sl-pause-btn');
    if (pb) { pb.textContent=running?'Pause':'Resume'; pb.classList.toggle('active',!running); }
  };

  window.slSetRes = function (n) {
    Nbase=n; init(S);
    const map={22:'sl-res-20',42:'sl-res-40',82:'sl-res-80'};
    ['sl-res-20','sl-res-40','sl-res-80'].forEach(id => {
      const b=document.getElementById(id); if(b) b.classList.toggle('active',map[n]===id);
    });
  };

  window.slSetLattice = function (type) {
    latticeType=type; init(S);
    document.getElementById('sl-sq-btn') .classList.toggle('active',type==='square');
    document.getElementById('sl-tri-btn').classList.toggle('active',type==='triangular');
  };

  window.slSetView = function (mode) {
    viewMode=mode;
    document.getElementById('sl-tile-btn')   .classList.toggle('active',mode==='tile');
    document.getElementById('sl-network-btn').classList.toggle('active',mode==='network');
  };

  document.getElementById('sl-k').addEventListener('input', function () {
    springK=parseFloat(this.value);
    document.getElementById('sl-k-val').textContent=springK.toFixed(0);
  });
  document.getElementById('sl-m').addEventListener('input', function () {
    massM=parseFloat(this.value);
    document.getElementById('sl-m-val').textContent=massM.toFixed(1);
  });
  document.getElementById('sl-damp').addEventListener('input', function () {
    damping=parseFloat(this.value);
    document.getElementById('sl-damp-val').textContent=damping.toFixed(2);
  });

})();
