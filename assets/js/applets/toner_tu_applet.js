(function () {
  const _cs = getComputedStyle(document.documentElement);
  const _c   = n => _cs.getPropertyValue(n).trim();
  const _rgb  = n => { const h = _c(n).replace('#',''); const v = parseInt(h,16); return [(v>>16)&0xFF,(v>>8)&0xFF,v&0xFF]; };

  /* ── Inject CSS: histogram sections share the leftover panel space ── */
  (function () {
    if (document.getElementById('tt-styles')) return;
    const s = document.createElement('style');
    s.id = 'tt-styles';
    s.textContent = `
      #tt-ctrl-panel { display:flex; flex-direction:column; overflow:hidden; }
      .tt-hist-section { flex:1; min-height:0; display:flex; flex-direction:column; padding:6px 12px 8px; border-top:1px solid var(--border-dark); }
      .tt-hist-section .applet-shell-ctrl-title { flex-shrink:0; margin-bottom:4px; }
      .tt-hist-section canvas { flex:1; min-height:0; display:block; width:100%; }
    `;
    document.head.appendChild(s);
  })();

  const N  = 128;
  const DT = 0.05;
  const DR = 1.0;

  let eta   = 0.2;
  let alpha = 0.3;
  let beta  = 1.0;
  let lam   = 0.5;
  let v0    = 0.5;

  let rho  = new Float32Array(N * N);
  let vx   = new Float32Array(N * N);
  let vy   = new Float32Array(N * N);
  let rho2 = new Float32Array(N * N);
  let vx2  = new Float32Array(N * N);
  let vy2  = new Float32Array(N * N);

  const D_rho = 0.2, D_v = 0.3;

  let _spare = null;
  function randn() {
    if (_spare !== null) { const s = _spare; _spare = null; return s; }
    let u, v, s;
    do { u = Math.random()*2-1; v = Math.random()*2-1; s = u*u+v*v; } while (s >= 1 || s === 0);
    const m = Math.sqrt(-2 * Math.log(s) / s);
    _spare = v * m;
    return u * m;
  }

  function idx(i, j) { return ((i + N) % N) * N + ((j + N) % N); }

  /* Precomputed wrapped neighbor offsets (avoids per-cell modulo) */
  const rowOff = new Int32Array(N), rowUp = new Int32Array(N), rowDn = new Int32Array(N);
  const colE = new Int32Array(N), colW = new Int32Array(N);
  for (let i = 0; i < N; i++) {
    rowOff[i] = i * N;
    rowUp[i]  = ((i - 1 + N) % N) * N;
    rowDn[i]  = ((i + 1) % N) * N;
    colE[i]   = (i + 1) % N;
    colW[i]   = (i - 1 + N) % N;
  }

  function init() {
    for (let k = 0; k < N * N; k++) {
      rho[k] = 1.0 + 0.1 * (Math.random() - 0.5);
      const angle = Math.random() * 2 * Math.PI;
      const speed = v0 * (0.5 + Math.random() * 0.5);
      vx[k] = speed * Math.cos(angle);
      vy[k] = speed * Math.sin(angle);
    }
  }

  function step() {
    const dr2 = DR * DR;
    const noiseMag = eta * Math.sqrt(DT);
    for (let i = 0; i < N; i++) {
      const r0 = rowOff[i], rU = rowUp[i], rD = rowDn[i];
      for (let j = 0; j < N; j++) {
        const jE = colE[j], jW = colW[j];
        const k = r0 + j;
        const r = rho[k], wx = vx[k], wy = vy[k];
        const rE=rho[r0+jE],rW=rho[r0+jW],rN=rho[rU+j],rS=rho[rD+j];
        const wE=vx[r0+jE],wW=vx[r0+jW],wN=vx[rU+j],wS=vx[rD+j];
        const hE=vy[r0+jE],hW=vy[r0+jW],hN=vy[rU+j],hS=vy[rD+j];
        const rNE=rho[rU+jE],rNW=rho[rU+jW],rSE=rho[rD+jE],rSW=rho[rD+jW];
        const wNE=vx[rU+jE],wNW=vx[rU+jW],wSE=vx[rD+jE],wSW=vx[rD+jW];
        const hNE=vy[rU+jE],hNW=vy[rU+jW],hSE=vy[rD+jE],hSW=vy[rD+jW];
        const lap_r  = ((2/3)*(rE+rW+rN+rS)+(1/6)*(rNE+rNW+rSE+rSW)-(10/3)*r) /dr2;
        const lap_wx = ((2/3)*(wE+wW+wN+wS)+(1/6)*(wNE+wNW+wSE+wSW)-(10/3)*wx)/dr2;
        const lap_wy = ((2/3)*(hE+hW+hN+hS)+(1/6)*(hNE+hNW+hSE+hSW)-(10/3)*wy)/dr2;
        const div_v = (wE-wW)/(2*DR) + (hS-hN)/(2*DR);
        const drdx = (rE-rW)/(2*DR), drdy = (rS-rN)/(2*DR);
        const advx = wx*(wE-wW)/(2*DR) + wy*(wS-wN)/(2*DR);
        const advy = wx*(hE-hW)/(2*DR) + wy*(hS-hN)/(2*DR);
        const v2 = wx*wx + wy*wy;
        const coeff = alpha - beta * v2;
        rho2[k] = Math.max(0.01, r + DT*(-(r*div_v + wx*drdx + wy*drdy) + D_rho*lap_r));
        vx2[k]  = wx + DT*(coeff*wx - lam*advx - drdx + D_v*lap_wx) + noiseMag*randn();
        vy2[k]  = wy + DT*(coeff*wy - lam*advy - drdy + D_v*lap_wy) + noiseMag*randn();
      }
    }
    let tmp;
    tmp=rho; rho=rho2; rho2=tmp;
    tmp=vx;  vx=vx2;   vx2=tmp;
    tmp=vy;  vy=vy2;   vy2=tmp;
  }

  const LUT_SIZE = 512;
  const lut = new Uint8Array(LUT_SIZE * 3);
  (function buildLUT() {
    const stops = [
      _rgb('--teal-dark'),
      _rgb('--teal-light'),
      _rgb('--cyan'),
      _rgb('--pink-light'),
      _rgb('--pink-dark'),
    ];
    for (let i = 0; i < LUT_SIZE; i++) {
      const t=i/(LUT_SIZE-1), ft=t*(stops.length-1);
      const lo=Math.floor(ft), hi=Math.min(lo+1,stops.length-1), f=ft-lo;
      lut[i*3]  =Math.round(stops[lo][0]+f*(stops[hi][0]-stops[lo][0]));
      lut[i*3+1]=Math.round(stops[lo][1]+f*(stops[hi][1]-stops[lo][1]));
      lut[i*3+2]=Math.round(stops[lo][2]+f*(stops[hi][2]-stops[lo][2]));
    }
  })();

  function rhoToLUT(r) {
    return Math.round(Math.max(0,Math.min(2,r))/2*(LUT_SIZE-1))*3;
  }

  let canvas, ctx, off, offCtx, imgData, buf;

  function render() {
    let q = 0;
    for (let k = 0; k < N * N; k++, q += 4) {
      const li = rhoToLUT(rho[k]);
      buf[q] = lut[li]; buf[q+1] = lut[li+1]; buf[q+2] = lut[li+2];
    }
    offCtx.putImageData(imgData, 0, 0);
    ctx.drawImage(off, 0, 0, canvas.width, canvas.height);
  }

  let _frameCount = 0;
  let running = false, frameId = null;
  let wasRunning = false;   // sim state stashed while the docs panel is open

  /* \u2500\u2500 |v| and \u03c1 histograms (control panel) \u2500\u2500 */
  let hctx = null, rctx = null;
  const N_BINS = 30;
  let smoothBins    = new Float64Array(N_BINS);
  let smoothRhoBins = new Float64Array(N_BINS);
  const HIST_ALPHA = 0.15;
  const [_TLR, _TLG, _TLB] = _rgb('--teal-light');
  const [_PLR, _PLG, _PLB] = _rgb('--pink-light');
  const [_TDR, _TDG, _TDB] = _rgb('--teal-dark');
  const [_PDR, _PDG, _PDB] = _rgb('--pink-dark');

  function renderHistogram() {
    const hc = document.getElementById('tt-hist-canvas');
    if (!hc || !hctx) return;
    const W = hc.clientWidth || 200, H = hc.clientHeight || 120;
    if (hc.width !== W || hc.height !== H) { hc.width = W; hc.height = H; }
    const vMaxPlot = Math.max(0.2, 2.5 * v0);
    const dv = vMaxPlot / N_BINS;
    const raw = new Float64Array(N_BINS);
    for (let k = 0; k < N * N; k++) {
      const spd = Math.sqrt(vx[k]*vx[k] + vy[k]*vy[k]);
      const b = (spd / dv) | 0;
      if (b < N_BINS) raw[b]++;   // out-of-range values fall in unprinted tail bins
    }
    for (let b = 0; b < N_BINS; b++) raw[b] /= (N * N);
    for (let b = 0; b < N_BINS; b++) smoothBins[b] += HIST_ALPHA * (raw[b] - smoothBins[b]);
    let ymax = 1e-9;
    for (let b = 0; b < N_BINS; b++) if (smoothBins[b] > ymax) ymax = smoothBins[b];
    ymax *= 1.1;
    const PL = 4, PR = 4, PT = 6, PB = 4;
    const pw = W - PL - PR, ph = H - PT - PB, bw = pw / N_BINS;
    hctx.clearRect(0, 0, W, H);
    hctx.save();
    hctx.shadowBlur = 6;
    for (let b = 0; b < N_BINS; b++) {
      const bh = Math.min(smoothBins[b] / ymax, 1) * ph;
      const hot = Math.min((b + 0.5) / N_BINS * vMaxPlot / Math.max(v0, 1e-3) * 0.5, 1);
      const lr = Math.round(_TLR + (_PLR - _TLR) * hot);
      const lg = Math.round(_TLG + (_PLG - _TLG) * hot);
      const lb = Math.round(_TLB + (_PLB - _TLB) * hot);
      const dr = Math.round(_TDR + (_PDR - _TDR) * hot);
      const dg = Math.round(_TDG + (_PDG - _TDG) * hot);
      const db = Math.round(_TDB + (_PDB - _TDB) * hot);
      hctx.shadowColor = `rgb(${dr},${dg},${db})`;
      hctx.fillStyle   = `rgba(${lr},${lg},${lb},0.75)`;
      hctx.fillRect(PL + b * bw, PT + ph - bh, bw - 1, bh);
    }
    hctx.restore();
  }

  /* ρ histogram — bars colored with the same LUT as the density field */
  function renderRhoHistogram() {
    const hc = document.getElementById('tt-rho-canvas');
    if (!hc || !rctx) return;
    const W = hc.clientWidth || 200, H = hc.clientHeight || 120;
    if (hc.width !== W || hc.height !== H) { hc.width = W; hc.height = H; }
    const RHO_MAX = 2.5;
    const dr_ = RHO_MAX / N_BINS;
    const raw = new Float64Array(N_BINS);
    for (let k = 0; k < N * N; k++) {
      const b = (rho[k] / dr_) | 0;
      if (b < N_BINS) raw[b]++;   // out-of-range values fall in unprinted tail bins
    }
    for (let b = 0; b < N_BINS; b++) raw[b] /= (N * N);
    for (let b = 0; b < N_BINS; b++) smoothRhoBins[b] += HIST_ALPHA * (raw[b] - smoothRhoBins[b]);
    let ymax = 1e-9;
    for (let b = 0; b < N_BINS; b++) if (smoothRhoBins[b] > ymax) ymax = smoothRhoBins[b];
    ymax *= 1.1;
    const PL = 4, PR = 4, PT = 6, PB = 4;
    const pw = W - PL - PR, ph = H - PT - PB, bw = pw / N_BINS;
    rctx.clearRect(0, 0, W, H);
    rctx.save();
    rctx.shadowBlur = 6;
    for (let b = 0; b < N_BINS; b++) {
      const bh = Math.min(smoothRhoBins[b] / ymax, 1) * ph;
      const li = rhoToLUT((b + 0.5) * dr_);
      const r = lut[li], g = lut[li+1], bb = lut[li+2];
      rctx.shadowColor = `rgb(${Math.round(r*0.55)},${Math.round(g*0.55)},${Math.round(bb*0.55)})`;
      rctx.fillStyle   = `rgba(${r},${g},${bb},0.8)`;
      rctx.fillRect(PL + b * bw, PT + ph - bh, bw - 1, bh);
    }
    rctx.restore();
  }

  function loop() {
    if (running) {
      for (let s = 0; s < 3; s++) step();
      render();
      if (++_frameCount % 3 === 0) { renderHistogram(); renderRhoHistogram(); }
    }
    frameId = requestAnimationFrame(loop);
  }

  /* ── Shell wiring ── */
  const shell = new AppletShell({
    id:    'tt',
    title: 'Toner&ndash;Tu',
    gap:   0,

    headerBtns: `<button class="applet-shell-header-btn" onclick="ttReset()">Initialise</button><button class="applet-shell-header-btn" id="tt-pause-btn" onclick="ttTogglePause()">Pause</button>`,


    ctrlHTML: `
      <div class="applet-shell-ctrl-section">
        <div class="applet-shell-ctrl-title">Noise &eta;</div>
        <div class="applet-shell-slider-row">
          <span class="applet-shell-side">Low</span>
          <input type="range" id="tt-noise" min="0.0" max="1.0" step="0.01" value="0.2">
          <span class="applet-shell-side">High</span>
        </div>
      </div>
      <div class="applet-shell-ctrl-section">
        <div class="applet-shell-ctrl-title">Activity &alpha;</div>
        <div class="applet-shell-slider-row">
          <span class="applet-shell-side">Disordered</span>
          <input type="range" id="tt-alpha" min="-1.0" max="1.0" step="0.05" value="0.3">
          <span class="applet-shell-side">Ordered</span>
        </div>
      </div>
      <div class="applet-shell-ctrl-section">
        <div class="applet-shell-ctrl-title">Self-advection &lambda;</div>
        <div class="applet-shell-slider-row">
          <span class="applet-shell-side">0</span>
          <input type="range" id="tt-lambda" min="0.0" max="2.0" step="0.05" value="0.5">
          <span class="applet-shell-side">2</span>
        </div>
      </div>
      <div class="applet-shell-ctrl-section">
        <div class="applet-shell-ctrl-title">Speed v&#8320;</div>
        <div class="applet-shell-slider-row">
          <span class="applet-shell-side">Slow</span>
          <input type="range" id="tt-speed" min="0.1" max="2.0" step="0.05" value="0.5">
          <span class="applet-shell-side">Fast</span>
        </div>
      </div>
      <div class="tt-hist-section">
        <div class="applet-shell-ctrl-title">Speed distribution |v|</div>
        <canvas id="tt-hist-canvas"></canvas>
      </div>
      <div class="tt-hist-section">
        <div class="applet-shell-ctrl-title">Density distribution &rho;</div>
        <canvas id="tt-rho-canvas"></canvas>
      </div>
    `,

    onOpen: function ({ canvas: c, S }) {
      canvas = c;
      ctx    = canvas.getContext('2d');
      canvas.width  = S;
      canvas.height = S;
      ctx.imageSmoothingEnabled = true;
      if (!off) {
        off        = document.createElement('canvas');
        off.width  = N;
        off.height = N;
        offCtx     = off.getContext('2d');
        imgData    = offCtx.createImageData(N, N);
        buf        = imgData.data;
        for (let i = 3; i < buf.length; i += 4) buf[i] = 255;
      }
      hctx = document.getElementById('tt-hist-canvas').getContext('2d');
      rctx = document.getElementById('tt-rho-canvas').getContext('2d');
      smoothBins.fill(0);
      smoothRhoBins.fill(0);
      init();
      running = true;
      const pb = document.getElementById('tt-pause-btn');
      if (pb) { pb.textContent = 'Pause'; pb.classList.remove('active'); }
      if (!frameId) loop();
    },

    onClose: function () {
      running = false;
      if (frameId) { cancelAnimationFrame(frameId); frameId = null; }
      const pb = document.getElementById('tt-pause-btn');
      if (pb) { pb.textContent = 'Pause'; pb.classList.remove('active'); }
    },

    onResize: function ({ canvas: c, S }) {
      canvas = c;
      canvas.width  = S;
      canvas.height = S;
      ctx.imageSmoothingEnabled = true;
    },
  });

  window.ttOpen  = () => shell.open();
  window.ttClose = () => shell.close();
  window.ttReset = init;
  window.ttTogglePause = function () {
    running = !running;
    const pb = document.getElementById('tt-pause-btn');
    if (pb) {
      pb.textContent = running ? 'Pause' : 'Resume';
      pb.classList.toggle('active', !running);
    }
  };

  /* Keep steady-state speed √(α/β) pinned to the v₀ slider */
  function updateBeta() {
    if (alpha > 0 && v0 > 0) beta = alpha / (v0 * v0);
  }

  document.getElementById('tt-noise').addEventListener('input', function () {
    eta = parseFloat(this.value);
  });
  document.getElementById('tt-alpha').addEventListener('input', function () {
    alpha = parseFloat(this.value);
    updateBeta();
  });
  document.getElementById('tt-lambda').addEventListener('input', function () {
    lam = parseFloat(this.value);
  });
  document.getElementById('tt-speed').addEventListener('input', function () {
    v0 = parseFloat(this.value);
    updateBeta();
  });
  updateBeta();

})();
