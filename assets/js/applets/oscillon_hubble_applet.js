(function () {
  'use strict';

  const _cs = getComputedStyle(document.documentElement);
  const _c   = n => _cs.getPropertyValue(n).trim();
  const _rgb  = n => { const h = _c(n).replace('#',''); const v = parseInt(h,16); return [(v>>16)&0xFF,(v>>8)&0xFF,v&0xFF]; };
  const _rgba = (n, a) => { const [r,g,b] = _rgb(n); return `rgba(${r},${g},${b},${a})`; };

  /* ── Grid parameters ── */
  const N        = 2048;
  const L        = 50.0;
  const dx       = L / N;
  const dt       = dx * 0.5;
  const LAMBDA   = 4.0;
  const T0       = 0.4;
  const H_DS     = 1.0;
  let SIGMA0   = 1.5;
  const SUBSTEPS = 4;

  /* ── Simulation state ── */
  let phi     = new Float64Array(N);
  let pi_f    = new Float64Array(N);
  let canvas, ctx, simW, simH;
  let simTime = T0;
  let running = false;
  let frameId = null;
  let frwMode = 'matter';

  /* ── FRW functions ── */
  function hubble(t) {
    if (frwMode === 'matter')    return (2.0 / 3.0) / t;
    if (frwMode === 'radiation') return 0.5 / t;
    return H_DS;
  }

  function scaleFactor(t) {
    if (frwMode === 'matter')    return Math.pow(t, 2.0 / 3.0);
    if (frwMode === 'radiation') return Math.sqrt(t);
    return Math.exp(H_DS * (t - T0));
  }

  let potMode = 'phi4';   // 'phi4' | 'dwell'

  function dV(phi_val) {
    if (potMode === 'dwell')
      return 2.0 * LAMBDA * phi_val * (phi_val - 1.0) * (2.0 * phi_val - 1.0);
    return 4.0 * LAMBDA * phi_val * phi_val * phi_val;
  }

  function syncPotButtons() {
    ['phi4', 'dwell'].forEach(m => {
      const btn = document.getElementById('oh-pot-' + m);
      if (btn) btn.classList.toggle('active', m === potMode);
    });
  }

  window.ohSetPot = function (mode) {
    potMode = mode;
    syncPotButtons();
    init();
  };

  /* ── 8th-order Laplacian with periodic BCs ── */
  function laplacian(arr, i) {
    const i1p = (i + 1 + N) % N, i1m = (i - 1 + N) % N;
    const i2p = (i + 2 + N) % N, i2m = (i - 2 + N) % N;
    const i3p = (i + 3 + N) % N, i3m = (i - 3 + N) % N;
    const i4p = (i + 4 + N) % N, i4m = (i - 4 + N) % N;
    return (
      -9.0     * (arr[i4p] + arr[i4m])
      + 128.0  * (arr[i3p] + arr[i3m])
      - 1008.0 * (arr[i2p] + arr[i2m])
      + 8064.0 * (arr[i1p] + arr[i1m])
      - 14350.0*  arr[i]
    ) / (5040.0 * dx * dx);
  }

  /* ── Leapfrog step ── */
  function step(ddt) {
    const t    = simTime;
    const H    = hubble(t);
    const a    = scaleFactor(t);
    const a2   = a * a;
    const damp1 = 1 - 0.5 * H * ddt;
    const damp2 = 1 + 0.5 * H * ddt;
    for (let i = 0; i < N; i++) {
      const acc = (laplacian(phi, i) / a2 - dV(phi[i])) / damp1;
      pi_f[i] = (pi_f[i] + ddt * acc) * damp1 / damp2;
    }
    for (let i = 0; i < N; i++) phi[i] += ddt * pi_f[i];
    simTime = t + ddt;
  }

  /* ── Energy density ── */
  function energyDensity(i, a) {
    const a2      = a * a;
    const dphi    = (phi[(i + 1) % N] - phi[i]) / dx;
    const phi_val = phi[i];
    const V = potMode === 'dwell'
      ? LAMBDA * phi_val * phi_val * (phi_val - 1.0) * (phi_val - 1.0)
      : LAMBDA * phi_val * phi_val * phi_val * phi_val;
    return 0.5 * pi_f[i] * pi_f[i] + 0.5 * dphi * dphi / a2 + V;
  }

  /* ── Gaussian RNG ── */
  function gauss() {
    let u;
    do { u = Math.random(); } while (u === 0);
    return Math.sqrt(-2.0 * Math.log(u)) * Math.cos(2.0 * Math.PI * Math.random());
  }

  /* ── Init ── */
  function init() {
    simTime = T0;
    for (let i = 0; i < N; i++) { phi[i] = SIGMA0 * gauss(); pi_f[i] = 0.0; }
    const SMOOTH_SIGMA = 2;
    const KR = Math.ceil(3 * SMOOTH_SIGMA);
    const kernel = new Float64Array(2 * KR + 1);
    let ksum = 0.0;
    for (let k = -KR; k <= KR; k++) { kernel[k + KR] = Math.exp(-0.5 * (k / SMOOTH_SIGMA) ** 2); ksum += kernel[k + KR]; }
    for (let k = 0; k < kernel.length; k++) kernel[k] /= ksum;
    const raw = phi.slice();
    for (let i = 0; i < N; i++) {
      let v = 0.0;
      for (let k = -KR; k <= KR; k++) v += kernel[k + KR] * raw[(i + k + N) % N];
      phi[i] = v;
    }
    updateTimeDisplay();
  }

  function updateTimeDisplay() {
    if (!ctx) return;
    const fs   = Math.round(simH * 0.055);
    const text = 't = ' + simTime.toFixed(2);
    ctx.save();
    ctx.font         = `${fs}px 'EB Garamond', Georgia, serif`;
    ctx.fillStyle    = _c('--teal-light');
    ctx.globalAlpha  = 0.85;
    ctx.textAlign    = 'right';
    ctx.textBaseline = 'top';
    ctx.fillText(text, simW - 12, 10);
    ctx.restore();
  }

  /* ── Render ── */
  function render() {
    if (!canvas || !ctx) return;
    ctx.fillStyle = _c('--bg-void');
    ctx.fillRect(0, 0, simW, simH);

    const a      = scaleFactor(simTime);
    const padY   = simH * 0.08;
    const plotH  = simH - 2 * padY;
    const midY   = simH * 0.5;
    const phiScale = plotH * 0.5 / 2.0;

    const rho = new Float32Array(N);
    let rhoTotal = 0.0;
    for (let i = 0; i < N; i++) { rho[i] = energyDensity(i, a); rhoTotal += rho[i]; }
    const rhoScale = (plotH * 0.005) / (rhoTotal / N);

    /* ρ glow */
    ctx.save();
    ctx.shadowColor = _c('--pink-light'); ctx.shadowBlur = 10;
    ctx.strokeStyle = _c('--pink-dark');  ctx.lineWidth  = 1.5; ctx.globalAlpha = 0.5;
    ctx.beginPath();
    for (let i = 0; i < N; i++) {
      const x = (i / N) * simW, y = simH - padY - rho[i] * rhoScale;
      i === 0 ? ctx.moveTo(x, y) : ctx.lineTo(x, y);
    }
    ctx.stroke(); ctx.restore();

    /* ρ main */
    ctx.save();
    ctx.shadowColor = _c('--pink-light'); ctx.shadowBlur = 18;
    ctx.strokeStyle = _c('--pink-light'); ctx.lineWidth  = 1.5;
    ctx.beginPath();
    for (let i = 0; i < N; i++) {
      const x = (i / N) * simW, y = simH - padY - rho[i] * rhoScale;
      i === 0 ? ctx.moveTo(x, y) : ctx.lineTo(x, y);
    }
    ctx.stroke(); ctx.restore();

    /* φ glow */
    ctx.save();
    ctx.shadowColor = _c('--teal-light'); ctx.shadowBlur = 10;
    ctx.strokeStyle = _c('--teal-dark');  ctx.lineWidth  = 1.5; ctx.globalAlpha = 0.5;
    ctx.beginPath();
    for (let i = 0; i < N; i++) {
      const x = (i / N) * simW, y = midY - phi[i] * phiScale;
      i === 0 ? ctx.moveTo(x, y) : ctx.lineTo(x, y);
    }
    ctx.stroke(); ctx.restore();

    /* φ main */
    ctx.save();
    ctx.shadowColor = _c('--teal-light'); ctx.shadowBlur = 18;
    ctx.strokeStyle = _c('--teal-light'); ctx.lineWidth  = 1.5;
    ctx.beginPath();
    for (let i = 0; i < N; i++) {
      const x = (i / N) * simW, y = midY - phi[i] * phiScale;
      i === 0 ? ctx.moveTo(x, y) : ctx.lineTo(x, y);
    }
    ctx.stroke(); ctx.restore();

    /* axis lines at φ=0 and φ=1 */
    ctx.save();
    ctx.strokeStyle = _rgba('--white', 0.3); ctx.lineWidth = 1; ctx.setLineDash([4, 6]);
    [0, 1].forEach(v => {
      const y = midY - v * phiScale;
      ctx.beginPath(); ctx.moveTo(0, y); ctx.lineTo(simW, y); ctx.stroke();
    });
    ctx.restore();

    updateTimeDisplay();
  }

  /* ── Animation loop ── */
  function loop() {
    if (running) { for (let s = 0; s < SUBSTEPS; s++) step(dt); render(); }
    frameId = requestAnimationFrame(loop);
  }

  /* ── FRW mode switching ── */
  function syncModeButtons() {
    ['matter', 'radiation', 'desitter'].forEach(m => {
      const btn = document.getElementById('oh-btn-' + m);
      if (btn) btn.classList.toggle('active', m === frwMode);
    });
  }

  window.ohSetMode = function (mode) {
    frwMode = mode;
    syncModeButtons();
    init();
  };

  /* ── Shell wiring ── */
  const shell = new AppletShell({
    id:     'oh',
    title:  'Oscillon Formation',
    gap:    0,
    layout: 'stacked',

    headerBtns: `<button class="applet-shell-header-btn" onclick="ohReset()">Reset</button><button class="applet-shell-header-btn" id="oh-pause-btn" onclick="ohTogglePause()">Pause</button>`,


    ctrlHTML: `
      <div class="applet-shell-ctrl-section" style="flex:0 0 auto;">
        <div class="applet-shell-btn-row" style="flex-wrap:nowrap;">
          <button class="applet-shell-btn active" id="oh-btn-matter"    onclick="ohSetMode('matter')">Matter</button>
          <button class="applet-shell-btn"         id="oh-btn-radiation" onclick="ohSetMode('radiation')">Radiation</button>
          <button class="applet-shell-btn"         id="oh-btn-desitter"  onclick="ohSetMode('desitter')">de Sitter</button>
        </div>
      </div>
      <div class="applet-shell-ctrl-section" style="flex:0 0 auto;">
        <div class="applet-shell-btn-row" style="flex-wrap:nowrap;">
          <button class="applet-shell-btn active" id="oh-pot-phi4"  onclick="ohSetPot('phi4')">&phi;<sup>4</sup></button>
          <button class="applet-shell-btn"         id="oh-pot-dwell" onclick="ohSetPot('dwell')">&phi;<sup>2</sup>(&phi;&minus;1)<sup>2</sup></button>
        </div>
      </div>
      <div class="applet-shell-ctrl-section" style="flex:1; min-width:140px;">
        <div class="applet-shell-ctrl-title">Fluctuation Amplitude</div>
        <div class="applet-shell-slider-row">
          <span class="applet-shell-side">Small</span>
          <input type="range" id="oh-sigma" min="0.1" max="4.0" step="0.1" value="1.5">
          <span class="applet-shell-side">Large</span>
          <span class="applet-shell-val" id="oh-sigma-val">1.5</span>
        </div>
      </div>
    `,

    onOpen: function ({ canvas: c, W, H, S }) {
      canvas = c;
      ctx    = canvas.getContext('2d');
      simW   = W || S;
      simH   = H || S;

      frwMode = 'matter'; potMode = 'phi4';
      syncModeButtons();
      syncPotButtons();
      const pb = document.getElementById('oh-pause-btn');
      if (pb) { pb.textContent = 'Pause'; pb.classList.remove('active'); }

      SIGMA0 = 1.5;
      const sigSl = document.getElementById('oh-sigma');
      const sigVl = document.getElementById('oh-sigma-val');
      if (sigSl) { sigSl.value = SIGMA0; }
      if (sigVl) { sigVl.textContent = SIGMA0.toFixed(1); }

      init();
      running = true;
      if (!frameId) frameId = requestAnimationFrame(loop);
    },

    onClose: function () {
      running = false;
      if (frameId) { cancelAnimationFrame(frameId); frameId = null; }
    },

    onResize: function ({ W, H, S }) {
      if (!canvas) return;
      simW = canvas.width  = W || S;
      simH = canvas.height = H || S;
      if (!running) render();
    },
  });

  window.ohOpen  = () => shell.open();
  window.ohClose = () => shell.close();

  window.ohReset = function () {
    init();
    if (!running) render();
  };

  document.getElementById('oh-sigma').addEventListener('input', function () {
    SIGMA0 = parseFloat(this.value);
    const vl = document.getElementById('oh-sigma-val');
    if (vl) vl.textContent = SIGMA0.toFixed(1);
    init();
  });

  window.ohTogglePause = function () {
    running = !running;
    const pb = document.getElementById('oh-pause-btn');
    if (pb) {
      pb.textContent = running ? 'Pause' : 'Resume';
      pb.classList.toggle('active', !running);
    }
  };

})();
