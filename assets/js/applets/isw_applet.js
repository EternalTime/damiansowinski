(function () {
  'use strict';

  /* ── Palette ── */
  const _cs = getComputedStyle(document.documentElement);
  const _c   = n => _cs.getPropertyValue(n).trim();
  const _rgb  = n => { const h = _c(n).replace('#',''); const v = parseInt(h,16); return [(v>>16)&0xFF,(v>>8)&0xFF,v&0xFF]; };
  const _rgba = (n, a) => { const [r,g,b] = _rgb(n); return `rgba(${r},${g},${b},${a})`; };

  const TEAL_DARK  = _c('--teal-dark');
  const TEAL_LIGHT = _c('--teal-light');
  const CYAN       = _c('--cyan');
  const PINK_LIGHT = _c('--pink-light');
  const PINK_DARK  = _c('--pink-dark');
  const BG         = _c('--bg-canvas');

  const [_TLR, _TLG, _TLB] = _rgb('--teal-light');
  const [_PLR, _PLG, _PLB] = _rgb('--pink-light');

  /* ── Physics constants (ℏ = m = L = 1) ── */
  const N_STATES = 30;
  const PDF_YMAX = 4.53;
  function energy(n) { return (n * n * Math.PI * Math.PI) / 2; }

  /* ── Simulation state ── */
  let canvas, ctx;
  let S = 500;
  let running  = false;
  let frameId  = null;
  let simTime  = 0;
  let simSpeed = 0.05;

  const amp   = new Float64Array(N_STATES);
  const phase = new Float64Array(N_STATES);

  /* ── Draw-PDF mode state ── */
  let drawMode    = false;
  let drawActive  = false;
  let drawPoints  = [];

  /* ── Gauge drag state ── */
  let gaugeDragging       = -1;
  let gaugeDragStartAngle = 0;
  let gaugeDragStartPhase = 0;

  /* ── Normalisation ── */
  function normFactor() {
    let s = 0;
    for (let i = 0; i < N_STATES; i++) s += amp[i] * amp[i];
    return s > 1e-12 ? 1 / Math.sqrt(s) : 0;
  }

  /* ── ψₙ(x) = √2·sin(nπx) ── */
  function psi_n(n, x) { return Math.SQRT2 * Math.sin(n * Math.PI * x); }

  /* ── Evaluate ψ(x,t) ── */
  function evalPsi(x, t) {
    const nf = normFactor();
    let re = 0, im = 0;
    for (let i = 0; i < N_STATES; i++) {
      if (amp[i] < 1e-9) continue;
      const n     = i + 1;
      const cn    = amp[i] * nf;
      const phi   = phase[i] - energy(n) * t;
      const basis = psi_n(n, x);
      re += cn * Math.cos(phi) * basis;
      im += cn * Math.sin(phi) * basis;
    }
    return { re, im };
  }

  /* ── Init ── */
  function init() {
    simTime = 0;
    for (let i = 0; i < N_STATES; i++) { amp[i] = (i === 0) ? 1.0 : 0.0; phase[i] = 0.0; }
    syncSliders();
    updateGaugeArrows(0);
  }

  function syncSliders() {
    for (let i = 0; i < N_STATES; i++) {
      const sl = document.getElementById('isw-amp-' + i);
      if (sl) sl.value = amp[i];
    }
  }

  /* ── Projection ── */
  const NQ = 500;
  function projectFunction(fVals) {
    for (let i = 0; i < N_STATES; i++) {
      const n = i + 1;
      let s = 0;
      for (let k = 0; k <= NQ; k++) {
        const x  = k / NQ;
        const w  = (k === 0 || k === NQ) ? 0.5 : 1.0;
        s += w * fVals[k] * psi_n(n, x);
      }
      s /= NQ;
      amp[i]   = Math.abs(s);
      phase[i] = s < 0 ? Math.PI : 0;
    }
    const mx = Math.max(...amp);
    if (mx > 1e-12) for (let i = 0; i < N_STATES; i++) amp[i] /= mx;
    syncSliders();
    updateGaugeArrows(simTime);
  }

  /* ── Preset: Dirac Delta ── */
  function presetDelta() {
    exitDrawMode();
    simTime = 0;
    const w  = 0.10;
    const hw = w / 2;
    const x0 = hw + Math.random() * (1 - w);
    const fVals = new Float64Array(NQ + 1);
    for (let k = 0; k <= NQ; k++) {
      const x  = k / NQ;
      const dx = x - x0;
      fVals[k] = Math.abs(dx) < hw ? Math.cos((Math.PI * dx) / w) : 0;
    }
    projectFunction(fVals);
  }

  /* ── Preset: Triangle PDF ── */
  function presetTriangle() {
    exitDrawMode();
    simTime = 0;
    const leftHalf = Math.random() < 0.5;
    const fVals = new Float64Array(NQ + 1);
    for (let k = 0; k <= NQ; k++) {
      const x = k / NQ;
      let pdf = 0;
      if (leftHalf) {
        if      (x <= 0.25) pdf = x / 0.25;
        else if (x <= 0.5)  pdf = (0.5 - x) / 0.25;
      } else {
        if      (x >= 0.5 && x <= 0.75) pdf = (x - 0.5) / 0.25;
        else if (x >= 0.75 && x <= 1.0) pdf = (1.0 - x) / 0.25;
      }
      fVals[k] = Math.sqrt(pdf);
    }
    projectFunction(fVals);
  }

  /* ── Draw PDF mode ── */
  function enterDrawMode() {
    drawMode   = true;
    drawPoints = [];
    if (canvas) canvas.style.cursor = 'crosshair';
    const btn = document.getElementById('isw-draw-btn');
    if (btn) btn.classList.add('active-pink');
  }

  function exitDrawMode() {
    drawMode   = false;
    drawActive = false;
    drawPoints = [];
    if (canvas) canvas.style.cursor = '';
    const btn = document.getElementById('isw-draw-btn');
    if (btn) btn.classList.remove('active-pink');
  }

  function toggleDrawMode() {
    if (drawMode) { exitDrawMode(); return; }
    enterDrawMode();
  }

  function canvasToProbCoords(cx, cy) {
    const L    = getLayout();
    const YMAX = PDF_YMAX;
    const x0 = L.mL, y0 = L.y2top, w = L.plotW, h = L.plotH;
    const xf = (cx - x0) / w;
    const yf = YMAX * (1 - (cy - y0) / h);
    return { xf: Math.max(0, Math.min(1, xf)), yf: Math.max(0, yf) };
  }

  function smoothDrawPoints() {
    if (drawPoints.length < 2) return null;

    const BINS = NQ + 1;
    const yBin = new Float64Array(BINS);
    const cnt  = new Float64Array(BINS);

    for (const p of drawPoints) {
      const k = Math.round(p.xf * NQ);
      if (k < 0 || k > NQ) continue;
      yBin[k] += p.yf;
      cnt[k]  += 1;
    }
    for (let k = 0; k <= NQ; k++) {
      if (cnt[k] > 0) yBin[k] /= cnt[k];
    }

    let lastHit = -1;
    for (let k = 0; k <= NQ; k++) {
      if (cnt[k] > 0) {
        if (lastHit >= 0 && k - lastHit > 1) {
          for (let j = lastHit + 1; j < k; j++) {
            const t = (j - lastHit) / (k - lastHit);
            yBin[j] = yBin[lastHit] * (1 - t) + yBin[k] * t;
          }
        }
        lastHit = k;
      }
    }

    const sigma = 0.025;
    const result = new Float64Array(BINS);
    const sigBins = sigma * NQ;
    const radius  = Math.ceil(3 * sigBins);
    for (let k = 0; k <= NQ; k++) {
      let num = 0, den = 0;
      for (let j = Math.max(0, k - radius); j <= Math.min(NQ, k + radius); j++) {
        const dx = k - j;
        const w  = Math.exp(-(dx * dx) / (2 * sigBins * sigBins));
        num += w * yBin[j];
        den += w;
      }
      result[k] = den > 1e-9 ? num / den : 0;
    }

    result[0] = 0; result[NQ] = 0;
    return result;
  }

  function applyDrawnPDF() {
    const density = smoothDrawPoints();
    if (!density) return;
    const fVals = new Float64Array(NQ + 1);
    for (let k = 0; k <= NQ; k++) fVals[k] = Math.sqrt(Math.max(0, density[k]));
    projectFunction(fVals);
  }

  /* ── Canvas pointer handlers ── */
  function getCanvasCoords(evt) {
    const rect = canvas.getBoundingClientRect();
    const scaleX = canvas.width  / rect.width;
    const scaleY = canvas.height / rect.height;
    const clientX = evt.touches ? evt.touches[0].clientX : evt.clientX;
    const clientY = evt.touches ? evt.touches[0].clientY : evt.clientY;
    return {
      cx: (clientX - rect.left) * scaleX,
      cy: (clientY - rect.top)  * scaleY,
    };
  }

  function isInProbPlot(cx, cy) {
    const L = getLayout();
    return cx >= L.mL && cx <= L.mL + L.plotW &&
           cy >= L.y2top && cy <= L.y2top + L.plotH;
  }

  function onCanvasDown(evt) {
    if (!drawMode) return;
    const { cx, cy } = getCanvasCoords(evt);
    if (!isInProbPlot(cx, cy)) return;
    evt.preventDefault();
    drawActive = true;
    drawPoints = [];
    const { xf, yf } = canvasToProbCoords(cx, cy);
    drawPoints.push({ xf, yf });
  }

  function onCanvasMove(evt) {
    if (!drawMode || !drawActive) return;
    evt.preventDefault();
    const { cx, cy } = getCanvasCoords(evt);
    if (!isInProbPlot(cx, cy)) {
      drawActive = false;
      applyDrawnPDF();
      exitDrawMode();
      return;
    }
    const { xf, yf } = canvasToProbCoords(cx, cy);
    drawPoints.push({ xf, yf });
  }

  function onCanvasUp() {
    if (!drawMode || !drawActive) return;
    drawActive = false;
    applyDrawnPDF();
    exitDrawMode();
  }

  /* ── Layout ── */
  function getLayout() {
    const mL   = Math.round(S * 0.06);
    const mR   = Math.round(S * 0.06);
    const mTop = Math.round(S * 0.04);
    const mBot = Math.round(S * 0.12);
    const mMid = Math.round(S * 0.06);
    const plotW = S - mL - mR;
    const plotH = Math.floor((S - mTop - mBot - mMid) / 2);
    const y1top = mTop;
    const y2top = mTop + plotH + mMid;
    return { mL, mR, mTop, mBot, mMid, plotW, plotH, y1top, y2top };
  }

  /* ── Rendering ── */
  const NX = 300;

  function render() {
    if (!canvas) return;
    ctx.fillStyle = BG;
    ctx.fillRect(0, 0, S, S);

    const L = getLayout();
    drawWell(L, L.y2top + L.plotH);

    if (drawMode) {
      drawDrawModeOverlay(L);
    } else {
      drawWavefunctionPlot(L.mL, L.y1top, L.plotW, L.plotH);
      drawProbabilityPlot (L.mL, L.y2top, L.plotW, L.plotH);
    }
  }

  /* ── Draw-mode overlay ── */
  function drawDrawModeOverlay(L) {
    const x0 = L.mL, y0 = L.y2top, w = L.plotW, h = L.plotH;

    ctx.save();
    ctx.fillStyle = _rgba('--pink-dark', 0.07);
    ctx.fillRect(x0, y0, w, h);
    ctx.restore();

    if (drawPoints.length > 1) {
      const YMAX = 4.53;
      function toY(v) { return y0 + h * (1 - v / YMAX); }
      function toX(f) { return x0 + f * w; }
      ctx.save();
      ctx.strokeStyle = PINK_LIGHT;
      ctx.lineWidth   = 2;
      ctx.globalAlpha = 0.7;
      ctx.beginPath();
      drawPoints.forEach((p, i) => {
        const px = toX(p.xf), py = toY(p.yf);
        if (i === 0) ctx.moveTo(px, py); else ctx.lineTo(px, py);
      });
      ctx.stroke();
      ctx.restore();
    }

    const fs = Math.max(10, Math.round(S * 0.028));
    ctx.save();
    ctx.font = `italic ${fs}px 'EB Garamond', Georgia, serif`;
    ctx.fillStyle = _rgba('--pink-light', 0.55);
    ctx.textAlign = 'center';
    ctx.textBaseline = 'middle';
    ctx.fillText('Draw a probability density', x0 + w / 2, y0 + h / 2);
    ctx.restore();
  }

  /* ── Well + crosshatch ── */
  function drawWell(L, yFloor) {
    const wallW  = Math.max(1, Math.round(S * 0.008));
    const xLeft  = L.mL;
    const xRight = L.mL + L.plotW;
    const yTop   = Math.round(S * 0.01);
    const step   = Math.round(S * 0.04);

    function hatchRect(rx, ry, rw, rh) {
      if (rw <= 0 || rh <= 0) return;
      ctx.save();
      ctx.beginPath();
      ctx.rect(rx, ry, rw, rh);
      ctx.clip();
      ctx.strokeStyle = _rgba('--text-dim', 0.30);
      ctx.lineWidth   = wallW;
      for (let d = -rh; d <= rw + rh; d += step) {
        ctx.beginPath();
        ctx.moveTo(rx + d,      ry);
        ctx.lineTo(rx + d - rh, ry + rh);
        ctx.stroke();
      }
      ctx.restore();
    }

    hatchRect(0, 0, xLeft, S);
    hatchRect(xRight, 0, S - xRight, S);
    hatchRect(xLeft, yFloor, xRight - xLeft, S - yFloor);

    ctx.save();
    ctx.strokeStyle = _rgba('--text-dim', 0.35);
    ctx.lineWidth   = wallW;
    ctx.lineCap     = 'square';
    ctx.beginPath();
    ctx.moveTo(xLeft,  yTop);
    ctx.lineTo(xLeft,  yFloor);
    ctx.lineTo(xRight, yFloor);
    ctx.lineTo(xRight, yTop);
    ctx.stroke();
    ctx.restore();
  }

  /* ── Dashed x-axis ── */
  function drawXAxis(x0, y0, w, h) {
    const zy = y0 + h * 0.5;
    ctx.save();
    ctx.setLineDash([4, 5]);
    ctx.strokeStyle = _rgba('--text-dim', 0.38);
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(x0, zy); ctx.lineTo(x0 + w, zy);
    ctx.stroke();
    ctx.restore();
  }

  /* ── Glow stroke helper ── */
  function strokeGlow(buildPath, glowColor, glowColorDim, mainColor) {
    ctx.save();
    ctx.shadowColor  = glowColor;
    ctx.shadowBlur   = 10;
    ctx.strokeStyle  = glowColorDim;
    ctx.lineWidth    = 2.5;
    ctx.globalAlpha  = 0.5;
    buildPath();
    ctx.stroke();
    ctx.restore();

    ctx.save();
    ctx.shadowColor  = glowColor;
    ctx.shadowBlur   = 18;
    ctx.strokeStyle  = mainColor;
    ctx.lineWidth    = 1.8;
    buildPath();
    ctx.stroke();
    ctx.restore();
  }

  function drawWavefunctionPlot(x0, y0, w, h) {
    const YMAX = 1.84;
    drawXAxis(x0, y0, w, h);

    const t = simTime;
    function toY(v) { return y0 + h * (1 - (v + YMAX) / (2 * YMAX)); }
    function toX(f) { return x0 + f * w; }

    const reArr = new Float32Array(NX + 1);
    const imArr = new Float32Array(NX + 1);
    for (let k = 0; k <= NX; k++) {
      const xf = k / NX;
      const { re, im } = evalPsi(xf, t);
      reArr[k] = re; imArr[k] = im;
    }

    strokeGlow(
      () => { ctx.beginPath(); for (let k = 0; k <= NX; k++) { const px = toX(k/NX), py = toY(imArr[k]); if(k===0) ctx.moveTo(px,py); else ctx.lineTo(px,py); } },
      PINK_LIGHT, PINK_DARK, PINK_LIGHT
    );
    strokeGlow(
      () => { ctx.beginPath(); for (let k = 0; k <= NX; k++) { const px = toX(k/NX), py = toY(reArr[k]); if(k===0) ctx.moveTo(px,py); else ctx.lineTo(px,py); } },
      TEAL_LIGHT, TEAL_DARK, TEAL_LIGHT
    );
  }

  function drawProbabilityPlot(x0, y0, w, h) {
    const YMAX = PDF_YMAX;
    const t = simTime;
    function toY(v) { return y0 + h * (1 - v / YMAX); }
    function toX(f) { return x0 + f * w; }

    const prob = new Float32Array(NX + 1);
    for (let k = 0; k <= NX; k++) {
      const xf = k / NX;
      const { re, im } = evalPsi(xf, t);
      prob[k] = re * re + im * im;
    }

    ctx.save();
    ctx.beginPath();
    ctx.moveTo(toX(0), y0 + h);
    for (let k = 0; k <= NX; k++) ctx.lineTo(toX(k/NX), toY(prob[k]));
    ctx.lineTo(toX(1), y0 + h);
    ctx.closePath();
    const grad = ctx.createLinearGradient(x0, y0, x0, y0 + h);
    grad.addColorStop(0,   _rgba('--cyan', 0.45));
    grad.addColorStop(0.6, _rgba('--cyan', 0.15));
    grad.addColorStop(1,   _rgba('--cyan', 0.03));
    ctx.fillStyle = grad;
    ctx.fill();
    ctx.restore();

    strokeGlow(
      () => { ctx.beginPath(); for (let k = 0; k <= NX; k++) ctx.lineTo(toX(k/NX), toY(prob[k])); },
      CYAN, TEAL_DARK, CYAN
    );
  }

  /* ── Animation loop ── */
  const DT = 1 / 60;

  function loop() {
    if (running) {
      simTime += DT * simSpeed;
      updateGaugeArrows(simTime);
    }
    render();
    frameId = requestAnimationFrame(loop);
  }

  /* ── Gauge arrows ── */
  function gaugeAngle(i, t) { return phase[i] - energy(i + 1) * t; }

  function updateGaugeArrows(t) {
    for (let i = 0; i < N_STATES; i++) {
      const arrow = document.getElementById('isw-gauge-arrow-' + i);
      if (!arrow) continue;
      const cx = 18, cy = 18, r = 12;
      if (amp[i] < 1e-9) {
        arrow.setAttribute('x2', (cx + r).toFixed(2));
        arrow.setAttribute('y2', cy.toFixed(2));
        continue;
      }
      const a  = gaugeAngle(i, t);
      arrow.setAttribute('x2', (cx + r * Math.cos(a)).toFixed(2));
      arrow.setAttribute('y2', (cy - r * Math.sin(a)).toFixed(2));
    }
  }

  /* ── Gauge drag ── */
  function gaugePointerDown(evt, i) {
    if (drawMode) return;
    evt.preventDefault();
    gaugeDragging = i;
    const svg  = document.getElementById('isw-gauge-svg-' + i);
    const rect = svg.getBoundingClientRect();
    const cx   = rect.left + rect.width  / 2;
    const cy   = rect.top  + rect.height / 2;
    const clientX = evt.touches ? evt.touches[0].clientX : evt.clientX;
    const clientY = evt.touches ? evt.touches[0].clientY : evt.clientY;
    gaugeDragStartAngle = Math.atan2(-(clientY - cy), clientX - cx);
    gaugeDragStartPhase = phase[i] + energy(i + 1) * simTime;
  }

  function gaugePointerMove(evt) {
    if (gaugeDragging < 0) return;
    evt.preventDefault();
    const i    = gaugeDragging;
    const svg  = document.getElementById('isw-gauge-svg-' + i);
    const rect = svg.getBoundingClientRect();
    const cx   = rect.left + rect.width  / 2;
    const cy   = rect.top  + rect.height / 2;
    const clientX = evt.touches ? evt.touches[0].clientX : evt.clientX;
    const clientY = evt.touches ? evt.touches[0].clientY : evt.clientY;
    const curAngle   = Math.atan2(-(clientY - cy), clientX - cx);
    const delta      = curAngle - gaugeDragStartAngle;
    phase[i] = (gaugeDragStartPhase + delta) - energy(i + 1) * simTime;
    updateGaugeArrows(simTime);
    if (!running) render();
  }

  function gaugePointerUp() { gaugeDragging = -1; }

  /* ── Build gauge SVG ── */
  function gaugeColor(i) {
    const t = i / Math.max(N_STATES - 1, 1);
    const r = Math.round(_TLR + (_PLR - _TLR) * t);
    const g = Math.round(_TLG + (_PLG - _TLG) * t);
    const b = Math.round(_TLB + (_PLB - _TLB) * t);
    return `rgb(${r},${g},${b})`;
  }

  function makeGaugeHTML(i) {
    const col = gaugeColor(i);
    return `
      <svg id="isw-gauge-svg-${i}" class="isw-gauge-svg" viewBox="0 0 36 36">
        <circle cx="18" cy="18" r="14" fill="none" stroke="${col}" stroke-width="1.5" opacity="0.45"/>
        <line x1="18" y1="4"  x2="18" y2="7"  stroke="${col}" stroke-width="1" opacity="0.4"/>
        <line x1="18" y1="29" x2="18" y2="32" stroke="${col}" stroke-width="1" opacity="0.4"/>
        <line x1="4"  y1="18" x2="7"  y2="18" stroke="${col}" stroke-width="1" opacity="0.4"/>
        <line x1="29" y1="18" x2="32" y2="18" stroke="${col}" stroke-width="1" opacity="0.4"/>
        <line id="isw-gauge-arrow-${i}" x1="18" y1="18" x2="30" y2="18"
              stroke="${col}" stroke-width="2" stroke-linecap="round" opacity="0.9"/>
        <circle cx="18" cy="18" r="2" fill="${col}" opacity="0.8"/>
      </svg>`;
  }

  /* ── Build eigenstate rows ── */
  function buildEigenRows() {
    const subs = ['₁','₂','₃','₄','₅','₆','₇','₈','₉'];
    let html = '';
    for (let i = 0; i < N_STATES; i++) {
      const n   = i + 1;
      const def = (i === 0) ? '1.0' : '0.0';
      const col = gaugeColor(i);
      const label = n <= 9
        ? `c${subs[n-1]}`
        : `c<sub style="font-size:0.7em">${n}</sub>`;
      html += `
        <div class="isw-eigen-row">
          <span class="isw-eigen-label" style="color:${col};">${label}</span>
          <input type="range" id="isw-amp-${i}" min="0" max="1" step="0.01" value="${def}" style="accent-color:${col};">
          <div class="isw-gauge-wrap" id="isw-gauge-wrap-${i}">
            ${makeGaugeHTML(i)}
          </div>
        </div>`;
    }
    return html;
  }

  /* ── Shell ── */
  const shell = new AppletShell({
    id:    'isw',
    title: 'Infinite Square Well &mdash; Quantum Mechanics',
    gap:   0,

    ctrlHTML: `
      <div class="applet-shell-ctrl-section" style="flex-shrink:0;">
        <div class="applet-shell-ctrl-title">Actions</div>
        <div class="applet-shell-btn-row">
          <button class="applet-shell-btn" onclick="iswReset()">Reset</button>
          <button class="applet-shell-btn active" id="isw-pause-btn" onclick="iswTogglePause()">Resume</button>
        </div>
        <div class="applet-shell-slider-row">
          <span class="applet-shell-side">Slow</span>
          <input type="range" id="isw-speed" min="0.002" max="1.0" step="0.002" value="0.05">
          <span class="applet-shell-side">Fast</span>
        </div>
      </div>

      <div class="applet-shell-ctrl-section" style="flex-shrink:0;">
        <div class="applet-shell-ctrl-title">Presets</div>
        <div class="applet-shell-btn-row">
          <button class="applet-shell-btn" onclick="iswPresetDelta()">Dirac Delta</button>
          <button class="applet-shell-btn" onclick="iswPresetTriangle()">Triangle</button>
        </div>
        <div class="applet-shell-btn-row">
          <button class="applet-shell-btn" id="isw-draw-btn" onclick="iswToggleDrawMode()">Draw PDF</button>
        </div>
      </div>

      <div class="applet-shell-ctrl-section" style="flex-shrink:0; padding: 4px 10px 2px;">
        <div class="applet-shell-ctrl-title">Eigenstates &nbsp;<span style="font-weight:normal;font-size:0.78em;color:var(--text-dim);">(amplitude | phase)</span></div>
      </div>

      <div id="isw-scrollable">
        ${buildEigenRows()}
      </div>
    `,

    onOpen: function ({ canvas: c, S: s }) {
      canvas = c;
      ctx    = canvas.getContext('2d');
      S      = s;

      init();

      for (let i = 0; i < N_STATES; i++) {
        const wrap = document.getElementById('isw-gauge-wrap-' + i);
        if (wrap) {
          wrap.addEventListener('mousedown',  (e) => gaugePointerDown(e, i));
          wrap.addEventListener('touchstart', (e) => gaugePointerDown(e, i), { passive: false });
        }
      }
      document.addEventListener('mousemove',  gaugePointerMove);
      document.addEventListener('mouseup',    gaugePointerUp);
      document.addEventListener('touchmove',  gaugePointerMove, { passive: false });
      document.addEventListener('touchend',   gaugePointerUp);

      canvas.addEventListener('mousedown',  onCanvasDown);
      canvas.addEventListener('mousemove',  onCanvasMove);
      canvas.addEventListener('mouseup',    onCanvasUp);
      canvas.addEventListener('touchstart', onCanvasDown, { passive: false });
      canvas.addEventListener('touchmove',  onCanvasMove, { passive: false });
      canvas.addEventListener('touchend',   onCanvasUp);

      for (let i = 0; i < N_STATES; i++) {
        (function(idx) {
          const sl = document.getElementById('isw-amp-' + idx);
          if (sl) sl.addEventListener('input', function () {
            amp[idx] = parseFloat(this.value);
            if (!running) render();
          });
        })(i);
      }

      const speedSl = document.getElementById('isw-speed');
      if (speedSl) {
        speedSl.value = simSpeed;
        speedSl.addEventListener('input', function () { simSpeed = parseFloat(this.value); });
      }

      running = false;
      const pb = document.getElementById('isw-pause-btn');
      if (pb) { pb.textContent = 'Resume'; pb.classList.add('active'); }

      if (!frameId) frameId = requestAnimationFrame(loop);
    },

    onClose: function () {
      running = false;
      exitDrawMode();
      if (frameId) { cancelAnimationFrame(frameId); frameId = null; }
      document.removeEventListener('mousemove',  gaugePointerMove);
      document.removeEventListener('mouseup',    gaugePointerUp);
      document.removeEventListener('touchmove',  gaugePointerMove);
      document.removeEventListener('touchend',   gaugePointerUp);
      const pb = document.getElementById('isw-pause-btn');
      if (pb) { pb.textContent = 'Resume'; pb.classList.add('active'); }
    },

    onResize: function ({ canvas: c, S: s }) {
      canvas = c;
      ctx    = canvas.getContext('2d');
      S      = s;
    },
  });

  /* ── Global entry points ── */
  window.iswOpen  = () => shell.open();
  window.iswClose = () => shell.close();

  window.iswReset = function () {
    simTime = 0;
    exitDrawMode();
    init();
    render();
  };

  window.iswTogglePause = function () {
    running = !running;
    const pb = document.getElementById('isw-pause-btn');
    if (pb) {
      pb.textContent = running ? 'Pause' : 'Resume';
      pb.classList.toggle('active', !running);
    }
  };

  window.iswPresetDelta    = presetDelta;
  window.iswPresetTriangle = presetTriangle;
  window.iswToggleDrawMode = toggleDrawMode;

})();
