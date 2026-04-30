(function () {
  'use strict';

  /* ── Palette ── */
  const TEAL_DARK  = '#23bbad';
  const TEAL_LIGHT = '#25d9c8';
  const CYAN       = '#2abed9';
  const PINK_LIGHT = '#ff6da2';
  const PINK_DARK  = '#f92672';
  const BG         = '#0d1117';

  /* ── Physics (ℏ=1, E=k²/2M) ── */
  const MASS      = 3.0;
  const K0        = 24.0;
  const SIGMA_K   = 2.5;

  /* ── Grid ── */
  const NX        = 2048;          /* must be power of 2 for FFT */
  const XMIN      = -18.0;
  const XMAX      =  22.0;
  const XSPAN     = XMAX - XMIN;
  const DX        = XSPAN / NX;    /* periodic grid: NX points, spacing DX */

  /* ── Time stepping ── */
  const DT        = 0.002;         /* dt: kinetic phase k²dt/2M < π/4 for k_max~33 */

  /* ── Absorbing boundary ── */
  const ABS_FRAC  = 0.20;          /* absorbing layer = 20% of domain on each side */
  const ABS_STR   = 20.0;          /* absorbing potential strength */

  /* ── Display window (interior region, excluding absorbing layers) ── */
  const XDISP_MIN = -10.0;
  const XDISP_MAX =  14.0;
  const XDISP_SPAN = XDISP_MAX - XDISP_MIN;

  /* ── Rendering ── */
  const BARRIER_MAX = 180.0;
  const PDF_SCALE   = 0.75;
  const YMAX_PSI    = 5.0 / 4.0;
  const YMAX_PROB   = 28.8 / PDF_SCALE / 48.0;
  const J_SCALE     = 0.01;        /* scale for flux diagnostic */

  /* ── State ── */
  let canvas, ctx, simW, simH;
  let frameId   = null;
  let firing    = false;
  let canFire   = true;
  let paused    = false;
  let probArmed = false;
  let probPeak  = 0;
  let speedVal = 24;    /* slider value 1–40; <5 → fractional steps, ≥5 → integer steps */
  let frameSkipCount = 0;
  let barrierH  = 60.0;
  let barrierW  = 0.8;

  /* wavefunction arrays (real + imag, Float64 for accuracy) */
  let psiRe, psiIm;

  /* flux integrals for R and T */
  let p0 = 1.0;    /* initial norm */
  let phiR = 0.0;  /* accumulated leftward flux through x=0 */
  let phiT = 0.0;  /* accumulated rightward flux through x=barrierW */

  /* precomputed arrays */
  let xArr;          /* physical x at each grid point */
  let vHalf;         /* half-step potential phase: exp(-i V(x) dt/2) */
  let vHalfRe, vHalfIm;
  let kPhase;        /* full-step kinetic phase: exp(-i k² dt / 2M) */
  let kPhaseRe, kPhaseIm;
  let absorb;        /* absorbing envelope, applied each half V step */

  /* ── FFT (Cooley-Tukey, in-place, radix-2) ── */
  function fft(re, im, inverse) {
    const n = re.length;
    /* bit-reversal permutation */
    for (let i = 1, j = 0; i < n; i++) {
      let bit = n >> 1;
      for (; j & bit; bit >>= 1) j ^= bit;
      j ^= bit;
      if (i < j) {
        let t = re[i]; re[i] = re[j]; re[j] = t;
            t = im[i]; im[i] = im[j]; im[j] = t;
      }
    }
    /* butterfly */
    for (let len = 2; len <= n; len <<= 1) {
      const ang = (inverse ? 1 : -1) * 2 * Math.PI / len;
      const wRe = Math.cos(ang), wIm = Math.sin(ang);
      for (let i = 0; i < n; i += len) {
        let uRe = 1, uIm = 0;
        for (let j = 0; j < len/2; j++) {
          const eRe = re[i+j+len/2]*uRe - im[i+j+len/2]*uIm;
          const eIm = re[i+j+len/2]*uIm + im[i+j+len/2]*uRe;
          re[i+j+len/2] = re[i+j] - eRe;
          im[i+j+len/2] = im[i+j] - eIm;
          re[i+j] += eRe;
          im[i+j] += eIm;
          const nuRe = uRe*wRe - uIm*wIm;
          uIm = uRe*wIm + uIm*wRe; uRe = nuRe;
        }
      }
    }
    if (inverse) { for (let i = 0; i < n; i++) { re[i] /= n; im[i] /= n; } }
  }

  /* ── Build V(x): rectangular barrier centred at x=0..barrierW ── */
  function buildV() {
    const v = new Float64Array(NX);
    for (let i = 0; i < NX; i++) {
      const x = xArr[i];
      if (x >= 0 && x <= barrierW) v[i] = barrierH;
    }
    return v;
  }

  /* ── Precompute propagator phases and absorbing layer ── */
  function precompute() {
    xArr    = new Float64Array(NX);
    vHalfRe = new Float64Array(NX);
    vHalfIm = new Float64Array(NX);
    kPhaseRe = new Float64Array(NX);
    kPhaseIm = new Float64Array(NX);
    absorb   = new Float64Array(NX);

    for (let i = 0; i < NX; i++) xArr[i] = XMIN + i * DX;

    /* potential half-step phase: e^{-i V(x) dt/2} */
    const v = buildV();
    for (let i = 0; i < NX; i++) {
      const ph = -v[i] * DT / 2;
      vHalfRe[i] = Math.cos(ph);
      vHalfIm[i] = Math.sin(ph);
    }

    /* kinetic full-step phase in k-space: e^{-i k² dt / 2M}
       k values for FFT: [0,1,...,N/2-1,-N/2,...,-1] * (2π/L) */
    const dk = 2 * Math.PI / XSPAN;
    for (let i = 0; i < NX; i++) {
      const ki = i <= NX/2 ? i : i - NX;
      const k  = ki * dk;
      const ph = -(k*k / (2*MASS)) * DT;
      kPhaseRe[i] = Math.cos(ph);
      kPhaseIm[i] = Math.sin(ph);
    }

    /* absorbing envelope: cos² ramp in boundary layers, 1 in interior
       applied as a multiplicative damping each half-V step */
    const absLen = Math.floor(ABS_FRAC * NX);
    for (let i = 0; i < NX; i++) {
      if (i < absLen) {
        const t = i / absLen;
        absorb[i] = Math.sin(0.5 * Math.PI * t);  /* 0→1 */
        absorb[i] *= absorb[i];
      } else if (i >= NX - absLen) {
        const t = (NX - 1 - i) / absLen;
        absorb[i] = Math.sin(0.5 * Math.PI * t);
        absorb[i] *= absorb[i];
      } else {
        absorb[i] = 1.0;
      }
    }
  }

  /* ── Initialise wavepacket directly in x-space ──
     ψ(x,0) = A · exp(-(x-x0)²/(4σ²)) · e^{iK0x}
     σ = 1/(2·SIGMA_K) so the k-space width matches SIGMA_K */
  function initPsi() {
    psiRe = new Float64Array(NX);
    psiIm = new Float64Array(NX);

    const x0    = -5.0;  /* fixed initial packet centre, independent of domain */
    const sigma = 1.0 / SIGMA_K;  /* spatial width */
    const sig2  = 2.0 * sigma * sigma;

    for (let i = 0; i < NX; i++) {
      const x   = xArr[i];
      const env = Math.exp(-(x-x0)*(x-x0) / sig2);
      psiRe[i]  = env * Math.cos(K0 * x) * absorb[i];
      psiIm[i]  = env * Math.sin(K0 * x) * absorb[i];
    }
  }

  /* ── Compute norm ── */
  function computeNorm() {
    let n = 0;
    for (let i = 0; i < NX; i++) n += (psiRe[i]*psiRe[i] + psiIm[i]*psiIm[i]) * DX;
    return n;
  }

  /* ── Flux at a single grid point ── */
  function fluxAt(i) {
    const i0 = Math.max(i-1,0), i1 = Math.min(i+1,NX-1);
    const dRe = (psiRe[i1]-psiRe[i0]) / ((i1-i0)*DX);
    const dIm = (psiIm[i1]-psiIm[i0]) / ((i1-i0)*DX);
    return (psiRe[i]*dIm - psiIm[i]*dRe) / MASS;
  }

  /* ── Single split-operator step ──
     ψ → e^{-iV dt/2} · IFFT[ e^{-ik²dt/2M} · FFT[ e^{-iV dt/2} · ψ ] ] */
  function step() {
    /* half V */
    for (let i = 0; i < NX; i++) {
      const re = psiRe[i], im = psiIm[i];
      psiRe[i] = re*vHalfRe[i] - im*vHalfIm[i];
      psiIm[i] = re*vHalfIm[i] + im*vHalfRe[i];
      /* absorb */
      psiRe[i] *= absorb[i];
      psiIm[i] *= absorb[i];
    }

    /* FFT */
    fft(psiRe, psiIm, false);

    /* kinetic */
    for (let i = 0; i < NX; i++) {
      const re = psiRe[i], im = psiIm[i];
      psiRe[i] = re*kPhaseRe[i] - im*kPhaseIm[i];
      psiIm[i] = re*kPhaseIm[i] + im*kPhaseRe[i];
    }

    /* IFFT */
    fft(psiRe, psiIm, true);

    /* half V */
    for (let i = 0; i < NX; i++) {
      const re = psiRe[i], im = psiIm[i];
      psiRe[i] = re*vHalfRe[i] - im*vHalfIm[i];
      psiIm[i] = re*vHalfIm[i] + im*vHalfRe[i];
      /* absorb */
      psiRe[i] *= absorb[i];
      psiIm[i] *= absorb[i];
    }
  }

  /* ── Render ── */
  function render() {
    if (!canvas) return;
    ctx.fillStyle = BG;
    ctx.fillRect(0, 0, simW, simH);

    const plotH  = Math.floor(simH * 0.38);
    const groundY = simH - Math.floor(simH * 0.18);
    const y1top  = Math.floor(simH * 0.04);

    /* barrier pixel coords — mapped through display window */
    const xBarL = Math.round(simW * (-XDISP_MIN / XDISP_SPAN));
    const xBarR = Math.round(simW * ((-XDISP_MIN + barrierW) / XDISP_SPAN));
    const yBarTop    = groundY - plotH * (barrierH / BARRIER_MAX);
    const barrierPixH = groundY - yBarTop;

    /* prob and psi arrays from current psiRe/psiIm */
    const prob   = new Float32Array(NX);
    const flux   = new Float32Array(NX);
    for (let i = 0; i < NX; i++) prob[i] = psiRe[i]*psiRe[i] + psiIm[i]*psiIm[i];
    /* flux: central differences */
    for (let i = 0; i < NX; i++) {
      const i0 = Math.max(i-1,0), i1 = Math.min(i+1,NX-1);
      const dRe = (psiRe[i1]-psiRe[i0]) / ((i1-i0)*DX);
      const dIm = (psiIm[i1]-psiIm[i0]) / ((i1-i0)*DX);
      flux[i] = (psiRe[i]*dIm - psiIm[i]*dRe) / MASS;
    }

    /* map helpers — x mapped through display window [XDISP_MIN, XDISP_MAX] */
    function toX(i)  { return ((xArr[i] - XDISP_MIN) / XDISP_SPAN) * simW; }
    function toY1(v) { return y1top + plotH * (1 - (v + YMAX_PSI) / (2*YMAX_PSI)); }
    function toY2(v) { return groundY - plotH * (v / YMAX_PROB); }
    /* index bounds for display window */
    const iDisp0 = Math.ceil((XDISP_MIN - XMIN) / DX);
    const iDisp1 = Math.floor((XDISP_MAX - XMIN) / DX);

    /* |ψ|² fill */
    ctx.save();
    ctx.beginPath();
    ctx.moveTo(toX(iDisp0), groundY);
    for (let i = iDisp0; i <= iDisp1; i++) ctx.lineTo(toX(i), toY2(prob[i]));
    ctx.lineTo(toX(iDisp1), groundY);
    ctx.closePath();
    const grad = ctx.createLinearGradient(0, groundY-plotH, 0, groundY);
    grad.addColorStop(0,   'rgba(42,190,217,0.85)');
    grad.addColorStop(0.6, 'rgba(42,190,217,0.55)');
    grad.addColorStop(1,   'rgba(42,190,217,0.20)');
    ctx.fillStyle = grad;
    ctx.fill();
    ctx.restore();

    /* |ψ|² outline */
    ctx.save();
    ctx.strokeStyle = CYAN; ctx.lineWidth = 1.5;
    ctx.shadowColor = CYAN; ctx.shadowBlur = 10;
    ctx.beginPath();
    for (let i = iDisp0; i <= iDisp1; i++) { const px=toX(i),py=toY2(prob[i]); i===iDisp0?ctx.moveTo(px,py):ctx.lineTo(px,py); }
    ctx.stroke();
    ctx.restore();

    /* incident energy dashed line */
    const E0  = K0*K0 / (2*MASS);
    const yE0 = groundY - plotH * (E0 / BARRIER_MAX);
    if (yE0 >= groundY - plotH && yE0 <= groundY) {
      ctx.save();
      ctx.setLineDash([4,6]);
      ctx.strokeStyle = 'rgba(255,109,162,0.55)';
      ctx.lineWidth = 2.5;
      ctx.beginPath(); ctx.moveTo(0,yE0); ctx.lineTo(simW,yE0);
      ctx.stroke();
      ctx.restore();
    }

    /* crosshatch */
    const hatchColor = 'rgba(100,130,150,0.30)';
    const hatchLW    = Math.max(1, Math.round(simW*0.003));
    const step       = Math.round(simW*0.035);
    function hatch(rx,ry,rw,rh) {
      if (rw<=0||rh<=0) return;
      ctx.save();
      ctx.beginPath(); ctx.rect(rx,ry,rw,rh); ctx.clip();
      ctx.strokeStyle = hatchColor; ctx.lineWidth = hatchLW;
      for (let d=-simH; d<=simW+simH; d+=step) {
        ctx.beginPath(); ctx.moveTo(d,groundY); ctx.lineTo(d-simH,groundY+simH); ctx.stroke();
        ctx.beginPath(); ctx.moveTo(d,groundY); ctx.lineTo(d+simH,groundY-simH); ctx.stroke();
      }
      ctx.restore();
    }
    ctx.save();
    ctx.globalAlpha = 0.7;
    hatch(0, groundY, simW, simH-groundY);
    hatch(xBarL, yBarTop, xBarR-xBarL, barrierPixH);
    const lineColor = 'rgba(100,130,150,0.65)';
    ctx.strokeStyle = lineColor; ctx.lineWidth = hatchLW*1.5;
    ctx.beginPath();
    ctx.moveTo(xBarL,groundY); ctx.lineTo(xBarL,yBarTop);
    ctx.lineTo(xBarR,yBarTop); ctx.lineTo(xBarR,groundY);
    ctx.stroke();
    ctx.restore();

    /* zero line for ψ plot */
    const zyFlux = y1top + plotH/2;
    ctx.save();
    ctx.setLineDash([4,5]);
    ctx.strokeStyle = 'rgba(168,192,208,0.28)';
    ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(0,zyFlux); ctx.lineTo(simW,zyFlux); ctx.stroke();
    ctx.restore();

    /* ±|ψ| envelope — thin cyan, no glow */
    ctx.save();
    ctx.shadowBlur = 0;
    ctx.shadowColor = 'transparent';
    ctx.strokeStyle = 'rgba(42,190,217,0.5)';
    ctx.lineWidth = 2.5;
    ctx.beginPath();
    for (let i=iDisp0;i<=iDisp1;i++) {
      const px=toX(i), py=toY1(Math.sqrt(prob[i]));
      i===iDisp0?ctx.moveTo(px,py):ctx.lineTo(px,py);
    }
    ctx.stroke();
    ctx.beginPath();
    for (let i=iDisp0;i<=iDisp1;i++) {
      const px=toX(i), py=toY1(-Math.sqrt(prob[i]));
      i===iDisp0?ctx.moveTo(px,py):ctx.lineTo(px,py);
    }
    ctx.stroke();
    ctx.restore();

    /* Im(ψ) — pink */
    function strokeGlow(buildPath, glowColor, dimColor, mainColor) {
      ctx.save();
      ctx.shadowColor=glowColor; ctx.shadowBlur=10;
      ctx.strokeStyle=dimColor;  ctx.lineWidth=2.5; ctx.globalAlpha=0.5;
      buildPath(); ctx.stroke(); ctx.restore();
      ctx.save();
      ctx.shadowColor=glowColor; ctx.shadowBlur=18;
      ctx.strokeStyle=mainColor; ctx.lineWidth=1.8;
      buildPath(); ctx.stroke(); ctx.restore();
    }
    strokeGlow(
      ()=>{ ctx.beginPath(); for(let i=iDisp0;i<=iDisp1;i++){const px=toX(i),py=toY1(psiIm[i]); i===iDisp0?ctx.moveTo(px,py):ctx.lineTo(px,py);} },
      PINK_LIGHT, PINK_DARK, PINK_LIGHT
    );
    /* Re(ψ) — teal */
    strokeGlow(
      ()=>{ ctx.beginPath(); for(let i=iDisp0;i<=iDisp1;i++){const px=toX(i),py=toY1(psiRe[i]); i===iDisp0?ctx.moveTo(px,py):ctx.lineTo(px,py);} },
      TEAL_LIGHT, TEAL_DARK, TEAL_LIGHT
    );


    /* probability diagnostics — integrate over display window only */
    let pL=0, pB=0, pR=0;
    for (let i=iDisp0;i<=iDisp1;i++) {
      const x = xArr[i];
      const w = (i===iDisp0||i===iDisp1) ? 0.5 : 1.0;
      const dP = prob[i]*w*DX;
      if      (x < 0)          pL += dP;
      else if (x <= barrierW)  pB += dP;
      else                     pR += dP;
    }
    const pTot = pL+pB+pR;
    const fs = Math.max(10, Math.round(simH*0.022));
    /* ── Diagnostics (commented out) ──
    ctx.save();
    ctx.font = `${fs}px monospace`;
    ctx.fillStyle = 'rgba(255,255,255,0.85)';
    ctx.textBaseline = 'bottom';
    const yLabel = groundY - 4;
    ctx.textAlign = 'center';
    ctx.fillText(pL.toFixed(4), xBarL/2, yLabel);
    ctx.fillText(pB.toFixed(4), (xBarL+xBarR)/2, yLabel);
    ctx.fillText(pR.toFixed(4), (xBarR+simW)/2, yLabel);
    ctx.textAlign = 'left';
    ctx.fillText('tot='+pTot.toFixed(4), 4, yLabel);
    ctx.fillText('p0='+p0.toFixed(4), 4, yLabel - fs - 2);
    const R = phiR / p0;
    const T = phiT / p0;
    const yRT = groundY - plotH * 0.85;
    ctx.textAlign = 'center';
    ctx.fillText('R=' + R.toFixed(3), simW / 6,       yRT);
    ctx.fillText('T=' + T.toFixed(3), simW * 5 / 6,   yRT);
    ctx.restore();
    ── end diagnostics ── */

    return pTot;
  }

  /* ── Fire button / slider lock ── */
  function setFireReady(ready) {
    canFire = ready;
    const btn = document.getElementById('qtun-fire-btn');
    if (!btn) return;
    btn.classList.toggle('qtun-active-pink',  ready);
    btn.classList.toggle('qtun-active-teal', !ready);
    const hSl = document.getElementById('qtun-height');
    const wSl = document.getElementById('qtun-width');
    if (hSl) hSl.disabled = !ready;
    if (wSl) wSl.disabled = !ready;
  }

  /* ── Grid indices for flux measurement ── */
  function iNearest(x) { return Math.round((x - XMIN) / DX); }

  /* ── Animation loop ── */
  function loop() {
    if (firing && !paused) {
      /* accumulate flux through measurement points each substep */
      const iL = iNearest(-9.0);
      const iR = iNearest(13.0);
      /* speedVal 1–4: step every (5-speedVal) frames; speedVal 5–40: 1–10 steps/frame */
      const nSteps = speedVal >= 5 ? Math.max(1, Math.round((speedVal - 4) / 4)) : 0;
      const framesPerStep = speedVal < 5 ? (5 - speedVal) : 1;
      frameSkipCount++;
      const doSteps = speedVal >= 5 || frameSkipCount >= framesPerStep;
      if (doSteps) {
        if (speedVal < 5) frameSkipCount = 0;
        const n = speedVal >= 5 ? nSteps : 1;
        for (let s=0; s<n; s++) {
          step();
          phiR -= fluxAt(iL) * DT;
          phiT += fluxAt(iR) * DT;
        }
      }
      const pTot = render();

      if (!probArmed && pTot > 0.5) { probArmed = true; probPeak = pTot; }
      if (probArmed && pTot > probPeak) probPeak = pTot;
      if (probArmed && pTot < 0.10 * probPeak) {
        firing = false; probArmed = false; probPeak = 0;
        setFireReady(true);
        const pb = document.getElementById('qtun-pause-btn');
        if (pb) pb.textContent = 'Pause';
      }
    }
    frameId = requestAnimationFrame(loop);
  }

  /* ── Layout ── */
  function layout() {
    const PAD = 20;
    const PHI = 1.6180339887;
    const vw  = window.innerWidth;
    const vh  = window.innerHeight;
    const W   = Math.min(vw - 2*PAD, 960);
    const Hsim  = Math.round(W / PHI);
    const Hhdr  = 48;
    const Hctrl = 90;
    const left  = Math.round((vw - W) / 2);
    const topHdr  = Math.round((vh - Hhdr - Hsim - Hctrl) / 2);
    const topBody = topHdr + Hhdr;
    const r = document.documentElement;
    r.style.setProperty('--qtun-W',       W+'px');
    r.style.setProperty('--qtun-left',    left+'px');
    r.style.setProperty('--qtun-top-hdr', topHdr+'px');
    r.style.setProperty('--qtun-top-body',topBody+'px');
    r.style.setProperty('--qtun-H-hdr',   Hhdr+'px');
    r.style.setProperty('--qtun-H-sim',   Hsim+'px');
    r.style.setProperty('--qtun-H-ctrl',  Hctrl+'px');

    if (canvas) {
      simW = canvas.width  = W;
      simH = canvas.height = Hsim;
      render();
    }
  }

  /* ── Open / Close ── */
  window.qtunOpen = function () {
    canvas = document.getElementById('qtun-canvas');
    ctx    = canvas.getContext('2d');

    /* set canvas dimensions before precompute/initPsi/render */
    const PAD = 20, PHI = 1.6180339887;
    const vw = window.innerWidth, vh = window.innerHeight;
    const W  = Math.min(vw - 2*PAD, 960);
    simW = canvas.width  = W;
    simH = canvas.height = Math.round(W / PHI);

    precompute();
    initPsi();
    p0 = computeNorm();
    phiR = 0; phiT = 0;
    layout();   /* sets CSS vars and calls render() with valid canvas */

    /* slider wiring */
    const hSl = document.getElementById('qtun-height');
    if (hSl) hSl.max = BARRIER_MAX;
    if (hSl) hSl.oninput = function () {
      barrierH = parseFloat(this.value);
      document.getElementById('qtun-height-val').textContent = barrierH.toFixed(1);
      if (!firing) { precompute(); initPsi(); p0 = computeNorm(); phiR = 0; phiT = 0; render(); }
    };
    const wSl = document.getElementById('qtun-width');
    if (wSl) wSl.oninput = function () {
      barrierW = parseFloat(this.value);
      document.getElementById('qtun-width-val').textContent = barrierW.toFixed(2);
      if (!firing) { precompute(); initPsi(); p0 = computeNorm(); phiR = 0; phiT = 0; render(); }
    };
    const spSl = document.getElementById('qtun-speed');
    if (spSl) spSl.oninput = function () {
      speedVal = parseInt(this.value);
      frameSkipCount = 0;
    };

    document.getElementById('qtun-overlay').classList.add('qtun-open');
    requestAnimationFrame(() => {
      document.getElementById('qtun-header').classList.add('qtun-open');
      document.getElementById('qtun-sim-panel').classList.add('qtun-open');
      document.getElementById('qtun-ctrl-panel').classList.add('qtun-open');
    });

    if (!frameId) frameId = requestAnimationFrame(loop);
  };

  window.qtunClose = function () {
    ['qtun-header','qtun-sim-panel','qtun-ctrl-panel'].forEach(id => {
      const el = document.getElementById(id);
      el.classList.remove('qtun-open');
      el.classList.add('qtun-closing');
    });
    setTimeout(() => {
      document.getElementById('qtun-overlay').classList.remove('qtun-open');
      ['qtun-header','qtun-sim-panel','qtun-ctrl-panel'].forEach(id =>
        document.getElementById(id).classList.remove('qtun-closing'));
    }, 600);
  };

  window.qtunFire = function () {
    if (!canFire) return;
    precompute();
    initPsi();
    p0 = computeNorm();
    phiR = 0; phiT = 0;
    probArmed = false; probPeak = 0; paused = false;
    firing = true;
    setFireReady(false);
    const pb = document.getElementById('qtun-pause-btn');
    if (pb) pb.textContent = 'Pause';
  };

  window.qtunPause = function () {
    if (!firing) return;
    paused = !paused;
    const btn = document.getElementById('qtun-pause-btn');
    if (btn) btn.textContent = paused ? 'Resume' : 'Pause';
  };

  window.qtunReset = function () {
    firing = false; paused = false;
    probArmed = false; probPeak = 0;
    phiR = 0; phiT = 0;
    setFireReady(true);
    const pb = document.getElementById('qtun-pause-btn');
    if (pb) pb.textContent = 'Pause';
    precompute(); initPsi(); p0 = computeNorm(); render();
  };

  window.addEventListener('resize', layout);

})();
