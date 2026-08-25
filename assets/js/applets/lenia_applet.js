(function () {
  'use strict';

  const _cs = getComputedStyle(document.documentElement);
  const _c   = n => _cs.getPropertyValue(n).trim();
  const _rgb  = n => { const h = _c(n).replace('#',''); const v = parseInt(h,16); return [(v>>16)&0xFF,(v>>8)&0xFF,v&0xFF]; };
  const _rgba = (n, a) => { const [r,g,b] = _rgb(n); return `rgba(${r},${g},${b},${a})`; };

  /* ── Creatures — params and RLE cell patterns from the Lenia species
        library (Chakazul/Lenia, animals.json), all kn=1, gn=1 ── */
  const CREATURES = [
    { label: 'Orbium', name: 'Orbium unicaudatus',
      R: 13, T: 10, b: [1], m: 0.15, s: 0.015,
      cells: '7.MD6.qL$6.pKqEqFURpApBRAqQ$5.VqTrSsBrOpXpWpTpWpUpCrQ$4.CQrQsTsWsApITNPpGqGvL$3.IpIpWrOsGsBqXpJ4.LsFrL$A.DpKpSpJpDqOqUqSqE5.ExD$qL.pBpTT2.qCrGrVrWqM5.sTpP$.pGpWpD3.qUsMtItQtJ6.tL$.uFqGH3.pXtOuR2vFsK5.sM$.tUqL4.GuNwAwVxBwNpC4.qXpA$2.uH5.vBxGyEyMyHtW4.qIpL$2.wV5.tIyG3yOxQqW2.FqHpJ$2.tUS4.rM2yOyJyOyHtVpPMpFqNV$2.HsR4.pUxAyOxLxDxEuVrMqBqGqKJ$3.sLpE3.pEuNxHwRwGvUuLsHrCqTpR$3.TrMS2.pFsLvDvPvEuPtNsGrGqIP$4.pRqRpNpFpTrNtGtVtStGsMrNqNpF$5.pMqKqLqRrIsCsLsIrTrFqJpHE$6.RpSqJqPqVqWqRqKpRXE$8.OpBpIpJpFTK!' },
    { label: 'Gyrorbium', name: 'Gyrorbium gyrans',
      R: 13, T: 10, b: [1], m: 0.156, s: 0.0224,
      cells: '10.EL2QLE$7.TpU2qHqCpXpUpNpFL$4.JrVtTuKuPuKtLrXqTqHqCpPpDG$3.qWtDqRpKqEsMuXvBtGrApXpUpSpIO$2.rQrN4.pAuAvRtTrIpUpIpKpFO$.pSsM6.tJwFuNsPsFrVpPpDL$.uFB6.tJ2yO2yLyOyDsKL$pDuC6.pFxW3yOwIwD2xPqH$rNtV5.EsMxCyIyOwXtJsMtJwFuX$sHuSV3.EpDvOwFxEwQsRqR2qHsFvWE$rQvJsWpPQpKpSqCvEvBuCpD3.BpDtGrQ$pXuKvMuPtLsWsCrIuCtBrS6.qWrQ$EsKvEwXyBwLtVrVsCrDqH6.pXrG$.qHtVxJyOwQrQpNqJpPV6.qJqE$.JsUxMyOrX10.pFqRJ$2.rQxPwIpI9.pKqJT$2.qJxEuPpKB7.qCpP$2.EvOvMpPO5.TrGqH$3.sCyOqEpIOEBOqHqEsRtG$4.xMsMqJqCpXqJqRpIqOuCtBsF$5.xPrAqTqMpSE.rSsMrLqRqHV.TpS$6.vErDE2.VpPB$7.pIrNqHpKQ!' },
    { label: 'Scutium', name: 'Scutium solidus',
      R: 13, T: 10, b: [1], m: 0.29, s: 0.045,
      cells: '5.pGQ$6.sUsDqRR$4.VpXrJwKvNtXrW$2.ApVrGrWsIwKyOyDwTuNO$2.qArOrIqIpPpTxH2yOyIvKqG$rVpIrNrJpH4.xP3yOvFqWA$sKvNsKqE5.pI4yOuHqPC$sNxBuMpD5.JuN3yOwXsUpN$sGxFyOT5.pCtIyF3yOuKqUH$rPwXyOxT5.qKtSxS3yOvIrQT$pQwK2yOtWXJXqHsGuWxW3yOvPsApC$.vQ3yOuVsLsGsXuMwOyK3yOvHrVpC$.rKyK3yOxEwBwCwXyE3yOxQuMrJU$.WvC5yOyM5yOwFtJqOL$2.rVwO9yOwWuMsApOC$2.pLtBwK6yOyJwQuTsQqLP$3.qHsWvEwUxTyCxRwWvRuGsNqSpA$3.HqBrUtHuGuPuLtVsXrSqIXB$4.DpDqEqW2rIqWqGpJN$6.ENTUPHA!' },
    { label: 'Paraptera', name: 'Paraptera arcus labens',
      R: 13, T: 10, b: [1], m: 0.347, s: 0.057,
      cells: '10.C2DC$7.KpBpOpXqDqCpUpMpAMD$5.JpOqVrWsQtAtBsVsMrXrGqMpPTE$4.TqTsUuLvLvXwAvSvIuUuFtOsSrTqRpOPA$3.SrPuKwSyH2yOyJxRxEwPwDvNuTtUsSrMqGpCH$3.rOvFyJ6yOyMyGxXxLwRvTuQtLsEqWpRQA$2.qMuU13yOxRwNvHuCsUrOqGpCG$.JtGxX6yO2yDyL5yOyGxAvTuPtKsDqVpON$.qOvQxQyM4yOwRwBwMxLyI5yOyHxBvXuWtWsSrKqBVD$.sMvLwXyI3yOvNtJtLuLvUxGyH4yOyLxHwLvVvHuKtFrWqMpFG$VsOuHwIyF2yOwKrFqJrGsOuCvRxAxTyCyBxWwUwCvXwDwEvSuXtSsHqUpMJ$pLrCtHvUyB2yOrL.BpIqUsLtVvCvRvTvNuNtTtUuOvSwRxBwOvLuCsPrDpRN$.pXsOvKxTyOvJK3.pEqQrWsRsVsPrSqVqTrOtCvFxGyJyExEvVuJsWrKpVP$.pDsAuVxIyDrSpF4.SpQpWpQW2.JpTsAuWxX3yOxNwCuPtDrOqAQ$.MrKuDwNsUrUqMpAD10.NrLvGxV4yOxNwBuRtFrQqBR$2.qRtEqUqMrJrNqWqDpKpBMC7.rQvHxR4yOyIxAvTuPtFrRqAQ$2.pSK.RpWqUrJrNrOrHqSpRN6.sFvPxRyKyN2yOxWwXwHvOuPtGrQqAQ$7.UpPqGqRrCrHrEqKpFH3.pXtXwRxSxUxOxHwUwAvUwBwDvRuRtGrPpXP$10.FUpOqMrCrAqKpVpNqTuTxPyJxPwGvDuLtMtJuDvIwKwRwBuTtFrOpVM$13.QpUqSrFrMsEwO3yOwQtKrUqOqIrBsLuOwRxSxKwHuUtDrLpRH$15.SpOqEuTxUyI2yOuPpW3.pDrHuDxF2yOxUwHuRtArDpHB$17.rCuWvRwXyHyGrT5.qOuFxT3yOxUwDuJsOqPS$17.rTsMtSvLxMvBpW5.qNuRxN4yOxKvQtVrWpSD$16.FpIqLsDuGtVrOqMG4.rAuKxF4yOyMwQuXtBqTQ$18.TqUrFpHqXrDpD4.qWuAwNyD4yOxLvUuBrTpI$19.K2.qArFqJN3.qNtHvJwXyC3yOyDwQuXsQpXD$22.TqLrBqDM2.pWsDtUvHwPxTyM2yOxLvTtKqNH$23.pCqQrDqBI.VqPsCtMuVwKxV2yOyHwPuArAM$24.pJrDrFpH2.VqHrQtFvExI3yOxJuMrEO$24.DqGrPqIF2.KpXrWuLxJ3yOxXuPrDM$25.pIrJrDV3.NrEvAyL3yOyGuLqSH$25.JqQrNqCJ2.FsQxX4yOxVtTqCB$26.pQrHrFpWUJrWyI5yOwTsOpE$26.KqIrMrNrDtN7yOuXqWF$27.pAqPrTuF7yOwQsMpA$28.pCrNxMxSxTxRxPxSyCxGtIpQ$29.uHvD3vJvOwFwOtFpQ$28.qDsBsIsKsLsStMuPrXV$28.WpHpKpLpRqKrQqA$33.L!' },
    { label: 'Helicium', name: 'Helicium solidus',
      R: 13, T: 10, b: [1], m: 0.35, s: 0.06,
      cells: '9.BG2LID$7.GXpSqIqOqMqFpRpEPE$5.DXqGrOsOtEtHtEsUsKrXrJqPpRTC$4.JpPrKtFuQvMvUvWvX2wAvVvMuRtLrUqCO$3.KqBsEuMwJxNyByI6yOyHwKuErUpRH$2.FpXsKvDxO11yOxIuLrUpRI$2.pKsBvCyE12yOwDtCqQWB9.O$.KrBuHxV4yOyM2yGyM4yOyCuRsNrLpUA9.pEqKpHF$.pKsSwO4yOyFwXwGwJxCyC8yOxOuHrPpWpARNLG2.tArSqLpKXB$.qMuI5yOwEuHtRtWuRwBxO9yOxHuMsTrVrKrBqRqKpWtFuFsSrJqOqM$.rJvO4yOxLtJrQrGrLsGtKvAwWyI7yOyGvSuAsSsBrQ2rOrNsEwNuWtJsHrVpA$BrXwK4yOvKqOpDQVpLqMrVtNvNxFyG2yOwMvGuOtQsLrGqEpLpFpIqArBrWyGwWvHuAtJrC$.sCwRyK3yOtLB5.JpNrBsOtS3.BOSMD5.QqKyCyHxAvPuTsXQ$.rSwJxVyL2yOsP25.JvQyOyGxBwEuMqC$.qXvMwVxT2yOtO26.tP2yOyBxHvQrC$.pTuFvPwNxVyOvWX14.pPpB9.tH3yOyCwHrQ$.HsNuFvBwLxXyBrEpLLB.CJXpRqFqEpOUpJvUwBuPsVrDpQL4.pCuN3yOyNwJrU$2.qNsRtIuQwHxSrVrMqVqKqJqTrPsUuKwGyG5yOyDwSuUtArMqKpTpPqDrSwO4yOvVrN$2.JrGrOsQuEvSrPrDrKrPrVsHtBuFwF9yOxQvVuHtCsIsDsRuPyH4yOuVqU$3.pUpWqQrWtJsNpApNpVqEqLqSrNtAvQ9yOyCwRvMuSuQvJxG4yOxOtMpV$4.GRpQqWsD8.pCrDtRwFyByOyHyG3yOyIxQ2xDxT5yOvLrVU$7.LpMP9.ASqJsPvN6yOyMyN5yOwPtEqEB$19.FpLrNuFxH12yOwUtUrAS$19.MqDsKvAxP10yOyDwGtTrGpE$19.GpPrOtSvRxGyFyOyLyDxTxNxKxDwLvCtDqXpEB$20.GpJqSsDtEtUuHuNuRuVuXuUuKtKrWqGUA$22.FUpNqFqUrJrVsFsHsCrJqJpFH$25.FPpBpL2pQpJUG$28.ACB!' },
    { label: 'Hydrogeminium', name: 'Hydrogeminium natans',
      R: 18, T: 10, b: [1/2, 1, 2/3], m: 0.26, s: 0.036,
      cells: '37.BEG3JGE$35.JpApPqCqJqOqMqEpUpKpAQG$33.BVqCrLsPtLtVuCtTtJsPrVrDqJpPTB$33.pIrAsRuIvWwXxMxPxHwQvTuStTsUrSqOpFG$32.pKrXuIwFxU6yOyLxJwFvBtVsPrDpNJ$31.TsCvT11yOyIwXvTuNsWrDpKG$31.rIwA14yOxHvWuKsRqTpDB$30.rQwS16yOxHvRuArXpXL$29.tG7yOyL4yOyL6yOwXvBsWqTpDB$28.tT7yOxUwSwIwVxHxJxMxUyI4yOyIwIuArNpNG$27.sW7yOxMvOtQsPtBuIvJwAwQxPyL4yOxPvBsHpXL$26.qTyL4yOxMvMuUtTrXqCpKqErStLuSvRwVyB5yOwFtDqMQ$25.pKuPyD3yOwVpS7.pIrQtJuSwDxPyL4yOxHuCrAV$24.QrXvOyB3yOqH9.pDrLtGvBxCyL4yOyBuSrLpA$23.rItDvMxU3yOtV11.pArGtLwIyI4yOyIvBrQpD$22.BuIwSyI2yOvE14.pNrXvMyG4yOyGvBrQpD$23.tJyD2yOpU15.VrNvEyB4yOxPuKrIpA$22.pItVyL2yO16.pSsRwAyD4yOwQtLqRQ$21.rLvJxW2yOsP16.rLvExPyL4yOvMsKpXG$20.tOxJ3yOwLsMJ15.sWxC5yOyBuCrIpF$20.wDyLyOyGuCsHrLpU15.tQyB5yOwLsPqEO$19.sCwS2yOtTpDLpNqHO14.uCyL5yOuFrApDB$16.JtTwNxCyG2yOtG2.BVG13.qRwF5yOuXrQpIG$16.rL6yOxUpD15.qJvRyB5yOuNrQpF$16.tJ7yOvO15.7yOtVqOV$14.pFrX3yOuXTtT3yOG13.rQ6yOvBpF$12.uI6yOtD2.xEyOxH14.6yOvJQ$9.JqEtB8yO2.JsKqR9.pFqErNuC3yOyLyGyLsU$9.tG4yOyD.rG3yOtO8.qCqE3.qT6yOxEtLB$8.wV5yO3.wS2yOtL7.qRwQuFE2.sU5yOxMtT$7.sCyL4yOwFJ3.xUxR8.vG2yOrS2.2yOuUuXyGxEtQ$4.tTvWuFvT2yO2.vWvMpU10.rA2.sR7yOsWsCvMqO$4.vGwNwLxRyOpI2.pIxWuC9.2yOrN.J7yOqM$4.uPwXxMyIyOV2.EyOuS8.uF2yOwV2.6yOuS$4.uSuUvB2yOqEB.JE9.vG3yOrArV2yOuIuPyOwV$4.sP2vRwXyOyGqJL8.O3.yG6yOrLqMvB$3.pFtT2xM3yOsMJ7.rGtQ3.xH6yO$3.rAvM4yOwIqW8.sFwDtQ2.xE5yOrG$.pPqTsHvOyL2yOwQsKV8.pAtT2uCwDyL4yOuP$rGtBsRsUuPxHyOwSuFqWO9.rIvJyB4yOvJuCpFJ$qJtQtOtBuAuFuUuPtTrQpIGEVvTyOV3.rIxC4yOuAsPrGpPB$.rGsRrGpUqWsCsWtQtGrQqCpPqWvWyOsKB2.uI4yOtJrAV$5.pKrDsFtQvMvTuUtDrQqTsW2uIvOyI4yOuUrL$5.OrLtOvGxWyLyOwVsHpPpNyD7yOvRpI$6.sKxC2yIyLyOxWuPrIuN7yOuUqJ$5.QuNyB2yOyLyIyL7yOyIpK$5.GsUxJ2yOyGwSvWwNxW4yOyLG$6.pKuAxHxUwNuCsRtBuSxCyG2yLsK$8.qJsRsWsCrApXqCsUqWG$10.TqHqMpN$11.QqHpK!' },
  ];

  /* Lenia RLE: digits = run count, '.' = 0, 'A'..'X' = 1..24, 'p'..'y' prefix
     extends to 25..255, 'o' = 255, '$' = end of row, '!' = end of pattern.
     Values are stored /255. Rows are returned bottom-up for GL texture order. */
  function rleDecode(st) {
    const rows = [[]];
    let count = '', last = '';
    st = st.replace(/!$/, '');
    for (const ch of st) {
      if (ch >= '0' && ch <= '9') { count += ch; continue; }
      if (ch >= 'p' && ch <= 'y') { last = ch; continue; }
      const n = count ? parseInt(count, 10) : 1;
      if (ch === '$') {
        for (let k = 0; k < n; k++) rows.push([]);
      } else {
        let v;
        if (ch === '.' || ch === 'b') v = 0;
        else if (ch === 'o') v = 255;
        else if (last) v = (last.charCodeAt(0) - 112) * 24 + (ch.charCodeAt(0) - 65 + 25);
        else v = ch.charCodeAt(0) - 64;
        const row = rows[rows.length - 1];
        for (let k = 0; k < n; k++) row.push(v / 255);
      }
      count = ''; last = '';
    }
    const h = rows.length;
    const w = Math.max(...rows.map(r => r.length));
    const data = new Float32Array(w * h);
    rows.forEach((r, y) => r.forEach((v, x) => { data[(h - 1 - y) * w + x] = v; }));
    return { w, h, data };
  }

  /* ── Simulation state ── */
  const N = 256;                 // grid cells per side (toroidal)
  let creature = 0;              // index into CREATURES
  let R = 13, T = 10, bRings = [1];
  let mu = 0.15, sigma = 0.015;
  let brushMode = 0;             // 0 = stamp, 1 = erase
  let brushRadius = 0.06;        // erase radius, uv units
  let running = false, frameId = null;
  let wasRunning = false;        // stashed while the docs panel is open
  let stamp = null;              // decoded cells of the selected creature
  let listenersAdded = false;
  const wall = new Uint8Array(N * N);   // 255 = wall cell (pinned to A = 0)
  let wallDirty = false;

  /* ── WebGL state ── */
  let canvas, gl;
  let simProg = null, eraseProg, dispProg;
  let tex = [null, null], fbo = [null, null], cur = 0;
  let kernelTex = null, kernelSize = 0;
  let wallTex = null;
  let rampTex, sampNearest, sampLinear, vao;
  let eraseQueue = [];           // {x0,y0,x1,y1}
  let pointerDown = false, lastUV = null;

  const WALL_RGB = _rgb('--gold');   // wall display color

  /* ── Palette ramp (256×1 RGBA8) ── */
  const RAMP_SIZE = 256;
  const rampPix = new Uint8Array(RAMP_SIZE * 4);
  (function () {
    const stops = [
      _rgb('--near-black'),
      _rgb('--teal-dark'),
      _rgb('--teal-light'),
      _rgb('--cyan'),
      _rgb('--pink-light'),
      _rgb('--pink-dark'),
    ];
    for (let i = 0; i < RAMP_SIZE; i++) {
      const t = i / (RAMP_SIZE - 1);
      const ft = t * (stops.length - 1);
      const lo = Math.floor(ft), hi = Math.min(lo + 1, stops.length - 1);
      const a = ft - lo;
      rampPix[i*4]   = Math.round(stops[lo][0] + a*(stops[hi][0]-stops[lo][0]));
      rampPix[i*4+1] = Math.round(stops[lo][1] + a*(stops[hi][1]-stops[lo][1]));
      rampPix[i*4+2] = Math.round(stops[lo][2] + a*(stops[hi][2]-stops[lo][2]));
      rampPix[i*4+3] = 255;
    }
  })();

  /* ── Kernel — ring shell K(r) = b[⌊Br⌋]·core(Br mod 1), core = exp(4 − 1/(r(1−r))) ── */
  function kernelCore(r) {
    if (r <= 0 || r >= 1) return 0;
    return Math.exp(4 - 1 / (r * (1 - r)));
  }
  function kernelShell(r, b) {          // unnormalized profile, for plotting too
    if (r >= 1) return 0;
    const B = b.length, Br = B * r;
    return b[Math.min(Math.floor(Br), B - 1)] * kernelCore(Math.min(Br % 1, 1));
  }
  function buildKernel(R_, b) {         // (2R+1)² window, normalized to sum 1
    const size = 2 * R_ + 1;
    const K = new Float32Array(size * size);
    let sum = 0;
    for (let dy = -R_; dy <= R_; dy++) {
      for (let dx = -R_; dx <= R_; dx++) {
        const v = kernelShell(Math.hypot(dx, dy) / R_, b);
        K[(dy + R_) * size + (dx + R_)] = v;
        sum += v;
      }
    }
    for (let i = 0; i < K.length; i++) K[i] /= sum;
    return { size, data: K };
  }

  /* ── Shaders ── */
  const VERT = `#version 300 es
  out vec2 v_uv;
  void main() {
    vec2 p = vec2(float((gl_VertexID << 1) & 2), float(gl_VertexID & 2));
    v_uv = p;
    gl_Position = vec4(p * 2.0 - 1.0, 0.0, 1.0);
  }`;

  /* Convolution radius is compiled in; the program is rebuilt per creature. */
  function makeSimFrag(R_) {
    return `#version 300 es
  precision highp float;
  uniform sampler2D uState;
  uniform sampler2D uKernel;
  uniform sampler2D uWall;
  uniform vec2  uTexel;
  uniform float uMu, uSigma, uDt;
  in vec2 v_uv;
  out vec4 frag;
  const int KR = ${R_};
  void main() {
    float u = 0.0;
    for (int dy = -KR; dy <= KR; dy++) {
      for (int dx = -KR; dx <= KR; dx++) {
        float w = texelFetch(uKernel, ivec2(dx + KR, dy + KR), 0).r;
        if (w > 0.0) {
          u += w * texture(uState, v_uv + vec2(float(dx), float(dy)) * uTexel).r;
        }
      }
    }
    float a = texture(uState, v_uv).r;
    float g = 2.0 * exp(-(u - uMu) * (u - uMu) / (2.0 * uSigma * uSigma)) - 1.0;
    float wl = texture(uWall, v_uv).r;
    frag = vec4(clamp(a + uDt * g, 0.0, 1.0) * (1.0 - wl), 0.0, 0.0, 1.0);
  }`;
  }

  const ERASE_FRAG = `#version 300 es
  precision highp float;
  uniform sampler2D uState;
  uniform vec2  uP0, uP1;
  uniform float uRadius;
  in vec2 v_uv;
  out vec4 frag;
  float distSeg(vec2 p, vec2 a, vec2 b) {
    vec2 ab = b - a;
    float t = clamp(dot(p - a, ab) / max(dot(ab, ab), 1e-12), 0.0, 1.0);
    return length(p - (a + t * ab));
  }
  void main() {
    float a = texture(uState, v_uv).r;
    float d = distSeg(v_uv, uP0, uP1);
    if (d < uRadius) a *= smoothstep(0.0, 1.0, d / uRadius);
    frag = vec4(a, 0.0, 0.0, 1.0);
  }`;

  const DISP_FRAG = `#version 300 es
  precision highp float;
  uniform sampler2D uState;
  uniform sampler2D uRamp;
  uniform sampler2D uWall;
  uniform vec3  uWallCol;
  in vec2 v_uv;
  out vec4 frag;
  void main() {
    float a = texture(uState, v_uv).r;
    vec3 col = texture(uRamp, vec2(clamp(a, 0.0, 1.0), 0.5)).rgb;
    float wl = texture(uWall, v_uv).r;
    frag = vec4(mix(col, uWallCol, wl), 1.0);
  }`;

  function compile(type, src) {
    const sh = gl.createShader(type);
    gl.shaderSource(sh, src);
    gl.compileShader(sh);
    if (!gl.getShaderParameter(sh, gl.COMPILE_STATUS)) {
      throw new Error('shader: ' + gl.getShaderInfoLog(sh));
    }
    return sh;
  }
  function makeProg(fragSrc) {
    const p = gl.createProgram();
    gl.attachShader(p, compile(gl.VERTEX_SHADER, VERT));
    gl.attachShader(p, compile(gl.FRAGMENT_SHADER, fragSrc));
    gl.linkProgram(p);
    if (!gl.getProgramParameter(p, gl.LINK_STATUS)) {
      throw new Error('link: ' + gl.getProgramInfoLog(p));
    }
    return p;
  }

  function makeStateTextures() {
    for (let i = 0; i < 2; i++) {
      if (tex[i]) { gl.deleteTexture(tex[i]); gl.deleteFramebuffer(fbo[i]); }
      tex[i] = gl.createTexture();
      gl.bindTexture(gl.TEXTURE_2D, tex[i]);
      gl.texStorage2D(gl.TEXTURE_2D, 1, gl.RG16F, N, N);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.REPEAT);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.REPEAT);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
      fbo[i] = gl.createFramebuffer();
      gl.bindFramebuffer(gl.FRAMEBUFFER, fbo[i]);
      gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, tex[i], 0);
    }
    gl.bindFramebuffer(gl.FRAMEBUFFER, null);
  }

  function uploadKernel(R_, b) {
    const { size, data } = buildKernel(R_, b);
    if (kernelTex) gl.deleteTexture(kernelTex);
    kernelTex = gl.createTexture();
    kernelSize = size;
    gl.bindTexture(gl.TEXTURE_2D, kernelTex);
    gl.pixelStorei(gl.UNPACK_ALIGNMENT, 1);
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.R32F, size, size, 0, gl.RED, gl.FLOAT, data);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
  }

  function clearField() {
    for (let i = 0; i < 2; i++) {
      gl.bindFramebuffer(gl.FRAMEBUFFER, fbo[i]);
      gl.viewport(0, 0, N, N);
      gl.clearColor(0.0, 0.0, 0.0, 1.0);
      gl.clear(gl.COLOR_BUFFER_BIT);
    }
    gl.bindFramebuffer(gl.FRAMEBUFFER, null);
    eraseQueue.length = 0;
  }

  /* Write a decoded pattern into the current state texture at uv center,
     split into up to four sub-rectangles across the periodic seams. */
  function blitPattern(pat, u, v) {
    const x0 = ((Math.round(u * N - pat.w / 2) % N) + N) % N;
    const y0 = ((Math.round(v * N - pat.h / 2) % N) + N) % N;
    gl.bindTexture(gl.TEXTURE_2D, tex[cur]);
    gl.pixelStorei(gl.UNPACK_ALIGNMENT, 1);
    const wA = Math.min(pat.w, N - x0), wB = pat.w - wA;
    const hA = Math.min(pat.h, N - y0), hB = pat.h - hA;
    const parts = [
      { dx: x0, dy: y0, sx: 0,  sy: 0,  w: wA, h: hA },
      { dx: 0,  dy: y0, sx: wA, sy: 0,  w: wB, h: hA },
      { dx: x0, dy: 0,  sx: 0,  sy: hA, w: wA, h: hB },
      { dx: 0,  dy: 0,  sx: wA, sy: hA, w: wB, h: hB },
    ];
    for (const p of parts) {
      if (p.w <= 0 || p.h <= 0) continue;
      const sub = new Float32Array(p.w * p.h * 2);
      for (let y = 0; y < p.h; y++) {
        for (let x = 0; x < p.w; x++) {
          const gi = (p.dy + y) * N + (p.dx + x);
          sub[(y * p.w + x) * 2] = wall[gi] ? 0 : pat.data[(p.sy + y) * pat.w + (p.sx + x)];
        }
      }
      gl.texSubImage2D(gl.TEXTURE_2D, 0, p.dx, p.dy, p.w, p.h, gl.RG, gl.FLOAT, sub);
    }
  }

  /* ── Walls — CPU mask, uploaded to an R8 texture when edited ── */
  function uploadWall() {
    gl.bindTexture(gl.TEXTURE_2D, wallTex);
    gl.pixelStorei(gl.UNPACK_ALIGNMENT, 1);
    gl.texSubImage2D(gl.TEXTURE_2D, 0, 0, 0, N, N, gl.RED, gl.UNSIGNED_BYTE, wall);
    wallDirty = false;
  }

  /* Rasterize a thick segment (uv coords, brushRadius wide) into the mask. */
  function paintWallSegment(x0, y0, x1, y1, on) {
    const r = brushRadius * N;
    const ax = x0 * N, ay = y0 * N, bx = x1 * N, by = y1 * N;
    const minx = Math.floor(Math.min(ax, bx) - r), maxx = Math.ceil(Math.max(ax, bx) + r);
    const miny = Math.floor(Math.min(ay, by) - r), maxy = Math.ceil(Math.max(ay, by) + r);
    const abx = bx - ax, aby = by - ay;
    const ab2 = Math.max(abx * abx + aby * aby, 1e-12);
    for (let y = miny; y <= maxy; y++) {
      for (let x = minx; x <= maxx; x++) {
        const px = x + 0.5, py = y + 0.5;
        let t = ((px - ax) * abx + (py - ay) * aby) / ab2;
        t = Math.max(0, Math.min(1, t));
        const dx = px - (ax + t * abx), dy = py - (ay + t * aby);
        if (dx * dx + dy * dy < r * r) {
          wall[(((y % N) + N) % N) * N + (((x % N) + N) % N)] = on ? 255 : 0;
        }
      }
    }
    wallDirty = true;
  }

  /* ── Random soup — smooth blobs, mean density ≈ 0.25 ── */
  function seedSoup() {
    if (!gl) return;
    const A = new Float32Array(N * N);
    const nb = Math.round((N * N) / (2.2 * R * R));
    for (let k = 0; k < nb; k++) {
      const cx = Math.random() * N, cy = Math.random() * N;
      const r0 = R * (0.4 + 0.6 * Math.random());
      const amp = 0.5 + 0.5 * Math.random();
      const s2 = 2 * (r0 / 2) * (r0 / 2);
      const span = Math.ceil(2 * r0);
      for (let dy = -span; dy <= span; dy++) {
        for (let dx = -span; dx <= span; dx++) {
          const x = ((Math.round(cx) + dx) % N + N) % N;
          const y = ((Math.round(cy) + dy) % N + N) % N;
          A[y * N + x] += amp * Math.exp(-(dx * dx + dy * dy) / s2);
        }
      }
    }
    const buf = new Float32Array(N * N * 2);
    for (let i = 0; i < N * N; i++) buf[i * 2] = wall[i] ? 0 : Math.min(1, A[i]);
    gl.bindTexture(gl.TEXTURE_2D, tex[cur]);
    gl.pixelStorei(gl.UNPACK_ALIGNMENT, 1);
    gl.texSubImage2D(gl.TEXTURE_2D, 0, 0, 0, N, N, gl.RG, gl.FLOAT, buf);
  }

  /* ── Creature selection ── */
  function setCreature(i) {
    creature = i;
    const c = CREATURES[i];
    R = c.R; T = c.T; bRings = c.b; mu = c.m; sigma = c.s;
    stamp = rleDecode(c.cells);
    syncSliders();
    if (gl) {
      uploadKernel(R, bRings);
      if (simProg) gl.deleteProgram(simProg);
      simProg = makeProg(makeSimFrag(R));
      clearField();
      blitPattern(stamp, 0.5, 0.5);
    }
    drawPlots();
  }

  function syncSliders() {
    const set = (id, val, txt) => {
      const el = document.getElementById('lenia-' + id);
      if (el) el.value = val;
      const vs = document.getElementById('lenia-' + id + '-val');
      if (vs) vs.textContent = txt;
    };
    set('mu', mu, mu.toFixed(3));
    set('sigma', sigma, sigma.toFixed(4));
    set('T', T, String(T));
  }

  /* ── GPU passes ── */
  function drawFullscreen() {
    gl.bindVertexArray(vao);
    gl.drawArrays(gl.TRIANGLES, 0, 3);
  }

  function applyErase(op) {
    const dst = 1 - cur;
    gl.bindFramebuffer(gl.FRAMEBUFFER, fbo[dst]);
    gl.viewport(0, 0, N, N);
    gl.useProgram(eraseProg);
    gl.activeTexture(gl.TEXTURE0);
    gl.bindTexture(gl.TEXTURE_2D, tex[cur]);
    gl.bindSampler(0, sampNearest);
    gl.uniform1i(gl.getUniformLocation(eraseProg, 'uState'), 0);
    gl.uniform2f(gl.getUniformLocation(eraseProg, 'uP0'), op.x0, op.y0);
    gl.uniform2f(gl.getUniformLocation(eraseProg, 'uP1'), op.x1, op.y1);
    gl.uniform1f(gl.getUniformLocation(eraseProg, 'uRadius'), brushRadius);
    drawFullscreen();
    cur = dst;
  }

  function simStep() {
    const dst = 1 - cur;
    gl.bindFramebuffer(gl.FRAMEBUFFER, fbo[dst]);
    gl.viewport(0, 0, N, N);
    gl.useProgram(simProg);
    gl.activeTexture(gl.TEXTURE0);
    gl.bindTexture(gl.TEXTURE_2D, tex[cur]);
    gl.bindSampler(0, sampNearest);
    gl.activeTexture(gl.TEXTURE1);
    gl.bindTexture(gl.TEXTURE_2D, kernelTex);
    gl.bindSampler(1, null);
    gl.activeTexture(gl.TEXTURE2);
    gl.bindTexture(gl.TEXTURE_2D, wallTex);
    gl.bindSampler(2, sampNearest);
    gl.uniform1i(gl.getUniformLocation(simProg, 'uState'), 0);
    gl.uniform1i(gl.getUniformLocation(simProg, 'uKernel'), 1);
    gl.uniform1i(gl.getUniformLocation(simProg, 'uWall'), 2);
    gl.uniform2f(gl.getUniformLocation(simProg, 'uTexel'), 1 / N, 1 / N);
    gl.uniform1f(gl.getUniformLocation(simProg, 'uMu'), mu);
    gl.uniform1f(gl.getUniformLocation(simProg, 'uSigma'), sigma);
    gl.uniform1f(gl.getUniformLocation(simProg, 'uDt'), 1 / T);
    drawFullscreen();
    cur = dst;
  }

  function render() {
    gl.bindFramebuffer(gl.FRAMEBUFFER, null);
    gl.viewport(0, 0, canvas.width, canvas.height);
    gl.useProgram(dispProg);
    gl.activeTexture(gl.TEXTURE0);
    gl.bindTexture(gl.TEXTURE_2D, tex[cur]);
    gl.bindSampler(0, sampLinear);
    gl.activeTexture(gl.TEXTURE1);
    gl.bindTexture(gl.TEXTURE_2D, rampTex);
    gl.bindSampler(1, null);
    gl.activeTexture(gl.TEXTURE2);
    gl.bindTexture(gl.TEXTURE_2D, wallTex);
    gl.bindSampler(2, sampLinear);
    gl.uniform1i(gl.getUniformLocation(dispProg, 'uState'), 0);
    gl.uniform1i(gl.getUniformLocation(dispProg, 'uRamp'), 1);
    gl.uniform1i(gl.getUniformLocation(dispProg, 'uWall'), 2);
    gl.uniform3f(gl.getUniformLocation(dispProg, 'uWallCol'),
                 WALL_RGB[0] / 255, WALL_RGB[1] / 255, WALL_RGB[2] / 255);
    drawFullscreen();
  }

  function initGL() {
    gl = canvas.getContext('webgl2', { antialias: false, alpha: false, preserveDrawingBuffer: false });
    if (!gl) return false;
    const ext = gl.getExtension('EXT_color_buffer_float')
             || gl.getExtension('EXT_color_buffer_half_float');
    if (!ext) { gl = null; return false; }

    eraseProg = makeProg(ERASE_FRAG);
    dispProg  = makeProg(DISP_FRAG);

    vao = gl.createVertexArray();

    sampNearest = gl.createSampler();
    gl.samplerParameteri(sampNearest, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
    gl.samplerParameteri(sampNearest, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
    gl.samplerParameteri(sampNearest, gl.TEXTURE_WRAP_S, gl.REPEAT);
    gl.samplerParameteri(sampNearest, gl.TEXTURE_WRAP_T, gl.REPEAT);
    sampLinear = gl.createSampler();
    gl.samplerParameteri(sampLinear, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
    gl.samplerParameteri(sampLinear, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
    gl.samplerParameteri(sampLinear, gl.TEXTURE_WRAP_S, gl.REPEAT);
    gl.samplerParameteri(sampLinear, gl.TEXTURE_WRAP_T, gl.REPEAT);

    rampTex = gl.createTexture();
    gl.bindTexture(gl.TEXTURE_2D, rampTex);
    gl.texStorage2D(gl.TEXTURE_2D, 1, gl.RGBA8, RAMP_SIZE, 1);
    gl.texSubImage2D(gl.TEXTURE_2D, 0, 0, 0, RAMP_SIZE, 1, gl.RGBA, gl.UNSIGNED_BYTE, rampPix);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);

    wallTex = gl.createTexture();
    gl.bindTexture(gl.TEXTURE_2D, wallTex);
    gl.texStorage2D(gl.TEXTURE_2D, 1, gl.R8, N, N);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.REPEAT);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.REPEAT);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);

    makeStateTextures();
    return true;
  }

  /* ── Animation loop ── */
  function loop() {
    if (gl) {
      if (wallDirty) uploadWall();
      for (let i = 0; i < eraseQueue.length; i++) applyErase(eraseQueue[i]);
      eraseQueue.length = 0;
      if (pointerDown && lastUV && brushMode === 1) {
        applyErase({ x0: lastUV[0], y0: lastUV[1], x1: lastUV[0], y1: lastUV[1] });
      }
      if (running) simStep();
      render();
    }
    frameId = requestAnimationFrame(loop);
  }

  /* ── Ctrl-panel plots — kernel profile K(r) and growth map G(u) ── */
  function sizePlots() {
    const dpr = window.devicePixelRatio || 1;
    ['lenia-kplot', 'lenia-gplot'].forEach(id => {
      const el = document.getElementById(id);
      if (!el) return;
      const w = el.clientWidth, h = el.clientHeight;
      if (w > 0) { el.width = Math.round(w * dpr); el.height = Math.round(h * dpr); }
    });
    drawPlots();
  }

  function plotFrame(ctx, w, h) {
    ctx.clearRect(0, 0, w, h);
    ctx.strokeStyle = _c('--border-mid');
    ctx.lineWidth = 1;
    ctx.strokeRect(0.5, 0.5, w - 1, h - 1);
  }

  function drawPlots() {
    const PADX = 6, PADY = 5;
    const kc = document.getElementById('lenia-kplot');
    if (kc && kc.width > 0) {
      const ctx = kc.getContext('2d'), w = kc.width, h = kc.height;
      plotFrame(ctx, w, h);
      let kmax = 0;
      const pts = [];
      for (let i = 0; i <= 220; i++) {
        const r = i / 220;
        const v = kernelShell(r, bRings);
        pts.push(v);
        if (v > kmax) kmax = v;
      }
      ctx.beginPath();
      pts.forEach((v, i) => {
        const x = PADX + (i / 220) * (w - 2 * PADX);
        const y = h - PADY - (v / kmax) * (h - 2 * PADY);
        i === 0 ? ctx.moveTo(x, y) : ctx.lineTo(x, y);
      });
      ctx.strokeStyle = _c('--teal-light');
      ctx.lineWidth = 1.5;
      ctx.stroke();
      ctx.lineTo(PADX + (w - 2 * PADX), h - PADY);
      ctx.lineTo(PADX, h - PADY);
      ctx.closePath();
      ctx.fillStyle = _rgba('--teal-light', 0.15);
      ctx.fill();
    }
    const gc = document.getElementById('lenia-gplot');
    if (gc && gc.width > 0) {
      const ctx = gc.getContext('2d'), w = gc.width, h = gc.height;
      plotFrame(ctx, w, h);
      const UMAX = 0.6;
      const y0 = h / 2;                       // G = 0 line
      ctx.strokeStyle = _rgba('--text-dim', 0.4);
      ctx.setLineDash([3, 4]);
      ctx.beginPath();
      ctx.moveTo(PADX, y0);
      ctx.lineTo(w - PADX, y0);
      ctx.stroke();
      ctx.setLineDash([]);
      const xMu = PADX + (mu / UMAX) * (w - 2 * PADX);
      ctx.strokeStyle = _c('--amber');
      ctx.beginPath();
      ctx.moveTo(xMu, PADY);
      ctx.lineTo(xMu, h - PADY);
      ctx.stroke();
      ctx.beginPath();
      for (let i = 0; i <= 220; i++) {
        const u = (i / 220) * UMAX;
        const g = 2 * Math.exp(-(u - mu) * (u - mu) / (2 * sigma * sigma)) - 1;
        const x = PADX + (i / 220) * (w - 2 * PADX);
        const y = y0 - g * (h / 2 - PADY);
        i === 0 ? ctx.moveTo(x, y) : ctx.lineTo(x, y);
      }
      ctx.strokeStyle = _c('--pink-light');
      ctx.lineWidth = 1.5;
      ctx.stroke();
    }
  }

  /* ── Shell wiring ── */
  const CREATURE_BTNS = CREATURES.map((c, i) =>
    `<button class="applet-shell-btn lenia-creature${i === 0 ? ' active' : ''}" data-i="${i}" title="${c.name}">${c.label}</button>`
  ).join('');

  const shell = new AppletShell({
    id:    'lenia',
    title: 'Lenia &mdash; Continuous Life',
    gap:   0,

    headerBtns: `<button class="applet-shell-header-btn" onclick="leniaReset()">Reset</button><button class="applet-shell-header-btn" id="lenia-pause-btn" onclick="leniaTogglePause()">Pause</button>`,

    docs: {
      whatis: `What is the smallest thing that deserves to be called an organism? Conway's Game of Life [gardner1970] made the question computational: on a lattice of dead-or-alive cells with a single update rule, gliders crawl, guns fire, and patterns consume one another. But its creatures are brittle — displace one cell and the glider shatters. Life of the wet kind is not like that; it is continuous, noisy, and forgiving.¶Lenia [chan2019], Bert Chan's generalization of the Game of Life, is what happens when the lattice is sent to the continuum the way a physicist would send it: states, space, and time all become smooth (SmoothLife [rafler2011] was an important waypoint). A field $A(\\mathbf{x}) \\in [0,1]$ senses its surroundings through a ring-shaped kernel $K$, and the sensed density $u = K * A$ feeds a growth map $G$ that pushes the field up where $u$ sits near an optimum $\\mu$ and down elsewhere:
$$A_{t+\\Delta t} = \\Big[\\, A_t + \\Delta t \\; G\\big(K * A_t\\big) \\Big]_0^1, \\qquad G(u) = 2\\, e^{-(u-\\mu)^2 / 2\\sigma^2} - 1,$$
where $[\\,\\cdot\\,]_0^1$ clips to the unit interval and $\\Delta t = 1/T$.¶One line of dynamics, and yet its attractors form a bestiary of over five hundred catalogued "species" [chan2020] — localized, self-repairing patterns that glide, spin, and swim, named with Linnaean binomials as if they were fauna. The six included here are stable solitons of their parameter sets: Orbium, the gliding mascot of the entire system; Gyrorbium, which orbits as it spins; Scutium, a marching shield; Paraptera, a broad two-winged flyer; Helicium, a rotating spiral; and Hydrogeminium, a three-ringed swimmer with visibly churning internals. Whether such patterns constitute a minimal model of life is an open question; that they make the question quantitative is Lenia's charm.`,

      howto: `The canvas shows the field $A$ on a periodic domain through the palette ramp, near-black at $A=0$ rising through teal into pink; the convolution runs on the GPU. Each creature button clears the world, loads that species' kernel rings, growth parameters, and time resolution, then stamps one specimen in the center. The two plots at the bottom of the panel show the world it now inhabits: the kernel profile $K(r)$ and the growth map $G(u)$, with the amber line marking $\\mu$.¶In Stamp mode a click deposits another specimen; in Erase mode the pointer carves holes, with the brush slider setting the width — creatures regrow small wounds, which is much of the fun. Wall raises gold barriers pinned to zero, which the rule reads as permanently empty space no pattern can grow into; they survive creature switches and Soup; Reset or the eraser removes them. Soup floods the world with random blobs and lets the current rule decide what survives.¶The sliders move the growth map underneath the creatures: nudge $\\mu$ or $\\sigma$ by a few thousandths and a glider fattens, stalls, or dissolves; $T$ sets the time resolution $\\Delta t = 1/T$, so higher $T$ is finer and slower. Reset restores a single centered specimen when an experiment ends badly, and Pause freezes time, though the brush still works.`,

      references: ['gardner1970', 'rafler2011', 'chan2019', 'chan2020'],
    },

    onDocsOpen:  function () { wasRunning = running; running = false; },
    onDocsClose: function () { running = wasRunning; },

    ctrlHTML: `
      <div class="applet-shell-ctrl-section">
        <div class="applet-shell-ctrl-title">Creature</div>
        <div class="applet-shell-btn-row">${CREATURE_BTNS}</div>
      </div>
      <div class="applet-shell-ctrl-section">
        <div class="applet-shell-ctrl-title">Brush</div>
        <div class="applet-shell-btn-row">
          <button class="applet-shell-btn lenia-brush active" data-mode="0">Stamp</button>
          <button class="applet-shell-btn lenia-brush" data-mode="1">Erase</button>
          <button class="applet-shell-btn lenia-brush" data-mode="2">Wall</button>
          <button class="applet-shell-btn" id="lenia-soup">Soup</button>
        </div>
        <div class="applet-shell-slider-row">
          <span class="applet-shell-side">Fine</span>
          <input type="range" id="lenia-brush-size" min="0.02" max="0.15" step="0.005" value="0.06">
          <span class="applet-shell-side">Broad</span>
        </div>
      </div>
      <div class="applet-shell-ctrl-section">
        <div class="applet-shell-ctrl-title">Growth Center &mu;</div>
        <div class="applet-shell-slider-row">
          <span class="applet-shell-side">Low</span>
          <input type="range" id="lenia-mu" min="0.05" max="0.45" step="0.001" value="0.15">
          <span class="applet-shell-side">High</span>
          <span class="applet-shell-val" id="lenia-mu-val">0.150</span>
        </div>
      </div>
      <div class="applet-shell-ctrl-section">
        <div class="applet-shell-ctrl-title">Growth Width &sigma;</div>
        <div class="applet-shell-slider-row">
          <span class="applet-shell-side">Narrow</span>
          <input type="range" id="lenia-sigma" min="0.004" max="0.09" step="0.0005" value="0.015">
          <span class="applet-shell-side">Wide</span>
          <span class="applet-shell-val" id="lenia-sigma-val">0.0150</span>
        </div>
      </div>
      <div class="applet-shell-ctrl-section">
        <div class="applet-shell-ctrl-title">Time Resolution T</div>
        <div class="applet-shell-slider-row">
          <span class="applet-shell-side">Coarse</span>
          <input type="range" id="lenia-T" min="2" max="60" step="1" value="10">
          <span class="applet-shell-side">Fine</span>
          <span class="applet-shell-val" id="lenia-T-val">10</span>
        </div>
      </div>
      <div class="applet-shell-ctrl-section">
        <div class="applet-shell-ctrl-title">Kernel K(r) &middot; Growth G(u)</div>
        <canvas id="lenia-kplot" style="width:100%;height:56px;display:block;"></canvas>
        <canvas id="lenia-gplot" style="width:100%;height:56px;display:block;"></canvas>
      </div>
    `,

    onOpen: function ({ canvas: c, S }) {
      canvas = c;
      canvas.style.cursor = 'crosshair';
      canvas.style.touchAction = 'none';

      if (!gl) {
        if (!initGL()) {
          const msg = document.createElement('div');
          msg.style.cssText = 'padding:24px;color:var(--text-bright);font-size:15px;';
          msg.textContent = 'WebGL2 with float render targets is not available in this browser.';
          canvas.parentNode.replaceChild(msg, canvas);
          return;
        }
        setCreature(0);
      }

      if (!listenersAdded) {
        listenersAdded = true;
        const uvCoords = function (e) {
          const rect = canvas.getBoundingClientRect();
          return [
            (e.clientX - rect.left) / rect.width,
            1.0 - (e.clientY - rect.top) / rect.height,
          ];
        };
        canvas.addEventListener('pointerdown', function (e) {
          canvas.setPointerCapture(e.pointerId);
          pointerDown = true;
          lastUV = uvCoords(e);
          if (brushMode === 0) {
            if (stamp && gl) blitPattern(stamp, lastUV[0], lastUV[1]);
          } else if (brushMode === 1) {
            eraseQueue.push({ x0: lastUV[0], y0: lastUV[1], x1: lastUV[0], y1: lastUV[1] });
            paintWallSegment(lastUV[0], lastUV[1], lastUV[0], lastUV[1], false);
          } else {
            paintWallSegment(lastUV[0], lastUV[1], lastUV[0], lastUV[1], true);
          }
        });
        canvas.addEventListener('pointermove', function (e) {
          if (!pointerDown || brushMode === 0) return;
          const uv = uvCoords(e);
          if (brushMode === 1) {
            eraseQueue.push({ x0: lastUV[0], y0: lastUV[1], x1: uv[0], y1: uv[1] });
            paintWallSegment(lastUV[0], lastUV[1], uv[0], uv[1], false);
          } else {
            paintWallSegment(lastUV[0], lastUV[1], uv[0], uv[1], true);
          }
          lastUV = uv;
        });
        const up = function () { pointerDown = false; lastUV = null; };
        canvas.addEventListener('pointerup', up);
        canvas.addEventListener('pointercancel', up);
      }

      sizePlots();
      running = true;
      const pb = document.getElementById('lenia-pause-btn');
      if (pb) { pb.textContent = 'Pause'; pb.classList.remove('active'); }
      if (!frameId) loop();
    },

    onClose: function () {
      running = false;
      pointerDown = false;
      if (frameId) { cancelAnimationFrame(frameId); frameId = null; }
      const pb = document.getElementById('lenia-pause-btn');
      if (pb) { pb.textContent = 'Pause'; pb.classList.remove('active'); }
    },

    onResize: function ({ canvas: c, S }) {
      canvas = c;   // shell has reset canvas.width/height; render() picks it up
      sizePlots();
    },
  });

  /* ── Global entry points ── */
  window.leniaOpen  = () => shell.open();
  window.leniaClose = () => shell.close();

  window.leniaReset = function () {
    if (!gl) return;
    wall.fill(0);
    wallDirty = true;
    clearField();
    if (stamp) blitPattern(stamp, 0.5, 0.5);
  };

  window.leniaTogglePause = function () {
    running = !running;
    const pb = document.getElementById('lenia-pause-btn');
    if (pb) {
      pb.textContent = running ? 'Pause' : 'Resume';
      pb.classList.toggle('active', !running);
    }
  };

  /* ── Control listeners ── */
  document.getElementById('lenia-mu').addEventListener('input', function () {
    mu = parseFloat(this.value);
    document.getElementById('lenia-mu-val').textContent = mu.toFixed(3);
    drawPlots();
  });
  document.getElementById('lenia-sigma').addEventListener('input', function () {
    sigma = parseFloat(this.value);
    document.getElementById('lenia-sigma-val').textContent = sigma.toFixed(4);
    drawPlots();
  });
  document.getElementById('lenia-T').addEventListener('input', function () {
    T = parseInt(this.value);
    document.getElementById('lenia-T-val').textContent = String(T);
  });
  document.getElementById('lenia-brush-size').addEventListener('input', function () {
    brushRadius = parseFloat(this.value);
  });
  document.getElementById('lenia-soup').addEventListener('click', function () {
    seedSoup();
    this.blur();
  });

  document.querySelectorAll('.lenia-brush').forEach(btn => {
    btn.addEventListener('click', function () {
      document.querySelectorAll('.lenia-brush').forEach(b => b.classList.remove('active'));
      this.classList.add('active');
      brushMode = parseInt(this.dataset.mode);
    });
  });

  document.querySelectorAll('.lenia-creature').forEach(btn => {
    btn.addEventListener('click', function () {
      document.querySelectorAll('.lenia-creature').forEach(b => b.classList.remove('active'));
      this.classList.add('active');
      setCreature(parseInt(this.dataset.i));
    });
  });

})();
