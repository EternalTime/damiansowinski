(function () {
'use strict';

const _cs  = getComputedStyle(document.documentElement);
const _c   = n => _cs.getPropertyValue(n).trim();

/* ── Physics (ℏ=1, E=k²/2M) ── */
const MASS    = 3.0;
const K0      = 24.0;
const SIGMA_K = 2.5;

/* ── Grid ── */
const NX   = 2048;
const XMIN = -18.0;
const XMAX =  22.0;
const XSPAN = XMAX - XMIN;
const DX    = XSPAN / NX;

/* ── Time stepping ── */
const DT = 0.002;

/* ── Absorbing boundary ── */
const ABS_FRAC = 0.20;

/* ── Display window ── */
const XDISP_MIN  = -10.0;
const XDISP_MAX  =  14.0;
const XDISP_SPAN = XDISP_MAX - XDISP_MIN;

/* ── Rendering ── */
const BARRIER_MAX = 180.0;
const YMAX_PSI  = 5.0 / 4.0;
const YMAX_PROB = 28.8 / 0.75 / 48.0;

/* ── State ── */
let canvas, ctx, simW, simH;
let frameId = null;
let firing  = false;
let canFire = true;
let paused  = false;
let probArmed = false, probPeak = 0;
let speedVal = 24;
let frameSkipCount = 0;
let barrierH = 60.0;
let barrierW = 0.8;
let running  = false;

let psiRe, psiIm;
let p0 = 1.0, phiR = 0.0, phiT = 0.0;

let xArr, vHalfRe, vHalfIm, kPhaseRe, kPhaseIm, absorb;

/* ── FFT (Cooley-Tukey in-place radix-2) ── */
function fft(re, im, inverse) {
  const n = re.length;
  for (let i=1,j=0; i<n; i++) {
    let bit = n>>1;
    for (; j&bit; bit>>=1) j^=bit;
    j^=bit;
    if (i<j) { let t=re[i];re[i]=re[j];re[j]=t; t=im[i];im[i]=im[j];im[j]=t; }
  }
  for (let len=2; len<=n; len<<=1) {
    const ang = (inverse?1:-1)*2*Math.PI/len;
    const wRe=Math.cos(ang), wIm=Math.sin(ang);
    for (let i=0; i<n; i+=len) {
      let uRe=1, uIm=0;
      for (let j=0; j<len/2; j++) {
        const eRe=re[i+j+len/2]*uRe-im[i+j+len/2]*uIm;
        const eIm=re[i+j+len/2]*uIm+im[i+j+len/2]*uRe;
        re[i+j+len/2]=re[i+j]-eRe; im[i+j+len/2]=im[i+j]-eIm;
        re[i+j]+=eRe; im[i+j]+=eIm;
        const nuRe=uRe*wRe-uIm*wIm; uIm=uRe*wIm+uIm*wRe; uRe=nuRe;
      }
    }
  }
  if (inverse) { for (let i=0;i<n;i++){re[i]/=n;im[i]/=n;} }
}

function buildV() {
  const v = new Float64Array(NX);
  for (let i=0; i<NX; i++) { const x=xArr[i]; if (x>=0&&x<=barrierW) v[i]=barrierH; }
  return v;
}

function precompute() {
  xArr    = new Float64Array(NX);
  vHalfRe = new Float64Array(NX);
  vHalfIm = new Float64Array(NX);
  kPhaseRe = new Float64Array(NX);
  kPhaseIm = new Float64Array(NX);
  absorb   = new Float64Array(NX);
  for (let i=0;i<NX;i++) xArr[i]=XMIN+i*DX;
  const v = buildV();
  for (let i=0;i<NX;i++) {
    const ph=-v[i]*DT/2;
    vHalfRe[i]=Math.cos(ph); vHalfIm[i]=Math.sin(ph);
  }
  const dk=2*Math.PI/XSPAN;
  for (let i=0;i<NX;i++) {
    const ki=i<=NX/2?i:i-NX, k=ki*dk;
    const ph=-(k*k/(2*MASS))*DT;
    kPhaseRe[i]=Math.cos(ph); kPhaseIm[i]=Math.sin(ph);
  }
  const absLen=Math.floor(ABS_FRAC*NX);
  for (let i=0;i<NX;i++) {
    if (i<absLen) { const t=i/absLen; absorb[i]=Math.sin(0.5*Math.PI*t)**2; }
    else if (i>=NX-absLen) { const t=(NX-1-i)/absLen; absorb[i]=Math.sin(0.5*Math.PI*t)**2; }
    else absorb[i]=1.0;
  }
}

function initPsi() {
  psiRe=new Float64Array(NX); psiIm=new Float64Array(NX);
  const x0=-5.0, sigma=1.0/SIGMA_K, sig2=2.0*sigma*sigma;
  for (let i=0;i<NX;i++) {
    const x=xArr[i], env=Math.exp(-(x-x0)*(x-x0)/sig2);
    psiRe[i]=env*Math.cos(K0*x)*absorb[i];
    psiIm[i]=env*Math.sin(K0*x)*absorb[i];
  }
}

function computeNorm() {
  let n=0; for (let i=0;i<NX;i++) n+=(psiRe[i]*psiRe[i]+psiIm[i]*psiIm[i])*DX; return n;
}

function fluxAt(i) {
  const i0=Math.max(i-1,0), i1=Math.min(i+1,NX-1);
  const dRe=(psiRe[i1]-psiRe[i0])/((i1-i0)*DX);
  const dIm=(psiIm[i1]-psiIm[i0])/((i1-i0)*DX);
  return (psiRe[i]*dIm-psiIm[i]*dRe)/MASS;
}

function step() {
  for (let i=0;i<NX;i++) {
    const re=psiRe[i],im=psiIm[i];
    psiRe[i]=(re*vHalfRe[i]-im*vHalfIm[i])*absorb[i];
    psiIm[i]=(re*vHalfIm[i]+im*vHalfRe[i])*absorb[i];
  }
  fft(psiRe,psiIm,false);
  for (let i=0;i<NX;i++) {
    const re=psiRe[i],im=psiIm[i];
    psiRe[i]=re*kPhaseRe[i]-im*kPhaseIm[i];
    psiIm[i]=re*kPhaseIm[i]+im*kPhaseRe[i];
  }
  fft(psiRe,psiIm,true);
  for (let i=0;i<NX;i++) {
    const re=psiRe[i],im=psiIm[i];
    psiRe[i]=(re*vHalfRe[i]-im*vHalfIm[i])*absorb[i];
    psiIm[i]=(re*vHalfIm[i]+im*vHalfRe[i])*absorb[i];
  }
}

/* ── Render ── */
function render() {
  if (!canvas||!ctx) return;
  ctx.fillStyle = _c('--bg-void');
  ctx.fillRect(0,0,simW,simH);

  const plotH  = Math.floor(simH*0.38);
  const groundY = simH - Math.floor(simH*0.18);
  const y1top   = Math.floor(simH*0.04);

  const xBarL = Math.round(simW*(-XDISP_MIN/XDISP_SPAN));
  const xBarR = Math.round(simW*((-XDISP_MIN+barrierW)/XDISP_SPAN));
  const yBarTop = groundY - plotH*(barrierH/BARRIER_MAX);
  const barrierPixH = groundY - yBarTop;

  const prob = new Float32Array(NX);
  for (let i=0;i<NX;i++) prob[i]=psiRe[i]*psiRe[i]+psiIm[i]*psiIm[i];

  function toX(i)  { return ((xArr[i]-XDISP_MIN)/XDISP_SPAN)*simW; }
  function toY1(v) { return y1top+plotH*(1-(v+YMAX_PSI)/(2*YMAX_PSI)); }
  function toY2(v) { return groundY-plotH*(v/YMAX_PROB); }
  const iDisp0=Math.ceil((XDISP_MIN-XMIN)/DX);
  const iDisp1=Math.floor((XDISP_MAX-XMIN)/DX);

  /* |ψ|² fill */
  const grad=ctx.createLinearGradient(0,groundY-plotH,0,groundY);
  grad.addColorStop(0,   'rgba(42,190,217,0.85)');
  grad.addColorStop(0.6, 'rgba(42,190,217,0.55)');
  grad.addColorStop(1,   'rgba(42,190,217,0.20)');
  ctx.save();
  ctx.beginPath(); ctx.moveTo(toX(iDisp0),groundY);
  for (let i=iDisp0;i<=iDisp1;i++) ctx.lineTo(toX(i),toY2(prob[i]));
  ctx.lineTo(toX(iDisp1),groundY); ctx.closePath();
  ctx.fillStyle=grad; ctx.fill(); ctx.restore();

  /* |ψ|² outline */
  ctx.save();
  ctx.strokeStyle=_c('--cyan'); ctx.lineWidth=1.5;
  ctx.shadowColor=_c('--cyan'); ctx.shadowBlur=10;
  ctx.beginPath();
  for (let i=iDisp0;i<=iDisp1;i++){const px=toX(i),py=toY2(prob[i]);i===iDisp0?ctx.moveTo(px,py):ctx.lineTo(px,py);}
  ctx.stroke(); ctx.restore();

  /* incident energy line */
  const E0=K0*K0/(2*MASS);
  const yE0=groundY-plotH*(E0/BARRIER_MAX);
  if (yE0>=groundY-plotH&&yE0<=groundY) {
    ctx.save(); ctx.setLineDash([4,6]);
    ctx.strokeStyle='rgba(255,109,162,0.55)'; ctx.lineWidth=2.5;
    ctx.beginPath(); ctx.moveTo(0,yE0); ctx.lineTo(simW,yE0); ctx.stroke(); ctx.restore();
  }

  /* barrier crosshatch */
  const hatchColor='rgba(100,130,150,0.30)';
  const hatchLW=Math.max(1,Math.round(simW*0.003));
  const hatchStep=Math.round(simW*0.035);
  function hatch(rx,ry,rw,rh) {
    if (rw<=0||rh<=0) return;
    ctx.save(); ctx.beginPath(); ctx.rect(rx,ry,rw,rh); ctx.clip();
    ctx.strokeStyle=hatchColor; ctx.lineWidth=hatchLW;
    for (let d=-simH;d<=simW+simH;d+=hatchStep){
      ctx.beginPath();ctx.moveTo(d,groundY);ctx.lineTo(d-simH,groundY+simH);ctx.stroke();
      ctx.beginPath();ctx.moveTo(d,groundY);ctx.lineTo(d+simH,groundY-simH);ctx.stroke();
    }
    ctx.restore();
  }
  ctx.save(); ctx.globalAlpha=0.7;
  hatch(0,groundY,simW,simH-groundY);
  hatch(xBarL,yBarTop,xBarR-xBarL,barrierPixH);
  ctx.strokeStyle='rgba(100,130,150,0.65)'; ctx.lineWidth=hatchLW*1.5;
  ctx.beginPath();
  ctx.moveTo(xBarL,groundY);ctx.lineTo(xBarL,yBarTop);
  ctx.lineTo(xBarR,yBarTop);ctx.lineTo(xBarR,groundY);
  ctx.stroke(); ctx.restore();

  /* zero line */
  ctx.save(); ctx.setLineDash([4,5]);
  ctx.strokeStyle='rgba(168,192,208,0.28)'; ctx.lineWidth=1;
  ctx.beginPath(); ctx.moveTo(0,y1top+plotH/2); ctx.lineTo(simW,y1top+plotH/2); ctx.stroke(); ctx.restore();

  /* ±|ψ| envelope */
  ctx.save(); ctx.strokeStyle='rgba(42,190,217,0.5)'; ctx.lineWidth=2.5;
  ctx.beginPath();
  for (let i=iDisp0;i<=iDisp1;i++){const px=toX(i),py=toY1(Math.sqrt(prob[i]));i===iDisp0?ctx.moveTo(px,py):ctx.lineTo(px,py);}
  ctx.stroke();
  ctx.beginPath();
  for (let i=iDisp0;i<=iDisp1;i++){const px=toX(i),py=toY1(-Math.sqrt(prob[i]));i===iDisp0?ctx.moveTo(px,py):ctx.lineTo(px,py);}
  ctx.stroke(); ctx.restore();

  /* Im(ψ) pink, Re(ψ) teal — double-pass glow */
  function strokeGlow(buildPath,glowCol,dimCol,mainCol) {
    ctx.save(); ctx.shadowColor=glowCol; ctx.shadowBlur=10;
    ctx.strokeStyle=dimCol; ctx.lineWidth=2.5; ctx.globalAlpha=0.5;
    buildPath(); ctx.stroke(); ctx.restore();
    ctx.save(); ctx.shadowColor=glowCol; ctx.shadowBlur=18;
    ctx.strokeStyle=mainCol; ctx.lineWidth=1.8;
    buildPath(); ctx.stroke(); ctx.restore();
  }
  strokeGlow(
    ()=>{ctx.beginPath();for(let i=iDisp0;i<=iDisp1;i++){const px=toX(i),py=toY1(psiIm[i]);i===iDisp0?ctx.moveTo(px,py):ctx.lineTo(px,py);}},
    _c('--pink-light'),_c('--pink-dark'),_c('--pink-light')
  );
  strokeGlow(
    ()=>{ctx.beginPath();for(let i=iDisp0;i<=iDisp1;i++){const px=toX(i),py=toY1(psiRe[i]);i===iDisp0?ctx.moveTo(px,py):ctx.lineTo(px,py);}},
    _c('--teal-light'),_c('--teal-dark'),_c('--teal-light')
  );

  /* prob sum for auto-stop */
  let pTot=0;
  for (let i=iDisp0;i<=iDisp1;i++) pTot+=prob[i]*DX;
  return pTot;
}

function iNearest(x) { return Math.round((x-XMIN)/DX); }

function setFireReady(ready) {
  canFire=ready;
  const btn=document.getElementById('qtun-fire-btn');
  if (btn) { btn.classList.toggle('active',ready); }
  const hSl=document.getElementById('qtun-height');
  const wSl=document.getElementById('qtun-width');
  if (hSl) hSl.disabled=!ready;
  if (wSl) wSl.disabled=!ready;
}

/* ── Animation loop ── */
function loop() {
  if (running) {
    if (firing && !paused) {
      const iL=iNearest(-9.0), iR=iNearest(13.0);
      const nSteps=speedVal>=5?Math.max(1,Math.round((speedVal-4)/4)):0;
      const framesPerStep=speedVal<5?(5-speedVal):1;
      frameSkipCount++;
      const doSteps=speedVal>=5||frameSkipCount>=framesPerStep;
      if (doSteps) {
        if (speedVal<5) frameSkipCount=0;
        const n=speedVal>=5?nSteps:1;
        for (let s=0;s<n;s++) {
          step();
          phiR-=fluxAt(iL)*DT;
          phiT+=fluxAt(iR)*DT;
        }
      }
      const pTot=render();
      if (!probArmed&&pTot>0.5) { probArmed=true; probPeak=pTot; }
      if (probArmed&&pTot>probPeak) probPeak=pTot;
      if (probArmed&&pTot<0.10*probPeak) {
        firing=false; probArmed=false; probPeak=0;
        setFireReady(true);
        const pb=document.getElementById('qtun-pause-btn');
        if (pb) { pb.textContent='Pause'; pb.classList.remove('active'); }
      }
    } else if (!firing) {
      render();
    }
  }
  frameId=requestAnimationFrame(loop);
}

/* ── Shell wiring ── */
const shell = new AppletShell({
  id:     'qtun',
  title:  'Quantum Tunneling',
  gap:    0,
  layout: 'stacked',

  headerBtns: `<button class="applet-shell-header-btn" onclick="qtunReset()">Reset</button><button class="applet-shell-header-btn" id="qtun-pause-btn" onclick="qtunPause()">Pause</button>`,


  ctrlHTML: `
    <div class="applet-shell-ctrl-section">
      <div class="applet-shell-btn-row">
        <button class="applet-shell-btn active" id="qtun-fire-btn" onclick="qtunFire()">Fire</button>
      </div>
    </div>
    <div class="applet-shell-ctrl-section" style="flex:1; min-width:120px;">
      <div class="applet-shell-ctrl-title">Barrier Height</div>
      <div class="applet-shell-slider-row">
        <span class="applet-shell-side">Low</span>
        <input type="range" id="qtun-height" min="0" max="180" step="1" value="60">
        <span class="applet-shell-side">High</span>
      </div>
    </div>
    <div class="applet-shell-ctrl-section" style="flex:1; min-width:120px;">
      <div class="applet-shell-ctrl-title">Barrier Width</div>
      <div class="applet-shell-slider-row">
        <span class="applet-shell-side">Thin</span>
        <input type="range" id="qtun-width" min="0.1" max="3.0" step="0.05" value="0.8">
        <span class="applet-shell-side">Thick</span>
      </div>
    </div>
    <div class="applet-shell-ctrl-section" style="flex:1; min-width:100px;">
      <div class="applet-shell-ctrl-title">Speed</div>
      <div class="applet-shell-slider-row">
        <span class="applet-shell-side">Slow</span>
        <input type="range" id="qtun-speed" min="1" max="40" step="1" value="24">
        <span class="applet-shell-side">Fast</span>
      </div>
    </div>
  `,

  onOpen: function ({ canvas: c, W, H, S }) {
    canvas = c;
    ctx    = canvas.getContext('2d');
    simW   = canvas.width  = W || S;
    simH   = canvas.height = H || S;

    barrierH=60; barrierW=0.8; speedVal=24;
    firing=false; paused=false; canFire=true;
    probArmed=false; probPeak=0; frameSkipCount=0;

    const hSl=document.getElementById('qtun-height');
    const wSl=document.getElementById('qtun-width');
    if (hSl) hSl.value=60;
    if (wSl) wSl.value=0.8;
    document.getElementById('qtun-speed').value=24;

    const fireBtn=document.getElementById('qtun-fire-btn');
    if (fireBtn) { fireBtn.classList.add('active'); fireBtn.disabled=false; }
    const pb=document.getElementById('qtun-pause-btn');
    if (pb) { pb.textContent='Pause'; pb.classList.remove('active'); }

    precompute();
    initPsi();
    p0=computeNorm(); phiR=0; phiT=0;
    render();

    running=true;
    if (!frameId) frameId=requestAnimationFrame(loop);
  },

  onClose: function () {
    running=false;
    firing=false;
    if (frameId) { cancelAnimationFrame(frameId); frameId=null; }
  },

  onResize: function ({ W, H, S }) {
    if (!canvas) return;
    simW = canvas.width  = W || S;
    simH = canvas.height = H || S;
    if (!firing) render();
  },
});

window.qtunOpen  = () => shell.open();
window.qtunClose = () => shell.close();

window.qtunFire = function () {
  if (!canFire) return;
  precompute(); initPsi(); p0=computeNorm(); phiR=0; phiT=0;
  probArmed=false; probPeak=0; paused=false; firing=true;
  setFireReady(false);
  const pb=document.getElementById('qtun-pause-btn');
  if (pb) { pb.textContent='Pause'; pb.classList.remove('active'); }
};

window.qtunPause = function () {
  if (!firing) return;
  paused=!paused;
  const btn=document.getElementById('qtun-pause-btn');
  if (btn) { btn.textContent=paused?'Resume':'Pause'; btn.classList.toggle('active',paused); }
};

window.qtunReset = function () {
  firing=false; paused=false; probArmed=false; probPeak=0; phiR=0; phiT=0;
  setFireReady(true);
  const pb=document.getElementById('qtun-pause-btn');
  if (pb) { pb.textContent='Pause'; pb.classList.remove('active'); }
  precompute(); initPsi(); p0=computeNorm(); render();
};

document.getElementById('qtun-height').addEventListener('input', function () {
  barrierH=parseFloat(this.value);
  if (!firing) { precompute(); initPsi(); p0=computeNorm(); phiR=0; phiT=0; render(); }
});
document.getElementById('qtun-width').addEventListener('input', function () {
  barrierW=parseFloat(this.value);
  if (!firing) { precompute(); initPsi(); p0=computeNorm(); phiR=0; phiT=0; render(); }
});
document.getElementById('qtun-speed').addEventListener('input', function () {
  speedVal=parseInt(this.value); frameSkipCount=0;
});

})();
