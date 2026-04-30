(function() {
'use strict';

const _cs = getComputedStyle(document.documentElement);
const _c   = n => _cs.getPropertyValue(n).trim();
const _rgb  = n => { const h = _c(n).replace('#',''); const v = parseInt(h,16); return [(v>>16)&0xFF,(v>>8)&0xFF,v&0xFF]; };

let leafsInitialized = false;
let leafsRunning = false;
let leafsExtinct = false;
let leafsAnimFrame = null;
let stepsPerFrame = 5;

// ── Geometry (golden ratio) ────────────────────────────────────────────────
const PHI = 1.6180339887;
const LGAP = 10; // px gap between panels

function leafsSetGeometry() {
  const vw = window.innerWidth;
  const vh = window.innerHeight;
  const margin = 24;
  const availW = vw - 2 * margin;
  const availH = vh - 2 * margin;

  // Body row is a golden rectangle: W × (W/φ)
  // Header = 10% of body height; one gap separates header from body
  // Total height = H_body * 1.1 + LGAP  ≤  availH
  const maxWfromH = (availH - LGAP) * PHI / 1.1;
  const W = Math.min(availW, maxWfromH);

  const H_body = W / PHI;
  const H_hdr  = H_body * 0.1;

  const left    = Math.round((vw - W) / 2);
  const topHdr  = Math.round((vh - (H_hdr + LGAP + H_body)) / 2);
  const topBody = topHdr + Math.round(H_hdr) + LGAP;

  // Horizontal split at golden ratio: sim = φ/(φ+1)·W, right = 1/(φ+1)·W
  const W_sim   = Math.round(W * PHI / (PHI + 1));
  const W_right = Math.round(W - W_sim - LGAP);

  const ov = document.getElementById('leafs-overlay');
  ov.style.setProperty('--leafs-left',     left     + 'px');
  ov.style.setProperty('--leafs-top-hdr',  topHdr   + 'px');
  ov.style.setProperty('--leafs-top-body', topBody  + 'px');
  ov.style.setProperty('--leafs-W',        Math.round(W)      + 'px');
  ov.style.setProperty('--leafs-H-hdr',    Math.round(H_hdr)  + 'px');
  ov.style.setProperty('--leafs-H-body',   Math.round(H_body) + 'px');
  ov.style.setProperty('--leafs-W-sim',    W_sim    + 'px');
  ov.style.setProperty('--leafs-W-right',  W_right  + 'px');
  ov.style.setProperty('--leafs-gap',      LGAP     + 'px');
}

// ── Toggle ─────────────────────────────────────────────────────────────────
window.leafsToggle = function() {
  const overlay = document.getElementById('leafs-overlay');
  const panels  = document.querySelectorAll('.leafs-panel');
  const btn     = document.querySelector('.leafs-launch-btn');
  const isOpen  = overlay.classList.contains('leafs-open');

  if (isOpen) {
    panels.forEach(p => p.classList.remove('leafs-open'));
    setTimeout(() => overlay.classList.remove('leafs-open'), 500);
    btn.innerHTML = '&#9654; LAUNCH LEAFS';
    leafsRunning = false;
    cancelAnimationFrame(leafsAnimFrame);
    document.getElementById('leafs-btn-run').textContent = 'RUN';
    document.getElementById('leafs-btn-run').classList.remove('active');
  } else {
    leafsSetGeometry();
    overlay.classList.add('leafs-open');
    const panelIds = ['leafs-header', 'leafs-sim-panel', 'leafs-right-panel'];
    const delays   = [0, 80, 160];
    panelIds.forEach((id, i) => {
      setTimeout(() => document.getElementById(id).classList.add('leafs-open'), delays[i]);
    });
    btn.innerHTML = '&#9646;&#9646; HIDE LEAFS';
    setTimeout(() => {
      if (!leafsInitialized) {
        leafsInit();
        leafsInitialized = true;
      } else {
        leafsOnResize();
      }
    }, 220);
  }
};

// ── Palette ────────────────────────────────────────────────────────────────
const PAL = {
  bg:            _c('--blue-mid'),
  resource:      _c('--green-light'),
  resourceInner: _c('--green-dark'),
  agentBody:     _c('--pink-muted'),
  agentBorder:   _c('--text-light'),
  cone:          _c('--blue-light'),
  coneFired:     _c('--pink-muted'),
  posEdge:       _c('--blue-lighter'),
  negEdge:       _c('--pink-deep'),
  hidNode:       _c('--pink-muted'),
  inNode:        _c('--blue-light'),
  outNode:       _c('--blue-muted'),
};

// ── Parameters ─────────────────────────────────────────────────────────────
const P = {
  Nx:5, Ny:5, L:10, dt:0.025,
  Gamma:0.001, gamma:0.1, epsilon:0.1,
  n_agents:5, n_sensors:3,
  sensor_phi: Math.PI/6,
  sensor_angles: [-Math.PI/3, 0, Math.PI/3],
  v:20.0, sigma:0.1, tau_A:1.0, dphi:Math.PI/8,
  r_col:1.0, R_sense:6.0, mu0:1.0, s_max:1.0,
  repro_threshold:0.8, repro_offset:1.5,
};
P.domainX = P.Nx*P.L; P.domainY = P.Ny*P.L;
P.Neq = P.Gamma*P.L*P.L/(P.epsilon*P.gamma);

// ── Genome ─────────────────────────────────────────────────────────────────
function makeGenome(scale) {
  const rn = () => (Math.random()*2-1)*scale;
  return {
    W_in:  Array.from({length:3}, () => Array.from({length:3}, rn)),
    W_rec: Array.from({length:3}, () => Array.from({length:3}, rn)),
    W_out: Array.from({length:2}, () => Array.from({length:3}, rn)),
    b_in:  Array.from({length:3}, rn),
    b_out: Array.from({length:2}, rn),
    sensor_angles: [...P.sensor_angles],
  };
}
function cloneGenome(gm) {
  return {
    W_in:  gm.W_in.map(r=>[...r]), W_rec: gm.W_rec.map(r=>[...r]),
    W_out: gm.W_out.map(r=>[...r]), b_in: [...gm.b_in],
    b_out: [...gm.b_out], sensor_angles: [...gm.sensor_angles],
  };
}
const PRESETS = {
  random: () => makeGenome(0.5),
  zero: () => ({
    W_in:  Array.from({length:3},()=>[0,0,0]),
    W_rec: Array.from({length:3},()=>[0,0,0]),
    W_out: Array.from({length:2},()=>[0,0,0]),
    b_in:[0,0,0], b_out:[0,0], sensor_angles:[...P.sensor_angles],
  }),
};
let genome = makeGenome(0.4);

// ── Environment ────────────────────────────────────────────────────────────
const NR = P.Nx*P.Ny;
const regionID = Array.from({length:P.Nx},(_,ix)=>Array.from({length:P.Ny},(_,iy)=>ix*P.Ny+iy));
const regionXY = [];
for (let ix=0;ix<P.Nx;ix++) for (let iy=0;iy<P.Ny;iy++) regionXY.push([ix,iy]);
const NBD_OFF = [[0,0],[-1,1],[0,1],[1,1],[-1,0],[1,0],[-1,-1],[0,-1],[1,-1]];
const neighbourhood = regionXY.map(([ix,iy])=>
  NBD_OFF.map(([dx,dy])=>{
    const nx=((ix+dx)%P.Nx+P.Nx)%P.Nx, ny=((iy+dy)%P.Ny+P.Ny)%P.Ny;
    return regionID[nx][ny];
  })
);
let regions=[], GammaField;

function gammaField() {
  const patches=[[0.2,0.2],[0.7,0.3],[0.4,0.7],[0.8,0.8],[0.1,0.6]];
  const G=Array.from({length:P.Nx},()=>new Array(P.Ny).fill(0));
  for (let ix=0;ix<P.Nx;ix++) for (let iy=0;iy<P.Ny;iy++) {
    const fx=(ix+0.5)/P.Nx, fy=(iy+0.5)/P.Ny;
    let patch=0;
    for (const [px,py] of patches) patch+=Math.exp(-((fx-px)**2+(fy-py)**2)/0.04);
    G[ix][iy]=P.Gamma*0.5+patch*0.5;
  }
  let sum=0;
  for (let ix=0;ix<P.Nx;ix++) for (let iy=0;iy<P.Ny;iy++) sum+=G[ix][iy];
  const mean=sum/NR;
  for (let ix=0;ix<P.Nx;ix++) for (let iy=0;iy<P.Ny;iy++) G[ix][iy]=G[ix][iy]/mean*P.Gamma;
  return G;
}
function regionLambda(ix,iy){return GammaField[ix][iy]*P.L*P.L*P.dt/P.epsilon;}
function regionPdecay(){return 1-Math.exp(-P.gamma*P.dt);}
function regionNeq(ix,iy){return GammaField[ix][iy]*P.L*P.L/(P.epsilon*P.gamma);}
function poissonSample(lam){
  if(lam<=0)return 0; const L=Math.exp(-lam); let k=0,p=1;
  do{k++;p*=Math.random();}while(p>L); return k-1;
}
function binomialSample(n,p){
  if(n<=0||p<=0)return 0; if(p>=1)return n;
  let k=0; for(let i=0;i<n;i++) if(Math.random()<p)k++; return k;
}
function initRegions(){
  regions=regionXY.map(([ix,iy])=>{
    const neq=regionNeq(ix,iy),N=poissonSample(neq);
    return{N,pos:Array.from({length:N},()=>[Math.random(),Math.random()])};
  });
}
function updateEnvironment(){
  const p=regionPdecay();
  for(let r=0;r<NR;r++){
    const [ix,iy]=regionXY[r];
    const kGrow=poissonSample(regionLambda(ix,iy)),kDecay=binomialSample(regions[r].N,p);
    if(Math.random()>0.5){growRegion(r,kGrow);decayRegion(r,kDecay);}
    else{decayRegion(r,kDecay);growRegion(r,kGrow);}
  }
}
function growRegion(r,k){for(let i=0;i<k;i++)regions[r].pos.push([Math.random(),Math.random()]);regions[r].N+=k;}
function decayRegion(r,k){
  k=Math.min(k,regions[r].N);
  for(let i=0;i<k;i++){
    const idx=Math.floor(Math.random()*regions[r].N);
    regions[r].pos[idx]=regions[r].pos[regions[r].N-1];regions[r].pos.pop();regions[r].N--;
  }
}
function harvestResource(r,idx){
  if(idx<0||idx>=regions[r].N)return false;
  regions[r].pos[idx]=regions[r].pos[regions[r].N-1];regions[r].pos.pop();regions[r].N--;return true;
}
function getResourcesInNeighbourhood(px,py){
  const ix=Math.min(P.Nx-1,Math.floor(px/P.L)),iy=Math.min(P.Ny-1,Math.floor(py/P.L));
  const normX=(px%P.L)/P.L,normY=(py%P.L)/P.L,r0=regionID[ix][iy],results=[];
  for(let nb=0;nb<9;nb++){
    const nid=neighbourhood[r0][nb],[dx,dy]=NBD_OFF[nb],reg=regions[nid];
    for(let i=0;i<reg.N;i++)
      results.push({disp:[(reg.pos[i][0]+dx-normX)*P.L,(reg.pos[i][1]+dy-normY)*P.L],r:nid,idx:i});
  }
  return results;
}

// ── Agent ──────────────────────────────────────────────────────────────────
function makeAgent(gm,pos,theta,health){
  return{
    pos:pos||[Math.random()*P.domainX,Math.random()*P.domainY],
    theta:theta!==undefined?theta:Math.random()*2*Math.PI,
    health:health!==undefined?health:P.s_max*(0.3+Math.random()*0.7),
    h:[0,0,0],genome:gm,sensorFired:[false,false,false],
  };
}
function softmax3(a,b,c){
  const m=Math.max(a,b,c),ea=Math.exp(a-m),eb=Math.exp(b-m),ec=Math.exp(c-m),s=ea+eb+ec;
  return[ea/s,eb/s,ec/s];
}
function randnJS(){return Math.sqrt(-2*Math.log(Math.random()))*Math.cos(2*Math.PI*Math.random());}

function agentStep(agent){
  const gm=agent.genome,resources=getResourcesInNeighbourhood(agent.pos[0],agent.pos[1]);
  const s=[0,0,0];
  for(const{disp}of resources){
    const d=Math.hypot(disp[0],disp[1]);
    if(d<=0||d>P.R_sense)continue;
    const ux=disp[0]/d,uy=disp[1]/d;
    for(let si=0;si<3;si++){
      const ang=agent.theta+gm.sensor_angles[si];
      if(ux*Math.cos(ang)+uy*Math.sin(ang)>=Math.cos(P.sensor_phi))s[si]=1;
    }
  }
  agent.sensorFired=[s[0]===1,s[1]===1,s[2]===1];
  const hPrev=[...agent.h];
  for(let hi=0;hi<3;hi++){
    let v=gm.b_in[hi];
    for(let xi=0;xi<3;xi++)v+=gm.W_in[hi][xi]*s[xi];
    for(let hj=0;hj<3;hj++)v+=gm.W_rec[hi][hj]*hPrev[hj];
    agent.h[hi]=Math.tanh(v);
  }
  let yL=gm.b_out[0],yR=gm.b_out[1];
  for(let hi=0;hi<3;hi++){yL+=gm.W_out[0][hi]*agent.h[hi];yR+=gm.W_out[1][hi]*agent.h[hi];}
  const[pL,,pR]=softmax3(yL,0,yR),rnd=Math.random();
  const action=rnd<pL?-1:rnd<pL+pR?1:0;
  agent.theta+=action*P.dphi+(P.sigma/Math.sqrt(P.tau_A))*randnJS()*Math.sqrt(P.dt);
  agent.pos[0]=((agent.pos[0]+P.v*Math.cos(agent.theta)*P.dt)%P.domainX+P.domainX)%P.domainX;
  agent.pos[1]=((agent.pos[1]+P.v*Math.sin(agent.theta)*P.dt)%P.domainY+P.domainY)%P.domainY;
  let hHarvest=0;
  if(agent.health<P.s_max){
    const cands=resources.filter(({disp})=>disp[0]**2+disp[1]**2<=P.r_col**2);
    cands.sort((a,b)=>a.disp[0]**2+a.disp[1]**2-b.disp[0]**2-b.disp[1]**2);
    for(const{r,idx}of cands)if(harvestResource(r,idx))hHarvest+=P.epsilon;
  }
  agent.health=Math.min(P.s_max,agent.health-P.mu0*P.dt+hHarvest);
}

// ── Simulation state ───────────────────────────────────────────────────────
let agents=[], simTime=0;

function leafsInitSim(){
  simTime=0; leafsExtinct=false;
  document.getElementById('leafs-extinction-stamp').style.display='none';
  initRegions();
  agents=Array.from({length:P.n_agents},()=>makeAgent(genome));
}
function simStep(){
  if(leafsExtinct)return;
  updateEnvironment();
  const newborns=[];
  for(let i=agents.length-1;i>=0;i--){
    const agent=agents[i]; agentStep(agent);
    if(agent.health>=P.repro_threshold){
      const angle=Math.random()*2*Math.PI;
      const ox=agent.pos[0]+Math.cos(angle)*P.repro_offset,oy=agent.pos[1]+Math.sin(angle)*P.repro_offset;
      const dpos=[((ox%P.domainX)+P.domainX)%P.domainX,((oy%P.domainY)+P.domainY)%P.domainY];
      agent.health/=2;
      newborns.push(makeAgent(cloneGenome(agent.genome),dpos,Math.random()*2*Math.PI,agent.health));
    }
    if(agent.health<=0)agents.splice(i,1);
  }
  agents.push(...newborns);
  if(agents.length===0){
    leafsExtinct=true; leafsRunning=false;
    document.getElementById('leafs-btn-run').textContent='RUN';
    document.getElementById('leafs-btn-run').classList.remove('active');
    document.getElementById('leafs-extinction-stamp').style.display='block';
    cancelAnimationFrame(leafsAnimFrame);
  }
  simTime++;
}

// ── Canvas drawing ─────────────────────────────────────────────────────────
const simCanvas = document.getElementById('leafs-sim-canvas');
const simCtx    = simCanvas.getContext('2d');
const genCanvas = document.getElementById('leafs-genome-canvas');
const genCtx    = genCanvas.getContext('2d');

function resizeSimCanvas(){
  const p=document.getElementById('leafs-sim-panel');
  const sb=document.getElementById('leafs-status-bar');
  simCanvas.width=p.clientWidth;
  simCanvas.height=p.clientHeight-sb.offsetHeight;
}
function drawSim(){
  const W=simCanvas.width,H=simCanvas.height;
  const scale=Math.min(W/P.domainX,H/P.domainY);
  const offX=(W-P.domainX*scale)/2,offY=(H-P.domainY*scale)/2;
  const tx=x=>offX+x*scale,ty=y=>offY+y*scale;
  simCtx.fillStyle=PAL.bg; simCtx.fillRect(0,0,W,H);
  for(let r=0;r<NR;r++){
    const[ix,iy]=regionXY[r],reg=regions[r];
    for(let i=0;i<reg.N;i++){
      const px=tx((ix+reg.pos[i][0])*P.L),py=ty((iy+reg.pos[i][1])*P.L);
      simCtx.beginPath();simCtx.arc(px,py,3.5,0,2*Math.PI);
      simCtx.fillStyle=PAL.resource;simCtx.fill();
      simCtx.beginPath();simCtx.arc(px,py,1.8,0,2*Math.PI);
      simCtx.fillStyle=PAL.resourceInner;simCtx.fill();
    }
  }
  for(const agent of agents){
    const ax=tx(agent.pos[0]),ay=ty(agent.pos[1]);
    const hf=Math.max(0,agent.health/P.s_max),agR=scale*1.2*(0.35+0.65*hf);
    const coneR=P.R_sense*scale*0.55;
    for(let si=0;si<3;si++){
      const ang=agent.theta+agent.genome.sensor_angles[si];
      simCtx.beginPath();simCtx.moveTo(ax,ay);
      simCtx.arc(ax,ay,coneR,ang-P.sensor_phi,ang+P.sensor_phi);simCtx.closePath();
      simCtx.fillStyle=agent.sensorFired[si]?PAL.coneFired+'a0':PAL.cone+'50';simCtx.fill();
    }
    simCtx.beginPath();simCtx.arc(ax,ay,agR,0,2*Math.PI);
    simCtx.fillStyle=PAL.agentBody;simCtx.fill();
    simCtx.strokeStyle=PAL.agentBorder;simCtx.lineWidth=1.5;simCtx.stroke();
    simCtx.beginPath();simCtx.moveTo(ax,ay);
    simCtx.lineTo(ax+Math.cos(agent.theta)*agR*1.4,ay+Math.sin(agent.theta)*agR*1.4);
    simCtx.strokeStyle=PAL.agentBorder;simCtx.lineWidth=1.5;simCtx.stroke();
  }
}

// ── Genome graph ───────────────────────────────────────────────────────────
let activeSliderKey=null,genomeHits=[];

function resizeGenomeCanvas(){
  genCanvas.width=genCanvas.offsetWidth;
  genCanvas.height=genCanvas.offsetHeight;
}
function getNodePositions(W,H){
  const ySensors=[0.25,0.50,0.75],arcCX=0.6,arcR=0.4,yMid=0.50;
  const inNodes=ySensors.map(y=>{const xOff=arcCX-Math.sqrt(Math.max(arcR**2-(y-yMid)**2,0));return[xOff*W,y*H];});
  const hidNodes=[[0.54*W,0.28*H],[0.54*W,0.50*H],[0.54*W,0.72*H]];
  const outNodes=[[0.82*W,0.40*H],[0.82*W,0.60*H]];
  return{inNodes,hidNodes,outNodes};
}
function hexToRgb(hex){return[parseInt(hex.slice(1,3),16),parseInt(hex.slice(3,5),16),parseInt(hex.slice(5,7),16)];}
function edgeCol(w,alpha){const[r,g,b]=hexToRgb(w>=0?PAL.posEdge:PAL.negEdge);return`rgba(${r},${g},${b},${alpha})`;}
function drawArrow(ctx,x1,y1,x2,y2,lw,col,NR){
  const dx=x2-x1,dy=y2-y1,d=Math.hypot(dx,dy);if(d<1e-6)return;
  const ux=dx/d,uy=dy/d,sx=x1+ux*NR,sy=y1+uy*NR,ex=x2-ux*NR,ey=y2-uy*NR;
  ctx.beginPath();ctx.moveTo(sx,sy);ctx.lineTo(ex,ey);ctx.strokeStyle=col;ctx.lineWidth=lw;ctx.stroke();
  const as=Math.max(4,lw*1.5),ang=Math.atan2(ey-sy,ex-sx);
  ctx.beginPath();ctx.moveTo(ex,ey);
  ctx.lineTo(ex-as*Math.cos(ang-0.4),ey-as*Math.sin(ang-0.4));
  ctx.lineTo(ex-as*Math.cos(ang+0.4),ey-as*Math.sin(ang+0.4));
  ctx.closePath();ctx.fillStyle=col;ctx.fill();
}
function drawBezierArc(ctx,x1,y1,x2,y2,ctrlOff,lw,col){
  const mx=(x1+x2)/2,my=(y1+y2)/2,dx=x2-x1,dy=y2-y1,d=Math.hypot(dx,dy);
  const px=-dy/Math.max(d,1e-6),py=dx/Math.max(d,1e-6);
  ctx.beginPath();ctx.moveTo(x1,y1);
  ctx.quadraticCurveTo(mx+ctrlOff*px,my+ctrlOff*py,x2,y2);
  ctx.strokeStyle=col;ctx.lineWidth=lw;ctx.stroke();
}
function drawSelfLoop(ctx,cx,cy,lw,col,sign){
  const loopCY=cy+sign*24;ctx.beginPath();ctx.ellipse(cx,loopCY,8,13,0,0,2*Math.PI);
  ctx.strokeStyle=col;ctx.lineWidth=lw;ctx.stroke();
}
function drawGenomeGraph(){
  const W=genCanvas.width,H=genCanvas.height;
  genCtx.clearRect(0,0,W,H);genCtx.fillStyle=_c('--bg-panel');genCtx.fillRect(0,0,W,H);
  const{inNodes,hidNodes,outNodes}=getNodePositions(W,H),gm=genome;
  const allW=[...gm.W_in.flat(),...gm.W_rec.flat(),...gm.W_out.flat(),...gm.b_in,...gm.b_out];
  const wMax=Math.max(...allW.map(Math.abs),1e-6);
  const MIN_LW=0.8,MAX_LW=5.5,W_THRESH=0.02,NODE_R=20,ALPHA=0.75;
  const lw=w=>MIN_LW+(MAX_LW-MIN_LW)*Math.abs(w)/wMax;
  genomeHits=[];
  const isActive=key=>key===activeSliderKey;

  // White glow helpers — set before drawing, clear after
  function glowOn()  { genCtx.shadowColor=_c('--white'); genCtx.shadowBlur=18; }
  function glowOff() { genCtx.shadowColor='transparent'; genCtx.shadowBlur=0; }

  for(let hi=0;hi<3;hi++)for(let ii=0;ii<3;ii++){
    const w=gm.W_in[hi][ii],key=`W_in_${hi}_${ii}`;
    if(Math.abs(w)<W_THRESH&&!isActive(key))continue;
    const[x1,y1]=inNodes[ii],[x2,y2]=hidNodes[hi];
    if(isActive(key))glowOn();
    drawArrow(genCtx,x1,y1,x2,y2,lw(w),edgeCol(w,isActive(key)?1.0:ALPHA),NODE_R);
    glowOff();
    genomeHits.push({x:(x1+x2)/2,y:(y1+y2)/2,r:10,key});
  }
  for(let oi=0;oi<2;oi++)for(let hi=0;hi<3;hi++){
    const w=gm.W_out[oi][hi],key=`W_out_${oi}_${hi}`;
    if(Math.abs(w)<W_THRESH&&!isActive(key))continue;
    const[x1,y1]=hidNodes[hi],[x2,y2]=outNodes[oi];
    if(isActive(key))glowOn();
    drawArrow(genCtx,x1,y1,x2,y2,lw(w),edgeCol(w,isActive(key)?1.0:ALPHA),NODE_R);
    glowOff();
    genomeHits.push({x:(x1+x2)/2,y:(y1+y2)/2,r:10,key});
  }
  const selfSign=[1,1,-1];
  for(let hi=0;hi<3;hi++)for(let hj=0;hj<3;hj++){
    const w=gm.W_rec[hi][hj],key=`W_rec_${hi}_${hj}`;
    if(Math.abs(w)<W_THRESH&&!isActive(key))continue;
    const col=edgeCol(w,isActive(key)?1.0:ALPHA),[x2,y2]=hidNodes[hi],[x1,y1]=hidNodes[hj];
    if(isActive(key))glowOn();
    if(hi===hj){drawSelfLoop(genCtx,x2,y2,lw(w),col,selfSign[hi]);genomeHits.push({x:x2+20,y:y2+selfSign[hi]*24,r:12,key});}
    else{drawBezierArc(genCtx,x1,y1,x2,y2,28*Math.sign(hi-hj),lw(w),col);genomeHits.push({x:(x1+x2)/2+10*Math.sign(hi-hj),y:(y1+y2)/2,r:10,key});}
    glowOff();
  }
  const biasLW=b=>1.0+4.5*Math.abs(b)/wMax,biasCol=b=>b>=0?PAL.posEdge:PAL.negEdge;
  inNodes.forEach(([x,y],i)=>{
    const key=`S${i+1}`;
    genCtx.beginPath();genCtx.arc(x,y,NODE_R,0,2*Math.PI);
    genCtx.fillStyle=PAL.inNode;genCtx.fill();
    genCtx.strokeStyle=_c('--blue');genCtx.lineWidth=1.5;genCtx.stroke();
    // S with subscript inside node
    genCtx.fillStyle=_c('--bg-panel');
    genCtx.font='bold 17px serif';genCtx.textAlign='center';
    genCtx.fillText('S',x-4,y+6);
    genCtx.font='bold 12px serif';
    genCtx.fillText(`${i+1}`,x+6,y+12);
  });
  hidNodes.forEach(([x,y],i)=>{
    const key=`b_in_${i}`,hl=isActive(key),b=gm.b_in[i];
    if(hl)glowOn();
    genCtx.beginPath();genCtx.arc(x,y,NODE_R,0,2*Math.PI);
    genCtx.fillStyle=PAL.hidNode;genCtx.fill();
    genCtx.strokeStyle=hl?_c('--white'):biasCol(b);genCtx.lineWidth=hl?3.5:biasLW(b);genCtx.stroke();
    glowOff();
    // H with subscript inside node
    genCtx.fillStyle=_c('--bg-panel');
    genCtx.font='bold 17px serif';genCtx.textAlign='center';
    genCtx.fillText('H',x-4,y+6);
    genCtx.font='bold 12px serif';
    genCtx.fillText(`${i+1}`,x+6,y+12);
    genomeHits.push({x,y,r:NODE_R,key});
  });
  outNodes.forEach(([x,y],i)=>{
    const key=`b_out_${i}`,hl=isActive(key),b=gm.b_out[i];
    if(hl)glowOn();
    genCtx.beginPath();genCtx.arc(x,y,NODE_R,0,2*Math.PI);
    genCtx.fillStyle=PAL.outNode;genCtx.fill();
    genCtx.strokeStyle=hl?_c('--white'):biasCol(b);genCtx.lineWidth=hl?3.5:biasLW(b);genCtx.stroke();
    glowOff();
    // A with subscript inside node
    genCtx.fillStyle=_c('--bg-panel');
    genCtx.font='bold 17px serif';genCtx.textAlign='center';
    genCtx.fillText('A',x-4,y+6);
    genCtx.font='bold 12px serif';
    genCtx.fillText(`${i+1}`,x+6,y+12);
    genomeHits.push({x,y,r:NODE_R,key});
  });
}
genCanvas.addEventListener('click',e=>{
  const rect=genCanvas.getBoundingClientRect();
  const mx=(e.clientX-rect.left)*(genCanvas.width/rect.width);
  const my=(e.clientY-rect.top)*(genCanvas.height/rect.height);
  for(const h of genomeHits){
    if((mx-h.x)**2+(my-h.y)**2<=h.r**2){
      activeSliderKey=h.key;
      const inp=leafsSliderInputEls[h.key];
      if(inp)inp.scrollIntoView({block:'nearest',behavior:'smooth'});
      drawGenomeGraph();return;
    }
  }
  activeSliderKey=null;leafsClearHighlights();drawGenomeGraph();
});

// ── Slider panel ───────────────────────────────────────────────────────────
function sliderDefs(){
  const gm=genome,defs=[];
  defs.push({section:'S→H'});
  for(let hi=0;hi<3;hi++)for(let ii=0;ii<3;ii++)
    defs.push({key:`W_in_${hi}_${ii}`,label:`W_in[h${hi+1},S${ii+1}]`,get:()=>gm.W_in[hi][ii],set:v=>{gm.W_in[hi][ii]=v;}});
  defs.push({section:'H→H'});
  for(let hi=0;hi<3;hi++)for(let hj=0;hj<3;hj++)
    defs.push({key:`W_rec_${hi}_${hj}`,label:`W_rec[h${hi+1},h${hj+1}]`,get:()=>gm.W_rec[hi][hj],set:v=>{gm.W_rec[hi][hj]=v;}});
  defs.push({section:'H→A'});
  const outNames=['←','→'];
  for(let oi=0;oi<2;oi++)for(let hi=0;hi<3;hi++)
    defs.push({key:`W_out_${oi}_${hi}`,label:`W_out[${outNames[oi]},h${hi+1}]`,get:()=>gm.W_out[oi][hi],set:v=>{gm.W_out[oi][hi]=v;}});
  defs.push({section:'Hidden Bias'});
  for(let hi=0;hi<3;hi++)
    defs.push({key:`b_in_${hi}`,label:`b_in[h${hi+1}]`,get:()=>gm.b_in[hi],set:v=>{gm.b_in[hi]=v;}});
  defs.push({section:'Action Bias'});
  for(let oi=0;oi<2;oi++)
    defs.push({key:`b_out_${oi}`,label:`b_out[${outNames[oi]}]`,get:()=>gm.b_out[oi],set:v=>{gm.b_out[oi]=v;}});
  return defs;
}
const leafsSliderRowEls={},leafsSliderInputEls={};
function leafsBuildSliders(){
  const panel=document.getElementById('leafs-slider-panel');panel.innerHTML='';
  const defs=sliderDefs();
  let currentCluster=null, currentRows=null;
  for(const d of defs){
    if(d.section){
      // New cluster block
      currentCluster=document.createElement('div');
      currentCluster.className='leafs-slider-cluster';
      const labelDiv=document.createElement('div');labelDiv.className='leafs-cluster-label';
      const labelSpan=document.createElement('span');labelSpan.textContent=d.section;
      labelDiv.appendChild(labelSpan);currentCluster.appendChild(labelDiv);
      currentRows=document.createElement('div');currentRows.className='leafs-cluster-rows';
      currentCluster.appendChild(currentRows);panel.appendChild(currentCluster);
      continue;
    }
    const row=document.createElement('div');row.className='leafs-slider-row';row.dataset.key=d.key;
    const inp=document.createElement('input');inp.type='range';inp.min=-1;inp.max=1;inp.step=0.01;inp.value=d.get().toFixed(2);
    inp.addEventListener('input',()=>{const v=parseFloat(inp.value);d.set(v);for(const a of agents)a.genome=genome;drawGenomeGraph();});
    inp.addEventListener('mousedown',()=>{activeSliderKey=d.key;drawGenomeGraph();});
    inp.addEventListener('touchstart',()=>{activeSliderKey=d.key;drawGenomeGraph();});
    inp.addEventListener('mouseup',()=>{activeSliderKey=null;drawGenomeGraph();});
    inp.addEventListener('touchend',()=>{activeSliderKey=null;drawGenomeGraph();});
    row.appendChild(inp);currentRows.appendChild(row);
    leafsSliderRowEls[d.key]=row;leafsSliderInputEls[d.key]=inp;
  }
}
function leafsRefreshSliders(){
  for(const d of sliderDefs()){
    if(d.section)continue;
    if(leafsSliderInputEls[d.key])leafsSliderInputEls[d.key].value=d.get().toFixed(2);
  }
}

// ── Animation loop ─────────────────────────────────────────────────────────
function loop(){
  if(!leafsRunning)return;
  for(let i=0;i<stepsPerFrame;i++)simStep();
  drawSim();drawGenomeGraph();
  const totalRes=regions.reduce((s,r)=>s+r.N,0);
  document.getElementById('leafs-status-bar').textContent=
    `t = ${(simTime*P.dt).toFixed(2)}  |  agents: ${agents.length}  |  resources: ${totalRes}`;
  leafsAnimFrame=requestAnimationFrame(loop);
}

// ── Controls ───────────────────────────────────────────────────────────────
document.getElementById('leafs-btn-run').addEventListener('click',()=>{
  if(leafsExtinct)return;
  leafsRunning=!leafsRunning;
  document.getElementById('leafs-btn-run').textContent=leafsRunning?'PAUSE':'RUN';
  document.getElementById('leafs-btn-run').classList.toggle('active',leafsRunning);
  if(leafsRunning)loop();
});
document.getElementById('leafs-btn-reset').addEventListener('click',()=>{
  leafsRunning=false;leafsExtinct=false;cancelAnimationFrame(leafsAnimFrame);
  document.getElementById('leafs-btn-run').textContent='RUN';
  document.getElementById('leafs-btn-run').classList.remove('active');
  leafsInitSim();drawSim();drawGenomeGraph();leafsRefreshSliders();
  document.getElementById('leafs-status-bar').textContent='t = 0  |  agents: 5  |  resources: —';
});
document.getElementById('leafs-speed-slider').addEventListener('input',function(){
  stepsPerFrame=parseInt(this.value);
  document.getElementById('leafs-speed-val').textContent=stepsPerFrame;
});
document.getElementById('leafs-preset-select').addEventListener('change',function(){
  if(PRESETS[this.value]){genome=PRESETS[this.value]();for(const a of agents)a.genome=genome;leafsBuildSliders();drawGenomeGraph();}
});

// ── Resize ─────────────────────────────────────────────────────────────────
function leafsOnResize(){
  const ov = document.getElementById('leafs-overlay');
  if (ov.classList.contains('leafs-open')) leafsSetGeometry();
  resizeSimCanvas();resizeGenomeCanvas();drawSim();drawGenomeGraph();
}
window.addEventListener('resize',leafsOnResize);

// ── Init ───────────────────────────────────────────────────────────────────
function leafsInit(){
  GammaField=gammaField();
  resizeSimCanvas();leafsBuildSliders();resizeGenomeCanvas();
  leafsInitSim();drawSim();drawGenomeGraph();
}

})();
