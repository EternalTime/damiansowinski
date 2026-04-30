(function () {
'use strict';

const _cs  = getComputedStyle(document.documentElement);
const _c   = n => _cs.getPropertyValue(n).trim();
const _rgb  = n => { const h = _c(n).replace('#',''); const v = parseInt(h,16); return [(v>>16)&0xFF,(v>>8)&0xFF,v&0xFF]; };
const _rgba = (n, a) => { const [r,g,b] = _rgb(n); return `rgba(${r},${g},${b},${a})`; };

// ── Mesh parameters (will iterate with user) ──────────────────────────────────
const NR = 4;   // elements in radial direction
const NZ = 12;  // elements in axial direction

// ── Physical parameters ───────────────────────────────────────────────────────
const R0 = 1.0;   // cylinder radius (normalised)
const H0 = 3.0;   // cylinder height (3:1 aspect ratio)

let E  = 1.0;     // Young's modulus (normalised)
let nu = 0.3;     // Poisson ratio
let delta = 0.0;  // prescribed axial displacement at top (+ = tension, - = compression)
const deltaMax = 0.75 * H0;  // deformed height ranges from 25% to 175% of H0

let stressComponent = 0; // 0=von Mises, 1=σ_r, 2=σ_z, 3=σ_θ, 4=τ_rz
let nonlinear = false;   // false=linear SRI, true=Neo-Hookean NR

// ── FEM state ─────────────────────────────────────────────────────────────────
let nodes = null;       // [nNodes][2] — r,z coords
let elems = null;       // [nElems][4] — node indices (CCW: bl, br, tr, tl)
let nNodes, nElems;
let K_fact = null;      // factorised stiffness (stored as banded Cholesky)
let u = null;           // displacement vector [2*nNodes]
let stressNodes = null; // [nNodes][5] — sr,sz,st,trz,vm at nodes

// ── Three.js state ────────────────────────────────────────────────────────────
let renderer, scene, camera;
let cylMesh = null;
let simCanvas;
const BASE_RADIUS = 6.0;
let orbit = { dragging: false, lastX: 0, lastY: 0, theta: 0.6 + Math.PI, phi: 1.1, radius: BASE_RADIUS };

let running = false, frameId = null;
let needsSolve = true;
let zCenter = 0, H_def = H0; // updated each solve

// ── Palette for stress coloring ───────────────────────────────────────────────
// teal-light → cyan → pink-dark → pink-light  (blue=low, pink=high)
function buildPalette() {
  const stops = [
    { u: 0.00, c: _rgb('--teal-light') },
    { u: 0.35, c: _rgb('--cyan') },
    { u: 0.65, c: _rgb('--pink-dark') },
    { u: 1.00, c: _rgb('--pink-light') },
  ];
  return stops;
}
let PALETTE = null;

function stressToColor(t) {
  // t in [0,1]
  t = Math.max(0, Math.min(1, t));
  const stops = PALETTE;
  let i = 0;
  while (i < stops.length - 2 && stops[i+1].u <= t) i++;
  const s0 = stops[i], s1 = stops[i+1];
  const f = (t - s0.u) / (s1.u - s0.u + 1e-12);
  const r = s0.c[0] + f*(s1.c[0]-s0.c[0]);
  const g = s0.c[1] + f*(s1.c[1]-s0.c[1]);
  const b = s0.c[2] + f*(s1.c[2]-s0.c[2]);
  return new THREE.Color(r/255, g/255, b/255);
}

// ── Build mesh ────────────────────────────────────────────────────────────────
function buildMesh() {
  const nrn = NR + 1, nzn = NZ + 1;
  nNodes = nrn * nzn;
  nElems = NR * NZ;

  nodes = new Float64Array(nNodes * 2);
  elems = new Int32Array(nElems * 4);

  for (let j = 0; j < nzn; j++) {
    for (let i = 0; i < nrn; i++) {
      const idx = j * nrn + i;
      nodes[idx*2]   = (i / NR) * R0;        // r
      nodes[idx*2+1] = (j / NZ) * H0;        // z
    }
  }

  let e = 0;
  for (let j = 0; j < NZ; j++) {
    for (let i = 0; i < NR; i++) {
      const bl = j*nrn + i,  br = j*nrn + i+1;
      const tl = (j+1)*nrn + i, tr = (j+1)*nrn + i+1;
      elems[e*4]   = bl;
      elems[e*4+1] = br;
      elems[e*4+2] = tr;
      elems[e*4+3] = tl;
      e++;
    }
  }
}

// ── Gauss quadrature points and weights ──────────────────────────────────────
const GP  = [-1/Math.sqrt(3), 1/Math.sqrt(3)]; // 2-point
const GW  = [1.0, 1.0];
const GP1 = [0.0];                              // 1-point (reduced)
const GW1 = [2.0];

// Shape functions for Q4 in [-1,1]^2
function shape(xi, eta) {
  return [
    0.25*(1-xi)*(1-eta),
    0.25*(1+xi)*(1-eta),
    0.25*(1+xi)*(1+eta),
    0.25*(1-xi)*(1+eta),
  ];
}
function shapeDeriv(xi, eta) {
  // dN/dxi, dN/deta  — rows: 0=dxi, 1=deta
  return [
    [-(1-eta)*0.25,  (1-eta)*0.25,  (1+eta)*0.25, -(1+eta)*0.25],
    [-(1-xi)*0.25,  -(1+xi)*0.25,  (1+xi)*0.25,   (1-xi)*0.25 ],
  ];
}

// ── Assemble global stiffness matrix — Selective Reduced Integration ──────────
// Deviatoric part: 2×2 Gauss (full integration)
// Volumetric part: 1×1 Gauss (reduced integration) — eliminates volumetric locking
// DOFs: node i → ur = 2i, uz = 2i+1
function assembleK() {
  const ndof = nNodes * 2;
  const Kg = new Float64Array(ndof * ndof);

  // Bulk modulus κ and shear modulus G
  const kappa = E / (3*(1-2*nu));   // volumetric
  const G     = E / (2*(1+nu));     // shear/deviatoric

  // Deviatoric D = D_full - D_vol
  // D_vol = κ * m*mᵀ  where m=[1,1,1,0]
  // D_full_ij for axisymmetric:
  function getDdev() {
    const c = E / ((1+nu)*(1-2*nu));
    const Df = [
      c*(1-nu), c*nu,    c*nu,    0,
      c*nu,     c*(1-nu),c*nu,    0,
      c*nu,     c*nu,    c*(1-nu),0,
      0,        0,       0,       c*(1-2*nu)/2,
    ];
    // subtract volumetric part κ*m*mᵀ
    const Dd = Df.slice();
    const mi = [0,1,2]; // indices where m=1
    for (const i of mi) for (const j of mi) Dd[i*4+j] -= kappa;
    return Dd;
  }

  const Ddev = getDdev();

  // Helper: compute B matrix and Jacobian det at a gauss point
  function computeB(xi, eta, rc, zc) {
    const N  = shape(xi, eta);
    const dN = shapeDeriv(xi, eta);
    let J00=0,J01=0,J10=0,J11=0;
    for (let k=0;k<4;k++){J00+=dN[0][k]*rc[k];J01+=dN[0][k]*zc[k];J10+=dN[1][k]*rc[k];J11+=dN[1][k]*zc[k];}
    const detJ = J00*J11 - J01*J10;
    const invJ00=J11/detJ, invJ01=-J01/detJ, invJ10=-J10/detJ, invJ11=J00/detJ;
    const dNr = new Float64Array(4), dNz = new Float64Array(4);
    for (let k=0;k<4;k++){
      dNr[k] = invJ00*dN[0][k] + invJ10*dN[1][k];
      dNz[k] = invJ01*dN[0][k] + invJ11*dN[1][k];
    }
    let rg = 0; for(let k=0;k<4;k++) rg += N[k]*rc[k];
    const B = new Float64Array(4*8);
    for (let k=0;k<4;k++){
      B[0*8 + 2*k]   = dNr[k];
      B[1*8 + 2*k+1] = dNz[k];
      B[2*8 + 2*k]   = N[k]/(rg+1e-14);
      B[3*8 + 2*k]   = dNz[k];
      B[3*8 + 2*k+1] = dNr[k];
    }
    return { B, detJ, rg };
  }

  // Helper: accumulate ke += fac * Bᵀ D B
  function addBtDB(ke, B, D, fac) {
    const DB = new Float64Array(4*8);
    for (let row=0;row<4;row++)
      for (let col=0;col<8;col++){
        let s=0; for(let m=0;m<4;m++) s+=D[row*4+m]*B[m*8+col];
        DB[row*8+col]=s;
      }
    for (let r2=0;r2<8;r2++)
      for (let c2=0;c2<8;c2++){
        let s=0; for(let m=0;m<4;m++) s+=B[m*8+r2]*DB[m*8+c2];
        ke[r2*8+c2]+=fac*s;
      }
  }

  for (let e = 0; e < nElems; e++) {
    const n0=elems[e*4], n1=elems[e*4+1], n2=elems[e*4+2], n3=elems[e*4+3];
    const econn = [n0,n1,n2,n3];
    const rc = [nodes[n0*2],nodes[n1*2],nodes[n2*2],nodes[n3*2]];
    const zc = [nodes[n0*2+1],nodes[n1*2+1],nodes[n2*2+1],nodes[n3*2+1]];
    const ke = new Float64Array(8*8);

    // ── Deviatoric part: 2×2 Gauss ────────────────────────────────────────
    for (let gi=0;gi<2;gi++) for (let gj=0;gj<2;gj++) {
      const { B, detJ, rg } = computeB(GP[gi], GP[gj], rc, zc);
      addBtDB(ke, B, Ddev, GW[gi]*GW[gj]*detJ*2*Math.PI*rg);
    }

    // ── Volumetric part: 1×1 Gauss ────────────────────────────────────────
    // D_vol = κ * m*mᵀ; Bᵀ D_vol B = κ * (m·B)ᵀ (m·B)
    {
      const { B, detJ, rg } = computeB(GP1[0], GP1[0], rc, zc);
      const fac = GW1[0]*GW1[0]*detJ*2*Math.PI*rg*kappa;
      // mB[col] = sum_i m[i]*B[i*8+col]  (m=[1,1,1,0])
      const mB = new Float64Array(8);
      for (let col=0;col<8;col++) mB[col] = B[0*8+col]+B[1*8+col]+B[2*8+col];
      for (let r2=0;r2<8;r2++)
        for (let c2=0;c2<8;c2++)
          ke[r2*8+c2] += fac*mB[r2]*mB[c2];
    }

    // Scatter ke into Kg
    const dofs = [econn[0]*2,econn[0]*2+1,econn[1]*2,econn[1]*2+1,
                  econn[2]*2,econn[2]*2+1,econn[3]*2,econn[3]*2+1];
    for (let r2=0;r2<8;r2++)
      for (let c2=0;c2<8;c2++)
        Kg[dofs[r2]*ndof + dofs[c2]] += ke[r2*8+c2];
  }
  return Kg;
}

// ── Apply BCs and solve ───────────────────────────────────────────────────────
// BCs: bottom nodes (j=0): uz=0, ur free (sliding)
//      top nodes (j=NZ): uz=delta, ur free (sliding)
//      axis nodes (i=0): ur=0 (symmetry)
function solveLinear() {
  const ndof = nNodes * 2;
  const Kg = assembleK();

  // Identify constrained DOFs
  const nrn = NR + 1;
  const constrained = new Map(); // dof -> prescribed value

  // Bottom face: sliding — uz=0, ur free
  for (let i = 0; i < nrn; i++) {
    const nd = i; // j=0
    constrained.set(nd*2+1, 0.0); // uz only
  }
  // Top face: sliding — uz=delta, ur free
  for (let i = 0; i < nrn; i++) {
    const nd = NZ*nrn + i; // j=NZ
    constrained.set(nd*2+1, delta); // uz only
  }
  // Axis nodes (i=0): ur=0 (symmetry)
  for (let j = 0; j < NZ+1; j++) {
    const nd = j*nrn + 0;
    constrained.set(nd*2, 0.0);
  }

  // Build reduced system via penalty method for simplicity
  const penalty = 1e14 * E;
  const f = new Float64Array(ndof);

  for (const [dof, val] of constrained) {
    Kg[dof*ndof + dof] += penalty;
    f[dof] += penalty * val;
  }

  // Solve with Cholesky (dense, simple — mesh is small)
  u = choleskysolve(Kg, f, ndof);
  computeStress();
}

function solve() {
  if (nonlinear) {
    solveNonlinear();
    computeStressNH();
  } else {
    solveLinear();
  }
}

// ── Dense Cholesky solve ──────────────────────────────────────────────────────
function choleskysolve(A, b, n) {
  // In-place Cholesky decomposition (A is modified to L)
  const L = new Float64Array(n*n);
  for (let i=0;i<n;i++){
    for (let j=0;j<=i;j++){
      let s = A[i*n+j];
      for (let k=0;k<j;k++) s -= L[i*n+k]*L[j*n+k];
      L[i*n+j] = (i===j) ? Math.sqrt(Math.max(s,1e-30)) : s/L[j*n+j];
    }
  }
  // Forward substitution Ly = b
  const y = new Float64Array(n);
  for (let i=0;i<n;i++){
    let s=b[i];
    for (let k=0;k<i;k++) s-=L[i*n+k]*y[k];
    y[i]=s/L[i*n+i];
  }
  // Back substitution L^T x = y
  const x = new Float64Array(n);
  for (let i=n-1;i>=0;i--){
    let s=y[i];
    for (let k=i+1;k<n;k++) s-=L[k*n+i]*x[k];
    x[i]=s/L[i*n+i];
  }
  return x;
}

// ── Neo-Hookean internal force vector ─────────────────────────────────────────
// Nearly-incompressible Neo-Hookean: Ψ = G/2*(Ī₁-3) + κ/2*(J-1)²
// Returns internal force vector fint of length ndof
function assembleFint(uCurr) {
  const ndof = nNodes * 2;
  const fint = new Float64Array(ndof);
  const kappa = E / (3*(1-2*nu));
  const G     = E / (2*(1+nu));
  const nrn = NR + 1;

  for (let e = 0; e < nElems; e++) {
    const n0=elems[e*4],n1=elems[e*4+1],n2=elems[e*4+2],n3=elems[e*4+3];
    const econn=[n0,n1,n2,n3];
    const rc=[nodes[n0*2],nodes[n1*2],nodes[n2*2],nodes[n3*2]];
    const zc=[nodes[n0*2+1],nodes[n1*2+1],nodes[n2*2+1],nodes[n3*2+1]];
    // element displacement
    const ue=[uCurr[n0*2],uCurr[n0*2+1],uCurr[n1*2],uCurr[n1*2+1],
              uCurr[n2*2],uCurr[n2*2+1],uCurr[n3*2],uCurr[n3*2+1]];

    // 2×2 Gauss
    for (let gi=0;gi<2;gi++) for (let gj=0;gj<2;gj++) {
      const xi=GP[gi], eta=GP[gj];
      const N=shape(xi,eta), dN=shapeDeriv(xi,eta);
      // Jacobian of reference map
      let J00=0,J01=0,J10=0,J11=0;
      for(let k=0;k<4;k++){J00+=dN[0][k]*rc[k];J01+=dN[0][k]*zc[k];J10+=dN[1][k]*rc[k];J11+=dN[1][k]*zc[k];}
      const detJ0 = J00*J11-J01*J10;
      const inv00=J11/detJ0, inv01=-J01/detJ0, inv10=-J10/detJ0, inv11=J00/detJ0;
      const dNr=new Float64Array(4),dNz=new Float64Array(4);
      for(let k=0;k<4;k++){dNr[k]=inv00*dN[0][k]+inv10*dN[1][k];dNz[k]=inv01*dN[0][k]+inv11*dN[1][k];}
      let rg=0; for(let k=0;k<4;k++) rg+=N[k]*rc[k];
      rg = Math.max(rg, 1e-12);

      // Deformation gradient F (axisymmetric: 4 independent components)
      // F = I + grad(u);  components: rr, rz, zr, zz, θθ
      let urr=0,urz=0,uzr=0,uzz=0,ur_r=0;
      for(let k=0;k<4;k++){
        urr += dNr[k]*ue[2*k];
        urz += dNz[k]*ue[2*k];
        uzr += dNr[k]*ue[2*k+1];
        uzz += dNz[k]*ue[2*k+1];
        // Near axis use ∂u_r/∂r (L'Hôpital), away from axis use u_r/r
        ur_r += (rg < 1e-6 ? dNr[k] : N[k]/rg) * ue[2*k];
      }
      const Frr=1+urr, Frz=urz, Fzr=uzr, Fzz=1+uzz, Ftt=1+ur_r;
      // J = det F (axisymmetric: θθ block decouples)
      const J = (Frr*Fzz - Frz*Fzr)*Ftt;
      if (J < 1e-10) continue; // degenerate

      // Left Cauchy-Green b (r-z 2×2 part + θθ)
      const brr = Frr*Frr + Frz*Frz;
      const bzz = Fzr*Fzr + Fzz*Fzz;
      const brz = Frr*Fzr + Frz*Fzz;
      const btt = Ftt*Ftt;

      // Isochoric part: J^(-2/3)*b
      const Jm23 = Math.pow(J, -2/3);
      const I1bar = Jm23*(brr+bzz+btt);
      const bbar_rr = Jm23*brr, bbar_zz = Jm23*bzz, bbar_tt = Jm23*btt, bbar_rz = Jm23*brz;

      // Kirchhoff stress τ = G*(bbar - I1bar/3*I) + κ*J*(J-1)*I
      const p_vol = kappa*J*(J-1);
      const I3 = I1bar/3;
      const tau_rr = G*(bbar_rr - I3) + p_vol;
      const tau_zz = G*(bbar_zz - I3) + p_vol;
      const tau_tt = G*(bbar_tt - I3) + p_vol;
      const tau_rz = G*bbar_rz;

      // First Piola-Kirchhoff: P = τ · F^{-T}  (τ = J·σ = Kirchhoff stress)
      // F^{-T} for axisymmetric 2×2 block: [Fzz,-Fzr; -Frz,Frr]/det(F_rz_block)
      // det of r-z block: d2 = Frr*Fzz - Frz*Fzr
      const d2 = Frr*Fzz - Frz*Fzr; // J = d2 * Ftt
      // F^{-T} 2×2: row r: [Fzz/d2, -Fzr/d2], row z: [-Frz/d2, Frr/d2]
      // P_ij = sum_k tau_ik * (F^{-T})_kj  (Voigt: rr,zz,rz,zr)
      // For axisymmetric we need P_rR, P_rZ, P_zR, P_zZ, P_tT
      // τ is symmetric: τ_rr, τ_zz, τ_rz=τ_zr, τ_tt
      const iFrr= Fzz/d2, iFrz=-Fzr/d2, iFzr=-Frz/d2, iFzz= Frr/d2, iFtt=1/Ftt;
      // P_rR = τ_rr*iFrr + τ_rz*iFzr
      const P_rR = tau_rr*iFrr + tau_rz*iFzr;
      // P_rZ = τ_rr*iFrz + τ_rz*iFzz
      const P_rZ = tau_rr*iFrz + tau_rz*iFzz;
      // P_zR = τ_rz*iFrr + τ_zz*iFzr
      const P_zR = tau_rz*iFrr + tau_zz*iFzr;
      // P_zZ = τ_rz*iFrz + τ_zz*iFzz
      const P_zZ = tau_rz*iFrz + tau_zz*iFzz;
      // P_tT = τ_tt * iFtt
      const P_tT = tau_tt * iFtt;

      // Virtual work in reference config: fint += ∫ P : δF dV_ref
      // dV_ref = 2π R detJ0 w  (R = reference radius = rg)
      // δF_rR = dNr_a, δF_rZ = dNz_a (for δu_r dof)
      // δF_zR = dNr_a, δF_zZ = dNz_a (for δu_z dof)
      // δF_tT = N_a/R               (for δu_r dof, circumferential)
      const wt = GW[gi]*GW[gj]*detJ0*2*Math.PI*rg;
      for(let k=0;k<4;k++){
        const NR_k = (rg < 1e-6) ? dNr[k] : N[k]/rg; // N/R with L'Hopital at axis
        const fi_r = (P_rR*dNr[k] + P_rZ*dNz[k] + P_tT*NR_k)*wt;
        const fi_z = (P_zR*dNr[k] + P_zZ*dNz[k])*wt;
        fint[econn[k]*2]   += fi_r;
        fint[econn[k]*2+1] += fi_z;
      }
    }
  }
  return fint;
}

// ── Newton-Raphson solve (Neo-Hookean) ────────────────────────────────────────
function solveNonlinear() {
  const ndof = nNodes * 2;
  const nrn = NR + 1;
  const penalty = 1e10 * E;
  const MAX_ITER = 20;
  const TOL = 1e-8;

  // Build constrained DOF map (same as linear)
  const constrained = new Map();
  for (let i=0;i<nrn;i++) constrained.set(i*2+1, 0.0);                  // bottom uz=0
  for (let i=0;i<nrn;i++) constrained.set((NZ*nrn+i)*2+1, delta);        // top uz=delta
  for (let j=0;j<=NZ;j++)  constrained.set(j*nrn*2, 0.0);                // axis ur=0

  // External force vector (only penalty contributions for prescribed DOFs)
  const fext = new Float64Array(ndof);
  for (const [dof, val] of constrained) fext[dof] = penalty * val;

  // Warm-start from linear solution
  solveLinear();
  let uCurr = new Float64Array(u);
  // Enforce BCs cleanly
  for (const [dof, val] of constrained) uCurr[dof] = val;

  for (let iter=0; iter<MAX_ITER; iter++) {
    // Internal force
    const fint = assembleFint(uCurr);
    // Add penalty for BCs to fint: penalty*(u_dof - val)
    for (const [dof, val] of constrained) fint[dof] += penalty*(uCurr[dof]-val);

    // Residual R = fint - fext
    const R = new Float64Array(ndof);
    for (let i=0;i<ndof;i++) R[i] = fint[i] - fext[i];

    // Check convergence: ||R|| / (||fext|| + 1)
    let normR=0, normF=0;
    for (let i=0;i<ndof;i++){normR+=R[i]*R[i];normF+=fext[i]*fext[i];}
    normR=Math.sqrt(normR); normF=Math.sqrt(normF);
    if (normR / (normF + 1e-30*E) < TOL) break;

    // Numerical tangent KT via column perturbation
    const eps_fd = 1e-7;
    const KT = new Float64Array(ndof*ndof);
    for (let j=0;j<ndof;j++){
      const uPlus = new Float64Array(uCurr);
      uPlus[j] += eps_fd;
      const fp = assembleFint(uPlus);
      for (const [dof, val] of constrained) fp[dof] += penalty*(uPlus[dof]-val);
      for (let i=0;i<ndof;i++) KT[i*ndof+j] = (fp[i]-fint[i])/eps_fd;
    }

    // Solve KT * du = -R
    const du = choleskysolve(KT, R.map(v=>-v), ndof);
    for (let i=0;i<ndof;i++) uCurr[i] += du[i];
  }

  u = uCurr;
}

// ── Compute Cauchy stress at nodes for Neo-Hookean ────────────────────────────
function computeStressNH() {
  const kappa = E / (3*(1-2*nu));
  const G     = E / (2*(1+nu));
  const nrn = NR + 1;
  stressNodes = new Float64Array(nNodes * 5);
  const count = new Float64Array(nNodes);

  for (let e=0; e<nElems; e++){
    const n0=elems[e*4],n1=elems[e*4+1],n2=elems[e*4+2],n3=elems[e*4+3];
    const econn=[n0,n1,n2,n3];
    const rc=[nodes[n0*2],nodes[n1*2],nodes[n2*2],nodes[n3*2]];
    const zc=[nodes[n0*2+1],nodes[n1*2+1],nodes[n2*2+1],nodes[n3*2+1]];
    const ue=[u[n0*2],u[n0*2+1],u[n1*2],u[n1*2+1],u[n2*2],u[n2*2+1],u[n3*2],u[n3*2+1]];

    const xi=0,eta=0;
    const N=shape(xi,eta),dN=shapeDeriv(xi,eta);
    let J00=0,J01=0,J10=0,J11=0;
    for(let k=0;k<4;k++){J00+=dN[0][k]*rc[k];J01+=dN[0][k]*zc[k];J10+=dN[1][k]*rc[k];J11+=dN[1][k]*zc[k];}
    const detJ0=J00*J11-J01*J10;
    const inv00=J11/detJ0,inv01=-J01/detJ0,inv10=-J10/detJ0,inv11=J00/detJ0;
    const dNr=new Float64Array(4),dNz=new Float64Array(4);
    for(let k=0;k<4;k++){dNr[k]=inv00*dN[0][k]+inv10*dN[1][k];dNz[k]=inv01*dN[0][k]+inv11*dN[1][k];}
    let rg=0; for(let k=0;k<4;k++) rg+=N[k]*rc[k];
    rg=Math.max(rg,1e-12);

    let urr=0,urz=0,uzr=0,uzz=0,ur_r=0;
    for(let k=0;k<4;k++){
      urr+=dNr[k]*ue[2*k]; urz+=dNz[k]*ue[2*k];
      uzr+=dNr[k]*ue[2*k+1]; uzz+=dNz[k]*ue[2*k+1];
      ur_r+=(rg < 1e-6 ? dNr[k] : N[k]/rg)*ue[2*k];
    }
    const Frr=1+urr,Frz=urz,Fzr=uzr,Fzz=1+uzz,Ftt=1+ur_r;
    const J=(Frr*Fzz-Frz*Fzr)*Ftt;
    if(J<1e-10){count[n0]++;count[n1]++;count[n2]++;count[n3]++;continue;}

    const brr=Frr*Frr+Frz*Frz, bzz=Fzr*Fzr+Fzz*Fzz, brz=Frr*Fzr+Frz*Fzz, btt=Ftt*Ftt;
    const Jm23=Math.pow(J,-2/3);
    const I1bar=Jm23*(brr+bzz+btt);
    const bbar_rr=Jm23*brr,bbar_zz=Jm23*bzz,bbar_tt=Jm23*btt,bbar_rz=Jm23*brz;
    const p_vol=kappa*J*(J-1), I3=I1bar/3, Jfac=1/J;

    const sr =(G*(bbar_rr-I3)+p_vol)*Jfac;
    const sz =(G*(bbar_zz-I3)+p_vol)*Jfac;
    const st =(G*(bbar_tt-I3)+p_vol)*Jfac;
    const trz=(G*bbar_rz)*Jfac;
    const vm =Math.sqrt(0.5*((sr-sz)**2+(sz-st)**2+(st-sr)**2+6*trz**2));

    for(let k=0;k<4;k++){
      const nd=econn[k];
      stressNodes[nd*5+0]+=sr; stressNodes[nd*5+1]+=sz;
      stressNodes[nd*5+2]+=st; stressNodes[nd*5+3]+=trz;
      stressNodes[nd*5+4]+=vm; count[nd]++;
    }
  }
  for(let nd=0;nd<nNodes;nd++) if(count[nd]>0) for(let s=0;s<5;s++) stressNodes[nd*5+s]/=count[nd];
}

// ── Compute stress at nodes (average from surrounding elements) ───────────────
function computeStress() {
  // stress at element centre, then average to nodes
  const nrn = NR + 1;
  stressNodes = new Float64Array(nNodes * 5); // sr,sz,st,trz,vm
  const count  = new Float64Array(nNodes);

  // Axisymmetric D
  const c = E / ((1+nu)*(1-2*nu));
  const D = [
    c*(1-nu), c*nu,    c*nu,    0,
    c*nu,     c*(1-nu),c*nu,    0,
    c*nu,     c*nu,    c*(1-nu),0,
    0,        0,       0,       c*(1-2*nu)/2,
  ];

  for (let e=0; e<nElems; e++){
    const n0=elems[e*4],n1=elems[e*4+1],n2=elems[e*4+2],n3=elems[e*4+3];
    const econn=[n0,n1,n2,n3];
    const rc=[nodes[n0*2],nodes[n1*2],nodes[n2*2],nodes[n3*2]];
    const zc=[nodes[n0*2+1],nodes[n1*2+1],nodes[n2*2+1],nodes[n3*2+1]];
    const ue=[u[n0*2],u[n0*2+1],u[n1*2],u[n1*2+1],u[n2*2],u[n2*2+1],u[n3*2],u[n3*2+1]];

    // Evaluate at element centroid (xi=eta=0)
    const xi=0, eta=0;
    const N=shape(xi,eta);
    const dN=shapeDeriv(xi,eta);
    let J00=0,J01=0,J10=0,J11=0;
    for(let k=0;k<4;k++){J00+=dN[0][k]*rc[k];J01+=dN[0][k]*zc[k];J10+=dN[1][k]*rc[k];J11+=dN[1][k]*zc[k];}
    const detJ=J00*J11-J01*J10;
    const invJ00=J11/detJ,invJ01=-J01/detJ,invJ10=-J10/detJ,invJ11=J00/detJ;
    const dNr=new Float64Array(4),dNz=new Float64Array(4);
    for(let k=0;k<4;k++){dNr[k]=invJ00*dN[0][k]+invJ10*dN[1][k];dNz[k]=invJ01*dN[0][k]+invJ11*dN[1][k];}
    let rg=0; for(let k=0;k<4;k++) rg+=N[k]*rc[k];

    // strain
    let er=0,ez=0,et=0,grz=0;
    for(let k=0;k<4;k++){
      er  += dNr[k]*ue[2*k];
      ez  += dNz[k]*ue[2*k+1];
      et  += N[k]*ue[2*k]/(rg+1e-14);
      grz += dNz[k]*ue[2*k] + dNr[k]*ue[2*k+1];
    }

    // stress = D * strain
    const sr  = D[0]*er + D[1]*ez + D[2]*et;
    const sz  = D[4]*er + D[5]*ez + D[6]*et;
    const st  = D[8]*er + D[9]*ez + D[10]*et;
    const trz = D[15]*grz;

    // von Mises for axisymmetric
    const vm = Math.sqrt(0.5*((sr-sz)**2+(sz-st)**2+(st-sr)**2+6*trz**2));

    // scatter to nodes
    for(let k=0;k<4;k++){
      const nd=econn[k];
      stressNodes[nd*5+0]+=sr; stressNodes[nd*5+1]+=sz;
      stressNodes[nd*5+2]+=st; stressNodes[nd*5+3]+=trz;
      stressNodes[nd*5+4]+=vm;
      count[nd]++;
    }
  }
  for(let nd=0;nd<nNodes;nd++){
    if(count[nd]>0){
      for(let s=0;s<5;s++) stressNodes[nd*5+s]/=count[nd];
    }
  }
}

// ── Build Three.js revolution surface ────────────────────────────────────────
// Revolve the deformed r-z cross-section around Z axis
const N_PHI = 32; // azimuthal subdivisions

function buildRevolutionMesh() {
  const nrn = NR + 1, nzn = NZ + 1;
  const HALF = N_PHI / 2; // opaque half: k=0..HALF-1, ghost half: k=HALF..N_PHI-1

  // Deformed node positions
  const rDef = new Float64Array(nNodes);
  const zDef = new Float64Array(nNodes);
  for (let nd=0; nd<nNodes; nd++){
    rDef[nd] = nodes[nd*2]   + u[nd*2];
    zDef[nd] = nodes[nd*2+1] + u[nd*2+1];
  }

  // Deformed extents — used to centre the mesh at y=0
  const zBot = zDef[0];
  const zTop = zDef[NZ * (NR+1)];
  zCenter = (zBot + zTop) * 0.5;
  H_def   = zTop - zBot;

  // Stress values
  const sVal = new Float64Array(nNodes);
  let sMin = Infinity, sMax = -Infinity;
  for (let nd=0; nd<nNodes; nd++){
    sVal[nd] = stressNodes[nd*5 + stressComponent];
    if (sVal[nd] < sMin) sMin = sVal[nd];
    if (sVal[nd] > sMax) sMax = sVal[nd];
  }
  const sRange = sMax - sMin + 1e-14;

  // Helper: build position/color arrays for a set of (nd, phi) pairs
  // Returns { positions, colors, indices } for a revolution surface
  // kStart..kEnd (exclusive) selects which phi slices to include
  function makeRevSurface(kStart, kEnd) {
    const nk = kEnd - kStart;
    const nVerts = nNodes * nk;
    const pos = new Float32Array(nVerts * 3);
    const col = new Float32Array(nVerts * 3);

    // vtx(nd, k_local) = nd*nk + k_local
    for (let j=0; j<nzn; j++){
      for (let i=0; i<nrn; i++){
        const nd = j*nrn + i;
        const r = rDef[nd], z = zDef[nd] - zCenter;
        const t = (sVal[nd] - sMin) / sRange;
        const c = stressToColor(t);
        for (let kl=0; kl<nk; kl++){
          const k = kStart + kl;
          const phi = 2*Math.PI*k/N_PHI;
          const vi = nd*nk + kl;
          pos[vi*3]   = r*Math.cos(phi);
          pos[vi*3+1] = z;
          pos[vi*3+2] = r*Math.sin(phi);
          col[vi*3]   = c.r; col[vi*3+1] = c.g; col[vi*3+2] = c.b;
        }
      }
    }

    const idx = [];
    const vt = (nd, kl) => nd*nk + kl;

    // Outer curved surface (i = NR only)
    for (let j=0; j<NZ; j++){
      const nb = j*nrn + NR, nt = (j+1)*nrn + NR;
      for (let kl=0; kl<nk-1; kl++){
        const b=vt(nb,kl), b1=vt(nb,kl+1), t=vt(nt,kl), t1=vt(nt,kl+1);
        idx.push(b,t,t1, b,t1,b1);
      }
    }

    // Top cap (j=NZ)
    for (let i=0; i<NR; i++){
      const nd0=NZ*nrn+i, nd1=NZ*nrn+i+1;
      for (let kl=0; kl<nk-1; kl++){
        const v0=vt(nd0,kl), v1=vt(nd1,kl), v0b=vt(nd0,kl+1), v1b=vt(nd1,kl+1);
        idx.push(v0,v1,v1b, v0,v1b,v0b);
      }
    }

    // Bottom cap (j=0)
    for (let i=0; i<NR; i++){
      const nd0=i, nd1=i+1;
      for (let kl=0; kl<nk-1; kl++){
        const v0=vt(nd0,kl), v1=vt(nd1,kl), v0b=vt(nd0,kl+1), v1b=vt(nd1,kl+1);
        idx.push(v0,v0b,v1b, v0,v1b,v1);
      }
    }

    return { pos, col, idx };
  }

  // ── Cross-section faces — two flat half-discs closing the cut at phi=0 and phi=π ──
  // We emit two copies of the r-z mesh: one at x=+r (phi=0) and one at x=-r (phi=π)
  function makeCrossSection() {
    const nVerts = nNodes * 2; // block 0: phi=0 (x=+r), block 1: phi=π (x=-r)
    const pos = new Float32Array(nVerts * 3);
    const col = new Float32Array(nVerts * 3);

    for (let nd=0; nd<nNodes; nd++){
      const r = rDef[nd], z = zDef[nd] - zCenter;
      const t = (sVal[nd] - sMin) / sRange;
      const c = stressToColor(t);
      // phi=0 block
      pos[nd*3]   =  r; pos[nd*3+1] = z; pos[nd*3+2] = 0;
      col[nd*3]   = c.r; col[nd*3+1] = c.g; col[nd*3+2] = c.b;
      // phi=π block
      const nd2 = nNodes + nd;
      pos[nd2*3]   = -r; pos[nd2*3+1] = z; pos[nd2*3+2] = 0;
      col[nd2*3]   = c.r; col[nd2*3+1] = c.g; col[nd2*3+2] = c.b;
    }

    const idx = [];
    for (let b=0; b<2; b++){
      const off = b * nNodes;
      for (let j=0; j<NZ; j++){
        for (let i=0; i<NR; i++){
          const bl=off+j*nrn+i, br=off+j*nrn+i+1;
          const tl=off+(j+1)*nrn+i, tr=off+(j+1)*nrn+i+1;
          idx.push(bl,br,tr, bl,tr,tl);
        }
      }
    }

    return { pos, col, idx };
  }

  function makeGeo({ pos, col, idx }) {
    const geo = new THREE.BufferGeometry();
    geo.setAttribute('position', new THREE.Float32BufferAttribute(pos, 3));
    geo.setAttribute('color',    new THREE.Float32BufferAttribute(col, 3));
    geo.setIndex(idx);
    geo.computeVertexNormals();
    return geo;
  }

  // ── Opaque half (k=0..HALF-1) ─────────────────────────────────────────────
  const opaqueMat = new THREE.MeshPhongMaterial({
    vertexColors: true,
    side: THREE.DoubleSide,
    shininess: 60,
  });
  const opaqueMesh = new THREE.Mesh(makeGeo(makeRevSurface(HALF, N_PHI+1)), opaqueMat);

  // ── Cross-section face ────────────────────────────────────────────────────
  const crossMat = new THREE.MeshPhongMaterial({
    vertexColors: true,
    side: THREE.DoubleSide,
    shininess: 30,
  });
  const crossMesh = new THREE.Mesh(makeGeo(makeCrossSection()), crossMat);

  // ── Ghost half (k=HALF..N_PHI-1) ─────────────────────────────────────────
  const [pgr, pgg, pgb] = _rgb('--text-gray');
  const ghostMat = new THREE.MeshPhongMaterial({
    color: new THREE.Color(pgr/255, pgg/255, pgb/255),
    transparent: true,
    opacity: 0.18,
    side: THREE.DoubleSide,
    depthWrite: false,
    shininess: 40,
  });
  const ghostMesh = new THREE.Mesh(makeGeo(makeRevSurface(0, HALF+1)), ghostMat);

  // ── Cross-section wireframe ───────────────────────────────────────────────
  // Draw element edges on both cut faces (phi=0 and phi=π)
  function makeCrossWire() {
    const lines = [];
    const nrn = NR + 1;
    function pushEdge(nd0, nd1, xSign) {
      const r0 = rDef[nd0], z0 = zDef[nd0] - zCenter;
      const r1 = rDef[nd1], z1 = zDef[nd1] - zCenter;
      lines.push(xSign*r0, z0, 0,  xSign*r1, z1, 0);
    }
    for (let s = -1; s <= 1; s += 2) { // s=-1 → phi=π (x=-r), s=+1 → phi=0 (x=+r)
      for (let j=0; j<=NZ; j++){
        for (let i=0; i<NR; i++){
          pushEdge(j*nrn+i, j*nrn+i+1, s); // horizontal edges
        }
      }
      for (let i=0; i<=NR; i++){
        for (let j=0; j<NZ; j++){
          pushEdge(j*nrn+i, (j+1)*nrn+i, s); // vertical edges
        }
      }
    }
    const geo = new THREE.BufferGeometry();
    geo.setAttribute('position', new THREE.Float32BufferAttribute(new Float32Array(lines), 3));
    return geo;
  }

  const wireMat = new THREE.LineBasicMaterial({
    color: new THREE.Color(_c('--bg-void')),
    transparent: true,
    opacity: 0.5,
  });
  const wireMesh = new THREE.LineSegments(makeCrossWire(), wireMat);

  const group = new THREE.Group();
  group.add(opaqueMesh, crossMesh, wireMesh, ghostMesh);
  return group;
}

// ── Plate meshes ──────────────────────────────────────────────────────────────
// PLATE_R is computed dynamically in buildPlates based on incompressible radius
const PLATE_T       = H0 * 0.04;   // plate thickness
let plateMeshes = [];

function buildPlates() {
  // Remove old plates
  plateMeshes.forEach(m => { scene.remove(m); m.geometry.dispose(); m.material.dispose(); });
  plateMeshes = [];

  const [pgr, pgg, pgb] = _rgb('--text-gray');
  const plateCol = new THREE.Color(pgr/255, pgg/255, pgb/255);
  const plateMat = new THREE.MeshPhongMaterial({
    color: plateCol,
    emissive: new THREE.Color(pgr/255*0.15, pgg/255*0.15, pgb/255*0.15),
    shininess: 120,
    side: THREE.DoubleSide,
  });

  // Place plates symmetrically around the centred deformed cylinder
  const eps = PLATE_T * 0.5;
  const botZ = -H_def/2 - PLATE_T/2 - eps;
  const topZ =  H_def/2 + PLATE_T/2 + eps;

  // Plate radius = incompressible radius at max compression (fixed)
  const PLATE_R = R0 * Math.sqrt(H0 / (H0 - deltaMax));
  const geo = new THREE.CylinderGeometry(PLATE_R, PLATE_R, PLATE_T, 48, 1);

  const botMesh = new THREE.Mesh(geo, plateMat);
  botMesh.position.set(0, botZ, 0);
  scene.add(botMesh); plateMeshes.push(botMesh);

  const topMesh = new THREE.Mesh(geo, plateMat.clone());
  topMesh.position.set(0, topZ, 0);
  scene.add(topMesh); plateMeshes.push(topMesh);
}

// ── Three.js init ─────────────────────────────────────────────────────────────
function initThree() {
  simCanvas = document.getElementById('fem-canvas');
  renderer = new THREE.WebGLRenderer({ canvas: simCanvas, antialias: true });
  renderer.setPixelRatio(window.devicePixelRatio);
  renderer.setClearColor(new THREE.Color(_c('--bg-void')), 1);

  camera = new THREE.PerspectiveCamera(45, 1, 0.01, 100);
  updateCamera();
  setupOrbit();

  scene = new THREE.Scene();
  scene.add(new THREE.AmbientLight(0xffffff, 0.5));
  const dl = new THREE.DirectionalLight(0xffffff, 0.8);
  dl.position.set(4, 6, 5); scene.add(dl);
  const dl2 = new THREE.DirectionalLight(0xffffff, 0.3);
  dl2.position.set(-3, -2, -4); scene.add(dl2);

  resizeRenderer();
}

function disposeGroup(g) {
  g.children.forEach(m => { m.geometry.dispose(); m.material.dispose(); });
  g.clear();
}

function updateScene() {
  if (cylMesh) { scene.remove(cylMesh); disposeGroup(cylMesh); }
  cylMesh = buildRevolutionMesh();
  scene.add(cylMesh);
  buildPlates();
  // Zoom camera to track deformed height (H_def set by buildRevolutionMesh)
  orbit.radius = BASE_RADIUS * Math.max(1, H_def / H0);
  updateCamera();
}

function updateCamera() {
  const { theta, phi, radius } = orbit;
  camera.position.set(
    radius * Math.sin(phi) * Math.cos(theta),
    radius * Math.cos(phi),
    radius * Math.sin(phi) * Math.sin(theta)
  );
  camera.lookAt(0, 0, 0);
  const el = document.getElementById('fem-cam-info');
  if (el) el.textContent = 'θ=' + theta.toFixed(2) + ' φ=' + phi.toFixed(2);
}

function resizeRenderer() {
  const w = simCanvas.clientWidth, h = simCanvas.clientHeight;
  renderer.setSize(w, h, false);
  camera.aspect = w / h;
  camera.updateProjectionMatrix();
}

function setupOrbit() {
  simCanvas.addEventListener('pointerdown', e => {
    orbit.dragging=true; orbit.lastX=e.clientX; orbit.lastY=e.clientY;
    simCanvas.setPointerCapture(e.pointerId);
  });
  simCanvas.addEventListener('pointermove', e => {
    if (!orbit.dragging) return;
    orbit.theta -= (e.clientX-orbit.lastX)*0.008;
    orbit.phi = Math.max(0.05, Math.min(Math.PI-0.05, orbit.phi+(e.clientY-orbit.lastY)*0.008));
    orbit.lastX=e.clientX; orbit.lastY=e.clientY; updateCamera();
  });
  simCanvas.addEventListener('pointerup', () => { orbit.dragging=false; });
  simCanvas.addEventListener('wheel', e => {
    e.preventDefault();
    orbit.radius = Math.max(1, Math.min(20, orbit.radius + e.deltaY*0.01));
    updateCamera();
  }, { passive: false });
}

// ── Animation loop ────────────────────────────────────────────────────────────
function loop() {
  if (running) {
    if (needsSolve) {
      solve();
      // Guard: don't update scene if solver produced NaN
      if (u && !isNaN(u[0])) {
        updateScene();
      }
      needsSolve = false;
    }
    renderer.render(scene, camera);
  }
  frameId = requestAnimationFrame(loop);
}

// ── Shell wiring ──────────────────────────────────────────────────────────────
const shell = new AppletShell({
  id:    'fem',
  title: 'Axisymmetric FEM &mdash; Elastic Cylinder',
  gap:   0,

  ctrlHTML: `
    <div class="applet-shell-ctrl-section">
      <div class="applet-shell-ctrl-title">Actions</div>
      <div class="applet-shell-btn-row">
        <button class="applet-shell-btn" onclick="femReset()">Reset</button>
        <button class="applet-shell-btn" id="fem-pause-btn" onclick="femTogglePause()">Pause</button>
        <!-- <button class="applet-shell-btn" id="fem-nl-btn" onclick="femToggleNonlinear()">Linear</button> -->
      </div>
    </div>
    <div class="applet-shell-ctrl-section">
      <div class="applet-shell-ctrl-title">Loading</div>
      <div class="applet-shell-slider-row">
        <span class="applet-shell-side">Compression</span>
        <div class="applet-shell-slider-wrap">
          <input type="range" id="fem-delta" min="-1" max="1" step="0.01" value="0">
          <div class="applet-shell-tick" style="left:50%"></div>
        </div>
        <span class="applet-shell-side">Tension</span>
      </div>
    </div>
    <div class="applet-shell-ctrl-section">
      <div class="applet-shell-ctrl-title">Young's Modulus</div>
      <div class="applet-shell-slider-row">
        <span class="applet-shell-side">Soft</span>
        <input type="range" id="fem-E" min="0.1" max="3" step="0.01" value="1">
        <span class="applet-shell-side">Stiff</span>
      </div>
    </div>
    <div class="applet-shell-ctrl-section">
      <div class="applet-shell-ctrl-title">Poisson Ratio</div>
      <div class="applet-shell-slider-row">
        <span class="applet-shell-side">Auxetic</span>
        <div class="applet-shell-slider-wrap">
          <input type="range" id="fem-nu" min="-0.25" max="0.499" step="0.001" value="0.3">
          <div class="applet-shell-tick" style="left:33.38%"></div>
        </div>
        <span class="applet-shell-side">Incompressible</span>
      </div>
    </div>
    <div class="applet-shell-ctrl-section">
      <div class="applet-shell-ctrl-title">Stress Component</div>
      <div class="applet-shell-btn-row" id="fem-stress-btns"></div>
    </div>
    <div class="applet-shell-ctrl-section">
      <div class="applet-shell-ctrl-title">Camera</div>
      <div class="applet-shell-val" id="fem-cam-info" style="width:100%;text-align:left;color:var(--text-dim);">θ=0.00 φ=0.00</div>
    </div>
  `,

  onOpen: function ({ canvas: c, S }) {
    PALETTE = buildPalette();
    buildMesh();
    delta = 0; E = 1; nu = 0.3; nonlinear = false;
    u = null;
    needsSolve = true;

    // stress buttons
    const labels = ['von Mises','σ_r','σ_z','σ_θ','τ_rz'];
    const row = document.getElementById('fem-stress-btns');
    row.innerHTML = '';
    labels.forEach((lb, i) => {
      const btn = document.createElement('button');
      btn.className = 'applet-shell-btn' + (i===0 ? ' active' : '');
      btn.textContent = lb;
      btn.addEventListener('click', () => {
        stressComponent = i;
        document.querySelectorAll('#fem-stress-btns .applet-shell-btn').forEach((b,j) => b.classList.toggle('active', j===i));
        needsSolve = true;
      });
      row.appendChild(btn);
    });

    const pb = document.getElementById('fem-pause-btn');
    if (pb) { pb.textContent = 'Pause'; pb.classList.remove('active'); }
    const nlb = document.getElementById('fem-nl-btn');
    if (nlb) { nlb.textContent = 'Linear'; nlb.classList.remove('active'); }

    function startThree() {
      setTimeout(() => {
        initThree();
        solve();
        updateScene();
        needsSolve = false;
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
  },

  onResize: function () {
    if (renderer) resizeRenderer();
  },
});

window.femOpen  = () => shell.open();
window.femClose = () => shell.close();

window.femReset = function () {
  delta = 0; E = 1; nu = 0.3;
  document.getElementById('fem-delta').value = 0;
  document.getElementById('fem-E').value = 1;
  document.getElementById('fem-nu').value = 0.3;
  u = null;
  needsSolve = true;
};

window.femToggleNonlinear = function () {
  nonlinear = !nonlinear;
  const btn = document.getElementById('fem-nl-btn');
  if (btn) {
    btn.textContent = nonlinear ? 'Neo-Hookean' : 'Linear';
    btn.classList.toggle('active', nonlinear);
  }
  u = null; // reset displacement
  needsSolve = true;
};

window.femTogglePause = function () {
  running = !running;
  const pb = document.getElementById('fem-pause-btn');
  if (pb) { pb.textContent = running ? 'Pause' : 'Resume'; pb.classList.toggle('active', !running); }
};

document.getElementById('fem-delta').addEventListener('input', function () {
  delta = parseFloat(this.value) * deltaMax;
  needsSolve = true;
});
document.getElementById('fem-E').addEventListener('input', function () {
  E = parseFloat(this.value);
  needsSolve = true;
});
document.getElementById('fem-nu').addEventListener('input', function () {
  nu = parseFloat(this.value);
  needsSolve = true;
});

})();
