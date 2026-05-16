(function () {
  'use strict';

  const _cs  = getComputedStyle(document.documentElement);
  const _c   = n => _cs.getPropertyValue(n).trim();
  const _rgb  = n => { const h = _c(n).replace('#',''); const v = parseInt(h,16); return [(v>>16)&0xFF,(v>>8)&0xFF,v&0xFF]; };
  const _rgba = (n, a) => { const [r,g,b] = _rgb(n); return `rgba(${r},${g},${b},${a})`; };

  /* ── Inject CSS ── */
  (function () {
    if (document.getElementById('pend-styles')) return;
    const s = document.createElement('style');
    s.id = 'pend-styles';
    s.textContent = `
      #pend-ctrl-panel { display: flex; flex-direction: column; overflow: hidden; }
      #pend-scrollable  { flex: 1; overflow-y: auto; min-height: 0; }
    `;
    document.head.appendChild(s);
  })();

  /* ── Canvas & shell state ── */
  let canvas, ctx;
  let S = 500, W = 500, H = 500;
  let running = false, frameId = null;

  /* ── Mode ── */
  let currentMode = 'double';

  const TOTAL_LENGTH_FRAC = 0.72;
  const G = 9.81;

  /* ── Zoom (shared) ── */
  let zoomT = 0.0;
  const ZOOM_MIN = 1.0;
  function totalPx()  { return H * TOTAL_LENGTH_FRAC; }
  function zoomMax()  { return (H * 0.5 * 0.88) / totalPx(); }
  function zoom()     { return ZOOM_MIN + (zoomMax() - ZOOM_MIN) * zoomT; }
  function pivotScreen() {
    return { px: W * 0.5, py: H * (0.12 + 0.38 * zoomT) };
  }

  /* ── Trail (shared) ── */
  const TRAIL_LEN      = 55;
  const TRAIL_INTERVAL = 2;
  let trailBuf  = new Array(TRAIL_LEN).fill(null);  // ring buffer
  let trailHead = 0, trailCount = 0;                // ring buffer pointers
  let trailTick = 0;
  let trailCvs  = null, trailCtx = null;

  /* ════════════════════════════════════════════
     DOUBLE PENDULUM
     ════════════════════════════════════════════ */
  let dp_lenRatio  = 0.5;
  let dp_massRatio = 0.5;
  let dp_damping   = 0.0;

  let dp_th1 = 0.3, dp_th2 = 0.0;
  let dp_w1  = 0,   dp_w2  = 0;

  let dp_grabbing      = false;
  let dp_grabRod       = 0;
  let dp_grabTh2Offset = 0;

  function dp_L1px() { return dp_lenRatio * totalPx(); }
  function dp_L2px() { return (1 - dp_lenRatio) * totalPx(); }
  function dp_l1()   { return dp_lenRatio; }
  function dp_l2()   { return 1 - dp_lenRatio; }
  function dp_m1()   { return dp_massRatio; }
  function dp_m2()   { return 1 - dp_massRatio; }

  function dp_derivs(th1_, th2_, w1_, w2_) {
    const _m1 = dp_m1(), _m2 = dp_m2();
    const _l1 = dp_l1(), _l2 = dp_l2();
    const dth = th1_ - th2_;
    const sd = Math.sin(dth), cd = Math.cos(dth);
    const M11 = (_m1 + _m2) * _l1 * _l1;
    const M12 = _m2 * _l1 * _l2 * cd;
    const M22 = _m2 * _l2 * _l2;
    const f1  = -(_m1 + _m2) * G * _l1 * Math.sin(th1_) - _m2 * _l1 * _l2 * w2_ * w2_ * sd - dp_damping * w1_;
    const f2  = -_m2 * G * _l2 * Math.sin(th2_) + _m2 * _l1 * _l2 * w1_ * w1_ * sd - dp_damping * w2_;
    const det = M11 * M22 - M12 * M12;
    return [w1_, w2_, (M22 * f1 - M12 * f2) / det, (M11 * f2 - M12 * f1) / det];
  }

  function dp_rk4step(dt_) {
    const [k1_w1, k1_w2, k1_a1, k1_a2] = dp_derivs(dp_th1, dp_th2, dp_w1, dp_w2);
    const [k2_w1, k2_w2, k2_a1, k2_a2] = dp_derivs(dp_th1+0.5*dt_*k1_w1, dp_th2+0.5*dt_*k1_w2, dp_w1+0.5*dt_*k1_a1, dp_w2+0.5*dt_*k1_a2);
    const [k3_w1, k3_w2, k3_a1, k3_a2] = dp_derivs(dp_th1+0.5*dt_*k2_w1, dp_th2+0.5*dt_*k2_w2, dp_w1+0.5*dt_*k2_a1, dp_w2+0.5*dt_*k2_a2);
    const [k4_w1, k4_w2, k4_a1, k4_a2] = dp_derivs(dp_th1+dt_*k3_w1, dp_th2+dt_*k3_w2, dp_w1+dt_*k3_a1, dp_w2+dt_*k3_a2);
    dp_th1 = wrapAngle(dp_th1 + (dt_/6)*(k1_w1+2*k2_w1+2*k3_w1+k4_w1));
    dp_th2 = wrapAngle(dp_th2 + (dt_/6)*(k1_w2+2*k2_w2+2*k3_w2+k4_w2));
    dp_w1  += (dt_/6)*(k1_a1+2*k2_a1+2*k3_a1+k4_a1);
    dp_w2  += (dt_/6)*(k1_a2+2*k2_a2+2*k3_a2+k4_a2);
  }

  const DT = 1/60, SUBSTEPS = 8, CH_SUBSTEPS = 4;

  function dp_stepPhysics() {
    const dt = DT / SUBSTEPS;
    for (let i = 0; i < SUBSTEPS; i++) dp_rk4step(dt);
  }

  function dp_joints() {
    const { px, py } = pivotScreen();
    const z = zoom();
    const l1s = dp_L1px() * z, l2s = dp_L2px() * z;
    const x2 = px + l1s * Math.sin(dp_th1), y2 = py + l1s * Math.cos(dp_th1);
    const x3 = x2 + l2s * Math.sin(dp_th2), y3 = y2 + l2s * Math.cos(dp_th2);
    return [{ x: px, y: py }, { x: x2, y: y2 }, { x: x3, y: y3 }];
  }

  function dp_recordTrail() {
    trailTick++;
    if (trailTick % TRAIL_INTERVAL !== 0) return;
    trailBuf[trailHead] = { th: [dp_th1, dp_th2] };
    trailHead = (trailHead + 1) % TRAIL_LEN;
    if (trailCount < TRAIL_LEN) trailCount++;
  }

  function dp_trailJoints(entry, _pivot, _zoom) {
    const { px, py } = _pivot || pivotScreen();
    const z = (_zoom !== undefined) ? _zoom : zoom();
    const l1s = dp_L1px() * z, l2s = dp_L2px() * z;
    const x2 = px + l1s * Math.sin(entry.th[0]), y2 = py + l1s * Math.cos(entry.th[0]);
    const x3 = x2 + l2s * Math.sin(entry.th[1]), y3 = y2 + l2s * Math.cos(entry.th[1]);
    return [{ x: px, y: py }, { x: x2, y: y2 }, { x: x3, y: y3 }];
  }

  function dp_init() {
    dp_th1 = 0.3; dp_th2 = 0.0;
    dp_w1  = 0;   dp_w2  = 0;
    trailHead = 0; trailCount = 0; trailTick = 0;
    dp_grabbing = false;
  }

  function dp_onPointerDown(x, y) {
    const j = dp_joints();
    const d2 = distToSegment(x, y, j[1].x, j[1].y, j[2].x, j[2].y);
    const d1 = distToSegment(x, y, j[0].x, j[0].y, j[1].x, j[1].y);
    if (d2 <= GRAB_THRESH) {
      dp_grabbing = true; dp_grabRod = 2; dp_w1 = 0; dp_w2 = 0; return true;
    } else if (d1 <= GRAB_THRESH) {
      dp_grabbing = true; dp_grabRod = 1; dp_grabTh2Offset = dp_th2 - dp_th1; dp_w1 = 0; dp_w2 = 0; return true;
    }
    return false;
  }

  function dp_onPointerMove(x, y) {
    if (!dp_grabbing) return;
    const { px, py } = pivotScreen();
    const z = zoom();
    if (dp_grabRod === 1) {
      const newTh = Math.atan2(x - px, y - py);
      if (!running) { dp_th1 = newTh; dp_th2 = dp_th1 + dp_grabTh2Offset; }
      else { dp_th1 = newTh; }
      dp_w1 = 0;
    } else {
      const x2 = px + dp_L1px() * z * Math.sin(dp_th1);
      const y2 = py + dp_L1px() * z * Math.cos(dp_th1);
      dp_th2 = Math.atan2(x - x2, y - y2);
      dp_w2  = 0;
    }
  }

  function dp_onPointerUp() { dp_grabbing = false; }

  /* ── Double render ── */
  function dp_render() {
    const j = dp_joints();
    drawRod(ctx, j[0].x, j[0].y, j[1].x, j[1].y, '--teal-light', 5, 18, 1.0);
    drawRod(ctx, j[1].x, j[1].y, j[2].x, j[2].y, '--pink-light',  5, 18, 1.0);
    drawDot(ctx, j[1].x, j[1].y, 4, '--teal-light');
    drawDot(ctx, j[0].x, j[0].y, 9, '--teal-light');
  }

  function dp_renderTrail() {
    const n = trailCount;
    if (n === 0) return;
    const ps = pivotScreen(), pz = zoom();
    trailCtx.save();
    trailCtx.lineCap = 'butt';
    trailCtx.lineWidth = 6;
    trailCtx.shadowBlur = 10;
    for (let i = 0; i < n; i++) {
      const idx = (trailHead - n + i + TRAIL_LEN) % TRAIL_LEN;
      const entry = trailBuf[idx];
      const alpha = Math.pow((i + 1) / n, 2) * 0.5;
      const j = dp_trailJoints(entry, ps, pz);
      trailCtx.globalAlpha = alpha;
      trailCtx.strokeStyle = _c('--teal-light');
      trailCtx.shadowColor = _c('--teal-dark');
      trailCtx.beginPath(); trailCtx.moveTo(j[0].x, j[0].y); trailCtx.lineTo(j[1].x, j[1].y); trailCtx.stroke();
      trailCtx.strokeStyle = _c('--pink-light');
      trailCtx.shadowColor = _c('--pink-dark');
      trailCtx.beginPath(); trailCtx.moveTo(j[1].x, j[1].y); trailCtx.lineTo(j[2].x, j[2].y); trailCtx.stroke();
    }
    trailCtx.restore();
  }

  /* ════════════════════════════════════════════
     CHAIN (N-pendulum)
     ════════════════════════════════════════════ */
  const CH_N       = 20;
  let CH_DAMPING = 0.05;  // mild dissipation (tunable via slider)
  const CH_G = G * CH_N;   // effective g for chain — compensates for short segment lengths
  // Each segment: dimensionless length = 1/N, mass = 1/N
  const ch_l = 1 / CH_N;
  const ch_m = 1 / CH_N;

  /* ── Scratch buffers (preallocated, reused each frame) ── */
  const _N  = CH_N;
  const _scrM   = new Float64Array(_N * _N);          // ch_buildSystem M (225)
  const _scrF   = new Float64Array(_N);               // ch_buildSystem f (15)
  const _scrA   = new Float64Array(_N * (_N + 1));    // ch_solveSystem augmented (240)
  const _scrX   = new Float64Array(_N);               // ch_solveSystem result (15)
  // RK4 stage buffers (shared by ch_rk4step, ch_rk4step_head, ch_rk4step_tail)
  const _rk_th2 = new Float64Array(_N);
  const _rk_w2  = new Float64Array(_N);
  const _rk_th3 = new Float64Array(_N);
  const _rk_w3  = new Float64Array(_N);
  const _rk_th4 = new Float64Array(_N);
  const _rk_w4  = new Float64Array(_N);
  // ch_derivs dth return buffer
  const _scrDth = new Float64Array(_N);
  /* ── Cached CSS values ── */
  let _white = null;

  let ch_th = new Float64Array(CH_N);   // absolute angles
  let ch_w  = new Float64Array(CH_N);   // angular velocities

  let ch_grabbing  = false;
  let ch_grabRod   = -1;   // 0-indexed rod being held
  let ch_cursorX   = 0;    // cursor world position (canvas px) during grab
  let ch_cursorY   = 0;

  /*
   * N-pendulum EOM using the standard Lagrangian (absolute angles).
   * The mass matrix M is symmetric with:
   *   M[i][j] = (sum of masses from max(i,j) to N-1) * l_i * l_j * cos(th_i - th_j)
   * Solved via Gaussian elimination each step.
   *
   * For equal masses and lengths (ch_m, ch_l):
   *   M[i][j] = ch_m * (N - max(i,j)) * ch_l^2 * cos(th[i]-th[j])
   *
   * RHS f[i] = -g * l_i * sin(th[i]) * (sum of masses from i to N-1)
   *            - l_i * sum_{j≠i} M_ij/l_i * w[j]^2 * sin(th[i]-th[j])
   *            - damping * w[i]
   *
   * We solve M * alpha = f for angular accelerations alpha.
   */
  function ch_buildSystem(th_, w_) {
    const N = CH_N, l = ch_l, m = ch_m;
    // Tail mass: tailM[i] = m*(N-i)
    // M[i][j] = tailM[max(i,j)] * l^2 * cos(th[i]-th[j])
    const M = _scrM;
    const f = _scrF;

    for (let i = 0; i < N; i++) {
      for (let j = 0; j < N; j++) {
        const tail = m * (N - Math.max(i, j));
        M[i*N+j] = tail * l * l * Math.cos(th_[i] - th_[j]);
      }
      // RHS: gravity term
      const tail_i = m * (N - i);
      f[i] = -CH_G * l * tail_i * Math.sin(th_[i]) - CH_DAMPING * w_[i];
      // Centripetal cross terms
      for (let j = 0; j < N; j++) {
        if (j === i) continue;
        const tail_ij = m * (N - Math.max(i, j));
        f[i] -= tail_ij * l * l * w_[j] * w_[j] * Math.sin(th_[i] - th_[j]);
      }
    }
    return { M, f };
  }

  // Gaussian elimination with partial pivoting, returns alpha array
  function ch_solveSystem(M, f) {
    const N = CH_N;
    // Augmented matrix [M | f]
    const A = _scrA;
    for (let i = 0; i < N; i++) {
      for (let j = 0; j < N; j++) A[i*(N+1)+j] = M[i*N+j];
      A[i*(N+1)+N] = f[i];
    }
    for (let col = 0; col < N; col++) {
      // Partial pivot
      let maxVal = Math.abs(A[col*(N+1)+col]), maxRow = col;
      for (let row = col+1; row < N; row++) {
        const v = Math.abs(A[row*(N+1)+col]);
        if (v > maxVal) { maxVal = v; maxRow = row; }
      }
      if (maxRow !== col) {
        for (let j = 0; j <= N; j++) {
          const tmp = A[col*(N+1)+j];
          A[col*(N+1)+j] = A[maxRow*(N+1)+j];
          A[maxRow*(N+1)+j] = tmp;
        }
      }
      const pivot = A[col*(N+1)+col];
      if (Math.abs(pivot) < 1e-14) continue;
      for (let row = col+1; row < N; row++) {
        const factor = A[row*(N+1)+col] / pivot;
        for (let j = col; j <= N; j++) A[row*(N+1)+j] -= factor * A[col*(N+1)+j];
      }
    }
    // Back substitution
    const x = _scrX;
    for (let i = N-1; i >= 0; i--) {
      x[i] = A[i*(N+1)+N];
      for (let j = i+1; j < N; j++) x[i] -= A[i*(N+1)+j] * x[j];
      x[i] /= A[i*(N+1)+i];
    }
    return x;
  }

  function ch_derivs(th_, w_) {
    const { M, f } = ch_buildSystem(th_, w_);
    const alpha = ch_solveSystem(M, f);
    // Copy w_ into preallocated buffer (avoids allocation per RK4 stage)
    _scrDth.set(w_);
    return { dth: _scrDth, dalpha: alpha };
  }

  function wrapAngle(a) {
    a = a % (2 * Math.PI);
    if (a >  Math.PI) a -= 2 * Math.PI;
    if (a < -Math.PI) a += 2 * Math.PI;
    return a;
  }

  function ch_rk4step(dt_) {
    const N = CH_N;
    const th = ch_th, w = ch_w;

    const k1 = ch_derivs(th, w);
    // Save k1 results before ch_derivs overwrites shared buffers on next call
    const k1dth = new Float64Array(k1.dth), k1da = new Float64Array(k1.dalpha);

    for (let i = 0; i < N; i++) { _rk_th2[i] = th[i]+0.5*dt_*k1dth[i]; _rk_w2[i] = w[i]+0.5*dt_*k1da[i]; }
    const k2 = ch_derivs(_rk_th2, _rk_w2);
    const k2dth = new Float64Array(k2.dth), k2da = new Float64Array(k2.dalpha);

    for (let i = 0; i < N; i++) { _rk_th3[i] = th[i]+0.5*dt_*k2dth[i]; _rk_w3[i] = w[i]+0.5*dt_*k2da[i]; }
    const k3 = ch_derivs(_rk_th3, _rk_w3);
    const k3dth = new Float64Array(k3.dth), k3da = new Float64Array(k3.dalpha);

    for (let i = 0; i < N; i++) { _rk_th4[i] = th[i]+dt_*k3dth[i]; _rk_w4[i] = w[i]+dt_*k3da[i]; }
    const k4 = ch_derivs(_rk_th4, _rk_w4);

    for (let i = 0; i < N; i++) {
      ch_th[i] = wrapAngle(ch_th[i] + (dt_/6)*(k1dth[i]+2*k2dth[i]+2*k3dth[i]+k4.dth[i]));
      ch_w[i]  += (dt_/6)*(k1da[i]+2*k2da[i]+2*k3da[i]+k4.dalpha[i]);
    }
  }

  /*
   * Grab physics: when rod k (0-indexed) is grabbed, split chain into:
   *   - "before" sub-chain: rods 0..k-1, with fixed top pivot and bottom
   *     boundary at the grabbed joint position — evolve freely.
   *   - rod k: pinned at cursor.
   *   - "after" sub-chain: rods k+1..N-1, with top pivot = joint k — evolve freely.
   *
   * Implementation: run the full N-body EOM but override rod k's angle to
   * track the cursor each sub-step. This naturally allows the neighbouring
   * rods to respond dynamically, with ch_w[k] zeroed each step.
   */
  /*
   * FABRIK IK for rods 0..k, with root fixed at pivot and tip pulled to cursor.
   * Returns joint positions as [{x,y}] of length k+2 (pivot + k+1 joints).
   * Rods k+1..N-1 are untouched — their angles evolve freely via EOM.
   *
   * After solving, we back-convert joint positions to angles and write
   * ch_th[0..k], zeroing ch_w[0..k].
   */
  function ch_applyGrabConstraint() {
    if (!ch_grabbing || ch_grabRod < 0) return;
    const k = ch_grabRod;
    const z = zoom();
    const lPx = (totalPx() / CH_N) * z;
    const { px, py } = pivotScreen();

    // Build current joint positions for rods 0..k (k+2 points)
    const pts = [];
    let cx = px, cy = py;
    pts.push({ x: cx, y: cy });
    for (let i = 0; i <= k; i++) {
      cx += lPx * Math.sin(ch_th[i]);
      cy += lPx * Math.cos(ch_th[i]);
      pts.push({ x: cx, y: cy });
    }

    // Target: cursor, clamped to reachable disk of radius (k+1)*lPx from pivot
    const maxReach = (k + 1) * lPx;
    let tx = ch_cursorX, ty = ch_cursorY;
    const distToPivot = Math.hypot(tx - px, ty - py);
    if (distToPivot > maxReach) {
      const angle = Math.atan2(ty - py, tx - px);
      tx = px + maxReach * Math.cos(angle);
      ty = py + maxReach * Math.sin(angle);
    }

    // FABRIK iterations
    const ITER = 10;
    const TOL  = 0.5; // px
    for (let iter = 0; iter < ITER; iter++) {
      // Forward pass: pull tip to target
      pts[k + 1].x = tx; pts[k + 1].y = ty;
      for (let i = k; i >= 0; i--) {
        const dx = pts[i].x - pts[i+1].x, dy = pts[i].y - pts[i+1].y;
        const d = Math.hypot(dx, dy);
        if (d < 1e-9) continue;
        pts[i].x = pts[i+1].x + (dx / d) * lPx;
        pts[i].y = pts[i+1].y + (dy / d) * lPx;
      }
      // Backward pass: re-anchor root at pivot
      pts[0].x = px; pts[0].y = py;
      for (let i = 0; i <= k; i++) {
        const dx = pts[i+1].x - pts[i].x, dy = pts[i+1].y - pts[i].y;
        const d = Math.hypot(dx, dy);
        if (d < 1e-9) continue;
        pts[i+1].x = pts[i].x + (dx / d) * lPx;
        pts[i+1].y = pts[i].y + (dy / d) * lPx;
      }
      if (Math.hypot(pts[k+1].x - tx, pts[k+1].y - ty) < TOL) break;
    }

    // Convert joint positions back to absolute angles, write into ch_th[0..k]
    for (let i = 0; i <= k; i++) {
      const dx = pts[i+1].x - pts[i].x, dy = pts[i+1].y - pts[i].y;
      ch_th[i] = Math.atan2(dx, dy);
      ch_w[i]  = 0;
    }
  }

  // RK4 step for rods 0..endIdx as independent sub-chain (top pivot fixed)
  function ch_rk4step_head(endIdx, dt_) {
    const M_ = endIdx + 1;
    if (M_ <= 0) return;

    // th0/w0 hold initial state for this step
    const th0 = new Float64Array(M_), w0 = new Float64Array(M_);
    for (let i = 0; i < M_; i++) { th0[i] = ch_th[i]; w0[i] = ch_w[i]; }

    function derivs(th_, w_) {
      const { M, f } = ch_buildSubSystem(th_, w_, M_);
      const alpha = ch_solveSubSystem(M, f, M_);
      // Copy w_ into _scrDth (reuse shared buffer)
      _scrDth.set(w_);
      return { dth: _scrDth, dalpha: alpha };
    }

    const k1 = derivs(th0, w0);
    const k1dth = new Float64Array(k1.dth.subarray(0, M_));
    const k1da  = new Float64Array(k1.dalpha.subarray(0, M_));

    for (let i = 0; i < M_; i++) { _rk_th2[i] = th0[i]+0.5*dt_*k1dth[i]; _rk_w2[i] = w0[i]+0.5*dt_*k1da[i]; }
    const k2 = derivs(_rk_th2, _rk_w2);
    const k2dth = new Float64Array(k2.dth.subarray(0, M_));
    const k2da  = new Float64Array(k2.dalpha.subarray(0, M_));

    for (let i = 0; i < M_; i++) { _rk_th3[i] = th0[i]+0.5*dt_*k2dth[i]; _rk_w3[i] = w0[i]+0.5*dt_*k2da[i]; }
    const k3 = derivs(_rk_th3, _rk_w3);
    const k3dth = new Float64Array(k3.dth.subarray(0, M_));
    const k3da  = new Float64Array(k3.dalpha.subarray(0, M_));

    for (let i = 0; i < M_; i++) { _rk_th4[i] = th0[i]+dt_*k3dth[i]; _rk_w4[i] = w0[i]+dt_*k3da[i]; }
    const k4 = derivs(_rk_th4, _rk_w4);

    for (let i = 0; i < M_; i++) {
      ch_th[i] = wrapAngle(ch_th[i] + (dt_/6)*(k1dth[i]+2*k2dth[i]+2*k3dth[i]+k4.dth[i]));
      ch_w[i]  += (dt_/6)*(k1da[i]+2*k2da[i]+2*k3da[i]+k4.dalpha[i]);
    }
  }

  function ch_stepPhysics() {
    const dt = DT / CH_SUBSTEPS;
    if (!ch_grabbing || ch_grabRod < 0) {
      // Free evolution: full N-body RK4
      for (let i = 0; i < CH_SUBSTEPS; i++) ch_rk4step(dt);
    } else {
      const k = ch_grabRod;
      for (let i = 0; i < CH_SUBSTEPS; i++) {
        // Inner chain 0..k: integrate freely, then FABRIK to enforce cursor constraint
        const thBefore = new Float64Array(k+1);
        for (let j = 0; j <= k; j++) thBefore[j] = ch_th[j];

        ch_rk4step_head(k, dt);

        // FABRIK to pull joint k+1 (tip of inner chain) back to cursor
        ch_applyGrabConstraint();

        // Velocity correction: finite difference from position change
        for (let j = 0; j <= k; j++) {
          ch_w[j] = (ch_th[j] - thBefore[j]) / dt;
        }

        // Outer chain k+1..N-1: integrate as independent sub-chain
        ch_rk4step_tail(k + 1, dt);
      }
    }
  }

  /*
   * Build EOM for a sub-chain of length M starting at startIdx,
   * treating joint startIdx-1 as a fixed pivot (no coupling to inner rods).
   * th_ and w_ are slices of length M (indices 0..M-1 = global startIdx..N-1).
   */
  function ch_buildSubSystem(th_, w_, M_) {
    const m = ch_m, l = ch_l;
    // Reuse _scrM and _scrF (safe: ch_buildSubSystem is never called concurrently with ch_buildSystem)
    const M = _scrM;
    const f = _scrF;
    for (let i = 0; i < M_; i++) {
      for (let j = 0; j < M_; j++) {
        const tail = m * (M_ - Math.max(i, j));
        M[i*M_+j] = tail * l * l * Math.cos(th_[i] - th_[j]);
      }
      const tail_i = m * (M_ - i);
      f[i] = -CH_G * l * tail_i * Math.sin(th_[i]) - CH_DAMPING * w_[i];
      for (let j = 0; j < M_; j++) {
        if (j === i) continue;
        const tail_ij = m * (M_ - Math.max(i, j));
        f[i] -= tail_ij * l * l * w_[j] * w_[j] * Math.sin(th_[i] - th_[j]);
      }
    }
    return { M, f };
  }

  function ch_solveSubSystem(M, f, M_) {
    // Reuse _scrA and _scrX (safe: ch_solveSubSystem is never called concurrently with ch_solveSystem)
    const A = _scrA;
    for (let i = 0; i < M_; i++) {
      for (let j = 0; j < M_; j++) A[i*(M_+1)+j] = M[i*M_+j];
      A[i*(M_+1)+M_] = f[i];
    }
    for (let col = 0; col < M_; col++) {
      let maxVal = Math.abs(A[col*(M_+1)+col]), maxRow = col;
      for (let row = col+1; row < M_; row++) {
        const v = Math.abs(A[row*(M_+1)+col]);
        if (v > maxVal) { maxVal = v; maxRow = row; }
      }
      if (maxRow !== col) {
        for (let j = 0; j <= M_; j++) {
          const tmp = A[col*(M_+1)+j]; A[col*(M_+1)+j] = A[maxRow*(M_+1)+j]; A[maxRow*(M_+1)+j] = tmp;
        }
      }
      const pivot = A[col*(M_+1)+col];
      if (Math.abs(pivot) < 1e-14) continue;
      for (let row = col+1; row < M_; row++) {
        const factor = A[row*(M_+1)+col] / pivot;
        for (let j = col; j <= M_; j++) A[row*(M_+1)+j] -= factor * A[col*(M_+1)+j];
      }
    }
    const x = _scrX;
    for (let i = M_-1; i >= 0; i--) {
      x[i] = A[i*(M_+1)+M_];
      for (let j = i+1; j < M_; j++) x[i] -= A[i*(M_+1)+j] * x[j];
      x[i] /= A[i*(M_+1)+i];
    }
    return x;
  }

  // RK4 step for rods startIdx..N-1 as an independent sub-chain
  function ch_rk4step_tail(startIdx, dt_) {
    const N = CH_N;
    const M_ = N - startIdx;
    if (M_ <= 0) return;

    // Extract sub-chain state
    const th0 = new Float64Array(M_), w0 = new Float64Array(M_);
    for (let i = 0; i < M_; i++) { th0[i] = ch_th[startIdx+i]; w0[i] = ch_w[startIdx+i]; }

    function derivs(th_, w_) {
      const { M, f } = ch_buildSubSystem(th_, w_, M_);
      const alpha = ch_solveSubSystem(M, f, M_);
      // Copy w_ into _scrDth (reuse shared buffer)
      _scrDth.set(w_);
      return { dth: _scrDth, dalpha: alpha };
    }

    const k1 = derivs(th0, w0);
    const k1dth = new Float64Array(k1.dth.subarray(0, M_));
    const k1da  = new Float64Array(k1.dalpha.subarray(0, M_));

    for (let i = 0; i < M_; i++) { _rk_th2[i] = th0[i]+0.5*dt_*k1dth[i]; _rk_w2[i] = w0[i]+0.5*dt_*k1da[i]; }
    const k2 = derivs(_rk_th2, _rk_w2);
    const k2dth = new Float64Array(k2.dth.subarray(0, M_));
    const k2da  = new Float64Array(k2.dalpha.subarray(0, M_));

    for (let i = 0; i < M_; i++) { _rk_th3[i] = th0[i]+0.5*dt_*k2dth[i]; _rk_w3[i] = w0[i]+0.5*dt_*k2da[i]; }
    const k3 = derivs(_rk_th3, _rk_w3);
    const k3dth = new Float64Array(k3.dth.subarray(0, M_));
    const k3da  = new Float64Array(k3.dalpha.subarray(0, M_));

    for (let i = 0; i < M_; i++) { _rk_th4[i] = th0[i]+dt_*k3dth[i]; _rk_w4[i] = w0[i]+dt_*k3da[i]; }
    const k4 = derivs(_rk_th4, _rk_w4);

    for (let i = 0; i < M_; i++) {
      ch_th[startIdx+i] = wrapAngle(ch_th[startIdx+i] + (dt_/6)*(k1dth[i]+2*k2dth[i]+2*k3dth[i]+k4.dth[i]));
      ch_w[startIdx+i]  += (dt_/6)*(k1da[i]+2*k2da[i]+2*k3da[i]+k4.dalpha[i]);
    }
  }

  // Compute joint screen positions from angles
  function ch_joints(th_, _pivot, _zoom) {
    const { px, py } = _pivot || pivotScreen();
    const z = (_zoom !== undefined) ? _zoom : zoom();
    const lPx = (totalPx() / CH_N) * z;
    const pts = [{ x: px, y: py }];
    let cx = px, cy = py;
    for (let i = 0; i < CH_N; i++) {
      cx += lPx * Math.sin(th_[i]);
      cy += lPx * Math.cos(th_[i]);
      pts.push({ x: cx, y: cy });
    }
    return pts;
  }

  function ch_recordTrail() {
    trailTick++;
    if (trailTick % TRAIL_INTERVAL !== 0) return;
    trailBuf[trailHead] = { th: Float64Array.from(ch_th) };
    trailHead = (trailHead + 1) % TRAIL_LEN;
    if (trailCount < TRAIL_LEN) trailCount++;
  }

  function ch_init() {
    ch_th = new Float64Array(CH_N);
    ch_w  = new Float64Array(CH_N);
    // Small initial offset on last rod only
    ch_th[CH_N - 1] = 0.25;
    trailHead = 0; trailCount = 0; trailTick = 0;
    ch_grabbing = false; ch_grabRod = -1;
    // Force color cache rebuild in case CH_N changed
    _ch_colors = null; _ch_colorLight = null; _ch_colorDark = null;
  }

  // Rod colour: interpolate teal→pink by index (resolved lazily on first use)
  let _ch_colors = null;
  let _ch_colorLight = null;  // array of 15 precomputed light color strings
  let _ch_colorDark  = null;  // array of 15 precomputed dark color strings
  function ch_ensureColors() {
    if (_ch_colors) return;
    _ch_colors = {
      TL: _rgb('--teal-light'), PL: _rgb('--pink-light'),
      TD: _rgb('--teal-dark'),  PD: _rgb('--pink-dark'),
    };
    const { TL, PL, TD, PD } = _ch_colors;
    _ch_colorLight = new Array(CH_N);
    _ch_colorDark  = new Array(CH_N);
    for (let i = 0; i < CH_N; i++) {
      const t = i / (CH_N - 1);
      _ch_colorLight[i] = `rgb(${Math.round(TL[0]+(PL[0]-TL[0])*t)},${Math.round(TL[1]+(PL[1]-TL[1])*t)},${Math.round(TL[2]+(PL[2]-TL[2])*t)})`;
      _ch_colorDark[i]  = `rgb(${Math.round(TD[0]+(PD[0]-TD[0])*t)},${Math.round(TD[1]+(PD[1]-TD[1])*t)},${Math.round(TD[2]+(PD[2]-TD[2])*t)})`;
    }
  }

  function ch_rodColor(i, dark) {
    ch_ensureColors();
    return dark ? _ch_colorDark[i] : _ch_colorLight[i];
  }

  function ch_render() {
    const j = ch_joints(ch_th);
    // Glow pass first (all rods), then white core pass on top
    ctx.save();
    ctx.lineCap = 'round';
    for (let i = 0; i < CH_N; i++) {
      ctx.shadowColor = ch_rodColor(i, true);
      ctx.shadowBlur  = 14;
      ctx.strokeStyle = ch_rodColor(i, false);
      ctx.lineWidth   = 4 * 2.5;
      ctx.beginPath(); ctx.moveTo(j[i].x, j[i].y); ctx.lineTo(j[i+1].x, j[i+1].y); ctx.stroke();
    }
    ctx.shadowBlur  = 0;
    ctx.strokeStyle = (_white || (_white = _c('--white')));
    ctx.lineWidth   = 4;
    for (let i = 0; i < CH_N; i++) {
      ctx.beginPath(); ctx.moveTo(j[i].x, j[i].y); ctx.lineTo(j[i+1].x, j[i+1].y); ctx.stroke();
    }
    ctx.restore();
    drawDot(ctx, j[0].x, j[0].y, 9, '--teal-light');
  }

  function ch_renderTrail() {
    const n = trailCount;
    if (n === 0) return;
    ch_ensureColors();
    const ps = pivotScreen(), pz = zoom();
    trailCtx.save();
    trailCtx.lineCap = 'butt';
    trailCtx.lineWidth = 6;
    trailCtx.shadowBlur = 10;
    for (let fi = 0; fi < n; fi++) {
      const idx = (trailHead - n + fi + TRAIL_LEN) % TRAIL_LEN;
      const entry = trailBuf[idx];
      const alpha = Math.pow((fi + 1) / n, 2) * 0.45;
      const j = ch_joints(entry.th, ps, pz);
      trailCtx.globalAlpha = alpha;
      for (let i = 0; i < CH_N; i++) {
        trailCtx.strokeStyle = _ch_colorLight[i];
        trailCtx.shadowColor = _ch_colorDark[i];
        trailCtx.beginPath(); trailCtx.moveTo(j[i].x, j[i].y); trailCtx.lineTo(j[i+1].x, j[i+1].y); trailCtx.stroke();
      }
    }
    trailCtx.restore();
  }

  function ch_onPointerDown(x, y) {
    const j = ch_joints(ch_th);
    for (let i = CH_N - 1; i >= 0; i--) {
      const d = distToSegment(x, y, j[i].x, j[i].y, j[i+1].x, j[i+1].y);
      if (d <= GRAB_THRESH) {
        ch_grabbing = true;
        ch_grabRod  = i;
        ch_cursorX  = x;
        ch_cursorY  = y;
        ch_w[i] = 0;
        return true;
      }
    }
    return false;
  }

  function ch_onPointerMove(x, y) {
    if (!ch_grabbing || ch_grabRod < 0) return;
    ch_cursorX = x;
    ch_cursorY = y;
    // Apply constraint immediately so paused view updates
    ch_applyGrabConstraint();
  }

  function ch_onPointerUp() { ch_grabbing = false; ch_grabRod = -1; }

  /* ════════════════════════════════════════════
     SHARED DRAWING HELPERS
     ════════════════════════════════════════════ */
  function drawRod(c, x1, y1, x2, y2, glowVar, lineW, glowW, alpha) {
    drawRodColored(c, x1, y1, x2, y2, _c(glowVar), _c(glowVar), lineW, glowW, alpha);
  }

  function drawRodColored(c, x1, y1, x2, y2, glowColor, glowDark, lineW, glowW, alpha) {
    c.save();
    c.globalAlpha = alpha;
    c.lineCap = 'round';
    c.shadowColor = glowDark;
    c.shadowBlur  = glowW;
    c.strokeStyle = glowColor;
    c.lineWidth   = lineW * 2.5;
    c.beginPath(); c.moveTo(x1, y1); c.lineTo(x2, y2); c.stroke();
    c.shadowBlur  = 0;
    c.strokeStyle = (_white || (_white = _c('--white')));
    c.lineWidth   = lineW;
    c.beginPath(); c.moveTo(x1, y1); c.lineTo(x2, y2); c.stroke();
    c.restore();
  }

  function drawDot(c, x, y, r, glowVar) {
    drawDotColored(c, x, y, r, _c(glowVar));
  }

  function drawDotColored(c, x, y, r, color) {
    c.save();
    c.shadowColor = color;
    c.shadowBlur  = 8;
    c.fillStyle   = (_white || (_white = _c('--white')));
    c.beginPath(); c.arc(x, y, r, 0, Math.PI * 2); c.fill();
    c.restore();
  }

  function trailStroke(c, x1, y1, x2, y2, colorVar, glowVar, alpha) {
    trailStrokeColored(c, x1, y1, x2, y2, _c(colorVar), _c(glowVar), alpha);
  }

  function trailStrokeColored(c, x1, y1, x2, y2, color, glow, alpha) {
    c.save();
    c.globalAlpha = alpha;
    c.shadowColor = glow;
    c.shadowBlur  = 10;
    c.strokeStyle = color;
    c.lineWidth   = 6;
    c.lineCap     = 'butt';
    c.beginPath(); c.moveTo(x1, y1); c.lineTo(x2, y2); c.stroke();
    c.restore();
  }

  function distToSegment(px_, py_, ax, ay, bx, by) {
    const dx = bx - ax, dy = by - ay;
    const lenSq = dx*dx + dy*dy;
    if (lenSq === 0) return Math.hypot(px_-ax, py_-ay);
    const t = Math.max(0, Math.min(1, ((px_-ax)*dx + (py_-ay)*dy) / lenSq));
    return Math.hypot(px_ - (ax + t*dx), py_ - (ay + t*dy));
  }

  const GRAB_THRESH = 20;

  /* ════════════════════════════════════════════
     SHARED LOOP & RENDER
     ════════════════════════════════════════════ */
  function renderTrail() {
    if (!trailCvs || trailCount === 0) return;
    trailCtx.clearRect(0, 0, W, H);
    if (currentMode === 'double') dp_renderTrail();
    else if (currentMode === 'chain') ch_renderTrail();
    ctx.drawImage(trailCvs, 0, 0);
  }

  function render() {
    ctx.fillStyle = _c('--black');
    ctx.fillRect(0, 0, W, H);
    renderTrail();
    if (currentMode === 'double') dp_render();
    else if (currentMode === 'chain') ch_render();
  }

  function loop() {
    const grabbing = (currentMode === 'double') ? dp_grabbing : ch_grabbing;
    if (running) {
      if (currentMode === 'double') { dp_stepPhysics(); dp_recordTrail(); }
      else if (currentMode === 'chain') { ch_stepPhysics(); ch_recordTrail(); }
    }
    render();
    frameId = requestAnimationFrame(loop);
  }

  /* ── Mode switching ── */
  function setMode(m) {
    currentMode = m;
    ['double', 'chain'].forEach(id => {
      const btn = document.getElementById('pend-btn-' + id);
      if (btn) btn.classList.toggle('active', id === m);
    });
    const dblEl = document.getElementById('pend-double-controls');
    const chEl  = document.getElementById('pend-chain-controls');
    if (dblEl) dblEl.style.display = (m === 'double') ? 'flex' : 'none';
    if (chEl)  chEl.style.display  = (m === 'chain')  ? 'flex' : 'none';
    if (m === 'double') dp_init();
    else if (m === 'chain') ch_init();
  }

  /* ── Pointer dispatch ── */
  function canvasCoords(evt) {
    const rect = canvas.getBoundingClientRect();
    const sx = canvas.width / rect.width, sy = canvas.height / rect.height;
    const src = evt.touches ? evt.touches[0] : evt;
    return { x: (src.clientX - rect.left) * sx, y: (src.clientY - rect.top) * sy };
  }

  function onPointerDown(evt) {
    const { x, y } = canvasCoords(evt);
    let hit = false;
    if (currentMode === 'double') hit = dp_onPointerDown(x, y);
    else if (currentMode === 'chain') hit = ch_onPointerDown(x, y);
    if (hit) evt.preventDefault();
  }

  function onPointerMove(evt) {
    const { x, y } = canvasCoords(evt);
    if (currentMode === 'double') dp_onPointerMove(x, y);
    else if (currentMode === 'chain') ch_onPointerMove(x, y);
    const grabbing = (currentMode === 'double') ? dp_grabbing : ch_grabbing;
    if (grabbing) { if (!running) render(); evt.preventDefault(); }
  }

  function onPointerUp() {
    if (currentMode === 'double') dp_onPointerUp();
    else if (currentMode === 'chain') ch_onPointerUp();
  }

  function onWheel(evt) {
    evt.preventDefault();
    zoomT = Math.max(0, Math.min(1, zoomT + (evt.deltaY > 0 ? 0.06 : -0.06)));
  }

  /* ── Shell ── */
  const shell = new AppletShell({
    id:    'pend',
    title: 'Pendulums',
    gap:   0,

    headerBtns: `<button class="applet-shell-header-btn" onclick="pendReset()">Reset</button><button class="applet-shell-header-btn" id="pend-pause-btn" onclick="pendTogglePause()">Pause</button>`,

    ctrlHTML: `
      <div id="pend-scrollable">

        <div class="applet-shell-ctrl-section">
          <div class="applet-shell-btn-row">
            <button class="applet-shell-btn active" id="pend-btn-double" onclick="pendSetMode('double')">Double</button>
            <button class="applet-shell-btn"        id="pend-btn-chain"  onclick="pendSetMode('chain')">Chain</button>
            <button class="applet-shell-btn" disabled style="opacity:0.4;">Wave</button>
          </div>
        </div>

        <div id="pend-double-controls" style="display:flex; flex-direction:column;">
          <div class="applet-shell-ctrl-section">
            <div class="applet-shell-ctrl-title">Length ratio</div>
            <div class="applet-shell-slider-row">
              <span class="applet-shell-side">Short/Long</span>
              <input type="range" id="pend-len-ratio" min="0.1" max="0.9" step="0.01" value="0.5">
              <span class="applet-shell-side">Long/Short</span>
              <span class="applet-shell-val" id="pend-len-ratio-val">0.50</span>
            </div>
          </div>
          <div class="applet-shell-ctrl-section">
            <div class="applet-shell-ctrl-title">Mass ratio</div>
            <div class="applet-shell-slider-row">
              <span class="applet-shell-side">Light/Heavy</span>
              <input type="range" id="pend-mass-ratio" min="0.1" max="0.9" step="0.01" value="0.5">
              <span class="applet-shell-side">Heavy/Light</span>
              <span class="applet-shell-val" id="pend-mass-ratio-val">0.50</span>
            </div>
          </div>
          <div class="applet-shell-ctrl-section">
            <div class="applet-shell-ctrl-title">Dissipation</div>
            <div class="applet-shell-slider-row">
              <span class="applet-shell-side">None</span>
              <input type="range" id="pend-damping" min="0" max="1" step="0.01" value="0">
              <span class="applet-shell-side">High</span>
              <span class="applet-shell-val" id="pend-damping-val">0.00</span>
            </div>
          </div>
        </div>

        <div id="pend-chain-controls" style="display:none; flex-direction:column;">
          <div class="applet-shell-ctrl-section">
            <div class="applet-shell-ctrl-title">Damping</div>
            <div class="applet-shell-slider-row">
              <span class="applet-shell-side">0</span>
              <input type="range" id="pend-ch-damping" min="0" max="0.2" step="0.001" value="0.05">
              <span class="applet-shell-side">0.2</span>
              <span class="applet-shell-val" id="pend-ch-damping-val">0.050</span>
            </div>
          </div>
        </div>

      </div>
    `,

    onOpen: function ({ canvas: c, S: s, W: ww, H: hh }) {
      canvas = c; ctx = canvas.getContext('2d');
      S = s; W = ww || s; H = hh || s;

      if (!trailCvs) { trailCvs = document.createElement('canvas'); trailCtx = trailCvs.getContext('2d'); }
      trailCvs.width = W; trailCvs.height = H;

      dp_lenRatio  = parseFloat(document.getElementById('pend-len-ratio').value);
      dp_massRatio = parseFloat(document.getElementById('pend-mass-ratio').value);
      dp_damping   = parseFloat(document.getElementById('pend-damping').value);

      setMode('double');

      canvas.addEventListener('mousedown',  onPointerDown);
      canvas.addEventListener('mousemove',  onPointerMove);
      canvas.addEventListener('mouseup',    onPointerUp);
      canvas.addEventListener('mouseleave', onPointerUp);
      canvas.addEventListener('touchstart', onPointerDown, { passive: false });
      canvas.addEventListener('touchmove',  onPointerMove, { passive: false });
      canvas.addEventListener('touchend',   onPointerUp);
      canvas.addEventListener('wheel',      onWheel,       { passive: false });

      running = true;
      const pb = document.getElementById('pend-pause-btn');
      if (pb) { pb.textContent = 'Pause'; pb.classList.remove('active'); }
      if (!frameId) frameId = requestAnimationFrame(loop);
    },

    onClose: function () {
      running = false;
      if (frameId) { cancelAnimationFrame(frameId); frameId = null; }
    },

    onResize: function ({ S: s, W: ww, H: hh }) {
      S = s; W = ww || s; H = hh || s;
      if (trailCvs) { trailCvs.width = W; trailCvs.height = H; }
      trailHead = 0; trailCount = 0;
    },
  });

  /* ── Global entry points ── */
  window.pendOpen  = () => shell.open();
  window.pendClose = () => shell.close();

  window.pendReset = function () {
    if (currentMode === 'double') dp_init();
    else if (currentMode === 'chain') ch_init();
  };

  window.pendTogglePause = function () {
    running = !running;
    if (!running) { trailHead = 0; trailCount = 0; }
    const pb = document.getElementById('pend-pause-btn');
    if (pb) { pb.textContent = running ? 'Pause' : 'Resume'; pb.classList.toggle('active', !running); }
  };

  window.pendSetMode = function (m) { setMode(m); };

  /* ── Slider listeners ── */
  document.getElementById('pend-len-ratio').addEventListener('input', function () {
    dp_lenRatio = parseFloat(this.value);
    document.getElementById('pend-len-ratio-val').textContent = dp_lenRatio.toFixed(2);
    trailHead = 0; trailCount = 0;
  });
  document.getElementById('pend-mass-ratio').addEventListener('input', function () {
    dp_massRatio = parseFloat(this.value);
    document.getElementById('pend-mass-ratio-val').textContent = dp_massRatio.toFixed(2);
  });
  document.getElementById('pend-damping').addEventListener('input', function () {
    dp_damping = parseFloat(this.value);
    document.getElementById('pend-damping-val').textContent = dp_damping.toFixed(2);
  });

  document.getElementById('pend-ch-damping').addEventListener('input', function () {
    CH_DAMPING = parseFloat(this.value);
    document.getElementById('pend-ch-damping-val').textContent = CH_DAMPING.toFixed(3);
  });

})();
