// Run with node _tools/test_swarmalator.cjs.
'use strict';
const assert = require('node:assert/strict');
const { createModel, PRESETS, SOFTENING } = require('../assets/js/applets/swarmalator_applet.js');
function random(seed) {
  return () => {
    seed = (Math.imul(1664525, seed) + 1013904223) >>> 0;
    return seed / 4294967296;
  };
}
function finite(model) {
  assert(model.state.every(Number.isFinite));
  const o = model.observe();
  assert(o.r >= 0 && o.r <= 1 + 1e-12 && o.S >= 0 && o.S <= 1 + 1e-12);
  return o;
}
function distance(a, b) {
  let sum = 0;
  for (let i = 0; i < a.length; i++) {
    const d = i % 3 === 2 ? Math.atan2(Math.sin(a[i] - b[i]), Math.cos(a[i] - b[i])) : a[i] - b[i];
    sum += d * d;
  }
  return Math.sqrt(sum / a.length);
}

// A two-particle case has an exact phase solution when J=0 and softened r=1.
const pair = createModel(2, random(1));
const separation = Math.sqrt(1 - SOFTENING ** 2);
pair.state.set([-separation / 2, 0, 0, separation / 2, 0, 1]);
pair.advance(2, 0, 1);
const delta = pair.state[5] - pair.state[2];
assert(Math.abs(delta - 2 * Math.atan(Math.tan(0.5) * Math.exp(-2))) < 1e-5);
assert(Math.abs(pair.state[3] - pair.state[0] - separation) < 1e-12);
console.log('PASS analytic two-particle phase relaxation and spatial equilibrium');

// Reciprocal pair interactions conserve the spatial centroid; K=0 freezes phases.
const invariant = createModel(80, random(2));
const initial = invariant.state.slice(), center = invariant.observe();
invariant.advance(10, 1, 0);
const after = finite(invariant);
assert(Math.hypot(after.x - center.x, after.y - center.y) < 1e-12);
for (let i = 2; i < initial.length; i += 3) assert(Math.abs(initial[i] - invariant.state[i]) < 1e-12);
console.log('PASS centroid conservation and frozen phases at K=0');

// Second-order convergence against the same initial data at three step sizes.
const solutions = [0.01, 0.005, 0.0025].map(dt => {
  const m = createModel(20, random(3));
  m.advance(1, 0.5, -0.5, dt);
  return m.state;
});
const coarseError = distance(solutions[0], solutions[1]);
const fineError = distance(solutions[1], solutions[2]);
assert(coarseError / fineError > 3.2 && coarseError / fineError < 4.8);
console.log('PASS Heun convergence ratio', (coarseError / fineError).toFixed(3));

// Long trajectories for all published presets; the small finite swarm is checked
// for boundedness and emergent structure, not equality to a continuum picture.
for (const [key, p] of Object.entries(PRESETS)) {
  const m = createModel(80, random(4));
  m.advance(key === 'splintered' || key === 'active' ? 500 : 100, p.J, p.K);
  const o = finite(m);
  let radius = 0;
  for (let i = 0; i < m.n; i++) radius = Math.max(radius, Math.hypot(m.state[3 * i] - o.x, m.state[3 * i + 1] - o.y));
  assert(radius < 3);
  if (key === 'sync') assert(o.r > 0.99);
  if (key === 'async') assert(o.r < 0.15);
  if (key === 'wave') assert(o.S > 0.9);
  console.log('PASS preset', key, JSON.stringify({ r: o.r, S: o.S, radius }));
}

// Supported parameter corners and maximum population, including the weakest
// confining attraction, plus a close pair that exercises the regularization.
for (const J of [-0.9, 0, 1]) for (const K of [-1, 0, 1]) {
  const m = createModel(240, random(5));
  m.advance(10, J, K);
  finite(m);
  assert(Math.max(...m.state.filter((_, i) => i % 3 !== 2).map(Math.abs)) < 4);
}
const close = createModel(80, random(6));
close.state[3] = close.state[0] + 1e-10; close.state[4] = close.state[1];
close.advance(2, 1, -1);
finite(close);
console.log('PASS supported parameter grid, maximum population and close encounter');
