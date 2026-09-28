#!/usr/bin/env node
/* Measure MFS/assets/embedding-turn.js, which the page turns an embedding diagram with, against
   every published embedding figure, for the tests in test_build_mfs_data.py to hold.

     node _tools/embedding_turn_check.cjs

   prints one JSON object. For every view it gives, at the figure's own camera, the length of
   every line class and the area of every fill class as the page would draw them beside the
   published ones, and how far each label lands from its published place; at a sweep of other
   cameras, the scale, how far any drawn point or label falls outside the box, and any number
   that is not one; and what share of the lines that a turn round the axis must leave in place
   moves under one. For Schwarzschild it also gives how far out on each sheet a circle is seen
   and how far out one is hidden, from straight above and from straight below. */
'use strict';
const fs = require('fs');
const path = require('path');

const root = path.resolve(__dirname, '..');
const turn = require(path.join(root, 'MFS', 'assets', 'embedding-turn.js'));
const dir = path.join(root, 'MFS', 'assets', 'data', 'embedding');

function length(lines) {
  let s = 0;
  for (const P of lines) for (let i = 0; i + 1 < P.length; i++) s += Math.hypot(P[i + 1][0] - P[i][0], P[i + 1][1] - P[i][1]);
  return s;
}
function inside(p, ring) {
  let c = false;
  for (let i = 0, j = ring.length - 1; i < ring.length; j = i++) {
    if ((ring[i][1] > p[1]) !== (ring[j][1] > p[1]) &&
        p[0] < (ring[j][0] - ring[i][0]) * (p[1] - ring[i][1]) / (ring[j][1] - ring[i][1]) + ring[i][0]) c = !c;
  }
  return c;
}
// The area the even odd rule paints: each ring counts in or out by how many others enclose it.
function area(rings) {
  return rings.reduce((sum, ring, i) => {
    let a = 0;
    for (let k = 0, j = ring.length - 1; k < ring.length; j = k++) a += ring[j][0] * ring[k][1] - ring[k][0] * ring[j][1];
    const depth = rings.filter((other, j) => j !== i && inside(ring[0], other)).length;
    return sum + (depth % 2 ? -1 : 1) * Math.abs(a) / 2;
  }, 0);
}
function byClass(layers) {
  const lines = {}, fills = {};
  for (const L of layers) {
    if (L.kind === 'line') (lines[L['class']] = lines[L['class']] || []).push(L.points);
    if (L.kind === 'fill') (fills[L['class']] = fills[L['class']] || []).push(L.points, ...(L.holes || []));
  }
  return { lines, fills };
}
function distance(p, lines) {
  let best = Infinity;
  for (const P of lines) {
    for (let i = 0; i + 1 < P.length; i++) {
      const a = P[i], b = P[i + 1], dx = b[0] - a[0], dy = b[1] - a[1], L2 = dx * dx + dy * dy;
      const t = L2 ? Math.max(0, Math.min(1, ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / L2)) : 0;
      best = Math.min(best, Math.hypot(p[0] - a[0] - t * dx, p[1] - a[1] - t * dy));
    }
  }
  return best;
}
// The length of polylines X lying further than tol from polylines Y, measured at the middle of
// each segment, in both directions.
function astray(A, B, tol) {
  let far = 0;
  for (const [X, Y] of [[A, B], [B, A]]) {
    for (const P of X) {
      for (let i = 0; i + 1 < P.length; i++) {
        const m = [(P[i][0] + P[i + 1][0]) / 2, (P[i][1] + P[i + 1][1]) / 2];
        if (distance(m, Y) > tol) far += Math.hypot(P[i + 1][0] - P[i][0], P[i + 1][1] - P[i][1]);
      }
    }
  }
  return far;
}

const out = { views: [] };
for (const file of fs.readdirSync(dir).filter(f => f.endsWith('.json')).sort()) {
  const data = JSON.parse(fs.readFileSync(path.join(dir, file), 'utf8'));
  for (const view of data.views) {
    const M = turn.prepare(view), fig = view.figure, a0 = fig.camera.azimuth, e0 = fig.camera.elevation;
    const [x0, x1, y0, y1] = fig.box, size = Math.max(x1 - x0, y1 - y0);
    const at = turn.draw(M, a0, e0, false), pub = byClass(fig.layers), mine = byClass(at.layers);
    const lines = {}, fills = {};
    for (const c of new Set([...Object.keys(pub.lines), ...Object.keys(mine.lines)])) {
      lines[c] = [length(pub.lines[c] || []), length(mine.lines[c] || [])];
    }
    for (const c of new Set([...Object.keys(pub.fills), ...Object.keys(mine.fills)])) {
      fills[c] = [area(pub.fills[c] || []), area(mine.fills[c] || [])];
    }
    const labels = Math.max(0, ...at.labels.map((L, i) => Math.hypot(L.at[0] - fig.labels[i].at[0], L.at[1] - fig.labels[i].at[1])));
    const shown = at.labels.every(L => L.shown);
    const sweep = [];
    for (const da of [0, 90, 200]) {
      for (const e of [-90, -45, 0, e0, 45, 90]) {
        const d = turn.draw(M, a0 + da, e, true);
        let outside = 0, bad = 0;
        const check = p => {
          if (!Number.isFinite(p[0]) || !Number.isFinite(p[1])) { bad++; return; }
          outside = Math.max(outside, x0 - p[0], p[0] - x1, y0 - p[1], p[1] - y1);
        };
        for (const L of d.layers) (L.kind === 'point' ? [[L.at]] : [L.points, ...(L.holes || [])]).forEach(r => r.forEach(check));
        d.labels.forEach(L => check(L.at));
        sweep.push({ azimuth: a0 + da, elevation: e, scale: d.scale, outside: outside / size, bad });
      }
    }
    // A surface of revolution looks the same from every side: turned round its axis, its outline,
    // its rims and its circles, seen and hidden, stand where they stood, all but slivers where a
    // circle crosses the outline, which fall between its points differently. A curve marked on
    // it that is no circle turns with it, as its meridians do.
    const curves = new Set(view.surfaces.flatMap(s => (s.curves || []).map(c => c['class'])));
    const still = c => !/^(meridian|cut|reference)(-far)?$/.test(c) && !curves.has(c.replace(/-far$/, ''));
    let moved = 0, total = 0;
    for (const e of [e0, 60, -30]) {
      const A = byClass(turn.draw(M, a0, e, false).layers).lines, B = byClass(turn.draw(M, a0 + 37, e, false).layers).lines;
      for (const c of new Set([...Object.keys(A), ...Object.keys(B)].filter(still))) {
        moved += astray(A[c] || [], B[c] || [], 2e-3 * size);
        total += length(A[c] || []) + length(B[c] || []);
      }
    }
    moved /= total;
    out.views.push({ metric: data.metric, view: view.id, size, boxArea: (x1 - x0) * (y1 - y0),
                     lines, fills, labels, shown, sweep, moved });
  }
}

const flamm = JSON.parse(fs.readFileSync(path.join(dir, 'schwarzschild.json'), 'utf8')).views[0];
const M = turn.prepare(flamm);
out.schwarzschild = {};
for (const [name, e] of [['above', 90], ['below', -90]]) {
  // Looking down the axis every point of a circle stands its radius from the axis on the page.
  const d = turn.draw(M, -90, e, false), reach = {};
  for (const L of d.layers) {
    if (L.kind !== 'line' || !/^r2?(-far)?$/.test(L['class'])) continue;
    for (const p of L.points) {
      const r = Math.hypot(p[0], p[1] - M.surfaces[0].at[1]) / d.scale;
      const k = L['class'];
      reach[k] = [Math.min(r, (reach[k] || [Infinity])[0]), Math.max(r, (reach[k] || [0, -Infinity])[1])];
    }
  }
  out.schwarzschild[name] = reach;
}
process.stdout.write(JSON.stringify(out) + '\n');
