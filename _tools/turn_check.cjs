#!/usr/bin/env node
/* Measure MFS/assets/turn.js, which the page turns a figure drawn in three dimensions with,
   against every published embedding figure and every published figure of light cones, for the
   tests in test_build_mfs_data.py to hold.

     node _tools/turn_check.cjs

   prints one JSON object. For every embedding view it gives, at the figure's own camera, the
   length of every line class and the area of every fill class as the page would draw them beside
   the published ones, and how far each label lands from its published place; at a sweep of other
   cameras, the scale, how far any drawn point or label falls outside the box, and any number
   that is not one; and what share of the lines that a turn round the axis must leave in place
   moves under one, for every view whose surfaces are surfaces of revolution. For Schwarzschild it
   also gives how far out on each sheet a circle is seen and how far out one is hidden, from
   straight above and from straight below. For a view that draws a height over a plane it gives,
   from straight above and straight below, how much of its lines is hidden and how much of the page
   its tint covers against its rim at the drawn scale; seen from just above the plane, how much of
   its grid is hidden; and, turned all the way round, how far it strays from its box and the
   least scale it is drawn at.

   For every figure of light cones it gives, at the figure's own camera, whether the layers come
   in the published order with the published kinds and classes, how far the drawn layers stray
   from the published ones, the cones' separately from the lines the generator thins, how far
   each label and slice lands from its published place and whether every label is shown; and at
   a sweep of other cameras the same as for an embedding view.

   For one figure of each kind, a surface of revolution, a height over a plane and a figure of
   light cones, it drags the drawing as the page does, with dragged(), keyed() and turned(): a
   quarter of its width across, a whole width up and down, a whole turn round and back home,
   giving each camera reached, whether the page counts it turned, how much of the drawing moved,
   how far the drawing strays from its box there, and for the light cones how far the cone
   facing the reader at the start moved across the page. */
'use strict';
const fs = require('fs');
const path = require('path');

const root = path.resolve(__dirname, '..');
const turn = require(path.join(root, 'MFS', 'assets', 'turn.js'));
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

/* A height over a plane looked at along the vertical hides nothing of itself, and its tint then
   covers its rim, drawn at the scale the turn chose; looked at from just above the plane, its
   relief hides part of its grid behind it; and turned all the way round it keeps to its box,
   drawn smaller where it would be wider or taller than it was. */
function heights(M, view, a0, e0) {
  const rim = view.surfaces[0].pieces.filter(p => p.grid).map(p => {
    const g = p.grid, polar = g.frame === 'polar', P = [];
    const node = (i, j) => polar ? [g.u[i] * Math.cos(g.v[j]), g.u[i] * Math.sin(g.v[j])] : [g.u[i], g.v[j]];
    if (polar) for (let j = 0; j < g.v.length; j++) P.push(node(g.u.length - 1, j));
    else P.push(node(0, 0), node(g.u.length - 1, 0), node(g.u.length - 1, g.v.length - 1), node(0, g.v.length - 1));
    return area([P]);
  }).reduce((a, b) => a + b, 0);
  const out = {};
  for (const [name, e] of [['above', 90], ['below', -90]]) {
    const d = turn.draw(M, a0, e, false), mine = byClass(d.layers);
    let far = 0, seen = 0;
    for (const [c, L] of Object.entries(mine.lines)) (/-far$/.test(c) ? (far += length(L)) : (seen += length(L)));
    out[name] = { far, seen, tint: area(mine.fills.cover || []), rim: rim * d.scale * d.scale };
  }
  const low = byClass(turn.draw(M, a0, 5, false).layers).lines;
  out.low = { far: length(low['grid-far'] || []), seen: length(low.grid || []) };
  const [x0, x1, y0, y1] = view.figure.box, size = Math.max(x1 - x0, y1 - y0);
  let outside = 0, least = 1;
  for (let a = a0; a < a0 + 360; a += 5) {
    const d = turn.draw(M, a, e0, true);
    least = Math.min(least, d.scale);
    for (const L of d.layers) {
      for (const p of (L.kind === 'point' ? [L.at] : L.points)) {
        outside = Math.max(outside, x0 - p[0], p[0] - x1, y0 - p[1], p[1] - y1);
      }
    }
  }
  out.turned = { outside: outside / size, scale: least };
  return out;
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
    // A height over a plane, and a stack of ellipses, which need not look the same from every side.
    const grids = view.surfaces.some(s => s.pieces.some(p => p.grid && p.grid.frame !== 'ellipses'));
    const stacks = [...view.surfaces, ...(view.movie ? view.movie.frames : [])]
      .some(s => s.pieces.some(p => p.grid && p.grid.frame === 'ellipses'));
    let moved = null;
    if (!grids && !stacks) {
      let total = 0;
      moved = 0;
      for (const e of [e0, 60, -30]) {
        const A = byClass(turn.draw(M, a0, e, false).layers).lines, B = byClass(turn.draw(M, a0 + 37, e, false).layers).lines;
        for (const c of new Set([...Object.keys(A), ...Object.keys(B)].filter(still))) {
          moved += astray(A[c] || [], B[c] || [], 2e-3 * size);
          total += length(A[c] || []) + length(B[c] || []);
        }
      }
      moved /= total;
    }
    const entry = { metric: data.metric, view: view.id, size, boxArea: (x1 - x0) * (y1 - y0),
                    lines, fills, labels, shown, sweep, moved };
    /* Each shade, drawn at the figure's own camera, covers what its published layers cover, class
       by class, and turned to the other side it still paints its class and keeps to the box; the
       figure without a shade is the published one again once the shade is taken away. */
    if (view.shades) {
      entry.shades = view.shades.map(sh => {
        turn.shade(M, sh);
        const drawn = byClass(turn.draw(M, a0, e0, false).layers).fills, published = byClass(sh.layers).fills, fills = {};
        for (const c of new Set([...Object.keys(published), ...Object.keys(drawn)])) {
          fills[c] = [area(published[c] || []), area(drawn[c] || [])];
        }
        const away = turn.draw(M, a0 + 180, -e0, false);
        return { view: sh.view, fills, away: area(byClass(away.layers).fills.shade || []), ...strays(away, fig.box) };
      });
      turn.shade(M, null);
      const back = byClass(turn.draw(M, a0, e0, false).layers).fills;
      entry.unshaded = Object.fromEntries(Object.keys(back).map(c => [c, area(back[c])]));
    }
    if (grids) entry.height = heights(M, view, a0, e0);
    if (stacks) entry.stack = true;
    /* A movie draws every frame at one scale about one place: at the figure's own camera every
       frame keeps to the box at scale 1, and turned, every frame keeps to it at the one scale the
       camera gives all of them. The label that names the frame shown stands at one place at each
       camera, whatever the frame: `label` is the farthest any frame puts it from the first's. */
    if (view.movie) {
      const frames = view.movie.frames.length, cams = [[a0, e0], [a0 + 90, 45], [a0 + 200, -45], [a0 + 30, 90], [a0, -90]];
      const named = view.figure.labels.map((L, i) => L.frame ? i : -1).filter(i => i >= 0);
      let outside = 0, bad = 0, scales = new Set(), home = 1, label = 0;
      for (const [a, e] of cams) {
        const seen = new Set();
        let first = null;
        for (let k = 0; k < frames; k++) {
          turn.frame(M, k);
          const d = turn.draw(M, a, e, true);
          seen.add(d.scale.toFixed(12));
          if (a === a0 && e === e0) home = Math.min(home, d.scale);
          const check = p => {
            if (!Number.isFinite(p[0]) || !Number.isFinite(p[1])) { bad++; return; }
            outside = Math.max(outside, x0 - p[0], p[0] - x1, y0 - p[1], p[1] - y1);
          };
          for (const L of d.layers) (L.kind === 'point' ? [[L.at]] : [L.points, ...(L.holes || [])]).forEach(r => r.forEach(check));
          d.labels.forEach(L => check(L.at));
          const at = named.map(i => d.labels[i].at);
          first = first || at;
          at.forEach((p, i) => { label = Math.max(label, Math.hypot(p[0] - first[i][0], p[1] - first[i][1])); });
        }
        scales.add(seen.size);
      }
      turn.frame(M, 0);
      entry.movie = { frames, outside: outside / size, bad, home, scalesPerCamera: Math.max(...scales), named: named.length, label };
    }
    out.views.push(entry);
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

// The largest distance from a point of either polyline to the other.
function apart(A, B) {
  let far = 0;
  for (const p of A) far = Math.max(far, distance(p, [B]));
  for (const p of B) far = Math.max(far, distance(p, [A]));
  return far;
}
// How far any drawn point, label or slice falls outside the box, as a share of its size, and how
// many numbers are not numbers.
function strays(d, box) {
  const [x0, x1, y0, y1] = box, size = Math.max(x1 - x0, y1 - y0);
  let outside = 0, bad = 0;
  const check = p => {
    if (!Number.isFinite(p[0]) || !Number.isFinite(p[1])) { bad++; return; }
    outside = Math.max(outside, x0 - p[0], p[0] - x1, y0 - p[1], p[1] - y1);
  };
  for (const L of d.layers) (L.kind === 'point' ? [[L.at]] : [L.points, ...(L.holes || [])]).forEach(r => r.forEach(check));
  for (const mark of d.slices) [...mark.lines, ...mark.fills.flat()].forEach(r => r.forEach(check));
  d.labels.forEach(L => check(L.at));
  return { outside: outside / size, bad };
}

const diagrams = path.join(root, 'MFS', 'assets', 'data', 'diagrams');
out.figures = [];
const figures = {};
for (const file of fs.readdirSync(diagrams).filter(f => f.endsWith('.json')).sort()) {
  const data = JSON.parse(fs.readFileSync(path.join(diagrams, file), 'utf8'));
  for (const [system, views] of Object.entries(data.projections || {})) {
    for (const fig of views) {
      if (!fig.turn) continue;
      figures[`${data.metric}/${fig.id}`] = fig;
      const M = turn.prepare(fig), a0 = fig.camera.azimuth, e0 = fig.camera.elevation;
      const at = turn.draw(M, a0, e0, false);
      const order = fig.layers.length === at.layers.length &&
        fig.layers.every((L, i) => L.kind === at.layers[i].kind && L['class'] === at.layers[i]['class']);
      // The generator thins only the lines it draws before the cones.
      const thinned = new Set(fig.turn.lines.map(L => L['class']));
      let cones = 0, lines = 0;
      if (order) {
        fig.layers.forEach((L, i) => {
          const D = at.layers[i];
          const far = L.kind === 'point' ? Math.hypot(L.at[0] - D.at[0], L.at[1] - D.at[1]) : apart(L.points, D.points);
          if (i < fig.turn.lines.length && thinned.has(L['class'])) lines = Math.max(lines, far);
          else cones = Math.max(cones, far);
        });
      }
      const labels = Math.max(0, ...at.labels.map((L, i) => Math.hypot(L.at[0] - fig.labels[i].at[0], L.at[1] - fig.labels[i].at[1])));
      let slices = 0;
      at.slices.forEach((mark, k) => {
        mark.lines.forEach((L, j) => { slices = Math.max(slices, apart(L, fig.slices[k].lines[j])); });
        mark.fills.forEach((rings, j) => rings.forEach((R, r) => { slices = Math.max(slices, apart(R, fig.slices[k].fills[j][r])); }));
      });
      const sweep = [];
      for (const da of [0, 90, 200]) {
        for (const e of [-90, -45, 0, e0, 45, 90]) {
          const d = turn.draw(M, a0 + da, e, true);
          sweep.push({ azimuth: a0 + da, elevation: e, scale: d.scale, ...strays(d, fig.box) });
        }
      }
      const [x0, x1, y0, y1] = fig.box;
      out.figures.push({ metric: data.metric, system, view: fig.id, size: Math.max(x1 - x0, y1 - y0), order,
                         cones, lines, labels, slices, shown: at.labels.every(L => L.shown), sweep });
    }
  }
}

/* The hand, on one figure of each kind, the drawing `width` pixels wide as a page shows it. A
   fifth of the width across turns the figure 36 degrees, which carries no meridian of Flamm's
   paraboloid onto another, and a whole width up or down would tilt it past straight down or
   straight up the axis, where it stops. The cone that faces the reader at the start, on the
   circle r_c at phi = -90 degrees, is followed across the page. */
out.hand = {};
const embedding = name => JSON.parse(fs.readFileSync(path.join(dir, name + '.json'), 'utf8')).views[0];
for (const [kind, view] of [['surface', embedding('schwarzschild')], ['height', embedding('krasnikov')],
                            ['cones', figures['godel/tipping']]]) {
  const fig = view.figure || view, M = turn.prepare(view), width = 600;
  const start = { azimuth: fig.camera.azimuth, elevation: fig.camera.elevation };
  const home = turn.draw(M, start.azimuth, start.elevation, false), before = byClass(home.layers).lines;
  const size = Math.max(fig.box[1] - fig.box[0], fig.box[3] - fig.box[2]);
  const steps = {
    across: turn.dragged(start, width / 5, 0, width),
    up: turn.dragged(start, 0, -width, width),
    down: turn.dragged(start, 0, width, width),
    round: turn.dragged(start, 2 * width, 0, width),
  };
  steps.home = turn.keyed(steps.across, 'Home', start);
  steps.left = turn.keyed(start, 'ArrowLeft', start);
  const entry = {};
  for (const [name, cam] of Object.entries(steps)) {
    const d = turn.draw(M, cam.azimuth, cam.elevation, false), after = byClass(d.layers).lines;
    let moved = 0, total = 0;
    for (const c of new Set([...Object.keys(before), ...Object.keys(after)])) {
      moved += astray(before[c] || [], after[c] || [], 2e-3 * size);
      total += length(before[c] || []) + length(after[c] || []);
    }
    entry[name] = { azimuth: cam.azimuth, elevation: cam.elevation, turned: turn.turned(start, cam),
                    moved: moved / total, ...strays(d, fig.box) };
    if (kind === 'cones') {
      const facing = fig.turn.cones.findIndex(c => Math.abs(c.apex[0]) < 1e-6 && c.apex[1] < -0.5);
      const apex = d => d.layers.find(L => L['class'] === 'cone-apex' && L.cone === facing).at[0];
      entry[name].facing = apex(d) - apex(home);
    }
  }
  out.hand[kind] = entry;
}
process.stdout.write(JSON.stringify(out) + '\n');
