/* The graph of the relations between the spacetimes, drawn behind the page.

   _layouts/mfs.html draws every spacetime as a point, every two spacetimes that list each other
   under "Related Spacetimes" as a line between their points, gathered in the room to the right of
   the list where the spacetime's panel rests and drawn across the whole window, and gathers the spacetimes a search finds into a view of their own,
   the rest standing back behind them. This file is its geometry and touches no page: where each
   point stands, how large it is, which names have the room to be written, what a press lands on
   and how a gathering moves. The page runs it, a worker runs it so that a layout is worked out
   off the page's main thread, and the tests run it in Node. _tools/README.md, "The graph behind
   the list", says what the page does with it.

   It is the timeline's graph, Scripts/timeline-page.html and Scripts/build-timeline.py in the
   application's repository, brought here with the numbers it has there: the stress majorization
   in three dimensions, the size of a point growing linearly with its relations, the spring a
   gathering moves on, the forty names at most and where each name is written. */
(function(root) {
  'use strict';

  // The eye stands CAMERA radii of the collection from its middle, the collection fills FILL of
  // the shorter side of its room, and a point of drag turns the figure TURN radians.
  var CAMERA = 3.2, FILL = 0.9, TURN = 0.008;
  // A name stands GAP clear of its own point and CLEAR of every other name and point.
  var GAP = 3, CLEAR = 2;
  // A press lands on a point from this far off, since a point is a few pixels across.
  var REACH = 9;
  // A spacetime related to nothing is THIN of the plain radius and the most related one WIDE.
  var THIN = 0.75, WIDE = 2.25;
  // A point is DOT in radius before its relations widen it and its nearness scales it, and
  // nearness draws nothing under SMALLEST of its plain size or over LARGEST. Two points stand
  // APART clear of each other, rim to rim.
  var DOT = 3.5, SMALLEST = 0.85, LARGEST = 1.3, APART = 2;
  // A name stands on a line LINE tall, one longer than LONG letters on two lines, and before the
  // page has measured it a letter is LETTER wide, its width at SIZE pixels in Source Code Pro.
  // No more than MOST names are written at once, since the room holds only so many.
  var SIZE = 13, LINE = 16, LONG = 20, LETTER = 7.83, MOST = 40;

  // The graph the page draws, from the list's index, in its order, and the relations the build
  // writes, MFS/assets/data/relations.json: each spacetime's id, its name in the list and how
  // many spacetimes it is related to, and each relation as the places of its two spacetimes.
  function model(index, relations) {
    var place = {}, spacetimes = index.map(function(m, i) {
      place[m.id] = i;
      return { id: m.id, name: m.name, degree: 0 };
    });
    var edges = [];
    (relations.edges || []).forEach(function(e) {
      var a = place[e[0]], b = place[e[1]];
      if (a === undefined || b === undefined || a === b) return;
      edges.push([a, b]);
      spacetimes[a].degree++;
      spacetimes[b].degree++;
    });
    return { spacetimes: spacetimes, edges: edges };
  }

  // The places in the model of the spacetimes with the given ids, in the model's order; ids it
  // does not have are passed over.
  function placesOf(graph, ids) {
    var wanted = {};
    ids.forEach(function(id) { wanted[id] = true; });
    var found = [];
    graph.spacetimes.forEach(function(s, i) { if (wanted[s.id]) found.push(i); });
    return found;
  }

  function start() {
    return { turn: [1, 0, 0, 0, 1, 0, 0, 0, 1] };
  }
  function isStart(camera) {
    var s = start().turn;
    return camera.turn.every(function(v, k) { return v === s[k]; });
  }

  // The figure turned about the screen's upright by `across` and about its level by `down`, so
  // whatever faces the reader follows the hand. The turn is made a rotation again every time,
  // since a thousand small turns multiplied together drift off being one.
  function turned(camera, across, down) {
    var m = camera.turn;
    var cb = Math.cos(across), sb = Math.sin(across), ca = Math.cos(down), sa = Math.sin(down);
    var r0 = [0, 0, 0], r1 = [0, 0, 0], k;
    for (k = 0; k < 3; k++) {
      var x = cb * m[k] + sb * m[6 + k], z = -sb * m[k] + cb * m[6 + k], y = m[3 + k];
      r0[k] = x;
      r1[k] = ca * y - sa * z;
    }
    var n0 = Math.sqrt(r0[0] * r0[0] + r0[1] * r0[1] + r0[2] * r0[2]);
    for (k = 0; k < 3; k++) r0[k] /= n0;
    var along = r1[0] * r0[0] + r1[1] * r0[1] + r1[2] * r0[2];
    for (k = 0; k < 3; k++) r1[k] -= along * r0[k];
    var n1 = Math.sqrt(r1[0] * r1[0] + r1[1] * r1[1] + r1[2] * r1[2]);
    for (k = 0; k < 3; k++) r1[k] /= n1;
    var r2 = [r0[1] * r1[2] - r0[2] * r1[1], r0[2] * r1[0] - r0[0] * r1[2], r0[0] * r1[1] - r0[1] * r1[0]];
    return { turn: r0.concat(r1, r2) };
  }

  // How many times the plain radius a point of `degree` relations is drawn at, where the most
  // related spacetime of the collection has `most`. It grows linearly, so every relation more
  // widens a point by the same step: the least related is still seen, and the most related is
  // WIDE / THIN across the least and no more, whatever the collection.
  function girth(degree, most) {
    if (!(most > 0)) return 1;
    return THIN + (WIDE - THIN) * Math.max(0, Math.min(degree, most)) / most;
  }

  // How much larger than plain a point and its name are drawn where nearness scales them by
  // `scale`, and the radius of a point drawn so.
  function drawnAt(scale) {
    return Math.max(SMALLEST, Math.min(LARGEST, scale));
  }
  function radius(f, degree, most) {
    return DOT * f * girth(degree, most);
  }

  // Where a point is seen in a room W by H: its place, how near the reader it stands from 0 at
  // the back of the collection to 1 at the front, and how much larger nearness draws it.
  function seen(camera, p, W, H) {
    var m = camera.turn;
    var X = m[0] * p.x + m[1] * p.y + m[2] * p.z;
    var Y = m[3] * p.x + m[4] * p.y + m[5] * p.z;
    var Z = m[6] * p.x + m[7] * p.y + m[8] * p.z;
    var scale = CAMERA / (CAMERA - Z);
    var reach = FILL * Math.min(W, H) / 2 * scale;
    return { x: W / 2 + reach * X, y: H / 2 - reach * Y, near: Math.max(0, Math.min(1, (Z + 1) / 2)), scale: scale };
  }

  // A name's box on one side of its point: after it, before it, above, below, and then the four
  // corners, which a name takes only where no side is free.
  var SIDES = 8;
  function box(item, side) {
    var g = item.r + GAP, c = g * 0.7;
    if (side === 0) return [item.x + g, item.y - item.h / 2, item.w, item.h];
    if (side === 1) return [item.x - g - item.w, item.y - item.h / 2, item.w, item.h];
    if (side === 2) return [item.x - item.w / 2, item.y - g - item.h, item.w, item.h];
    if (side === 3) return [item.x - item.w / 2, item.y + g, item.w, item.h];
    return [side % 2 ? item.x - c - item.w : item.x + c, side < 6 ? item.y - c - item.h : item.y + c, item.w, item.h];
  }
  function dot(item) {
    return [item.x - item.r, item.y - item.r, 2 * item.r, 2 * item.r];
  }
  function meet(a, b) {
    return a[0] < b[0] + b[2] + CLEAR && b[0] < a[0] + a[2] + CLEAR &&
           a[1] < b[1] + b[3] + CLEAR && b[1] < a[1] + a[3] + CLEAR;
  }
  function within(b, edge) {
    return b[0] >= edge[0] && b[1] >= edge[1] && b[0] + b[2] <= edge[2] && b[1] + b[3] <= edge[3];
  }

  // A long name is set on two lines, broken at the hyphen or space nearest its middle.
  function twoLines(name) {
    if (name.length <= LONG) return [name];
    var best = -1;
    for (var k = 1; k < name.length - 1; k++) {
      if ((name[k] === '-' || name[k] === ' ') && (best < 0 || Math.abs(k - name.length / 2) < Math.abs(best - name.length / 2))) best = k;
    }
    if (best < 0) return [name];
    return name[best] === '-' ? [name.slice(0, best + 1), name.slice(best + 1)] : [name.slice(0, best), name.slice(best + 1)];
  }

  // The order names are given room in: the lower rank first, then a name already written, so a
  // name does not come and go as the figure turns, then the spacetime related to more of the
  // collection, then the nearer the reader.
  function inOrder(items) {
    return items.map(function(item, i) { return i; }).sort(function(a, b) {
      var A = items[a], B = items[b];
      return A.rank - B.rank || (B.side >= 0) - (A.side >= 0) || B.weight - A.weight || B.near - A.near || a - b;
    });
  }

  // The side each name is written on, or -1 for a name with no room. A name is written whole
  // inside the room W by H, or inside `edge`, [left, top, right, bottom] in the room's own
  // pixels, where the drawing reaches past the room, clear of every point and of every name given room before it, and keeps the
  // side it had (item.side) while that side is still free; no more than MOST are written. The
  // name that `must` be written, the one under the pointer, gives up the room's edge and the
  // points before it gives up being written, and still never lies on another name. A point
  // standing `behind`, sent back while a search gathers others, has no name unless it must and
  // stands in the way of none.
  function labels(items, W, H, edge) {
    var dots = items.map(dot), placed = [], sides = items.map(function() { return -1; }), written = 0;
    edge = edge || [0, 0, W, H];
    inOrder(items).forEach(function(i) {
      var item = items[i];
      if (item.x < edge[0] || item.y < edge[1] || item.x > edge[2] || item.y > edge[3]) return;
      if (item.behind && !item.must) return;
      if (written >= MOST && !item.must) return;
      var tries = item.side >= 0 ? [item.side] : [];
      for (var side = 0; side < SIDES; side++) tries.push(side);
      function free(strict) {
        for (var t = 0; t < tries.length; t++) {
          var b = box(item, tries[t]), ok = true, k;
          if (strict) {
            ok = within(b, edge);
            for (k = 0; ok && k < dots.length; k++) if (k !== i && !items[k].behind && meet(b, dots[k])) ok = false;
          }
          for (k = 0; ok && k < placed.length; k++) if (meet(b, placed[k])) ok = false;
          if (ok) { placed.push(b); return tries[t]; }
        }
        return -1;
      }
      sides[i] = free(true);
      if (sides[i] < 0 && item.must) sides[i] = free(false);
      if (sides[i] >= 0 && !item.must) written++;
    });
    return sides;
  }

  // The point a press at (px, py) lands on: the point it is on, and where it is on two the one
  // whose middle is nearer, failing that the name it lies in, and failing that the nearest point
  // within reach. -1 where it lands on none. Only a spacetime the list shows can be pointed at or
  // pressed, so a point a search sent behind is passed over, as if it were not there.
  function hit(items, sides, px, py) {
    var nearest = -1, least = Infinity, on = -1, onLeast = Infinity;
    items.forEach(function(item, i) {
      if (item.behind) return;
      var d = Math.sqrt((item.x - px) * (item.x - px) + (item.y - py) * (item.y - py));
      if (d < least) { nearest = i; least = d; }
      if (d < item.r + CLEAR && d < onLeast) { on = i; onLeast = d; }
    });
    if (on >= 0) return on;
    for (var i = 0; i < items.length; i++) {
      if (sides[i] < 0 || items[i].behind) continue;
      var b = box(items[i], sides[i]);
      if (px >= b[0] && px <= b[0] + b[2] && py >= b[1] && py <= b[1] + b[3]) return i;
    }
    return least <= REACH ? nearest : -1;
  }

  // Spacetimes placed by stress majorization over the relations between two of them and no
  // other: an edge is LENGTH long and STRETCH longer the fewer neighbours among them its two ends
  // share, so a family stands together and a spacetime related to half the collection does not
  // pull the whole of it into one knot, two spacetimes no path joins stand farther apart than any
  // two it does, the places start on a spiral over the sphere and are swept SWEEPS times over,
  // the middle is the origin, the widest spread lies along x and the narrowest along z, the way
  // the figure is first seen, and every place is drawn EVENNESS of the way to where it would stand
  // were the ball evenly filled, so names find room away from the middle.
  var LENGTH = 1, STRETCH = 2, SWEEPS = 400, EVENNESS = 0.75;

  // The places of `members`, spacetimes by their place in the model, in the unit ball, in their
  // order: the same for the same members and relations every time.
  function arranged(graph, members) {
    var n = members.length, place = graph.spacetimes.map(function() { return -1; }), i, j, k;
    members.forEach(function(m, k) { place[m] = k; });
    if (n === 0) return [];
    if (n === 1) return [[0, 0, 0]];
    var neighbours = members.map(function() { return []; }), pairs = [];
    graph.edges.forEach(function(e) {
      var a = place[e[0]], b = place[e[1]];
      if (a < 0 || b < 0) return;
      pairs.push([a, b]);
      neighbours[a].push(b);
      neighbours[b].push(a);
    });
    var d = [];
    for (i = 0; i < n; i++) {
      d.push([]);
      for (j = 0; j < n; j++) d[i].push(i === j ? 0 : Infinity);
    }
    pairs.forEach(function(p) {
      var a = p[0], b = p[1], own = {}, shared = 0, both = 0;
      neighbours[a].concat([a]).forEach(function(x) { own[x] = 1; both++; });
      neighbours[b].concat([b]).forEach(function(x) { if (own[x]) shared++; else both++; });
      d[a][b] = d[b][a] = LENGTH + STRETCH * (1 - shared / both);
    });
    for (k = 0; k < n; k++) {
      var dk = d[k];
      for (i = 0; i < n; i++) {
        var dik = d[i][k];
        if (dik === Infinity) continue;
        var di = d[i];
        for (j = 0; j < n; j++) if (dik + dk[j] < di[j]) di[j] = dik + dk[j];
      }
    }
    var farthest = 0;
    for (i = 0; i < n; i++) for (j = 0; j < n; j++) if (d[i][j] !== Infinity && d[i][j] > farthest) farthest = d[i][j];
    for (i = 0; i < n; i++) for (j = 0; j < n; j++) if (d[i][j] === Infinity) d[i][j] = farthest + LENGTH + STRETCH;
    var golden = Math.PI * (3 - Math.sqrt(5)), at = [];
    for (i = 0; i < n; i++) {
      var z = 1 - 2 * (i + 0.5) / n, r = Math.sqrt(1 - z * z);
      at.push([r * Math.cos(golden * i), r * Math.sin(golden * i), z]);
    }
    for (var sweep = 0; sweep < SWEEPS; sweep++) {
      for (i = 0; i < n; i++) {
        var pi = at[i], sx = 0, sy = 0, sz = 0, total = 0;
        for (j = 0; j < n; j++) {
          if (i === j) continue;
          var pj = at[j], want = d[i][j], weight = 1 / (want * want);
          var dx = pi[0] - pj[0], dy = pi[1] - pj[1], dz = pi[2] - pj[2];
          var pull = want / (Math.sqrt(dx * dx + dy * dy + dz * dz) || 1e-9);
          sx += weight * (pj[0] + pull * dx);
          sy += weight * (pj[1] + pull * dy);
          sz += weight * (pj[2] + pull * dz);
          total += weight;
        }
        pi[0] = sx / total; pi[1] = sy / total; pi[2] = sz / total;
      }
    }
    var middle = [0, 1, 2].map(function(c) { var s = 0; at.forEach(function(p) { s += p[c]; }); return s / n; });
    at = at.map(function(p) { return [p[0] - middle[0], p[1] - middle[1], p[2] - middle[2]]; });
    var frame = axes(at);
    at = at.map(function(p) { return frame.map(function(axis) { return p[0] * axis[0] + p[1] * axis[1] + p[2] * axis[2]; }); });
    var out = at.map(function(p) { return Math.sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]); });
    var reach = Math.max.apply(null, out) || 1, rank = [];
    members.map(function(m, k) { return k; }).sort(function(a, b) { return out[a] - out[b] || a - b; })
      .forEach(function(k, place) { rank[k] = place; });
    return at.map(function(p, k) {
      var even = Math.pow((rank[k] + 1) / n, 1 / 3);
      var scale = out[k] ? ((1 - EVENNESS) * out[k] / reach + EVENNESS * even) / out[k] : 0;
      return [p[0] * scale, p[1] * scale, p[2] * scale];
    });
  }

  // The three directions a cloud of points spreads along, widest first, as rows, by Jacobi's
  // rotations on its covariance, each signed so the point lying farthest along it lies on its
  // positive side, and never mirrored.
  function axes(points) {
    var n = points.length, a = [], v = [[1, 0, 0], [0, 1, 0], [0, 0, 1]], i, j, k;
    for (i = 0; i < 3; i++) {
      a.push([]);
      for (j = 0; j < 3; j++) { var s = 0; points.forEach(function(p) { s += p[i] * p[j]; }); a[i].push(s / n); }
    }
    var offs = [[0, 1], [0, 2], [1, 2]];
    for (var step = 0; step < 60; step++) {
      var best = offs[0];
      offs.forEach(function(pq) { if (Math.abs(a[pq[0]][pq[1]]) > Math.abs(a[best[0]][best[1]])) best = pq; });
      var p = best[0], q = best[1];
      if (Math.abs(a[p][q]) < 1e-15) break;
      var theta = 0.5 * Math.atan2(2 * a[p][q], a[q][q] - a[p][p]), c = Math.cos(theta), sn = Math.sin(theta), t;
      for (k = 0; k < 3; k++) { t = a[k][p]; a[k][p] = c * t - sn * a[k][q]; a[k][q] = sn * t + c * a[k][q]; }
      for (k = 0; k < 3; k++) { t = a[p][k]; a[p][k] = c * t - sn * a[q][k]; a[q][k] = sn * t + c * a[q][k]; }
      for (k = 0; k < 3; k++) { t = v[k][p]; v[k][p] = c * t - sn * v[k][q]; v[k][q] = sn * t + c * v[k][q]; }
    }
    var order = [0, 1, 2].sort(function(x, y) { return a[y][y] - a[x][x] || x - y; });
    var rows = order.slice(0, 2).map(function(i) { return [v[0][i], v[1][i], v[2][i]]; });
    rows.forEach(function(axis) {
      var far = 0;
      points.forEach(function(p) {
        var along = p[0] * axis[0] + p[1] * axis[1] + p[2] * axis[2];
        if (Math.abs(along) > Math.abs(far)) far = along;
      });
      if (far < 0) for (var c = 0; c < 3; c++) axis[c] = -axis[c];
    });
    var x = rows[0], y = rows[1];
    rows.push([x[1] * y[2] - x[2] * y[1], x[2] * y[0] - x[0] * y[2], x[0] * y[1] - x[1] * y[0]]);
    return rows;
  }

  // The members gathered into a room W by H, as the figure is first seen: filling the room's
  // middle FILL, their widest spread along its longer side, no point drawn on another, and as
  // many of their names as there is room for written beside them by labels(): every one, never
  // more than MOST, and fewer only where the room has none, a quarter fewer a time. A name has its
  // room after its point, the first side labels() tries, kept clear by ROOM past CLEAR, since a
  // name is measured a little wider than its letters, and every push goes PAST what it has to by
  // a little, so two things never come to rest on the line between touching and not. The place
  // of each member in the ball, in their order, and how many of them, first in the order names
  // are given room, have it. Only the members and the relations between them decide it, and the
  // same room gives the same places every time. The whole collection is gathered so too, which
  // is the graph as the page first shows it.
  var ROOM = 2, PAST = 0.1, ROUNDS = 1000, CROWD = 0.5;
  function gathered(graph, members, W, H, arrangement) {
    var base = arrangement || arranged(graph, members), most = 0;
    graph.spacetimes.forEach(function(s) { most = Math.max(most, s.degree); });
    if (!members.length) return { places: [], named: 0 };
    var reach = FILL * Math.min(W, H) / 2, margin = (1 - FILL) * Math.min(W, H) / 2, upright = H > W;
    var items = base.map(function(p, k) {
      var s = graph.spacetimes[members[k]], x = upright ? -p[1] : p[0], y = upright ? p[0] : p[1], z = p[2];
      var scale = CAMERA / (CAMERA - z), f = drawnAt(scale), parts = twoLines(s.name), letters = 0;
      parts.forEach(function(part) { letters = Math.max(letters, part.length); });
      return { u: scale * x, v: scale * y, z: z, scale: scale, x: W / 2, y: H / 2, near: Math.max(0, Math.min(1, (z + 1) / 2)),
               weight: s.degree, f: f, r: radius(f, s.degree, most), w: letters * LETTER * f, h: parts.length * LINE * f,
               rank: 2, must: false, side: -1 };
    });
    // A room whose middle the footprints would cover past CROWD of has never had the room for
    // them, so no more names are tried than leave it that free.
    var order = inOrder(items), frame = (W - 2 * margin) * (H - 2 * margin), covered = 0, count = 0;
    items.forEach(function(item) { covered += (2 * item.r + APART) * (2 * item.r + APART); });
    while (count < Math.min(items.length, MOST)) {
      var item = items[order[count]], more = (item.w + GAP + CLEAR) * (item.h + CLEAR);
      if ((covered + more) / frame > CROWD) break;
      covered += more;
      count++;
    }
    var counts = [count];
    while (counts[counts.length - 1] > 0) counts.push(Math.floor(counts[counts.length - 1] * 3 / 4));
    var tried, named;
    for (var c = 0; c < counts.length; c++) {
      named = items.map(function() { return false; });
      order.slice(0, counts[c]).forEach(function(k) { named[k] = true; });
      tried = roomFor(items, named, W, H, margin);
      if (tried.clear && writes(tried.items, order.slice(0, counts[c]), W, H)) break;
    }
    return {
      places: tried.items.map(function(item) {
        return [(item.x - W / 2) / (reach * item.scale), (H / 2 - item.y) / (reach * item.scale), item.z];
      }),
      named: counts[Math.min(c, counts.length - 1)]
    };
  }
  // Whether labels() writes the name of every one of `owed`.
  function writes(items, owed, W, H) {
    var sides = labels(items, W, H);
    return owed.every(function(k) { return sides[k] >= 0; });
  }
  // The members spread to fill the room's middle with the footprint of each, its point and,
  // where `named`, the name after it, inside, and then moved apart, a pair at a time and each
  // half the way, until no two points are nearer than APART rim to rim and no name lies on
  // another or on a point but its own. `clear` says whether that was reached.
  function roomFor(start, named, W, H, margin) {
    var n = start.length, i, j;
    var items = start.map(function(item) { var copy = {}; for (var key in item) copy[key] = item[key]; return copy; });
    var left = items.map(function(item) { return item.r; });
    var right = items.map(function(item, k) { return named[k] ? item.r + GAP + item.w : item.r; });
    var half = items.map(function(item, k) { return named[k] ? Math.max(item.r, item.h / 2) : item.r; });
    var across = W - 2 * margin, tall = H - 2 * margin, fill = Infinity;
    for (i = 0; i < n; i++) for (j = 0; j < n; j++) {
      if (items[i].u > items[j].u) fill = Math.min(fill, (across - right[i] - left[j]) / (items[i].u - items[j].u));
      if (items[i].v > items[j].v) fill = Math.min(fill, (tall - half[i] - half[j]) / (items[i].v - items[j].v));
    }
    fill = isFinite(fill) ? Math.max(0, fill) : 0;
    var x0 = -Infinity, x1 = Infinity, y0 = -Infinity, y1 = Infinity;
    items.forEach(function(item, k) {
      x0 = Math.max(x0, margin + left[k] - fill * item.u); x1 = Math.min(x1, W - margin - right[k] - fill * item.u);
      y0 = Math.max(y0, margin + half[k] + fill * item.v); y1 = Math.min(y1, H - margin - half[k] + fill * item.v);
    });
    var cx = (x0 + x1) / 2, cy = (y0 + y1) / 2;
    items.forEach(function(item) { item.x = cx + fill * item.u; item.y = cy - fill * item.v; });
    return settled(items, named, W, H, margin);
  }
  // The members moved apart where they stand, as roomFor() moves them.
  function settled(items, named, W, H, margin) {
    var n = items.length, i, j, k, apart = CLEAR + ROOM;
    var left = items.map(function(item) { return item.r; });
    var right = items.map(function(item, k) { return named[k] ? item.r + GAP + item.w : item.r; });
    var half = items.map(function(item, k) { return named[k] ? Math.max(item.r, item.h / 2) : item.r; });
    var extent = items.map(function(item, k) { return Math.max(left[k], right[k], half[k]) + apart + APART; });
    function shapes(k) {
      var item = items[k], own = [item.x - item.r, item.y - item.r, 2 * item.r, 2 * item.r];
      return named[k] ? [own, box(item, 0)] : [own];
    }
    function held(k) {
      var item = items[k];
      item.x = Math.max(margin + left[k], Math.min(W - margin - right[k], item.x));
      item.y = Math.max(margin + half[k], Math.min(H - margin - half[k], item.y));
    }
    var golden = Math.PI * (3 - Math.sqrt(5)), clear = false, pad = apart + APART + ROOM + PAST;
    var lined = items.map(function(item, k) { return k; });
    function from(k) { return items[k].x - left[k] - pad; }
    // Each round goes along the room from its leading edge, and two whose footprints, with room
    // round them, do not meet across it are no pair.
    for (var round = 0; round < ROUNDS && !clear; round++) {
      clear = true;
      lined.sort(function(a, b) { return from(a) - from(b) || a - b; });
      for (var one = 0; one < n; one++) for (var other = one + 1; other < n; other++) {
        if (from(lined[other]) > items[lined[one]].x + right[lined[one]] + pad) break;
        i = Math.min(lined[one], lined[other]); j = Math.max(lined[one], lined[other]);
        var a = items[i], b = items[j], dx = b.x - a.x, dy = b.y - a.y, reachOf = extent[i] + extent[j];
        if (Math.abs(dy) > reachOf) continue;
        var d = Math.sqrt(dx * dx + dy * dy), need = a.r + b.r + APART + ROOM;
        if (d < need) {
          var ux = d > 1e-9 ? dx / d : Math.cos(golden * (i * n + j)), uy = d > 1e-9 ? dy / d : Math.sin(golden * (i * n + j));
          var push = (need - d + PAST) / 2;
          a.x -= ux * push; a.y -= uy * push; b.x += ux * push; b.y += uy * push;
          clear = false;
        }
        var A = shapes(i), B = shapes(j);
        for (var s = 0; s < A.length; s++) for (var t = 0; t < B.length; t++) {
          if (s === 0 && t === 0) continue;
          var p = A[s], q = B[t];
          var px = Math.min(p[0] + p[2] + apart - q[0], q[0] + q[2] + apart - p[0]);
          var py = Math.min(p[1] + p[3] + apart - q[1], q[1] + q[3] + apart - p[1]);
          if (px <= 0 || py <= 0) continue;
          if (px <= py) {
            var sx = p[0] + p[2] / 2 <= q[0] + q[2] / 2 ? 1 : -1;
            a.x -= sx * (px + PAST) / 2; b.x += sx * (px + PAST) / 2;
          } else {
            var sy = p[1] + p[3] / 2 <= q[1] + q[3] / 2 ? 1 : -1;
            a.y -= sy * (py + PAST) / 2; b.y += sy * (py + PAST) / 2;
          }
          clear = false;
          A = shapes(i); B = shapes(j);
        }
      }
      for (k = 0; k < n; k++) held(k);
    }
    return { items: items, clear: clear && apartAll(items) };
  }
  // Whether every two points stand APART clear, rim to rim.
  function apartAll(items) {
    for (var i = 0; i < items.length; i++) for (var j = i + 1; j < items.length; j++) {
      if (Math.hypot(items[i].x - items[j].x, items[i].y - items[j].y) < items[i].r + items[j].r + APART) return false;
    }
    return true;
  }

  // Where a spacetime stands while a search gathers others: RECEDE of the way to the back of the
  // ball and drawn that much smaller, so the whole collection still stands behind what the
  // search found, for the reader to see what it left out.
  var RECEDE = 0.5;
  function receded(p) {
    return { x: p.x * (1 - RECEDE), y: p.y * (1 - RECEDE), z: p.z * (1 - RECEDE) - RECEDE };
  }

  // Where every spacetime is going and whether it stands behind. `home` is the whole collection
  // gathered, each spacetime's place in the ball; `found` is the places of the spacetimes a
  // search found, null while nothing is searched, and `gathering` their places gathered, in
  // their order. With nothing searched every spacetime goes home and none stands behind; with a
  // search, what it found goes to its gathered place and everything else stands back behind.
  function goal(home, found, gathering) {
    if (!found) return { to: home.map(copy), behind: home.map(function() { return false; }) };
    var to = home.map(receded), behind = home.map(function() { return true; });
    found.forEach(function(i, k) {
      var p = gathering[k];
      to[i] = { x: p[0], y: p[1], z: p[2] };
      behind[i] = false;
    });
    return { to: to, behind: behind };
  }
  function copy(p) {
    return { x: p.x, y: p.y, z: p.z };
  }

  // A gathering moves on a spring: SPRING is its angular frequency in radians a second and
  // DAMPING its share of critical damping, so it passes its mark by under two parts in a hundred
  // and has come within two parts in a thousand of it by SETTLE seconds, where it stands still.
  // How far it has gone, from 0 to 1, `seconds` after it set off.
  var SPRING = 7, DAMPING = 0.8, SETTLE = 1.2;
  function sprung(seconds) {
    if (!(seconds > 0)) return 0;
    if (seconds >= SETTLE) return 1;
    var ring = SPRING * Math.sqrt(1 - DAMPING * DAMPING), decay = DAMPING * SPRING;
    return 1 - Math.exp(-decay * seconds) * (Math.cos(ring * seconds) + decay / ring * Math.sin(ring * seconds));
  }

  // A turn part of the way, `p`, from the turn `a` to the turn `b`, the shortest way round.
  function quaternion(m) {
    var trace = m[0] + m[4] + m[8], s;
    if (trace > 0) { s = 2 * Math.sqrt(trace + 1); return [s / 4, (m[7] - m[5]) / s, (m[2] - m[6]) / s, (m[3] - m[1]) / s]; }
    if (m[0] > m[4] && m[0] > m[8]) { s = 2 * Math.sqrt(1 + m[0] - m[4] - m[8]); return [(m[7] - m[5]) / s, s / 4, (m[1] + m[3]) / s, (m[2] + m[6]) / s]; }
    if (m[4] > m[8]) { s = 2 * Math.sqrt(1 + m[4] - m[0] - m[8]); return [(m[2] - m[6]) / s, (m[1] + m[3]) / s, s / 4, (m[5] + m[7]) / s]; }
    s = 2 * Math.sqrt(1 + m[8] - m[0] - m[4]);
    return [(m[3] - m[1]) / s, (m[2] + m[6]) / s, (m[5] + m[7]) / s, s / 4];
  }
  function rotation(q) {
    var l = Math.sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
    var w = q[0] / l, x = q[1] / l, y = q[2] / l, z = q[3] / l;
    return [1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y),
            2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x),
            2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)];
  }
  function blend(a, b, p) {
    if (p === 0) return a;
    if (p === 1) return b;
    var qa = quaternion(a.turn), qb = quaternion(b.turn), cos = 0, k;
    for (k = 0; k < 4; k++) cos += qa[k] * qb[k];
    if (cos < 0) { qb = qb.map(function(c) { return -c; }); cos = -cos; }
    var angle = Math.acos(Math.min(1, cos)), wa = 1 - p, wb = p;
    if (angle > 1e-6) { wa = Math.sin((1 - p) * angle) / Math.sin(angle); wb = Math.sin(p * angle) / Math.sin(angle); }
    return { turn: rotation([0, 1, 2, 3].map(function(k) { return wa * qa[k] + wb * qb[k]; })) };
  }

  /* The layouts the page asks for, worked out where they are asked: in a worker, so that the
     page's main thread, which the list and the spacetime need, is never held by them. A layout is
     kept once worked out, for the members and the room it was worked out for. The page sends the
     graph once, as `graph`, and then a request at a time, `members` and the room W by H, and is
     answered with the request's `seq`, the members' places gathered and how many are named. */
  function layouts() {
    var graph = null, arrangements = {}, gatherings = {};
    return function(message) {
      if (message.graph) { graph = message.graph; arrangements = {}; gatherings = {}; return null; }
      var key = message.members.join(','), sized = key + ' ' + message.W + ' ' + message.H;
      if (!arrangements[key]) arrangements[key] = arranged(graph, message.members);
      if (!gatherings[sized]) gatherings[sized] = gathered(graph, message.members, message.W, message.H, arrangements[key]);
      return { seq: message.seq, places: gatherings[sized].places, named: gatherings[sized].named };
    };
  }

  var MfsGraph = {
    CAMERA: CAMERA, FILL: FILL, TURN: TURN, GAP: GAP, CLEAR: CLEAR, REACH: REACH, THIN: THIN, WIDE: WIDE,
    DOT: DOT, SMALLEST: SMALLEST, LARGEST: LARGEST, APART: APART, SIZE: SIZE, LINE: LINE, LONG: LONG,
    LETTER: LETTER, MOST: MOST, SIDES: SIDES, LENGTH: LENGTH, STRETCH: STRETCH, SWEEPS: SWEEPS,
    EVENNESS: EVENNESS, RECEDE: RECEDE, SPRING: SPRING, DAMPING: DAMPING, SETTLE: SETTLE,
    model: model, placesOf: placesOf, start: start, isStart: isStart, turned: turned, girth: girth,
    drawnAt: drawnAt, radius: radius, seen: seen, box: box, dot: dot, meet: meet, within: within,
    twoLines: twoLines, inOrder: inOrder, labels: labels, hit: hit, arranged: arranged, axes: axes,
    gathered: gathered, apartAll: apartAll, receded: receded, goal: goal, sprung: sprung, blend: blend,
    layouts: layouts
  };

  if (typeof module !== 'undefined' && module.exports) module.exports = MfsGraph;
  else root.MfsGraph = MfsGraph;
  // Run as a worker, it answers the page's requests for layouts.
  if (typeof WorkerGlobalScope !== 'undefined' && root instanceof WorkerGlobalScope) {
    var answer = layouts();
    root.onmessage = function(event) {
      var reply = answer(event.data);
      if (reply) root.postMessage(reply);
    };
  }
})(typeof self !== 'undefined' ? self : this);
