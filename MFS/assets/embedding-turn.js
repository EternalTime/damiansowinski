/* Turning an embedding diagram.

   _layouts/mfs.html draws an embedding diagram's published figure, which
   _tools/derivations/embedding.py projected once from a fixed camera, until the reader drags it,
   and from then on draws the figure again from the view's surfaces at the camera the drag has
   reached, with draw() below. Everything here is geometry and touches no page, so the tests run
   it in Node against the published figures; _tools/README.md, "Turning the figure", is the
   definition the application follows as well.

   The figure is drawn by the rules the generator draws it by: each piece's meridians, its
   outline where it turns edge on to the camera, the rim at an end where the drawing stops, its
   marked circles and marked meridians, every line split where a surface hides it into the part
   seen and the part hidden, which carries -far, and each tinted piece filled where it is the
   surface nearest the camera. What hides a line is found as the generator finds it, by casting a
   ray from each of its points toward the camera through the truncated cones between neighbouring
   circles of every profile, which is exact for the surface the points describe, and a point of
   the outline, where the line of sight only grazes the surface, is judged a little off it on
   either side. The tint is found on a grid of the page, as the generator's is, from the same
   cones cut into facets and cast through exactly wherever two regions meet.

   Each surface turns about the point of its axis halfway up its drawing, which stays where it
   stood on the page. Its width on the page never changes as it turns, a surface of revolution
   being as wide as its widest circle from every side, but its height does as it tilts, so the
   surfaces are drawn smaller, all at one scale, only as far as a tilt needs more height than a
   surface had at the start, and moved up or down only as far as that height needs. So the
   figure keeps its box and its labels their room, and at the start it is the published figure. */
(function(root) {
  'use strict';

  var RAD = Math.PI / 180;

  /* The generator's camera, projections.Camera: an orthographic view from `azimuth`, measured
     round the axis from phi = 0 toward phi = pi/2, and `elevation` above the plane z = 0, both
     in degrees. A point P is drawn at (P . right, P . up) on the page, y up, and a larger
     P . toward is nearer the reader. */
  function camera(azimuth, elevation) {
    var a = azimuth * RAD, e = elevation * RAD;
    var ca = Math.cos(a), sa = Math.sin(a), ce = Math.cos(e), se = Math.sin(e);
    return { azimuth: azimuth, elevation: elevation,
             toward: [ce * ca, ce * sa, se], right: [-sa, ca, 0], up: [-se * ca, -se * sa, ce] };
  }

  // What draw() needs of a view that does not depend on the camera.
  function prepare(view) {
    var figure = view.figure, turn = figure.turn, box = figure.box;
    var size = Math.max(box[1] - box[0], box[3] - box[2]);
    var e0 = figure.camera.elevation * RAD, classes = [];
    Object.keys(turn.tint).forEach(function(k) { if (classes.indexOf(turn.tint[k]) < 0) classes.push(turn.tint[k]); });
    classes.sort();
    var surfaces = view.surfaces.map(function(s, k) {
      var low = Infinity, high = -Infinity, byId = {};
      var pieces = s.pieces.map(function(p) {
        var n = p.points.length, rho = new Float64Array(n), z = new Float64Array(n);
        for (var i = 0; i < n; i++) {
          rho[i] = p.points[i][1];
          z[i] = p.points[i][2];
          if (z[i] < low) low = z[i];
          if (z[i] > high) high = z[i];
        }
        var piece = { id: p.id, cls: p['class'], reference: !!p.reference, rho: rho, z: z,
                      ends: [p.start.kind, p.end.kind] };
        byId[p.id] = piece;
        return piece;
      });
      // The scene: every cone between two neighbouring circles of a piece of the slice, with
      // what its piece is tinted, 1 for a piece left clear and 2 on for the fill classes in
      // order. A reference piece is not part of the slice and hides nothing.
      var r1 = [], dr = [], z1 = [], dz = [], code = [];
      pieces.forEach(function(p) {
        if (p.reference) return;
        p.code = turn.tint[p.cls] ? 2 + classes.indexOf(turn.tint[p.cls]) : 1;
        for (var i = 0; i + 1 < p.rho.length; i++) {
          r1.push(p.rho[i]); dr.push(p.rho[i + 1] - p.rho[i]);
          z1.push(p.z[i]); dz.push(p.z[i + 1] - p.z[i]);
          code.push(p.code);
        }
      });
      var S = { pieces: pieces, byId: byId, rings: s.rings, zc: (low + high) / 2,
                r1: new Float64Array(r1), dr: new Float64Array(dr), z1: new Float64Array(z1), dz: new Float64Array(dz),
                code: new Uint8Array(code) };
      var o = turn.origins[k];
      S.at = [o[0], o[1] + S.zc * Math.cos(e0)];
      var h = height(S, e0);
      S.top0 = h[0];
      S.bottom0 = h[1];
      return S;
    });
    var order = {};
    figure.legend.forEach(function(item, i) { order[item[1]] = i; });
    return { figure: figure, turn: turn, box: box, size: size, surfaces: surfaces, order: order, classes: classes,
             eps: 1e-7 * size, unit: (box[1] - box[0]) / 560 };
  }

  /* How far a surface reaches above and below its centre on the page at elevation e: a point
     of its profile draws a circle whose height on the page runs over (z - zc) cos e +- rho sin e,
     and along a cone between two points both are linear, so the points bound the whole. */
  function height(S, e) {
    var ce = Math.cos(e), se = Math.abs(Math.sin(e)), top = -Infinity, bottom = Infinity;
    S.pieces.forEach(function(p) {
      for (var i = 0; i < p.rho.length; i++) {
        var y = (p.z[i] - S.zc) * ce, w = p.rho[i] * se;
        if (y + w > top) top = y + w;
        if (y - w < bottom) bottom = y - w;
      }
    });
    return [top, bottom];
  }

  // The one scale, never above 1, and each surface's move up or down, that keep every surface
  // within the height it had at the start.
  function fitting(M, cam) {
    var e = cam.elevation * RAD, scale = 1;
    var reach = M.surfaces.map(function(S) { return height(S, e); });
    M.surfaces.forEach(function(S, k) {
      var need = reach[k][0] - reach[k][1], have = S.top0 - S.bottom0;
      if (need > have) scale = Math.min(scale, have / need);
    });
    var shift = M.surfaces.map(function(S, k) {
      return Math.min(Math.max(0, S.bottom0 - scale * reach[k][1]), S.top0 - scale * reach[k][0]);
    });
    return { scale: scale, shift: shift };
  }

  /* The cone of the surface's scene that a ray from (qx, qy, qz) toward the camera meets
     farthest along it, further than eps, or with `any` the first found, or -1 where it meets
     none. Along the ray P + t v a cone's circle at height z1 + u dz has radius r1 + u dr, which
     gives a quadratic in u, the generator's Scene.reach(); when the ray runs level it keeps its
     height, and u is that height's place on the cone. */
  var reached = 0;   // how far along its ray cast() last met its cone
  function cast(S, qx, qy, qz, v, eps, any) {
    var R1 = S.r1, DR = S.dr, Z1 = S.z1, DZ = S.dz, m = R1.length, i, u, t, best = -1, far = eps;
    var vx = v[0], vy = v[1], vz = v[2], A = vx * vx + vy * vy;
    var Bp = 2 * (qx * vx + qy * vy), C = qx * qx + qy * qy;
    if (Math.abs(vz) > 1e-6) {
      var ivz = 1 / vz;
      for (i = 0; i < m; i++) {
        var alpha = (Z1[i] - qz) * ivz, beta = DZ[i] * ivz;
        if (alpha <= far && alpha + beta <= far) continue;
        var r1 = R1[i], dr = DR[i];
        var qa = A * beta * beta - dr * dr;
        var qb = Bp * beta + 2 * A * alpha * beta - 2 * r1 * dr;
        var qc = C + Bp * alpha + A * alpha * alpha - r1 * r1;
        var disc = qb * qb - 4 * qa * qc;
        if (disc < 0) continue;
        var u1, u2;
        if (Math.abs(qa) < 1e-14 * Math.max(1, Math.abs(qb))) {
          u1 = u2 = -qc / qb;
        } else {
          var root = Math.sqrt(disc);
          u1 = (-qb - root) / (2 * qa);
          u2 = (-qb + root) / (2 * qa);
        }
        for (var k = 0; k < 2; k++) {
          u = k ? u2 : u1;
          t = alpha + beta * u;
          if (u >= 0 && u <= 1 && t > far) {
            if (any) return i;
            best = i;
            far = reached = t;
          }
        }
      }
      return best;
    }
    for (i = 0; i < m; i++) {
      if (DZ[i] === 0) continue;
      u = (qz - Z1[i]) / DZ[i];
      if (u < 0 || u > 1) continue;
      var rho = R1[i] + u * DR[i], d = Bp * Bp - 4 * A * (C - rho * rho);
      if (d < 0) continue;
      t = (-Bp + Math.sqrt(d)) / (2 * A);
      if (t > far) {
        if (any) return i;
        best = i;
        far = reached = t;
      }
    }
    return best;
  }

  function hidden(S, qx, qy, qz, v, eps) { return cast(S, qx, qy, qz, v, eps, true) >= 0; }

  /* Whether a point of the outline is hidden. The line of sight only grazes the surface
     there, and along it the neighbouring cones lie as near the line as the rounding of the
     published points, a part in 10^7 of a piece, which for a straight profile such as the
     cone's is enough to hide every other point of its outline. So the point is taken a
     hundred thousandth of the drawing off the surface on either side, along its normal, and is
     hidden only if both are: on the side the surface folds toward, the fold hides it, and on
     the other only what truly lies in front does. */
  function outlined(S, P, v, M) {
    var d = 1e-5 * M.size, n = P.normal;
    return hidden(S, P[0] + d * n[0], P[1] + d * n[1], P[2] + d * n[2], v, M.eps) &&
           hidden(S, P[0] - d * n[0], P[1] - d * n[1], P[2] - d * n[2], v, M.eps);
  }

  // Ramer-Douglas-Peucker, as null_rays.thin(): the fewest points within tol of the polyline.
  function thin(P, tol) {
    var n = P.length;
    if (n < 3) return P;
    var keep = new Uint8Array(n), stack = [[0, n - 1]];
    keep[0] = keep[n - 1] = 1;
    while (stack.length) {
      var ab = stack.pop(), a = ab[0], b = ab[1];
      if (b <= a + 1) continue;
      var cx = P[b][0] - P[a][0], cy = P[b][1] - P[a][1], length = Math.sqrt(cx * cx + cy * cy);
      var best = -1, far = -1;
      for (var i = a + 1; i < b; i++) {
        var rx = P[i][0] - P[a][0], ry = P[i][1] - P[a][1];
        var d = length < 1e-12 ? Math.sqrt(rx * rx + ry * ry) : Math.abs(rx * cy - ry * cx) / length;
        if (d > far) { far = d; best = i; }
      }
      if (far > tol) {
        keep[best] = 1;
        stack.push([a, best], [best, b]);
      }
    }
    return P.filter(function(_, i) { return keep[i]; });
  }

  /* The polylines where F = 0 on a grid of ni rows and nj columns, F[i * nj + j], found by
     marching squares with a zero counted as positive and a saddle settled by the mean of its
     corners. Each polyline is a list of [i, j] in grid units, and a closed one repeats its first
     point at its end. */
  function isolines(F, ni, nj) {
    var H = ni * nj, at = new Map(), links = new Map();
    function place(key, i0, j0, i1, j1, p, q) {
      if (!at.has(key)) {
        var w = p / (p - q);
        at.set(key, [i0 + w * (i1 - i0), j0 + w * (j1 - j0)]);
      }
      return key;
    }
    function link(a, b) {
      if (!links.has(a)) links.set(a, []);
      if (!links.has(b)) links.set(b, []);
      links.get(a).push(b);
      links.get(b).push(a);
    }
    for (var i = 0; i + 1 < ni; i++) {
      for (var j = 0; j + 1 < nj; j++) {
        var a = F[i * nj + j], b = F[i * nj + j + 1], c = F[(i + 1) * nj + j + 1], d = F[(i + 1) * nj + j];
        var ia = a >= 0, ib = b >= 0, ic = c >= 0, id = d >= 0;
        if (ia === ib && ib === ic && ic === id) continue;
        var T = ia !== ib ? place(i * nj + j, i, j, i, j + 1, a, b) : -1;
        var R = ib !== ic ? place(H + i * nj + j + 1, i, j + 1, i + 1, j + 1, b, c) : -1;
        var B = id !== ic ? place((i + 1) * nj + j, i + 1, j, i + 1, j + 1, d, c) : -1;
        var L = ia !== id ? place(H + i * nj + j, i, j, i + 1, j, a, d) : -1;
        if (T >= 0 && R >= 0 && B >= 0 && L >= 0) {
          if ((a + b + c + d >= 0) === ia) { link(T, R); link(B, L); } else { link(T, L); link(R, B); }
        } else {
          var ends = [T, R, B, L].filter(function(e) { return e >= 0; });
          link(ends[0], ends[1]);
        }
      }
    }
    var out = [], done = new Set();
    function walk(start) {
      var line = [start], prev = -1, cur = start;
      done.add(start);
      for (;;) {
        var next = links.get(cur).filter(function(n) { return n !== prev && !done.has(n); })[0];
        if (next === undefined) break;
        line.push(next);
        done.add(next);
        prev = cur;
        cur = next;
      }
      if (line.length > 2 && links.get(cur).indexOf(start) >= 0) line.push(start);
      return line.map(function(key) { return at.get(key); });
    }
    links.forEach(function(n, key) { if (n.length === 1 && !done.has(key)) out.push(walk(key)); });
    links.forEach(function(n, key) { if (!done.has(key)) out.push(walk(key)); });
    return out;
  }

  /* Where a piece turns edge on to the camera, as the generator's outline() finds it: on the
     circle of each point of the profile, the angles at which the normal (dz cos phi,
     dz sin phi, -drho) is square to the line of sight, cos(phi - a) = (drho/dz) tan e, with the
     tangent taken as numpy's gradient takes it, joined from point to point. A run of such
     points inside the piece turns at both ends, where its two sides meet; one that reaches an
     end of the piece runs on to it on each side. */
  function outline(p, cam) {
    var m = p.rho.length, a = cam.azimuth * RAD, e = cam.elevation * RAD, runs = [];
    var c = new Float64Array(m), ok = new Uint8Array(m), nr = new Float64Array(m), nz = new Float64Array(m);
    for (var i = 0; i < m; i++) {
      var lo = Math.max(0, i - 1), hi = Math.min(m - 1, i + 1), span = hi - lo;
      var drho = (p.rho[hi] - p.rho[lo]) / span, dz = (p.z[hi] - p.z[lo]) / span, length = Math.sqrt(drho * drho + dz * dz);
      c[i] = drho * Math.sin(e) / (dz * Math.cos(e));
      ok[i] = isFinite(c[i]) && Math.abs(c[i]) <= 1 && p.rho[i] > 0 ? 1 : 0;
      nr[i] = dz / length;
      nz[i] = -drho / length;
    }
    // Each point carries the unit normal of the surface there, which outlined() needs.
    function side(i, j, s) {
      var out = [];
      for (var k = i; k < j; k++) {
        var phi = a + s * Math.acos(Math.max(-1, Math.min(1, c[k]))), cp = Math.cos(phi), sp = Math.sin(phi);
        var P = [p.rho[k] * cp, p.rho[k] * sp, p.z[k]];
        P.normal = [nr[k] * cp, nr[k] * sp, nz[k]];
        out.push(P);
      }
      return out;
    }
    for (i = 0; i < m;) {
      if (!ok[i]) { i++; continue; }
      var j = i;
      while (j < m && ok[j]) j++;
      if (j - i > 1) {
        if (i > 0 && j < m) runs.push(side(i, j, 1).reverse().concat(side(i, j, -1)));
        else runs.push(side(i, j, 1), side(i, j, -1));
      }
      i = j;
    }
    return runs;
  }

  // A polyline of n points on each segment of P, as the generator's densify().
  function densify(P, n) {
    var out = [];
    for (var i = 0; i + 1 < P.length; i++) {
      for (var k = 0; k < n; k++) {
        var s = k / n;
        out.push([P[i][0] + s * (P[i + 1][0] - P[i][0]), P[i][1] + s * (P[i + 1][1] - P[i][1]),
                  P[i][2] + s * (P[i + 1][2] - P[i][2])]);
      }
    }
    out.push(P[P.length - 1]);
    return out;
  }

  function circle(rho, z, n) {
    var out = [];
    for (var k = 0; k <= n; k++) {
      var phi = 2 * Math.PI * k / n;
      out.push([rho * Math.cos(phi), rho * Math.sin(phi), z]);
    }
    return out;
  }

  function meridian(p, phi) {
    var c = Math.cos(phi), s = Math.sin(phi), out = [];
    for (var i = 0; i < p.rho.length; i++) out.push([p.rho[i] * c, p.rho[i] * s, p.z[i]]);
    return out;
  }

  // How finely a figure is drawn: as the generator draws it, or coarser while a drag goes on.
  var FINE = { densify: 4, ring: 720, facets: 128, grid: 420 };
  var QUICK = { densify: 2, ring: 240, facets: 72, grid: 240 };

  /* The figure's layers and labels at the camera (azimuth, elevation), in the published
     figure's form: layers painted in order, fills first, then every line hidden and every line
     seen, each in the order of the legend, then points, and for each published label where it
     stands and whether it is shown. */
  function draw(M, azimuth, elevation, quick, sizes) {
    var Q = quick ? QUICK : FINE, cam = camera(azimuth, elevation), fit = fitting(M, cam);
    var s = fit.scale, ce = Math.cos(elevation * RAD), R = cam.right, U = cam.up, V = cam.toward;
    var lines = [], dots = [], outlines = M.surfaces.map(function() { return []; });
    var tol = 4e-4 * M.size;

    M.surfaces.forEach(function(S, k) {
      var X0 = S.at[0], Y0 = S.at[1] + fit.shift[k] - s * S.zc * ce;
      function page(P) { return [X0 + s * (P[0] * R[0] + P[1] * R[1]), Y0 + s * (P[0] * U[0] + P[1] * U[1] + P[2] * U[2])]; }
      // A line on the surface, split into the parts seen and the parts hidden, cut halfway
      // between neighbouring points of opposite kinds, as Figure.line() splits it.
      function line(cls, P) {
        var n = P.length, hid = new Uint8Array(n), XY = new Array(n);
        for (var i = 0; i < n; i++) {
          hid[i] = (P[i].normal ? outlined(S, P[i], V, M) : hidden(S, P[i][0], P[i][1], P[i][2], V, M.eps)) ? 1 : 0;
          XY[i] = page(P[i]);
        }
        var start = 0;
        for (var c = 0; c < n; c++) {
          if (c < n - 1 && hid[c + 1] === hid[c]) continue;
          var run = XY.slice(start, c + 1);
          if (c < n - 1) run.push([(XY[c][0] + XY[c + 1][0]) / 2, (XY[c][1] + XY[c + 1][1]) / 2]);
          if (start > 0) run.unshift([(XY[start - 1][0] + XY[start][0]) / 2, (XY[start - 1][1] + XY[start][1]) / 2]);
          if (run.length > 1) {
            var L = { cls: cls + (hid[start] ? '-far' : ''), points: thin(run, tol) };
            lines.push(L);
            if (L.cls === 'outline') outlines[k].push(L.points);
          }
          start = c + 1;
        }
      }
      S.pieces.forEach(function(p) {
        for (var j = 0; j < M.turn.meridians; j += p.reference ? 2 : 1) {
          line(p.reference ? 'reference' : 'meridian', densify(meridian(p, 2 * Math.PI * j / M.turn.meridians), Q.densify));
        }
        outline(p, cam).forEach(function(run) { line(p.reference ? 'reference' : 'outline', run); });
        p.ends.forEach(function(kind, end) {
          var i = end ? p.rho.length - 1 : 0;
          if ((kind === 'edge' || kind === 'stops') && !p.reference) line('outline', circle(p.rho[i], p.z[i], Q.ring));
        });
      });
      S.rings.forEach(function(ring) {
        if (ring.rho > 0) line(ring['class'], circle(ring.rho, ring.z, Q.ring));
        else dots.push({ kind: 'point', 'class': ring['class'], at: page([0, 0, ring.z]) });
      });
      M.turn.marks.forEach(function(mark) {
        if (mark.surface === k) line(mark['class'], densify(meridian(S.byId[mark.piece], mark.phi), Q.densify));
      });
    });

    var layers = M.figure.layers.filter(function(L) { return L.flat && L.kind === 'fill'; });
    layers = layers.concat(fills(M, cam, fit, Q));
    M.figure.layers.forEach(function(L) { if (L.flat && L.kind === 'line') lines.push({ cls: L['class'], points: L.points }); });
    function rank(L) { var r = M.order[L.cls.replace(/-far$/, '')]; return r === undefined ? 99 : r; }
    [true, false].forEach(function(far) {
      lines.filter(function(L) { return /-far$/.test(L.cls) === far; })
        .sort(function(a, b) { return rank(a) - rank(b); })
        .forEach(function(L) { layers.push({ kind: 'line', 'class': L.cls, points: L.points }); });
    });
    layers = layers.concat(dots);
    return { layers: layers, labels: labels(M, cam, fit, outlines, sizes), scale: s };
  }

  /* Where each fill class is the surface nearest the camera, as polygons on a grid of the
     page, as the generator's Figure.fills_seen() finds them. Every cone of every piece of the
     slice is cut into facets and rasterized with its depth, the nearest kept at each point of
     the grid; then wherever a point's nearest piece differs from a neighbour's, where facets
     could have it wrong, it is found again by casting a ray through the cones themselves, and
     so is every neighbour of a point that changes, so that the regions are the generator's.
     Each class's region is traced round, and a piece whose class is not tinted is left clear
     but hides what lies behind it. */
  function fills(M, cam, fit, Q) {
    var box = M.box, step = M.size / Q.grid, x0 = box[0] - step, y0 = box[2] - step;
    var nx = Math.ceil((box[1] - box[0]) / step) + 3, ny = Math.ceil((box[3] - box[2]) / step) + 3;
    var depth = new Float64Array(nx * ny).fill(-Infinity), code = new Uint8Array(nx * ny);
    var n = Q.facets, cs = new Float64Array(n + 1), sn = new Float64Array(n + 1);
    for (var j = 0; j <= n; j++) { cs[j] = Math.cos(2 * Math.PI * j / n); sn[j] = Math.sin(2 * Math.PI * j / n); }
    var s = fit.scale, ce = Math.cos(cam.elevation * RAD), R = cam.right, U = cam.up, V = cam.toward;
    var places = M.surfaces.map(function(S, k) { return [S.at[0], S.at[1] + fit.shift[k] - s * S.zc * ce]; });

    function triangle(ax, ay, ad, bx, by, bd, cx, cy, cd, c) {
      var area = (bx - ax) * (cy - ay) - (cx - ax) * (by - ay);
      if (area === 0) return;
      var minx = Math.max(0, Math.ceil(Math.min(ax, bx, cx))), maxx = Math.min(nx - 1, Math.floor(Math.max(ax, bx, cx)));
      var miny = Math.max(0, Math.ceil(Math.min(ay, by, cy))), maxy = Math.min(ny - 1, Math.floor(Math.max(ay, by, cy)));
      for (var py = miny; py <= maxy; py++) {
        for (var px = minx; px <= maxx; px++) {
          var wa = (bx - px) * (cy - py) - (cx - px) * (by - py);
          var wb = (cx - px) * (ay - py) - (ax - px) * (cy - py);
          var wc = area - wa - wb;
          if (area > 0 ? (wa < 0 || wb < 0 || wc < 0) : (wa > 0 || wb > 0 || wc > 0)) continue;
          var d = (wa * ad + wb * bd + wc * cd) / area, at = py * nx + px;
          if (d > depth[at]) { depth[at] = d; code[at] = c; }
        }
      }
    }

    M.surfaces.forEach(function(S, k) {
      var X0 = places[k][0], Y0 = places[k][1];
      S.pieces.forEach(function(p) {
        if (p.reference) return;
        var m = p.rho.length, X = new Float64Array(m * (n + 1)), Y = new Float64Array(m * (n + 1)), D = new Float64Array(m * (n + 1));
        for (var i = 0; i < m; i++) {
          for (var j = 0; j <= n; j++) {
            var px = p.rho[i] * cs[j], py = p.rho[i] * sn[j], pz = p.z[i], at = i * (n + 1) + j;
            X[at] = (X0 + s * (px * R[0] + py * R[1]) - x0) / step;
            Y[at] = (Y0 + s * (px * U[0] + py * U[1] + pz * U[2]) - y0) / step;
            D[at] = px * V[0] + py * V[1] + pz * V[2];
          }
        }
        for (i = 0; i + 1 < m; i++) {
          for (j = 0; j < n; j++) {
            var a = i * (n + 1) + j, b = a + 1, d = a + n + 1, e = d + 1;
            triangle(X[a], Y[a], D[a], X[b], Y[b], D[b], X[d], Y[d], D[d], p.code);
            triangle(X[b], Y[b], D[b], X[e], Y[e], D[e], X[d], Y[d], D[d], p.code);
          }
        }
      });
    });

    // The nearest piece at a point of the grid, by a ray from far behind every surface.
    function exact(at) {
      var X = x0 + (at % nx) * step, Y = y0 + Math.floor(at / nx) * step, best = 0, near = -Infinity;
      M.surfaces.forEach(function(S, k) {
        var a = (X - places[k][0]) / s, b = (Y - places[k][1]) / s, back = 4 * M.size / s;
        var i = cast(S, a * R[0] + b * U[0] - back * V[0], a * R[1] + b * U[1] - back * V[1], b * U[2] - back * V[2],
                     V, M.eps, false);
        if (i >= 0 && reached - back > near) { near = reached - back; best = S.code[i]; }
      });
      return best;
    }
    var known = new Uint8Array(nx * ny), queue = [];
    function around(at, f) {
      var x = at % nx;
      if (x > 0) f(at - 1);
      if (x < nx - 1) f(at + 1);
      if (at >= nx) f(at - nx);
      if (at < nx * (ny - 1)) f(at + nx);
    }
    for (var at = 0; at < nx * ny; at++) {
      around(at, function(b) { if (code[b] !== code[at]) queue.push(at); });
    }
    while (queue.length) {
      var at2 = queue.pop();
      if (known[at2]) continue;
      known[at2] = 1;
      var c = exact(at2);
      if (c !== code[at2]) {
        code[at2] = c;
        around(at2, function(b) { if (!known[b]) queue.push(b); });
      }
    }

    var out = [], tol = 5e-4 * M.size;
    M.classes.forEach(function(cls, index) {
      var F = new Float64Array(nx * ny), any = false;
      for (var at = 0; at < nx * ny; at++) {
        var on = code[at] === 2 + index;
        F[at] = on ? 0.5 : -0.5;
        if (on) any = true;
      }
      if (!any) return;
      var rings = isolines(F, ny, nx).map(function(line) {
        return thin(line.map(function(g) { return [x0 + g[1] * step, y0 + g[0] * step]; }), tol);
      }).filter(function(ring) { return ring.length > 3; });
      if (!rings.length) return;
      var layer = { kind: 'fill', 'class': cls, points: rings[0] };
      if (rings.length > 1) layer.holes = rings.slice(1);
      out.push(layer);
    });
    return out;
  }

  /* A generous width of a TeX label in ems, and the box it takes, as the generator's
     label_width() and label_box() give them, in the units of the figure's box, or with `size`
     the box of that width and height, as a page measures a label it has set. */
  var LAB = { lab: 15, small: 13 };
  var ANCHOR = { l: [0, -0.5], r: [-1, -0.5], t: [-0.5, 0], b: [-0.5, -1], c: [-0.5, -0.5],
                 tl: [0, 0], tr: [-1, 0], bl: [0, -1], br: [-1, -1] };
  var COMMANDS = ['\\chi', '\\ell', '\\pi', '\\phi', '\\eta', '\\theta', '\\sqrt', '\\frac', '\\infty', '\\mu', '\\delta'];
  function labelWidth(text) {
    var plain = text.replace(/\$/g, '').split('\\,').join(' ');
    COMMANDS.forEach(function(cmd) { plain = plain.split(cmd).join('x'); });
    plain = plain.replace(/[{}_^]/g, '');
    return 0.55 * Array.from(plain).length + 0.3;
  }
  function labelBox(L, at, unit, measured) {
    var size = (LAB[L['class']] || 13) * unit, w = labelWidth(L.text) * size, h = 1.25 * size;
    if (measured) { w = measured[0]; h = measured[1]; }
    var a = ANCHOR[L.anchor || 'c'];
    var x0 = at[0] + (L.dx || 0) * unit + a[0] * w, y0 = at[1] - (L.dy || 0) * unit - (a[1] + 1) * h;
    return [x0, y0, x0 + w, y0 + h];
  }

  /* A label that names a circle stands beside the end of the circle on the side the file
     gives, the point of it farthest right or left on the page, which is where the published
     label stands, and with `clear` past the outline where the outline crosses within `clear`
     above or below that end, as the generator's ring_label() sets it. It names its circle
     wherever the circle runs, in front of the surface or behind it, as the published labels
     do, and gives way only where it would overlap a label that stays put or a label before it
     in the file, measured by `sizes` where the page gives each label's width and height as it
     set it. Every other label stays where it is and is always shown. */
  function labels(M, cam, fit, outlines, sizes) {
    var s = fit.scale, ce = Math.cos(cam.elevation * RAD);
    var out = M.figure.labels.map(function(L) {
      var at = L.at;
      if (L.ring) {
        var k = L.ring.surface, S = M.surfaces[k], ring = S.rings[L.ring.ring], side = L.ring.side;
        at = [S.at[0] + s * side * ring.rho, S.at[1] + fit.shift[k] + s * (ring.z - S.zc) * ce];
        if (L.clear) {
          var band = L.clear;
          outlines[k].forEach(function(P) {
            [at[1] - band, at[1], at[1] + band].forEach(function(y) {
              for (var i = 0; i + 1 < P.length; i++) {
                var p = P[i], q = P[i + 1];
                if ((p[1] - y) * (q[1] - y) > 0 || p[1] === q[1]) continue;
                var x = p[0] + (y - p[1]) * (q[0] - p[0]) / (q[1] - p[1]);
                at = [side > 0 ? Math.max(at[0], x) : Math.min(at[0], x), at[1]];
              }
            });
          });
        }
      }
      return { at: at, shown: true };
    });
    var boxes = [];
    function boxOf(i) { return labelBox(M.figure.labels[i], out[i].at, M.unit, sizes && sizes[i]); }
    M.figure.labels.forEach(function(L, i) { if (!L.ring) boxes.push(boxOf(i)); });
    M.figure.labels.forEach(function(L, i) {
      if (!L.ring || !out[i].shown) return;
      var a = boxOf(i);
      out[i].shown = boxes.every(function(b) { return !(a[0] < b[2] && b[0] < a[2] && a[1] < b[3] && b[1] < a[3]); });
      if (out[i].shown) boxes.push(a);
    });
    return out;
  }

  var api = { camera: camera, prepare: prepare, draw: draw, fitting: fitting, hidden: hidden, isolines: isolines, labelBox: labelBox };
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  else root.MfsTurn = api;
})(this);
