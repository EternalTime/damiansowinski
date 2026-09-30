/* Turning a figure drawn in three dimensions.

   _layouts/mfs.html draws a figure's published drawing, which a generator projected once from a
   fixed camera, until the reader drags it, and from then on draws the figure again from its
   pieces in three dimensions at the camera the drag has reached, with prepare() and draw()
   below. Two kinds of figure turn: an embedding diagram, from the surfaces of its view, which
   _tools/derivations/embedding.py writes, and a figure of light cones in three dimensions, from
   the `turn` of a spacetime diagram's figure, which _tools/derivations/projections.py writes.
   The same hand turns both, as dragged(), keyed() and turned() say, and both keep to the box
   they were published in. Everything here is geometry and touches no page, so the tests run it
   in Node against the published figures; _tools/README.md, "Turning the figure" and "Turning a
   figure of light cones", is the definition the application follows as well.

   An embedding diagram is drawn by the rules the generator draws it by: each piece's meridians,
   its outline where it turns edge on to the camera, the rim at an end where the drawing stops,
   its marked circles and marked meridians, the curves and points marked on a surface that are no
   circles, as a ring of free particles on a flat plane, every line split where a surface hides
   it into the part seen and the part hidden, which carries -far, and each tinted piece filled
   where it is the surface nearest the camera. What hides a line is found as the generator finds
   it, by casting a ray from each of its points toward the camera through the truncated cones
   between neighbouring circles of every profile, which is exact for the surface the points
   describe, and a point of the outline, where the line of sight only grazes the surface, is
   judged a little off it on either side. The tint is found on a grid of the page, as the
   generator's is, from the same cones cut into facets and cast through exactly wherever two
   regions meet.

   Each surface turns about the point of its axis halfway up its drawing, which stays where it
   stood on the page. Its width on the page never changes as it turns, a surface of revolution
   being as wide as its widest circle from every side, but its height does as it tilts, so the
   surfaces are drawn smaller, all at one scale, only as far as a tilt needs more height than a
   surface had at the start, and moved up or down only as far as that height needs. So the
   figure keeps its box and its labels their room, and at the start it is the published figure.

   A grid piece, a quantity drawn as a height over a plane, is drawn by the same rules from its
   flat triangles: the lines of its grid that the figure names, its rim, and its outline where the
   normal at its nodes turns square to the line of sight. What hides a point is found from the
   triangles that cover the point of the page it falls on, which is exact for triangles seen along
   one direction, and a point on a height is judged from just above and just below it, the two
   sides of a height. A grid need not look the same from every side, so it is kept within the
   width it had at the start as well as the height, a rectangle turned toward its diagonal drawn
   smaller.

   A figure of light cones is drawn as drawCones() says. */
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

  // What drawSurfaces() needs of an embedding view that does not depend on the camera.
  function prepareSurfaces(view) {
    var figure = view.figure, turn = figure.turn, box = figure.box;
    var size = Math.max(box[1] - box[0], box[3] - box[2]);
    var e0 = figure.camera.elevation * RAD, a0 = figure.camera.azimuth * RAD, classes = [];
    Object.keys(turn.tint).forEach(function(k) { if (classes.indexOf(turn.tint[k]) < 0) classes.push(turn.tint[k]); });
    classes.sort();
    // A movie's frames are drawn one at a time, each on its axis at the one place.
    var movie = view.movie || null;
    var surfaces = (movie ? movie.frames : view.surfaces).map(function(s, k) {
      var low = Infinity, high = -Infinity, byId = {}, grids = [];
      var pieces = s.pieces.filter(function(p) {
        if (!p.grid) return true;
        var g = gridPiece(p);
        for (var i = 0; i < g.Z.length; i++) {
          if (g.Z[i] < low) low = g.Z[i];
          if (g.Z[i] > high) high = g.Z[i];
        }
        g.code = turn.tint[g.cls] ? 2 + classes.indexOf(turn.tint[g.cls]) : 1;
        byId[p.id] = g;
        grids.push(g);
        return false;
      }).map(function(p) {
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
      var S = { pieces: pieces, grids: grids, byId: byId, rings: s.rings, curves: s.curves || [], dots: s.dots || [],
                axis: s.axis || null, low: low, high: high, zc: (low + high) / 2,
                r1: new Float64Array(r1), dr: new Float64Array(dr), z1: new Float64Array(z1), dz: new Float64Array(dz),
                code: new Uint8Array(code) };
      S.origin = turn.origins[movie ? 0 : k];
      return S;
    });
    // Every frame of a movie turns about one centre, halfway between the lowest and the highest
    // point any frame reaches, and keeps within the height and width every frame takes together,
    // so the frames keep one size and one place as they change and as the reader turns them.
    if (movie) {
      var low = Math.min.apply(null, surfaces.map(function(S) { return S.low; }));
      var high = Math.max.apply(null, surfaces.map(function(S) { return S.high; }));
      surfaces.forEach(function(S) { S.zc = (low + high) / 2; });
    }
    surfaces.forEach(function(S) {
      S.at = [S.origin[0], S.origin[1] + S.zc * Math.cos(e0)];
      var h = height(S, e0, a0);
      S.top0 = h[0];
      S.bottom0 = h[1];
      if (S.grids.length) {
        var w = width(S, a0);
        S.left0 = w[0];
        S.right0 = w[1];
      }
    });
    if (movie) {
      var top0 = Math.max.apply(null, surfaces.map(function(S) { return S.top0; }));
      var bottom0 = Math.min.apply(null, surfaces.map(function(S) { return S.bottom0; }));
      surfaces.forEach(function(S) { S.top0 = top0; S.bottom0 = bottom0; });
    }
    var order = {};
    figure.legend.forEach(function(item, i) { order[item[1]] = i; });
    return { kind: 'surfaces', figure: figure, turn: turn, box: box, size: size, surfaces: surfaces, order: order, classes: classes,
             eps: 1e-7 * size, unit: (box[1] - box[0]) / 560, movie: movie, frame: 0 };
  }

  // Which frame of a movie draw() draws, the first until another is chosen.
  function frame(M, k) { M.frame = k; }

  /* A grid piece: a quantity drawn as a height over a plane, sampled on a grid and drawn as flat
     triangles, as _tools/README.md defines it. Its nodes in the surface's own frame, the node of
     row i of u and column j of v at i * n + j, a polar grid's X and Y being u cos v and u sin v;
     each cell cut into two triangles along its diagonal from (i, j) to (i + 1, j + 1), listed
     first every triangle (i, j), (i + 1, j), (i + 1, j + 1) and then every (i, j), (i + 1, j + 1),
     (i, j + 1), a polar grid's last column joining its first, as the generator lists them; and at
     each node the normal pointing up, the sum of the normals of the triangles that meet there,
     each as long as twice its area, as the generator's grid_normals() takes it. */
  function gridPiece(p) {
    var G = p.grid, m = G.u.length, n = G.v.length, N = m * n, ellipses = G.frame === 'ellipses';
    var wrap = G.frame === 'polar' || (ellipses && !G.open);
    var X = new Float64Array(N), Y = new Float64Array(N), Z = new Float64Array(N), i, j, k;
    for (i = 0; i < m; i++) {
      for (j = 0; j < n; j++) {
        k = i * n + j;
        if (ellipses) {
          X[k] = G.a[i] * Math.cos(G.v[j]);
          Y[k] = G.b[i] * Math.sin(G.v[j]);
          Z[k] = G.z[i];
        } else {
          X[k] = G.frame === 'polar' ? G.u[i] * Math.cos(G.v[j]) : G.u[i];
          Y[k] = G.frame === 'polar' ? G.u[i] * Math.sin(G.v[j]) : G.v[j];
          Z[k] = G.z[i][j];
        }
      }
    }
    var cols = wrap ? n : n - 1, half = (m - 1) * cols, T = 2 * half;
    var A = new Int32Array(T), B = new Int32Array(T), C = new Int32Array(T);
    for (i = 0; i + 1 < m; i++) {
      for (j = 0; j < cols; j++) {
        var a = i * n + j, b = (i + 1) * n + j, c = (i + 1) * n + (j + 1) % n, d = i * n + (j + 1) % n, t = i * cols + j;
        A[t] = a; B[t] = b; C[t] = c;
        A[half + t] = a; B[half + t] = c; C[half + t] = d;
      }
    }
    var nx = new Float64Array(N), ny = new Float64Array(N), nz = new Float64Array(N);
    for (t = 0; t < T; t++) {
      var ex = X[B[t]] - X[A[t]], ey = Y[B[t]] - Y[A[t]], ez = Z[B[t]] - Z[A[t]];
      var fx = X[C[t]] - X[A[t]], fy = Y[C[t]] - Y[A[t]], fz = Z[C[t]] - Z[A[t]];
      var cx = ey * fz - ez * fy, cy = ez * fx - ex * fz, cz = ex * fy - ey * fx;
      // A height's normal points up; a stack of ellipses turns every triangle one way, rows up
      // and columns round, and its normal points away from the axis.
      if (ellipses ? true : cz < 0) { cx = -cx; cy = -cy; cz = -cz; }
      [A[t], B[t], C[t]].forEach(function(q) { nx[q] += cx; ny[q] += cy; nz[q] += cz; });
    }
    // A node that no triangle with any area meets, as a cone's apex at the end of its cut, takes
    // the normal of the node before it in its row, or after it at the row's start.
    for (i = 0; i < m; i++) {
      for (j = 0; j < n; j++) {
        k = i * n + j;
        if (nx[k] === 0 && ny[k] === 0 && nz[k] === 0) {
          var o = j > 0 ? k - 1 : k + 1;
          nx[k] = nx[o]; ny[k] = ny[o]; nz[k] = nz[o];
        }
      }
    }
    for (k = 0; k < N; k++) {
      var L = Math.sqrt(nx[k] * nx[k] + ny[k] * ny[k] + nz[k] * nz[k]);
      nx[k] /= L; ny[k] /= L; nz[k] /= L;
    }
    // What hides a point on it is measured in the grid's own width across the plane, as the
    // generator measures it.
    var lo = [Infinity, Infinity], hi = [-Infinity, -Infinity];
    for (k = 0; k < N; k++) {
      lo[0] = Math.min(lo[0], X[k]); hi[0] = Math.max(hi[0], X[k]);
      lo[1] = Math.min(lo[1], Y[k]); hi[1] = Math.max(hi[1], Y[k]);
    }
    var extent = Math.max(hi[0] - lo[0], hi[1] - lo[1]);
    return { id: p.id, cls: p['class'], grid: true, ellipses: ellipses, lines: G.lines, m: m, n: n, wrap: wrap, X: X, Y: Y, Z: Z,
             A: A, B: B, C: C, nx: nx, ny: ny, nz: nz, eps: 1e-7 * extent, lift: 1e-5 * extent };
  }

  function node(g, i, j) { var k = i * g.n + j; return [g.X[k], g.Y[k], g.Z[k]]; }

  // Row i of a grid piece, round and back to its first node where the grid runs round.
  function gridRow(g, i) {
    var out = [];
    for (var j = 0; j < g.n; j++) out.push(node(g, i, j));
    if (g.wrap) out.push(node(g, i, 0));
    return out;
  }
  function gridColumn(g, j) {
    var out = [];
    for (var i = 0; i < g.m; i++) out.push(node(g, i, j));
    return out;
  }
  // The edge of a grid piece, closed: a polar grid's outermost row, or round a rectangle.
  function gridRim(g) {
    if (g.wrap) return gridRow(g, g.m - 1);
    var out = [], i, j;
    for (j = 0; j + 1 < g.n; j++) out.push(node(g, 0, j));
    for (i = 0; i + 1 < g.m; i++) out.push(node(g, i, g.n - 1));
    for (j = g.n - 1; j > 0; j--) out.push(node(g, g.m - 1, j));
    for (i = g.m - 1; i > 0; i--) out.push(node(g, i, 0));
    out.push(out[0]);
    return out;
  }

  /* A grid piece seen from a camera: each node on the page and its depth, and the triangles that
     cover any area of the page sorted into a grid of bins over the page by their boxes, about as
     many bins as triangles, as the generator's Facets sorts them. A ray toward the camera from a
     point meets a triangle exactly where the triangle covers the point of the page it falls on,
     after the difference of their depths. */
  function facetFrame(g, cam) {
    var R = cam.right, U = cam.up, V = cam.toward, N = g.X.length, T = g.A.length, k, t;
    var px = new Float64Array(N), py = new Float64Array(N), pd = new Float64Array(N);
    for (k = 0; k < N; k++) {
      px[k] = g.X[k] * R[0] + g.Y[k] * R[1];
      py[k] = g.X[k] * U[0] + g.Y[k] * U[1] + g.Z[k] * U[2];
      pd[k] = g.X[k] * V[0] + g.Y[k] * V[1] + g.Z[k] * V[2];
    }
    var area = new Float64Array(T), biggest = 0;
    for (t = 0; t < T; t++) {
      var a = g.A[t], b = g.B[t], c = g.C[t];
      area[t] = (px[b] - px[a]) * (py[c] - py[a]) - (px[c] - px[a]) * (py[b] - py[a]);
      biggest = Math.max(biggest, Math.abs(area[t]));
    }
    var floor = 1e-15 * Math.max(1, biggest), kept = [];
    for (t = 0; t < T; t++) if (Math.abs(area[t]) > floor) kept.push(t);
    var nb = Math.max(1, Math.floor(Math.sqrt(Math.max(kept.length, 1))));
    var lo = kept.map(function(t) {
      return [Math.min(px[g.A[t]], px[g.B[t]], px[g.C[t]]), Math.min(py[g.A[t]], py[g.B[t]], py[g.C[t]])];
    });
    var hi = kept.map(function(t) {
      return [Math.max(px[g.A[t]], px[g.B[t]], px[g.C[t]]), Math.max(py[g.A[t]], py[g.B[t]], py[g.C[t]])];
    });
    var gx = Infinity, gy = Infinity, hx = -Infinity, hy = -Infinity;
    lo.forEach(function(q) { gx = Math.min(gx, q[0]); gy = Math.min(gy, q[1]); });
    hi.forEach(function(q) { hx = Math.max(hx, q[0]); hy = Math.max(hy, q[1]); });
    if (!kept.length) { gx = gy = 0; hx = hy = 1; }
    var sx = Math.max((hx - gx) / nb, 1e-12), sy = Math.max((hy - gy) / nb, 1e-12);
    function bin(v, g0, s) { return Math.min(nb - 1, Math.max(0, Math.floor((v - g0) / s))); }
    var count = new Int32Array(nb * nb + 1), span = kept.map(function(t, q) {
      return [bin(lo[q][0], gx, sx), bin(hi[q][0], gx, sx), bin(lo[q][1], gy, sy), bin(hi[q][1], gy, sy)];
    });
    span.forEach(function(r) {
      for (var y = r[2]; y <= r[3]; y++) for (var x = r[0]; x <= r[1]; x++) count[y * nb + x + 1]++;
    });
    for (k = 0; k < nb * nb; k++) count[k + 1] += count[k];
    var list = new Int32Array(count[nb * nb]), fill = count.slice(0, nb * nb);
    span.forEach(function(r, q) {
      for (var y = r[2]; y <= r[3]; y++) for (var x = r[0]; x <= r[1]; x++) list[fill[y * nb + x]++] = kept[q];
    });
    return { px: px, py: py, pd: pd, area: area, nb: nb, gx: gx, gy: gy, hx: hx, hy: hy, sx: sx, sy: sy, start: count, list: list };
  }

  /* The deepest reach of a grid piece's triangles over the point (sx, sy) of the page, their depth
     there less dq, above eps, or with `any` the first found; -Infinity where none lies over it. */
  function facetReach(g, F, sx, sy, dq, eps, any) {
    var best = -Infinity;
    if (!(sx >= F.gx && sx <= F.hx && sy >= F.gy && sy <= F.hy)) return best;
    // A point on the far edge of the bins falls in the last of them, as a triangle's box does.
    var x = Math.min(F.nb - 1, Math.floor((sx - F.gx) / F.sx)), y = Math.min(F.nb - 1, Math.floor((sy - F.gy) / F.sy));
    var b = y * F.nb + x, px = F.px, py = F.py, pd = F.pd;
    for (var q = F.start[b]; q < F.start[b + 1]; q++) {
      var t = F.list[q], a = g.A[t], bb = g.B[t], c = g.C[t], area = F.area[t];
      var wa = ((px[bb] - sx) * (py[c] - sy) - (px[c] - sx) * (py[bb] - sy)) / area;
      var wb = ((px[c] - sx) * (py[a] - sy) - (px[a] - sx) * (py[c] - sy)) / area;
      var wc = 1 - wa - wb;
      if (wa < -1e-12 || wb < -1e-12 || wc < -1e-12) continue;
      var d = wa * pd[a] + wb * pd[bb] + wc * pd[c] - dq;
      if (d > eps) {
        if (any) return d;
        if (d > best) best = d;
      }
    }
    return best;
  }

  /* Where a grid piece turns edge on to the camera, as the generator's grid_outline() finds it:
     the lines where the normal at the nodes, carried along each edge of the grid, is square to
     the line of sight, by marching squares over the grid, a polar grid's first column repeated
     after its last. Each point lies on an edge of the grid, and so on the triangles. */
  function gridOutline(g, cam) {
    var V = cam.toward, m = g.m, nj = g.wrap ? g.n + 1 : g.n, F = new Float64Array(m * nj);
    for (var i = 0; i < m; i++) {
      for (var jj = 0; jj < nj; jj++) {
        var k = i * g.n + jj % g.n;
        F[i * nj + jj] = g.nx[k] * V[0] + g.ny[k] * V[1] + g.nz[k] * V[2];
      }
    }
    return isolines(F, m, nj).map(function(line) {
      return line.map(function(at) {
        var i0 = Math.floor(at[0]), j0 = Math.floor(at[1]), i1 = Math.min(i0 + 1, m - 1), j1 = Math.min(j0 + 1, nj - 1);
        var wi = at[0] - i0, wj = at[1] - j0, P = [0, 0, 0];
        [[i0, j0, (1 - wi) * (1 - wj)], [i1, j0, wi * (1 - wj)], [i0, j1, (1 - wi) * wj], [i1, j1, wi * wj]].forEach(function(c) {
          var q = node(g, c[0], c[1] % g.n);
          P[0] += c[2] * q[0]; P[1] += c[2] * q[1]; P[2] += c[2] * q[2];
        });
        return P;
      });
    });
  }

  /* How far a surface reaches above and below its centre on the page at elevation e: a point
     of its profile draws a circle whose height on the page runs over (z - zc) cos e +- rho sin e,
     and along a cone between two points both are linear, so the points bound the whole. A grid
     piece is flat between its nodes, so its nodes bound it, seen from the azimuth a as well. */
  function height(S, e, a) {
    var ce = Math.cos(e), se = Math.abs(Math.sin(e)), top = -Infinity, bottom = Infinity;
    S.pieces.forEach(function(p) {
      for (var i = 0; i < p.rho.length; i++) {
        var y = (p.z[i] - S.zc) * ce, w = p.rho[i] * se;
        if (y + w > top) top = y + w;
        if (y - w < bottom) bottom = y - w;
      }
    });
    if (S.grids.length) {
      var u0 = -Math.sin(e) * Math.cos(a), u1 = -Math.sin(e) * Math.sin(a);
      S.grids.forEach(function(g) {
        for (var k = 0; k < g.Z.length; k++) {
          var y = g.X[k] * u0 + g.Y[k] * u1 + (g.Z[k] - S.zc) * ce;
          if (y > top) top = y;
          if (y < bottom) bottom = y;
        }
      });
    }
    return [top, bottom];
  }

  /* How far a surface with a grid piece reaches left and right of its axis on the page from the
     azimuth a, [left, right], both positive. A surface of revolution is as wide from every side;
     a grid over a rectangle is not. */
  function width(S, a) {
    var r0 = -Math.sin(a), r1 = Math.cos(a), left = 0, right = 0;
    S.grids.forEach(function(g) {
      for (var k = 0; k < g.X.length; k++) {
        var x = g.X[k] * r0 + g.Y[k] * r1;
        if (x > right) right = x;
        if (-x > left) left = -x;
      }
    });
    return [left, right];
  }

  // The one scale, never above 1, and each surface's move up or down, that keep every surface
  // within the height it had at the start.
  function fitting(M, cam) {
    var e = cam.elevation * RAD, a = cam.azimuth * RAD, scale = 1;
    var reach = M.surfaces.map(function(S) { return height(S, e, a); });
    if (M.movie) {
      // A movie's frames reach as far as all of them together.
      var all = [Math.max.apply(null, reach.map(function(r) { return r[0]; })), Math.min.apply(null, reach.map(function(r) { return r[1]; }))];
      reach = reach.map(function() { return all; });
      if (M.surfaces.some(function(S) { return S.grids.length; })) {
        var w = M.surfaces.map(function(S) { return S.grids.length ? width(S, a) : [0, 0]; });
        var left0 = Math.max.apply(null, M.surfaces.map(function(S) { return S.left0 || 0; }));
        var right0 = Math.max.apply(null, M.surfaces.map(function(S) { return S.right0 || 0; }));
        var left = Math.max.apply(null, w.map(function(x) { return x[0]; })), right = Math.max.apply(null, w.map(function(x) { return x[1]; }));
        if (left > left0) scale = Math.min(scale, left0 / left);
        if (right > right0) scale = Math.min(scale, right0 / right);
      }
      var need = all[0] - all[1], have = M.surfaces[0].top0 - M.surfaces[0].bottom0;
      if (need > have) scale = Math.min(scale, have / need);
      var move = Math.min(Math.max(0, M.surfaces[0].bottom0 - scale * all[1]), M.surfaces[0].top0 - scale * all[0]);
      return { scale: scale, shift: M.surfaces.map(function() { return move; }) };
    }
    M.surfaces.forEach(function(S, k) {
      var need = reach[k][0] - reach[k][1], have = S.top0 - S.bottom0;
      if (need > have) scale = Math.min(scale, have / need);
      if (S.grids.length) {
        var w = width(S, a);
        if (w[0] > S.left0) scale = Math.min(scale, S.left0 / w[0]);
        if (w[1] > S.right0) scale = Math.min(scale, S.right0 / w[1]);
      }
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

  function hidden(S, qx, qy, qz, v, eps) {
    if (cast(S, qx, qy, qz, v, eps, true) >= 0) return true;
    for (var k = 0; k < S.grids.length; k++) {
      var g = S.grids[k], F = g.frame, cam = S.camera;
      var sx = qx * cam.right[0] + qy * cam.right[1], sy = qx * cam.up[0] + qy * cam.up[1] + qz * cam.up[2];
      if (facetReach(g, F, sx, sy, qx * v[0] + qy * v[1] + qz * v[2], g.eps, true) > g.eps) return true;
    }
    return false;
  }

  /* Whether a point of the outline is hidden. The line of sight only grazes the surface
     there, and along it the neighbouring cones lie as near the line as the rounding of the
     published points, a part in 10^7 of a piece, which for a straight profile such as the
     cone's is enough to hide every other point of its outline. So the point is taken a
     hundred thousandth of the drawing off the surface on either side, along its normal, and is
     hidden only if both are: on the side the surface folds toward, the fold hides it, and on
     the other only what truly lies in front does. A point on a grid piece is taken a hundred
     thousandth of the grid's width above and below it, as the generator takes it. */
  function outlined(S, P, v, M) {
    var d = P.lift || 1e-5 * M.size, n = P.normal;
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

  /* The points of a line on a grid piece, each to be judged from a hundred thousandth of the
     grid's width to either side of it: above and below a height, and straight out from the axis
     and back for a stack of ellipses, which that crosses wherever it is not flat, or above and
     below on the axis itself, as the generator's sides() takes them. */
  var UP = [0, 0, 1];
  function upright(P, g) {
    P.forEach(function(q) {
      var r = Math.sqrt(q[0] * q[0] + q[1] * q[1]);
      q.normal = g.ellipses && r > 0 ? [q[0] / r, q[1] / r, 0] : UP;
      q.lift = g.lift;
    });
    return P;
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

  // A polyline with each segment cut into as few equal pieces as keep every piece shorter than
  // `most`, as the generator's finely().
  function finely(P, most) {
    var out = [P[0]];
    for (var i = 0; i + 1 < P.length; i++) {
      var a = P[i], b = P[i + 1], d = Math.sqrt(Math.pow(b[0] - a[0], 2) + Math.pow(b[1] - a[1], 2) + Math.pow(b[2] - a[2], 2));
      var n = Math.max(1, Math.ceil(d / most));
      for (var k = 1; k <= n; k++) {
        var t = k / n;
        out.push([a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1]), a[2] + t * (b[2] - a[2])]);
      }
    }
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
  function drawSurfaces(M, azimuth, elevation, quick, sizes) {
    var Q = quick ? QUICK : FINE, cam = camera(azimuth, elevation), fit = fitting(M, cam);
    var s = fit.scale, ce = Math.cos(elevation * RAD), R = cam.right, U = cam.up, V = cam.toward;
    var lines = [], dots = [], outlines = M.surfaces.map(function() { return []; });
    var tol = 4e-4 * M.size;

    M.surfaces.forEach(function(S, k) {
      if (M.movie && k !== M.frame) return;
      var X0 = S.at[0], Y0 = S.at[1] + fit.shift[k] - s * S.zc * ce;
      S.camera = cam;
      S.grids.forEach(function(g) { g.frame = facetFrame(g, cam); });
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
      // A grid piece: the lines of its grid that the figure names, its rim and its outline, every
      // point of each judged by the two points just above and below it, as the generator's
      // draw_grid() judges them, the two sides of a height being above and below it.
      (M.turn.grid || []).forEach(function(G) {
        if (G.surface !== k) return;
        var g = S.byId[G.piece];
        G.u.forEach(function(i) { line(G['class'], upright(densify(gridRow(g, i), Q.densify), g)); });
        G.v.forEach(function(j) { line(G['class'], upright(densify(gridColumn(g, j), Q.densify), g)); });
      });
      S.grids.forEach(function(g) {
        // A stack of ellipses draws only the lines it names itself, and its outline.
        (g.lines || []).forEach(function(G) {
          G.u.forEach(function(i) { line(G['class'], upright(densify(gridRow(g, i), Q.densify), g)); });
          G.v.forEach(function(j) { line(G['class'], upright(densify(gridColumn(g, j), Q.densify), g)); });
        });
        if (!g.ellipses) line('outline', upright(densify(gridRim(g), Q.densify), g));
        gridOutline(g, cam).forEach(function(run) { line('outline', upright(run, g)); });
      });
      // A curve marked on the surface is drawn through its own points, and back to the first
      // where it closes; a point marked on it is drawn wherever it stands, as the generator's are.
      S.curves.forEach(function(c) {
        var P = c.closed ? c.points.concat([c.points[0]]) : c.points;
        var on = S.byId[c.piece];
        if (on && on.grid) P = upright(P.map(function(q) { return q.slice(); }), on);
        line(c['class'], P);
      });
      S.dots.forEach(function(d) { dots.push({ kind: 'point', 'class': d['class'], at: page(d.at) }); });
      // The axis of time a stack of moments stands on, hidden where the surface lies in front.
      if (S.axis) line('axis', finely([[0, 0, S.axis.from], [0, 0, S.axis.to]], M.size / 360));
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
    return { layers: layers, slices: [], labels: labels(M, cam, fit, outlines, sizes), scale: s };
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
      if (M.movie && k !== M.frame) return;
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
      // A grid piece is its triangles already.
      S.grids.forEach(function(g) {
        var N = g.X.length, X = new Float64Array(N), Y = new Float64Array(N), D = new Float64Array(N);
        for (var q = 0; q < N; q++) {
          var px = g.X[q], py = g.Y[q], pz = g.Z[q];
          X[q] = (X0 + s * (px * R[0] + py * R[1]) - x0) / step;
          Y[q] = (Y0 + s * (px * U[0] + py * U[1] + pz * U[2]) - y0) / step;
          D[q] = px * V[0] + py * V[1] + pz * V[2];
        }
        for (var t = 0; t < g.A.length; t++) {
          var a = g.A[t], b = g.B[t], c = g.C[t];
          triangle(X[a], Y[a], D[a], X[b], Y[b], D[b], X[c], Y[c], D[c], g.code);
        }
      });
    });

    // The nearest piece at a point of the grid, by a ray from far behind every surface.
    function exact(at) {
      var X = x0 + (at % nx) * step, Y = y0 + Math.floor(at / nx) * step, best = 0, near = -Infinity;
      M.surfaces.forEach(function(S, k) {
        if (M.movie && k !== M.frame) return;
        var a = (X - places[k][0]) / s, b = (Y - places[k][1]) / s, back = 4 * M.size / s;
        var i = cast(S, a * R[0] + b * U[0] - back * V[0], a * R[1] + b * U[1] - back * V[1], b * U[2] - back * V[2],
                     V, M.eps, false);
        if (i >= 0 && reached - back > near) { near = reached - back; best = S.code[i]; }
        S.grids.forEach(function(g) {
          var d = facetReach(g, g.frame, a, b, 0, -Infinity, false);
          if (d > near) { near = d; best = g.code; }
        });
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
    function place(k) { return M.movie ? M.frame : k; }
    var out = M.figure.labels.map(function(L) {
      var at = L.at, S, X0, Y0, k;
      function page(P) {
        return [X0 + s * (P[0] * cam.right[0] + P[1] * cam.right[1]),
                Y0 + s * (P[0] * cam.up[0] + P[1] * cam.up[1] + P[2] * cam.up[2])];
      }
      if (L.curve || L.axis) {
        k = place((L.curve || L.axis).surface);
        S = M.surfaces[k];
        X0 = S.at[0];
        Y0 = S.at[1] + fit.shift[k] - s * S.zc * ce;
      }
      if (L.curve) {
        // Beside the end of the curve on its side, its point farthest right or left on the page,
        // the first of them where two tie.
        var P = S.curves[L.curve.curve].points, best = -Infinity;
        P.forEach(function(q) {
          var p = page(q), x = L.curve.side * p[0];
          if (x > best) { best = x; at = p; }
        });
      } else if (L.axis) {
        at = page([0, 0, S.axis.to]);
      }
      if (L.ring) {
        var k = place(L.ring.surface), S = M.surfaces[k], ring = S.rings[L.ring.ring], side = L.ring.side;
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
    return giveWay(M.figure.labels, out, M.unit, sizes, function(i) {
      var L = M.figure.labels[i];
      return !!(L.ring || L.curve || L.axis);
    });
  }

  /* Labels at their places `out`, each shown unless `moves` says it follows the figure round and
     it would overlap a label that stays put or a label before it that is shown, measured by
     `sizes` where the page gives each label's width and height as it set it. */
  function giveWay(labels, out, unit, sizes, moves) {
    var boxes = [];
    function boxOf(i) { return labelBox(labels[i], out[i].at, unit, sizes && sizes[i]); }
    labels.forEach(function(L, i) { if (!moves(i)) boxes.push(boxOf(i)); });
    labels.forEach(function(L, i) {
      if (!moves(i) || !out[i].shown) return;
      var a = boxOf(i);
      out[i].shown = boxes.every(function(b) { return !(a[0] < b[2] && b[0] < a[2] && a[1] < b[3] && b[1] < a[3]); });
      if (out[i].shown) boxes.push(a);
    });
    return out;
  }

  /* ── Light cones in three dimensions ──

     A figure of _tools/derivations/projections.py that carries `turn` is drawn from its pieces
     in the drawing's (X, Y, T), T up, by the rules the generator's Figure projects them by:
     every line first, in the order published, projected and thinned by Ramer-Douglas-Peucker to
     0.0005 of the page's units; then every cone, farthest first by the depth of its apex, as the
     convex hull of its projected apex and rim, filled, its rim closed, every generator a `ribs`th
     of the way round the rim drawn from the apex, the two generators that bound the hull where
     the apex lies on it, and the apex. Cones whose apexes stand at depths closer than a hundred
     thousandth of the drawing keep the order they were published in, the order the generator
     painted them in at the figure's own camera. A slice of the embedding diagram on the floor
     is projected and thinned as a line is, and the page paints it under everything else.

     A label stands at the projection of its point, or where it names a circle about the axis
     at the point of that circle as far round from the camera's azimuth as it stood at the
     figure's own camera, so it keeps to the side of the circle nearest the reader, and gives
     way where it would overlap a label before it, as a label naming an embedding's circle does.

     The figure turns about `centre`, on the axis halfway up the drawing, which stays where it
     stood on the page. It is drawn at one scale, never above 1, the largest that keeps it within
     the width it had at the start on each side of the axis and within the height it had at the
     start, and moved up or down only as far as that height needs, as an embedding's height over
     a plane is; so at the figure's own camera it is the published figure. */
  function prepareCones(figure) {
    var box = figure.box, turn = figure.turn, c = turn.centre;
    var start = camera(figure.camera.azimuth, figure.camera.elevation);
    var M = { kind: 'cones', figure: figure, turn: turn, box: box, centre: c,
              size: Math.max(box[1] - box[0], box[3] - box[2]), unit: (box[1] - box[0]) / 560,
              origin: [c[0] * start.right[0] + c[1] * start.right[1] + c[2] * start.right[2],
                       c[0] * start.up[0] + c[1] * start.up[1] + c[2] * start.up[2]] };
    M.reach0 = coneReach(M, start);
    return M;
  }

  // Where a label of a figure of light cones stands in the drawing, seen from the camera.
  function labelPoint(place, cam) {
    if (place.at) return place.at;
    var phi = (cam.azimuth + place.angle) * RAD, rho = place.circle[0];
    return [rho * Math.cos(phi), rho * Math.sin(phi), place.circle[1]];
  }

  // Every point a figure of light cones draws, seen from the camera.
  function eachPoint(M, cam, f) {
    var turn = M.turn;
    turn.lines.forEach(function(L) { L.points.forEach(f); });
    turn.cones.forEach(function(cone) { f(cone.apex); cone.rim.forEach(f); });
    turn.labels.forEach(function(place) { f(labelPoint(place, cam)); });
    turn.slices.forEach(function(mark) {
      mark.lines.forEach(function(P) { P.forEach(f); });
      mark.fills.forEach(function(rings) { rings.forEach(function(P) { P.forEach(f); }); });
    });
  }

  // How far the figure reaches left, right, down and up of its centre on the page at scale 1.
  function coneReach(M, cam) {
    var c = M.centre, R = cam.right, U = cam.up, out = [Infinity, -Infinity, Infinity, -Infinity];
    eachPoint(M, cam, function(P) {
      var x = P[0] - c[0], y = P[1] - c[1], z = P[2] - c[2];
      var u = x * R[0] + y * R[1] + z * R[2], v = x * U[0] + y * U[1] + z * U[2];
      if (u < out[0]) out[0] = u;
      if (u > out[1]) out[1] = u;
      if (v < out[2]) out[2] = v;
      if (v > out[3]) out[3] = v;
    });
    return out;
  }

  // The one scale, never above 1, and the move up or down, that keep the figure within the
  // width it had at the start on each side of its axis and within the height it had.
  function coneFitting(M, cam) {
    var r0 = M.reach0, r = coneReach(M, cam), scale = 1;
    if (r[3] - r[2] > r0[3] - r0[2]) scale = Math.min(scale, (r0[3] - r0[2]) / (r[3] - r[2]));
    if (r[0] < r0[0]) scale = Math.min(scale, r0[0] / r[0]);
    if (r[1] > r0[1]) scale = Math.min(scale, r0[1] / r[1]);
    return { scale: scale, shift: Math.min(Math.max(0, r0[2] - scale * r[2]), r0[3] - scale * r[3]) };
  }

  /* The convex hull of points of the page, counterclockwise, by Andrew's monotone chain, each
     point taken to twelve decimals, as the generator's hull() takes it. */
  function hull(points) {
    var P = points.map(function(p) { return [Math.round(p[0] * 1e12) / 1e12, Math.round(p[1] * 1e12) / 1e12]; });
    P.sort(function(a, b) { return a[0] - b[0] || a[1] - b[1]; });
    if (P.length < 3) return P;
    function cross(o, a, b) { return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0]); }
    var lower = [], upper = [], i;
    for (i = 0; i < P.length; i++) {
      while (lower.length >= 2 && cross(lower[lower.length - 2], lower[lower.length - 1], P[i]) <= 0) lower.pop();
      lower.push(P[i]);
    }
    for (i = P.length - 1; i >= 0; i--) {
      while (upper.length >= 2 && cross(upper[upper.length - 2], upper[upper.length - 1], P[i]) <= 0) upper.pop();
      upper.push(P[i]);
    }
    return lower.slice(0, -1).concat(upper.slice(0, -1));
  }

  // Where the apex a stands on the hull H, as numpy's isclose() finds it, or -1.
  function onHull(H, a) {
    var x = Math.round(a[0] * 1e12) / 1e12, y = Math.round(a[1] * 1e12) / 1e12;
    for (var k = 0; k < H.length; k++) {
      if (Math.abs(H[k][0] - x) <= 1e-8 + 1e-5 * Math.abs(x) && Math.abs(H[k][1] - y) <= 1e-8 + 1e-5 * Math.abs(y)) return k;
    }
    return -1;
  }

  /* A figure of light cones at the camera (azimuth, elevation), in the published figure's form:
     layers painted in order, each cone's apex naming the cone's place in `turn.cones`, the
     slices of the embedding diagram, which the page paints under them, and for each published
     label where it stands and whether it is shown. */
  function drawCones(M, azimuth, elevation, quick, sizes) {
    var cam = camera(azimuth, elevation), fit = coneFitting(M, cam), s = fit.scale, c = M.centre, turn = M.turn;
    var R = cam.right, U = cam.up, V = cam.toward, X0 = M.origin[0], Y0 = M.origin[1] + fit.shift;
    function page(P) {
      var x = P[0] - c[0], y = P[1] - c[1], z = P[2] - c[2];
      return [X0 + s * (x * R[0] + y * R[1] + z * R[2]), Y0 + s * (x * U[0] + y * U[1] + z * U[2])];
    }
    function line(P) { return thin(P.map(page), 5e-4); }
    var layers = turn.lines.map(function(L) { return { kind: 'line', 'class': L['class'], points: line(L.points) }; });
    var tie = 1e-5 * M.size;
    turn.cones.map(function(cone, i) {
      var a = cone.apex;
      return { cone: cone, i: i, depth: Math.round((a[0] * V[0] + a[1] * V[1] + a[2] * V[2]) / tie) };
    }).sort(function(p, q) { return p.depth - q.depth || p.i - q.i; }).forEach(function(o) {
      var cls = o.cone['class'], a = page(o.cone.apex), rim = o.cone.rim.map(page), n = rim.length, H = hull([a].concat(rim));
      layers.push({ kind: 'fill', 'class': cls, points: H });
      layers.push({ kind: 'line', 'class': cls + '-rim', points: rim.concat([rim[0]]) });
      for (var i = 0; i < n; i += n / turn.ribs) layers.push({ kind: 'line', 'class': cls + '-rib', points: [a, rim[i]] });
      var k = onHull(H, a);
      if (k >= 0) layers.push({ kind: 'line', 'class': cls, points: [H[(k + H.length - 1) % H.length], a, H[(k + 1) % H.length]] });
      layers.push({ kind: 'point', 'class': cls + '-apex', at: a, cone: o.i });
    });
    var slices = turn.slices.map(function(mark) {
      return { lines: mark.lines.map(line), fills: mark.fills.map(function(rings) { return rings.map(line); }) };
    });
    var labels = giveWay(M.figure.labels, turn.labels.map(function(place) { return { at: page(labelPoint(place, cam)), shown: true }; }),
                         M.unit, sizes, function(i) { return !!turn.labels[i].circle; });
    return { layers: layers, slices: slices, labels: labels, scale: s };
  }

  /* ── The hand ──
     How the reader turns every figure that turns: a drag across the drawing's whole width, of
     `width` pixels, turns it half a turn round its axis, the way the hand moves, as far as the
     reader likes, and a drag up or down as far tilts it as much, from looking straight down the
     axis, elevation 90, to looking straight up it, -90, and never past, so the axis always
     stands up the page. The arrow keys turn it by 15 degrees, as a drag their way would, and
     Home and Escape bring back the figure's own camera, `start`, as a double click and a double
     tap do. */
  function aimed(azimuth, elevation) {
    return { azimuth: azimuth, elevation: Math.max(-90, Math.min(90, elevation)) };
  }
  function dragged(from, dx, dy, width) {
    var k = 180 / width;
    return aimed(from.azimuth - k * dx, from.elevation + k * dy);
  }
  function keyed(at, key, start) {
    if (key === 'ArrowLeft') return aimed(at.azimuth + 15, at.elevation);
    if (key === 'ArrowRight') return aimed(at.azimuth - 15, at.elevation);
    if (key === 'ArrowUp') return aimed(at.azimuth, at.elevation - 15);
    if (key === 'ArrowDown') return aimed(at.azimuth, at.elevation + 15);
    if (key === 'Home' || key === 'Escape') return aimed(start.azimuth, start.elevation);
    return null;
  }
  // Whether the camera is anywhere but the figure's own, a whole number of turns apart being the same.
  function turned(start, at) {
    return at.elevation !== start.elevation || Math.abs((((at.azimuth - start.azimuth) % 360) + 540) % 360 - 180) > 1e-9;
  }

  // What draw() needs of a figure that does not depend on the camera: an embedding view, which
  // carries its surfaces, or a figure of light cones, which carries `turn`.
  function prepare(view) { return view.surfaces ? prepareSurfaces(view) : prepareCones(view); }

  /* The figure at the camera (azimuth, elevation), coarser where `quick` while a drag goes on,
     with the labels measured by `sizes`: its layers, the slices of the embedding diagram a
     figure of light cones carries, its labels and the scale it is drawn at. */
  function draw(M, azimuth, elevation, quick, sizes) {
    return (M.kind === 'cones' ? drawCones : drawSurfaces)(M, azimuth, elevation, quick, sizes);
  }

  var api = { camera: camera, prepare: prepare, draw: draw, frame: frame, dragged: dragged, keyed: keyed, turned: turned,
              fitting: fitting, hidden: hidden, isolines: isolines, labelBox: labelBox, hull: hull };
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  else root.MfsTurn = api;
})(this);
