#!/usr/bin/env node
/* Walk every spacetime with more than one chart on the spacetimes page in headless Chrome, as
   a reader does, and hold the Conventions section to the chart shown in every state.

     bundle exec jekyll serve
     node _tools/chart_conventions.mjs http://127.0.0.1:4000
     node _tools/chart_conventions.mjs http://127.0.0.1:4000 --phone

   Each spacetime is opened from the list, and its charts are chosen along a walk that goes
   from every chart to every other, once by pointer, a mouse on a desktop and a finger on a
   phone, and once from the keyboard, Enter and Space in turn on the chart's button. At every
   chart every other row of choices, the views of each diagram and the index placements, is
   pressed. Then a chart is chosen and another spacetime opened and this one again, a chart is
   chosen and the page reloaded and this spacetime opened again, and two charts are chosen one
   straight after the other, before the first is drawn.

   A spacetime keeps its chart while the page is open, so it opens again in the chart chosen
   in it last, and a reloaded page opens every spacetime in its first chart. Every file of the
   collection the page fetched has to carry its build's stamp in its address.

   After each step the page is left to settle, and the Conventions paragraph, the pressed chart
   button, the chart's domains, the line element, the metric, every tensor's components and
   the views of the spacetime diagram are each read back, every formula as the TeX MathJax set
   it from, and held to one chart of the metric file. The paragraph has to read exactly that
   chart's `convention` followed by the spacetime's shared `convention`, every other part has
   to be that same chart, and the chart has to be the one the reader chose, or the spacetime's
   chart last chosen where it was opened again. Any other state counts as a mismatch, named
   with the spacetime, the control, what was shown and what should have been.

   --metric <id> walks one spacetime and is repeatable; --phone lays the page out as an iPhone
   held upright, 390 by 844, and taps with a finger. It prints the pairs of spacetime and chart
   it held, then every mismatch and page error, and exits non-zero if there is any. */
import { launch, sleep } from './chrome.mjs';

const argv = process.argv.slice(2);
const base = (argv[0] && !argv[0].startsWith('--') ? argv[0] : 'http://127.0.0.1:4000').replace(/\/$/, '');
const only = argv.flatMap((a, i) => (a === '--metric' ? [argv[i + 1]] : []));
const phone = argv.includes('--phone');
const width = phone ? 390 : 1440, height = phone ? 844 : 900;
const LIMIT_MS = 120000;

const { send, evaluate, finish, onError } = await launch({ width, height, phone });
const mismatches = [], errors = [];
let opened = 'the list';
onError(message => errors.push(`${opened}: ${message}`));

await send('Page.addScriptToEvaluateOnNewDocument', { source: `
  window.__cc = { busy: 0, errors: [] };
  window.addEventListener('error', function (e) { window.__cc.errors.push(e.message); });
  window.addEventListener('unhandledrejection', function (e) { window.__cc.errors.push('unhandled rejection: ' + (e.reason && e.reason.stack || e.reason)); });
  (function wrap() {
    if (!(window.MathJax && MathJax.typesetPromise && MathJax.startup && MathJax.startup.document)) { setTimeout(wrap, 20); return; }
    var typeset = MathJax.typesetPromise;
    MathJax.typesetPromise = function () {
      window.__cc.busy++;
      function done() { window.__cc.busy--; }
      return typeset.apply(this, arguments).then(function (v) { done(); return v; }, function (e) { done(); throw e; });
    };
    window.__cc.ready = true;
  })();
  /* Resolves once the panel stands in place with a header other than the one given, all its
     mathematics set, no typesetting outstanding and its fonts in, and two frames have passed. */
  window.__cc.settle = function (before, limit) {
    var panel = document.getElementById('mfs-content-panel'), start = performance.now();
    return new Promise(function (resolve) {
      (function poll() {
        var header = panel.querySelector('.mfs-header');
        var ready = header && header !== before && window.__cc.busy === 0 && !panel.classList.contains('mfs-fonts-wait') &&
          /^translateX\\(0(px)?\\)$/.test(panel.style.transform) && !/\\\\[(\\[]/.test(panel.textContent);
        if (!ready) {
          if (performance.now() - start > limit) { resolve(false); return; }
          setTimeout(poll, 25); return;
        }
        requestAnimationFrame(function () { requestAnimationFrame(function () { setTimeout(function () { resolve(true); }, 0); }); });
      })();
    });
  };
  /* Every part of the panel that depends on the chart, each formula as the TeX it was set from. */
  window.__cc.read = function () {
    var panel = document.getElementById('mfs-content-panel');
    var roots = new Map();
    MathJax.startup.document.getMathItemsWithin(panel).forEach(function (m) { roots.set(m.typesetRoot, m.math); });
    function read(node) {
      if (roots.has(node)) return '$' + roots.get(node) + '$';
      if (node.nodeType === 3) return node.data;
      return [].map.call(node.childNodes, read).join('');
    }
    function texOf(el) {
      var out = [];
      (function walk(node) {
        if (roots.has(node)) { out.push(roots.get(node)); return; }
        [].forEach.call(node.childNodes, walk);
      })(el);
      return out;
    }
    function section(label) {
      return [].find.call(panel.querySelectorAll('.mfs-section'), function (s) {
        var l = s.querySelector(':scope > .mfs-section-label, :scope > .mfs-section-head > .mfs-section-label');
        return l && l.textContent.trim().toLowerCase() === label;
      });
    }
    var box = panel.querySelector('.mfs-convention');
    var paragraphs = box ? [].filter.call(box.querySelectorAll('p'), function (p) { return !p.querySelector('.mfs-convention-k'); }) : [];
    var pressed = [].filter.call(panel.querySelectorAll('.mfs-charts .mfs-choice'), function (b) {
      return b.getAttribute('aria-pressed') === 'true' || b.classList.contains('mfs-choice-on');
    }).map(function (b) { return +b.dataset.chart; });
    var line = section('line element'), metric = section('metric'), diagram = panel.querySelector('.mfs-diagram:not(.mfs-conformal):not(.mfs-embedding)');
    var tensors = {};
    ['christoffel', 'riemann', 'ricci', 'einstein', 'weyl'].forEach(function (t) {
      var comp = document.getElementById('mfs-comp-' + t);
      if (!comp) return;
      var on = comp.parentNode.querySelector('.mfs-choice-on[data-placement]');
      tensors[t] = { placement: on ? on.dataset.placement : null,
        lines: [].map.call(comp.querySelectorAll('.mfs-line-body'), function (b) { return texOf(b).join(''); }) };
    });
    return {
      name: (panel.querySelector('.mfs-title') || {}).textContent || '',
      conventions: paragraphs.map(read).join(' ').replace(/\\u00AD/g, '').replace(/\\s+/g, ' ').trim(),
      paragraphs: paragraphs.length,
      pressed: pressed,
      domains: [].map.call(panel.querySelectorAll('.mfs-domain-body'), function (d) { return texOf(d).join(''); }),
      line: line ? texOf(line).join('') : null,
      metric: metric ? texOf(metric).join('') : null,
      tensors: tensors,
      views: diagram ? diagram.querySelectorAll('.mfs-nr-view').length : 0,
      viewLabels: diagram ? [].map.call(diagram.querySelectorAll('.mfs-section-head .mfs-choice'), function (b) {
        return read(b).replace(/\\s+/g, ' ').trim();
      }) : []
    };
  };
` });

await send('Page.navigate', { url: `${base}/MFS/` });
async function listed() {
  for (let i = 0; i < 600; i++) {
    if (await evaluate(`!!(window.__cc && window.__cc.ready && document.querySelector('.mfs-result'))`)) return;
    await sleep(100);
  }
  console.error(`No spacetimes listed at ${base}/MFS/.`);
  finish(2);
}
await listed();

const index = await (await fetch(`${base}/MFS/assets/data/metrics_index.json`)).json();
const allIds = (Array.isArray(index) ? index : index.metrics || []).map(e => (typeof e === 'string' ? e : e.id));
const files = {}, diagrams = {};
for (const id of allIds) {
  files[id] = await (await fetch(`${base}/MFS/assets/data/metrics/${id}.json`)).json();
}
const ids = allIds.filter(id => (files[id].coordinates || []).length > 1 && (!only.length || only.includes(id)));
for (const id of only) if (!ids.includes(id)) { console.error(`No spacetime ${id} with more than one chart.`); finish(2); }
for (const id of allIds) {
  const r = await fetch(`${base}/MFS/assets/data/diagrams/${id}.json`);
  diagrams[id] = r.ok ? await r.json() : null;
}

/* What each part of the panel reads in chart c of spacetime `data`, built as the page builds it. */
function trimmed(t) { return t && t.trim() ? t.trim() : ''; }
function conventionOf(data, c) {
  return [data.coordinates[c].convention, data.convention].map(trimmed).filter(Boolean).join(' ');
}
function metricOf(coord) {
  if (!coord.metric_components) return null;
  const n = coord.coords.length, grid = Array.from({ length: n }, () => Array(n).fill('0'));
  for (const m of coord.metric_components) {
    const i = coord.coords.indexOf(m.indices[0]), j = coord.coords.indexOf(m.indices[1]);
    if (i >= 0 && j >= 0) grid[i][j] = m.value;
  }
  return 'g_{\\mu\\nu} = \\begin{bmatrix}' + grid.map(r => r.join(' & ')).join(' \\\\ ') + '\\end{bmatrix}';
}
const TENSOR_FIELD = { christoffel: 'christoffel', riemann: 'riemann', ricci: 'ricci_tensor', einstein: 'einstein_tensor', weyl: 'weyl_tensor' };
function viewsOf(id, coord) {
  const d = diagrams[id];
  if (!d) return [];
  return ((d.systems || {})[coord.id] || []).concat(((d.projections || {})[coord.id]) || []);
}
/* The charts of `data` that each part read on the page could be. */
function chartsOf(id, data, shown) {
  const all = data.coordinates.map((_, c) => c);
  const parts = {
    conventions: all.filter(c => shown.conventions === conventionOf(data, c)),
    button: shown.pressed.length === 1 ? [shown.pressed[0]] : [],
    domains: all.filter(c => JSON.stringify(shown.domains) === JSON.stringify(data.coordinates[c].domains || [])),
    'line element': all.filter(c => shown.line === (data.coordinates[c].line_element || null)),
    metric: all.filter(c => shown.metric === metricOf(data.coordinates[c])),
    components: all.filter(c => {
      const coord = data.coordinates[c];
      return Object.keys(TENSOR_FIELD).every(t => {
        const field = coord[TENSOR_FIELD[t]], got = shown.tensors[t];
        if (!field || !got) return !field === !got;
        const placement = got.placement || field.default, variant = field.variants[placement];
        if (!variant) return false;
        const values = new Set(variant.nonzero.map(v => v.value));
        if (!variant.nonzero.length) return got.lines.length === 0;
        return got.lines.length > 0 && got.lines.every(l => [...values].some(v => l.endsWith(' = ' + v) || l === '{} = ' + v));
      });
    }),
    'spacetime diagram': all.filter(c => {
      const views = viewsOf(id, data.coordinates[c]);
      if (views.length !== shown.views) return false;
      return views.length < 2 || JSON.stringify(views.map(v => v.label.replace(/\s+/g, ' ').trim())) === JSON.stringify(shown.viewLabels);
    })
  };
  return parts;
}
const held = new Set();
let steps = 0;
function name(data, c) { return `${data.coordinates[c].name} (${data.coordinates[c].id})`; }
/* Hold the page to chart `want` of spacetime `id`, reached by `how`. */
async function check(id, want, how) {
  const data = files[id], shown = await evaluate('window.__cc.read()');
  steps++;
  if (ids.includes(id)) held.add(`${id} ${want}`);
  const where = `${id} / ${how}`;
  if (shown.name !== data.name) { mismatches.push(`${where}: the panel shows ${JSON.stringify(shown.name)}, not ${JSON.stringify(data.name)}`); return; }
  const parts = chartsOf(id, data, shown);
  const wrong = Object.keys(parts).filter(p => !parts[p].includes(want));
  if (shown.paragraphs !== 1) wrong.push(`conventions in ${shown.paragraphs} paragraphs`);
  if (!wrong.length) return;
  const convAs = parts.conventions.length ? parts.conventions.map(c => name(data, c)).join(' or ') : `no chart: ${JSON.stringify(shown.conventions)}`;
  const lineAs = parts['line element'].length ? parts['line element'].map(c => name(data, c)).join(' or ') : 'no chart';
  mismatches.push(`${where}: chose ${name(data, want)}; conventions read as ${convAs}, line element as ${lineAs}, ` +
    `button ${shown.pressed.map(c => name(data, c)).join(', ') || 'none'}; wrong: ${wrong.join(', ')}`);
}

/* Mark the panel as it stands before an action, and wait after it for the panel to be drawn
   again (or, with `again` false, only for it to settle). */
async function mark() {
  await evaluate(`window.__cc.before = document.querySelector('#mfs-content-panel .mfs-header')`);
}
async function settle(again, limit = LIMIT_MS) {
  return evaluate(`window.__cc.settle(${again ? 'window.__cc.before' : 'null'}, ${limit})`);
}
/* The centre of the element `selector` picks, scrolled into the middle of the view. */
async function centre(selector) {
  return evaluate(`(function () {
    var el = document.querySelector(${JSON.stringify(selector)});
    if (!el) return null;
    el.scrollIntoView({ block: 'center', inline: 'center' });
    var r = el.getBoundingClientRect();
    var x = r.left + r.width / 2, y = r.top + r.height / 2, hit = document.elementFromPoint(x, y);
    return { x: x, y: y, hit: !!hit && (hit === el || el.contains(hit)) };
  })()`);
}
/* Press the element `selector` picks where it stands still, since the list slides in as the
   page opens, as a reader's finger or mouse would. */
async function point(selector) {
  let at = await centre(selector);
  for (let i = 0; i < 100 && at; i++) {
    await sleep(100);
    const now = await centre(selector);
    const still = now && now.hit && Math.abs(now.x - at.x) < 0.5 && Math.abs(now.y - at.y) < 0.5;
    at = now;
    if (still) break;
  }
  if (!at || !at.hit) throw new Error(`nothing to press at ${selector}`);
  if (phone) {
    await send('Input.dispatchTouchEvent', { type: 'touchStart', touchPoints: [{ x: at.x, y: at.y }] });
    await send('Input.dispatchTouchEvent', { type: 'touchEnd', touchPoints: [] });
  } else {
    await send('Input.dispatchMouseEvent', { type: 'mouseMoved', x: at.x, y: at.y });
    await send('Input.dispatchMouseEvent', { type: 'mousePressed', x: at.x, y: at.y, button: 'left', clickCount: 1 });
    await send('Input.dispatchMouseEvent', { type: 'mouseReleased', x: at.x, y: at.y, button: 'left', clickCount: 1 });
  }
}
let presses = 0;
async function key(selector) {
  await evaluate(`document.querySelector(${JSON.stringify(selector)}).focus()`);
  const enter = presses++ % 2 === 0;
  const k = enter ? { key: 'Enter', code: 'Enter', windowsVirtualKeyCode: 13, text: '\r' } : { key: ' ', code: 'Space', windowsVirtualKeyCode: 32, text: ' ' };
  await send('Input.dispatchKeyEvent', { type: 'keyDown', ...k });
  await send('Input.dispatchKeyEvent', { type: 'keyUp', key: k.key, code: k.code, windowsVirtualKeyCode: k.windowsVirtualKeyCode });
  return enter ? 'Enter' : 'Space';
}
const chartButton = c => `#mfs-content-panel .mfs-charts [data-chart="${c}"]`;

async function open(id) {
  opened = id;
  await mark();
  await point(`.mfs-result[data-id="${id}"]`);
  if (!await settle(true)) mismatches.push(`${id}: never drawn once opened from the list`);
}
/* Choose chart c with the pointer or the keyboard, and wait for it to be drawn. */
async function choose(id, c, by) {
  let how = by === 'keyboard' ? '' : phone ? 'tap' : 'click';
  await mark();
  if (by === 'keyboard') how = await key(chartButton(c)); else await point(chartButton(c));
  if (!await settle(true, 20000)) mismatches.push(`${id} / ${how} on ${files[id].coordinates[c].id}: the page was not drawn again`);
  return how;
}
/* A walk through the complete directed graph on n charts from chart 0 that takes every step
   from one chart to another once, Hierholzer's circuit. */
function walk(n) {
  const out = Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => j).filter(j => j !== i).reverse());
  const stack = [0], path = [];
  while (stack.length) {
    const v = stack[stack.length - 1];
    if (out[v].length) stack.push(out[v].pop()); else path.push(stack.pop());
  }
  return path.reverse();
}
/* Press every other row of choices once, the views of each diagram and the placements. */
async function others(id, c) {
  const rows = await evaluate(`[].map.call(document.querySelectorAll('#mfs-content-panel .mfs-choices:not(.mfs-charts)'), function (row, r) {
    row.setAttribute('data-cc-row', r);
    return row.querySelectorAll('.mfs-choice').length;
  })`);
  for (let r = 0; r < rows.length; r++) {
    for (let b = rows[r] - 1; b >= 0; b--) {
      const sel = `#mfs-content-panel [data-cc-row="${r}"] .mfs-choice:nth-child(${b + 1})`;
      if (!await evaluate(`!!document.querySelector(${JSON.stringify(sel)})`)) continue;
      await point(sel);
      await settle(false, 20000);
      await check(id, c, `${phone ? 'tap' : 'click'} on ${await evaluate(`document.querySelector(${JSON.stringify(sel)}).closest('.mfs-section').querySelector('.mfs-section-label').textContent.trim()`)} choice ${b + 1} in ${files[id].coordinates[c].id}`);
    }
  }
}

/* Every file of the collection the page fetched carries the stamp of its build: a metric,
   diagram, conformal or embedding file its index entry's stamp, the index and the
   bibliography the time the site was built. */
async function stamps(id) {
  const index = await evaluate('window._mfsIndex');
  const urls = await evaluate(`performance.getEntriesByType('resource').map(function (e) { return e.name; })
    .filter(function (u) { return u.indexOf('/assets/data/') >= 0; })`);
  const field = { metrics: 'version', diagrams: 'diagrams', conformal: 'conformal', embedding: 'embedding' };
  for (const u of urls) {
    const url = new URL(u), v = url.searchParams.get('v'), m = url.pathname.match(/\/MFS\/assets\/data\/(metrics|diagrams|conformal|embedding)\/([^/]+)\.json$/);
    if (!v) { errors.push(`${id}: fetched ${url.pathname} with no stamp`); continue; }
    if (!m) continue;
    const entry = index.find(e => e.id === m[2]);
    if (!entry || entry[field[m[1]]] !== v) errors.push(`${id}: fetched ${url.pathname} at ${v}, not its stamp ${entry && entry[field[m[1]]]}`);
  }
}

const expected = {};   // the chart of each spacetime the reader chose last
for (const id of ids) {
  const data = files[id], n = data.coordinates.length, order = walk(n);
  const other = ids[(ids.indexOf(id) + 1) % ids.length] === id ? allIds.find(x => x !== id) : ids[(ids.indexOf(id) + 1) % ids.length];

  await open(id);
  let want = expected[id] ?? 0;
  await check(id, want, 'opened from the list');
  for (const by of ['pointer', 'keyboard']) {
    if (want !== 0) { await choose(id, 0, by); want = 0; expected[id] = 0; }
    for (const c of order.slice(1)) {
      const how = await choose(id, c, by);
      want = c; expected[id] = c;
      await check(id, c, `${how} on ${data.coordinates[c].id}`);
      if (by === 'pointer') await others(id, c);
    }
  }

  // Away to another spacetime and back, with the last chart of the walk chosen.
  const last = n - 1;
  if (want !== last) { await choose(id, last, 'pointer'); want = last; expected[id] = last; }
  await open(other);
  await check(other, expected[other] ?? 0, `opened after ${id}`);
  await open(id);
  await check(id, expected[id], `opened again after ${other}, ${data.coordinates[last].id} chosen before`);

  // Two charts chosen one straight after the other, the first not yet drawn.
  const a = (expected[id] + 1) % n, b = (expected[id] + 2) % n;
  await mark();
  await evaluate(`document.querySelector(${JSON.stringify(chartButton(a))}).click();
    document.querySelector(${JSON.stringify(chartButton(b))}) && document.querySelector(${JSON.stringify(chartButton(b))}).click();`);
  await settle(true, 20000);
  expected[id] = b;
  await check(id, b, `click on ${data.coordinates[a].id} and at once ${data.coordinates[b].id}`);

  // Reloaded with a chart chosen, and opened again: a new page, so the first chart.
  await stamps(id);
  opened = `${id} reloaded`;
  await send('Page.reload', {});
  await sleep(300);
  await listed();
  await open(id);
  const before = expected[id];
  for (const k of Object.keys(expected)) delete expected[k];
  await check(id, 0, `reloaded with ${data.coordinates[before].id} chosen, then opened`);
  for (const e of await evaluate(`window.__cc.errors.splice(0)`)) errors.push(`${id}: page error: ${e}`);
}

const pairs = ids.reduce((s, id) => s + files[id].coordinates.length, 0);
for (const m of mismatches) console.log(`MISMATCH ${m}`);
for (const e of errors) console.log(`ERROR ${e}`);
console.log(`${ids.length} spacetimes, ${held.size} of ${pairs} spacetime-chart pairs held in ${steps} states at ${width} by ${height}` +
  `${phone ? ' as a phone' : ''}: ${mismatches.length} mismatches, ${errors.length} errors`);
finish(mismatches.length || errors.length || held.size !== pairs ? 1 : 0);
