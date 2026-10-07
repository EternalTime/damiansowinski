#!/usr/bin/env node
/* Open the spacetimes page in headless Chrome as a reader does and hold the graph of relations
   behind the list to what _tools/README.md, "The graph behind the list", says it does:

     bundle exec jekyll serve
     node _tools/background_graph.mjs http://127.0.0.1:4000

   On a desktop the graph is drawn in the room to the right of the list and covers no panel,
   every spacetime is in it, and no more than forty names are written. Typing a keyword into
   the search, a key at a time, gathers exactly the spacetimes the list then shows, which are the
   ones carrying it, with the rest behind; clearing the search brings the whole graph back.
   The list and the graph read as one thing: the pointer or the keyboard on a name in the list
   lights its spacetime in the graph, that one alone and with its name written, and the pointer
   on a spacetime in the graph turns its name in the list the page's pink, with no background,
   and scrolls the list to it, and both hold under a search, where a spacetime the search left out
   is never lit and never opened. A press on a spacetime in the graph opens it, and the graph
   stays where it is under the spacetime's panel, which covers it, standing still.
   The search field keeps its own colours in every state, empty, focused, typed in, with names
   found and with none, with keywords offered and one gone to, chosen and cleared, and the browser
   offers nothing of its own in it. From the first letter typed the keywords the search can be
   finished with are offered under it, one letter and two offering every keyword with a word they
   start, set as the names in the list are; the arrow keys and Enter choose one, and so does the
   pointer, and the list and the graph then show exactly the spacetimes carrying it. With --phone
   the page is laid out as an iPhone held upright, the graph is not drawn and its relations are
   never fetched, and the field and its keywords are held to the same. Every console error and
   page error is an error.
   It exits non-zero on any of them. _tools/chrome.mjs starts Chrome. */
import { launch, sleep } from './chrome.mjs';

const argv = process.argv.slice(2);
const base = (argv[0] && !argv[0].startsWith('--') ? argv[0] : 'http://127.0.0.1:4000').replace(/\/$/, '');
const phone = argv.includes('--phone');
const keywords = ['conformally flat', 'wormhole', 'vacuum', 'Petrov type D'];

const { send, evaluate, finish, onError } = await launch({ width: phone ? 390 : 1440, height: phone ? 844 : 900, phone });
const errors = [], failures = [];
onError(message => errors.push(message));
await send('Page.addScriptToEvaluateOnNewDocument', { source: `
  window.__mfsErrors = [];
  window.addEventListener('error', function (e) { window.__mfsErrors.push(String(e.message)); });
  window.addEventListener('unhandledrejection', function (e) { window.__mfsErrors.push(String(e.reason)); });
` });

function check(ok, what) {
  console.log(`${ok ? 'ok  ' : 'FAIL'} ${what}`);
  if (!ok) failures.push(what);
}
async function until(expression, what, seconds = 20) {
  for (let t = 0; t < seconds * 20; t++) {
    if (await evaluate(expression)) return true;
    await sleep(50);
  }
  check(false, `${what} within ${seconds}s`);
  return false;
}
const shows = () => evaluate('window._mfsGraph.shows()');
const listed = () => evaluate("[].map.call(document.querySelectorAll('#mfs-search-results .mfs-result'), function (r) { return r.dataset.id; })");
const index = () => evaluate('window._mfsIndex');
// Types a key at a time and then puts away the keywords offered, which lie over the top of the list.
async function type(text) {
  await evaluate("document.getElementById('mfs-search-input').focus()");
  for (const key of text) await send('Input.insertText', { text: key });
  await key('Escape', 'Escape', 27);
}
async function clear() {
  await evaluate("document.getElementById('mfs-search-input').select()");
  await send('Input.dispatchKeyEvent', { type: 'keyDown', key: 'Backspace', code: 'Backspace', windowsVirtualKeyCode: 8 });
  await send('Input.dispatchKeyEvent', { type: 'keyUp', key: 'Backspace', code: 'Backspace', windowsVirtualKeyCode: 8 });
}
async function point(x, y) {
  await send('Input.dispatchMouseEvent', { type: 'mouseMoved', x, y, button: 'none' });
}
async function key(name, code, keyCode) {
  await send('Input.dispatchKeyEvent', { type: 'keyDown', key: name, code, windowsVirtualKeyCode: keyCode });
  await send('Input.dispatchKeyEvent', { type: 'keyUp', key: name, code, windowsVirtualKeyCode: keyCode });
}
// Somewhere on the page that is neither the list nor the graph: the grid above the panels.
const away = () => point(700, 8);
// Where a name in the list is, whether all of it is in sight in the list, and how it is drawn.
const nameIn = id => evaluate(`(function () {
  var b = document.querySelector('#mfs-search-results .mfs-result[data-id="${id}"]'), list = document.getElementById('mfs-search-results');
  if (!b) return null;
  var r = b.getBoundingClientRect(), l = list.getBoundingClientRect(), cs = getComputedStyle(b);
  return { x: r.left + r.width / 2, y: r.top + r.height / 2, inSight: r.top >= l.top - 0.5 && r.bottom <= l.bottom + 0.5,
           background: cs.backgroundColor, color: cs.color };
})()`);
const litNames = () => evaluate("[].map.call(document.querySelectorAll('#mfs-search-results .mfs-result-lit'), function (r) { return r.dataset.id; })");
// The spacetimes whose points are drawn in the page's pink, read from the canvas at the middle
// of each point in front.
const pinkPoints = () => evaluate(`(function () {
  var c = document.getElementById('mfs-graph'), box = c.getBoundingClientRect(), k = c.width / box.width, pen = c.getContext('2d');
  var hex = getComputedStyle(document.documentElement).getPropertyValue('--pink-dark').trim();
  var pink = [1, 3, 5].map(function (i) { return parseInt(hex.slice(i, i + 2), 16); });
  return window._mfsGraph.shows().front.filter(function (s) {
    var d = pen.getImageData(Math.floor((s.x - box.left) * k), Math.floor((s.y - box.top) * k), 1, 1).data;
    return Math.abs(d[0] - pink[0]) + Math.abs(d[1] - pink[1]) + Math.abs(d[2] - pink[2]) <= 9 && d[3] === 255;
  }).map(function (s) { return s.id; });
})()`);
// The graph has drawn what was last asked of it.
const drawn = "!window._mfsGraph.shows().drawing";

/* The pointer and then the keyboard on a name in the list light its spacetime in the graph, that
   one alone, with its name written, and nothing is lit once they leave; the pointer on a
   spacetime in the graph turns its name in the list the page's pink, with nothing else about it
   changed, and brings it into sight, and nothing is lit once it leaves. `ids` are the spacetimes the list shows. */
async function asOneThing(ids, what) {
  await away();
  await evaluate("document.getElementById('mfs-search-results').scrollTop = 0");
  await until(`window._mfsGraph.shows().lit === null && ${drawn}`, `${what}: nothing is lit at first`);
  check((await pinkPoints()).length === 0 && (await litNames()).length === 0, `${what}: no point and no name is lit at first`);

  const first = ids[0], at = await nameIn(first);
  await point(at.x, at.y);
  await until(`window._mfsGraph.shows().lit === ${JSON.stringify(first)} && ${drawn}`, `${what}: the pointer on ${first} in the list lights it in the graph`);
  let state = await shows();
  check(state.litNamed, `${what}: ${first} lit from the list has its name written`);
  check(JSON.stringify(await pinkPoints()) === JSON.stringify([first]), `${what}: ${first} alone is drawn lit (${await pinkPoints()})`);
  await away();
  await until(`window._mfsGraph.shows().lit === null && ${drawn}`, `${what}: leaving ${first} in the list lets it go in the graph`);
  check((await pinkPoints()).length === 0, `${what}: no point is lit once the pointer leaves the list`);

  // The keyboard: Tab from the search field reaches the names in their order.
  await evaluate("document.getElementById('mfs-search-input').focus()");
  for (let k = 0; k < Math.min(2, ids.length); k++) {
    await key('Tab', 'Tab', 9);
    const focused = await evaluate("document.activeElement.dataset ? document.activeElement.dataset.id || null : null");
    check(focused === ids[k], `${what}: Tab ${k + 1} reaches ${ids[k]} in the list (${focused})`);
    await until(`window._mfsGraph.shows().lit === ${JSON.stringify(ids[k])} && ${drawn}`, `${what}: the keyboard on ${ids[k]} lights it in the graph`);
    check(JSON.stringify(await pinkPoints()) === JSON.stringify([ids[k]]), `${what}: the keyboard lights ${ids[k]} alone`);
  }
  await evaluate('document.activeElement.blur()');
  await until(`window._mfsGraph.shows().lit === null && ${drawn}`, `${what}: the keyboard leaving the list lets go in the graph`);

  // The graph: a named spacetime in front whose name in the list is out of sight where it can be.
  await evaluate("document.getElementById('mfs-search-results').scrollTop = 0");
  state = await shows();
  const named = state.front.filter(s => s.named && ids.includes(s.id));
  let target = null;
  for (const s of named) if (!(await nameIn(s.id)).inSight) { target = s; break; }
  target = target || named[named.length - 1];
  const bare = await nameIn(target.id);
  await point(target.x, target.y);
  await until(`window._mfsGraph.shows().lit === ${JSON.stringify(target.id)} && ${drawn}`, `${what}: the pointer on ${target.id} in the graph lights it`);
  await until(`(function () { var b = document.querySelector('#mfs-search-results .mfs-result[data-id="${target.id}"]'), l = document.getElementById('mfs-search-results');
    var r = b.getBoundingClientRect(), m = l.getBoundingClientRect(); return r.top >= m.top - 0.5 && r.bottom <= m.bottom + 0.5; })()`,
    `${what}: the list scrolls ${target.id} into sight`);
  check(JSON.stringify(await litNames()) === JSON.stringify([target.id]), `${what}: ${target.id} alone is lit in the list (${await litNames()})`);
  const pink = await evaluate(`(function () {
    var hex = getComputedStyle(document.documentElement).getPropertyValue('--pink-light').trim();
    return 'rgb(' + [1, 3, 5].map(function (i) { return parseInt(hex.slice(i, i + 2), 16); }).join(', ') + ')';
  })()`);
  const lit = await nameIn(target.id);
  check(lit.color === pink, `${what}: lit from the graph, ${target.id}'s name turns the page's pink (${lit.color}, pink ${pink})`);
  check(lit.background === bare.background, `${what}: lit from the graph, ${target.id}'s button keeps its background (${lit.background}, ${bare.background} unlit)`);
  check(JSON.stringify(await pinkPoints()) === JSON.stringify([target.id]), `${what}: ${target.id} alone is drawn lit in the graph`);
  await away();
  await until(`window._mfsGraph.shows().lit === null && ${drawn}`, `${what}: leaving ${target.id} in the graph lets it go`);
  check((await litNames()).length === 0, `${what}: no name is lit in the list once the pointer leaves the graph`);
  return target;
}

/* The search field's colours, every one it is drawn in, which are the same in every state as at
   rest; and the keywords, which the page offers for what is typed. */
const look = () => evaluate(`(function () {
  var cs = getComputedStyle(document.getElementById('mfs-search-input')), out = {};
  ['backgroundColor', 'backgroundImage', 'color', 'webkitTextFillColor', 'caretColor', 'borderTopColor', 'borderRightColor',
   'borderBottomColor', 'borderLeftColor', 'borderTopWidth', 'borderBottomWidth', 'borderTopStyle', 'outlineStyle',
   'outlineColor', 'outlineWidth', 'boxShadow', 'opacity', 'filter'].forEach(function (k) { out[k] = cs[k]; });
  return out;
})()`);
const offered = () => evaluate(`(function () {
  var k = document.getElementById('mfs-keywords');
  return k.hidden ? null : [].map.call(k.children, function (o) { return o.textContent; });
})()`);
const goneTo = () => evaluate("[].map.call(document.querySelectorAll('#mfs-keywords .mfs-keyword-on'), function (o) { return o.textContent; })");
function keywordsFor(idx, typed) {
  const q = typed.trim().toLowerCase(), tags = [...new Set(idx.flatMap(m => m.tags))].filter(t => wordMatches(t, q));
  const order = (a, b) => a.toLowerCase() < b.toLowerCase() ? -1 : a.toLowerCase() > b.toLowerCase() ? 1 : a < b ? -1 : a > b ? 1 : 0;
  return q ? [...tags.filter(t => t.toLowerCase().startsWith(q)).sort(order), ...tags.filter(t => !t.toLowerCase().startsWith(q)).sort(order)] : [];
}
async function clearField() {
  await evaluate("document.getElementById('mfs-search-input').select()");
  await key('Backspace', 'Backspace', 8);
}

/* The search field keeps its own look in every state, and the keywords are offered from the first
   letter, chosen by the keyboard and by the pointer, and select exactly the spacetimes carrying
   them, in the list and, where it is drawn, in the graph. */
async function searchField() {
  const idx = await index(), field = "document.getElementById('mfs-search-input')";
  const attributes = await evaluate(`['autocomplete', 'autocorrect', 'autocapitalize', 'spellcheck', 'role', 'aria-expanded'].map(function (a) { return ${field}.getAttribute(a); })`);
  check(JSON.stringify(attributes) === JSON.stringify(['off', 'off', 'off', 'false', 'combobox', 'false']),
        `the browser offers nothing of its own in the search field (${attributes})`);
  await evaluate(`${field}.blur()`);
  const rest = await look();
  const own = await evaluate(`(function () {
    var root = getComputedStyle(document.documentElement);
    function rgb(name) { var hex = root.getPropertyValue(name).trim(); return 'rgb(' + [1, 3, 5].map(function (i) { return parseInt(hex.slice(i, i + 2), 16); }).join(', ') + ')'; }
    return { cyan: rgb('--cyan'), text: rgb('--text-bright') };
  })()`);
  check(rest.backgroundColor === 'rgba(255, 255, 255, 0.05)' && rest.color === own.text && rest.borderTopColor === own.cyan &&
        rest.borderBottomColor === own.cyan && rest.outlineStyle === 'none' && rest.boxShadow === 'none',
        `at rest the search field has the site's own look (${JSON.stringify(rest)})`);
  const states = [];
  async function state(name) { states.push([name, await look()]); }

  await evaluate(`${field}.focus()`);
  await state('focused');
  await send('Input.insertText', { text: 'w' });
  await state('one letter typed');
  check(JSON.stringify(await offered()) === JSON.stringify(keywordsFor(idx, 'w')),
        `one letter, "w", offers every keyword with a word it starts (${(await offered() || []).length})`);
  await send('Input.insertText', { text: 'o' });
  await state('two letters typed');
  check(JSON.stringify(await offered()) === JSON.stringify(keywordsFor(idx, 'wo')), `two letters, "wo", offer ${keywordsFor(idx, 'wo')}`);
  await clearField();
  await send('Input.insertText', { text: 'c' });
  check(JSON.stringify(await offered()) === JSON.stringify(keywordsFor(idx, 'c')), `one letter, "c", offers its ${keywordsFor(idx, 'c').length} keywords`);
  await send('Input.insertText', { text: 'o' });
  const co = keywordsFor(idx, 'co');
  check(JSON.stringify(await offered()) === JSON.stringify(co), `two letters, "co", offer its ${co.length} keywords`);

  // The keywords lie under the field, inside the panel, set as the names in the list are.
  const box = await evaluate(`(function () {
    var k = document.getElementById('mfs-keywords').getBoundingClientRect(), f = ${field}.getBoundingClientRect();
    var pane = document.getElementById('mfs-list-pane').getBoundingClientRect();
    var o = getComputedStyle(document.querySelector('#mfs-keywords .mfs-keyword')), n = getComputedStyle(document.querySelector('#mfs-search-results .mfs-result:not(.mfs-result-active):not(.mfs-result-lit)'));
    var b = getComputedStyle(document.getElementById('mfs-keywords'));
    return { under: k.top >= f.bottom && k.left === f.left && k.right === f.right, inside: k.bottom <= pane.bottom + 0.5,
             same: ['fontFamily', 'fontSize', 'color', 'backgroundColor', 'letterSpacing', 'paddingLeft'].every(function (p) { return o[p] === n[p]; }),
             plain: b.boxShadow === 'none' && o.boxShadow === 'none' && o.textShadow === 'none' && b.filter === 'none' };
  })()`);
  check(box.under, 'the keywords lie under the search field, as wide as it');
  check(box.inside, 'the keywords end inside the panel');
  check(box.same, 'a keyword is set as a name in the list is');
  check(box.plain, 'nothing about the keywords glows');

  // The keyboard goes through them, the field keeping it, and Enter chooses one.
  await key('ArrowDown', 'ArrowDown', 40);
  await key('ArrowDown', 'ArrowDown', 40);
  check(JSON.stringify(await goneTo()) === JSON.stringify([co[1]]), `the arrow keys go to "${co[1]}" (${await goneTo()})`);
  check(await evaluate(`document.activeElement === ${field} && ${field}.getAttribute('aria-activedescendant') === 'mfs-keyword-1'`),
        'the field keeps the keyboard while the arrow keys go through the keywords');
  const lit = await evaluate(`(function () { var o = getComputedStyle(document.querySelector('#mfs-keywords .mfs-keyword-on')); return [o.backgroundColor, o.color]; })()`);
  check(JSON.stringify(lit) === JSON.stringify(['rgba(0, 229, 255, 0.1)', own.cyan]), `the keyword gone to is set as a name under the pointer is (${lit})`);
  await state('a keyword gone to');
  await key('ArrowUp', 'ArrowUp', 38);
  check(JSON.stringify(await goneTo()) === JSON.stringify([co[0]]), `ArrowUp goes back to "${co[0]}"`);
  await key('ArrowDown', 'ArrowDown', 40);
  await key('Enter', 'Enter', 13);
  await chosen(co[1], 'by the keyboard');
  await state(`"${co[1]}" chosen`);

  // Escape and leaving the field put them away.
  await clearField();
  await send('Input.insertText', { text: 'de' });
  check((await offered() || []).includes('de Sitter'), 'two letters, "de", offer "de Sitter"');
  await key('Escape', 'Escape', 27);
  check((await offered()) === null && await evaluate(`${field}.value === 'de' && ${field}.getAttribute('aria-expanded') === 'false'`),
        'Escape puts the keywords away and leaves what was typed');
  await key('ArrowDown', 'ArrowDown', 40);
  check(JSON.stringify(await offered()) === JSON.stringify(keywordsFor(idx, 'de')), 'ArrowDown offers them again');
  // The pointer: a press on one chooses it.
  const at = await evaluate(`(function () {
    var o = [].find.call(document.querySelectorAll('#mfs-keywords .mfs-keyword'), function (o) { return o.textContent === 'de Sitter'; });
    var r = o.getBoundingClientRect(); return [r.left + r.width / 2, r.top + r.height / 2];
  })()`);
  await point(at[0], at[1]);
  check(JSON.stringify(await goneTo()) === JSON.stringify(['de Sitter']), 'the pointer on "de Sitter" goes to it');
  await press(at[0], at[1]);
  await chosen('de Sitter', 'by the pointer');
  check(await evaluate(`document.activeElement === ${field}`), 'the field keeps the keyboard after a press on a keyword');
  await state('"de Sitter" chosen');

  await clearField();
  await state('cleared');
  await send('Input.insertText', { text: 'zzqq' });
  check((await listed()).length === 0 && (await offered()) === null, '"zzqq" finds nothing and offers nothing');
  await state('nothing found');
  await clearField();
  check((await listed()).length === idx.length, 'clearing the field lists every spacetime');
  await send('Input.insertText', { text: 'vacuum' });
  await state('names found');
  await evaluate(`${field}.blur()`);
  check((await offered()) === null, 'leaving the field puts the keywords away');
  await state('left with names found');
  await evaluate(`${field}.focus()`);
  await clearField();
  await evaluate(`${field}.blur()`);
  await state('cleared and left');
  for (const [name, seen] of states) {
    const changed = Object.keys(rest).filter(k => seen[k] !== rest[k]);
    check(changed.length === 0, `the search field keeps its colours ${name}` +
          (changed.length ? ` (${changed.map(k => `${k} ${seen[k]} for ${rest[k]}`).join(', ')})` : ''));
  }
  if (!phone) await until('window._mfsGraph.shows().settled && window._mfsGraph.shows().whole', 'the whole graph is back after the search field\'s checks');
}
// A keyword chosen is what the field holds, and the list and the graph show exactly its carriers.
async function chosen(keyword, how) {
  const carriers = (await index()).filter(m => m.tags.includes(keyword)).map(m => m.id);
  check(await evaluate(`document.getElementById('mfs-search-input').value === ${JSON.stringify(keyword)}`) && (await offered()) === null,
        `"${keyword}" chosen ${how} is put in the field and the keywords are put away`);
  check(JSON.stringify(await listed()) === JSON.stringify(carriers), `"${keyword}" chosen ${how}: the list shows exactly the ${carriers.length} spacetimes carrying it`);
  if (phone) return;
  await until(`window._mfsGraph.shows().settled && !window._mfsGraph.shows().whole && ${drawn}`, `"${keyword}" chosen ${how} is gathered`);
  const front = (await shows()).front.map(s => s.id);
  check(JSON.stringify(front) === JSON.stringify(carriers), `"${keyword}" chosen ${how}: the graph gathers exactly the spacetimes carrying it`);
}

async function press(x, y) {
  for (const type of ['mouseMoved', 'mousePressed', 'mouseReleased']) {
    await send('Input.dispatchMouseEvent', { type, x, y, button: type === 'mouseMoved' ? 'none' : 'left', clickCount: 1 });
  }
}
function wordMatches(tag, q) {
  const t = tag.toLowerCase();
  for (let i = t.indexOf(q); i !== -1; i = t.indexOf(q, i + 1)) if (i === 0 || ' -/'.includes(t[i - 1])) return true;
  return false;
}

await send('Page.navigate', { url: `${base}/MFS/` });
await until("!!window._mfsIndex && document.querySelectorAll('#mfs-search-results .mfs-result').length > 0", 'the list is drawn');
const all = (await index()).map(m => m.id);

if (phone) {
  await sleep(2500);
  const state = await evaluate(`(function () {
    var c = document.getElementById('mfs-graph'), cs = getComputedStyle(c);
    return { display: cs.display, fetched: performance.getEntriesByType('resource').some(function (e) { return /relations\\.json/.test(e.name); }) };
  })()`);
  check(state.display === 'none', 'a phone draws no graph');
  check(!state.fetched, 'a phone never fetches the relations');
  const hits = await evaluate(`(function () {
    var found = [];
    for (var y = 5; y < innerHeight; y += 40) for (var x = 5; x < innerWidth; x += 40) {
      var e = document.elementFromPoint(x, y);
      if (e && e.id === 'mfs-graph') found.push([x, y]);
    }
    return found;
  })()`);
  check(hits.length === 0, 'no tap on a phone lands on the graph');
  await searchField();
} else {
  await until(`window._mfsGraph.shows().on && window._mfsGraph.shows().settled && ${drawn}`, 'the graph is drawn');
  let state = await shows();
  check(state.whole && state.front.length === all.length, `the whole graph holds all ${all.length} spacetimes`);
  const named = state.front.filter(s => s.named).length;
  check(named > 0 && named <= 40, `forty names at most are written (${named})`);

  const room = await evaluate(`(function () {
    var g = document.getElementById('mfs-graph-room').getBoundingClientRect();
    var panels = ['mfs-search-panel', 'mfs-coffee-panel'].map(function (id) { return document.getElementById(id).getBoundingClientRect(); });
    return { left: g.left, right: g.right, top: g.top, bottom: g.bottom, panelRight: Math.max(panels[0].right, panels[1].right) };
  })()`);
  check(room.left >= room.panelRight, 'the graph is gathered to the right of the list');

  // The canvas reaches every edge of the window at every desktop size, so nothing of the graph
  // is cut off short of it, and it is drawn at the screen's own pixels.
  // Each size is a resize, since headless Chrome changes the pixel ratio without telling the page.
  for (const [width, height] of [[1280, 720], [1920, 1080], [1024, 768], [1440, 900]]) {
    await send('Emulation.setDeviceMetricsOverride', { width, height, deviceScaleFactor: 2, mobile: false });
    await until(`window.innerWidth === ${width} && window.innerHeight === ${height}`, `the window is ${width} by ${height}`);
    await until(`window._mfsGraph.shows().settled && ${drawn}`, `the graph is drawn at ${width} by ${height}`);
    const c = (await shows()).canvas, view = await evaluate('[window.innerWidth, window.innerHeight]');
    check(c.left === 0 && c.top === 0 && c.right === view[0] && c.bottom === view[1],
          `at ${width} by ${height} the graph's canvas spans the whole window (${c.left},${c.top} to ${c.right},${c.bottom})`);
    await until(`window._mfsGraph.shows().canvas.width === ${width * 2} && window._mfsGraph.shows().canvas.height === ${height * 2}`,
                `at ${width} by ${height} the graph is drawn at the screen's pixels`, 5);
  }
  await send('Emulation.setDeviceMetricsOverride', { width: 1440, height: 900, deviceScaleFactor: 1, mobile: false });
  await until('window.innerWidth === 1440 && window.innerHeight === 900', 'the window is 1440 by 900 again');
  await until(`window._mfsGraph.shows().settled && ${drawn}`, 'the graph is drawn again at 1440 by 900');
  state = await shows();
  const covered = await evaluate(`(function () {
    var found = [];
    ['mfs-search-panel', 'mfs-coffee-panel', 'mfs-exit', 'mfs-studio-sign'].forEach(function (id) {
      var r = document.getElementById(id).getBoundingClientRect();
      for (var y = r.top + 3; y < r.bottom - 3; y += 15) for (var x = r.left + 3; x < r.right - 3; x += 15) {
        var e = document.elementFromPoint(x, y);
        if (e && e.id === 'mfs-graph') found.push(id);
      }
    });
    return found;
  })()`);
  check(covered.length === 0, 'the graph covers no panel');
  // The title takes no press, so a press on it reaches the graph, which is drawn under it.
  check(await evaluate("+getComputedStyle(document.getElementById('mfs-title')).zIndex > +getComputedStyle(document.getElementById('mfs-graph')).zIndex"),
        'the title is drawn over the graph');
  check(state.front.every(s => s.x >= room.left && s.x <= room.right && s.y >= room.top && s.y <= room.bottom),
        'every spacetime is gathered in the graph\'s room');

  await asOneThing(all, 'the whole graph');
  await searchField();

  for (const keyword of keywords) {
    await type(keyword);
    await until('window._mfsGraph.shows().settled && !window._mfsGraph.shows().whole', `"${keyword}" is gathered`);
    state = await shows();
    const list = await listed();
    const carriers = (await index()).filter(m => m.name.toLowerCase().includes(keyword.toLowerCase()) ||
                                                 m.tags.some(t => wordMatches(t, keyword.toLowerCase()))).map(m => m.id);
    const front = state.front.map(s => s.id);
    check(JSON.stringify(front) === JSON.stringify(list), `"${keyword}": the graph gathers what the list shows (${front.length})`);
    check(JSON.stringify(front) === JSON.stringify(carriers), `"${keyword}": those are the spacetimes carrying it`);
    check(state.front.filter(s => s.named).length <= 40, `"${keyword}": forty names at most`);
    if (keyword === keywords[0]) {
      await asOneThing(list, `"${keyword}"`);
      // A spacetime the search left out is never lit, from the list or the graph, nor opened.
      const out = all.find(id => !list.includes(id));
      await evaluate(`window._mfsGraph.light(${JSON.stringify(out)})`);
      await sleep(100);
      check((await shows()).lit === null && (await pinkPoints()).length === 0, `"${keyword}": ${out}, left out, is never lit`);
      await evaluate('window._mfsGraph.light(null)');
      const back = await evaluate(`window._mfsGraph.shows().back`);
      // One the pointer reaches on the graph, standing clear of every point in front and every name written.
      const open = await evaluate(`window._mfsGraph.shows().back.map(function (b) { var e = document.elementFromPoint(b.x, b.y); return !!e && e.id === 'mfs-graph'; })`);
      const lone = back.find((b, k) => open[k] && state.front.every(f => Math.hypot(f.x - b.x, f.y - b.y) > 30 &&
        !(f.name && b.x > f.name[0] - 15 && b.x < f.name[0] + f.name[2] + 15 && b.y > f.name[1] - 15 && b.y < f.name[1] + f.name[3] + 15)));
      if (lone) {
        await point(lone.x, lone.y);
        await sleep(150);
        check((await shows()).lit === null && (await litNames()).length === 0, `"${keyword}": the pointer on ${lone.id}, left out, lights nothing`);
        await press(lone.x, lone.y);
        await sleep(300);
        check(!new URLSearchParams(await evaluate('location.search')).get('spacetime'), `"${keyword}": a press on ${lone.id}, left out, opens nothing`);
        await away();
      } else check(false, `"${keyword}": a spacetime left out stands clear of those in front`);
    }
    await clear();
    await until('window._mfsGraph.shows().settled && window._mfsGraph.shows().whole', `clearing "${keyword}" brings back the whole graph`);
    state = await shows();
    check(state.front.length === all.length, `cleared "${keyword}": every spacetime is back in front`);
  }

  // A press on a named spacetime opens it as its name in the list does, and the graph stays
  // where it is under the spacetime's panel, covered by it and standing still.
  const target = state.front.find(s => s.named);
  await press(target.x, target.y);
  await until(`new URLSearchParams(location.search).get('spacetime') === ${JSON.stringify(target.id)}`, `a press opens ${target.id}`);
  await until("document.getElementById('mfs-content-panel').style.pointerEvents === 'all'", 'the spacetime slides in', 30);
  check(await evaluate(`document.querySelector('#mfs-search-results .mfs-result-active').dataset.id === ${JSON.stringify(target.id)}`),
        `${target.id} is marked open in the list, as pressing its name marks it`);
  await sleep(500);
  state = await shows();
  check(state.on && state.settled && !state.drawing, 'the graph stays drawn and still under the open spacetime');
  check(JSON.stringify(state.front.map(s => [s.id, s.x, s.y])) === JSON.stringify((await shows()).front.map(s => [s.id, s.x, s.y])),
        'the graph stands where it was under the open spacetime');
  const under = await evaluate(`(function () {
    var g = document.getElementById('mfs-graph-room').getBoundingClientRect(), found = [];
    for (var y = g.top + 3; y < g.bottom - 3; y += 25) for (var x = g.left + 3; x < g.right - 3; x += 25) {
      var e = document.elementFromPoint(x, y);
      if (!e || !e.closest('#mfs-content-panel')) found.push([x, y]);
    }
    return found;
  })()`);
  check(under.length === 0, `the spacetime's panel covers the graph's room (${under.length} points not covered)`);
  await evaluate("document.getElementById('mfs-toc-back').click()");
  await until("document.getElementById('mfs-content-panel').style.pointerEvents !== 'all' && !document.getElementById('mfs-list-pane').inert", 'the list comes back');
  check((await shows()).on, 'the graph is there with the list');
}

const pageErrors = await evaluate('window.__mfsErrors');
[...errors, ...pageErrors].forEach(e => check(false, e));
console.log(failures.length ? `${failures.length} failed` : 'all held');
finish(failures.length ? 1 : 0);
