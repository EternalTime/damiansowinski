#!/usr/bin/env node
/* Open the spacetimes page in headless Chrome as a reader does and hold the graph of relations
   behind the list to what _tools/README.md, "The graph behind the list", says it does:

     bundle exec jekyll serve
     node _tools/background_graph.mjs http://127.0.0.1:4000

   On a desktop the graph is drawn in the room to the right of the list and covers no panel,
   every spacetime is in it, and no more than forty names are written. Typing a keyword into
   the search, a key at a time, gathers exactly the spacetimes the list then shows, which are the
   ones carrying it, with the rest behind; clearing the search brings the whole graph back; a
   press on a spacetime in the graph opens it and the graph steps aside, and "‹ All spacetimes"
   brings it back. With --phone the page is laid out as an iPhone held upright, the graph is not
   drawn and its relations are never fetched. Every console error and page error is an error.
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
async function type(text) {
  await evaluate("document.getElementById('mfs-search-input').focus()");
  for (const key of text) await send('Input.insertText', { text: key });
}
async function clear() {
  await evaluate("document.getElementById('mfs-search-input').select()");
  await send('Input.dispatchKeyEvent', { type: 'keyDown', key: 'Backspace', code: 'Backspace', windowsVirtualKeyCode: 8 });
  await send('Input.dispatchKeyEvent', { type: 'keyUp', key: 'Backspace', code: 'Backspace', windowsVirtualKeyCode: 8 });
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
} else {
  await until('window._mfsGraph.shows().on && window._mfsGraph.shows().settled', 'the graph is drawn');
  let state = await shows();
  check(state.whole && state.front.length === all.length, `the whole graph holds all ${all.length} spacetimes`);
  const named = state.front.filter(s => s.named).length;
  check(named > 0 && named <= 40, `forty names at most are written (${named})`);

  const room = await evaluate(`(function () {
    var g = document.getElementById('mfs-graph').getBoundingClientRect();
    var panels = ['mfs-search-panel', 'mfs-coffee-panel'].map(function (id) { return document.getElementById(id).getBoundingClientRect(); });
    return { left: g.left, right: g.right, top: g.top, bottom: g.bottom, panelRight: Math.max(panels[0].right, panels[1].right) };
  })()`);
  check(room.left >= room.panelRight, 'the graph lies to the right of the list');
  const covered = await evaluate(`(function () {
    var found = [];
    ['mfs-search-panel', 'mfs-coffee-panel', 'mfs-title', 'mfs-exit'].forEach(function (id) {
      var r = document.getElementById(id).getBoundingClientRect();
      for (var y = r.top + 3; y < r.bottom - 3; y += 15) for (var x = r.left + 3; x < r.right - 3; x += 15) {
        var e = document.elementFromPoint(x, y);
        if (e && e.id === 'mfs-graph') found.push(id);
      }
    });
    return found;
  })()`);
  check(covered.length === 0, 'the graph covers no panel');
  check(state.front.every(s => s.x >= room.left && s.x <= room.right && s.y >= room.top && s.y <= room.bottom),
        'every spacetime is drawn in the graph\'s room');

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
    await clear();
    await until('window._mfsGraph.shows().settled && window._mfsGraph.shows().whole', `clearing "${keyword}" brings back the whole graph`);
    state = await shows();
    check(state.front.length === all.length, `cleared "${keyword}": every spacetime is back in front`);
  }

  // A press on a named spacetime opens it, and the graph steps aside until the list is back.
  const target = state.front.find(s => s.named);
  await press(target.x, target.y);
  await until(`new URLSearchParams(location.search).get('spacetime') === ${JSON.stringify(target.id)}`, `a press opens ${target.id}`);
  await until('!window._mfsGraph.shows().on', 'the graph steps aside while a spacetime is open');
  await until("document.getElementById('mfs-content-panel').style.pointerEvents === 'all'", 'the spacetime slides in', 30);
  await evaluate("document.getElementById('mfs-toc-back').click()");
  await until('window._mfsGraph.shows().on', 'the graph comes back with the list');
}

const pageErrors = await evaluate('window.__mfsErrors');
[...errors, ...pageErrors].forEach(e => check(false, e));
console.log(failures.length ? `${failures.length} failed` : 'all held');
finish(failures.length ? 1 : 0);
