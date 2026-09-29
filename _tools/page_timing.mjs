#!/usr/bin/env node
/* Open every spacetime on the spacetimes page in headless Chrome, as a reader does, and time
   how long each takes to be ready: from the click on its name in the list, or the choice of a
   chart in the selector, until its mathematics is typeset and the page answers again.

   It measures in the page itself, so the time includes whatever keeps the tab busy after
   MathJax is done, which is where the Natário warp drive's minute went until 25 September
   2026, and it reports the longest single task that blocked the page, which is the thing a
   reader feels as a frozen tab. It needs Chrome and a served copy of the site, and nothing
   else; Node 22 or later carries the WebSocket it talks to Chrome over:

     bundle exec jekyll serve
     node _tools/page_timing.mjs http://127.0.0.1:4000

   --metric <id> times one spacetime and is repeatable; --phone lays the page out as an
   iPhone held upright, 390 by 844; --width and --height set another screen; --text <px>
   sets the browser's default font size, 16 by default, as a reader who enlarges text does;
   --cpu <n> slows the processor n times over, as Chrome's own tools do to stand in for a
   slower phone; and --budget <seconds>, 10 by default, is the time past which a chart
   counts as failing. CHROME names the browser when it is not where macOS keeps it.

   It also watches every number and label of every drawing shown, and the plot of each
   spacetime diagram, from the first frame the reader sees them until two frames after the
   last font has loaded, as each spacetime opens, as each chart is chosen and as each view of
   a diagram is chosen, and counts any that moves or is shown or hidden after that first frame
   as an error: until 29 September 2026 the first spacetime to use one of MathJax's fonts,
   Bianchi's fraktur or Malament-Hogarth's small capitals, fitted its numbers in the face
   drawn in its place and moved them by up to 4 pixels when the font arrived.

   In each of those states it holds every number, name of an axis and label of a spacetime
   diagram shown to the margin of 1.5 em of the diagram's words that the page gives them, and
   counts any that stands nearer the edge of the diagram's ground as an error. --print does the
   same for the print copy of each chart, as the page lays it out for paper once its print
   button is pressed, the browser's print dialog left closed.

   It prints one line per chart, slowest last, then every page error and console error the
   page raised and every label that moved, and exits non-zero if any chart missed the budget
   or never became ready, or if the page raised any error at all. */
import { spawn } from 'node:child_process';
import { mkdtempSync, readFileSync, rmSync, existsSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';

const argv = process.argv.slice(2);
function option(name, fallback) {
  const i = argv.indexOf(name);
  return i >= 0 ? argv[i + 1] : fallback;
}
const base = (argv[0] && !argv[0].startsWith('--') ? argv[0] : 'http://127.0.0.1:4000').replace(/\/$/, '');
const only = argv.flatMap((a, i) => (a === '--metric' ? [argv[i + 1]] : []));
const phone = argv.includes('--phone');
const width = +option('--width', phone ? 390 : 1440);
const height = +option('--height', phone ? 844 : 900);
const text = +option('--text', 16);
const cpu = +option('--cpu', 1);
const budget = +option('--budget', 10);
const print = argv.includes('--print');
const LIMIT_S = Math.max(300, budget * 10);

const chromePath = process.env.CHROME ||
  ['/Applications/Google Chrome.app/Contents/MacOS/Google Chrome', '/usr/bin/google-chrome', '/usr/bin/chromium']
    .find(existsSync);
if (!chromePath) { console.error('No Chrome found; set CHROME to its path.'); process.exit(2); }

const profile = mkdtempSync(join(tmpdir(), 'mfs-page-timing-'));
const chrome = spawn(chromePath, ['--headless=new', '--remote-debugging-port=0', `--user-data-dir=${profile}`,
  '--no-first-run', '--no-default-browser-check', '--hide-scrollbars', 'about:blank'], { stdio: 'ignore' });
function finish(code) {
  chrome.kill('SIGKILL');
  try { rmSync(profile, { recursive: true, force: true }); } catch {}
  process.exit(code);
}
// Whatever stops the run, Chrome goes with it.
process.on('uncaughtException', e => { console.error(e); finish(2); });
process.on('unhandledRejection', e => { console.error(e); finish(2); });
const sleep = ms => new Promise(r => setTimeout(r, ms));

// Chrome writes the port it chose into the profile once it listens.
let port = null;
for (let i = 0; i < 200 && !port; i++) {
  try { port = readFileSync(join(profile, 'DevToolsActivePort'), 'utf8').split('\n')[0]; } catch { await sleep(50); }
}
if (!port) { console.error('Chrome did not start.'); finish(2); }
const targets = await (await fetch(`http://127.0.0.1:${port}/json`)).json();
const socket = new WebSocket(targets.find(t => t.type === 'page').webSocketDebuggerUrl);
await new Promise(r => socket.addEventListener('open', r));
let seq = 0;
const pending = new Map();
/* Every console.error and every error the browser logs itself, such as a file that failed to
   load, with the spacetime open when it came. The page's own errors are gathered in the page,
   below, since Chrome reports "ResizeObserver loop completed with undelivered notifications"
   to the page's error handlers alone. */
const errors = [];
let opened = 'the list';
socket.addEventListener('message', e => {
  const m = JSON.parse(e.data);
  if (pending.has(m.id)) { pending.get(m.id)(m); pending.delete(m.id); }
  if (m.method === 'Runtime.consoleAPICalled' && m.params.type === 'error') {
    errors.push(`${opened}: console error: ${m.params.args.map(a => a.value ?? a.description).join(' ')}`);
  } else if (m.method === 'Log.entryAdded' && m.params.entry.level === 'error') {
    errors.push(`${opened}: console error: ${m.params.entry.text}${m.params.entry.url ? ` (${m.params.entry.url})` : ''}`);
  }
});
function send(method, params = {}) {
  return new Promise(r => { const id = ++seq; pending.set(id, r); socket.send(JSON.stringify({ id, method, params })); });
}
async function evaluate(expression) {
  const m = await send('Runtime.evaluate', { expression, returnByValue: true, awaitPromise: true });
  const r = m.result || {};
  if (r.exceptionDetails) throw new Error(r.exceptionDetails.exception?.description || r.exceptionDetails.text);
  return r.result?.value;
}

await send('Emulation.setDeviceMetricsOverride', { width, height, deviceScaleFactor: 1, mobile: phone });
if (phone) await send('Emulation.setTouchEmulationEnabled', { enabled: true, maxTouchPoints: 5 });
if (cpu > 1) await send('Emulation.setCPUThrottlingRate', { rate: cpu });
if (text !== 16) await send('Page.setFontSizes', { fontSizes: { standard: text } });
await send('Runtime.enable');
await send('Log.enable');

/* Before the page's own scripts run: count the typesetting MathJax has in hand, keep every
   task that blocked the page for more than 50ms, which is what the browser calls a long task,
   and keep every error the page raises. */
await send('Page.enable');
await send('Page.addScriptToEvaluateOnNewDocument', { source: `
  window.__mfsTiming = { busy: 0, long: [], errors: [] };
  window.addEventListener('error', function (e) {
    window.__mfsTiming.errors.push(e.message + (e.filename ? ' (' + e.filename + ':' + e.lineno + ')' : ''));
  });
  window.addEventListener('unhandledrejection', function (e) {
    window.__mfsTiming.errors.push('unhandled rejection: ' + (e.reason && e.reason.stack || e.reason));
  });
  try {
    new PerformanceObserver(function (list) {
      list.getEntries().forEach(function (e) { window.__mfsTiming.long.push([e.startTime, e.duration]); });
    }).observe({ entryTypes: ['longtask'] });
  } catch (e) {}
  (function wrap() {
    if (!(window.MathJax && MathJax.typesetPromise && MathJax.startup && MathJax.startup.document)) { setTimeout(wrap, 20); return; }
    var typeset = MathJax.typesetPromise;
    MathJax.typesetPromise = function () {
      window.__mfsTiming.busy++;
      function done() { window.__mfsTiming.busy--; }
      return typeset.apply(this, arguments).then(function (v) { done(); return v; }, function (e) { done(); throw e; });
    };
    window.__mfsTiming.ready = true;
  })();
  /* Where every number and label of every drawing shown stands, and each spacetime diagram's
     plot, from the top left of its figure, since the panel slides in, and whether it is shown. */
  window.__mfsLabels = function (panel) {
    var out = {};
    [].forEach.call(panel.querySelectorAll('.mfs-nr-figure, .mfs-cd-figure'), function (fig, i) {
      if (!fig.getClientRects().length) return;
      var b = fig.getBoundingClientRect();
      [].forEach.call(fig.querySelectorAll('.nr-plot, .nr-tick, .mfs-slice-label, .cd-at'), function (el, j) {
        if (!el.getClientRects().length) return;
        var r = el.getBoundingClientRect();
        out[i + ' ' + j] = { name: el.className.split(' ')[0] + ' "' + el.textContent.trim().slice(0, 24) + '"',
          at: [r.left - b.left, r.top - b.top, r.right - b.left, r.bottom - b.top],
          shown: getComputedStyle(el).visibility !== 'hidden' };
      });
    });
    return out;
  };
  /* From the first frame at which shown() holds, every frame is held against that first one
     until stop() is called, and each label that moved by more than a twentieth of a pixel, the
     noise of the panel's slide, or was shown or hidden, is kept once with how far it went. */
  window.__mfsWatch = function (panel, shown) {
    var first = null, start = 0, moved = {}, on = true;
    function compare() {
      var now = window.__mfsLabels(panel);
      Object.keys(now).forEach(function (k) {
        var a = first[k], b = now[k];
        if (!a || moved[k]) return;
        var d = Math.max.apply(null, a.at.map(function (v, i) { return Math.abs(v - b.at[i]); }));
        if (d > 0.05 || a.shown !== b.shown) {
          moved[k] = a.name + (a.shown !== b.shown ? (b.shown ? ' shown' : ' hidden') : ' moved ' + d.toFixed(2) + 'px') +
            ' at ' + Math.round(performance.now() - start) + 'ms after first shown';
        }
      });
    }
    (function frame() {
      if (!on) return;
      if (first) compare();
      else if (shown()) { first = window.__mfsLabels(panel); start = performance.now(); }
      requestAnimationFrame(frame);
    })();
    return function stop() {
      on = false;
      if (first) compare();
      return Object.keys(moved).map(function (k) { return moved[k]; });
    };
  };
  /* Every number, name of an axis and label of each spacetime diagram shown in root that
     stands nearer the edge of the diagram's ground than 1.5 em of its words, the margin
     .mfs-nr-figure gives them, with how near it stands, in those ems. */
  window.__mfsMargins = function (root) {
    var out = [];
    [].forEach.call(root.querySelectorAll('.mfs-nr-figure'), function (fig) {
      if (!fig.getClientRects().length) return;
      var b = fig.getBoundingClientRect(), em = parseFloat(getComputedStyle(fig).fontSize);
      [].forEach.call(fig.querySelectorAll('.nr-tick, .nr-name, .nr-ref, .mfs-slice-label'), function (el) {
        if (!el.getClientRects().length || getComputedStyle(el).visibility === 'hidden') return;
        var r = el.getBoundingClientRect();
        var near = Math.min(r.left - b.left, r.top - b.top, b.right - r.right, b.bottom - r.bottom);
        // Less half a pixel, to which a browser may round where a turned or shifted box stands.
        if (near < 1.5 * em - 0.5) {
          out.push(el.className.split(' ')[0] + ' "' + el.textContent.trim().slice(0, 24) + '" stands ' +
            (near / em).toFixed(2) + ' em from the edge of its diagram, inside the margin of 1.5 em');
        }
      });
    });
    return out;
  };
  window.__mfsSettled = function () {
    return document.fonts.ready.then(function () {
      return new Promise(function (r) { requestAnimationFrame(function () { requestAnimationFrame(r); }); });
    });
  };` });
await send('Page.navigate', { url: `${base}/MFS/` });
for (let i = 0; i < 600; i++) {
  if (await evaluate(`!!(window.__mfsTiming && window.__mfsTiming.ready && document.querySelector('.mfs-result'))`)) break;
  await sleep(100);
}
const ids = await evaluate(`[].map.call(document.querySelectorAll('.mfs-result'), function (e) { return e.dataset.id; })`);
if (!ids || !ids.length) { console.error(`No spacetimes listed at ${base}/MFS/.`); finish(2); }
for (const id of only) if (!ids.includes(id)) { console.error(`No spacetime ${id} in the list.`); finish(2); }

/* Act, then resolve once the panel has been drawn again, holds the chart asked for with every
   formula typeset and no typesetting outstanding, and two frames have passed, so the
   observers that run once MathJax is done have run too. Every wait is a timer inside the
   page, so while a task holds the page nothing here moves, and the time taken includes that
   task. The labels are watched from the first frame at which the panel stands in place with
   its mathematics set and its figures in sight, until the page has settled and the fonts are
   in; that wait is not timed. */
function timed(action, chart) {
  return evaluate(`new Promise(function (resolve) {
    var t = window.__mfsTiming, panel = document.getElementById('mfs-content-panel');
    var before = panel.querySelector('.mfs-header');
    var stop = window.__mfsWatch(panel, function () {
      var header = panel.querySelector('.mfs-header');
      return header && header !== before && t.busy === 0 && !panel.classList.contains('mfs-fonts-wait') &&
        /^translateX\\(0(px)?\\)$/.test(panel.style.transform);
    });
    var start = performance.now();
    ${action}
    (function poll() {
      var on = panel.querySelector('.mfs-charts .mfs-choice-on'), header = panel.querySelector('.mfs-header');
      var ready = header && header !== before && t.busy === 0 &&
        (!on || on.dataset.chart === '${chart}') && !/\\\\[(\\[]/.test(panel.textContent);
      if (!ready) {
        if (performance.now() - start > ${LIMIT_S * 1000}) { stop(); resolve(null); return; }
        setTimeout(poll, 25); return;
      }
      requestAnimationFrame(function () { requestAnimationFrame(function () { setTimeout(function () {
        var end = performance.now();
        var longest = t.long.filter(function (l) { return l[0] + l[1] > start; })
                            .reduce(function (m, l) { return Math.max(m, l[1]); }, 0);
        var result = { seconds: (end - start) / 1000, longest: longest / 1000,
                       chart: on ? on.textContent : '',
                       lines: panel.querySelectorAll('.mfs-line').length,
                       elements: panel.getElementsByTagName('*').length };
        window.__mfsSettled().then(function () {
          result.moved = stop();
          result.margins = window.__mfsMargins(panel);
          resolve(result);
        });
      }, 0); }); });
    })();
  })`);
}

/* Choose each view of each diagram but the one first shown, and watch its labels from the
   frame after the choice until the page has settled. */
async function views() {
  const moved = await evaluate(`(async function () {
    var panel = document.getElementById('mfs-content-panel'), moved = [];
    var buttons = [].slice.call(panel.querySelectorAll('.mfs-diagram > .mfs-section-head .mfs-choice'));
    for (var i = 0; i < buttons.length; i++) {
      if (buttons[i].classList.contains('mfs-choice-on')) continue;
      var button = buttons[i], clicked = false;
      var stop = window.__mfsWatch(panel, function () { return clicked; });
      button.click();
      clicked = true;
      await window.__mfsSettled();
      await new Promise(function (r) { requestAnimationFrame(r); });
      stop().forEach(function (m) { moved.push(button.textContent.trim() + ': ' + m); });
      window.__mfsMargins(panel).forEach(function (m) { moved.push(button.textContent.trim() + ': ' + m); });
    }
    return moved;
  })()`);
  for (const m of moved) errors.push(`${opened}: view ${m}`);
}

/* Press the print button with the browser's print dialog held closed, wait until the print
   copy is set, and hold its spacetime diagrams to their margin as paper lays them out. */
async function printCopy() {
  const ready = await evaluate(`new Promise(function (resolve) {
    var start = performance.now(), printed = false;
    window.print = function () { printed = true; resolve(true); };
    document.getElementById('mfs-print-btn').click();
    (function wait() {
      if (printed) return;
      if (performance.now() - start > ${LIMIT_S * 1000}) { resolve(false); return; }
      setTimeout(wait, 50);
    })();
  })`);
  if (!ready) { errors.push(`${opened}: the print copy was never set`); return; }
  await send('Emulation.setEmulatedMedia', { media: 'print' });
  for (const m of await evaluate(`window.__mfsMargins(document.getElementById('mfs-print-root'))`)) {
    errors.push(`${opened}: print: ${m}`);
  }
  await send('Emulation.setEmulatedMedia', { media: '' });
}

async function pageErrors() {
  for (const e of await evaluate(`window.__mfsTiming.errors.splice(0)`)) errors.push(`${opened}: page error: ${e}`);
}
await pageErrors();

const results = [];
for (const id of ids) {
  if (only.length && !only.includes(id)) continue;
  opened = id;
  const first = await timed(`document.querySelector('.mfs-result[data-id="${id}"]').click();`, 0);
  results.push({ id, ...(first || { chart: '?' }), failed: !first });
  if (!first) continue;
  for (const m of first.moved.concat(first.margins)) errors.push(`${opened} / ${first.chart}: ${m}`);
  await views();
  if (print) await printCopy();
  const charts = await evaluate(`document.querySelectorAll('#mfs-content-panel .mfs-charts .mfs-choice').length || 1`);
  for (let c = 1; c < charts; c++) {
    const r = await timed(`document.querySelector('#mfs-content-panel .mfs-charts [data-chart="${c}"]').click();`, c);
    results.push({ id, ...(r || { chart: String(c) }), failed: !r });
    if (!r) continue;
    for (const m of r.moved.concat(r.margins)) errors.push(`${opened} / ${r.chart}: ${m}`);
    await views();
    if (print) await printCopy();
  }
  await pageErrors();
}

results.sort((a, b) => (a.failed - b.failed) || (a.seconds - b.seconds));
let bad = 0;
for (const r of results) {
  const over = r.failed || r.seconds > budget;
  if (over) bad++;
  const name = `${r.id} / ${r.chart}`.padEnd(44);
  console.log(r.failed
    ? `${name}  never ready within ${LIMIT_S}s`
    : `${name} ${r.seconds.toFixed(2).padStart(6)}s  longest task ${r.longest.toFixed(2).padStart(5)}s  ` +
      `${String(r.lines).padStart(4)} lines ${String(r.elements).padStart(7)} elements${over ? '  OVER BUDGET' : ''}`);
}
for (const e of errors) console.log(e);
console.log(`${results.length} charts at ${width} by ${height}${phone ? ' as a phone' : ''}` +
  `${text !== 16 ? ` with ${text}px text` : ''}${print ? ' and their print copies' : ''}${cpu > 1 ? `, processor slowed ${cpu} times` : ''}: ` +
  (bad ? `${bad} over the ${budget}s budget` : `all within the ${budget}s budget`) +
  `, ${errors.length ? `${errors.length} errors` : 'no errors'}`);
finish(bad || errors.length ? 1 : 0);
