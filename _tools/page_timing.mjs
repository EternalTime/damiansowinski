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
   iPhone held upright, 390 by 844; --width and --height set another screen; --cpu <n> slows
   the processor n times over, as Chrome's own tools do to stand in for a slower phone; and
   --budget <seconds>, 10 by default, is the time past which a chart counts as failing.
   CHROME names the browser when it is not where macOS keeps it.

   It prints one line per chart, slowest last, and exits non-zero if any chart missed the
   budget or never became ready. */
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
const cpu = +option('--cpu', 1);
const budget = +option('--budget', 10);
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
socket.addEventListener('message', e => {
  const m = JSON.parse(e.data);
  if (pending.has(m.id)) { pending.get(m.id)(m); pending.delete(m.id); }
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

/* Before the page's own scripts run: count the typesetting MathJax has in hand and keep every
   task that blocked the page for more than 50ms, which is what the browser calls a long task. */
await send('Page.enable');
await send('Page.addScriptToEvaluateOnNewDocument', { source: `
  window.__mfsTiming = { busy: 0, long: [] };
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
  })();` });
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
   task. */
function timed(action, chart) {
  return evaluate(`new Promise(function (resolve) {
    var t = window.__mfsTiming, panel = document.getElementById('mfs-content-panel');
    var before = panel.querySelector('.mfs-header');
    var start = performance.now();
    ${action}
    (function poll() {
      var select = document.getElementById('mfs-coord-select'), header = panel.querySelector('.mfs-header');
      var ready = header && header !== before && t.busy === 0 &&
        (!select || select.value === '${chart}') && !/\\\\[(\\[]/.test(panel.textContent);
      if (!ready) {
        if (performance.now() - start > ${LIMIT_S * 1000}) { resolve(null); return; }
        setTimeout(poll, 25); return;
      }
      requestAnimationFrame(function () { requestAnimationFrame(function () { setTimeout(function () {
        var end = performance.now();
        var longest = t.long.filter(function (l) { return l[0] + l[1] > start; })
                            .reduce(function (m, l) { return Math.max(m, l[1]); }, 0);
        resolve({ seconds: (end - start) / 1000, longest: longest / 1000,
                  chart: select ? select.options[select.selectedIndex].text : '',
                  lines: panel.querySelectorAll('.mfs-line').length,
                  elements: panel.getElementsByTagName('*').length });
      }, 0); }); });
    })();
  })`);
}

const results = [];
for (const id of ids) {
  if (only.length && !only.includes(id)) continue;
  const first = await timed(`document.querySelector('.mfs-result[data-id="${id}"]').click();`, 0);
  results.push({ id, ...(first || { chart: '?' }), failed: !first });
  if (!first) continue;
  const charts = await evaluate(`(document.getElementById('mfs-coord-select') || { options: [1] }).options.length`);
  for (let c = 1; c < charts; c++) {
    const r = await timed(`var s = document.getElementById('mfs-coord-select'); s.value = '${c}';
                           s.dispatchEvent(new Event('change'));`, c);
    results.push({ id, ...(r || { chart: String(c) }), failed: !r });
  }
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
console.log(`${results.length} charts at ${width} by ${height}${phone ? ' as a phone' : ''}` +
  `${cpu > 1 ? `, processor slowed ${cpu} times` : ''}: ` +
  (bad ? `${bad} over the ${budget}s budget` : `all within the ${budget}s budget`));
finish(bad ? 1 : 0);
