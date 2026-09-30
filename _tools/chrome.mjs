/* Headless Chrome for the checks that open the spacetimes page as a reader does,
   _tools/page_timing.mjs and _tools/turn_drag.mjs, driven over the DevTools protocol with the
   WebSocket Node 22 and later carry, so they need Chrome and nothing else. CHROME names the
   browser when it is not where macOS keeps it.

   launch() starts Chrome on a profile of its own, opens one tab at `width` by `height`, as a
   phone held upright with touch where `phone` is set, with `text` as the browser's default font
   size and the processor slowed `cpu` times, and gives back send() and evaluate() for that tab,
   onError() to hear every console.error and every error the browser logs itself, such as a file
   that failed to load, and finish(), which quits Chrome, removes its profile and exits. Whatever
   stops the run, Chrome goes with it. */
import { spawn } from 'node:child_process';
import { mkdtempSync, mkdirSync, readFileSync, rmSync, existsSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';

export const sleep = ms => new Promise(r => setTimeout(r, ms));

export async function launch({ width = 1440, height = 900, phone = false, text = 16, cpu = 1 } = {}) {
  const chromePath = process.env.CHROME ||
    ['/Applications/Google Chrome.app/Contents/MacOS/Google Chrome', '/usr/bin/google-chrome', '/usr/bin/chromium']
      .find(existsSync);
  if (!chromePath) { console.error('No Chrome found; set CHROME to its path.'); process.exit(2); }

  // The profile's folder ends in .noindex, which keeps Spotlight from indexing what Chrome
  // writes there while it runs.
  const home = mkdtempSync(join(tmpdir(), 'mfs-chrome-')), profile = join(home, 'profile.noindex');
  mkdirSync(profile);
  const chrome = spawn(chromePath, ['--headless=new', '--remote-debugging-port=0', `--user-data-dir=${profile}`,
    '--no-first-run', '--no-default-browser-check', '--hide-scrollbars', 'about:blank'], { stdio: 'ignore' });
  function finish(code) {
    chrome.kill('SIGKILL');
    try { rmSync(home, { recursive: true, force: true }); } catch {}
    process.exit(code);
  }
  process.on('uncaughtException', e => { console.error(e); finish(2); });
  process.on('unhandledRejection', e => { console.error(e); finish(2); });

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
  const pending = new Map(), heard = [];
  socket.addEventListener('message', e => {
    const m = JSON.parse(e.data);
    if (pending.has(m.id)) { pending.get(m.id)(m); pending.delete(m.id); }
    if (m.method === 'Runtime.consoleAPICalled' && m.params.type === 'error') {
      heard.forEach(f => f(`console error: ${m.params.args.map(a => a.value ?? a.description).join(' ')}`));
    } else if (m.method === 'Log.entryAdded' && m.params.entry.level === 'error') {
      heard.forEach(f => f(`console error: ${m.params.entry.text}${m.params.entry.url ? ` (${m.params.entry.url})` : ''}`));
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
  await send('Page.enable');
  return { send, evaluate, finish, onError: f => heard.push(f) };
}
