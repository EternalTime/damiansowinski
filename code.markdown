---
layout: default
title: Code
---

<a href="{{ '/' | relative_url }}" style="float: right;">← Back to Home</a>
<br>

<h1 class="gradient-text2" style="font-size: 2em; text-align: left; filter: drop-shadow(1px 1px 0px rgba(0,0,0,0.9));">
  Code Repositories
</h1>

<div>
Open-source libraries behind the physics. Pick a repository and its description
and last update load into the terminal, with links to the documentation
<svg viewBox="0 0 24 24" width="14" height="14" style="fill: var(--black); vertical-align: -2px;" aria-hidden="true"><path d="M10 20v-6h4v6h5v-8h3L12 3 2 12h3v8z"/></svg>
and its GitHub
<svg viewBox="0 0 16 16" width="14" height="14" style="fill: var(--black); vertical-align: -2px;" aria-hidden="true"><path d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.01 8.01 0 0016 8c0-4.42-3.58-8-8-8z"/></svg>
source.
</div><br>

<style>
.repo-explorer { display: flex; gap: 20px; margin-top: 20px; align-items: flex-start; }
.repo-buttons { display: flex; flex-direction: column; flex: 0 0 auto; min-width: 230px; }
.repo-buttons .applet-launch-btn { width: 100%; margin: 6px 0; text-align: center; }
.repo-buttons .applet-launch-btn.active {
  background: var(--bg-active-teal);
  border-color: var(--teal-light);
  box-shadow: 0 0 16px rgba(var(--teal-light-rgb), 0.8);
}
.repo-terminal {
  --neon: var(--teal-light);
  flex: 1 1 auto;
  min-width: 0;
  border: 1px solid var(--neon);
  border-radius: 8px;
  overflow: hidden;
  background: var(--bg-void);
  box-shadow: 0 0 16px rgba(0,0,0,0.5), 0 0 16px color-mix(in srgb, var(--neon) 45%, transparent);
  font-family: 'Source Code Pro', monospace;
}
.rt-bar { display: flex; align-items: center; gap: 8px; padding: 8px 14px; background: var(--bg-panel); border-bottom: 1px solid var(--neon); }
.rt-glyph { color: var(--neon); font-weight: 700; letter-spacing: 1px; }
.rt-title { color: var(--neon); font-size: 0.85em; letter-spacing: 0.03em; }
.rt-body { padding: 16px 18px; color: var(--text-light); font-size: 0.9em; line-height: 1.6; min-height: 250px; }
.rt-prompt, .rt-key { color: var(--neon); }
.rt-out { display: block; color: var(--text-dim); white-space: pre-wrap; }
.rt-btns { display: flex; justify-content: center; gap: 12px; padding: 0 18px 18px; }
.rt-btns:empty { display: none; }
.rt-btns a {
  display: inline-flex; align-items: center; gap: 6px;
  text-decoration: none !important; text-shadow: none;
  color: var(--neon) !important;
  border: 1px solid var(--neon); border-radius: 6px;
  padding: 5px 12px; font-size: 0.85em;
  transition: background 0.15s, color 0.15s;
}
.rt-btns a:hover { background: var(--neon); color: var(--bg-void) !important; }
.rt-btns a svg { display: block; }
@media (max-width: 620px) {
  .repo-explorer { flex-direction: column; }
  .repo-buttons { flex-direction: row; flex-wrap: wrap; min-width: 0; }
}
</style>

<div class="repo-explorer">
  <div class="repo-buttons" id="repo-buttons"></div>
  <div class="repo-terminal" id="repo-terminal">
    <div class="rt-bar">
      <span class="rt-glyph">&gt;_</span>
      <span class="rt-title" id="rt-title">repositories</span>
    </div>
    <div class="rt-body" id="rt-body"><span class="rt-out">select a repository →</span></div>
    <div class="rt-btns" id="rt-btns"></div>
  </div>
</div>

<script>
(function () {
  const GH_USER = 'EternalTime';
  const GH_PATH = 'M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.01 8.01 0 0016 8c0-4.42-3.58-8-8-8z';

  // label = display name you control; neon = a site palette.css variable.
  const REPOS = [
    { repo: 'CellularAutomata', label: 'Cellular Automata',       docs: '/pyCA/' },
    { repo: 'pyCE',             label: 'Configurational Entropy', docs: '/pyCE/' },
    { repo: 'pyCoop',           label: 'Cooperation Games',       docs: '/pyCoop/' },
    { repo: 'pyEDW',            label: 'Exo-Daisy World',         docs: '/pyEDW/' },
    { repo: 'pyLEAFS',          label: 'Foraging Simulator',      docs: '/LEAFS/', gh: 'LEAFS' },
    { repo: 'pyGD',             label: 'Graph Dynamics',          docs: '/pyGD/' }
  ];

  const btnWrap = document.getElementById('repo-buttons');
  const term  = document.getElementById('repo-terminal');
  const title = document.getElementById('rt-title');
  const body  = document.getElementById('rt-body');
  const btns  = document.getElementById('rt-btns');
  const cache = {};
  let typer = null;

  function esc(s) { return (s || '—').replace(/[&<>]/g, function (c) { return { '&': '&amp;', '<': '&lt;', '>': '&gt;' }[c]; }); }
  function fmtDate(iso) { return iso ? new Date(iso).toLocaleDateString('en-US', { year: 'numeric', month: 'short', day: 'numeric' }) : 'unknown'; }

  function typeLines(lines) {
    if (typer) clearTimeout(typer);
    body.innerHTML = '';
    let i = 0;
    (function tick() {
      if (i >= lines.length) return;
      body.insertAdjacentHTML('beforeend', lines[i] + '\n');
      i++; typer = setTimeout(tick, 150);
    })();
  }

  function render(cfg, d) {
    typeLines([
      '<span class="rt-out"><span class="rt-key">description</span> : ' + esc(d.description) + '</span>',
      '<span class="rt-out"><span class="rt-key">updated</span>     : ' + fmtDate(d.pushed_at) + '</span>'
    ]);
  }

  function select(cfg, btn) {
    document.querySelectorAll('#repo-buttons .applet-launch-btn').forEach(function (b) { b.classList.remove('active'); });
    btn.classList.add('active');
    title.textContent = cfg.repo;
    const gh = cfg.gh || cfg.repo;
    btns.innerHTML =
      '<a href="' + cfg.docs + '" target="_blank" rel="noopener" aria-label="Documentation"><svg viewBox="0 0 24 24" width="15" height="15" fill="currentColor" aria-hidden="true"><path d="M10 20v-6h4v6h5v-8h3L12 3 2 12h3v8z"/></svg></a>' +
      '<a href="https://github.com/' + GH_USER + '/' + gh + '" target="_blank" rel="noopener" aria-label="GitHub"><svg viewBox="0 0 16 16" width="15" height="15" fill="currentColor" aria-hidden="true"><path d="' + GH_PATH + '"/></svg></a>';
    body.innerHTML = '<span class="rt-out">loading…</span>';
    if (cache[gh]) { render(cfg, cache[gh]); return; }
    fetch('https://api.github.com/repos/' + GH_USER + '/' + gh)
      .then(function (r) { if (!r.ok) throw new Error(r.status); return r.json(); })
      .then(function (d) { cache[gh] = d; render(cfg, d); })
      .catch(function () { body.innerHTML = '<span class="rt-out" style="color:var(--red)">could not reach github — use the links below.</span>'; });
  }

  REPOS.forEach(function (cfg) {
    const b = document.createElement('button');
    b.className = 'applet-launch-btn';
    b.textContent = cfg.label;
    b.addEventListener('click', function () { select(cfg, b); });
    btnWrap.appendChild(b);
  });
})();
</script>
