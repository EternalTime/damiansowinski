---
layout: default
title: Code
---

<link rel="stylesheet" href="{{ '/assets/css/repo-palette.css' | relative_url }}">

<a href="{{ '/' | relative_url }}" style="float: right;">← Back to Home</a>
<br>

<h1 class="gradient-text2" style="font-size: 2em; text-align: left; filter: drop-shadow(1px 1px 0px rgba(0,0,0,0.9));">
  Code Repositories
</h1>

<div>
Open-source libraries behind the physics. Each panel pulls its description and
last-updated date from GitHub; the buttons open the docs and the source in a
new tab.
</div><br>

<style>
@font-face {
  font-family: 'LatoDocs';
  src: url('{{ "/pyCE/_static/css/fonts/lato-normal.woff2" | relative_url }}') format('woff2');
  font-weight: 400; font-style: normal; font-display: swap;
}
@font-face {
  font-family: 'LatoDocs';
  src: url('{{ "/pyCE/_static/css/fonts/lato-bold.woff2" | relative_url }}') format('woff2');
  font-weight: 700; font-style: normal; font-display: swap;
}
@font-face {
  font-family: 'RobotoSlabDocs';
  src: url('{{ "/pyCE/_static/css/fonts/Roboto-Slab-Regular.woff2" | relative_url }}') format('woff2');
  font-weight: 400; font-display: swap;
}
@font-face {
  font-family: 'RobotoSlabDocs';
  src: url('{{ "/pyCE/_static/css/fonts/Roboto-Slab-Bold.woff2" | relative_url }}') format('woff2');
  font-weight: 700; font-display: swap;
}
.doc-cards { display: flex; flex-direction: column; gap: 22px; margin-top: 20px; }
.doc-card {
  --upper: var(--blue);
  --lower: var(--blue-mid);
  display: flex;
  background: var(--white);
  border: 1px solid var(--border-dark);
  border-radius: 8px;
  overflow: hidden;
  box-shadow: 0 1px 5px rgba(0,0,0,0.18);
}
.doc-spine { flex: 0 0 14px; background: var(--lower); }
.doc-body { flex: 1; min-width: 0; font-family: 'LatoDocs', 'Lato', 'Helvetica Neue', Arial, sans-serif; }
.doc-header { background: var(--upper); padding: 12px 20px; }
.doc-title { font-family: 'RobotoSlabDocs', 'Roboto Slab', Georgia, serif; font-size: 1.3em; font-weight: 700; color: var(--white); line-height: 1.2; margin: 0; }
.doc-content { padding: 14px 20px 16px; }
.doc-desc { color: var(--text-body); font-size: 0.98em; line-height: 1.5; margin-bottom: 10px; }
.doc-updated { color: var(--text-muted); font-size: 0.8em; margin-bottom: 14px; }
.doc-btns { display: flex; gap: 10px; }
.doc-btns a {
  text-decoration: none;
  font-family: 'LatoDocs', 'Lato', 'Helvetica Neue', Arial, sans-serif;
  font-size: 0.85em;
  color: var(--upper);
  border: 1px solid var(--upper);
  border-radius: 5px;
  padding: 5px 14px;
  transition: background 0.15s, color 0.15s;
}
.doc-btns a:hover { background: var(--upper); color: var(--white); }
</style>

<div class="doc-cards" id="doc-cards"></div>

<script>
(function () {
  const GH_USER = 'EternalTime';
  // label = the display name you control (not pulled from the repo).
  const REPOS = [
    { repo: 'pyCE',             label: 'Configurational Entropy', docs: '/pyCE/',   upper: '--pyce-primary',   lower: '--pyce-nav' },
    { repo: 'CellularAutomata', label: 'Cellular Automata',       docs: '/pyCA/',   upper: '--pyca-primary',   lower: '--pyca-nav' },
    { repo: 'pyGD',             label: 'Graph Dynamics',          docs: '/pyGD/',   upper: '--pygd-primary',   lower: '--pygd-nav' },
    { repo: 'pyCoop',           label: 'Cooperation Games',       docs: '/pyCoop/', upper: '--pycoop-primary', lower: '--pycoop-nav' },
    { repo: 'pyLEAFS',          label: 'Foraging Simulator',      docs: '/LEAFS/',  upper: '--pyleafs-primary', lower: '--pyleafs-nav', gh: 'LEAFS' }
  ];

  const container = document.getElementById('doc-cards');

  function fmtDate(iso) {
    if (!iso) return '';
    const d = new Date(iso);
    return 'Updated ' + d.toLocaleDateString('en-US', { year: 'numeric', month: 'short', day: 'numeric' });
  }

  REPOS.forEach(function (cfg) {
    const gh = cfg.gh || cfg.repo;
    const card = document.createElement('div');
    card.className = 'doc-card';
    card.style.setProperty('--upper', 'var(' + cfg.upper + ')');
    card.style.setProperty('--lower', 'var(' + cfg.lower + ')');
    card.innerHTML =
      '<div class="doc-spine"></div>' +
      '<div class="doc-body">' +
        '<div class="doc-header"><div class="doc-title"></div></div>' +
        '<div class="doc-content">' +
          '<div class="doc-desc" data-desc>Loading…</div>' +
          '<div class="doc-updated" data-updated></div>' +
          '<div class="doc-btns">' +
            '<a href="' + cfg.docs + '" target="_blank" rel="noopener">Docs</a>' +
            '<a href="https://github.com/' + GH_USER + '/' + gh + '" target="_blank" rel="noopener">GitHub</a>' +
          '</div>' +
        '</div>' +
      '</div>';
    card.querySelector('.doc-title').textContent = cfg.label;
    container.appendChild(card);

    const desc = card.querySelector('[data-desc]');
    const updated = card.querySelector('[data-updated]');
    fetch('https://api.github.com/repos/' + GH_USER + '/' + gh)
      .then(function (r) { if (!r.ok) throw new Error(r.status); return r.json(); })
      .then(function (d) {
        desc.textContent = d.description || '—';
        updated.textContent = fmtDate(d.pushed_at);
      })
      .catch(function () {
        desc.textContent = 'Could not reach GitHub — see the Docs or GitHub links below.';
      });
  });
})();
</script>
