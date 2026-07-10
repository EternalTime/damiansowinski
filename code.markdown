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
Open-source libraries behind the physics. 
Links to both the github <svg viewBox="0 0 16 16" width="14" height="14" style="fill: var(--black); vertical-align: -2px;" aria-hidden="true"><path d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.01 8.01 0 0016 8c0-4.42-3.58-8-8-8z"/></svg> repositories and the documentation <svg viewBox="0 0 24 24" width="14" height="14" style="fill: var(--black); vertical-align: -2px;" aria-hidden="true"><path d="M10 20v-6h4v6h5v-8h3L12 3 2 12h3v8z"/></svg> are in the drop down menus.
Feel free to fork any project and make it your own. 
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
  border: 2px solid var(--lower);
  border-radius: 8px;
  overflow: hidden;
  box-shadow:
    0 0 24px color-mix(in srgb, var(--upper) 60%, transparent),
    0 0 10px color-mix(in srgb, var(--upper) 45%, transparent);
}
.doc-spine { flex: 0 0 14px; background: var(--lower); }
.doc-body { flex: 1; min-width: 0; font-family: 'LatoDocs', 'Lato', 'Helvetica Neue', Arial, sans-serif; }
.doc-header { background: var(--upper); padding: 12px 20px; display: flex; align-items: center; gap: 12px; cursor: pointer; user-select: none; }
.doc-title { flex: 1; font-family: 'RobotoSlabDocs', 'Roboto Slab', Georgia, serif; font-size: 1.3em; font-weight: 700; color: var(--white); line-height: 1.2; margin: 0; }
.doc-toggle {
  flex: 0 0 auto;
  width: 0; height: 0;
  border-top: 7px solid transparent;
  border-bottom: 7px solid transparent;
  border-left: 10px solid var(--white);
  transition: transform 0.25s ease;
}
.doc-card.open .doc-toggle { transform: rotate(90deg); }
.doc-content { max-height: 0; overflow: hidden; transition: max-height 0.3s ease; }
.doc-content-inner { padding: 14px 20px 16px; }
.doc-desc { color: var(--text-body); font-size: 0.98em; line-height: 1.5; margin-bottom: 10px; }
.doc-updated { color: var(--text-muted); font-size: 0.8em; margin-bottom: 14px; }
.doc-btns { display: flex; gap: 10px; }
.doc-btns a {
  display: inline-flex;
  align-items: center;
  gap: 6px;
  text-decoration: none;
  font-family: 'LatoDocs', 'Lato', 'Helvetica Neue', Arial, sans-serif;
  font-size: 0.85em;
  color: var(--white) !important;
  text-shadow: none;
  background: var(--upper);
  border: 1px solid var(--upper);
  border-radius: 5px;
  padding: 5px 14px;
  transition: filter 0.15s;
}
.doc-btns a:hover { color: var(--white) !important; filter: brightness(1.15); }
.doc-btns a svg { display: block; }
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
    { repo: 'pyLEAFS',          label: 'Foraging Simulator',      docs: '/LEAFS/',  upper: '--pyleafs-primary', lower: '--pyleafs-nav', gh: 'LEAFS' },
    { repo: 'pyEDW',            label: 'Exo-Daisy World',         docs: '/pyEDW/',  upper: '--pyedw-primary',   lower: '--pyedw-nav' }
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
        '<div class="doc-header">' +
          '<div class="doc-title"></div>' +
          '<span class="doc-toggle" aria-hidden="true"></span>' +
        '</div>' +
        '<div class="doc-content">' +
          '<div class="doc-content-inner">' +
            '<div class="doc-desc" data-desc>Loading…</div>' +
            '<div class="doc-updated" data-updated></div>' +
            '<div class="doc-btns">' +
              '<a href="' + cfg.docs + '" target="_blank" rel="noopener" aria-label="Documentation">' +
                '<svg viewBox="0 0 24 24" width="15" height="15" fill="currentColor" aria-hidden="true"><path d="M10 20v-6h4v6h5v-8h3L12 3 2 12h3v8z"/></svg>' +
              '</a>' +
              '<a href="https://github.com/' + GH_USER + '/' + gh + '" target="_blank" rel="noopener" aria-label="GitHub">' +
                '<svg viewBox="0 0 16 16" width="15" height="15" fill="currentColor" aria-hidden="true"><path d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.01 8.01 0 0016 8c0-4.42-3.58-8-8-8z"/></svg>' +
              '</a>' +
            '</div>' +
          '</div>' +
        '</div>' +
      '</div>';
    card.querySelector('.doc-title').textContent = cfg.label;
    container.appendChild(card);

    const header = card.querySelector('.doc-header');
    const content = card.querySelector('.doc-content');
    header.setAttribute('role', 'button');
    header.setAttribute('tabindex', '0');
    header.setAttribute('aria-expanded', 'false');
    function toggle() {
      const open = card.classList.toggle('open');
      header.setAttribute('aria-expanded', open ? 'true' : 'false');
      content.style.maxHeight = open ? content.scrollHeight + 'px' : '0';
    }
    header.addEventListener('click', toggle);
    header.addEventListener('keydown', function (e) {
      if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); toggle(); }
    });

    const desc = card.querySelector('[data-desc]');
    const updated = card.querySelector('[data-updated]');
    fetch('https://api.github.com/repos/' + GH_USER + '/' + gh)
      .then(function (r) { if (!r.ok) throw new Error(r.status); return r.json(); })
      .then(function (d) {
        desc.textContent = d.description || '—';
        updated.textContent = fmtDate(d.pushed_at);
        if (card.classList.contains('open')) content.style.maxHeight = content.scrollHeight + 'px';
      })
      .catch(function () {
        desc.textContent = 'Could not reach GitHub — see the Docs or GitHub links below.';
        if (card.classList.contains('open')) content.style.maxHeight = content.scrollHeight + 'px';
      });
  });
})();
</script>
