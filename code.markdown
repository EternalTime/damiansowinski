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
Open-source libraries behind the physics. Each window pulls live metadata from
GitHub; the buttons open the docs and the source in a new tab.
</div><br>

<style>
.code-terminals {
  display: flex;
  flex-direction: column;
  gap: 34px;
  margin-top: 20px;
}
.term {
  --accent: var(--cyan);
  border: 2px solid var(--accent);
  border-radius: 10px;
  overflow: hidden;
  background: var(--bg-void);
  box-shadow: 0 0 18px rgba(0,0,0,0.55), 0 0 22px var(--accent);
  font-family: 'Source Code Pro', monospace;
}
.term-bar {
  display: flex;
  align-items: center;
  gap: 8px;
  padding: 8px 12px;
  background: var(--bg-panel);
  border-bottom: 1px solid var(--accent);
}
.term-dot {
  width: 12px; height: 12px;
  border-radius: 50%;
  display: inline-block;
}
.term-dot.r { background: var(--red); }
.term-dot.y { background: var(--amber); }
.term-dot.g { background: var(--green-light); }
.term-title {
  margin-left: 8px;
  color: var(--accent);
  font-size: 0.9em;
  letter-spacing: 0.03em;
}
.term-body {
  padding: 16px 18px 18px;
  color: var(--text-light);
  font-size: 0.92em;
  line-height: 1.55;
  min-height: 96px;
}
.term-prompt { color: var(--accent); }
.term-cmd    { color: var(--text-bright); }
.term-out    { color: var(--text-dim); white-space: pre-wrap; }
.term-key    { color: var(--accent); }
.term-caret::after {
  content: "▊";
  color: var(--accent);
  animation: term-blink 1s steps(1) infinite;
}
@keyframes term-blink { 50% { opacity: 0; } }
.term-btns {
  display: flex;
  gap: 12px;
  padding: 0 18px 18px;
}
.term-btns a {
  font-family: 'Source Code Pro', monospace;
  font-size: 0.86em;
  text-decoration: none;
  color: var(--accent);
  border: 1px solid var(--accent);
  border-radius: 6px;
  padding: 6px 14px;
  background: transparent;
  transition: background 0.15s, color 0.15s;
}
.term-btns a:hover {
  background: var(--accent);
  color: var(--bg-void);
}
</style>

<div class="code-terminals" id="code-terminals"></div>

<script>
(function () {
  const GH_USER = 'EternalTime';
  const REPOS = [
    { repo: 'pyCE',             docs: '/pyCE/',   accent: '--pyce-accent',   prompt: 'damian@cosmos:~/pyCE$' },
    { repo: 'CellularAutomata', docs: '/pyCA/',   accent: '--pyca-accent',   prompt: 'damian@lattice:~/pyCA$' },
    { repo: 'pyGD',             docs: '/pyGD/',   accent: '--pygd-accent',   prompt: 'damian@sync:~/pyGD$' },
    { repo: 'pyCoop',           docs: '/pyCoop/', accent: '--pycoop-accent', prompt: 'damian@game:~/pyCoop$' },
    { repo: 'pyLEAFS',          docs: '/LEAFS/',  accent: '--pyleafs-accent', prompt: 'damian@canopy:~/pyLEAFS$', gh: 'LEAFS' }
  ];

  const container = document.getElementById('code-terminals');

  function fmtDate(iso) {
    if (!iso) return 'unknown';
    const d = new Date(iso);
    return d.toLocaleDateString('en-US', { year: 'numeric', month: 'short', day: 'numeric' });
  }

  function typeOut(el, lines, done) {
    let li = 0;
    (function tick() {
      if (li >= lines.length) { el.classList.remove('term-caret'); if (done) done(); return; }
      el.insertAdjacentHTML('beforeend', '<span class="term-out">' + lines[li] + '</span>\n');
      li++;
      setTimeout(tick, 180);
    })();
  }

  REPOS.forEach(function (cfg) {
    const gh = cfg.gh || cfg.repo;
    const card = document.createElement('div');
    card.className = 'term';
    card.style.setProperty('--accent', 'var(' + cfg.accent + ')');
    card.innerHTML =
      '<div class="term-bar">' +
        '<span class="term-dot r"></span><span class="term-dot y"></span><span class="term-dot g"></span>' +
        '<span class="term-title">' + cfg.repo + '</span>' +
      '</div>' +
      '<div class="term-body">' +
        '<span class="term-prompt">' + cfg.prompt + '</span> ' +
        '<span class="term-cmd">cat info</span>\n' +
        '<span class="term-feed term-caret"></span>' +
      '</div>' +
      '<div class="term-btns">' +
        '<a href="' + cfg.docs + '" target="_blank" rel="noopener">./docs</a>' +
        '<a href="https://github.com/' + GH_USER + '/' + gh + '" target="_blank" rel="noopener">git remote</a>' +
      '</div>';
    container.appendChild(card);

    const feed = card.querySelector('.term-feed');
    fetch('https://api.github.com/repos/' + GH_USER + '/' + gh)
      .then(function (r) { if (!r.ok) throw new Error(r.status); return r.json(); })
      .then(function (d) {
        const lines = [
          '<span class="term-key">description</span> : ' + (d.description || '—'),
          '<span class="term-key">language</span>    : ' + (d.language || '—'),
          '<span class="term-key">updated</span>     : ' + fmtDate(d.pushed_at)
        ];
        typeOut(feed, lines);
      })
      .catch(function () {
        feed.classList.remove('term-caret');
        feed.innerHTML = '<span class="term-out" style="color:var(--red)">could not reach github — see ./docs or git remote below.</span>';
      });
  });
})();
</script>
