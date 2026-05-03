---
layout: default
title: Applets
---

<script src="{{ '/assets/js/applet-shell.js' | relative_url }}"></script>
<script src="{{ '/assets/js/applet-shell-mobile.js' | relative_url }}"></script>
{% for f in site.static_files %}{% if f.path contains '/assets/js/applets/' and f.extname == '.js' %}<script src="{{ f.path | relative_url }}" defer></script>
{% endif %}{% endfor %}

<style>
.wip-wrap { display: inline-flex; flex-direction: column; align-items: center; gap: 3px; }
.wip-wrap .applet-launch-btn { outline: 2px solid var(--amber); outline-offset: 2px; box-shadow: 0px 1px 1px rgba(0,0,0,0.95); }
.wip-label { font-size: 0.6em; color: var(--amber); letter-spacing: 0.04em; line-height: 1; text-shadow: 0px 1px 1px rgba(0,0,0,0.95); position: relative; z-index: 10; }
</style>

<a href="{{ '/' | relative_url }}" style="float: right;">← Back to Home</a>
<br>


<h1 class="gradient-text2" style="font-size: 2em; text-align: left; filter: drop-shadow(1px 1px 0px rgba(0,0,0,0.9));">
  Computational Physics+ Applets
</h1>

<div>
Computational physics is awesome. 
It's the art of telling a computer how to simulate phenomena we see in the world. 
But it's not limited to just those phenomena we know about. 
We can use explore unknown worlds obeying strange rules, limited solely by our imaginations. 
</div><br>

---

<h2 class="gradient-text2" style="font-size: 1.5em;">Statistical Physics</h2>
<button class="applet-launch-btn" onclick="isingOpen(); this.blur()">Ising Model</button>
<button class="applet-launch-btn" onclick="xyOpen(); this.blur()">XY Model</button>
<button class="applet-launch-btn" onclick="glOpen(); this.blur()">Ginzburg&ndash;Landau</button>
<button class="applet-launch-btn" onclick="bmOpen(); this.blur()">Brownian Motion</button>
<span class="wip-wrap"><button class="applet-launch-btn" onclick="gasOpen(); this.blur()">Hard Spheres</button><span class="wip-label">Under Construction</span></span>
<button class="applet-launch-btn" onclick="demonOpen(); this.blur()">Maxwell's Demon</button>

---

<h2 class="gradient-text2" style="font-size: 1.5em;">Pattern Formation</h2>
<button class="applet-launch-btn" onclick="gsOpen(); this.blur()">Gray&ndash;Scott</button>

---

<h2 class="gradient-text2" style="font-size: 1.5em;">Ecological Physics</h2>
<button class="applet-launch-btn" onclick="vicsekOpen(); this.blur()">Vicsek Model</button>
<span class="wip-wrap"><button class="applet-launch-btn" onclick="ttOpen(); this.blur()">Toner&ndash;Tu Model</button><span class="wip-label">Under Construction</span></span>
<button class="applet-launch-btn" onclick="faOpen(); this.blur()">Forager</button>
<button class="applet-launch-btn" onclick="fa3Open(); this.blur()">Forager 3D</button>
<span class="wip-wrap"><button class="applet-launch-btn" onclick="antsOpen(); this.blur()">Ant Colony</button><span class="wip-label">Under Construction</span></span>
<button class="applet-launch-btn" onclick="kmOpen(); this.blur()">Kuramoto Model</button>

---

<h2 class="gradient-text2" style="font-size: 1.5em;">Special Relativity</h2>
<button class="applet-launch-btn" onclick="ltOpen(); this.blur()">Relativistic Lattice</button>
<button class="applet-launch-btn" onclick="ltzOpen(); this.blur()">Lorentz Transformation</button>

---

<h2 class="gradient-text2" style="font-size: 1.5em;">Continuum &amp; Classical Mechanics</h2>
<button class="applet-launch-btn" onclick="slOpen(); this.blur()">Spring Lattice</button>
<button class="applet-launch-btn" onclick="shoOpen(); this.blur()">SHO</button>
<span class="wip-wrap"><button class="applet-launch-btn" onclick="femOpen(); this.blur()">Elastic Cylinder</button><span class="wip-label">Under Construction</span></span>
<button class="applet-launch-btn" onclick="tbOpen(); this.blur()">Three-Body Problem</button>

---

<h2 class="gradient-text2" style="font-size: 1.5em;">Quantum Mechanics</h2>
<button class="applet-launch-btn" onclick="iswOpen(); this.blur()">Infinite Square Well</button>
<button class="applet-launch-btn" onclick="qhoOpen(); this.blur()">Quantum Harmonic Oscillator</button>
<button class="applet-launch-btn" onclick="qtunOpen(); this.blur()">Quantum Tunneling</button>

---

<h2 class="gradient-text2" style="font-size: 1.5em;">Field Theory</h2>
<button class="applet-launch-btn" onclick="ohOpen(); this.blur()">Oscillon Formation</button>
