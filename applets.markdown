---
layout: default
title: Applets
---

<script src="{{ '/assets/js/applet-shell.js' | relative_url }}"></script>
<script src="{{ '/assets/js/applet-shell-mobile.js' | relative_url }}"></script>

<a href="{{ '/' | relative_url }}" style="float: right;">← Back to Home</a>
<br>


<h1 class="gradient-text2" style="font-size: 2em; text-align: left; filter: drop-shadow(1px 1px 0px rgba(0,0,0,0.9));">
  Computational Physics+ Applets
</h1>

<div>
Computational physics is awesome. It's the art of telling a computer how to simulate phenomena we see in the world. But it's not limited to just those phenomena we know about. We can use 
</div>

<img class="research-img-placeholder" src="{{ '/assets/images/Under-Construction-Sign.png' | relative_url }}" alt="Under Construction">

---

<h1 class="gradient-text2" style="font-size: 2em; text-align: left; filter: drop-shadow(1px 1px 0px rgba(0,0,0,0.9));">
  Statistical Physics
</h1>
<button class="applet-launch-btn" onclick="isingOpen(); this.blur()">Ising Model</button>
<button class="applet-launch-btn" onclick="xyOpen(); this.blur()">XY Model</button>
<button class="applet-launch-btn" onclick="glOpen(); this.blur()">Ginzburg&ndash;Landau</button>
<button class="applet-launch-btn" onclick="bmOpen(); this.blur()">Brownian Motion</button>
<button class="applet-launch-btn" onclick="gasOpen(); this.blur()">Hard Spheres</button>
<button class="applet-launch-btn" onclick="demonOpen(); this.blur()">Maxwell's Demon</button>
{% include applets/ising_applet.html %}
{% include applets/xy_applet.html %}
{% include applets/gl_applet.html %}
{% include applets/brownian_applet.html %}
{% include applets/gas_applet.html %}
{% include applets/maxwell_demon_applet.html %}

---

<h1 class="gradient-text2" style="font-size: 2em; text-align: left; filter: drop-shadow(1px 1px 0px rgba(0,0,0,0.9));">
  Pattern Formation
</h1>
<button class="applet-launch-btn" onclick="gsOpen(); this.blur()">Gray&ndash;Scott</button>
{% include applets/greyscott_applet.html %}

---

<h1 class="gradient-text2" style="font-size: 2em; text-align: left; filter: drop-shadow(1px 1px 0px rgba(0,0,0,0.9));">
  Ecological Physics
</h1>
<button class="applet-launch-btn" onclick="vicsekOpen(); this.blur()">Vicsek Model</button>
<button class="applet-launch-btn" onclick="ttOpen(); this.blur()">Toner&ndash;Tu Model</button>
<button class="applet-launch-btn" onclick="faOpen(); this.blur()">Forager</button>
<button class="applet-launch-btn" onclick="kmOpen(); this.blur()">Kuramoto Model</button>
{% include applets/vicsek_applet.html %}
{% include applets/forager_applet.html %}
{% include applets/toner_tu_applet.html %}
{% include applets/kuramoto_applet.html %}

---

<h1 class="gradient-text2" style="font-size: 2em; text-align: left; filter: drop-shadow(1px 1px 0px rgba(0,0,0,0.9));">
  Special Relativity
</h1>
<button class="applet-launch-btn" onclick="ltOpen(); this.blur()">Relativistic Lattice</button>
<button class="applet-launch-btn" onclick="ltzOpen(); this.blur()">Lorentz Transformation</button>
{% include applets/lattice_applet.html %}
{% include applets/lorentz_applet.html %}

---

<h1 class="gradient-text2" style="font-size: 2em; text-align: left; filter: drop-shadow(1px 1px 0px rgba(0,0,0,0.9));">
  Continuum & Classical Mechanics
</h1>
<button class="applet-launch-btn" onclick="slOpen(); this.blur()">Spring Lattice</button>
<button class="applet-launch-btn" onclick="shoOpen(); this.blur()">SHO</button>
{% include applets/spring_lattice_applet.html %}
{% include applets/sho_applet.html %}

---

<h1 class="gradient-text2" style="font-size: 2em; text-align: left; filter: drop-shadow(1px 1px 0px rgba(0,0,0,0.9));">
  Quantum Mechanics
</h1>

<button class="applet-launch-btn" onclick="iswOpen(); this.blur()">Infinite Square Well</button>
<button class="applet-launch-btn" onclick="qhoOpen(); this.blur()">Quantum Harmonic Oscillator</button>
{% include applets/isw_applet.html %}
{% include applets/qho_applet.html %}

---

<h1 class="gradient-text2" style="font-size: 2em; text-align: left; filter: drop-shadow(1px 1px 0px rgba(0,0,0,0.9));">
  Field Theory
</h1>
<button class="applet-launch-btn" onclick="ohOpen(); this.blur()">Oscillon Formation</button>
{% include applets/oscillon_hubble_applet.html %}

{% comment %}
<button class="applet-launch-btn" onclick="qtunOpen(); this.blur()">Quantum Tunneling</button>

{% include qtun_applet.html %}
{% endcomment %}

