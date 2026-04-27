---
layout: default
title: Applets
---

<script src="{{ '/assets/js/applet-shell.js' | relative_url }}"></script>

<a href="{{ '/' | relative_url }}" style="float: right;">← Back to Home</a>
<br>

<img src="{{ '/assets/images/Under-Construction-Sign.png' | relative_url }}"
     style="float: center; width: 55%; margin: 0 20px 20px 20px; border-radius: 15px; border: 5px solid #ff6da2; filter: drop-shadow(2px 2px 2px rgba(0,0,0,0.9));"
     alt="Description"> 

---

<button class="applet-launch-btn" onclick="isingOpen(); this.blur()">Ising Model</button>

{% include applets/ising_applet.html %}

<button class="applet-launch-btn" onclick="vicsekOpen(); this.blur()">Vicsek Model</button>

{% include applets/vicsek_applet.html %}

<button class="applet-launch-btn" onclick="xyOpen(); this.blur()">XY Model</button>

{% include applets/xy_applet.html %}

<button class="applet-launch-btn" onclick="glOpen(); this.blur()">Ginzburg&ndash;Landau</button>

{% include applets/gl_applet.html %}

<button class="applet-launch-btn" onclick="gasOpen(); this.blur()">Hard Spheres</button>

{% include applets/gas_applet.html %}

<button class="applet-launch-btn" onclick="bmOpen(); this.blur()">Brownian Motion</button>

{% include applets/brownian_applet.html %}

<button class="applet-launch-btn" onclick="demonOpen(); this.blur()">Maxwell's Demon</button>

{% include applets/maxwell_demon_applet.html %}

<button class="applet-launch-btn" onclick="gsOpen(); this.blur()">Gray&ndash;Scott</button>

{% include applets/greyscott_applet.html %}

<button class="applet-launch-btn" onclick="kmOpen(); this.blur()">Kuramoto Model</button>

{% include applets/kuramoto_applet.html %}

<button class="applet-launch-btn" onclick="ttOpen(); this.blur()">Toner&ndash;Tu Model</button>

{% include applets/toner_tu_applet.html %}

<button class="applet-launch-btn" onclick="faOpen(); this.blur()">Forager</button>

{% include applets/forager_applet.html %}

<button class="applet-launch-btn" onclick="ltOpen(); this.blur()">Relativistic Lattice</button>

{% include applets/lattice_applet.html %}

<button class="applet-launch-btn" onclick="slOpen(); this.blur()">Spring Lattice</button>

{% include applets/spring_lattice_applet.html %}

<button class="applet-launch-btn" onclick="ltzOpen(); this.blur()">Lorentz Transformation</button>

{% include applets/lorentz_applet.html %}

<button class="applet-launch-btn" onclick="ohOpen(); this.blur()">Oscillon &mdash; FRW</button>

{% include applets/oscillon_hubble_applet.html %}

<button class="applet-launch-btn" onclick="shoOpen(); this.blur()">SHO</button>

{% include applets/sho_applet.html %}
