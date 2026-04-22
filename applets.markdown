---
layout: default
title: Applets
---

<script src="{{ '/assets/js/applet-shell.js' | relative_url }}"></script>

<a href="{{ '/' | relative_url }}" style="float: right;">← Back to Home</a>
<br>

<button class="applet-launch-btn" onclick="isingOpen(); this.blur()">Ising Model</button>

{% include ising_applet.html %}

<button class="applet-launch-btn" onclick="vicsekOpen(); this.blur()">Vicsek Model</button>

{% include vicsek_applet.html %}

<button class="applet-launch-btn" onclick="xyOpen(); this.blur()">XY Model</button>

{% include xy_applet.html %}

<button class="applet-launch-btn" onclick="glOpen(); this.blur()">Ginzburg&ndash;Landau</button>

{% include gl_applet.html %}

<button class="applet-launch-btn" onclick="gasOpen(); this.blur()">Hard Spheres</button>

{% include gas_applet.html %}

<button class="applet-launch-btn" onclick="bmOpen(); this.blur()">Brownian Motion</button>

{% include brownian_applet.html %}

<button class="applet-launch-btn" onclick="demonOpen(); this.blur()">Maxwell's Demon</button>

{% include maxwell_demon_applet.html %}

<button class="applet-launch-btn" onclick="gsOpen(); this.blur()">Gray&ndash;Scott</button>

{% include greyscott_applet.html %}

<button class="applet-launch-btn" onclick="kmOpen(); this.blur()">Kuramoto Model</button>

{% include kuramoto_applet.html %}

<button class="applet-launch-btn" onclick="ttOpen(); this.blur()">Toner&ndash;Tu Model</button>

{% include toner_tu_applet.html %}

<button class="applet-launch-btn" onclick="faOpen(); this.blur()">Forager</button>

{% include forager_applet.html %}

<button class="applet-launch-btn" onclick="ltOpen(); this.blur()">Relativistic Lattice</button>

{% include lattice_applet.html %}

<button class="applet-launch-btn" onclick="slOpen(); this.blur()">Spring Lattice</button>

{% include spring_lattice_applet.html %}

<button class="applet-launch-btn" onclick="ltzOpen(); this.blur()">Lorentz Transformation</button>

{% include lorentz_applet.html %}

<button class="applet-launch-btn" onclick="ohOpen(); this.blur()">Oscillon &mdash; FRW</button>

{% include oscillon_hubble_applet.html %}

<button class="applet-launch-btn" onclick="shoOpen(); this.blur()">SHO</button>

{% include sho_applet.html %}
