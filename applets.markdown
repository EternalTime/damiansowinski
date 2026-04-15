---
layout: default
title: Applets
---


<button class="applet-launch-btn" onclick="isingOpen()">Ising Model</button>

{% include ising_applet.html %}

<button class="applet-launch-btn" onclick="vicsekOpen()">Vicsek Model</button>

{% include vicsek_applet.html %}

<button class="applet-launch-btn" onclick="xyOpen()">XY Model</button>

{% include xy_applet.html %}

<button class="applet-launch-btn" onclick="glOpen()">Ginzburg&ndash;Landau</button>

{% include gl_applet.html %}

<button class="applet-launch-btn" onclick="gasOpen()">Hard Spheres</button>

{% include gas_applet.html %}

<button class="applet-launch-btn" onclick="demonOpen()">Maxwell's Demon</button>

{% include maxwell_demon_applet.html %}
