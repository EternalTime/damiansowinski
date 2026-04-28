---
layout: default
mathjax: true
---

<h1 class="gradient-text home-title">
  <strong>Damian R. Sowinski</strong>
</h1>
<div class="gradient-text2" style="font-size: 2em; text-align: center; filter: drop-shadow(1px 1px 0px rgba(0,0,0,0.9));">
Theoretical Physicist & Applied Mathemagician
</div>

<div class="home-photo-wrap">
  <img class="home-photo-img" src="{{ '/assets/images/Damian_in_Key_West.jpg' | relative_url }}" alt="Damian Sowinski">
</div>

<div class="gradient-text2 home-subtitle" style="text-align: center; filter: drop-shadow(1px 1px 0px rgba(0,0,0,0.9));">
Artificial Life & Intelligence | Complexity | Information | Simulation |Visualization
</div>
<br>

<div class="home-videos">
  <video class="home-video" autoplay muted loop playsinline>
    <source src="{{ '/assets/videos/flockMF_0p45.mp4' | relative_url }}" type="video/mp4">
  </video>
  <video class="home-video" autoplay muted loop playsinline>
    <source src="{{ '/assets/videos/flockMF_0p15.mp4' | relative_url }}" type="video/mp4">
  </video>
  <video class="home-video" autoplay muted loop playsinline>
    <source src="{{ '/assets/videos/flockMF_0p75.mp4' | relative_url }}" type="video/mp4">
  </video>
</div>

<script>
  window.addEventListener('load', function() {
    document.querySelectorAll('video').forEach(function(v) {
      v.loop = true;
      v.muted = true;
      v.play().catch(function() {});
      v.addEventListener('ended', function() { v.currentTime = 0; v.play(); });
      v.addEventListener('pause', function() { if (!document.hidden) { v.currentTime = 0; v.play(); } });
    });
  });
</script>

<h3 class="gradient-text2" style="font-size: 1.5em; text-align: left; filter: drop-shadow(1px 1px 0px rgba(0,0,0,0.9));">
  What is Agency?
</h3>

Our description of the physical stuff of the world, in the language of modern physics, paints a picture of Nature driven by extremum principles.
Energy minimization and entropy maximization guide this soup of particles, the *physical substrate*, towards thermodynamic equilibrium where, once there, nothing <span title="...at least not on any relevant timescales. See Freeman Dyson's work on eschatology.">interesting ever happens</span>.
Yet, in this inexorable slide towards inevitable boredom, matter organizes itself from atoms to crystals and catalytic networks, patterns forming and growing in complexity and scale.
And, if the conditions are just right, macroscopic structure supervenes over microscopic details giving birth to agency within the substrate.
Matter, once subject solely to the whims of energy and entropy, gains the ability to scream **No! Not today Equilibrium!**, to act to suppress its slide towards that inevitable end.
Our description of the physical stuff of the world, in the language of modern physics, does not know how matter does that.

As AI gets more sophisticated, the ethical question of whether it is alive will grow in importance over the coming decades.
My research is a cross-disciplinary approach towards trying to understand how agency emerges from the physical substrate of the world, a question of profound importance to humanity's understanding of itself, artificial intelligence, the biosphere both are a part of.
Living systems act teleologically; since prediction is the providence of information theory, my approach towards tackling this enormous question draws not just from the physical and life sciences, but computer science and machine learning, neuroscience, the philosophy of language, and epistemology.
A large component of this program is using the techniques of information theory to create informational narratives of systems complementing their traditional physical descriptions.
This means characterizing how information is generated, stored, and flows in biological, ecological, and sociological systems, and identifying relationships between subsystems that lead to these informational structures.



---

## Links

<div style="text-align: center; margin-bottom: 40px;" markdown="1">

[Email](mailto:DRSowinski@gmail.com) | [Curriculum Vitae]({{ '/assets/docs/Sowinski_CV.pdf' | relative_url }})|  [GitHub](https://github.com/eternaltime) | [LinkedIn](https://www.linkedin.com/in/damian-sowinski-413b4214)

</div>
