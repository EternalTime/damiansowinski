---
layout: default
title: Publications
permalink: /publications/
---

<h1>Publications</h1>

<div id="publications-list"></div>

<script src="https://cdn.jsdelivr.net/npm/citation-js"></script>

<script>
const Cite = require('citation-js');

const bibtex = `

@article{sowinski2025exo,
  title     = {Exo-Daisy World: Revisiting Gaia theory through an informational architecture perspective},
  author    = {Sowinski, Damian R and Ghoshal, Gourab and Frank, Adam},
  journal   = {The Planetary Science Journal},
  volume    = {6},
  number    = {7},
  pages     = {176},
  year      = {2025},
  publisher = {The American Astronomical Society},
  url       = {/damiansowinski/assets/docs/ExoDaisy_Sowinski_2025.pdf}
}

@inproceedings{sowinski2025causal,
  title     = {Causal Leverage Density: A Universal Framework for Semantic Information},
  author    = {Sowinski, Damian R and Holler, Silvia and Wong, Mike and Segal, Gary and Parkinson, David and Vidal, Clément and Henderson Cleaves, James and Sinhadc, Pritvik and Prabhu, Anirudh and Bartlett, Stuart},
  booktitle = {Artificial Life Conference Proceedings 37},
  volume    = {2025},
  number    = {1},
  pages     = {43},
  year      = {2025},
  organization = {MIT Press},
  url       = {/damiansowinski/assets/docs/CLD_Sowinski_2025.pdf}
}

@article{sowinski2025configurational,
  title     = {Configurational information measures, phase transitions, and an upper bound on complexity},
  author    = {Sowinski, Damian R and Kelty, Sean and Ghoshal, Gourab},
  journal   = {arXiv preprint arXiv:2503.02980},
  year      = {2025},
  url       = {/damiansowinski/assets/docs/Ising_Sowinski_2025.pdf}
}

@article{sowinski2024information,
  title     = {Information-theoretic description of a feedback-control Kuramoto model},
  author    = {Sowinski, Damian R and Frank, Adam and Ghoshal, Gourab},
  journal   = {Physical Review Research},
  volume    = {6},
  number    = {4},
  pages     = {043188},
  year      = {2024},
  publisher = {APS},
  url       = {/damiansowinski/assets/docs/Kuramoto_Sowinski_2024.pdf}
}

@article{sowinski2023semantic,
  title     = {Semantic information in a model of resource gathering agents},
  author    = {Sowinski, Damian R and Carroll-Nellenback, Jonathan and Markwick, Robert N and Pinero, Jordi and Gleiser, Marcelo and Kolchinsky, Artemy and Ghoshal, Gourab and Frank, Adam},
  journal   = {PRX Life},
  volume    = {1},
  number    = {2},
  pages     = {023003},
  year      = {2023},
  publisher = {APS},
  url       = {/damiansowinski/assets/docs/Foraging_Sowinski_2023.pdf}
}

@article{sowinski2022consensus,
  title     = {The consensus problem in polities of agents with dissimilar cognitive architectures},
  author    = {Sowinski, Damian R. and Carroll-Nellenback, Jonathan and DeSilva, Jeremy and Frank, Adam and Ghoshal, Gourab and Gleiser, Marcelo},
  journal   = {Entropy},
  volume    = {24},
  number    = {10},
  pages     = {1378},
  year      = {2022},
  publisher = {MDPI},
  url       = {/damiansowinski/assets/docs/Consensus_Sowinski_2022.pdf}
}

@article{sowinski2021poroelasticity,
  title     = {Poroelasticity as a model of soft tissue structure: Hydraulic permeability reconstruction for magnetic resonance elastography in silico},
  author    = {Sowinski, Damian R and McGarry, Matthew DJ and Van Houten, Elijah EW and Gordon-Wylie, Scott and Weaver, John B and Paulsen, Keith D},
  journal   = {Frontiers in physics},
  volume    = {8},
  pages     = {617582},
  year      = {2021},
  publisher = {Frontiers Media SA},
  url       = {/damiansowinski/assets/docs/Poro_Sowinski_2020.pdf}
}

@article{gleiser2018configurational,
  title     = {Configurational information approach to instantons and false vacuum decay in D-dimensional spacetime},
  author    = {Sowinski, Damian and Gleiser, Marcelo},
  journal   = {Physical Review D},
  volume    = {98},
  number    = {5},
  pages     = {056026},
  year      = {2018},
  publisher = {APS},
  url       = {/damiansowinski/assets/docs/Instantons_Sowinski_2018.pdf}
}

@article{sowinski2017information,
  title     = {Information dynamics at a phase transition},
  author    = {Sowinski, Damian and Gleiser, Marcelo},
  journal   = {Journal of Statistical Physics},
  volume    = {167},
  number    = {5},
  pages     = {1221--1232},
  year      = {2017},
  publisher = {Springer},
  url       = {/damiansowinski/assets/docs/PhaseDynamics_Sowinski_2017.pdf}
}

@article{sowinski2016complexity,
  title     = {Complexity and Stability for Epistemic Agents: The Foundations and Phenomenology of Configurational Entropy},
  author    = {Sowinski, Damian Radoslaw},
  journal   = {Ph. D. Thesis},
  year      = {2016}
}

@article{gleiser2013information,
  title     = {Information-entropic stability bound for compact objects: Application to Q-balls and the Chandrasekhar limit of polytropes},
  author    = {Gleiser, Marcelo and Sowinski, Damian},
  journal   = {Physics Letters B},
  volume    = {727},
  number    = {1-3},
  pages     = {272--275},
  year      = {2013},
  publisher = {Elsevier},
  url       = {/damiansowinski/assets/docs/Polytropes_Sowinski_2013.pdf}
}

@article{pinero2026information,
  title     = {Information bounds production in replicator systems},
  author    = {Pinero, Jordi and Sowinski, Damian R and Ghoshal, Gourab and Frank, Adam and Kolchinsky, Artemy},
  journal   = {Communications Physics},
  year      = {2026},
  publisher = {Nature Publishing Group UK London}
}

@article{quillen2025notions,
  title     = {Notions of adiabatic drift in the quantized Harper model},
  author    = {Quillen, Alice C and Skerrett, Nathan and Sowinski, Damian R and Miakhel, Abobakar Sediq},
  journal   = {Physical Review A},
  volume    = {112},
  number    = {4},
  pages     = {042226},
  year      = {2025},
  publisher = {APS}
}

@article{mcgarry2025vivo,
  title     = {In vivo magnetic resonance imaging of the interstitial pressure gradients (pgMRI) using a pulsatile poroelastic computational model},
  author    = {McGarry, Matthew and Sowinski, Damian and Tan, Likun and Weaver, John and Zwanenburg, Jacobus JM and Paulsen, Keith},
  journal   = {Interface focus},
  volume    = {15},
  number    = {1},
  year      = {2025},
  publisher = {The Royal Society}
}

@article{mcgarry2022mapping,
  title     = {Mapping heterogenous anisotropic tissue mechanical properties with transverse isotropic nonlinear inversion MR elastography},
  author    = {McGarry, Matthew and Van Houten, Elijah and Sowinski, Damian and Jyoti, Dhrubo and Smith, Daniel R and Caban-Rivera, Diego A and McIlvain, Grace and Bayly, Philip and Johnson, Curtis L and Weaver, John and others},
  journal   = {Medical Image Analysis},
  volume    = {78},
  year      = {2022},
  publisher = {Elsevier}
}

@article{mcgarry2021heterogenous,
  title     = {A heterogenous, time harmonic, nearly incompressible transverse isotropic finite element brain simulation platform for MR elastography},
  author    = {McGarry, Matthew and Houten, Elijah Van and Guertler, Charlotte and Okamoto, Ruth and Smith, Daniel and Sowinski, Damian and Johnson, Curtis and Bayly, Philip and Weaver, John and Paulsen, Keith},
  journal   = {Physics in Medicine \& Biology},
  volume    = {66},
  number    = {5},
  pages     = {055029},
  year      = {2021},
  publisher = {IOP Publishing}
}

@article{hannum2022correlated,
  title     = {Correlated noise in brain magnetic resonance elastography},
  author    = {Hannum, Ariel J and McIlvain, Grace and Sowinski, Damian and McGarry, Matthew DJ and Johnson, Curtis L},
  journal   = {Magnetic resonance in medicine},
  volume    = {87},
  number    = {3},
  pages     = {1313--1328},
  year      = {2022},
  publisher = {Wiley Online Library}
}

@article{jyoti2022quantifying,
  title     = {Quantifying stability of parameter estimates for in vivo nearly incompressible transversely-isotropic brain MR elastography},
  author    = {Jyoti, Dhrubo and McGarry, Matthew and Van Houten, Elijah and Sowinski, Damian and Bayly, Philip V and Johnson, Curtis L and Paulsen, Keith},
  journal   = {Biomedical physics and engineering express},
  volume    = {8},
  number    = {3},
  pages     = {035015},
  year      = {2022},
  publisher = {IOP Publishing}
}

@article{bowen2022visual,
  title     = {Visual stream connectivity predicts assessments of image quality},
  author    = {Bowen, Elijah FW and Rodriguez, Antonio M and Sowinski, Damian R and Granger, Richard},
  journal   = {Journal of vision},
  volume    = {22},
  number    = {11},
  pages     = {4--4},
  year      = {2022},
  publisher = {The Association for Research in Vision and Ophthalmology}
}

@article{gleiser2018configurational,
  title     = {Configurational entropy as a lifetime predictor and pattern discriminator for oscillons},
  author    = {Gleiser, Marcelo and Stephens, Michelle and Sowinski, Damian},
  journal   = {Physical Review D},
  volume    = {97},
  number    = {9},
  pages     = {096007},
  year      = {2018},
  publisher = {APS}
}

@incollection{key,
  author    = "Marcelo Gleiser and Damian Sowinski",
  title     = "How We Make Sense of the World: Information, Map-Making, and The Scientiﬁc Narrative",
  booktitle = "The map and the territory: Exploring the foundations of science, thought and reality",
  publisher = "Springer",
  year      = "2018",
  editor    = "Wuppuluri, Shyam and Doria, Francisco Antonio",
  pages     = "100--120"
}

@article{gleiser2015information,
  title     = {Information-entropic signature of the critical point},
  author    = {Gleiser, Marcelo and Sowinski, Damian},
  journal   = {Physics Letters B},
  volume    = {747},
  pages     = {125--128},
  year      = {2015},
  publisher = {Elsevier}
}

@article{gleiser2013information,
  title     = {Information-entropic stability bound for compact objects: Application to Q-balls and the Chandrasekhar limit of polytropes},
  author    = {Gleiser, Marcelo and Sowinski, Damian},
  journal   = {Physics Letters B},
  volume    = {727},
  number    = {1-3},
  pages     = {272--275},
  year      = {2013},
  publisher = {Elsevier}
}

`;

const cite = new Cite(bibtex);
const entries = cite.data;
let output = '';
entries.forEach(entry => {
  const single = new Cite(entry);
  output += single.format('bibliography', {
    format: 'html',
    template: 'acm',
    lang: 'en-US'
  });
});

output = output.replace(
  /(https?:\/\/[^\s<>"]+|\/damiansowinski\/assets\/docs\/[^\s<>"]+)/g,
  '<a href="$1" target="_blank" style="font-size: 0.7em; font-weight: bold;">[PDF]</a>'
);

document.getElementById('publications-list').innerHTML = output;
</script>

<style>
  .csl-entry {
    margin-bottom: 1.2em;
    line-height: 1.1;
  }
</style>
