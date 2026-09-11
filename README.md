# GlissADe.jl

[![Stable Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://thealanjason.github.io/GlissADe.jl/stable)
[![Development documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://thealanjason.github.io/GlissADe.jl/dev)
[![Test workflow status](https://github.com/thealanjason/GlissADe.jl/actions/workflows/Test.yml/badge.svg?branch=main)](https://github.com/thealanjason/GlissADe.jl/actions/workflows/Test.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/thealanjason/GlissADe.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/thealanjason/GlissADe.jl)
[![Lint workflow Status](https://github.com/thealanjason/GlissADe.jl/actions/workflows/Lint.yml/badge.svg?branch=main)](https://github.com/thealanjason/GlissADe.jl/actions/workflows/Lint.yml?query=branch%3Amain)
[![Docs workflow Status](https://github.com/thealanjason/GlissADe.jl/actions/workflows/Docs.yml/badge.svg?branch=main)](https://github.com/thealanjason/GlissADe.jl/actions/workflows/Docs.yml?query=branch%3Amain)
[![JuliaCon26 Presentation](https://img.shields.io/badge/JuliaCon26-Presentation-9558B2)](https://thealanjason.github.io/GlissADe.jl/presentation/)
[![DOI](https://zenodo.org/badge/DOI/FIXME)](https://doi.org/FIXME)
[![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-2.1-4baaaa.svg)](CODE_OF_CONDUCT.md)
[![All Contributors](https://img.shields.io/github/all-contributors/thealanjason/GlissADe.jl?labelColor=5e1ec7&color=c0ffee&style=flat-square)](#contributors)
[![BestieTemplate](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/JuliaBesties/BestieTemplate.jl/main/docs/src/assets/badge.json)](https://github.com/JuliaBesties/BestieTemplate.jl)

GlissADe.jl is a differentiable simulator for surface flow over complex geometry, developed at the Chair of Methods for Model-based Development in Computational Engineering, RWTH Aachen University. The goal is an efficient simulator that is fully differentiable with respect to all model inputs, from low-dimensional parameters, such as constitutive model parameters, to high-dimensional fields, such as the underlying geometry.

## Getting Started

1. [Install Julia](https://julialang.org/downloads/)
2. Clone this repository
3. Change directory to the cloned folder
4. Setup the environment using `julia setup.jl`

## How to Cite

If you use GlissADe.jl in your work, please cite using the reference given in [CITATION.cff](https://github.com/thealanjason/GlissADe.jl/blob/main/CITATION.cff), or the BibTeX entry below.

```bibtex
@software{Tanish_Jain_and_Alan_Correa_GlissADe_jl,
  author = {Tanish Jain and Alan Correa},
  title = {{GlissADe.jl: Differentiable Simulator for Surface Flow over Complex Geometry}},
  year = {2024},
}
```

## Contributing

If you want to make contributions of any kind, please first take a look into our [contributing guide directly on GitHub](docs/src/90-contributing.md) or the [contributing page on the website](https://thealanjason.github.io/GlissADe.jl/dev/90-contributing/)

---

### Contributors

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->
<table>
  <tbody>
    <tr>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/thealanjason"><img src="https://avatars.githubusercontent.com/u/40885942?v=4?s=100" width="100px;" alt="Alan Correa"/><br /><sub><b>Alan Correa</b></sub></a><br /><a href="https://github.com/thealanjason/GlissADe.jl/commits?author=thealanjason" title="Code">💻</a> <a href="https://github.com/thealanjason/GlissADe.jl/commits?author=thealanjason" title="Documentation">📖</a> <a href="#infra-thealanjason" title="Infrastructure (Hosting, Build-Tools, etc)">🚇</a> <a href="#projectManagement-thealanjason" title="Project Management">📆</a> <a href="#maintenance-thealanjason" title="Maintenance">🚧</a> <a href="#ideas-thealanjason" title="Ideas, Planning, & Feedback">🤔</a></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/reckylurker"><img src="https://avatars.githubusercontent.com/u/102153509?v=4?s=100" width="100px;" alt="Tanish Jain"/><br /><sub><b>Tanish Jain</b></sub></a><br /><a href="https://github.com/thealanjason/GlissADe.jl/commits?author=reckylurker" title="Code">💻</a> <a href="https://github.com/thealanjason/GlissADe.jl/commits?author=reckylurker" title="Documentation">📖</a> <a href="https://github.com/thealanjason/GlissADe.jl/commits?author=reckylurker" title="Tests">⚠️</a></td>
    </tr>
  </tbody>
</table>

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->
