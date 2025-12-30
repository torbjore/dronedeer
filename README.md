# **Supporting material for** 'Estimating red deer (*Cervus elaphus*) population density using drones in a steep and rugged terrain' (DOI: 10.1002/wlb3.01533)


This repository contains all data and computer code used for the paper, in addition to the [Supporting Information Document (html)](https://torbjore.github.io/dronedeer/).


## Directories

- `R`: R code for running models and some utility functions
  - `R/nimble_models`: NIMBLE models (sourced in from R code)
- `data`: Data used by R code and scripts for formatting data for the NIMBLE models. See [`SiteDataInfo.html`](https://torbjore.github.io/dronedeer/SiteDataInfo.html) for description of the main data file.
- `posterior_samples`: Posterior samples form NIMBLE models. Only posterior samples from the `gamma_2` model are included in the online repository. Run models in the R directory to generate posterior samples from other models.
- `Supplemental_Info`: Quarto source code for [Supporting Information Document (html)](https://torbjore.github.io/dronedeer/).
