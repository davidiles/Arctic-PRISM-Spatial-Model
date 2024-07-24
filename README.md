# Arctic PRISM Spatial Model

 First draft of a Bayesian spatial analysis of shorebird abundance for Arctic PRISM.  Goal is to develop a develop a model-based analysis of shorebird density (and change over time), confirm it works as intended through spatially explicit simulations, apply the model to empirical data, and compare to design-based analysis.
 
 
## Survey locations

![PRISM survey locations](https://github.com/davidiles/Arctic-PRISM-Spatial-Model/blob/main/output/PRISM_survey_locations.png)
 
 
## Proposed model structure

The Bayesian model has several key features, but most of these can be changed/adjusted easily:

- Negative binomial or overdispersed Poisson error
    - Can also include zero-inflation if necessary
- Spatial autocorrelation process
    - Modeled using a Matern covariance structure / Gaussian Markov random field (widely used, easily implemented, spatially continuous)
- Plot-level random effects (to account for the fact that many plots are repeatedly surveyed)
- Survey-level random effect (to account for overdispersion)
- Habitat covariates (extracted from relevant raster layers)
- PRISM strata can be included as discrete covariates
- Integrated analysis of "rapid surveys" and "intensive plots" which allows for detectability corrections to be embedded in the model and uncertainty to be propagated
- Informative priors where required to help constrain estimates, improve precision, improve parameter identifiability, etc

Models will be fit using Bayesian methods, [using INLA](https://www.r-inla.org/)

## Simulation

The first step of this analysis is to conduct simulations to confirm the Bayesian model (and if possible, the design-based analysis) produces unbiased estimates with appropriate credible interval coverage (i.e., the true population total should be within the 90% credible interval for 90% of simulations).

The approach entails several steps:
1) Download eBird species distribution rasters to use as "truth" in simulations.  These are just used as examples of how species *might be* distributed, and we want to confirm that if the species is distributed in *any* manner, the model can accurately estimate the spatial pattern of abundance as well as the total abundance (i.e., the sum of all pixels)
2) Sample from the "true" rasters at existing PRISM survey locations.  These represent the observed data.  Also include some degree of observation error.
3) Analyze simulated data using Bayesian model and design-based approach.  Evaluate bias in estimates and credible interval coverage.  Repeat this for many species.

- Script "PRISM_Simulation.R" conducts simulations for 46 species.  Results for each species are contained in the [output sub-folder](https://github.com/davidiles/Arctic-PRISM-Spatial-Model/tree/main/output). 
- [An example of a simulation for Horned Lark can be found here](https://github.com/davidiles/Arctic-PRISM-Spatial-Model/blob/main/output/simulation_Horned%20Lark_nbinomial.png).
- [Summary of results across 46 species can be viewed here](https://github.com/davidiles/Arctic-PRISM-Spatial-Model/blob/main/output/simulation_species_estimates_nbinomial.png).

For estimates of overall population size, credible interval coverage is slightly lower than nominal.  The 90% credible interval contains the true population total for *76%* of species.
On average, bias in total population size is about 4% (so very small).  56% of species have estimates that are negatively biased, which is excellent (a perfect estimator would underestimate 50% of the time).

- need to replicate this analysis using the design-based approach.

## Empirical analysis

I applied the Bayesian model to empirical data for as many species as possible.  Results are contained in the [output sub-folder](https://github.com/davidiles/Arctic-PRISM-Spatial-Model/tree/main/output), where .

An example of model results for SESA is shown below:
![SESA](https://github.com/davidiles/Arctic-PRISM-Spatial-Model/blob/main/output/empirical_SESA_nbinomial.png)

Across all species for which I fit models, here is a comparison of estimated population sizes (whiskers are 90% credible intervals).  Black dots are the point estimates from the design-based analysis reported in Paul Smith's document.

![Population estimates](https://github.com/davidiles/Arctic-PRISM-Spatial-Model/blob/main/output/species_estimates.png)

## Issues / Ideas

1) Need to obtain the original sampling frame (grid across the survey area). Currently, I am fitting the model to an eye-balled portion of BCR 3, which means I am probably predicting into areas that were not part of the original design.  This would reduce precision of the estimates (and also result in over-estimated population size, compared to a design-based analysis of a smaller subset of the region).

2) Need to obtain covariate layers that are appropriate for Arctic Canada.  Currently using a sub-optimal covariate layer (land cover of Canada 2020).  Ideally, we should start by using identical covariate definitions and stratum categories as in the original PRISM design.

3) Precision of Bayesian estimates is generally lower (and largely uncorrelated) with precision reported in Smith et al.'s design-based analysis.  Need to confirm that design-based precision is accurate (through simulation), and ensure same data and covariates/strata are being used.

4) Reconnaissance surveys on spaghetti transects - do we have an effort measurement to associate with these?  If so, can include in analysis.  Can also consider integrating additional sources of information, if any exist (e.g., from eBird?).  However, this is a rabbit hole that would need to be carefully considered, so likely want to avoid this in the near-term.
 
