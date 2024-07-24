# Arctic PRISM Spatial Model

 First draft of a Bayesian spatial analysis of shorebird abundance for Arctic PRISM.  Goal is to develop a compare design-based analysis to model-based analysis.
 
 
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

- Script "PRISM_1_Simulation.R" conducts simulations, confirms that confidence intervals are appropriate (95% credible intervals contain the true population total for 95% of species).  Bias is minimal.

## Empirical analysis

Currently fitting simple model for each species.  Data is from rapid surveys.



Model output:
- Maps of relative density
- Estimates of total abundance

![SESA](https://github.com/davidiles/Arctic-PRISM-Spatial-Model/blob/main/output/empirical_SESA_nbinomial.png)

Estimates of population size from Bayesian analysis (blue), compared to estimates reported in Smith et al.'s design-based analysis.
![Population estimates](https://github.com/davidiles/Arctic-PRISM-Spatial-Model/blob/main/output/species_estimates.png)

*Note: Also eventually analyze the same data using design-based analysis* 

## Issues / Ideas

1) Currently only using Rapid Surveys, and removing all observations outside plot.
 - Include Intensive Surveys in same analysis (improve precision?).  Set prior of "intensive survey effect" to 1.13 as per Paul's previous analysis (mean across species), but only fit effect for each species separately.
 - Reconnaissance surveys on spaghetti transects - do we have an effort measurement to associate with these?  If so, can include in analysis.
 

2) Need to use same covariate definitions as in Paul's original analysis, along with same strata that were used to select survey locations.

3) Precision of Bayesian estimates is generally lower (and largely uncorrelated) with precision reported in Smith et al.'s design-based analysis.  Need to confirm that design-based precision is accurate (through simulation), and ensure same data and covariates/strata are being used.

4) Evaluate change analysis through simulation.  Start by assuming same plots are revisited.  2 levels of spatial autocorrelation in change pattern (highly patchy vs highly continuous).  For each simulation scenario, evaluate correlation (x-y plot) between true change and estimate change of overall population.

5) Random year effects?  These will be confounded with space, but may be somewhat (?) estimable since many sites have been visited twice?  Or perhaps use 5-year windows?



- No LESA detected in intensive plots?