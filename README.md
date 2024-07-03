# Arctic PRISM Spatial Model

 Bayesian spatial analysis of shorebird abundance for Arctic PRISM.  Goal is to compare design-based analysis to model-based analysis.
 
![PRISM survey locations](https://github.com/davidiles/Arctic-PRISM-Spatial-Model/blob/main/output/PRISM_survey_locations.png)
 
## Model

Currently using [Zero-Adjusted Poisson models](https://inlabru-org.github.io/inlabru/articles/zip_zap_models.html) because there are large portions of the study area that are not suitable, which could be better modeled by separating the zero process from the non-zero process.

Model components:

- spatial random fields for the zero (presence/absence) and non-zero (counts) process.
- covariate effects (currently Ecozones as categorical variables; i.e., strata)
- barrier effects caused by large waterbodies.

I also defined a (fairly arbitrary) "in range / out of range" cutoff for each pixel on the landscape. Dropped any pixels where the Bayesian model estimated there was a >50% chance the relative abundance was less than 0.1.  This stops the model from summing up a bunch of pixels that have (on average) a low probability of being suitable for the species but have extremely high uncertainty which can contribute substantially to the population total.

## Simulation

- Use eBird species distribution rasters as "truth" and sample from those rasters at existing PRISM survey locations.

- Analyze simulated data using Bayesian model.  

- Script "PRISM_1_Simulation.R" conducts simulations, confirms that confidence intervals are appropriate (95% credible intervals contain the true population total for 95% of species).  Bias is minimal.

## Empirical analysis

Currently fitting simple model for each species.  Data is from rapid surveys.

Model structure:
- Poisson error
- Gaussian random field for spatial autocorrelation
- Plot-level random effect (several plots were surveyed multiple times)
- Survey-level random effect (to account for overdispersion)
- Habitat covariates (proportion of landscape within 5 km of survey comprised by different land cover classes from LCC2020)

Model output:
- Maps of relative density
- Estimates of total abundance, uncorrected for detection

![SESA](https://github.com/davidiles/Arctic-PRISM-Spatial-Model/blob/main/output/empirical_SESA.png)


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