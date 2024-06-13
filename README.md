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

#### Example results below

Simulated surface, sample data, and fitted density surface for SESA:

![SESA simulation](https://github.com/davidiles/Arctic-PRISM-Spatial-Model/blob/main/output/Semipalmated%20Sandpiper.png)
 
The resulting species distribution raster looks pretty good. Other species are shown in the [output folder](https://github.com/davidiles/Arctic-PRISM-Spatial-Model/tree/main/output) .

To evaluate how well the model performs for estimating abundance over the entire species arctic range, sum up the estimates and compare to the "true" sum.  Results are shown for 18 species below.

![Population estimates](https://github.com/davidiles/Arctic-PRISM-Spatial-Model/blob/main/output/species_estimates.png)
 
*Note: Also eventually analyze the same data using design-based analysis* 

## Issues

1) Covariates need to be refined (currently only using ecozones from EBAR shapefile).  Additional covariates would be helpful.

2) Precision in simulations is likely overestimated.  Real data will be noisier and will require additional variance components.

3) Need to get design-based analysis working to compare bias/precision to model-based estimates.  Even if model-based estimates aren't perfect, they may be better than design-based ones.  
   - Adapt Brandon Edwards' code.
   - Contact Jon Bart.