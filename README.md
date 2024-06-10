# Arctic PRISM Spatial Model
 Bayesian spatial analysis of shorebird abundance for Arctic PRISM.  Goal is to compare design-based analysis to model-based analysis.
 
 
![PRISM survey locations](https://github.com/davidiles/Arctic-PRISM-Spatial-Model/blob/main/output/PRISM_survey_locations.png)
 
## Model

Currently using [Zero-Adjusted Poisson models](https://inlabru-org.github.io/inlabru/articles/zip_zap_models.html) because there are large portions of the study area that are not suitable, which could be better modeled by separating the "zero" process from the non-zero process

## Simulation


## Issues

1) The model seems to overestimate abundance when there are large portions of the landscape that are unsuitable (true abundance = 0). This is probably because the model cannot estimate an abundance of exactly zero (or negative), and any uncertainty will only add positive values.

2) Covariates need to be refined (currently using land cover of canada 2020)

3) Consider implementing [barrier models] (https://haakonbakkagit.github.io/btopic128.html)

4) Need to get design-based analysis working to compare bias/precision to model-based estimates.  Even if model-based estimates aren't perfect, they may be better than design-based ones.