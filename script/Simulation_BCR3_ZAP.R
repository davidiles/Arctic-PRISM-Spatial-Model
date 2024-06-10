# -----------------------------------------------------
# Load / prepare study area
# -----------------------------------------------------

library(tidyverse)
library(dplyr)
library(sf)
library(ggspatial)
library(ebirdst)
library(terra)
library(tidyterra)
library(cowplot)
library(INLA)
library(inlabru)
library(spsurvey)
library(exactextractr)
library(factoextra)
library(viridis)
library(stars)
library(ggpubr)

rm(list=ls())

# Timeout INLA after 10 minutes (if it has not fit by then, it has likely stalled)
inla.setOption(inla.timeout = 60*10)

`%!in%` <- Negate(`%in%`)

sf_use_s2(FALSE)

setwd("C:/Users/IlesD/OneDrive - EC-EC/Iles/Projects/Shorebirds/Arctic_PRISM/script")

arctic_proj <- "+proj=lcc +lat_0=79 +lon_0=-101 +lat_1=80 +lat_2=77 +x_0=29500000 +y_0=3500000 +ellps=GRS80 +units=km +no_defs "

#load BCR boundaries
allBCR <- read_sf("../data/Spatial_Covariates/BCR/BCR_Terrestrial_master.shp") %>%
  filter(COUNTRY %in% c("CANADA","USA")) %>%
  st_transform(crs = arctic_proj) %>%
  subset(WATER == 3)

study_region <- subset(allBCR, BCR == 3)
study_region$BCR_PROV <- paste0(study_region$PROVINCE_S,"_",study_region$BCR) %>% as.factor()

# -----------------------------------------------------
# Arctic PRISM survey locations
# -----------------------------------------------------

surveys <- readRDS("../data/fromChristine/GIS_shapefiles.RDS") %>%
  st_transform(st_crs(arctic_proj)) %>%
  st_centroid()

# Strata at each location
surveys <- st_intersection(surveys,study_region)

# Remove strata that have fewer than 50 surveys
n_surveys_per_stratum <- table(surveys$BCR_PROV)
strata_to_keep <- n_surveys_per_stratum[n_surveys_per_stratum>=50]

surveys <- subset(surveys, BCR_PROV %in% names(strata_to_keep))
study_region <- subset(study_region, BCR_PROV %in% names(strata_to_keep))
study_region$BCR_PROV <- factor(study_region$BCR_PROV, levels = unique(study_region$BCR_PROV))
surveys$BCR_PROV <- factor(surveys$BCR_PROV, levels = unique(study_region$BCR_PROV))

dummy_vars <- model.matrix(~BCR_PROV, data = surveys)[,-1] %>% as.data.frame()
colnames(dummy_vars) <- paste0("stratum_",2:(ncol(dummy_vars)+1))
surveys <- cbind(surveys,dummy_vars)

# -----------------------------------------------------
# Prepare sampling frame
# -----------------------------------------------------

name <- "Semipalmated Sandpiper"
code <- subset(ebirdst_runs,common_name == name)$species_code
ebirdst_download_status(code,pattern = "_mean_9km_")

# Raster used as template
ebirdSDM <- load_raster(code, product = "abundance", 
                        period = "seasonal", metric = "mean", 
                        resolution = "9km")

if ("resident" %in% names(ebirdSDM)){
  ebirdSDM <- ebirdSDM$resident
} else{
  ebirdSDM <- ebirdSDM$breeding
}

ebirdSDM <- terra::project(ebirdSDM, arctic_proj) 
sampling_frame <- ebirdSDM
values(sampling_frame) <- 1

ebirdSDM <- ebirdSDM %>% crop(vect(study_region), mask = TRUE)
sampling_frame <- sampling_frame %>% crop(vect(study_region), mask = TRUE)
values(ebirdSDM)[is.na(values(ebirdSDM)) & !is.na(values(sampling_frame))] <- 0

# -----------------------------------------------------
# Extract strata at each location
# -----------------------------------------------------

# Dummy variables for strata
sampling_frame_sf <- as.points(sampling_frame) %>% st_as_sf() %>% st_intersection(study_region)
sampling_frame_sf <- sampling_frame_sf[,-1]

dummy_vars <- model.matrix(~BCR_PROV, data = sampling_frame_sf)[,-1] %>% as.data.frame()
colnames(dummy_vars) <- paste0("stratum_",2:(ncol(dummy_vars)+1))
sampling_frame_sf <- cbind(sampling_frame_sf,dummy_vars)

# Remove strata that have fewer than 50 surveys
sampling_frame_sf <- subset(sampling_frame_sf, BCR_PROV %in% names(strata_to_keep))

# # -----------------------------------------------------
# # Use GRTS to select spatially balanced random sample of locations
# # -----------------------------------------------------
# 
# surveys <- grts(sampling_frame_sf,n_base = nrow(surveys))$sites_base

# # Plot surveys
bbox <- st_bbox(study_region) %>% st_as_sfc()
xlim <- range(as.data.frame(st_coordinates(bbox))$X)
ylim <- range(as.data.frame(st_coordinates(bbox))$Y)

map1 <- ggplot()+
  
  geom_sf(data=study_region,colour="gray80", fill = "transparent")+
  geom_sf(data=surveys,size = 1)+
  
  coord_sf(xlim = xlim, ylim = ylim, crs = arctic_proj)+
  
  annotation_scale(style = "ticks",
                   text_face = "bold")+
  
  annotation_north_arrow(which_north = "true",
                         location = "tr",
                         pad_x = unit(0.25, "cm"), pad_y = unit(0.25, "cm"),
                         height = unit(1, "cm"),
                         width = unit(1, "cm"),
                         style = north_arrow_fancy_orienteering(text_col = 'black',
                                                                line_col = 'gray20',
                                                                text_face = "bold",
                                                                fill = 'gray80'))+
  
  theme_bw()+
  
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))+
  scale_fill_gradient(low = 'white', high = 'red', na.value=NA, name = "incl. prob.")+
  ggtitle("Survey locations")

map1

png("../output/PRISM_survey_locations.png", height=4, width=6, units="in", res = 600)
print(map1)
dev.off()

# -----------------------------------------------------
# Load spatial covariates
# -----------------------------------------------------

study_region_buffer <- st_union(study_region) %>% st_buffer(10)

# Mean Annual Temperature
MAT <- rast("../data/Spatial_Covariates/Bioclimate/Normal_1991_2020_MAT.tif") %>% 
  crop(st_transform(study_region_buffer,crs(.))) %>%
  project(ebirdSDM, align = TRUE, method = "bilinear")%>%
  resample(y = ebirdSDM,method="bilinear")

# Mean Annual Precipitation
MAP <- rast("../data/Spatial_Covariates/Bioclimate/Normal_1991_2020_MAP.tif") %>% 
  crop(st_transform(study_region_buffer,crs(.))) %>%
  project(ebirdSDM, align = TRUE, method = "bilinear")%>%
  resample(y = ebirdSDM,method="bilinear")

raster_list <- list(MAT,MAP)
raster_stack <- rast(raster_list)
names(raster_stack) = c("MAT","MAP")

# LCC2020
lcc2020 <- rast("../data/Spatial_Covariates/LandCoverCanada2020/landcover-2020-classification.tif") %>% 
  crop(st_transform(study_region_buffer,crs(.)))

# ---------------------------------------------------
# For each element in survey dataset, extract land cover within 10 km
# ---------------------------------------------------

surveys = terra::extract(raster_stack,vect(surveys) , bind = TRUE) %>% st_as_sf()

# Proportion of each land cover class from lcc 2020
surveys_10km <- surveys %>% st_buffer(10)

prop_LCC_10km <- exact_extract(lcc2020,st_transform(surveys_10km, st_crs(lcc2020)),"frac") %>% suppressWarnings()
names(prop_LCC_10km) <- paste0(str_replace(names(prop_LCC_10km),"frac","LCC"),"_10km")
prop_LCC_10km[setdiff(paste0("LCC_",seq(1,19),"_10km"),names(prop_LCC_10km))] <- 0
prop_LCC_10km <-  prop_LCC_10km %>% dplyr::select(sort(names(.)))

# Combine land cover classes
prop_LCC_10km <- prop_LCC_10km %>%
  mutate(
    Needleleaf_forest_10km = LCC_1_10km + LCC_2_10km,
    Mixed_forest_10km = LCC_5_10km + LCC_6_10km,
    Shrub_10km = LCC_8_10km + LCC_11_10km,
    Grass_10km = LCC_10_10km + LCC_12_10km,
    Barren_10km = LCC_13_10km + LCC_16_10km,
    Wetland_10km = LCC_14_10km,
    Crop_10km = LCC_15_10km,
    Urban_10km = LCC_17_10km,
    Water_10km = LCC_18_10km,
    Snow_10km = LCC_19_10km) %>%
  dplyr::select(Needleleaf_forest_10km:Snow_10km)

surveys <- bind_cols(surveys,prop_LCC_10km)
surveys <- na.omit(surveys)

# ---------------------------------------------------
# For each element in sampling frame, extract land cover within 10 km
# ---------------------------------------------------

sampling_frame_sf = terra::extract(raster_stack,vect(sampling_frame_sf) , bind = TRUE) %>% st_as_sf()

# Proportion of each land cover class from lcc 2020
sampling_frame_sf_10km <- sampling_frame_sf %>% st_buffer(10)

prop_LCC_10km <- exact_extract(lcc2020,st_transform(sampling_frame_sf_10km, st_crs(lcc2020)),"frac") %>% suppressWarnings()
names(prop_LCC_10km) <- paste0(str_replace(names(prop_LCC_10km),"frac","LCC"),"_10km")
prop_LCC_10km[setdiff(paste0("LCC_",seq(1,19),"_10km"),names(prop_LCC_10km))] <- 0
prop_LCC_10km <-  prop_LCC_10km %>% dplyr::select(sort(names(.)))

# Combine land cover classes
prop_LCC_10km <- prop_LCC_10km %>%
  mutate(
    Needleleaf_forest_10km = LCC_1_10km + LCC_2_10km,
    Mixed_forest_10km = LCC_5_10km + LCC_6_10km,
    Shrub_10km = LCC_8_10km + LCC_11_10km,
    Grass_10km = LCC_10_10km + LCC_12_10km,
    Barren_10km = LCC_13_10km + LCC_16_10km,
    Wetland_10km = LCC_14_10km,
    Crop_10km = LCC_15_10km,
    Urban_10km = LCC_17_10km,
    Water_10km = LCC_18_10km,
    Snow_10km = LCC_19_10km) %>%
  dplyr::select(Needleleaf_forest_10km:Snow_10km)

sampling_frame_sf <- bind_cols(sampling_frame_sf,prop_LCC_10km)

# *******************************************************************
# *******************************************************************
# Conduct principal components analysis
# *******************************************************************
# *******************************************************************

covars_for_PCA <- c("Barren_10km",
                    "Grass_10km",
                    "Water_10km",
                    "Shrub_10km",
                    "Snow_10km",
                    "Needleleaf_forest_10km",
                    "Wetland_10km")


dat_for_PCA <- surveys %>%
  as.data.frame() %>%
  dplyr::select(covars_for_PCA)

pca <- prcomp(dat_for_PCA, scale = TRUE)

# ------------------------------------------
# Interpretation of specific axes (e.g., axes 1 and 2)
# ------------------------------------------

summary(pca)   # Proportion variance explaind by axes
fviz_eig(pca)  # Scree plot (first 5 axes explain 85% of variation in habitat between sites)
pca            # Variable loadings

fviz_pca_var(pca,
             axes = c(1,2),
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = viridis(10),
             repel = TRUE     # Avoid text overlapping
)

# ------------------------------------------
# Predict PCA values for each survey location and standardize (mean = 0, sd = 1)
# ------------------------------------------

surveys_PCA <- predict(pca, newdata = as.data.frame(surveys)[,names(pca$center)])
sampling_frame_PCA <- predict(pca, newdata = as.data.frame(sampling_frame_sf)[,names(pca$center)])

for (covar in colnames(surveys_PCA)){
  
  covar_mean <- mean(surveys_PCA[,covar],na.rm = TRUE)
  covar_sd <- sd(surveys_PCA[,covar],na.rm = TRUE)
  
  surveys_PCA[,covar] <- (as.data.frame(surveys_PCA)[,covar] - covar_mean)/covar_sd
  sampling_frame_PCA[,covar] <- (as.data.frame(sampling_frame_PCA)[,covar] - covar_mean)/covar_sd
  
}

surveys <- surveys %>% bind_cols(surveys_PCA)
sampling_frame_sf <- bind_cols(sampling_frame_sf,sampling_frame_PCA)

# ------------------------------------------
# Plot maps of covariates
# ------------------------------------------

VarRast <- sampling_frame_sf %>% 
  dplyr::select(covars_for_PCA,PC1:PC7) %>%
  stars::st_rasterize()

covar_to_plot <- names(VarRast)

covar_plotlist <- list()

for (covar in covar_to_plot){
  cplot <- ggplot() + 
    geom_stars(data = VarRast, aes(fill = !!sym(covar)))+
    scale_fill_gradientn(colours = viridis(10), name = covar,na.value="transparent")+
    geom_sf(data = study_region,colour="gray50",fill=NA,lwd=0.3,show.legend = F)+
    ggtitle(covar)+
    theme_bw()+
    xlab("")+ylab("")
  
  covar_plotlist[[covar]] <- cplot
}

covar_plots <- ggarrange(plotlist = covar_plotlist,nrow=2,ncol=length(covars_for_PCA))

png("../output/Covariate_Maps/Covariate_Maps.png", width=60, height=10, units="in", res=300, type="cairo")
print(covar_plots)
dev.off()

# -----------------------------------------------------
# Conduct simulations
# -----------------------------------------------------

sr_for_mesh <- st_union(study_region) %>% st_buffer(100)

# For spatial analysis in INLA
mesh_spatial <- fm_mesh_2d_inla(
  boundary = sr_for_mesh, 
  max.edge = c(100, 500), # km inside and outside
  cutoff = 100, 
  crs = st_crs(arctic_proj)
)

mesh_locs <- mesh_spatial$loc[,c(1,2)] %>% as.data.frame()
dim(mesh_locs)
plot(mesh_spatial)

# Need to set access key (only once) before downloading ranges
#ebirdst::set_ebirdst_access_key()
#usethis::edit_r_environ()

set.seed(111)

species_to_run <- subset(ebirdst_runs, common_name %in% c("Black-bellied Plover",
                                                          "American Golden Plover",
                                                          "Pacific Golden Plover",
                                                          "Semipalmated Plover",
                                                          "Whimbrel",
                                                          "Hudsonian Godwit",
                                                          "Bar-tailed Godwit",
                                                          "Ruddy Turnstone",
                                                          "Black Turnstone",
                                                          "Rock Sandpiper",
                                                          "Purple Sandpiper",
                                                          "Red Knot",
                                                          "Sanderling",
                                                          "Dunlin",
                                                          "Semipalmated Sandpiper",
                                                          "Western Sandpiper",
                                                          "Least Sandpiper",
                                                          "White-rumped Sandpiper",
                                                          "Baird's Sandpiper",
                                                          "Pectoral Sandpiper",
                                                          "Buff-breasted Sandpiper",
                                                          "Long-billed Dowitcher",
                                                          "Stilt Sandpiper",
                                                          "Wilson's Snipe",
                                                          "Red-necked Phalarope",
                                                          "Red Phalarope"))

results <- data.frame()
maps <- list()

for (i in 1:nrow(species_to_run)){
  
  species <- species_to_run[i,]
  
  if (species$common_name %in% results$common_name){
    species_to_run <- subset(species_to_run, !(common_name %in% results$common_name))
    next
  }
  
  ebirdst_download_status(species$species_code,pattern = "_mean_9km_")
  
  # Raster used as template
  ebirdSDM <- load_raster(species$species_code, product = "abundance", 
                          period = "seasonal", metric = "mean", 
                          resolution = "9km")
  
  if ("resident" %in% names(ebirdSDM)){
    ebirdSDM <- ebirdSDM$resident
  } else{
    ebirdSDM <- ebirdSDM$breeding
  }
  
  ebirdSDM <- terra::project(ebirdSDM, arctic_proj) 
  ebirdSDM <- ebirdSDM %>% crop(vect(study_region), mask = TRUE)
  
  # Assume NA values have zero abundance
  values(ebirdSDM)[is.na(values(ebirdSDM)) & !is.na(values(sampling_frame))] <- 0
  values(ebirdSDM) <- values(ebirdSDM) * 5
  #values(ebirdSDM) <- values(ebirdSDM)*0 + 1
  ebirdSDM_sf <- as.points(ebirdSDM) %>% st_as_sf() 
  
  # -----------------------------------------------------
  # Plot species distribution, overlaid with true species distribution
  # -----------------------------------------------------
  
  map2 <- ggplot()+
    
    geom_spatraster(data = ebirdSDM) +
    geom_sf(data=study_region,colour="gray80", fill = "transparent")+
    #geom_sf(data=surveys,size = 0.1)+
    
    #coord_sf(xlim = xlim, ylim = ylim, crs = arctic_proj)+
    
    annotation_scale(style = "ticks",
                     text_face = "bold")+
    
    annotation_north_arrow(which_north = "true",
                           location = "tr",
                           pad_x = unit(0.25, "cm"), pad_y = unit(0.25, "cm"),
                           height = unit(1, "cm"),
                           width = unit(1, "cm"),
                           style = north_arrow_fancy_orienteering(text_col = 'black',
                                                                  line_col = 'gray20',
                                                                  text_face = "bold",
                                                                  fill = 'gray80'))+
    
    theme_bw()+
    
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))+
    scale_fill_gradient(low = 'white', high = 'blue', na.value=NA)+
    ggtitle(paste0(species$common_name," - simulated density surface"))
  
  map2
  
  # -----------------------------------------------------
  # Extract abundance at survey locations
  # -----------------------------------------------------
  
  surveys$count <- extract(ebirdSDM,vect(surveys))[,2]
  surveys$count <- rpois(length(surveys$count),surveys$count)
  
  # Ensure the species is abundant enough to fit model
  if (sum(surveys$count>0, na.rm = TRUE) < 100){
    next
  }
  
  simcount_map <- ggplot()+
    
    
    geom_sf(data=study_region,colour="gray80", fill = "transparent")+
    
    geom_sf(data=surveys,col = "gray")+
    geom_sf(data=subset(surveys, count > 0),aes(size = count))+
    
    #coord_sf(xlim = xlim, ylim = ylim, crs = arctic_proj)+
    
    annotation_scale(style = "ticks",
                     text_face = "bold")+
    
    annotation_north_arrow(which_north = "true",
                           location = "tr",
                           pad_x = unit(0.25, "cm"), pad_y = unit(0.25, "cm"),
                           height = unit(1, "cm"),
                           width = unit(1, "cm"),
                           style = north_arrow_fancy_orienteering(text_col = 'black',
                                                                  line_col = 'gray20',
                                                                  text_face = "bold",
                                                                  fill = 'gray80'))+
    
    theme_bw()+
    
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))+
    scale_fill_gradient(low = 'white', high = 'blue', na.value=NA)+
    ggtitle(paste0(species$common_name," - simulated density surface"))
  
  simcount_map
  
  # --------------------------------
  # ZAP
  # --------------------------------
  
  # Controls the 'residual spatial field'.  This can be adjusted to create smoother surfaces.
  prior_range <- c(250, 0.99)       # 1% chance range is smaller than 1000km
  prior_sigma <- c(0.1,0.01)         # 1% chance sd is larger than 0.1
  matern_count <- inla.spde2.pcmatern(mesh_spatial,
                                      prior.range = prior_range, 
                                      prior.sigma = prior_sigma,
                                      constr = TRUE
  )
  
  # Controls the 'residual spatial field'.  This can be adjusted to create smoother surfaces.
  prior_range <- c(1000, 0.01)        # 1% chance range is smaller than 1000km
  prior_sigma <- c(1,0.5)          # 50% chance sd is larger than 1
  matern_present <- inla.spde2.pcmatern(mesh_spatial,
                                        prior.range = prior_range, 
                                        prior.sigma = prior_sigma,
                                        constr = TRUE
  )
  comps <- ~
    Intercept_count(1) +
    spde_count(geometry, model = matern_count) +
    
    Intercept_present(1) +
    spde_present(geometry, model = matern_present)# +
    #count_PC2(PC2, model = "linear")+
    #present_PC2(PC2, model = "linear")
    
  
  surveys$present <- as.numeric(surveys$count>0)
  
  truncated_poisson_like <-
    if (package_version(getNamespaceVersion("INLA")) < "23.10.19-1") {
      like(
        family = "zeroinflatedpoisson0",
        data = subset(surveys,present==1),
        formula = count ~ Intercept_count + spde_count ,
        control.family = list(hyper = list(theta = list(
          initial = -20, fixed = TRUE
        )))
      )
    } else {
      like(
        family = "nzpoisson",
        data = subset(surveys,present==1),
        formula = count ~ Intercept_count + spde_count  
      )
    }
  
  present_like <- like(
    family = "binomial",
    data = surveys,
    formula = present ~ Intercept_present + spde_present
  )
  
  fit_zap <- bru(
    comps,
    present_like,
    truncated_poisson_like,
    options = list(bru_verbose = 4)
  )
  
  # ****************************************************************************
  # ****************************************************************************
  # GENERATE PREDICTIONS FOR EVERY PIXEL ON LANDSCAPE
  # ****************************************************************************
  # ****************************************************************************
  
  pred_df <- sampling_frame_sf
  pred_df$true_count <- extract(ebirdSDM,vect(pred_df))[,2]
  
  pred <- generate(
    fit_zap,
    pred_df,
    ~ {
      presence_prob <- plogis(Intercept_present + spde_present)
      lambda <- exp(Intercept_count + spde_count)
      expect_param <- presence_prob * lambda
      expect <- expect_param / (1 - exp(-lambda))
      
      #variance <- expect * (1 - exp(-lambda) * expect)
      # list(
      #   presence = presence_prob,
      #   lambda = lambda,
      #   expect = expect,
      #   variance = variance#,
      #   #obs_prob = (1 - presence_prob) * (count == 0) + (count > 0) * presence_prob * dpois(count, expect_param) / (1 - dpois(0, expect_param))
      # )
      
    },
    n.samples = 2500
  )
  
  pred_df <- cbind(as.data.frame(pred_df),st_coordinates(pred_df))
  
  # Median predictions
  pred_df$pred_mean  <- apply(pred,1,function(x) mean(x,na.rm = TRUE))
  pred_df$pred_q50  <- apply(pred,1,function(x) quantile(x,0.5,na.rm = TRUE))
  pred_df$flag <- apply(pred,1,function(x) sum(x>100)>0)
  
  rast_to_plot <- rast(pred_df[,c("X","Y","pred_mean")], type="xyz", crs = arctic_proj)
  map3 <- ggplot() +
    
    geom_spatraster(data = rast_to_plot)+
    geom_sf(data=study_region,colour="gray80", fill = "transparent")+
    coord_sf(xlim = xlim, ylim = ylim, crs = arctic_proj)+
    
    annotation_scale(style = "ticks",
                     text_face = "bold")+
    
    annotation_north_arrow(which_north = "true",
                           location = "tr",
                           pad_x = unit(0.25, "cm"), pad_y = unit(0.25, "cm"),
                           height = unit(1, "cm"),
                           width = unit(1, "cm"),
                           style = north_arrow_fancy_orienteering(text_col = 'black',
                                                                  line_col = 'gray20',
                                                                  text_face = "bold",
                                                                  fill = 'gray80'))+
    
    theme_bw()+
    
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))+
    scale_fill_gradient(low = 'white', high = 'blue', na.value=NA, limits = c(0,max(values(rast_to_plot))))+
    ggtitle("Estimated density")
  
  map3
  
  # 
  # # Which cells have unrealistic predictions?
  # rast_to_plot <- rast(pred_df[,c("X","Y","flag")], type="xyz", crs = arctic_proj)
  # 
  # map3 <- ggplot() +
  #   
  #   geom_spatraster(data = rast_to_plot)+
  #   geom_sf(data=study_region,colour="gray80", fill = "transparent")+
  #   coord_sf(xlim = xlim, ylim = ylim, crs = arctic_proj)+
  #   
  #   annotation_scale(style = "ticks",
  #                    text_face = "bold")+
  #   
  #   annotation_north_arrow(which_north = "true",
  #                          location = "tr",
  #                          pad_x = unit(0.25, "cm"), pad_y = unit(0.25, "cm"),
  #                          height = unit(1, "cm"),
  #                          width = unit(1, "cm"),
  #                          style = north_arrow_fancy_orienteering(text_col = 'black',
  #                                                                 line_col = 'gray20',
  #                                                                 text_face = "bold",
  #                                                                 fill = 'gray80'))+
  #   
  #   theme_bw()+
  #   
  #   theme(panel.grid.major = element_blank(),
  #         panel.grid.minor = element_blank(),
  #         panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))+
  #   scale_fill_gradientn(colours = c("red","white","blue"), na.value=NA)+
  #   ggtitle("Estimated density")
  # 
  # map3
  # 
  # Proportion of estimated total population in pixels where species is truly absent
  # j = which(pred_df$true_count > 0)
  # prop <- apply(pred,2,function(x) sum(x[j])/sum(x))
  # sum_est_present <- apply(pred[j,],2,sum, na.rm = TRUE)
  # sum_true_present <- sum(pred_df[j,"true_count"])
  # hist(sum_est_present)
  # abline(v = sum_true_present, lwd = 2)
  # 
  # j = which(pred_df$true_count == 0)
  # sum_est_absent <- apply(pred[j,],2,sum, na.rm = TRUE)
  # sum_true_absent <- sum(pred_df[j,"true_count"])
  # hist(sum_est_absent)
  # abline(v = sum_true_absent, lwd = 2)

  # -----------------------------------------------------
  # Estimate of sum
  # -----------------------------------------------------
  
  sum_est <- apply(pred,2,sum, na.rm = TRUE)
  
  # True sum
  sum_true <- sum(pred_df$true_count,na.rm = TRUE)
  
  ggplot()+
    geom_histogram(aes(x = sum_est), fill = "dodgerblue")+
    geom_vline(aes(xintercept = sum_true), size = 2)+
    theme_bw()+
    xlab("Total")+
    ylab("")+
    ggtitle("Sum of all pixels")
  
  # -----------------------------------------------------
  # Plot
  # -----------------------------------------------------
  
  species_maps <- plot_grid(map2, map3, nrow = 1, align = "hv")
  maps[[species$common_name]] <- species_maps
  
  # --------------------------------
  # Design-based analysis
  # --------------------------------
  
  # N <- nrow(sampling_frame_sf)
  # n <- nrow(surveys)
  # surveys$wgt <- N/n
  # 
  # design_analysis <- spsurvey::cont_analysis(surveys, 
  #                                           vars = "count", 
  #                                           weight = "wgt")
  # 
  
  # --------------------------------
  # Save results
  # --------------------------------
  
  results <- rbind(results, data.frame(common_name = species$common_name,
                                       species_code = species$species_code,
                                       species_number = NA,
                                       sum_true = sum_true,
                                       sum_est = median(sum_est),
                                       sum_lcl = quantile(sum_est,0.025),
                                       sum_ucl = quantile(sum_est,0.975)))
  
  results <- results %>% arrange(sum_true)
  results$species_number <- 1:nrow(results)
  
  # --------------------------------
  # Plot results
  # --------------------------------
  # 
  # lim <- range(c(results$sum_est,results$sum_true))
  # lim[1] <- lim[1] / 1.2
  # lim[2] <- lim[2] *1.2
  # 
  # result_plot <- ggplot(data = results)+
  #   geom_errorbar(aes(x = sum_true, ymin = sum_lcl, ymax = sum_ucl), width = 0, col = "dodgerblue")+
  #   geom_point(aes(x = sum_true, y = sum_est), col = "dodgerblue")+
  #   geom_text(aes(x = sum_true, y = sum_est, label = common_name), col = "dodgerblue", hjust = -0.1)+
  #   geom_abline(slope=1,intercept=0,linetype=2)+
  #   coord_cartesian(xlim = lim, ylim = lim)+
  #   scale_y_continuous(trans = "log10")+
  #   scale_x_continuous(trans = "log10")+
  #   
  #   theme_bw()
  # 
  # print(result_plot)
  # 
  # 
  # --------------------------------
  # Plot results
  # --------------------------------
  
  result_plot <- ggplot(data = results)+
    
    geom_errorbarh(aes(y = species_number, xmin = sum_lcl, xmax = sum_ucl), height = 0.01, col = "dodgerblue", size = 2)+
    geom_point(aes(y = species_number, x = sum_est), col = "dodgerblue", size = 5)+
    
    geom_point(aes(y = species_number, x = sum_true), fill = "white", size = 4, pch = 23)+
    geom_point(aes(y = species_number, x = sum_true), fill = "black", size = 2, pch = 23)+
    
    scale_x_continuous(trans = "log10")+
    scale_y_continuous(labels = results$common_name, name = "", breaks = 1:nrow(results))+
    
    theme_bw()+
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank())
  
  print(result_plot)
  
}
