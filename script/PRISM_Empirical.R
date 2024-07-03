# -----------------------------------------------------
# Load / prepare study area
# -----------------------------------------------------

library(tidyverse)
library(INLA)
library(inlabru)
library(sf)
library(ebirdst)
library(tidyterra)
library(terra)
library(ggspatial)
library(exactextractr)
library(factoextra)
library(viridis)
library(stars)
library(ggpubr)
library(scales)
library(INLAspacetime)
library(patchwork)
library(lubridate)
library(ggrepel)

rm(list=ls())

# Timeout INLA after 10 minutes (if it has not fit by then, it has likely stalled)
inla.setOption(inla.timeout = 60*10)
sf_use_s2(FALSE)
setwd("C:/Users/IlesD/OneDrive - EC-EC/Iles/Projects/Shorebirds/Arctic-PRISM-Spatial-Model/script")

`%!in%` <- Negate(`%in%`)

# ****************************************************************
# ****************************************************************
# Load data
# ****************************************************************
# ****************************************************************

# -----------------------------------------------------
# Load BCR Shapefile
# -----------------------------------------------------

arctic_proj <- "+proj=lcc +lat_0=79 +lon_0=-101 +lat_1=80 +lat_2=77 +x_0=29500000 +y_0=3500000 +ellps=GRS80 +units=km +no_defs "

#load BCR boundaries
study_region <- read_sf("../data/Spatial_Covariates/BCR/BCR_Terrestrial_master.shp") %>%
  filter(COUNTRY %in% c("CANADA","USA")) %>%
  st_transform(crs = arctic_proj) %>%
  subset(WATER == 3)

study_region$BCR_PROV <- paste0(study_region$PROVINCE_S,"_",study_region$BCR) %>% as.factor()

# -----------------------------------------------------
# Load species counts at each survey location (count data)
# -----------------------------------------------------

rawdat <- survey_data <- read.csv("../data/fromChristine/PRISM files for Dave/PRISM_Contacts_and_Survey_Counts_20220110.csv")

# -----------------------------------------------------
# Load survey locations
# -----------------------------------------------------

# Location information
survey_locations <- read.csv("../data/fromChristine/PRISM files for Dave/PRISM_Plot_Coordinates_20220110.csv")

survey_coords <- survey_locations %>%
  subset(Coordinate_Type %in% c("NW Corner","NE Corner","SE Corner","SW Corner")) %>%
  group_by(PlotID) %>%
  summarize(lat = mean(LL_Coordinate1,na.rm = TRUE),
            lon = mean(LL_Coordinate2,na.rm = TRUE)) %>%
  na.omit() %>%
  st_as_sf(coords = c("lon", "lat"),crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") %>%
  st_transform(arctic_proj)

surveys <- unique(na.omit(survey_data[,c("PlotID","Plot_type","Survey_Method","Year","Month","Day","Plot_Area_km2","Proportion_Surveyed")]))
surveys <- left_join(survey_coords,surveys) %>% na.omit() # Drops surveys that have no location information

# BCR/PROV at each location
surveys <- st_intersection(surveys,study_region)

# Remove strata that have fewer than 10 surveys
n_surveys_per_stratum <- table(surveys$BCR_PROV)
strata_to_keep <- n_surveys_per_stratum[n_surveys_per_stratum>=10]

surveys <- subset(surveys, BCR_PROV %in% names(strata_to_keep))
study_region <- subset(study_region, BCR_PROV %in% names(strata_to_keep))
study_region$BCR_PROV <- factor(study_region$BCR_PROV, levels = unique(study_region$BCR_PROV))
surveys$BCR_PROV <- factor(surveys$BCR_PROV, levels = unique(study_region$BCR_PROV))
surveys$Year <- factor(surveys$Year, levels = seq(min(surveys$Year),max(surveys$Year)))

study_region_outline <- st_union(study_region)

# ****************************************************************
# ****************************************************************
# Exploratory analyses
# ****************************************************************
# ****************************************************************

# -----------------------------------------------------
# Plot survey locations
# -----------------------------------------------------

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
# 
# # -----------------------------------------------------
# # Evaluate survey intensity through time
# # -----------------------------------------------------
# 
# map2 <- ggplot()+
# 
#   geom_sf(data=study_region,colour="gray80", fill = "transparent")+
#   geom_sf(data=surveys,size = 1)+
# 
#   coord_sf(xlim = xlim, ylim = ylim, crs = arctic_proj)+
# 
#   theme_bw()+
# 
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))+
#   scale_fill_gradient(low = 'white', high = 'red', na.value=NA, name = "incl. prob.")+
#   facet_wrap(Year~., drop = FALSE)+
#   ggtitle("Survey locations")
# 
# map2
# 
# png("../output/PRISM_survey_locations_year.png", height=12, width=12, units="in", res = 600)
# print(map2)
# dev.off()
# 
# # # -----------------------------------------------------
# # # How many sites were visited more than once?
# # # -----------------------------------------------------
# #
# # nsurveys <- surveys$PlotID %>% table() %>% sort(decreasing = TRUE)
# # sum(nsurveys>1) # 239 plots surveyed more than once
# #
# # map3 <- ggplot()+
# #
# #   geom_sf(data=study_region,colour="gray80", fill = "transparent")+
# #   geom_sf(data=subset(surveys, PlotID %in% names(nsurveys[nsurveys>1])),size = 1)+
# #
# #   coord_sf(xlim = xlim, ylim = ylim, crs = arctic_proj)+
# #
# #   theme_bw()+
# #
# #   theme(panel.grid.major = element_blank(),
# #         panel.grid.minor = element_blank(),
# #         panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))+
# #   scale_fill_gradient(low = 'white', high = 'red', na.value=NA, name = "incl. prob.")+
# #   ggtitle("Plots that have been surveyed more than once")
# #
# # map3
# 
# # -----------------------------------------------------
# # Prepare matrix object containing species counts associated with each survey
# # -----------------------------------------------------
# 
# # Restrict survey_data to surveys that will be included in the analysis
# surveys$survey_ID <- paste(surveys$PlotID,surveys$Year,surveys$Month,surveys$Day,sep="-")
# survey_data$survey_ID <- paste(survey_data$PlotID,survey_data$Year,survey_data$Month,survey_data$Day,sep="-")
# survey_data <- subset(survey_data, survey_ID %in% surveys$survey_ID)
# 
# # Only include relevant species groups
# survey_data <- subset(survey_data, Species_Group %!in% c("not recorded","not applicable","Mammal"))
# 
# # Only include relevant count type
# survey_data <- subset(survey_data, Count_Type %in% c("Unknown Sex","Pair","Male","Female"))
# 
# # Only include observations collected inside plots
# survey_data <- subset(survey_data, Sighting_Type %in% c("inside plot"))
# 
# # Multiplier
# survey_data$Count_adjusted <- survey_data$Count * survey_data$Multiplier_Individuals
# 
# # Sum up counts of each species on each survey
# survey_count_df <- survey_data %>%
#   group_by(survey_ID,Species_Code) %>% summarize(Count = sum(Count_adjusted))
# 
# # Relative commonness
# relabund_raw <- survey_count_df %>% group_by(Species_Code) %>% summarize(Count_Sum = sum(Count>0)) %>% arrange(desc(Count_Sum))
# 
# # Only include species with at least 50 unique detections
# species_list <- subset(relabund_raw, Count_Sum >= 50 & Species_Code %!in% c("UNSH","UNRE"))$Species_Code
# length(species_list) # 34 species
# 
# # Create matrix of species counts associated with each survey
# surveys <- surveys %>% arrange(Year,Month,Day,PlotID)
# count_matrix <- matrix(0,nrow=nrow(surveys), ncol = length(species_list), dimnames = list(surveyID = surveys$survey_ID,species = species_list))
# 
# for (i in 1:nrow(count_matrix)){
# 
#   # Count data for this survey
#   survID <- rownames(count_matrix)[i]
#   surv <- subset(survey_count_df, survey_ID == survID & Species_Code %in% species_list)
# 
#   # Fill in cells assocaited with each species count
#   if (nrow(surv)>0){
#     for (j in 1:nrow(surv)){
#       count_matrix[i,surv$Species_Code[j]] <- count_matrix[i,surv$Species_Code[j]] + surv$Count[j]
#     }
#   }
# }
# 
# # -----------------------------------------------------
# # Some additional visualizations
# # -----------------------------------------------------
# 
# # Create hexagon grid across study area
# hexgrid <- st_make_grid(study_region_outline, cellsize = 50, square=FALSE, what = "polygons") %>%
#   st_as_sf() %>% mutate(hexid = 1:nrow(.)) %>% dplyr::rename(geometry = x)
# 
# hexcent <- st_centroid(hexgrid)
# 
# # Intersect with survey information
# surveys_hex <- surveys %>% st_intersection(hexgrid) %>% arrange(Year,Month,Day,PlotID)
# 
# # For each species in each hexagon, calculate mean count
# species_plots <- list()
# for (species in species_list){
# 
#   species_hex <- surveys_hex %>% mutate(count = count_matrix[,species])
# 
#   # Summarize total effort and mean count in each hexagon
#   species_hex <- species_hex %>%
#     as.data.frame() %>%
#     group_by(PlotID,hexid) %>%
#     summarize(count = mean(count),
#               effort = mean(Plot_Area_km2 * Proportion_Surveyed)) %>%
#     group_by(hexid) %>%
#     summarize(count = sum(count)/sum(effort),
#               effort = sum(effort))
# 
#   species_hex <- as.data.frame(species_hex) %>%
#     left_join(hexgrid,.)
# 
#   col_lim <- c(0,max(species_hex$count,na.rm = TRUE))
# 
#   species_plot <- ggplot() +
#     geom_sf(data=study_region_outline,colour="gray70", fill = "gray80")+
#     geom_sf(data=species_hex, aes(fill= count), col = "transparent")+
#     geom_sf(data=subset(species_hex,count>0), col = "black", fill = "transparent",size = 0.1)+
# 
#     annotation_scale(style = "ticks",
#                      text_face = "bold")+
# 
#     annotation_north_arrow(which_north = "true",
#                            location = "tr",
#                            pad_x = unit(0.25, "cm"), pad_y = unit(0.25, "cm"),
#                            height = unit(1, "cm"),
#                            width = unit(1, "cm"),
#                            style = north_arrow_fancy_orienteering(text_col = 'black',
#                                                                   line_col = 'gray20',
#                                                                   text_face = "bold",
#                                                                   fill = 'gray80'))+
# 
#     theme_bw()+
# 
#     theme(panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank(),
#           panel.background = element_rect(fill = 'gray80', colour = 'black'))+
#     scale_fill_gradient(low = 'white', high = 'darkred', na.value=NA,
#                         label = comma)+
#     ggtitle(paste0(species," - Observed count per km^2"))
#   species_plot
# 
#   species_plots[[species]] <- species_plot
# }
# 
# # -----------------------------------------------------
# # Prepare sampling frame
# # -----------------------------------------------------
# 
# # Create sampling frame
# sampling_frame_cell_area = 100 # km2
# sampling_frame_sf = st_make_grid(study_region, cellsize = sqrt(sampling_frame_cell_area), square = TRUE, what = "centers")
# sampling_frame_sf = sampling_frame_sf[study_region] # only select cells that intersect with study region
# sampling_frame_sf <- st_as_sf(sampling_frame_sf)
# sampling_frame_sf <- mutate(sampling_frame_sf, log_offset = log(sampling_frame_cell_area)) %>%
#   dplyr::rename(geometry = x)
# 
# # -----------------------------------------------------
# # Load spatial covariates
# # -----------------------------------------------------
# 
# # LCC2020
# lcc2020 <- rast("../data/Spatial_Covariates/LandCoverCanada2020/landcover-2020-classification.tif")
# 
# # ---------------------------------------------------
# # For each element in survey dataset, extract land cover within 5 km radius (10 km diameter)
# # ---------------------------------------------------
# 
# surveys_5km <- surveys %>% st_buffer(sqrt(sampling_frame_cell_area)/2)
# 
# prop_LCC_5km <- exact_extract(lcc2020,st_transform(surveys_5km, st_crs(lcc2020)),"frac") %>% suppressWarnings()
# names(prop_LCC_5km) <- paste0(str_replace(names(prop_LCC_5km),"frac","LCC"),"_5km")
# prop_LCC_5km[setdiff(paste0("LCC_",seq(1,19),"_5km"),names(prop_LCC_5km))] <- 0
# prop_LCC_5km <-  prop_LCC_5km %>% dplyr::select(sort(names(.)))
# 
# # Combine land cover classes
# prop_LCC_5km <- prop_LCC_5km %>%
#   mutate(
#     Needleleaf_forest_5km = LCC_1_5km + LCC_2_5km,
#     Mixed_forest_5km = LCC_5_5km + LCC_6_5km,
#     Shrub_5km = LCC_8_5km + LCC_11_5km,
#     Grass_5km = LCC_10_5km + LCC_12_5km,
#     Barren_5km = LCC_13_5km + LCC_16_5km,
#     Wetland_5km = LCC_14_5km,
#     Crop_5km = LCC_15_5km,
#     Urban_5km = LCC_17_5km,
#     Water_5km = LCC_18_5km,
#     Snow_5km = LCC_19_5km) %>%
#   dplyr::select(Needleleaf_forest_5km:Snow_5km)
# 
# surveys <- bind_cols(surveys,prop_LCC_5km)
# surveys <- na.omit(surveys)
# 
# # ---------------------------------------------------
# # For each element in sampling frame, extract land cover within 10 km
# # ---------------------------------------------------
# 
# # Proportion of each land cover class from lcc 2020
# sampling_frame_sf_5km <- sampling_frame_sf %>% st_buffer(5)
# 
# prop_LCC_5km <- exact_extract(lcc2020,st_transform(sampling_frame_sf_5km, st_crs(lcc2020)),"frac") %>% suppressWarnings()
# names(prop_LCC_5km) <- paste0(str_replace(names(prop_LCC_5km),"frac","LCC"),"_5km")
# prop_LCC_5km[setdiff(paste0("LCC_",seq(1,19),"_5km"),names(prop_LCC_5km))] <- 0
# prop_LCC_5km <-  prop_LCC_5km %>% dplyr::select(sort(names(.)))
# 
# # Combine land cover classes
# prop_LCC_5km <- prop_LCC_5km %>%
#   mutate(
#     Needleleaf_forest_5km = LCC_1_5km + LCC_2_5km,
#     Mixed_forest_5km = LCC_5_5km + LCC_6_5km,
#     Shrub_5km = LCC_8_5km + LCC_11_5km,
#     Grass_5km = LCC_10_5km + LCC_12_5km,
#     Barren_5km = LCC_13_5km + LCC_16_5km,
#     Wetland_5km = LCC_14_5km,
#     Crop_5km = LCC_15_5km,
#     Urban_5km = LCC_17_5km,
#     Water_5km = LCC_18_5km,
#     Snow_5km = LCC_19_5km) %>%
#   dplyr::select(Needleleaf_forest_5km:Snow_5km)
# 
# sampling_frame_sf <- bind_cols(sampling_frame_sf,prop_LCC_5km)
# 
# # *******************************************************************
# # *******************************************************************
# # DISCRETIZE HABITAT?
# # *******************************************************************
# # *******************************************************************
# 
# # *******************************************************************
# # *******************************************************************
# # Conduct principal components analysis
# # *******************************************************************
# # *******************************************************************
# 
# covars_for_PCA <- c("Barren_5km",
#                     "Grass_5km",
#                     "Water_5km",
#                     "Shrub_5km",
#                     "Wetland_5km")
# 
# dat_for_PCA <- surveys %>%
#   as.data.frame() %>%
#   dplyr::select(covars_for_PCA)
# 
# pca <- prcomp(dat_for_PCA, scale = TRUE)
# 
# # ------------------------------------------
# # Interpretation of specific axes (e.g., axes 1 and 2)
# # ------------------------------------------
# 
# summary(pca)   # Proportion variance explaind by axes
# fviz_eig(pca)  # Scree plot (first 5 axes explain 99.99% of variation in habitat between sites)
# pca            # Variable loadings
# 
# fviz_pca_var(pca,
#              axes = c(1,2),
#              col.var = "contrib", # Color by contributions to the PC
#              gradient.cols = viridis(10),
#              repel = TRUE     # Avoid text overlapping
# )
# 
# # ------------------------------------------
# # Predict PCA values for each survey location and standardize (mean = 0, sd = 1)
# # ------------------------------------------
# 
# surveys_PCA <- predict(pca, newdata = as.data.frame(surveys)[,names(pca$center)])
# sampling_frame_PCA <- predict(pca, newdata = as.data.frame(sampling_frame_sf)[,names(pca$center)])
# 
# for (covar in colnames(surveys_PCA)){
# 
#   covar_mean <- mean(surveys_PCA[,covar],na.rm = TRUE)
#   covar_sd <- sd(surveys_PCA[,covar],na.rm = TRUE)
# 
#   surveys_PCA[,covar] <- (as.data.frame(surveys_PCA)[,covar] - covar_mean)/covar_sd
#   sampling_frame_PCA[,covar] <- (as.data.frame(sampling_frame_PCA)[,covar] - covar_mean)/covar_sd
# 
# }
# 
# surveys <- surveys %>% bind_cols(surveys_PCA)
# sampling_frame_sf <- bind_cols(sampling_frame_sf,sampling_frame_PCA)
# 
# # ------------------------------------------
# # Plot maps of covariates
# # ------------------------------------------
# 
# VarRast <- sampling_frame_sf %>%
#   dplyr::select(covars_for_PCA,PC1:PC5) %>%
#   stars::st_rasterize()
# 
# covar_to_plot <- names(VarRast)
# 
# covar_plotlist <- list()
# 
# for (covar in covar_to_plot){
#   cplot <- ggplot() +
#     geom_stars(data = VarRast, aes(fill = !!sym(covar)))+
#     scale_fill_gradientn(colours = viridis(10), name = covar,na.value="transparent")+
#     geom_sf(data = study_region,colour="gray50",fill=NA,lwd=0.3,show.legend = F)+
#     ggtitle(covar)+
#     theme_bw()+
#     xlab("")+ylab("")
# 
#   covar_plotlist[[covar]] <- cplot
# }
# 
# covar_plots <- ggarrange(plotlist = covar_plotlist,nrow=2,ncol=length(covars_for_PCA))
# 
# png("../output/Covariate_Maps/Covariate_Maps.png", width=60, height=10, units="in", res=300, type="cairo")
# print(covar_plots)
# dev.off()
# 
# # ---------------------------
# # Create mesh and barrier feature (water)
# # Barrier example taken from: https://eliaskrainski.github.io/INLAspacetime/articles/web/barrierExample.html
# # ---------------------------
# 
# # For spatial analysis in INLA
# mesh_spatial <- fm_mesh_2d_inla(
#   boundary = st_buffer(study_region_outline,100),
#   max.edge = c(50, 200), # km inside and outside
#   cutoff = 50,
#   crs = st_crs(arctic_proj)
# )
# 
# mesh_locs <- mesh_spatial$loc[,c(1,2)] %>% as.data.frame()
# dim(mesh_locs)
# # Create barrier object for analysis
# barriers <- st_difference(bbox %>% st_buffer(10000),st_union(study_region)) %>% st_buffer(-5)
# 
# triBarrier <- unlist(fm_contains(
#   x = barriers,
#   y = mesh_spatial,
#   type = "centroid"))
# 
# triCenters.xy <- cbind(
#   mesh_spatial$loc[mesh_spatial$graph$tv[,1], 1:2] +
#     mesh_spatial$loc[mesh_spatial$graph$tv[,2], 1:2] +
#     mesh_spatial$loc[mesh_spatial$graph$tv[,3], 1:2])/3
# 
# ggplot() + theme_bw() +
#   gg(mesh_spatial) +
#   geom_point(aes(
#     x = triCenters.xy[triBarrier, 1],
#     y = triCenters.xy[triBarrier, 2]))
# 
#save.image("C:/Users/IlesD/OneDrive - EC-EC/Iles/Projects/Shorebirds/Arctic-PRISM-Spatial-Model/output/PRISM_wksp.RData")

load("C:/Users/IlesD/OneDrive - EC-EC/Iles/Projects/Shorebirds/Arctic-PRISM-Spatial-Model/output/PRISM_wksp.RData")

# -----------------------------------------------------
# Data analysis
# -----------------------------------------------------

results <- data.frame()
maps <- list()

for (species in (species_list)){
  
  #if (species %!in% species_list) next
  print(species)
  set.seed(111)
  
  if (file.exists("../output/empirical_results.RDS")) results <- readRDS("../output/empirical_results.RDS")
  if (nrow(results) > 0 & species %in% results$species) next
  sdat <- surveys %>% st_intersection(hexgrid) %>% arrange(Year,Month,Day,PlotID)
  sdat <- sdat %>% mutate(count = count_matrix[,species])
  sdat$present <- as.numeric(sdat$count > 0)
  sdat$log_offset <- log(sdat$Proportion_Surveyed * sdat$Plot_Area_km2)
  sdat <- subset(sdat,Proportion_Surveyed >0.5)
  sdat$Plot_Year <- paste0(sdat$PlotID,"-",sdat$Year) %>% as.factor() %>% as.numeric()
  
  # --------------------------------
  # How many plots were surveyed in multiple ways?
  # --------------------------------
  
  rapid_counts <- subset(sdat, Survey_Method == "rapid") %>%
    as.data.frame() %>%
    group_by(Plot_Year) %>%
    summarize(count_rapid = mean(count/exp(log_offset)))
  
  intensive_counts <- subset(sdat, Survey_Method == "intensive") %>%
    as.data.frame() %>%
    group_by(Plot_Year) %>%
    summarize(count_intensive = mean(count/exp(log_offset)))
  
  # Compare
  cmpr <- full_join(rapid_counts,intensive_counts) %>% na.omit()
  
  ggplot(data = na.omit(cmpr))+
    geom_jitter(aes(x = count_rapid, y = count_intensive))+
    theme_bw()
  
  # --------------------------------
  # Prepare priors + model components
  # --------------------------------
  
  prior_range <- c(200,0.1)         # 10% chance range is smaller than 200km
  prior_sigma <- c(0.1,0.1)         # 10% chance sd is larger than 0.1
  
  matern_count <- inla.spde2.pcmatern(mesh_spatial,
                                      prior.range = prior_range, 
                                      prior.sigma = prior_sigma,
                                      constr = TRUE
  )
  
  pc_prec <- list(prior = "pcprec", param = c(0.1, 0.01))
  
  sdat$plot_idx <- as.numeric(as.factor(sdat$Plot_Year))
  
  comps <- ~
    spde_count(geometry, model = matern_count) +
    intercept_rapid(rep(1, nrow(.data.))) +
    intensive(rep(1, nrow(.data.)), mean.linear = log(1/1.15), prec.linear = 16) +
    plot(plot_idx, model = "iid", constr = TRUE, hyper = list(prec = pc_prec))+
    PC1_beta1(PC1, model = "linear", mean.linear = 0, prec.linear = 100)+
    PC1_beta2(PC1^2, model = "linear", mean.linear = 0, prec.linear = 100)+
    PC2_beta1(PC2, model = "linear", mean.linear = 0, prec.linear = 100)+
    PC2_beta2(PC1^2, model = "linear", mean.linear = 0, prec.linear = 100)+
    PC3_beta1(PC3, model = "linear", mean.linear = 0, prec.linear = 100)+
    PC3_beta2(PC1^2, model = "linear", mean.linear = 0, prec.linear = 100)
  
  # --------------------------------
  # Model formulas
  # --------------------------------
  
  like_rapid <- like(
    family = "nbinomial",
    data = subset(sdat, Survey_Method == "rapid"),
    
    formula = count ~ 
      intercept_rapid + 
      spde_count + 
      plot + 
      log_offset +
      PC1_beta1+
      PC2_beta1+
      PC3_beta1+
      PC1_beta2+
      PC2_beta2+
      PC3_beta2)
  
  like_intensive <- like(
    family = "nbinomial",
    data = subset(sdat, Survey_Method == "intensive"),
    
    formula = count ~ 
      intercept_rapid +
      intensive +
      spde_count + 
      plot + 
      log_offset +
      PC1_beta1+
      PC2_beta1+
      PC3_beta1+
      PC1_beta2+
      PC2_beta2+
      PC3_beta2)
  
  fit <- bru(
    comps,
    like_rapid,
    like_intensive,
    options = list(bru_verbose = 4)
  )
  
  summary(fit)
  
  # ****************************************************************************
  # ****************************************************************************
  # MODEL ASSESSMENT / GOODNESS-OF-FIT
  # ****************************************************************************
  # ****************************************************************************
  
  
  # ****************************************************************************
  # ****************************************************************************
  # GENERATE PREDICTIONS FOR EVERY PIXEL ON LANDSCAPE
  # ****************************************************************************
  # ****************************************************************************
  
  # Random effect estimates
  #var_svy <- 1/fit$summary.hyperpar$'0.5quant'[3]
  var_plt <-  1/fit$summary.hyperpar$'0.5quant'[5]
  
  pred <- generate(
    fit,
    sampling_frame_sf,
    ~ intercept_rapid + 
      spde_count + 
      PC1_beta1+
      PC2_beta1+
      PC3_beta1,
    n.samples = 500)
  
  pred <- exp(pred + 0.5*var_plt)
  
  pred_df <- sampling_frame_sf
  pred_df <- cbind(as.data.frame(pred_df),st_coordinates(pred_df))
  
  # Median predictions
  pred_df$pred_med  <- apply(pred,1,function(x) median(x,na.rm = TRUE))
  
  rast_to_plot <- rast(pred_df[,c("X","Y","pred_med")], type="xyz", crs = arctic_proj)
  
  map3 <- ggplot() +
    
    geom_spatraster(data = rast_to_plot)+
    geom_sf(data=study_region_outline,colour="gray80", fill = "transparent")+
    
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
    scale_fill_gradient(low = 'white', high = 'darkred', na.value=NA, 
                        label = comma)+
    ggtitle(paste0(species," - expected count per km^2"))
  
  map3
  
  # -----------------------------------------------------
  # Estimate of sum across sampling frame
  # -----------------------------------------------------
  
  pred <- pred * sampling_frame_cell_area
  sum_est <- apply(pred,2,sum, na.rm = TRUE)
  
  xlim <- c(0,max(sum_est))
  estimate_posterior <- ggplot()+
    geom_histogram(aes(x = sum_est), fill = "dodgerblue", alpha = 0.5)+
    geom_vline(aes(xintercept = median(sum_est)), size = 2, col = "dodgerblue")+
    geom_vline(aes(xintercept = quantile(sum_est,c(0.025,0.975))), size = 1, col = "dodgerblue", linetype = 2)+
    
    theme_bw()+
    xlab("Total")+
    ylab("")+
    ggtitle("Total estimated abundance")+
    scale_x_continuous(labels = comma, limits = xlim)
  estimate_posterior
  sd(sum_est)/mean(sum_est)
  
  # -----------------------------------------------------
  # Plot
  # -----------------------------------------------------
  
  species_maps <- ggarrange(species_plots[[species]],map3,nrow = 2, ncol=1,align = "hv")
  
  png(paste0("../output/empirical_",species,"_intensive.png"), height=8, width=6, units="in", res = 600)
  print(species_maps)
  dev.off()
  
  # maps[[species]] <- species_maps
  # 
  # if (file.exists("../output/empirical_maps.RDS")) maps <- readRDS("../output/empirical_maps.RDS")
  # 
  # saveRDS(maps,"../output/empirical_maps.RDS")
  # 
  # --------------------------------
  # Save results
  # --------------------------------
  if (file.exists("../output/empirical_results.RDS")) results <- readRDS("../output/empirical_results.RDS")
  
  results <- rbind(results, data.frame(species = species,
                                       species_number = NA,
                                       n_detections = sum(sdat$count>0),
                                       sum_count = sum(sdat$count),
                                       sum_est = median(sum_est),
                                       sum_lcl = quantile(sum_est,0.025),
                                       sum_ucl = quantile(sum_est,0.975),
                                       sum_CV = sd(sum_est)/mean(sum_est),
                                       nbin_size_rapid = fit$summary.hyperpar$mean[1],
                                       nbin_size_intensive = fit$summary.hyperpar$mean[2],
                                       spde_range = fit$summary.hyperpar$mean[3],
                                       spde_sd = fit$summary.hyperpar$mean[4],
                                       plotRE_sd = sqrt(var_plt),
                                       intensive_effect = exp(fit$summary.fixed["intensive","mean"])))
  
  results <- results %>% arrange(sum_est)
  results$species_number <- 1:nrow(results)
  
  saveRDS(results,"../output/empirical_results.RDS")
  
  # --------------------------------
  # Plot results
  # --------------------------------
  
  result_plot <- ggplot(data = results)+
    
    geom_errorbarh(aes(y = species_number, xmin = sum_lcl, xmax = sum_ucl), height = 0.01, col = "dodgerblue", size = 2)+
    geom_point(aes(y = species_number, x = sum_est), col = "dodgerblue", size = 5)+
    
    scale_x_continuous(trans = "log10", name = "Estimate (sum of pixels)", labels = comma)+
    scale_y_continuous(labels = results$species, name = "", breaks = 1:nrow(results))+
    
    theme_bw()+
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank())+
    ggtitle("Estimated population sizes\n\n(uncorrected for detection)")
  
  print(result_plot)
  
  # # -----------------------------------------------------
  # # Estimate of sum within each Ecozone
  # # -----------------------------------------------------
  # 
  # # Estimates of proportion of population in each ecozone
  # stratum_proportions <- matrix(NA, nrow=length(unique(sampling_frame_sf$Ecozone)), ncol = ncol(pred))
  # 
  # for (j in 1:ncol(pred)){
  #   overall_sum = sum(pred[,j])
  #   estimates <- data.frame(Ecozone = sampling_frame_sf$Ecozone,pred = pred[,j]) %>%
  #     group_by(Ecozone) %>%
  #     summarize(prop = sum(pred)/overall_sum)
  #   stratum_proportions[,j] <- estimates$prop
  # }
  # rownames(stratum_proportions) = estimates$Ecozone
  # 
  # stratum_summary <- data.frame()
  # for (j in unique(sampling_frame_sf$Ecozone)){
  #   
  #   rows <- which(sampling_frame_sf$Ecozone == j)
  #   
  #   props <- stratum_proportions[j,]
  #   
  #   true <- sum(pred_df$true_count[rows])/sum(pred_df$true_count)
  #   mean_count <- mean(subset(sdat, Ecozone == j)$count)
  #   n_sdat <- nrow(subset(sdat, Ecozone == j))
  #   stratum_summary <- rbind(stratum_summary, data.frame(stratum = j,
  #                                                        n_sdat = n_sdat,
  #                                                        true = true,
  #                                                        est_mean = mean(props),
  #                                                        est_lcl = quantile(props,0.025),
  #                                                        est_ucl = quantile(props,0.975),
  #                                                        mean_count = mean_count))
  # }
  # 
  # ggplot()+
  #   geom_point(data = stratum_summary, aes(y = stratum, x = est_mean), col = "dodgerblue")+
  #   geom_errorbarh(data = stratum_summary, aes(y = stratum, xmin = est_lcl, xmax = est_ucl), col = "dodgerblue", height = 0)+
  #   geom_point(data = stratum_summary, aes(y = stratum, x = true))+
  #   theme_bw()
  # 
  rm(pred)
  rm(pred_df)
}

# ------------------------------------
# Comparison to existing PRISM estimates
# ------------------------------------

results <- readRDS("../output/empirical_results.RDS")

oldPrism <- read.csv("../data/fromPaul/prism_estimates_designbased.csv")
results$species[results$species %!in% oldPrism$Species_Code]

result_comparison <- full_join(oldPrism, results, by = c("Species_Code" = "species")) %>%
  subset(!is.na(sum_est)) %>%
  arrange(sum_est)

result_comparison_plot <- ggplot(data = result_comparison)+
  
  geom_errorbarh(aes(y = species_number, xmin = sum_lcl, xmax = sum_ucl), height = 0.01, col = "dodgerblue", size = 2)+
  geom_point(aes(y = species_number, x = sum_est), col = "dodgerblue", size = 5)+
  geom_point(aes(y = species_number, x = Est_Uncorrected), col = "black", size = 2)+
  
  scale_x_continuous(name = "Population Estimate", labels = comma, trans = "log10")+
  scale_y_continuous(labels = result_comparison$Species_Code, name = "", breaks = 1:nrow(result_comparison))+
  
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank())+
  ggtitle("Estimated population sizes\n\n(uncorrected for detection)")

result_comparison_plot

png("../output/species_estimates.png", height=8, width=6, units="in", res = 600)
print(result_comparison_plot)
dev.off()

lim <- range(result_comparison[,c("sum_est","Est_Uncorrected")],na.rm = TRUE)
est_comparison_plot <- ggplot(data = result_comparison)+
  
  geom_point(aes(y = sum_est, x = Est_Uncorrected))+
  geom_label_repel(aes(y = sum_est, x = Est_Uncorrected, label = Species_Code))+
  geom_abline(intercept=0,slope=1)+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank())+
  scale_x_continuous(limits=lim,trans="log10")+
  scale_y_continuous(limits=lim,trans="log10")+
  
  ylab("Bayesian")+
  xlab("Design-based")+
  ggtitle("Comparison of Population Estimate")

est_comparison_plot

lim <- range(result_comparison[,c("CV","sum_CV")],na.rm = TRUE)
CV_comparison_plot <- ggplot(data = result_comparison)+
  
  geom_point(aes(y = sum_CV, x = CV))+
  geom_label_repel(aes(y = sum_CV, x = CV, label = Species_Code))+
  geom_abline(intercept=0,slope=1)+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank())+
  coord_cartesian(xlim=lim,ylim=lim)+
  ylab("Bayesian")+
  xlab("Design-based")+
  ggtitle("Comparison of CV (lower is better)")

CV_comparison_plot

cor(result_comparison$CV,result_comparison$sum_CV,use = "complete.obs")


shorebirds <- c("SESA","REPH","WRSA","AMGP","DUNL","PESA","LESA","RNPH","BASA","SEPL","STSA","BBPL","RUTU","WISN","WHIM","HUGO")

median(1/subset(results, species %in% shorebirds)$intensive_effect)
hist(1/subset(results, species %in% shorebirds)$intensive_effect)
