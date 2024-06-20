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

survey_data <- read.csv("../data/fromChristine/PRISM files for Dave/PRISM_Contacts_and_Survey_Counts_20220110.csv")
survey_data <- survey_data %>% subset(Plot_type %in% c("Rapid"))

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

surveys <- unique(na.omit(survey_data[,c("PlotID","Year","Month","Day","Plot_Area_km2","Proportion_Surveyed")])) # n = 3106
surveys <- left_join(survey_coords,surveys) %>% na.omit() # Drops surveys that have no location information; n = 2686 remain

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

# -----------------------------------------------------
# Evaluate survey intensity through time
# -----------------------------------------------------

map2 <- ggplot()+
  
  geom_sf(data=study_region,colour="gray80", fill = "transparent")+
  geom_sf(data=surveys,size = 1)+
  
  coord_sf(xlim = xlim, ylim = ylim, crs = arctic_proj)+
  
  theme_bw()+
  
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))+
  scale_fill_gradient(low = 'white', high = 'red', na.value=NA, name = "incl. prob.")+
  facet_wrap(Year~., drop = FALSE)+
  ggtitle("Survey locations")

map2

png("../output/PRISM_survey_locations_year.png", height=12, width=12, units="in", res = 600)
print(map2)
dev.off()

# -----------------------------------------------------
# How many sites were visited more than once?
# -----------------------------------------------------

nsurveys <- surveys$PlotID %>% table() %>% sort(decreasing = TRUE)
sum(nsurveys>1) # 239 plots surveyed more than once

map3 <- ggplot()+
  
  geom_sf(data=study_region,colour="gray80", fill = "transparent")+
  geom_sf(data=subset(surveys, PlotID %in% names(nsurveys[nsurveys>1])),size = 1)+
  
  coord_sf(xlim = xlim, ylim = ylim, crs = arctic_proj)+
  
  theme_bw()+
  
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))+
  scale_fill_gradient(low = 'white', high = 'red', na.value=NA, name = "incl. prob.")+
  ggtitle("Plots that have been surveyed more than once")

map3

# -----------------------------------------------------
# Prepare matrix object containing species counts associated with each survey
# -----------------------------------------------------

# Restrict survey_data to surveys that will be included in the analysis
surveys$survey_ID <- paste(surveys$PlotID,surveys$Year,surveys$Month,surveys$Day,sep="-")
survey_data$survey_ID <- paste(survey_data$PlotID,survey_data$Year,survey_data$Month,survey_data$Day,sep="-")
survey_data <- subset(survey_data, survey_ID %in% surveys$survey_ID)

# Only include relevant species groups
survey_data <- subset(survey_data, Species_Group %!in% c("not recorded","not applicable","Mammal"))

# Only include relevant count type
survey_data <- subset(survey_data, Count_Type %in% c("Unknown Sex","Pair","Male","Female"))

# Multiplier
survey_data$Count_adjusted <- survey_data$Count * survey_data$Multiplier_Individuals

# Sum up counts of each species on each survey
survey_count_df <- survey_data %>% group_by(survey_ID,Species_Code) %>% summarize(Count = sum(Count_adjusted))

# Relative commonness
relabund_raw <- survey_count_df %>% group_by(Species_Code) %>% summarize(Count_Sum = sum(Count>0)) %>% arrange(desc(Count_Sum))

# Only include species with at least 50 unique detections
species_list <- subset(relabund_raw, Count_Sum >= 50 & Species_Code %!in% c("UNSH"))$Species_Code
length(species_list) # 46 species

# Create matrix of species counts associated with each survey
surveys <- surveys %>% arrange(Year,Month,Day,PlotID)
count_matrix <- matrix(0,nrow=nrow(surveys), ncol = length(species_list), dimnames = list(surveyID = surveys$survey_ID,species = species_list))

for (i in 1:nrow(count_matrix)){
  
  # Count data for this survey
  survID <- rownames(count_matrix)[i]
  surv <- subset(survey_count_df, survey_ID == survID & Species_Code %in% species_list)
  
  # Fill in cells assocaited with each species count
  if (nrow(surv)>0){
    for (j in 1:nrow(surv)){
      count_matrix[i,surv$Species_Code[j]] <- count_matrix[i,surv$Species_Code[j]] + surv$Count[j]
    }
  }
}


# -----------------------------------------------------
# Prepare sampling frame
# -----------------------------------------------------

# Create sampling frame
sampling_frame_cell_area = 100 # km2
sampling_frame_sf = st_make_grid(study_region, cellsize = sqrt(sampling_frame_cell_area), square = TRUE, what = "centers")
sampling_frame_sf = sampling_frame_sf[study_region] # only select cells that intersect with study region
sampling_frame_sf <- st_as_sf(sampling_frame_sf)
sampling_frame_sf <- mutate(sampling_frame_sf, log_offset = log(100)) %>%
  dplyr::rename(geometry = x)

# ---------------------------
# Create mesh and barrier feature (water)
# Barrier example taken from: https://eliaskrainski.github.io/INLAspacetime/articles/web/barrierExample.html
# ---------------------------

# For spatial analysis in INLA
mesh_spatial <- fm_mesh_2d_inla(
  boundary = st_buffer(study_region_outline,100),
  max.edge = c(50, 200), # km inside and outside
  cutoff = 50,
  crs = st_crs(arctic_proj)
)

mesh_locs <- mesh_spatial$loc[,c(1,2)] %>% as.data.frame()

# Create barrier object for analysis
barriers <- st_difference(bbox %>% st_buffer(10000),st_union(study_region)) %>% st_buffer(-5)

triBarrier <- unlist(fm_contains(
  x = barriers, 
  y = mesh_spatial, 
  type = "centroid"))

triCenters.xy <- cbind(
  mesh_spatial$loc[mesh_spatial$graph$tv[,1], 1:2] +
    mesh_spatial$loc[mesh_spatial$graph$tv[,2], 1:2] +
    mesh_spatial$loc[mesh_spatial$graph$tv[,3], 1:2])/3

ggplot() + theme_bw() +
  gg(mesh_spatial) +
  geom_point(aes(
    x = triCenters.xy[triBarrier, 1],
    y = triCenters.xy[triBarrier, 2])) 

# -----------------------------------------------------
# Data analysis
# -----------------------------------------------------

species_list <-  c("Black-bellied Plover",
                   "American Golden Plover",
                   "Pacific Golden Plover",
                   "Semipalmated Plover",
                   "Whimbrel",
                   "Hudsonian Godwit",
                   "Ruddy Turnstone",
                   "Black Turnstone",
                   "Rock Sandpiper",
                   "Purple Sandpiper",
                   "Spotted Sandpiper",
                   "Lesser Yellowlegs",
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
                   "Red Phalarope",
                   
                   "Canada Goose",
                   "Snow Goose",
                   "Cackling Goose",
                   "Tundra Swan",
                   "Northern Shoveler",
                   "Northern Pintail",
                   "Green-winged Teal",
                   "Greater Scaup",
                   "King Eider",
                   "Long-tailed Duck",
                   "Rock Ptarmigan",
                   "Sandhill Crane",
                   
                   "Lapland Longspur",
                   "Willow Ptarmigan",
                   "Long-tailed Jaeger",
                   "Pomarine Jaeger",
                   "Sabine's Gull",
                   "Herring Gull",
                   "Glaucous Gull",
                   "Arctic Tern",
                   "Red-throated Loon",
                   "Pacific Loon",
                   
                   "Northern Harrier",
                   "Rough-legged Hawk",
                   "Short-eared Owl",
                   "Peregrine Falcon",
                   "Common Raven",
                   
                   "Horned Lark",
                   "American Pipit",
                   "Redpoll",
                   "Lapland Longspur",
                   "Smith's Longspur",
                   "Snow Bunting",
                   "White-crowned Sparrow",
                   "Harris's Sparrow",
                   "Yellow Warbler",
                   "Savannah Sparrow",
                   "American Tree Sparrow")

species_to_run <- subset(ebirdst_runs, common_name %in% species_list)
species_list[species_list %!in% species_to_run$common_name]

results <- data.frame()
maps <- list()

for (i in 1:nrow(species_to_run)){
  
  set.seed(111)
  
  species <- species_to_run[i,]
  
  if (file.exists("../output/results_simulation.RDS")) results <- readRDS("../output/results_simulation.RDS")
  if (nrow(results) > 0 & species$common_name %in% results$common_name) next
  
  # -----------------------------------------------------
  # Load and prepare ebird raster
  # -----------------------------------------------------
  
  ebirdst_download_status(species$species_code,pattern = "_mean_9km_")
  
  ebirdSDM <- load_raster(species$species_code, product = "abundance", 
                          period = "seasonal", metric = "mean", 
                          resolution = "9km")
  
  if ("resident" %in% names(ebirdSDM)){
    ebirdSDM <- ebirdSDM$resident
  } else{
    ebirdSDM <- ebirdSDM$breeding
  }
  values(ebirdSDM)[is.na(values(ebirdSDM))] <- 0
  values(ebirdSDM) <- values(ebirdSDM)*10
  
  # -----------------------------------------------------
  # Extract expected counts at survey location
  # Simulate poisson observations
  # -----------------------------------------------------
  
  sdat <- surveys
  sdat$true <- extract(ebirdSDM,vect(sdat %>% st_transform(crs(ebirdSDM))))[,2]
  sdat$count <- rpois(nrow(sdat),sdat$true * sdat$Plot_Area_km2 * sdat$Proportion_Surveyed)
  sdat$present <- as.numeric(sdat$count > 0)
  sdat$log_offset <- log(sdat$Proportion_Surveyed * sdat$Plot_Area_km2)
  
  # If the species was detected in fewer than 10% of plots, skip it
  nplots <- sdat %>% 
    as.data.frame() %>%
    group_by(PlotID) %>%
    summarize(det = as.numeric(sum(count)>0))
  
  if (mean(nplots$det)<0.1) next
  
  # -----------------------------------------------------
  # Reproject ebird raster
  # -----------------------------------------------------
  
  ebirdSDM <- terra::project(ebirdSDM, arctic_proj) 
  ebirdSDM <- ebirdSDM %>% crop(vect(study_region), mask = TRUE)
  ebirdSDM_pixel_area <- res(ebirdSDM)[1] * res(ebirdSDM)[2] # km2
  
  # --------------------------------
  # Prepare priors + model components
  # --------------------------------
  
  # Controls the 'residual spatial field'.  This can be adjusted to create smoother surfaces.
  prior_range <- c(200,0.1)         # 50% chance range is smaller than 100km
  prior_sigma <- c(0.2,0.1)         # 1% chance sd is larger than 1
  matern_count <- barrierModel.define(mesh_spatial,
                                      barrier.triangles = triBarrier,
                                      prior.range = prior_range, 
                                      prior.sigma = prior_sigma,
                                      constr = TRUE,
                                      range.fraction = 0.05
  )
  
  matern_count <- inla.spde2.pcmatern(mesh_spatial,
                                       prior.range = prior_range, 
                                       prior.sigma = prior_sigma,
                                       constr = TRUE
  )

  
  pc_prec <- list(prior = "pcprec", param = c(1, 0.1))
  
  sdat$surv_idx <- as.numeric(as.factor(sdat$survey_ID))
  sdat$plot_idx <- as.numeric(as.factor(sdat$PlotID))
  
  comps <- ~
    spde_count(geometry, model = matern_count) +
    intercept_count(rep(1, nrow(.data.))) +
    surv_effect(surv_idx, model = "iid", constr = TRUE, hyper = list(prec = pc_prec)) +
    plot_effect(plot_idx, model = "iid", constr = TRUE, hyper = list(prec = pc_prec))
  
  # --------------------------------
  # Model formulas
  # --------------------------------
  
  model_formula_count = count ~ intercept_count + spde_count + log_offset + surv_effect + plot_effect
  
  count_like <- like(
    family = "poisson",
    data = sdat,
    formula = model_formula_count)
  
  fit <- bru(
    comps,
    count_like,
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
  
  pred <- generate(
      fit,
      sampling_frame_sf,
      ~ exp(intercept_count + spde_count),
      n.samples = 2500)
   
  pred_df <- sampling_frame_sf
  pred_df <- cbind(as.data.frame(pred_df),st_coordinates(pred_df))
  
  pred_df$pred_med  <- apply(pred,1,function(x) median(x,na.rm = TRUE))
  
  # -----------------------------------------------------
  # Plot true species distribution
  # -----------------------------------------------------
  
  vals <- na.omit(values(ebirdSDM)) %>% as.numeric()
  max <- max(c(vals,pred_df$pred_med))
  lim <- c(0,max) 
  
  map2 <- ggplot()+
    
    geom_spatraster(data = ebirdSDM) +
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
                        limits = lim, 
                        label = comma)+
    ggtitle("Simulated density (count / km^2)")
  
  map2
  
  # ------------------------------
  # Median predictions
  # ------------------------------
  
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
    scale_fill_gradient(low = 'white', 
                        high = 'darkred', 
                        na.value=NA, 
                        label = comma, 
                        limits = lim)+
    ggtitle(paste0(species$common_name," - expected count per km^2"))
  
  map3
  
  # -----------------------------------------------------
  # Estimate of sum across sampling frame
  # -----------------------------------------------------
  
  pred <- pred * sampling_frame_cell_area
  sum_est <- apply(pred,2,sum, na.rm = TRUE)
  sum_true <- sum(values(ebirdSDM)*ebirdSDM_pixel_area,na.rm = TRUE)
  
  xlim <- c(0,max(c(sum_true,quantile(sum_est,0.975))))
  estimate_posterior <- ggplot()+
    geom_histogram(aes(x = sum_est), fill = "dodgerblue", alpha = 0.5)+
    geom_vline(aes(xintercept = sum_true), size = 2)+
    geom_vline(aes(xintercept = median(sum_est)), size = 2, col = "dodgerblue")+
    geom_vline(aes(xintercept = quantile(sum_est,c(0.025,0.975))), size = 1, col = "dodgerblue", linetype = 2)+
    
    theme_bw()+
    xlab("Total")+
    ylab("")+
    ggtitle("Sum of all pixels")+
    coord_cartesian(xlim=xlim)+
    scale_x_continuous(labels = comma)
  estimate_posterior
  
  # -----------------------------------------------------
  # Plot
  # -----------------------------------------------------
  
  species_maps <- ggarrange(map2,map3,estimate_posterior,nrow = 3, ncol=1,align = "hv")
  
  png(paste0("../output/simulation_",species$common_name,".png"), height=9, width=4, units="in", res = 600)
  print(species_maps)
  dev.off()
  
  maps[[species$common_name]] <- species_maps
  
  if (file.exists("../output/maps_simulation.RDS")) maps <- readRDS("../output/maps_simulation.RDS")
  
  saveRDS(maps,"../output/maps_simulation.RDS")
  
  # --------------------------------
  # Save results
  # --------------------------------
  if (file.exists("../output/results_simulation.RDS")) results <- readRDS("../output/results_simulation.RDS")
  
  results <- rbind(results, data.frame(common_name = species$common_name,
                                       species_code = species$species_code,
                                       n_detections = sum(sdat$count>0),
                                       sum_count = sum(sdat$count),
                                       sum_true = sum_true,
                                       sum_est = median(sum_est),
                                       sum_lcl = quantile(sum_est,0.025),
                                       sum_ucl = quantile(sum_est,0.975),
                                       sum_CV = sd(sum_est)/mean(sum_est),
                                       species_number = NA))
  
  results <- results %>% arrange(sum_true)
  results$species_number <- 1:nrow(results)
  
  saveRDS(results,"../output/results_simulation.RDS")
  
  # --------------------------------
  # Plot results
  # --------------------------------
  
  result_plot <- ggplot(data = results)+
    
    geom_errorbarh(aes(y = species_number, xmin = sum_lcl, xmax = sum_ucl), height = 0.01, col = "dodgerblue", size = 2)+
    geom_point(aes(y = species_number, x = sum_est), col = "dodgerblue", size = 5)+
    
    geom_point(aes(y = species_number, x = sum_true), fill = "white", size = 4, pch = 23)+
    geom_point(aes(y = species_number, x = sum_true), fill = "black", size = 2, pch = 23)+
    
    scale_x_continuous(name = "Estimate (sum of pixels)", labels = comma,trans="log10")+
    scale_y_continuous(labels = results$common_name, name = "", breaks = 1:nrow(results))+
    
    theme_bw()+
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank())+
    ggtitle("True vs estimated population sizes")
  
  print(result_plot)
  
}

png("../output/simulation_species_estimates.png", height=6, width=8, units="in", res = 600)
print(result_plot)
dev.off()


# Credible interval coverage
mean(results$sum_lcl < results$sum_true & results$sum_ucl > results$sum_true) # 59%

# Bias
exp(mean(log(results$sum_est) - log(results$sum_true))) # about 3%

# Proportion of over and under-estimates
mean(results$sum_est > results$sum_true) # 61% of estimates are over-estimates
