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

# Timeout INLA after 3 minutes (if it has not fit by then, it has likely stalled)
inla.setOption(inla.timeout = 60*3)
sf_use_s2(FALSE)
setwd("C:/Users/IlesD/OneDrive - EC-EC/Iles/Projects/Shorebirds/Arctic-PRISM-Spatial-Model/script")

`%!in%` <- Negate(`%in%`)

# -----------------------------------------------------
# Load data package prepared by script "PRISM_Empirical"
# -----------------------------------------------------

load("C:/Users/IlesD/OneDrive - EC-EC/Iles/Projects/Shorebirds/Arctic-PRISM-Spatial-Model/output/PRISM_wksp.RData")

# -----------------------------------------------------
# Use ebird relative abundance maps as "truth" in simulations
# - assess if they can be recovered based on PRISM sampling
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

for (i in (1:nrow(species_to_run))){
  
  set.seed(111)
  
  species <- species_to_run[i,]
  print(species$common_name)
  
  if (file.exists("../output/simulation_results_nbinomial.RDS")) results <- readRDS("../output/simulation_results_nbinomial.RDS")
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
  sdat$log_offset <- log(sdat$Proportion_Surveyed * sdat$Plot_Area_km2)
  sdat <- subset(sdat,Proportion_Surveyed > 0.5)
  sdat$Plot_Year <- paste0(sdat$PlotID,"-",sdat$Year) %>% as.factor() %>% as.numeric()
  
  sdat$count <- rpois(nrow(sdat),lambda = sdat$true * sdat$Plot_Area_km2 * sdat$Proportion_Surveyed)
  sdat$present <- as.numeric(sdat$count > 0)
  
  range(sdat$count)
  
  # If the species was detected in fewer than 50 of plots, skip it
  nplots <- sdat %>% 
    as.data.frame() %>%
    group_by(PlotID) %>%
    summarize(det = as.numeric(sum(count)>0))
  
  if (sum(nplots$det)<50) next
  
  # -----------------------------------------------------
  # Reproject ebird raster
  # -----------------------------------------------------
  
  ebirdSDM <- terra::project(ebirdSDM, arctic_proj) 
  ebirdSDM <- ebirdSDM %>% crop(vect(study_region), mask = TRUE)
  ebirdSDM_pixel_area <- res(ebirdSDM)[1] * res(ebirdSDM)[2] # km2
  
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
    ropedrag_effect(rep(1, nrow(.data.)), mean.linear = 0, prec.linear = 16) +
    intensive_effect(rep(1, nrow(.data.)), mean.linear = 0, prec.linear = 16) +
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
    family = "poisson",
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
    family = "poisson",
    data = subset(sdat, Survey_Method == "intensive"),
    
    formula = count ~ 
      intercept_rapid +
      intensive_effect +
      spde_count + 
      plot + 
      log_offset +
      PC1_beta1+
      PC2_beta1+
      PC3_beta1+
      PC1_beta2+
      PC2_beta2+
      PC3_beta2)
  
  like_ropedrag <- like(
    family = "poisson",
    data = subset(sdat, Survey_Method == "rope drag"),
    
    formula = count ~ 
      intercept_rapid +
      ropedrag_effect +
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
    like_ropedrag,
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
  var_plt <-  1/fit$summary.hyperpar$'0.5quant'[5]
  
  pred <- generate(
    fit,
    sampling_frame_sf,
    ~ intercept_rapid + 
      spde_count + 
      PC1_beta1+
      PC2_beta1+
      PC3_beta1+
      PC1_beta2+
      PC2_beta2+
      PC3_beta2,
    n.samples = 2000)
  
  pred <- exp(pred + 0.5*var_plt)
  
  pred_df <- sampling_frame_sf
  pred_df <- cbind(as.data.frame(pred_df),st_coordinates(pred_df))
  
  # Median predictions
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
  
  #map2
  
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
  
  #map3
  
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
    xlab("Population Estimate")+
    ylab("")+
    ggtitle("Population Estimate")+
    coord_cartesian(xlim=xlim)+
    scale_x_continuous(labels = comma)
  estimate_posterior
  
  # -----------------------------------------------------
  # Plot
  # -----------------------------------------------------
  
  species_maps <- ggarrange(map2,map3,estimate_posterior,nrow = 3, ncol=1,align = "hv")
  
  png(paste0("../output/simulation_",species$common_name,"_nbinomial.png"), height=9, width=4, units="in", res = 600)
  print(species_maps)
  dev.off()
  
  #maps[[species$common_name]] <- species_maps
  
  #if (file.exists("../output/simulation_maps.RDS")) maps <- readRDS("../output/simulation_maps.RDS")
  
  #saveRDS(maps,"../output/simulation_maps.RDS")
  
  # --------------------------------
  # Save results
  # --------------------------------
  if (file.exists("../output/simulation_results_nbinomial.RDS")) results <- readRDS("../output/simulation_results_nbinomial.RDS")
  
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
  
  saveRDS(results,"../output/simulation_results_nbinomial.RDS")
  
  # --------------------------------
  # Plot results
  # --------------------------------
  
  result_plot <- ggplot(data = results)+
    
    geom_errorbarh(aes(y = species_number, xmin = sum_lcl, xmax = sum_ucl), height = 0.01, col = "dodgerblue", size = 2)+
    geom_point(aes(y = species_number, x = sum_est), col = "dodgerblue", size = 5)+
    
    geom_point(aes(y = species_number, x = sum_true), fill = "white", size = 4, pch = 23)+
    geom_point(aes(y = species_number, x = sum_true), fill = "black", size = 2, pch = 23)+
    
    scale_x_continuous(name = "Population Estimate", labels = comma,trans="log10")+
    scale_y_continuous(labels = results$common_name, name = "", breaks = 1:nrow(results))+
    
    theme_bw()+
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank())+
    ggtitle("True vs estimated population sizes")
  
  print(result_plot)
  
  rm(pred)
  rm(pred_df)
  
}

png("../output/simulation_species_estimates_nbinomial.png", height=6, width=8, units="in", res = 600)
print(result_plot)
dev.off()

# Credible interval coverage
mean(results$sum_lcl < results$sum_true & results$sum_ucl > results$sum_true) # 0.75

# Bias
exp(mean(log(results$sum_est) - log(results$sum_true))) # 1.15

# Proportion of over and under-estimates
mean(results$sum_est > results$sum_true) # 0.75 of estimates are over-estimates
