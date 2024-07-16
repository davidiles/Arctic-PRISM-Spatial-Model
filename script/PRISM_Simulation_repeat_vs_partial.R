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

species_list <-  c("Baird's Sandpiper",
                   "Least Sandpiper",
                   "Ruddy Turnstone")

species_to_run <- subset(ebirdst_runs, common_name %in% species_list)
species_list[species_list %!in% species_to_run$common_name]

results <- data.frame()
maps <- list()

for (sim_run in (50:100)){
  for (sd in c(0,0.3)){
    for (i in (1:nrow(species_to_run))){
      
      start <- Sys.time()
      
      set.seed(sim_run)
      
      species <- species_to_run[i,]
      print(species$common_name)
      
      if (file.exists("../output/simulation_results_repeat_vs_partial.RDS")) results <- readRDS("../output/simulation_results_repeat_vs_partial.RDS")
      
      # -----------------------------------------------------
      # if this combination of settings has already been run, skip it
      # -----------------------------------------------------
      
      if (nrow(results)>0){
        if (nrow(subset(results, common_name == species$common_name & sd_field == sd & run == sim_run)) > 0 ) next
      }
      
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
      
      values(ebirdSDM) <- values(ebirdSDM)*10
      values(ebirdSDM)[is.na(values(ebirdSDM))] <- 0
      
      # -----------------------------------------------------
      # Simulate a change surface across the landscape
      # -----------------------------------------------------
      
      # Grid on which to simulate change (assume change will be a coarse process)
      change_grid <- st_make_grid(
        study_region_outline,
        cellsize = units::set_units(50*50,km^2),
        what = "polygons",
        square = TRUE,
        flat_topped = FALSE)%>%
        st_as_sf() %>%
        rename(geometry = x)
      
      coords <- change_grid %>% st_centroid() %>% st_coordinates()
      change_field_sim <- geoR::grf(grid = coords,
                                    cov.pars = c(0.5,1000))
      
      change_sim <- scale(change_field_sim$data)*sd + log(0.7)
      change_grid$change <- change_sim
      sampling_frame_sf$change <- c(as.data.frame(st_intersection(change_grid,sampling_frame_sf))$change)
      change_df <- cbind(as.data.frame(sampling_frame_sf),st_coordinates(sampling_frame_sf))
      
      change_rast <- rast(change_df[,c("X","Y","change")], type="xyz", crs = arctic_proj)
      
      lims <- max(abs(sampling_frame_sf$change))
      change_map <- ggplot()+
        
        geom_spatraster(data = change_rast) +
        geom_sf(data=surveys,colour="black", size = 0.1)+
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
        scale_fill_gradient2(low=muted("red"),mid="white",high=muted("blue"),midpoint=0, 
                             na.value=NA, 
                             limits = c(-lims,lims),
                             label = comma)+
        ggtitle("Simulated change surface")
      
      change_map
      
      sampling_frame_sf$lambda_1 = extract(ebirdSDM,vect(sampling_frame_sf %>% st_transform(crs(ebirdSDM))))[,2]
      sampling_frame_sf$lambda_2 = exp(log(sampling_frame_sf$lambda_1) + sampling_frame_sf$change)
      
      # -----------------------------------------------------
      # Extract expected counts and population change at survey location
      # Simulate nbinomial observations
      # -----------------------------------------------------
      
      sdat <- surveys
      sdat <- st_intersection(sdat,change_grid)
      sdat$true_1 <- extract(ebirdSDM,vect(sdat %>% st_transform(crs(ebirdSDM))))[,2]
      sdat$true_2 <- exp(log(sdat$true_1) + sdat$change)
      
      sdat$log_offset <- log(sdat$Proportion_Surveyed * sdat$Plot_Area_km2)
      sdat <- subset(sdat,Proportion_Surveyed > 0.5)
      sdat$Plot_Year <- paste0(sdat$PlotID,"-",sdat$Year) %>% as.factor() %>% as.numeric()
      
      sdat$count_1 <- rpois(nrow(sdat), lambda = sdat$true_1 * sdat$Plot_Area_km2 * sdat$Proportion_Surveyed)
      sdat$count_2 <- rpois(nrow(sdat), lambda = sdat$true_2 * sdat$Plot_Area_km2 * sdat$Proportion_Surveyed)
      
      sdat$present_1 <- as.numeric(sdat$count_1 > 0)
      sdat$present_2 <- as.numeric(sdat$count_2 > 0)
      
      # If the species was detected in fewer than 50 of plots, skip it
      nplots <- sdat %>% 
        as.data.frame() %>%
        group_by(PlotID) %>%
        summarize(det = as.numeric(sum(count_1)>0))
      
      if (sum(nplots$det)<50) next
      
      # -----------------------------------------------------
      # Create large hexagon grid across study area to represent "sampling regions"
      # for different sampling scenarios
      # -----------------------------------------------------
      
      hexgrid <- st_make_grid(study_region_outline, cellsize = 200, square=FALSE, what = "polygons") %>%
        st_as_sf() %>% mutate(hexid = 1:nrow(.)) %>% dplyr::rename(geometry = x)
      
      sdat <- sdat %>% st_intersection(hexgrid) %>% arrange(Year,Month,Day,PlotID)
      
      # -----------------------------------------------------
      # Survey scenario: half of hexagons surveyed during round 1, half during round 2
      # -----------------------------------------------------
      
      all_hexagons <- unique(sdat$hexid)
      
      # Hexagons to be surveyed during first cycle
      hexagons_cycle1 <- sample(unique(sdat$hexid),round(length(all_hexagons)*(2/3)))
      
      # Hexagons to be resurveyed during second cycle
      hexagons_resurveyed <- sample(hexagons_cycle1,round(length(hexagons_cycle1)/2))
      
      # New hexagons to be surveyed during cycle 2
      hexagons_new <- sample(all_hexagons[-which(all_hexagons %in% hexagons_cycle1)],length(hexagons_cycle1)-length(hexagons_resurveyed))
      
      hexagons_cycle2 <- c(hexagons_resurveyed,hexagons_new)
      
      # -----------------------------------------------------
      # Reproject ebird raster
      # -----------------------------------------------------
      
      ebirdSDM <- terra::project(ebirdSDM, arctic_proj) 
      ebirdSDM <- ebirdSDM %>% crop(vect(study_region), mask = TRUE)
      ebirdSDM_pixel_area <- res(ebirdSDM)[1] * res(ebirdSDM)[2] # km2
      
      # --------------------------------
      # Prepare priors + model components
      # --------------------------------
      
      matern_count <- inla.spde2.pcmatern(mesh_spatial,
                                          prior.range = c(200,0.1)  , 
                                          prior.sigma = c(0.1,0.1)  ,
                                          constr = TRUE
      )
      
      pc_prec <- list(prior = "pcprec", param = c(0.1, 0.1))
      
      sdat$plot_idx <- as.numeric(as.factor(sdat$Plot_Year))
      
      comps <- ~
        spde_count(geometry, model = matern_count) +
        spde_change(geometry, model = matern_count) +
        
        intercept_rapid(rep(1, nrow(.data.))) +
        mean_change(rep(1, nrow(.data.)))+
        plot(plot_idx, model = "iid", constr = TRUE, hyper = list(prec = pc_prec))+
        PC1_beta1(PC1, model = "linear", mean.linear = 0, prec.linear = 100)+
        PC1_beta2(PC1^2, model = "linear", mean.linear = 0, prec.linear = 100)+
        PC2_beta1(PC2, model = "linear", mean.linear = 0, prec.linear = 100)+
        PC2_beta2(PC1^2, model = "linear", mean.linear = 0, prec.linear = 100)+
        PC3_beta1(PC3, model = "linear", mean.linear = 0, prec.linear = 100)+
        PC3_beta2(PC1^2, model = "linear", mean.linear = 0, prec.linear = 100)
      
      
      # ****************************************************************************
      # ****************************************************************************
      # ANALYSIS WHEN RESURVEYING SAME LOCATIONS
      # ****************************************************************************
      # ****************************************************************************
      
      # --------------------------------
      # Model formulas
      # --------------------------------
      
      like_rapid_1 <- like(
        family = "poisson",
        data = subset(sdat, Survey_Method == "rapid" & hexid %in% hexagons_cycle1) ,
        
        formula = count_1 ~ 
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
      
      like_rapid_2 <- like(
        family = "poisson",
        data = subset(sdat, Survey_Method == "rapid" & hexid %in% hexagons_cycle1),
        
        formula = count_2 ~ 
          intercept_rapid + 
          spde_count + 
          mean_change +
          spde_change +
          plot + 
          log_offset +
          PC1_beta1+
          PC2_beta1+
          PC3_beta1+
          PC1_beta2+
          PC2_beta2+
          PC3_beta2)
      
      fit_repeat <- bru(
        comps,
        like_rapid_1,
        like_rapid_2,
        options = list(bru_verbose = 4)
      )
      
      # Estimate of population total during phase 1 of sampling
      N_1_repeat <- generate(
        fit_repeat,
        sampling_frame_sf,
        ~ exp(intercept_rapid + 
                spde_count + 
                PC1_beta1+
                PC2_beta1+
                PC3_beta1+
                PC1_beta2+
                PC2_beta2+
                PC3_beta2),
        n.samples = 1000)
      
      # Estimate of population total during phase 2 of sampling
      N_2_repeat <- generate(
        fit_repeat,
        sampling_frame_sf,
        ~ exp(intercept_rapid + 
                spde_count + 
                mean_change +
                spde_change +
                PC1_beta1+
                PC2_beta1+
                PC3_beta1+
                PC1_beta2+
                PC2_beta2+
                PC3_beta2),
        n.samples = 1000)
      
      # Estimate of population change between cycles
      change_repeat <- generate(
        fit_repeat,
        sampling_frame_sf,
        ~ log(sum(exp(intercept_rapid + 
                        spde_count + 
                        mean_change +
                        spde_change +
                        PC1_beta1+
                        PC2_beta1+
                        PC3_beta1+
                        PC1_beta2+
                        PC2_beta2+
                        PC3_beta2))/
                sum(exp(intercept_rapid + 
                          spde_count + 
                          PC1_beta1+
                          PC2_beta1+
                          PC3_beta1+
                          PC1_beta2+
                          PC2_beta2+
                          PC3_beta2)))
        ,
        n.samples = 1000)
      
      # ****************************************************************************
      # ****************************************************************************
      # ANALYSIS WHEN SURVEYING 50% NEW SITES AND RE-SURVEYED 50% OF OLD SITES
      # ****************************************************************************
      # ****************************************************************************
      
      # --------------------------------
      # Model formulas
      # --------------------------------
      
      like_rapid_1 <- like(
        family = "poisson",
        data = subset(sdat, Survey_Method == "rapid" & hexid %in% hexagons_cycle1) ,
        
        formula = count_1 ~ 
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
      
      like_rapid_2 <- like(
        family = "poisson",
        data = subset(sdat, Survey_Method == "rapid" & hexid %in% hexagons_cycle2),
        
        formula = count_2 ~ 
          intercept_rapid + 
          spde_count + 
          mean_change +
          spde_change +
          plot + 
          log_offset +
          PC1_beta1+
          PC2_beta1+
          PC3_beta1+
          PC1_beta2+
          PC2_beta2+
          PC3_beta2)
      
      fit_partial <- bru(
        comps,
        like_rapid_1,
        like_rapid_2,
        options = list(bru_verbose = 4)
      )
      
      # Estimate of population total during phase 1 of sampling
      N_1_partial <- generate(
        fit_partial,
        sampling_frame_sf,
        ~ exp(intercept_rapid + 
                spde_count + 
                PC1_beta1+
                PC2_beta1+
                PC3_beta1+
                PC1_beta2+
                PC2_beta2+
                PC3_beta2),
        n.samples = 1000)
      
      # Estimate of population total during phase 2 of sampling
      N_2_partial <- generate(
        fit_partial,
        sampling_frame_sf,
        ~ exp(intercept_rapid + 
                spde_count + 
                mean_change +
                spde_change +
                PC1_beta1+
                PC2_beta1+
                PC3_beta1+
                PC1_beta2+
                PC2_beta2+
                PC3_beta2),
        n.samples = 1000)
      
      # Estimate of population change between cycles
      change_partial <- generate(
        fit_partial,
        sampling_frame_sf,
        ~ log(sum(exp(intercept_rapid + 
                        spde_count + 
                        mean_change +
                        spde_change +
                        PC1_beta1+
                        PC2_beta1+
                        PC3_beta1+
                        PC1_beta2+
                        PC2_beta2+
                        PC3_beta2))/
                sum(exp(intercept_rapid + 
                          spde_count + 
                          PC1_beta1+
                          PC2_beta1+
                          PC3_beta1+
                          PC1_beta2+
                          PC2_beta2+
                          PC3_beta2)))
        ,
        n.samples = 1000)
      
      end <- Sys.time() 
      
      end-start # 2 min to produce predictions
      
      # --------------------------------
      # Estimates of total population size
      # --------------------------------
      
      N1_total_repeat <- apply(N_1_repeat,2,sum)
      N2_total_repeat <- apply(N_2_repeat,2,sum)
      N1_total_partial <- apply(N_1_partial,2,sum)
      N2_total_partial <- apply(N_2_partial,2,sum)
      
      # --------------------------------
      # Save results
      # --------------------------------
      
      if (file.exists("../output/simulation_results_repeat_vs_partial.RDS")) results <- readRDS("../output/simulation_results_repeat_vs_partial.RDS")
      
      results <- rbind(results, data.frame(common_name = species$common_name,
                                           species_code = species$species_code,
                                           run = sim_run,
                                           sd_field = sd,
                                           
                                           # True total abundance during sampling period 1
                                           true_sum_N1 = sum(sampling_frame_sf$lambda_1),
                                           
                                           # Estimate of N1 under full revisit scenario
                                           est_repeat_sum_N1_q05 = quantile(N1_total_repeat,0.05),
                                           est_repeat_sum_N1_q50 = quantile(N1_total_repeat,0.50),
                                           est_repeat_sum_N1_q95 = quantile(N1_total_repeat,0.95),
                                           
                                           # Estimate of N1 under partial revisit scenario
                                           est_partial_sum_N1_q05 = quantile(N1_total_partial,0.05),
                                           est_partial_sum_N1_q50 = quantile(N1_total_partial,0.50),
                                           est_partial_sum_N1_q95 = quantile(N1_total_partial,0.95),
                                           
                                           # True total abundance during sampling period 2
                                           true_sum_N2 = sum(sampling_frame_sf$lambda_2),
                                           
                                           # Estimate of N2 under full revisit scenario
                                           est_repeat_sum_N2_q05 = quantile(N2_total_repeat,0.05),
                                           est_repeat_sum_N2_q50 = quantile(N2_total_repeat,0.50),
                                           est_repeat_sum_N2_q95 = quantile(N2_total_repeat,0.95),
                                           
                                           # Estimate of N2 under partial revisit scenario
                                           est_partial_sum_N2_q05 = quantile(N2_total_partial,0.05),
                                           est_partial_sum_N2_q50 = quantile(N2_total_partial,0.50),
                                           est_partial_sum_N2_q95 = quantile(N2_total_partial,0.95),
                                           
                                           # True log change
                                           true_log_change =  log(sum(sampling_frame_sf$lambda_2)/sum(sampling_frame_sf$lambda_1)),
                                           
                                           # Estimate of N2 under full revisit scenario
                                           est_repeat_change_q05 = quantile(change_repeat,0.05),
                                           est_repeat_change_q50 = quantile(change_repeat,0.50),
                                           est_repeat_change_q95 = quantile(change_repeat,0.95),
                                           
                                           # Estimate of N2 under partial revisit scenario
                                           est_partial_change_q05 = quantile(change_partial,0.05),
                                           est_partial_change_q50 = quantile(change_partial,0.50),
                                           est_partial_change_q95 = quantile(change_partial,0.95),
                                           
                                           species_number = NA))
      
      
      results <- results %>% arrange(true_sum_N1)
      results$species_number <- as.numeric(as.factor(results$species_code))
      
      saveRDS(results,"../output/simulation_results_repeat_vs_partial.RDS")
      
      # --------------------------------
      # Plot results
      # --------------------------------
      
      species_names <- unique(results[,c("species_number","common_name")]) %>% arrange(species_number)
      result_plot_N1 <- ggplot(data = results)+
        
        geom_errorbarh(aes(y = run, xmin = est_repeat_sum_N1_q05, xmax = est_repeat_sum_N1_q95, col = "Full repeat"), height = 0, size = 2)+
        geom_point(aes(y = run, x = est_repeat_sum_N1_q50, col = "Full repeat"), size = 5)+

        geom_errorbarh(aes(y = run+0.5, xmin = est_partial_sum_N1_q05, xmax = est_partial_sum_N1_q95, col = "Partial partial"), height = 0, size = 2)+
        geom_point(aes(y = run+0.5, x = est_partial_sum_N1_q50, col = "Partial partial"), size = 5)+
        
        geom_point(aes(y = run+0.25, x = true_sum_N1), fill = "white", size = 4, pch = 23)+
        geom_point(aes(y = run+0.25, x = true_sum_N1), fill = "black", size = 2, pch = 23)+
        
        scale_x_continuous(name = "Total abundance during period 1", labels = comma)+
        scale_y_continuous(labels = 1:max(results$run), name = "", breaks = 1:max(results$run))+
        theme_bw()+
        theme(panel.grid.minor = element_blank(),
              panel.grid.major.y = element_blank())+
        ggtitle("Total abundance during period 1")+
        facet_grid(sd_field~common_name, scales = "free")
      
      print(result_plot_N1)
      
      result_plot_N2 <- ggplot(data = results)+
        
        geom_errorbarh(aes(y = run, xmin = est_repeat_sum_N2_q05, xmax = est_repeat_sum_N2_q95, col = "Full repeat"), height = 0, size = 2)+
        geom_point(aes(y = run, x = est_repeat_sum_N2_q50, col = "Full repeat"), size = 5)+

        geom_errorbarh(aes(y = run+0.5, xmin = est_partial_sum_N2_q05, xmax = est_partial_sum_N2_q95, col = "Partial repeat"), height = 0, size = 2)+
        geom_point(aes(y = run+0.5, x = est_partial_sum_N2_q50, col = "Partial repeat"), size = 5)+
        
        geom_point(aes(y = run+0.25, x = true_sum_N2), fill = "white", size = 4, pch = 23)+
        geom_point(aes(y = run+0.25, x = true_sum_N2), fill = "black", size = 2, pch = 23)+
        
        scale_x_continuous(name = "Total abundance during period 2", labels = comma)+
        scale_y_continuous(labels = 1:max(results$run), name = "", breaks = 1:max(results$run))+
        theme_bw()+
        theme(panel.grid.minor = element_blank(),
              panel.grid.major.y = element_blank())+
        ggtitle("Total abundance during period 2")+
        facet_grid(sd_field~common_name, scales = "free")
      
      print(result_plot_N2)
      
      result_plot_change <- ggplot(data = results)+
        
        geom_errorbarh(aes(y = run, xmin = est_repeat_change_q05, xmax = est_repeat_change_q95, col = "Full repeat"), height = 0, size = 2)+
        geom_point(aes(y = run, x = est_repeat_change_q50, col = "Full repeat"), size = 5)+
        
        geom_errorbarh(aes(y = run+0.5, xmin = est_partial_change_q05, xmax = est_repeat_change_q95, col = "Partial repeat"), height = 0, size = 2)+
        geom_point(aes(y = run+0.5, x = est_partial_change_q50, col = "Partial repeat"), size = 5)+
        
        geom_point(aes(y = run+0.25, x = true_log_change), fill = "white", size = 4, pch = 23)+
        geom_point(aes(y = run+0.25, x = true_log_change), fill = "black", size = 2, pch = 23)+
        
        scale_x_continuous(name = "Estimate of log-scale population change", labels = comma)+
        scale_y_continuous(labels = 1:max(results$run), name = "", breaks = 1:max(results$run))+
        theme_bw()+
        theme(panel.grid.minor = element_blank(),
              panel.grid.major.y = element_blank())+
        ggtitle("Estimate of log-scale population change")+
        facet_grid(sd_field~common_name, scales = "free")
      
      print(result_plot_change)
      
    }
  }
}

# -------------------------------------------------
# Estimates of total population size in first "cycle" of surveys
#   - which data structure leads to better predictions?
# -------------------------------------------------

if (file.exists("../output/simulation_results_repeat_vs_partial.RDS")) results <- readRDS("../output/simulation_results_repeat_vs_partial.RDS")

result_summary <- results %>%
  #subset(common_name == "Baird's Sandpiper") %>%
  group_by(common_name,sd_field) %>%
  summarize(
    
    
    # -------------------------------------------------
    # Estimates of total population size in first "cycle" of surveys
    #   - which data structure leads to better predictions?
    # -------------------------------------------------
    
    median_bias_N1_repeat = median(log(est_repeat_sum_N1_q50/true_sum_N1)),
    median_bias_N1_partial = median(log(est_partial_sum_N1_q50/true_sum_N1)),
    prob_repeat_better_N1 = mean((abs(est_repeat_sum_N1_q50) - true_sum_N1) < (abs(est_partial_sum_N1_q50) - true_sum_N1)),
    mean_coverage_N1_repeat = mean(est_repeat_sum_N1_q05 < true_sum_N1 & est_repeat_sum_N1_q95 > true_sum_N1),
    mean_coverage_N1_partial = mean(est_partial_sum_N1_q05 < true_sum_N1 & est_partial_sum_N1_q95 > true_sum_N1),
    mean_CI_N1_repeat = mean(est_repeat_sum_N1_q95 - est_repeat_sum_N1_q05),
    mean_CI_N1_partial = mean(est_partial_sum_N1_q95 - est_partial_sum_N1_q05),
    
    # -------------------------------------------------
    # Estimates of total population size in second "cycle" of surveys
    #   - which data structure leads to better predictions?
    # -------------------------------------------------
    
    median_bias_N2_repeat = median(log(est_repeat_sum_N2_q50/true_sum_N2)),
    median_bias_N2_partial = median(log(est_partial_sum_N2_q50/true_sum_N2)),
    prob_repeat_better_N2 = mean((abs(est_repeat_sum_N2_q50) - true_sum_N2) < (abs(est_partial_sum_N2_q50) - true_sum_N2)),
    mean_coverage_N2_repeat = mean(est_repeat_sum_N2_q05 < true_sum_N2 & est_repeat_sum_N2_q95 > true_sum_N2),
    mean_coverage_N2_partial = mean(est_partial_sum_N2_q05 < true_sum_N2 & est_partial_sum_N2_q95 > true_sum_N2),
    mean_CI_N2_repeat = mean(est_repeat_sum_N2_q95 - est_repeat_sum_N2_q05),
    mean_CI_N2_partial = mean(est_partial_sum_N2_q95 - est_partial_sum_N2_q05),
    
    # -------------------------------------------------
    # Estimates of change between each cycle of surveys
    #   - which data structure leads to better predictions?
    # -------------------------------------------------
    
    median_bias_change_repeat = median(est_repeat_change_q50 - true_log_change),
    median_bias_change_partial = median(est_partial_change_q50 - true_log_change),
    prob_repeat_better_change = mean((abs(est_repeat_change_q50) - true_log_change) < (abs(est_partial_change_q50) - true_log_change)),
    mean_coverage_change_repeat = mean(est_repeat_change_q05 < true_log_change & est_repeat_change_q95 > true_log_change),
    mean_coverage_change_partial = mean(est_partial_change_q05 < true_log_change & est_partial_change_q95 > true_log_change),
    mean_CI_change_repeat = mean(est_repeat_change_q95 - est_repeat_change_q05),
    mean_CI_change_partial = mean(est_partial_change_q95 - est_partial_change_q05)
    
  ) %>%
  as.data.frame()
result_summary

ggplot() +
  geom_bar(data = result_summary, aes(x = "Full Repeat", y = mean_CI_N1_repeat), stat = "identity")+
  geom_bar(data = result_summary, aes(x = "Partial Repeat", y = mean_CI_N1_partial), stat = "identity")+
  
  facet_grid(common_name~sd_field, scales = "free")+
  ylab("Width of 90% CI")+
  xlab("")+
  theme_bw()

ggplot() +
  geom_bar(data = result_summary, aes(x = "Full Repeat", y = mean_CI_N2_repeat), stat = "identity")+
  geom_bar(data = result_summary, aes(x = "Partial Repeat", y = mean_CI_N2_partial), stat = "identity")+
  
  facet_grid(common_name~sd_field, scales = "free")+
  xlab("")+
  theme_bw()

ggplot() +
  geom_bar(data = result_summary, aes(x = "Full Repeat", y = mean_CI_change_repeat), stat = "identity")+
  geom_bar(data = result_summary, aes(x = "Partial Repeat", y = mean_CI_change_partial), stat = "identity")+
  facet_grid(common_name~sd_field, scales = "free")+
  xlab("")+
  theme_bw()

# Predictions:

# Estimates of species abundance should be more precise using partial revisits
#         - this effect possibly less pronounced when population growth is variable across space?

# Estimates of population change should be more precise when using full revisits
