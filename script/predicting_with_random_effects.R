library(tidyverse)
library(INLA)
library(inlabru)

rm(list=ls())

results <- data.frame()
for (reps in 1:1000){
  
  dat <- expand.grid(site = 1:10000,
                     survey = 1:100)
  
  site_effect <- rnorm(max(dat$site),0,1)
  
  dat$site_effect <- site_effect[dat$site]
  dat$lambda <- exp(log(10) + dat$site_effect)
  
  dat$count <- rnbinom(nrow(dat),mu = dat$lambda, size = 0.8)
  range(dat$count)
  
  sdat <- sample_n(dat,1000)
  sdat$site <- as.numeric(as.factor(sdat$site))
  
  # ------------------------------------------------------
  # Inlabru
  # ------------------------------------------------------
  
  pc_prec <- list(prior = "pcprec", param = c(2, 0.1))
  
  comps <- ~
    intercept_count(rep(1, nrow(.data.))) +
    st(survey, constr = TRUE, hyper = list(prec = pc_prec),model = "iid", mapper = bru_mapper_index(max(sdat$site)))
  
  model_formula_count = count ~ intercept_count + st
  
  count_like <- like(
    family = "nbinomial",
    data = sdat,
    formula = model_formula_count)
  
  fit <- bru(
    comps,
    count_like,
    options = list(bru_verbose = 4)
  )
  
  summary(fit)
  
  # ****************************************************************************
  # GENERATE PREDICTIONS FOR EVERY PIXEL ON LANDSCAPE
  # ****************************************************************************
  
  sampling_frame <- data.frame(site = 1:max(dat$site))
  
  # True sum
  true_sum <- dat %>%
    group_by(site) %>%
    summarize(mean = mean(count)) %>%
    ungroup() %>%
    summarize(sum = sum(mean))
  
  # Predictions
  # accounting for random effects
  pred <- generate(fit,
                   data = sampling_frame,
                   formula = ~ intercept_count,
                   n.samples = 1000)
  
  var_st <-  1/fit$summary.hyperpar$'0.5quant'[2]
  
  pred <- exp(pred + 0.5*var_st)
  
  # Estimate of sum
  sum <- apply(pred,2,sum)
  
  # Predictions_naive
  pred_naive <- generate(fit,
                         data = sampling_frame,
                         formula = ~ exp(intercept_count),
                         n.samples = 1000)
  
  # Estimate of sum
  sum_naive <- apply(pred_naive,2,sum)
  
  results <- rbind(results, data.frame(rep = reps,
                                       true_sum = true_sum$sum,
                                       
                                       lcl = quantile(sum,0.025),
                                       med = median(sum),
                                       ucl = quantile(sum,0.975),
                                       
                                       lcl_naive  = quantile(sum_naive ,0.025),
                                       med_naive  = median(sum_naive ),
                                       ucl_naive  = quantile(sum_naive ,0.975)
  ))
  
  result_plot <- ggplot(data = results) +
    geom_errorbar(aes(x = rep, ymin = lcl, ymax = ucl), width = 0, size = 2, col = "dodgerblue")+
    geom_point(aes(x = rep, y = med),size = 4, col = "dodgerblue")+
    
    geom_errorbar(aes(x = rep, ymin = lcl_naive, ymax = ucl_naive), width = 0, size = 2, col = "orangered")+
    geom_point(aes(x = rep, y = med_naive),size = 4, col = "orangered")+
    
    geom_point(aes(x = rep, y = true_sum),size = 2, col = "black")+
    theme_bw()
  
  print(result_plot)
}

# Credible interval coverage when including random effect variance in summation
mean(results$lcl<results$true_sum & results$ucl > results$true_sum) # 87%

# Credible interval coverage when omitting random effect variance in summation
mean(results$lcl_naive<results$true_sum & results$ucl_naive > results$true_sum) # 64%
