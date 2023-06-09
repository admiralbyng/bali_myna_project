########################################################################
# PROJECT: Bali myna population viability analysis
# CONTENTS: 
#  - run matrix population model in file 'bali-myna_MPM_penalised_supplementation.R' and generate outputs used to populate Results section
#  DEPENDENCIES:
#  - model code file (bali-myna_MPM_penalised_supplementation.R)
# AUTHOR: tom squires
# LAST UPDATED: 27/01/2023
########################################################################

# set seed to give same results every time script is run
set.seed(101)

####---------------------------####
# create folder to hold results and figures
####---------------------------####

# check working directory is correct - if not, change
getwd()
# setwd(...)

dir()

# create a container for the results in working directory
res_fold <- file.path(paste0("./bm_pva_results_",Sys.Date()))
if(!dir.exists(res_fold)) dir.create(res_fold, recursive = TRUE) 

# now use "res_fold" to save all outputs into for each run


####---------------------------####
# read in the baseline MPM function (needs to in your working directory)
####---------------------------####

source("./bali-myna_MPM_penalised_supplementation.R")


####------------------------------------------------####
# read in packages required for script. the package will then be installed if it is not available
####------------------------------------------------####

install.load.package <- function(x){
  if(!require(x, character.only = TRUE)) 
    install.packages(x, repos = 'http://cran.us.r-project.org')
  require(x, character.only = TRUE)
}

package_vector <- c("ggplot2","egg","popbio","mvtnorm","msm") 

sapply(package_vector, install.load.package)

####-------------------------------------------------####
# function for basic plot of the entire simulation, with all simulation runs shown in grey and a mean trend line in black
####-------------------------------------------------####
plot_MPM <- function(NTot = model_name, 
                     title_name = "", 
                     y_axis = "t", 
                     list_element = 1, # this is needed because function output is a list - the first element of the basic function is the NTot matrix with the model data
                     plot_future_only = FALSE, # if true you get rid of first 10 years
                     confidence_intervals = TRUE,
                     y_upper_limit = 4000,
                     y_axis_label = "Total population"
                     ) {
    
  NTot <- NTot[[list_element]] # first element of results list is the NTot matrix
  
  if(plot_future_only == TRUE){
    NTot <- NTot[ ,11:length(NTot[1,])] # this is to filter out the historical years
  }
  
  plot(NTot[1,], 
         xaxt = "n",
         xaxs = "r",
         yaxs = "i",
         yaxt = y_axis,
         type = "l", 
         lwd = 1, 
         lty = 1, 
         xlab = expression(bold("")),
         ylab = y_axis_label, 
         ylim = c(0, y_upper_limit), 
         col="grey72")
    
    for (i in 2:nrow(NTot)){
      lines(NTot[i,], type = "l", lwd = 1, lty = 1,col="grey72")
    }
    # plots the results of each of the simulations through a loop. The x axis is purposefully blanked to fill with the appropriate years instead of the model years
    
    axis(1, at=1:length(NTot[1,]), 
         labels = c(as.numeric(names(NTot[1,])[1]):as.numeric(names(NTot[1,])[length(NTot[1,])])))
    
    trend <- colMeans(NTot)
    # calculates the mean of totals from each simulation for each given year i.e. mean trend of all population trajectories
    
    lines(trend, type = "l", lwd = 2, lty = 1)
    # plots the mean trend, the 0.5 is to get it out of being a female-only total
    
    if(confidence_intervals == TRUE){
      q5 <- apply(NTot, MARGIN = 2, function(x) quantile(x, probs = 0.025))
      q95 <- apply(NTot, MARGIN = 2, function(x) quantile(x, probs = 0.975))    
      lines(q5, type = "l", lwd = 1, lty = 2)
      lines(q95, type = "l", lwd = 1, lty = 2) 
      }
    
    title(title_name)
  }
  
####-------------------------------------------------####
# function for plotting the historic and model data (MODEL VALIDATION)
####-------------------------------------------------####

plot_model_validation <- function(NTot = model_name) {
    
  NTot <- NTot[[1]] # the first element of the results list is the NTot matrix
  validation <- NTot[,1:10]
    
    years <- names(validation[1,])
    nyears <- length(validation[1,])
    
    censusData <- matrix(c(1:10, 15,32,48,57,81,109,184,256,341,420),ncol = 2)
  
    plot(validation[1,],
         xaxt = "n",
         xaxs = "r",
         yaxs = "i",
         type = "l",
         lwd = 1,
         lty = 1,
         xlab = expression(bold("")),ylab = "Total population",
         ylim = c(0,max(c(max(validation), max(censusData)))*1.05),col="grey72")
    
    for(i in 2:length(validation[,1])){
      lines(validation[i,], type = "l", lwd = 1, lty = 1,col="grey72")
         }
    
    axis(1, at = 1:nyears, labels=c(years[1]:years[nyears]))
    
    trend <- colMeans(validation)
    
    q5 <- apply(validation, MARGIN = 2, function(x) quantile(x, probs = 0.025))
    
    q95 <- apply(validation, MARGIN = 2, function(x) quantile(x, probs = 0.975))
    
    lines(trend, type = "l", lwd = 2, lty = 1)

    lines(q5, type = "l", lwd = 1, lty = 2)
    
    lines(q95, type = "l", lwd = 1, lty = 2)
    
    points(x = censusData[,1], y = censusData[,2], type = "p", col = "red", pch = 16)
   
    legend(x = "topleft",
           inset = 0.02,
           legend = "Population census data",  # Legend texts
           pch = 16,           
           col = "red",
           bty = "n",
           lty = 0)     
    
    # title("Model validation")
    
    }
  
  # show the options available in the function
  args(balimynaMPM)
  
####-------------------------------------------####
#### Model 1 - baseline model using the established default settings 
####-------------------------------------------####

mod1_baseline <- balimynaMPM() 
  
  plot_MPM(mod1_baseline, plot_future_only = TRUE)
  plot_model_validation(NTot = mod1_baseline)
  
  # what is the final mean and median population size? 
  mean(mod1_baseline[[1]][,"2032"])
  median(mod1_baseline[[1]][,21])
  
  ## compare simulated vs. validated population sizes for validation period 
  
  simulation_est <- round(as.numeric(colMeans(mod1_baseline[[1]])),digits = 0)
  simulation_est <- simulation_est[1:10]
  census <- c(15,32,48,57,81,109,184,256,341,420)
  
 proportional_difference <- (1-simulation_est/census) * 100
  proportional_difference <- round(abs(proportional_difference),digits = 1)
  proportional_difference
  
  mean(proportional_difference)
  
  num_difference <- abs(simulation_est - census)
  num_difference
  # use table() to check how many population simulations reached certain status by x date
  table(mod1_baseline[[1]][,"2032"] < 1000)
  
  ## how many birds are juvenile? 
  mod1_meanPropJuvenile <- apply(mod1_baseline[[2]][,11:20], 
                                 MARGIN = 2, 
                                 function(x) mean(x))
  mod1_meanPropJuvenile
  mean(mod1_meanPropJuvenile)
  # result shows % of 1st years in the future years modelled

####-------------------------------------------####
#### FIGURE 1. Model validation
####-------------------------------------------####
  
  png(filename = paste0(file.path(res_fold, 
                        "fig1_model_validation.png")),
      width = 700,
      height = 500,
      pointsize = 14)

  par(mfrow = c(1,1),
      mar = c(3,4,3,3))

  plot_model_validation(NTot = mod1_baseline)

  dev.off()

####-------------------------------------------####
## model 2a > SUPPLEMENTATION CEASES 
####-------------------------------------------####

  mod2a_supplementation_stops <- balimynaMPM(future.supplementation = 0)
  
  ## how many birds are juvenile? 
  mod2a_meanPropJuvenile <- apply(mod2a_supplementation_stops[[2]][,11:21], 
                                 MARGIN = 1, 
                                 function(x) mean(x))
  mean(mod2a_meanPropJuvenile)
  # result shows % of birds juvenile in the future years modelled
  
  table(mod2a_supplementation_stops[[1]][,"2032"] < 1000)[2]/10# % chance of being below 1000 individuals
  
  mean(mod2a_supplementation_stops[[1]][,"2032"])
  
####-------------------------------------------####
## model 2b > SUPPLEMENTATION CEASES AFTER 5 YEARS
####-------------------------------------------####
  
  mod2b_supplementation_stops_delayed <- balimynaMPM(delayedSupplementation = TRUE)
  
  table(mod2b_supplementation_stops_delayed[[1]][,"2032"] < 1000)[2]/10 # % chance of population being below 1000 at the end of the modelling period
  
  mean(mod2b_supplementation_stops_delayed[[1]][,"2032"])
  

####-------------------------------------------####
#### FIGURE 2. Baseline model with supplementation vs. without supplementation vs. delayed end to supplementation 
####-------------------------------------------####

  png(filename = paste0(file.path(res_fold, "/fig2_baseline_model_95CI_v2.png")),
      width = 750,
      height = 300,
      pointsize = 12)
  par(oma=c(0,3,0,0),mar=c(3,2,3,1.5),mfrow=c(1,3))
  plot_MPM(NTot = mod1_baseline,
           title = "Continue supplementation",
           plot_future_only = TRUE,
           confidence_intervals = TRUE,
           y_upper_limit = 5000,
           y_axis_label = ""
           )
  # abline(h = mean(mod1_baseline[[1]][,21]), col = "red", lty = "dashed")

   plot_MPM(NTot = mod2b_supplementation_stops_delayed,
           title = "Cease supplementation after 2026",
           plot_future_only = TRUE,
           confidence_intervals = TRUE,
           y_upper_limit = 5000,
           y_axis_label = ""
  )

   plot_MPM(NTot = mod2a_supplementation_stops,
            title = "Cease supplementation after 2021 \n (final baseline year)",
            plot_future_only = TRUE,
            confidence_intervals = TRUE,
            y_upper_limit = 5000,
            y_axis_label = ""
   )
  #
  # abline(h = mean(mod2_supplementation_stops[[1]][,21]), col = "red", lty = "dashed")

  mtext(text="Total population",side=2,line=1,outer=TRUE, cex = 1)

  dev.off()
  
####-------------------------------------------####
#### SCENARIOS: 
####  - trapping birds (above background rate in baseline model)
####  - supplementing the population
####-------------------------------------------####
  
  # two types of trapping: 
  # 1. nest trapping - only first year birds removed from pop
  # 2. mist-netting - adults and first years taken in random proportions
  
  # three trapping volumes
  # low (40 birds/year), med (80 birds/year), high (120 birds/year)
  
  # supplementation
  # continues
  # stops
  # ceases after 5 years
 
  # run models using all the different scenario settings
  
  # Scenario 1: Nest trapping + supplementation continues
  # Scenario 2: Nest trapping + supplementation stops immediately
  # Scenario 3: Nest trapping + supplementation stops after 5 years
  # Scenario 4: Mist netting + supplementation continues
  # Scenario 5: Mist netting + supplementation stops
  # Scenario 6: Mist netting + supplementation stops after 5 years
  
 #### Scenario 1. nest trapping + supplementation continues ####
  nest_trapping_sim <- c(40, 80, 120)
  nest_trapping_list <- list()
  nest_trapping_out <- matrix(0, 
                              nrow = length(nest_trapping_sim),        
                              ncol = length(mod1_baseline[[1]][1,]))
  
  for(scenario in 1:length(nest_trapping_sim)) {
    
    res <- balimynaMPM(trapping = TRUE,
                       trappingType = "nest",
                       trappedN = nest_trapping_sim[scenario])
    
    nest_trapping_list[[scenario]] <- res[[1]]
    nest_trapping_out[scenario,] <- colMeans(nest_trapping_list[[scenario]])
  }
  nest_trapping_fig <- nest_trapping_out[,11:21]
  
  png(filename = paste0(file.path(res_fold,
                        "nest_poaching_supp_continues.png")))
  plot(nest_trapping_fig[1,], 
       xaxt = "n", 
       xaxs = "r",
       yaxs = "i",
       type = "l",lwd = 1,
       lty = 1,
       xlab = expression(bold("")),
       ylab = "Total population",
       ylim = c(0,5000),
       main = "nest poaching, immature birds taken,\n supplementation continues")
  axis(1, 
       at = c(1,6,11),
       labels = c("2022","2027","2032"))
  for(i in 2:length(nest_trapping_fig[,1])){
    lines(nest_trapping_fig[i,], type = "l", lwd = 1, lty = 1)
  }
  dev.off()
  
  #### Scenario 2. Nest trapping + supplementation stops ####
  nest_trapping_noSupp_list <- list()
  nest_trapping_noSupp_out <- matrix(0, 
                                     nrow = length(nest_trapping_sim),        
                                     ncol = length(mod1_baseline[[1]][1,]))
  
  for(scenario in 1:length(nest_trapping_sim)) {
    
    res <- balimynaMPM(trapping = TRUE,
                       trappingType = "nest",
                       trappedN = nest_trapping_sim[scenario],
                       future.supplementation = 0)
    
    nest_trapping_noSupp_list[[scenario]] <- res[[1]]
    nest_trapping_noSupp_out[scenario,] <- colMeans(nest_trapping_noSupp_list[[scenario]])
    
  }
  
  nest_trapping_noSupp_fig <- nest_trapping_noSupp_out[,11:21]
  
  png(filename = paste0(file.path(res_fold,
                        "nest_poaching_supplementation_stops.png")))
  
  plot(nest_trapping_noSupp_fig[1,], 
       xaxt = "n", 
       xaxs = "r",
       yaxs = "i",
       type = "l",lwd = 1,
       lty = 1,
       xlab = expression(bold("")),
       ylab = "Total population",
       ylim = c(0,4000),
       main = "nest poaching, immature birds taken,\nsupplementation stops")
  axis(1, 
       at = c(1,6,11),
       labels = c("2022","2027","2032"))
  for(i in 2:length(nest_trapping_noSupp_fig[,1])){
    lines(nest_trapping_noSupp_fig[i,], type = "l", lwd = 1, lty = 1)
  }
    dev.off()
  
  #### Scenario 3. nest trapping + supplementation stops after 5 years ####
  nest_trapping_supp5_list <- list()
  nest_trapping_supp5_out <- matrix(0, 
                                     nrow = length(nest_trapping_sim),        
                                     ncol = length(mod1_baseline[[1]][1,]))
  
  for(scenario in 1:length(nest_trapping_sim)) {
    
    res <- balimynaMPM(trapping = TRUE,
                       trappingType = "nest",
                       trappedN = nest_trapping_sim[scenario],
                       delayedSupplementation = TRUE)
    
    nest_trapping_supp5_list[[scenario]] <- res[[1]]
    nest_trapping_supp5_out[scenario,] <- colMeans(nest_trapping_supp5_list[[scenario]])
    
  }
  
  nest_trapping_supp5_fig <- nest_trapping_supp5_out[,11:21]
  
  png(filename = paste0(file.path(res_fold, "nest_poaching_supp_stops_5y.png")))
  
  plot(nest_trapping_supp5_fig[1,], 
       xaxt = "n", 
       xaxs = "r",
       yaxs = "i",
       type = "l",lwd = 1,
       lty = 1,
       xlab = expression(bold("")),
       ylab = "Total population",
       ylim = c(0,4000),
       main = "nest poaching, immature birds taken,\nsupplementation ceases after 5 years")
  axis(1, 
       at = c(1,6,11),
       labels = c("2022","2027","2032"))
    for(i in 2:length(nest_trapping_supp5_fig[,1])){
    lines(nest_trapping_supp5_fig[i,], type = "l", lwd = 1, lty = 1)
  }
    dev.off()
  
  #### Scenario 4. Mist netting + supplementation continues ####
    
  mistNet_trapping_sim <- c(40,80,120)
  mistNet_trapping_list <- list()
  mistNet_trapping_out <- matrix(0, 
                                 nrow = length(mistNet_trapping_sim),        
                                 ncol = length(mod1_baseline[[1]][1,]))
  
  for(scenario in 1:length(mistNet_trapping_sim)) {
    
    res <- balimynaMPM(trapping = TRUE,
                       trappingType = "mist-net",
                       trappedN = mistNet_trapping_sim[scenario])
    
    mistNet_trapping_list[[scenario]] <- res[[1]]
    mistNet_trapping_out[scenario,] <- colMeans(mistNet_trapping_list[[scenario]])
    
  }
  
  mistNet_trapping_fig <- mistNet_trapping_out[,11:21]
  
  png(filename = paste0(file.path(res_fold,
                        "mist-net-trapping_scenario_suppCont.png")))
  
  plot(mistNet_trapping_fig[1,], 
       xaxt = "n", 
       xaxs = "r",
       yaxs = "i",
       type = "l",lwd = 1,
       lty = 1,
       xlab = expression(bold("")),
       ylab = "Total population",
       ylim = c(0,4000),
       main = "Trapped with mist-nets, random take of juvs and adults,\nsupplementation continues (n = 59/year)")
  axis(1, 
       at = c(1,6,11),
       labels = c("2022","2027","2032"))
  for(i in 1:length(mistNet_trapping_fig[,1])){
    lines(mistNet_trapping_fig[i,], type = "l", lwd = 1, lty = 1)
  }
    dev.off()
  
  #### Scenario 5. Mist netting + supplementation stops ####
  mistNet_trapping_sim <- c(40,80,120)
  mistNet_trapping_noSupp_list <- list()
  mistNet_trapping_noSupp_out <- matrix(0, 
                                        nrow = length(mistNet_trapping_sim),        
                                        ncol = length(mod1_baseline[[1]][1,]))
  
  for(scenario in 1:length(mistNet_trapping_sim)) {
    
    res <- balimynaMPM(trapping = TRUE,
                       trappingType = "mist-net",
                       trappedN = mistNet_trapping_sim[scenario],
                       future.supplementation = 0)
    
    mistNet_trapping_noSupp_list[[scenario]] <- res[[1]]
    mistNet_trapping_noSupp_out[scenario,] <- colMeans(mistNet_trapping_noSupp_list[[scenario]])
    
  }
  
  mistNet_trapping_noSupp_fig <- mistNet_trapping_noSupp_out[,11:21]
  
  png(filename = paste0(file.path(res_fold,
                        "mist-net-trapping_scenario_suppStops.png")))
  
  plot(mistNet_trapping_noSupp_fig[1,], 
       xaxt = "n", 
       xaxs = "r",
       yaxs = "i",
       type = "l",lwd = 1,
       lty = 1,
       xlab = expression(bold("")),
       ylab = "Total population",
       ylim = c(0,4000),
       main = "Trapped with mist-nets, random take juvs and adults,\nno supplementation")
  axis(1, 
       at = c(1,6,11),
       labels = c("2022","2027","2032"))
  for(i in 2:length(mistNet_trapping_noSupp_fig[,1])){
    lines(mistNet_trapping_noSupp_fig[i,], type = "l", lwd = 1, lty = 1)
  }
    dev.off()
  
#### Scenario 6. Mist netting + supplementation stops after 5 years ####
    mistNet_trapping_sim <- c(40,80,120)
    mistNet_trapping_suppStops5_list <- list()
    mistNet_trapping_suppStops5_out <- matrix(0, 
                                          nrow = length(mistNet_trapping_sim),
                                          ncol = length(mod1_baseline[[1]][1,]))
    
    for(scenario in 1:length(mistNet_trapping_sim)) {
      
      res <- balimynaMPM(trapping = TRUE,
                         trappingType = "mist-net",
                         trappedN = mistNet_trapping_sim[scenario],
                         delayedSupplementation = TRUE)
      
      mistNet_trapping_suppStops5_list[[scenario]] <- res[[1]]
      mistNet_trapping_suppStops5_out[scenario,] <- colMeans(mistNet_trapping_suppStops5_list[[scenario]])
      
    }
    
    mistNet_trapping_suppStops5_fig <- mistNet_trapping_suppStops5_out[,11:21]
    
    png(filename = paste0(file.path(res_fold,
                          "mist-net-trapping_scenario_suppStops5.png")))
    
    plot(mistNet_trapping_suppStops5_fig[1,], 
         xaxt = "n", 
         xaxs = "r",
         yaxs = "i",
         type = "l",lwd = 1,
         lty = 1,
         xlab = expression(bold("")),
         ylab = "Total population",
         ylim = c(0,4000),
         main = "Trapped with mist-nets, random take juvs and adults,\nsupplementation ceases after 5 years")
    axis(1, 
         at = c(1,6,11),
         labels = c("2022","2027","2032"))
    for(i in 2:length(mistNet_trapping_suppStops5_fig[,1])){
      lines(mistNet_trapping_suppStops5_fig[i,], type = "l", lwd = 1, lty = 1)
    }
    dev.off()
    
####-------------------------------------------####
#### TABLE 2 summary stats to show mean population and risk of population being below 1000 individuals after ten years 
####-------------------------------------------####
  # Number of individuals after 10 years with nest poaching at 1 = low, 2 = medium and 3 = high levels 
  # supp continues
  mean(nest_trapping_list[[1]][,"2032"]) # 2689
  mean(nest_trapping_list[[2]][,"2032"]) # 2237
  mean(nest_trapping_list[[3]][,"2032"]) # 1771
  #supp stops after 5 years
  mean(nest_trapping_supp5_list[[1]][,"2032"]) # 2202
  mean(nest_trapping_supp5_list[[2]][,"2032"]) # 1738
  mean(nest_trapping_supp5_list[[3]][,"2032"]) # 1300
  
  # chance of population remaining under 1,000 individuals
  # supp continues
  table(nest_trapping_list[[1]][,"2032"] < 1000) # <1%
  table(nest_trapping_list[[2]][,"2032"] < 1000) # 2.3%
  table(nest_trapping_list[[3]][,"2032"] < 1000) # 10.9%
  # supp stops after 5 years
  table(nest_trapping_supp5_list[[1]][,"2032"] < 1000) # 3.7%
  table(nest_trapping_supp5_list[[2]][,"2032"] < 1000) # 17.9%
  table(nest_trapping_supp5_list[[3]][,"2032"] < 1000) # 38.8%
  
  # Number of individuals after 10 years with mist-netting at 1 = low, 2 = medium and 3 = high levels 
  
  # supp continues
  mean(mistNet_trapping_list[[1]][,"2032"]) # 2440
  mean(mistNet_trapping_list[[2]][,"2032"]) # 1744
  mean(mistNet_trapping_list[[3]][,"2032"]) # 1066
  #supp stops after 5 years
  mean(mistNet_trapping_suppStops5_list[[1]][,"2032"]) # 1951
  mean(mistNet_trapping_suppStops5_list[[2]][,"2032"]) # 1254
  mean(mistNet_trapping_suppStops5_list[[3]][,"2032"]) # 624
  
  # chance of population remaining under 1,000 individuals
  # supp continues
  table(mistNet_trapping_list[[1]][,"2032"] < 1000) # <1%
  table(mistNet_trapping_list[[2]][,"2032"] < 1000) # 12.7%
  table(mistNet_trapping_list[[3]][,"2032"] < 1000) # 54.2%
  # supp stops after 5 years
  table(mistNet_trapping_suppStops5_list[[1]][,"2032"] < 1000) # 8.8%
  table(mistNet_trapping_suppStops5_list[[2]][,"2032"] < 1000) # 40.6%
  table(mistNet_trapping_suppStops5_list[[3]][,"2032"] < 1000) # 78.7%
  
#####-------------------#####
#### FIGURE 3: Plot trapping scenarios side by side. Columns are trapping type, rows are trapping intensity, and the trend and 5/95 percentiles are shown for supplemented and non-supplemented population in the same plot. Result should be a 3 by 2 grid of plots
#####-------------------#####

  # function to create the dataframes needed to generate plots
  filterSim <- function(name1,
                        name2,
                        name3,
                        list.element = 1, 
                        conf.width = 0.025
  ){
    
    # establish continued supplementation DF
    res1 <- name1[[list.element]]
    res1 <- res1[,11:length(res1[1,])]
    
    df1 <- data.frame(
      meanPop = colMeans(res1),
      lower = apply(res1,
                    MARGIN = 2,
                    function(x) quantile(x,
                                         probs = conf.width)),
      upper = apply(res1,
                    MARGIN = 2,
                    function(x) quantile(x,
                                         probs = 1-conf.width)),
      Year = 2022:2032,
      supp = "Supplementation continues"
    )
    
    # establish supplementation stops DF
    res2 <- name2[[list.element]]
    res2 <- res2[,11:length(res2[1,])]
    
    df2 <- data.frame(
      meanPop = colMeans(res2),
      lower = apply(res2,
                    MARGIN = 2,
                    function(x) quantile(x, probs = conf.width)),
      upper = apply(res2,
                    MARGIN = 2,
                    function(x) quantile(x, probs = 1-conf.width)),
      Year = 2022:2032,
      supp = "Supplementation ceases after 2021"
    )
    
    # establish supplementation stops DF
    res3 <- name3[[list.element]]
    res3 <- res3[,11:length(res3[1,])]
    
    df3 <- data.frame(
      meanPop = colMeans(res3),
      lower = apply(res3,
                    MARGIN = 2,
                    function(x) quantile(x, probs = conf.width)),
      upper = apply(res3,
                    MARGIN = 2,
                    function(x) quantile(x, probs = 1-conf.width)),
      Year = 2022:2032,
      supp = "Supplementation ceases after 2026"
    )
    
    b <- list(df1,df2,df3)
    df <- do.call(rbind, b)
    
    return(df)
  }
  
  ## generate dataframes for plotting
  nest_low <- filterSim(name1 = nest_trapping_list,
                        name2 = nest_trapping_noSupp_list,
                        name3 = nest_trapping_supp5_list,
                        list.element = 1,
                        conf.width = 0.05)
  mist_low <- filterSim(name1 = mistNet_trapping_list,
                        name2 = mistNet_trapping_noSupp_list,
                        name3 = mistNet_trapping_suppStops5_list,
                        list.element = 1,
                        conf.width = 0.05)
  nest_mid <- filterSim(name1 = nest_trapping_list,
                        name2 = nest_trapping_noSupp_list,
                        name3 = nest_trapping_supp5_list,
                        list.element = 2,
                        conf.width = 0.05)
  mist_med <- filterSim(name1 = mistNet_trapping_list,
                        name2 = mistNet_trapping_noSupp_list,
                        name3 = mistNet_trapping_suppStops5_list,
                        list.element = 2,
                        conf.width = 0.05)
  nest_high <- filterSim(name1 = nest_trapping_list,
                         name2 = nest_trapping_noSupp_list,
                         name3 = nest_trapping_supp5_list,
                         list.element = 3,
                         conf.width = 0.05)
  mist_high <- filterSim(name1 = mistNet_trapping_list,
                         name2 = mistNet_trapping_noSupp_list,
                         name3 = mistNet_trapping_suppStops5_list,
                         list.element = 3,
                         conf.width = 0.05)
  
  #### plot 1 ####
  p1 <-   ggplot() +
    geom_ribbon(data = nest_low[nest_low$supp == "Supplementation continues", ], 
                aes(x = Year, ymin = lower, ymax = upper, fill = supp, colour = supp),
                alpha = 0.4,
                linetype = 2,
                linewidth = 0.8) +
    
    geom_ribbon(data = nest_low[nest_low$supp == "Supplementation ceases after 2026", ],                
                aes(x = Year, ymin = lower, ymax = upper, fill = supp, colour = factor(supp)),
                alpha = 0.4,
                linetype = 2,
                linewidth = 0.8) +
    
    geom_line(data = nest_low[nest_low$supp == "Supplementation continues", ],
              aes(x = Year, y = meanPop, colour = supp),
              linewidth = 1.1) +
    geom_line(data = nest_low[nest_low$supp == "Supplementation ceases after 2026", ],
              aes(x = Year, y = meanPop, colour = supp),
              linewidth = 1.1) +
    
    
    scale_y_continuous(breaks = seq(1000,5000,1000),
                       expand = c(0,0),
                       name = "Total population") +
    scale_x_continuous(breaks = seq(2022,2032,1),
                       name = NULL,
                       expand = c(0,0)) +
    coord_cartesian(ylim = c(0,5000)) +
    ggtitle("Nest poaching",
            subtitle = expression("Offtake: 40 first-year birds year"^"–1")) +
    
    theme_bw() +
    theme(plot.margin = margin(0.4,0.5,0.2,0.2, "cm"),
          panel.grid.minor = element_blank(),
          axis.ticks.length=unit(.15, "cm"),
          legend.justification = c(0,1),
          legend.position = c(0.05,0.95),
          legend.title = element_blank(),
          legend.key.size = unit(0.5, "cm"),
          legend.text = element_text(size=12),
          axis.text = element_text(size = 10),
          axis.title.y = element_text(size = 12, 
                                      margin = margin(0,5,0,0)),
          title = element_text(size = 14),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(size = 12)) +
    
    scale_colour_manual(breaks = c("Supplementation continues", 
                                   "Supplementation ceases after 2026"),
                        values = c("#67a9cf","#ef8a62")) +
    
    scale_fill_manual(breaks = c("Supplementation continues", 
                                 "Supplementation ceases after 2026"),
                      values = c("#67a9cf","#ef8a62"),
                      labels = c("Supplementation continues",
                                 "Ceases after further 5 years")) +
    
    guides(colour = FALSE)
  
  plot(p1)
  #### end p1 
  
  #### plot 2 ####
  p2 <- ggplot() +
    geom_ribbon(data = mist_low[mist_low$supp == "Supplementation continues", ], 
                aes(x = Year, ymin = lower, ymax = upper, fill = supp, colour = supp),
                alpha = 0.4,
                linetype = 2,
                linewidth = 0.8) +
    geom_ribbon(data = mist_low[mist_low$supp == "Supplementation ceases after 2026", ], 
                aes(x = Year, ymin = lower, ymax = upper, fill = supp, colour = supp),
                alpha = 0.4,
                linetype = 2,
                linewidth = 0.8) +
    
    geom_line(data = mist_low[mist_low$supp == "Supplementation ceases after 2026", ],
              aes(x = Year, y = meanPop, colour = supp),
              linewidth = 1.1) +
    geom_line(data = mist_low[mist_low$supp == "Supplementation continues", ],
              aes(x = Year, y = meanPop, colour = supp),
              linewidth = 1.1) +
    
    scale_y_continuous(breaks = seq(1000,5000,1000),
                       expand = c(0,0),
                       name = NULL) +
    scale_x_continuous(breaks = seq(2022,2032,1),
                       name = NULL,
                       expand = c(0,0)) +
    coord_cartesian(ylim = c(0,5000)) +
    
    ggtitle("Mist-netting",
            subtitle = expression("Offtake: 40 birds year"^"–1"~"(random life stage taken)")) +
    
    theme_bw() +
    theme(plot.margin = margin(0.4,0.5,0.2,0.2, "cm"),
          axis.ticks.length=unit(.15, "cm"),
          legend.position = "none",
          panel.grid.minor = element_blank(),
          title = element_text(size = 14),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(size = 12)) +
    
    scale_colour_manual(breaks = c("Supplementation continues", 
                                   "Supplementation ceases after 2026"),
                        values = c("#67a9cf","#ef8a62")) +
    
    scale_fill_manual(breaks = c("Supplementation continues", 
                                 "Supplementation ceases after 2026"),
                      values = c("#67a9cf","#ef8a62"),
                      labels = c("Supplementation continues",
                                 "Ceases after further 5 years"))
  
  
  plot(p2)  
  #### end p2 
  
  #### plot 3 ####
  p3 <-   ggplot() +
    
    geom_ribbon(data = nest_mid[nest_mid$supp == "Supplementation continues", ], 
                aes(x = Year, ymin = lower, ymax = upper, fill = supp, colour = supp),
                alpha = 0.4,
                linetype = 2,
                linewidth = 0.8) +
    geom_ribbon(data = nest_mid[nest_mid$supp == "Supplementation ceases after 2026", ], 
                aes(x = Year, ymin = lower, ymax = upper, fill = supp, colour = supp),
                alpha = 0.4,
                linetype = 2,
                linewidth = 0.8) +
    geom_line(data = nest_mid[nest_mid$supp == "Supplementation ceases after 2026", ],
              aes(x = Year, y = meanPop, colour = supp),
              linewidth = 1.1) +
    geom_line(data = nest_mid[nest_mid$supp == "Supplementation continues", ],
              aes(x = Year, y = meanPop, colour = supp),
              linewidth = 1.1) +
    
    scale_y_continuous(breaks = seq(1000,5000,1000),
                       expand = c(0,0),
                       name = "Total population") +
    scale_x_continuous(breaks = seq(2022,2032,1),
                       name = NULL,
                       expand = c(0,0)) +
    coord_cartesian(ylim = c(0,5000)) +
    
    ggtitle(expression("Offtake: 80 first-year birds year"^"–1")) +
    
    theme_bw() +
    theme(plot.margin = margin(0.4,0.5,0.2,0.2, "cm"),
          axis.ticks.length=unit(.15, "cm"),
          legend.position = "none",
          axis.text = element_text(size = 10),
          axis.title.y = element_text(size = 12, 
                                      margin = margin(0,5,0,0)),
          plot.title = element_text(size = 12),
          panel.grid.minor = element_blank()) +
    
    scale_colour_manual(breaks = c("Supplementation continues", 
                                   "Supplementation ceases after 2026"),
                        values = c("#67a9cf","#ef8a62")) +
    
    scale_fill_manual(breaks = c("Supplementation continues", 
                                 "Supplementation ceases after 2026"),
                      values = c("#67a9cf","#ef8a62"),
                      labels = c("Supplementation continues",
                                 "Ceases after further 5 years"))
  plot(p3)
  ##### end plot p3
  
  #### plot 4 #####
  p4 <- ggplot() +
    geom_ribbon(data = mist_med[mist_med$supp == "Supplementation continues", ], 
                aes(x = Year, ymin = lower, ymax = upper, fill = supp, colour = supp),
                alpha = 0.4,
                linetype = 2,
                linewidth = 0.8) +
    geom_ribbon(data = mist_med[mist_med$supp == "Supplementation ceases after 2026", ], 
                aes(x = Year, ymin = lower, ymax = upper, fill = supp, colour = supp),
                alpha = 0.4,
                linetype = 2,
                linewidth = 0.8) +
    geom_line(data = mist_med[mist_med$supp == "Supplementation ceases after 2026", ],
              aes(x = Year, y = meanPop, colour = supp),
              linewidth = 1.1) +
    geom_line(data = mist_med[mist_med$supp == "Supplementation continues", ],
              aes(x = Year, y = meanPop, colour = supp),
              linewidth = 1.1) +
    
    scale_y_continuous(breaks = seq(1000,5000,1000),
                       expand = c(0,0),
                       name = NULL) +
    scale_x_continuous(breaks = seq(2022,2032,1),
                       name = NULL,
                       expand = c(0,0)) +
    coord_cartesian(ylim = c(0,5000)) +
    
    ggtitle(expression("Offtake: 80 birds year"^"–1"~"(random life stage taken)")) +
    
    theme_bw() +
    theme(plot.margin = margin(0.4,0.5,0.2,0.2, "cm"),
          axis.ticks.length=unit(.15, "cm"),
          legend.position = "none",
          panel.grid.minor = element_blank(),
          plot.title = element_text(size = 12)) +
    
    scale_colour_manual(breaks = c("Supplementation continues", 
                                   "Supplementation ceases after 2026"),
                        values = c("#67a9cf","#ef8a62")) +
    
    scale_fill_manual(breaks = c("Supplementation continues", 
                                 "Supplementation ceases after 2026"),
                      values = c("#67a9cf","#ef8a62"),
                      labels = c("Supplementation continues",
                                 "Ceases after further 5 years"))
  plot(p4)
  
  ##### end plot p4
  
  #### plot 5 ####
  p5 <-   ggplot() +
    
    geom_ribbon(data = nest_high[nest_high$supp == "Supplementation continues", ], 
                aes(x = Year, ymin = lower, ymax = upper, fill = supp, colour = supp),
                alpha = 0.4,
                linetype = 2,
                linewidth = 0.8) +
    geom_ribbon(data = nest_high[nest_high$supp == "Supplementation ceases after 2026", ], 
                aes(x = Year, ymin = lower, ymax = upper, fill = supp, colour = supp),
                alpha = 0.4,
                linetype = 2,
                linewidth = 0.8) +
    geom_line(data = nest_high[nest_high$supp == "Supplementation ceases after 2026", ],
              aes(x = Year, y = meanPop, colour = supp),
              linewidth = 1.1) +
    geom_line(data = nest_high[nest_high$supp == "Supplementation continues", ],
              aes(x = Year, y = meanPop, colour = supp),
              linewidth = 1.1) +
    
    scale_y_continuous(breaks = seq(1000,5000,1000),
                       expand = c(0,0),
                       name = "Total population") +
    scale_x_continuous(breaks = seq(2022,2032,1),
                       name = "Year",
                       expand = c(0,0)) +
    coord_cartesian(ylim = c(0,5000)) +
    
    ggtitle(expression("Offtake: 120 first-year birds year"^"–1")) +
    
    theme_bw() +
    theme(plot.margin = margin(0.4,0.5,0.2,0.2, "cm"),
          panel.grid.minor = element_blank(),
          axis.ticks.length=unit(.15, "cm"),
          legend.position = "none",
          axis.text = element_text(size = 10),
          axis.title.y = element_text(size = 12, 
                                      margin = margin(0,5,0,0)),
          axis.title.x = element_text(size = 12,
                                      margin = margin(5,0,0,0)),
          plot.title = element_text(size = 12)) +
    
    scale_colour_manual(breaks = c("Supplementation continues", 
                                   "Supplementation ceases after 2026"),
                        values = c("#67a9cf","#ef8a62")) +
    
    scale_fill_manual(breaks = c("Supplementation continues", 
                                 "Supplementation ceases after 2026"),
                      values = c("#67a9cf","#ef8a62"),
                      labels = c("Supplementation continues",
                                 "Ceases after further 5 years"))
  plot(p5)
  #### end p5
  
  #### plot 6 ####
  p6 <- ggplot() +
    geom_ribbon(data = mist_high[mist_high$supp == "Supplementation continues", ], 
                aes(x = Year, ymin = lower, ymax = upper, fill = supp, colour = supp),
                alpha = 0.4,
                linetype = 2,
                linewidth = 0.8) +
    geom_ribbon(data = mist_high[mist_high$supp == "Supplementation ceases after 2026", ], 
                aes(x = Year, ymin = lower, ymax = upper, fill = supp, colour = supp),
                alpha = 0.4,
                linetype = 2,
                linewidth = 0.8) +
    geom_line(data = mist_high[mist_high$supp == "Supplementation ceases after 2026", ],
              aes(x = Year, y = meanPop, colour = supp),
              linewidth = 1.1) +
    geom_line(data = mist_high[mist_high$supp == "Supplementation continues", ],
              aes(x = Year, y = meanPop, colour = supp),
              linewidth = 1.1) +
    
    scale_y_continuous(breaks = seq(1000,5000,1000),
                       expand = c(0,0),
                       name = NULL) +
    scale_x_continuous(breaks = seq(2022,2032,1),
                       name = "Year",
                       expand = c(0,0)) +
    coord_cartesian(ylim = c(0,5000)) +
    
    ggtitle(expression("Offtake: 120 birds year"^"–1"~"(random life stage taken)")) +
    
    theme_bw() +
    theme(plot.margin = margin(0.4,0.5,0.2,0.2, "cm"),
          axis.ticks.length=unit(.15, "cm"),
          legend.position = "none",
          panel.grid.minor = element_blank(),
          axis.title.x = element_text(size = 12,
                                      margin = margin(5,0,0,0)),
          plot.title = element_text(size = 12)) +
    
    scale_colour_manual(breaks = c("Supplementation continues", 
                                   "Supplementation ceases after 2026"),
                        values = c("#67a9cf","#ef8a62")) +
    
    scale_fill_manual(breaks = c("Supplementation continues", 
                                 "Supplementation ceases after 2026"),
                      values = c("#67a9cf","#ef8a62"),
                      labels = c("Supplementation continues",
                                 "Ceases after further 5 years"))
  plot(p6)
  
  #### arrange plots and save ####
  Fig3 <- ggarrange(p1,p2,p3,p4,p5,p6,
                    ncol = 2)
  plot(Fig3)
  
  ggsave(paste0(file.path(res_fold,"Fig3_scenarios.png")),
         plot = Fig3,
         width = 9, height = 12, units = "in" 
  )
  
