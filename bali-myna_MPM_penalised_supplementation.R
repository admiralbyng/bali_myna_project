########################################################################
# PROJECT: Bali myna population viability analysis
# CONTENTS: 
#  - building a matrix population model for validation and projection of future scenarios of population demographics
#  DEPENDENCIES:
#  - Code documents
#  - Data files
# AUTHOR: tom squires
########################################################################

# PREAMBLE ===============================================================

# rm(list=ls())

## Directories ------------------------------------------------------------
### Define directories in relation to project directory
Dir.Base <- getwd()
Dir.Data <- file.path(Dir.Base, "data")
Dir.Exports <- file.path(Dir.Base, "exports")
### Create directories which aren't present yet
Dirs <- c(Dir.Data, Dir.Exports)
sapply(Dirs, function(x) if(!dir.exists(x)) dir.create(x))


## FUNCTION ==================================================

# implementing a stage-based female-only matrix population model (MPM) to simulate the Bali Myna population trajectory over the modelling period for different scenarios

balimynaMPM <- function(
  # default settings for what will be the 'baseline' Bali myna demographic model
  
  lambda = 1.374, 
  # = initial estimate of population increase i.e. 10 root of (N in 2021 / N in 2012) = (420/15)^(1/10)
  
  stages = c("Juvenile","Adult"), # life stages
  T1 = 1, # length of lifestage 1 (immature)
  T2 = 15, # length of lifestage 2 (adult) max lifespan in Ross et al. 2021 "Population Analysis & Breeding and Transfer Plan: Bali Myna (Leucopsar rothschildi) AZA Species Survival Plan® Green Program"
  firstYear = 2012, # first year of modelling period
  historicalYears = 10, # number of years with historical data
  projectionYears = 10, # number of years to project into future
  nJuvsInit = 0, # initial number of immature birds in the population at t = 0
  nAdultsInit = 15, # initial number of adults in the population at t = 0
  nsimul = 1000, # number of times to run the simulation
  
  Sjuv_mean = 0.55, # mean immature survival
  Sjuv_sd = 0.05, # variation in immature survival
  Sadult_mean = 0.8, # mean adult survival
  Sadult_sd = 0.05, # variation in adult survival
  
  female.ratio = 0.5, # birth sex ratio
  female.breeders_mean = 0.61, # mean proportion of females breeding based on Yuni et al 2022 breeding females = 30, BBNP 2021 census = 256 with the 108 birds released since Oct 2018 removed = 148. 148 * 0.66 (R. Martin's 2/3 rule for adult/immature birds in IUCN Red List calcs) = 97.6 * 0.5 (sex ratio) = 48.8. 30/48.8 females breeding = 0.61
  female.breeders_sd = 0.1, # variation in proportion of females breeding
  F.penalty.juv = 0.7, # penalty to first year birds during first breeding season because their fecundity is probably lower than adults
  
  chicksPerAttempt = 1.6,
  attemptsMin = 1,
  attemptsMax = 4,
  attemptsMean = 1.6,
  attemptsSD = 0.9,
  # Fmean = 1.6, # mean number of chicks per brood
  # Fsd = 1.1, # variation in number of chicks per brood
  supplementation = ceiling(c(10,24,14,12,22,28,42,76,80,68)), # known values from BBNP (2021)
  future.supplementation = 59, # continuing supplementation at the average number of birds released over previous five years
  supplementation_penalty_S = 0.3,  # penalty to adult survival of birds released previous yr; lower values = higher pen
  supplementation_penalty_F = 0.3, # penalty to fecundity of birds released previous yr; lower values = higher pen
  delayedSupplementation = FALSE,
  supplementationDelay = 5,
  
  densityDependence = FALSE, # if TRUE,  nest site limitation is used to impose D-D
  N_nestSites = 150, # current number of nest sites, used as baseline from which to increase the nest sites using the growth rate
  nestSiteGrowthRate = 1.05, # rate of nest site increase per year
  max_nest_sites = 600, # arbitrary maximum number of nest sites estimated for the study region to be able to put a density dependent term in the model
  
  # trapping scenario
  trapping = FALSE,
  trappingType = "nest", # options = "nest" or "mist-net"
  trappedN = 40 # total number trapped (males + females)
)
{
  
  #####--------------------------------#####  
  ## load packages ---------------------------------------------------------------
  #####--------------------------------#####  
  install.load.package <- function(x){
    if(!require(x, character.only = TRUE)) 
      install.packages(x, repos = 'http://cran.us.r-project.org')
    require(x, character.only = TRUE)
  }
  
  package_vector <- c("popbio","mvtnorm","msm")
  sapply(package_vector, install.load.package)
  
  #####--------------------------------#####  
  ## store some of the values for parameters based on the selected model settings
  #####--------------------------------#####  
  
  ## define length of modelling period
  yearsModelled <- historicalYears + projectionYears
  lastYear <- firstYear + yearsModelled
  
  ## create vector of number of nest sites available in each year
  nest_sites <- round(N_nestSites * nestSiteGrowthRate^(0:yearsModelled), digits = 0) # number of nest sites increases each year by x% (installation of more nest boxes)
  
  ## create empty matrix to be filled during each ITERATION of the simulation (n = 1000)
  N <- matrix(0, nrow = length(stages), ncol = yearsModelled + 1) 
  rownames(N) <- stages # name the rows
  colnames(N) <- firstYear:lastYear # name the cols with actual years
  
  # the first year of the model contains the known population figures
  N[,1] <- ceiling(c(nJuvsInit,nAdultsInit) * female.ratio) 
  
  ## create a matrix to contain the population size in each year and run of simulation
  NTot <- matrix(0, nrow = nsimul, ncol = yearsModelled + 1)
  rownames(NTot) <- 1:nsimul # name rows
  colnames(NTot) <- firstYear:lastYear # name cols with actual years
  
  ## create a matrix to contain the proportion of birds that were immature in each year
  propJuv <- matrix(0, nrow = nsimul, ncol = yearsModelled + 1)
  rownames(propJuv) <- 1:nsimul # name rows
  colnames(propJuv) <- firstYear:lastYear # name cols with actual years
  
  ## numbers supplemented needs to be converted into females only:
  
  supplementationFemales <- supplementation * female.ratio
  
  future.supplementationFemales <- future.supplementation * female.ratio
  
  #####--------------------------------#####  
  
  ## FUNCTIONS ----------------------------------------------------------------
  
  #####--------------------------------#####  
  
  
  ## immature survival (growth into adult stage) ----------------------------
  
  juvenile.survival.growth <- function(){
    
    # immature survival = 0 because this life stage lasts for t = 1 (they either become adults or die)
    s1 <- round(rtnorm(n = 1, mean = Sjuv_mean, sd = Sjuv_sd,lower = 0, upper = 1), 
                digits = 3) 
    # generate a random annual survival rate from the survival rate for immatures (mean and sd)
    
    p1 <- round(s1 * (1 - ((((s1/lambda)^T1) - ((s1/lambda)^(T1-1))) / (((s1/lambda)^T1)-1))), digits = 3)
    
    g1 <- round(s1 * ((((s1/lambda) ^T1) - ((s1/lambda) ^(T1-1))) / (((s1/lambda)^T1)-1)),
                digits = 3)
    
    result <- matrix(c(s1, p1, g1), 
                     nrow = 1, 
                     dimnames = list("",c("s1","p1","g1")))
  }
  
  ## adult survival ---------------------------------------------------------
  
  adult.survival <- function(){ 
    
    s2 <- round(rtnorm(n = 1, mean = Sadult_mean, sd = Sadult_sd, lower = 0, upper = 1), digits = 3)
    
    p2 <- round(s2 * (1 - ((((s2 / lambda) ^ T2) - ((s2 / lambda) ^ (T2 - 1))) / (((s2 / lambda) ^  T2) - 1))), digits = 3) # this is a correction to the adult survival rate to account for the loss of individuals that reach the maximum life expectancy of the species
    
    # modify survival for adults that were released in t - 1
    p2_released <- round(p2 * supplementation_penalty_S, digits = 3) # penalty to survival of last year's released birds
    
    p2_adjusted <- (p2_released * prop.released) + (p2 * (1 - prop.released)) # combine the survival rates for wild and released adults based on relative proportions
  }
  
  
  ## adult fecundity -------------------------------------------------------
 
  fecundity <- function(){ 
    
    ### not the simplest function. There are a few important points here:
    # 1. number of attempts is drawn for the number of breeding females, which is then multiplied by the fixed parameter of chicks per attempt taken from our other paper (Yuni et al 2022)
    # 2. our census data are pre-breeding. So, birds with age 0 in the current timestep do not breed until the next timestep, because the next breeding season follows the next census. This means immature females should be excluded from fecundity calculations, unlike if it was a post-breeding census. 
    # 3. as it is a pre-breeding census, all of our fecundity rates should NOT be adjusted down to account for mortality (birds breed immediately after census - presumed that there is almost no mortality). Fecundity rates should NOT be subjected to stage-specific survival rate. 
    # 4. it is likely that birds breeding for the first time are not as productive as full adults, partly because a lot of the 1st year males may not even breed. So a penalty should be applied to the fecundity rate of 1st year birds (in addition to the survival rate penalty)
    
    if(female.breeders_mean < 1){
      female_breeders <- rtnorm(n = 1, 
                                mean = female.breeders_mean,
                                sd = female.breeders_sd, 
                                lower = 0, 
                                upper = 1) # randomly select % of females that breed
    } else{female_breeders <- 1}
    
    N_female_breeders <- round((N_adults_female + N_juvs_female) * female_breeders, digits = 0)
    
    if(N_female_breeders > 0){
      
      N_attempts <- sum(floor(rtnorm(n = N_female_breeders, 
                                     mean = attemptsMean, 
                                     sd = attemptsSD, 
                                     lower = attemptsMin, 
                                     upper = attemptsMax))) # total number of breeding attempts
    } else{N_attempts <- 0}
    
    N_chicks <- N_attempts * chicksPerAttempt * female.ratio # female chicks per attempt is fixed based on our productivity data from 2019
    
    if(N_chicks > 0 & N_female_breeders > 0){
      
    f2_basic <- N_chicks / N_female_breeders
      
# fecundity of released birds is penalised and subject to the survival rate of first year birds (PREBREEDING CENSUS - chicks must survive the year first)
    f2_released <- f2_basic * supplementation_penalty_F * g1
      
# fecundity of all 'wild' birds is boosted because our baseline data come from a supplemented population
    f2_boosted  <- f2_basic * 1.5 * g1 ### PREBREEDING so include first year survival in the fertility coefficient
    
    f2_year1    <- f2_basic * 1.5 * F.penalty.juv * g1 # PREBREEDING so include IMMATURE survival rate
      
    f2_adjusted <- (f2_released * prop.released) + (f2_year1 * prop.juv) + (f2_boosted * prop.wildAdFemales)
        
    f2 <- f2_adjusted 
    } else{
        
      f2 <- 0
    }

      if(densityDependence == TRUE){
      
      if(N_female_breeders > max_nest_sites) {N_female_breeders <- max_nest_sites} # number of breeding females cannot be higher than the max number of nest sites
      
      dd_f2 <- f2 * (1 - (N_female_breeders/max_nest_sites) ^ 10)
      # if you decrease the power here, the density dependent effect is strengthened, i.e. females are affected by the lack of nest sites sooner (due to increased intraspecific competition for nest sites - a likely occurrence for Bali Myna)
      
      f2 <- dd_f2
    }
    
    return(f2)
  }
  
  
  ## trapping scenario ---------------------------------------------------------
  
  trappingScenario <- function(){
    
    # one of these two options should be used, otherwise throw an error
    if(trappingType != "mist-net" & trappingType != "nest"){
      stop("trappingType should be either nest or mist-net")
    }
    
    # convert trappedN into females only (female-only model)
    trappedFemales <- trappedN * female.ratio
    
    # calculate number of immature and adult females to remove from the population. There is a check in here to make sure the population can't go below zero  
    
    if(trappingType == "mist-net"){
      
      # for mist-net trapping, a random proportion of birds trapped are assumed to be immature, and it changes each year
      proportionJuv <- round(runif(1, min = 0.2, max = 0.8), digits = 1)
      
      # number of immature birds trapped
      nJuvsTrapped <- round(trappedFemales * proportionJuv, digits = 0)
      nJuvs <- if(nJuvsTrapped > N_juvs_female){N_juvs_female}else{nJuvsTrapped}
      
      #number of adults trapped
      nAdultsTrapped <- round(trappedFemales * (1 - proportionJuv), digits = 1)
      
      nAdults <- if(nAdultsTrapped > N_adults_female) {
        N_adults_female
      } else{
        nAdultsTrapped
      }
      
      result <- c(nJuvs, nAdults)
    }
    
    if(trappingType == "nest"){
      nJuvs <- if(N_juvs_female < trappedFemales) {N_juvs_female
      } else{ 
        trappedFemales
      }
      nAdults <- 0
      result <- c(nJuvs, nAdults)
    }
    
    return(result)
  }
  
  ## run model using FOR loops =================================================
  
  for (j in 1:nsimul) { 
    
    for (t in 1:historicalYears) {
      # this loop is for the first 10 historical years (2012-2021) of data taken from BBNP's 2021 Bali Myna report
      
      # establish the number of juvs, females, and released females in last year
    N_juvs_female   <- N[1, t] # number of immature birds at the start of this year
    N_adults_female <- N[2, t] # number of adults at the start of this year
    N_released      <- if(t > 1) {supplementationFemales[t-1]} else {N_adults_female * 0.5}       # for the first year of modelling half of the adult females are considered released individuals (actually, in 2011 4 adults were released, so that is about right )
    
    # establish proportion of all birds that are each type:
prop.wildAdFemales <- (N_adults_female - N_released) / (N_adults_female + N_juvs_female)
prop.released <- N_released/(N_adults_female + N_juvs_female)  
prop.juv <- N_juvs_female/(N_juvs_female + N_adults_female)
      
      ## immature survival and growth into adult life stage --------------------
      
      juvenile_survival_metrics <- juvenile.survival.growth()
      
      p1 <- juvenile_survival_metrics[1, "p1"]    
      g1 <- juvenile_survival_metrics[1, "g1"]  
      
      ## adult survival --------------------------------------------------------
      
      p2 <- adult.survival()
      
      ## fecundity ------------------------------------------------------------
      
      f2 <- fecundity()
      
      ## construct the transition matrix for period t --------------------------
      
      TMat <- matrix(c(p1, g1, f2, p2), nrow = 2, ncol = 2) # create the transition matrix
      rownames(TMat) <- stages
      colnames(TMat) <- stages
      
      # matrix of supplementation
      supMat <- matrix(0, nrow = 2, ncol = 1) # matrix to store number of birds supplemented
      rownames(supMat) <- stages
      
      supMat[,1] <- c(0,supplementationFemales[t]) 
      
      N[ ,t+1] <- (TMat %*% N[,t]) + supMat[ ,1]
      
      NTot[j, t+1] <- round(sum(N[, t+1]) / female.ratio, digits=0) # convert the number of females back to total number of adults based on female:male ratio
      
      propJuv[j, t+1] <- round(N[1, t+1] / sum(N[ , t+1]), digits=2)
      
    } 
    
    for(t in (historicalYears + 1):yearsModelled){
      
      # future years of the modelling period
      
      # establish the number of immature, adult, and released females were in the modelled population in the last year
      N_juvs_female <- N[1, t] # number of immature birds at the start of this year
      N_adults_female <- N[2, t] # number of adult birds at the start of this year
      N_released <- future.supplementationFemales
      
      prop.wildAdFemales <- (N_adults_female - N_released) / (N_adults_female + N_juvs_female)
      prop.released <- N_released/N_adults_female # proportion of all birds that are released
      prop.juv <- N_juvs_female / (N_juvs_female + N_adults_female)
      
      ## immature metrics -----------------------------------------------------
      juvenile_survival_metrics <- juvenile.survival.growth()
      
      p1 <- juvenile_survival_metrics[1, "p1"]    
      g1 <- juvenile_survival_metrics[1, "g1"]  
      
      ## adult survival -------------------------------------------------------
      
      s2 <- adult.survival()
      
      ## fecundity ------------------------------------------------------------
      
      f2 <- fecundity()
      
      ## trapped birds
      trapped_birds <- if(trapping == TRUE){trappingScenario()} else{c(0,0)}
      
      
      ## construct the transition matrix for period t ------------------------- 
      
      TMat <- matrix(c(p1, g1, f2, p2), nrow = 2, ncol = 2) 
      
      rownames(TMat) <- stages
      colnames(TMat) <- stages
      
      # add supplemented birds and include delay if the option is selected in function input
      if(delayedSupplementation == TRUE){
        if(t <= historicalYears + supplementationDelay){
          futureSup <- matrix(c(0, future.supplementationFemales), nrow = 2)
        } else{
          futureSup <- matrix(c(0, 0), nrow = 2)
        }
      } else{
        
        futureSup <- matrix(c(0, future.supplementationFemales), nrow = 2)
      }
      ## add year of data to N and then to NTot ----------------------------------
      
      N[, t + 1] <- (TMat %*% N[ , t]) + futureSup[,1] - trapped_birds
      
      NTot[j, t+1] <- round(sum(N[ ,t+1]) / female.ratio, digits = 0)
      
      propJuv[j, t+1] <- round(N[1, t + 1] / sum(N[ , t + 1]), digits = 2)
    }
    
    # this does the sum to fill in the first year because that didn't get covered by the loop
    NTot[j, 1] <- sum(N[,1]) / female.ratio
    propJuv[j, 1] <- round(N[1,1] / sum(N[, 1]), digits = 2)
    
  }
  
  result <- list(NTot, propJuv)
  names(result)[1] <- "NTot"
  names(result)[2] <- "proportionJuveniles"
  return(result)
}


