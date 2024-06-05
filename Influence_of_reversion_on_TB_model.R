#----------------------------------------------------------------------------------------------------#
### Exploring the influence of reversion on a simple dynamic tuberculosis (TB) compartmental model ###
#----------------------------------------------------------------------------------------------------#
#' Model assumptions (among others)
#'  - Stable population 
#'  - No Age structure
#'  - No risk groups (e.g. people living with HIV)
#'  - No drug resistant TB transmission
#'  
#'  Model compartment definitions:
#'    - S - susceptible - exposure naive
#'    - SPE - susceptible - prior exposure
#'    - L - latent/immunoreactive
#'    - I - active TB disease
#'    - R - recovered
#'  
# Author: Katie Dale 
#'        adapted from https://github.com/juanvesga/TBpractical/tree/master, with thanks. 
#'        The model structure is roughly based on Model E from Menzies et al 2018
# Date: August 2023
#----------------------------------------------------------------------------------------------------#
#                                   #### Load packages ####
#----------------------------------------------------------------------------------------------------#
library(data.table)
library(deSolve)
library(gridExtra)
library(ggplot2)
library(viridis)

#----------------------------------------------------------------------------------------------------#
#                        #### Set objects and model parameters ####
#----------------------------------------------------------------------------------------------------#
# Define the path for saving output, if you wish.
strOutputPath <- "???"

# Define the main parameters
intStartYear  <- 2019
numBirthR     <- 0.017     # Birth rate (the overall death rate is the same, after taking TB 
                           # mortality into account, i.e. the population is stable)
intBeta       <- 10        # Per capita rate at which two specific individuals come into 
                           # effective contact per unit time (per year).
intInfDurYrs  <- 2         # Duration of infectious period for untreated TB, in years
numFastP      <- 0.065     # Fraction of infected that fast-progress to active TB
numReactR     <- 0.001     # Reactivation rate from latent/immunoreactive into active disease
numSusReactR  <- 0.00005   # Reactivation rate from susceptible - previously exposed to active disease 
numTBMortP    <- 0.4       # TB mortality probability, untreated
numTBMortTxP  <- 0.03      # TB mortality probability, treated
numImm        <- 0.79      # Protection against reinfection for those in L compartment
numRelaps     <- 0.002857  # Relapse rate per year for those with prior untreated TB
numRelapsTx   <- 0.000525  # Relapse rate per year for those with prior treated TB
numReverP     <- 0.00      # Reversion probability

# Altered parameters for the models with reversion
numReverP2 <- 0.02         # Reversion probability
numReverP3 <- 0.04         # Reversion probability, model 2

# Parameters associated with the 'passive treatment' intervention
numAccessCareP    <- 0.70   # Probability of accessing care
numSympToTreatYrs <- 0.3    # Time delay (yrs) between developing symptoms and receiving treatment
numDxP            <- 0.95   # Probability of being diagnosed once sought care
numTxPAfterDxP     <- 0.95   # probability of receiving correct Tx if diagnosed
numTxPEffectYrs    <- 0.04   # Time until infectivity is eliminated, yrs 

numTxP <- numAccessCareP * numDxP * numTxPAfterDxP # the proportion treated

numTxPInfDur <- (intInfDurYrs * (1 - numTxP)) + # Duration of infectiousness for untreated
  (numTxP * (numSympToTreatYrs + numTxPEffectYrs)) # Duration of infectiousness for treated
numTxPRelaps <- (numRelaps * (1 - numTxP)) + # Relapse rate reduces for the treated
  (numTxP * numRelapsTx) 
numTxPMuTB <- (numTBMortP * (1 - numTxP)) + # TB mortality reduces for those treated
  (numTxP * numTBMortTxP) 

# Parameters associated with the 'some active screening' intervention
numAccessCareP     <- 0.9    # Probability of accessing care
numSympToTreatYrs <- 0.1    # Time delay (yrs) between developing symptoms and seeking for care

numTxP <- numAccessCareP * numDxP * numTxPAfterDxP # the proportion treated

numACFTxInfDur <- (intInfDurYrs * (1 - numTxP)) + # Duration of infectiousness for untreated
  (numTxP * (numSympToTreatYrs + numTxPEffectYrs)) # Duration of infectiousness for treated
numACFTxRelaps <- (numRelaps * (1 - numTxP)) + # Relapse rate reduces for the treated
  (numTxP * numRelapsTx) 
numACFTxMuTB <- (numTBMortP * (1 - numTxP)) + # TB mortality reduces for those treated
  (numTxP * numTBMortTxP) 

# Parameters associated with the 'latent Mtb infection screening and treatment' intervention
numScreenP         <- 0.6   # Probability of being screened
numTSTIGRASens     <- 0.7   # Sensitivity of screening tool 
numTxComplP        <- 0.8   # Likelihood of completing Latent Mtb treatment 
numReactBeforeTxP  <- 0.3   # Likelihood of reactivating before treatment
numTxEfficP        <- 0.69  # Latent Mtb infection treatment efficacy

numLatentTx <- numTSTIGRASens * numTxComplP * numTxEfficP * 
  numScreenP * (1 - numReactBeforeTxP) # the proportion treated

numReactRPrev  <- numReactR * (1 - numLatentTx) + 
  (numReactR * numTxEfficP) * numLatentTx

numFastPR <- numFastP * 0.95 # reduction in fast progression can't be much given
                             # so few receive treatment, and given that so few could be
                             # identified between exposure and reactivation anyway.

# Simulation descriptions
strName0 <- "Baseline - no interventions"
strName1 <- "Passive treatment"
strName4 <- "Passive treatment + TST/IGRA screening and latent Mtb infection treatment strategy"
strName2 <- "Passive treatment + active screening strategy"
strName5 <- "Hypothetical - no disease caused by remote infection (latent Mtb infection reactivation)"
strName3 <- "Hypothetical - no disease caused by recent infection"

lstLabels <- c(strName0, strName1,
               expression(paste("Passive treatment + TST/IGRA screening and latent",
                                italic(" Mtb "), "infection treatment strategy")),
               strName2,
               expression(paste("Hypothetical - no disease caused by remote infection (latent",
                                italic(" Mtb "), "infection reactivation)")),
               strName3)

lstSimNames <- c(strName0, strName1, strName4, 
                 strName2, strName5, strName3)

# Plot characteristics
numMaxYaxis <- 200
numTextSize <- 15
numLegendTextSize <- 15
numDist <- 6
numGeomTextSize <- 5
lstLegendPosition <- c(0.7, 0.8)
strTextColour <- "black"
numLineSize <- 1.5 
lstColourCode <- c("#90d743", "#35b779", "#21918c", 
                   "#31688e",  "#443983", "#440154")
names(lstColourCode) <- lstSimNames
lstLineStyleCode <- c("solid", "solid", "solid",
                      "solid", "dashed", "twodash")
names(lstLineStyleCode) <- lstSimNames

# Set simulation parameters 
intSimYears   <- 400      # years of simulation
intScaleYears <- 3        # Scaling up time of interventions
lstModelRun = seq(0, intSimYears, by = 1) # time scale

# Prepare list of parameters
lstParams <- c(numReactR = numReactR, numSusReactR = numSusReactR,
               numTBMortP = numTBMortP, intBeta = intBeta, numFastP = numFastP,
               numImm = numImm, numRelaps = numRelaps, numReverP = numReverP,
               intInfDurYrs = intInfDurYrs)

lstParamsR <- c(numReactR = numReactR, numSusReactR = numSusReactR,
               numTBMortP = numTBMortP, intBeta = intBeta, numFastP = numFastP,
               numImm = numImm, numRelaps = numRelaps, numReverP = numReverP2,
               intInfDurYrs = intInfDurYrs)

lstParamsRTwo <- c(numReactR = numReactR, numSusReactR = numSusReactR,
                numTBMortP = numTBMortP, intBeta = intBeta, numFastP = numFastP,
                numImm = numImm, numRelaps = numRelaps, numReverP = numReverP3,
                intInfDurYrs = intInfDurYrs)

# Initial Conditions
numTotalPop <- 1     # Total population equal to unity 
numZeroSeed <- 1e-6  # TB seed at time 0

dfStart <- data.frame(S = numTotalPop - numZeroSeed,
                      SPE = 0,
                      L = 0,
                      I = numZeroSeed,  
                      R = 0,
                      Incidence = 0, 
                      Irecent = 0, 
                      Iremote = 0, 
                      dRelap = 0, 
                      dTBmort = 0, 
                      dTBreact = 0, 
                      dTBreactSPE = 0)   

#----------------------------------------------------------------------------------------------------#
#                                      #### Set functions ####
#----------------------------------------------------------------------------------------------------#
# Model
fxBasicModel <- function (t, state, parameters) {
  
  with(as.list(c(state, parameters)),    
       
       {
         
         # Adjust risk/probability parameters to yearly risks (rough, but close enough)
         numTBMortP   <- ifelse(intInfDurYrs <= 1, numTBMortP,
                             numTBMortP * (1 / intInfDurYrs)) # TB mortality probability distributed  
                                                        # over period of illness
         numCureP   <- 1 - numTBMortP 
         numCureP   <- ifelse(intInfDurYrs <= 1, numCureP,
                             numCureP * (1 / intInfDurYrs)) # Recovery probability distributed over 
                                                        # period of illness
         
         # The duration of illness also needs to modify the beta value if it is less than one
         intBeta <- ifelse(intInfDurYrs >= 1, intBeta,
                           intBeta * intInfDurYrs)
         
         # Adjust the death rate so that it and the TB mortality equal the birth rate
         numMu <- numBirthR - (numTBMortP * I)  
         
         numTotalPop  <- S + SPE + L + I + R  # total population
         lambda <- intBeta * I/numTotalPop    # force of Infection - the risk that a susceptible
                                              # person is infected in a year
         
         # Change in 'susceptible - Mtb naive' compartment
         dS <- numBirthR -                     # new births (population that results from the deaths) 
           S * lambda * numFastP -             # minus those that become infected and progress quickly
           S * lambda * (1 - numFastP) -       # minus those that head to the latent/imm compartment
           S * numMu                           # minus those that die  
         
         # Change in 'susceptible - previously exposed' compartment
         dSPE <- L * numReverP +               # add those from the latent/imm compartment who revert
           R * numReverP -                     # add the those from the recovered compartment who revert
           SPE * lambda * numFastP -           # minus those that become infected and progress quickly
           SPE * lambda * (1 - numFastP) -     # minus those that head to the latent/imm compartment
           SPE * numMu -                       # minus those that die 
           SPE * numSusReactR                  # minus those that reactivate (assuming that TST/IGRA 
                                               # sill have a low reactivation risk)
         
         # Change in latent/imm compartment
         dL <- S * lambda * (1 - numFastP) +   # add 'susceptible - naive' multiplied by FOI and the fraction 
                                               # that don't progress to active disease
           SPE * lambda * (1 - numFastP) -     # add 'susceptible - prior exposure' multiplied by FOI 
                                               # and the fraction that don't progress to active disease
           L * numMu -                         # minus those that die 
           L * numReactR -                     # minus those that reactivate
           L * numReverP -                     # minus those that revert to negative
           L * (lambda * numFastP * (1 - numImm)) # minus those that get reinfected...despite being
                                               # somewhat protected
                
         # Change in active TB
         dI <-  S * lambda * numFastP +       # add 'susceptible - naive' multiplied by FOI and the fraction 
                                              # that progress to active disease
           SPE * lambda * numFastP +          # add 'susceptible - prior' exposure multiplied by FOI and 
                                              # the fraction that progress to active disease
           R * (lambda * numFastP * (1 - numImm)) + # add those in recovered compartment that get reinfected 
                                              # and progress to active disease
           L * (lambda * numFastP * (1 - numImm)) + # add those in the latent/imm compartment
                                              # that get reinfected and progress to active disease
           L * numReactR +                    # add those in the latent/imm compartment that reactivate 
           R * numRelaps -                    # add those in the recovered compartment that relapse 
           I * numMu -                        # minus those that die 
           I * numTBMortP -                      # minus those that die from TB
           I * numCureP +                      # minus those that recover/self-cure
           SPE * numSusReactR                 # add those in the 'Susceptible - prior exposure' 
                                              # compartment that reactivate
         
         # Change in recovered compartment
         dR <-  I * numCureP -                 # add those with active TB that recover
           R * (lambda * numFastP * (1 - numImm)) - #  minus those that have recovered that 
                                              # get reinfected and progress to active disease
           R * numRelaps -                    # minus those that relapse  
           R * numMu -                        # minus those that die 
           R * numReverP                      # minus those that revert
         
         # Model outcomes - Incidence of active disease as a result of recent infection
         dIrecent <- S * (lambda * numFastP) + # those who are exposure naive who progress quickly to disease
           SPE * (lambda * numFastP) +         # those in 'Susceptible - prior exposure' compartment
                                               # who progress quickly to disease
           L * (lambda * numFastP * (1 - numImm)) + # those in the latent/imm compartment that  
                                               # progress quickly to disease
           R * (lambda * numFastP * (1 - numImm)) # those in the recovered compartment who progress 
         
         # Model outcomes - Incidence of active disease as a result of remote infection
         dIremote <- L * numReactR +  # those in the latent/imm compartment that reactivate 
           R * numRelaps +            # those in the recovered compartment that relapse 
           SPE * numSusReactR         # those in the 'Susceptible - prior exposure' compartment
                                      # that reactivate   
         
         # Model outcomes - Total Incidence of active disease
         dIncidence <- dIrecent + dIremote

         # calculate some useful output measures
         dRelap   <- R * numRelaps 
         dTBmort  <- I * numTBMortP   
         dTBreact <- L * numReactR
         dTBreactSPE <- SPE * numSusReactR 
         
         # Define output
         dx <- c(dS, dSPE, dL, dI, dR, dIncidence, dIrecent, dIremote, dRelap, 
                 dTBmort, dTBreact, dTBreactSPE)
         
         list(dx)
         
       }
  )
}

# Run model
fxRunModel <- function (dfStartParams, params_new, params_old, 
                        lstModelRun_new, lstTimeInterv, 
                        fx_scale, fx_basic, strName, PriorData) {
  
  # Starting conditions
  lstStart <- c(S = dfStartParams$S, 
                SPE = dfStartParams$SPE, 
                L = dfStartParams$L,  
                I = dfStartParams$I,  
                R = dfStartParams$R,
                Incidence = dfStartParams$Incidence,
                Irecent = dfStartParams$Irecent,  
                Iremote = dfStartParams$Iremote,  
                dRelap = dfStartParams$dRelap,  
                dTBmort = dfStartParams$dTBmort,  
                dTBreact = dfStartParams$dTBreact,  
                dTBreactSPE = dfStartParams$dTBreactSPE) 
  
  # Select type of function
  if (length(lstTimeInterv) == 1) {
    
    fx <- fx_basic
    
  }  else {
    
    fx <- fx_scale
    
  }
  
  # Run the model
  dfOut <- as.data.frame(ode(y = lstStart, times = lstModelRun_new, 
                             func = fx, parms = params_new))  
  
  # Model output
  numTotalPop  <- dfOut$S + dfOut$SPE + dfOut$L + dfOut$I + dfOut$R  
  numRateInc   <- 1e5 * (diff(dfOut$Incidence) / 
                          numTotalPop[1 : length(numTotalPop) - 1])
  numRemotFrac <- diff(dfOut$Iremote) / diff(dfOut$Incidence)
  time         <- dfOut$time[1 : length(dfOut$time) - 1]
  dfData       <- data.frame(Years = time + (intStartYear - intSimYears), 
                            Incidence = numRateInc)
  dfData$Sim   <- strName
  
  # If it is a first run, nothing to append 
  if (length(PriorData) == 1) {
    
    dfData <- dfData
    
  } else { # Append previous runs
    
    dfData  <- rbind(PriorData, dfData)
    
  }
  
  lstOutput <- list("out" = dfOut, 
                    "data" = dfData,
                    "params" = params_new)
  
  return(lstOutput)
  
}

# Intervention scaling function
fxScaleUp <- function (t, state, parameters, lstTimeInterv, parameters_old, fx) {
  
  scale <- min((t - lstTimeInterv[1]) / (lstTimeInterv[2] - lstTimeInterv[1]), 1);
  
  if (scale < 0) {
    
    scale = 0
    
  }
  
  pars_scaled <- parameters;
  
  pars_scaled <- parameters_old + scale * (parameters - parameters_old)
  
  return(fx(t, state, pars_scaled))
  
}

# Function handles to pass as arguments.
fx_basic <- fxBasicModel
fx_scale <- function(t, state, parameters) fxScaleUp(t, state, parameters, 
                                                     lstTimeInterv, lstParams, fx_basic)
fx_scaleR <- function(t, state, parameters) fxScaleUp(t, state, parameters, 
                                                      lstTimeInterv, lstParamsR, fx_basic)
fx_scaleRTwo <- function(t, state, parameters) fxScaleUp(t, state, parameters, 
                                                      lstTimeInterv, lstParamsRTwo, fx_basic)
#----------------------------------------------------------------------------------------------------#
#                                     #### Initialise the model ####
#----------------------------------------------------------------------------------------------------#
strScenarioName   <- "Initialise"

lstSim  <- fxRunModel(dfStart, lstParams, NA, lstModelRun, NA, NA, 
                      fxBasicModel, strScenarioName, NA) 

lstSimR  <- fxRunModel(dfStart, lstParamsR, NA, lstModelRun, NA, NA, 
                       fxBasicModel, strScenarioName, NA) 

lstSimRTwo  <- fxRunModel(dfStart, lstParamsRTwo, NA, lstModelRun, NA, NA, 
                       fxBasicModel, strScenarioName, NA) 

#----------------------------------------------------------------------------------------------------#
#                        #### Simulation 0 - Baseline - no intervention ####
#----------------------------------------------------------------------------------------------------#
# Initial conditions 
dfStartParams   <- tail(lstSim$out, 1)
dfStartParamsR  <- tail(lstSimR$out, 1)
dfStartParamsRTwo <- tail(lstSimRTwo$out, 1)

lstParamsBase0   <- lstParams
lstParamsBaseR0  <- lstParamsR
lstParamsBaseRTwo0 <- lstParamsRTwo

lstModelRun_new <- seq(intSimYears, intSimYears + 25 , by = 1)
lstTimeInterv   <- c(lstModelRun_new[2], lstModelRun_new[2] + intScaleYears)

# Run model
lstSim0 <- fxRunModel(dfStartParams, lstParamsBase0, lstParamsBase0, 
                      lstModelRun_new, NA, fx_scale, fx_basic, 
                      strName0, NA)

lstSimR0 <- fxRunModel(dfStartParamsR, lstParamsBaseR0, lstParamsBaseR0, 
                       lstModelRun_new, NA,  fx_scaleR, fx_basic, 
                       strName0, NA)

lstSimRTwo0 <- fxRunModel(dfStartParamsRTwo, lstParamsBaseRTwo0, lstParamsBaseRTwo0, 
                          lstModelRun_new, NA, fx_scaleRTwo, fx_basic, 
                          strName0, NA)

#----------------------------------------------------------------------------------------------------#
#                           #### Simulation 1 - Passive treatment ####
#----------------------------------------------------------------------------------------------------#
# An Intervention simulating introduction of treatment
# Update parameter results data frame to append new results
lstParams1  <- lstParamsBase0
lstParamsR1  <- lstParamsBaseR0
lstParamsRTwo1  <- lstParamsBaseRTwo0
dfDataPrior1 <- lstSim0$data
dfDataPriorR1 <- lstSimR0$data
dfDataPriorRTwo1 <- lstSimRTwo0$data

lstParams1["numRelaps"] <- numTxPRelaps
lstParamsR1["numRelaps"] <- numTxPRelaps
lstParamsRTwo1["numRelaps"] <- numTxPRelaps

lstParams1["numTBMortP"] <- numTxPMuTB
lstParamsR1["numTBMortP"] <- numTxPMuTB
lstParamsRTwo1["numTBMortP"] <- numTxPMuTB

lstParams1["intInfDurYrs"] <- numTxPInfDur
lstParamsR1["intInfDurYrs"] <- numTxPInfDur
lstParamsRTwo1["intInfDurYrs"] <- numTxPInfDur

lstSim1 <- fxRunModel(dfStartParams, lstParams1, lstParamsBase0, 
                      lstModelRun_new, lstTimeInterv, fx_scale, fx_basic, 
                      strName1, dfDataPrior1)

lstSimR1 <- fxRunModel(dfStartParamsR, lstParamsR1, lstParamsBaseR0, 
                      lstModelRun_new, lstTimeInterv, fx_scaleR, fx_basic, 
                      strName1, dfDataPriorR1)

lstSimRTwo1 <- fxRunModel(dfStartParamsRTwo, lstParamsRTwo1, lstParamsBaseRTwo0, 
                       lstModelRun_new, lstTimeInterv, fx_scaleRTwo, fx_basic, 
                       strName1, dfDataPriorRTwo1)

#----------------------------------------------------------------------------------------------------#
#                         #### Simulation 2 - Passive treatment + ACF ####
#----------------------------------------------------------------------------------------------------#
# An Intervention simulating introduction of treatment

# Update parameter results data frame to append new results
lstParams2  <- lstParams1
lstParamsR2  <- lstParamsR1
lstParamsRTwo2  <- lstParamsRTwo1

dfDataPrior2 <- lstSim1$data
dfDataPriorR2 <- lstSimR1$data
dfDataPriorRTwo2 <- lstSimRTwo1$data

lstParams2["numRelaps"] <- numACFTxRelaps
lstParamsR2["numRelaps"] <- numACFTxRelaps
lstParamsRTwo2["numRelaps"] <- numACFTxRelaps

lstParams2["numTBMortP"] <- numACFTxMuTB
lstParamsR2["numTBMortP"] <- numACFTxMuTB
lstParamsRTwo2["numTBMortP"] <- numACFTxMuTB

lstParams2["intInfDurYrs"] <- numACFTxInfDur
lstParamsR2["intInfDurYrs"] <- numACFTxInfDur
lstParamsRTwo2["intInfDurYrs"] <- numACFTxInfDur

lstSim2 <- fxRunModel(dfStartParams, lstParams2, lstParams1, 
                      lstModelRun_new, lstTimeInterv, fx_scale, fx_basic, 
                      strName2, dfDataPrior2)

lstSimR2 <- fxRunModel(dfStartParamsR, lstParamsR2, lstParamsR1, 
                       lstModelRun_new, lstTimeInterv, fx_scaleR, fx_basic, 
                       strName2, dfDataPriorR2)

lstSimRTwo2 <- fxRunModel(dfStartParamsRTwo, lstParamsRTwo2, lstParamsRTwo1, 
                        lstModelRun_new, lstTimeInterv, fx_scaleRTwo, fx_basic, 
                        strName2, dfDataPriorRTwo2)

#----------------------------------------------------------------------------------------------------#
#                            #### Simulation 3 - Transmission stop ####
#----------------------------------------------------------------------------------------------------#
# A hypothetical ceasing of transmission (and nothing else, i.e. no treatment)
# Update parameter results data frame to append new results
lstParams3  <- lstParamsBase0
lstParamsR3  <- lstParamsBaseR0
lstParamsRTwo3  <- lstParamsBaseRTwo0

dfDataPrior3 <- lstSim0$data
dfDataPriorR3 <- lstSimR0$data
dfDataPriorRTwo3 <- lstSimRTwo0$data

# Adjust the beta value to zero to simulate no transmission
lstParams3["intBeta"] <- 0
lstParamsR3["intBeta"] <- 0
lstParamsRTwo3["intBeta"] <- 0

lstSim3 <- fxRunModel(dfStartParams, lstParams3, lstParamsBase0, 
                      lstModelRun_new, lstTimeInterv, fx_scale, fx_basic, 
                      strName3, dfDataPrior3)

lstSimR3 <- fxRunModel(dfStartParamsR, lstParamsR3, lstParamsBaseR0, 
                       lstModelRun_new, lstTimeInterv, fx_scaleR, fx_basic, 
                       strName3, dfDataPriorR3) 

lstSimRTwo3 <- fxRunModel(dfStartParamsRTwo, lstParamsRTwo3, lstParamsBaseRTwo0, 
                       lstModelRun_new, lstTimeInterv, fx_scaleRTwo, fx_basic, 
                       strName3, dfDataPriorRTwo3) 

#----------------------------------------------------------------------------------------------------#
#                 #### Simulation 4 - Latent Mtb infection treatment strategy ####
#----------------------------------------------------------------------------------------------------#
# An Intervention simulating LTBI screening and treatment

# Update parameter results data frame to append new results
lstParams4  <- lstParams1
lstParamsR4  <- lstParamsR1
lstParamsRTwo4  <- lstParamsRTwo1
dfDataPrior4 <- lstSim0$data
dfDataPriorR4 <- lstSimR0$data
dfDataPriorRTwo4  <- lstSimRTwo0$data

# Modify the reactivation and fast progression rates
lstParams4["numReactR"] <- numReactRPrev 
lstParamsR4["numReactR"] <- numReactRPrev 
lstParamsRTwo4["numReactR"] <- numReactRPrev

lstParams4["numFastP"] <- numFastPR 
lstParamsR4["numFastP"] <- numFastPR 
lstParamsRTwo4["numFastP"] <- numFastPR 

lstSim4 <- fxRunModel(dfStartParams, lstParams4, lstParamsBase0, 
                      lstModelRun_new, lstTimeInterv, fx_scale, fx_basic, 
                      strName4, dfDataPrior4)


lstSimR4 <- fxRunModel(dfStartParamsR, lstParamsR4, lstParamsBaseR0, 
                       lstModelRun_new, lstTimeInterv, fx_scaleR, fx_basic, 
                       strName4, dfDataPriorR4)


lstSimRTwo4 <- fxRunModel(dfStartParamsRTwo, lstParamsRTwo4, lstParamsBaseRTwo0, 
                          lstModelRun_new, lstTimeInterv, fx_scaleRTwo, fx_basic, 
                          strName4, dfDataPriorRTwo4)

#----------------------------------------------------------------------------------------------------#
#                              #### Simulation 5 - Zero reactivation ####
#----------------------------------------------------------------------------------------------------#
# An Intervention simulating zero TB reactivation (and nothing else, i.e. no treatment)
# Update parameter results data frame to append new results
lstParams5  <- lstParamsBase0
lstParamsR5  <- lstParamsBaseR0
lstParamsRTwo5  <- lstParamsBaseRTwo0
dfDataPrior5 <- lstSim0$data
dfDataPriorR5 <- lstSimR0$data
dfDataPriorRTwo5 <- lstSimRTwo0$data

# Modify the reactivation and fast progression rates
lstParams5["numReactR"] <- 0 
lstParamsR5["numReactR"] <- 0 
lstParamsRTwo5["numReactR"] <- 0 

lstParams5["numSusReactR"] <- 0 
lstParamsR5["numSusReactR"] <- 0 
lstParamsRTwo5["numSusReactR"] <- 0 

lstSim5 <- fxRunModel(dfStartParams, lstParams5, lstParamsBase0, 
                      lstModelRun_new, lstTimeInterv, fx_scale, fx_basic, 
                      strName5, dfDataPrior5)

lstSimR5 <- fxRunModel(dfStartParamsR, lstParamsR5, lstParamsBaseR0, 
                       lstModelRun_new, lstTimeInterv, fx_scaleR, fx_basic, 
                       strName5, dfDataPriorR5)

lstSimRTwo5 <- fxRunModel(dfStartParamsRTwo, lstParamsRTwo5, lstParamsBaseRTwo0, 
                          lstModelRun_new, lstTimeInterv, fx_scaleRTwo, fx_basic, 
                          strName5, dfDataPriorRTwo5)

# ================================================================================================= #
#                              #### Create plots ####
# ================================================================================================= #
# Plot the findings from the model with no reversion
plot <- unique(rbind(lstSim3$data, lstSim2$data, 
                     lstSim4$data, lstSim5$data))

# Creating the geom_text for the plot giving the proportion of disease due to fast progression
strRecentPcnt <- paste0(round(sum(lstSim1$out$Irecent) / 
                                (sum(lstSim1$out$Irecent) + 
                                   sum(lstSim1$out$Iremote)) * 100, 1), "% fast progression")

# Position this text properly in the plot
numXTextCoord <- 2034
setDT(lstSim3$data)
numYTextCoord <- lstSim3$data[Years == intStartYear & Sim %like% "Baseline", Incidence] + 10

plot <- 
  ggplot(data = plot, mapping = aes(x = Years, y = Incidence, 
                                    col = Sim, linetype = Sim)) + 
  geom_line(linewidth = numLineSize) +
  labs(y = "Rate per 100,000 population",
       col = "Simulation", linetype = "Simulation",
       tag = "a)") +  
  geom_text(x = numXTextCoord, y = numYTextCoord, 
            label = strRecentPcnt, colour = strTextColour,
            size = numGeomTextSize) +
  coord_cartesian(ylim = c(0,
                           numMaxYaxis
  )) +
  scale_colour_manual(values = lstColourCode) +
  scale_linetype_manual(values = lstLineStyleCode) +
  theme_bw() + 
  theme(text = element_text(size = numTextSize, colour = strTextColour),
        legend.position = "none",
        axis.title.y = element_text(margin = margin(t = 0, r = numDist, b = 0, l = 0),
                                    size = numTextSize, colour = strTextColour),
        axis.title.x = element_text(margin = margin(t = numDist, r = 0, b = 0, l = 0),
                                    size = numTextSize, colour = strTextColour),
        axis.text = element_text(size = numTextSize, colour = strTextColour),
        legend.text = element_text(size = numLegendTextSize, colour = strTextColour))

# Plot the findings from the model with reversion (first example, low levels)
plotR <- unique(rbind(lstSimR3$data, lstSimR2$data, 
                      lstSimR4$data, lstSimR5$data))

# Creating the geom_text for the plot giving proportion of disease due to fast progression
strRecentPcnt <- paste0(round(sum(lstSimR1$out$Irecent) / 
                                (sum(lstSimR1$out$Irecent) + 
                                   sum(lstSimR1$out$Iremote)) * 100, 1), "% fast progression")

setDT(lstSimR3$data)
numYTextCoord <- lstSimR3$data[Years == intStartYear & Sim %like% "Baseline", Incidence] + 10

plotR <- 
  ggplot(data = plotR, mapping = aes(x = Years, y = Incidence, 
                                     col = Sim, linetype = Sim)) + 
  geom_line(size = numLineSize) +
  labs(y = "Rate per 100,000 population",
       col = "Simulation", linetype = "Simulation",
       tag = "b)") +
  scale_colour_manual(values = lstColourCode) +
  geom_text(x = numXTextCoord, y = numYTextCoord, 
            label = strRecentPcnt, colour = strTextColour,
            size = numGeomTextSize) +
  scale_linetype_manual(values = lstLineStyleCode) +
  coord_cartesian(ylim = c(0,
                           numMaxYaxis
  )) +
  theme_bw() + 
  theme(text = element_text(size = numTextSize, colour = strTextColour),
        legend.position = "none",
        axis.title.y = element_text(margin = margin(t = 0, r = numDist, b = 0, l = 0),
                                    size = numTextSize, colour = strTextColour),
        axis.title.x = element_text(margin = margin(t = numDist, r = 0, b = 0, l = 0),
                                    size = numTextSize, colour = strTextColour),
        axis.text = element_text(size = numTextSize, colour = strTextColour),
        legend.text = element_text(size = numLegendTextSize, colour = strTextColour))

# Plot the findings from the model with reversion (second example, higher levels)
plotRTwo <- unique(rbind(lstSimRTwo3$data, lstSimRTwo2$data, 
                         lstSimRTwo4$data, lstSimRTwo5$data))

# Creating the geom_text for the plot giving proportion of disease due to fast progression
strRecentPcnt <- paste0(round(sum(lstSimRTwo1$out$Irecent) / 
                                (sum(lstSimRTwo1$out$Irecent) + 
                                   sum(lstSimRTwo1$out$Iremote)) * 100, 1), "% fast progression")

setDT(lstSimRTwo3$data)
numYTextCoord <- lstSimRTwo3$data[Years == intStartYear & Sim %like% "Baseline", Incidence] + 10

plotRTwo <- 
  ggplot(data = plotRTwo, mapping = aes(x = Years, y = Incidence, 
                                     col = Sim, linetype = Sim)) + 
  geom_line(linewidth = numLineSize) +
  labs(y = "Rate per 100,000 population",
       col = "Simulation", linetype = "Simulation",
       tag = "c)") +
  geom_text(x = numXTextCoord, y = numYTextCoord, 
            label = strRecentPcnt, colour = strTextColour,
            size = numGeomTextSize) +
  scale_colour_manual(values = lstColourCode) +
  scale_linetype_manual(values = lstLineStyleCode) +
  coord_cartesian(ylim = c(0,
                           numMaxYaxis
  )) +
  theme_bw() + 
  theme(text = element_text(size = numTextSize, colour = strTextColour),
        legend.position = "none",
        axis.title.y = element_text(margin = margin(t = 0, r = numDist, b = 0, l = 0),
                                    size = numTextSize, colour = strTextColour),
        axis.title.x = element_text(margin = margin(t = numDist, r = 0, b = 0, l = 0),
                                    size = numTextSize, colour = strTextColour),
        axis.text = element_text(size = numTextSize, colour = strTextColour),
        legend.text = element_text(size = numLegendTextSize, colour = strTextColour))

# Legend plot
plotLegend <- unique(rbind(lstSimR3$data, lstSimR2$data,
                           lstSimR4$data, lstSimR5$data))

plotLegend$Sim <- factor(plotLegend$Sim, levels = lstSimNames)


plotLegend <-
  ggplot(data = plotLegend, mapping = aes(x = Years, y = Incidence,
                                          col = Sim, linetype = Sim)) +
  geom_line(size = numLineSize) +
  labs(col = "Simulation", linetype = "Simulation") +
  scale_colour_manual(values = lstColourCode,
                      labels = lstLabels) +
  scale_linetype_manual(values = lstLineStyleCode,
                        labels = lstLabels) +
  theme_bw() +
  theme(text = element_text(size = numTextSize, colour = strTextColour),
        legend.position = lstLegendPosition,
        legend.key.width = unit(6, "line"),
        legend.text.align = 0,
        legend.text = element_text(size = numLegendTextSize, colour = strTextColour))
plotLegend

get_legend <- function(myggplot){
  
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  
  return(legend)
  
}

legend <- get_legend(plotLegend)

# View plot
grid.arrange(plot, plotR, plotRTwo, ncol = 3)

# ================================================================================================= #
#       #### Create some simple table outputs for examining model outputs/differences ####
# ================================================================================================= #
dfNone <- lstSim1$out
dfR <- lstSimR1$out
dfR2 <- lstSimRTwo1$out
setDT(dfNone)
setDT(dfR)
setDT(dfR2)
dfNone[, SumSLIR := S + SPE + L + I + R]
dfNone[, SumRecentRemote := Irecent + Iremote]
dfR[, SumSLIR := S + SPE + L + I + R]
dfR[, SumRecentRemote := Irecent + Iremote]
dfR[, SumSusc := S + SPE]
dfR2[, SumSLIR := S + SPE + L + I + R]
dfR2[, SumRecentRemote := Irecent + Iremote]
dfR2[, SumSusc := S + SPE]

# ================================================================================================= #
#                                      #### Save plots ####
# ================================================================================================= #
tiff(paste0(strOutputPath, 'Figure6.tiff'), units = "in",
     width = 14, height = 10, res = 300)

grid.arrange(plot, plotR, plotRTwo, legend,
             layout_matrix = rbind(c(1, 2, 3),
                                   c(6, 4, 6)))

dev.off()
