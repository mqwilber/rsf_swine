## Description
## -----------
##
## For all of the formatted studies, fit four different models to the movement data
## for each pig individually
##
## 1. A main effects model where there is no effect of day or season on effects
## 2. A daily model in which certain effects can vary cyclically with the hour 
##    of the day
## 3. A seasonal and daily model in which effects can vary both seasonally and 
##    over the a day
## 4. A temperature/precipitation model where the effect of season is replaced
##    by an interactive effect of temperature and precipitation.
##    Note that for this model, collaring needs to longer than 1 month for the
##    the temperature and precipitation variable to make sense 
##    (i.e. you need variation in these variables)
##
## TODOS:  
## 1. A few of the effects of crop locations are huge as there are not enough
## crop locations to estimate.  Need to limit the crop effect analysis to pigs that
## have spent only a certain amount of time in crop fields.  I would have thought 
## regularizing would have accounted for this...
##
## 2. Clean up this script to make it easier to understand
## 3. Save all of the movement results. 
## 4. 


library(data.table)
library(raster)
library(yaml)
library(glmnet)
library(ggplot2)
require(doMC)
source("pigfxns.R")

anal_params = yaml.load_file("analysis_parameters.yml")
allstudies = Sys.glob("~/Repos/rsf_swine/results/glmdata_by_study/*.csv")
studynames = sapply(strsplit(basename(allstudies), ".", fixed=T), function(x) x[1])

seasons = c("winter", "spring", "summer", "fall")

model1params = list()
allnongam = list()
model2params = list()
alldailyeffects = list()
alldailyloceffects = list()
model3params = list()
allseasoneffects = list()
allseasonbb = list()
alltpffects = list()

# Loop through analysis for each study
for(studynm in "florida"){

  cat("Beginning analysis of", studynm, "\n")
  dat = fread(paste0("~/Repos/rsf_swine/results/glmdata_by_study/", studynm, ".csv"))
  dat[ , datetime:=as.POSIXct(t, origin = '1970-01-01', tz = 'GMT')]
  dat[ , hourofday:=hour(datetime)]
  dat[ , monthofyear:=month(datetime)]

  # Extract the crop-using pigs
  allpigs = unique(dat$pigID)
  incroppigs = dat[z == 1][, list(incrop=any(crop_loc == 1)), by=pigID][incrop == TRUE, pigID]
  croppigs = dat[pigID %in% incroppigs]

  timeincrops = croppigs[(z == 1) & (crop_loc == 1), list(tottime = sum(tau)), by=pigID]

  # Pigs must spend longer than 1 hour in crops to be crop users
  # longcroppigs = timeincrops[tottime > 3600]$pigID
  longcroppigs = incroppigs

  ctypes = anal_params$croptypes
  croptypes = croppigs[z == 1, lapply(.SD, function(x) any(x == 1)),
                          .SDcols=paste0(ctypes, "_loc")]

  if(nrow(croptypes) != 0){

    cropsused = colnames(croptypes)[(croptypes[1, ] == TRUE) & 
                                        (colnames(croptypes) != "crop_loc")]

    cropdistgrad = sapply(strsplit(cropsused, "_loc"), function(x) paste0(x[1], "dists_grad"))
  } else{

    # If no crops are used just set as a dummy variable
    cropsused = "fruit_and_nuts_loc"
    cropdistgrad = "fruit_and_nutsdists_grad"
  }

  # Define columns for analysis below
  stdcols = c(c("ndvi_loc", "ndvi_grad", "masting_loc",
              "masting_grad", "water_grad",
              "elevation_loc", "elevation_grad",
              "canopycover_loc", "canopycover_grad", "crw", "developed_grad"), 
              cropdistgrad)
  nonstdcols = cropsused

  # Interaction between directional movement and presence in crop field
  interactions = paste0("crw:", cropsused) 
  splinecols = c(c("ndvi_grad", "masting_grad", "elevation_grad",
                 "canopycover_grad", "water_grad", "crw"), cropdistgrad)
  seasonalcols = c(c("masting_grad", "masting_loc", "ndvi_grad", 
                       "ndvi_loc", "water_grad", "crw"), cropdistgrad, cropsused,
                        interactions)

  # Pattern list for plotting
  plist = c(c("^base_hour_", "^ndvi_grad_", 
              "^masting_grad_", "^elevation_grad_", 
              "^canopycover_grad_", "^water_grad_*", "^crw_"),
              paste0("^", cropdistgrad, "_"), 
              paste0("^", cropsused, "_"))

  plistseas = c(c("^base_hour_", "^ndvi_grad_", 
              "^masting_grad_", "^water_grad_", "^crw_"),
              paste0("^", cropdistgrad, "_"),
              paste0("^", cropsused, "_"))

  # Flip the coefficients for interpretability
  flippattern = c("^water_grad_", paste0("^", cropdistgrad, "_"))
  flippatternbase = c("water_grad", cropdistgrad, "developed_grad")

  # Append at the end of each file
  saveappend = "crop_day_factor"

  ###################################################### 
  ### Step 1: Fit model with no time-varying effects ###
  ######################################################

  # Only fitting to pigs that use crops
  cat("Model 1: Main effects", "\n")

  allres = list()
  nongam = list()
  allbestbetas = list()

  for(pig in allpigs){
    
    # Extract a format data
    tdat = dat[pigID == pig]
    cat("Fitting pig", pig, "\n")
    
    # Set fngrad to 0 if it is NA, otherwise keep it the same
    tdat = remove_nas_from_crop(tdat, cropdistgrad)

    # If pigs are not in crops, don't fit the interactions
    if(pig %in% incroppigs){
      Xfull = build_design_matrix2(tdat, stdcols, cropsused, 
                                            interactions=interactions)
    } else{
      Xfull = build_design_matrix2(tdat, stdcols, cropsused, 
                                            interactions=NULL)
    }

    registerDoMC(cores=4)
    tfit = cv.glmnet(Xfull, y=tdat$z, offset=log(tdat$tau), 
                         family="poisson", 
                         alpha=1, nlambda=20, nfolds=5, parallel=TRUE)
    
    nongam[[pig]] = tfit$glmnet.fit$beta[, which(tfit$lambda == tfit$lambda.min), drop=T]

  }

  nongam_dt = do.call(rbind, lapply(names(nongam), function(x) {
                                                      t = melt(nongam[[x]])
                                                      t$coef = rownames(t)
                                                      t$pigID = x
                                                      return(t)}
                                                                  ))
  nongam_dt$study = studynm
  model1params[[studynm]] = nongam_dt

  nongam_dt$cropuser = nongam_dt$pigID %in% longcroppigs
  flipind = nongam_dt$coef %in% flippatternbase
  nongam_dt[flipind, "value"] = nongam_dt[flipind, "value"] * -1

  allnongam[[studynm]] = nongam_dt

  tplot = ggplot(nongam_dt) + geom_boxplot(aes(x=coef, y=value)) +
                      geom_point(aes(x=coef, y=value, color=pigID)) + 
                      geom_hline(aes(yintercept=0)) + theme_bw() + 
                      xlab("Coefficient Name") + ylab("Effect size") + 
                      theme(axis.text.x = element_text(angle = 90, hjust = 1))

  ggsave(file.path("../results/effect_plots",
                paste0("model1maineffects_", studynm, saveappend, ".pdf")), 
                tplot, width=10, height=7)


  ###################################################### 
  ### Step 2: Fit model with daily effects           ###
  ######################################################

  cat("Model 2: Daily effects", "\n")

  dailybestbetas = list()
  dailyeffects = list()

  for(pig in allpigs){
    
    # Extract a format data
    tdat = dat[pigID == pig]
    cat("Fitting pig", pig, "\n")

    # Set fngrad to 0 if it is NA, otherwise keep it the same
    tdat = remove_nas_from_crop(tdat, cropdistgrad)
    
    if(pig %in% incroppigs){
      Xgam = build_daily_design_matrix(tdat, stdcols, 
                              nonstdcols, splinecols, interactions=interactions)
    } else{

      Xgam = build_daily_design_matrix(tdat, stdcols, 
                              nonstdcols, 
                              splinecols, 
                              interactions=NULL)
    }
    
    registerDoMC(cores=4)
    fitgam = cv.glmnet(Xgam, y=tdat$z, offset=log(tdat$tau), family="poisson", 
                       alpha=1, nlambda=20, nfolds=5, parallel=TRUE)
    
    # Save the best fit betas
    bestbetas = fitgam$glmnet.fit$beta[, which(fitgam$lambda == fitgam$lambda.min),
                                                           drop=T]
    dailybestbetas[[pig]] = bestbetas
    
    dailyres = list()
    for(pattern in plist) {
      
      # Make sure model has this coef
      if(length(grep(pattern, names(bestbetas))) > 0){
        # Flip the sign for interpretability
        flip = ifelse(pattern %in% flippattern, -1, 1)

        eff = flip * bestbetas[grep(pattern, names(bestbetas))]
        ind = order(names(eff))

        tod = c("evening", "midday", "morning")
        effdt = data.frame(hour=tod, beta=eff[ind], coef=pattern)
        dailyres[[pattern]] = effdt
      }

    }
    
    dailyres = do.call(rbind, dailyres)
    dailyres$pigID = pig
    dailyeffects[[pig]] = dailyres
    
  }

  # Plot the daily results
  dailyeffects_dt = as.data.table(do.call(rbind, dailyeffects))



  dailyeffects_dt$hour = factor(dailyeffects_dt$hour, 
                        levels=c("morning", "midday", "evening"))

  dailyeffects_dt$study = studynm
  dailyeffects_dt$cropuser = dailyeffects_dt$pigID %in% longcroppigs
  alldailyeffects[[studynm]] = dailyeffects_dt

  meandailyeffects_dt = dailyeffects_dt[, list(meanbeta=mean(beta)), by=list(coef, hour)]
  tplot = ggplot(data=NULL) + geom_boxplot(data=dailyeffects_dt, 
                                  aes(x=hour, y=beta), alpha=0.3) +  
                             geom_point(data=dailyeffects_dt, 
                                  aes(x=hour, y=beta, color=pigID)) + 
                      geom_point(data=meandailyeffects_dt, aes(x=hour, y=meanbeta), 
                                  size=1, color="red") +
                      xlab("Time of day") + 
                      ylab("Effect on daily movement rate") + 
                      theme_bw() + facet_wrap(~coef, scales="free")

  ggsave(file.path("../results/effect_plots", 
                paste0("model2dailyeffects_", studynm, saveappend, ".pdf")), 
                tplot, width=10, height=7)


  dailybestbetas_dt = do.call(rbind, lapply(names(dailybestbetas), function(x) {
                                                      t = melt(dailybestbetas[[x]])
                                                      t$coef = rownames(t)
                                                      t$pigID = x
                                                      return(t)}
                                                                  ))
  dailybestbetas_dt = as.data.table(dailybestbetas_dt)
  dailybestbetas_dt$study = studynm
  model2params[[studynm]] = dailybestbetas_dt

  dailyloceffects = rbind(dailybestbetas_dt[coef %like% "*_loc"], 
                          dailybestbetas_dt[coef %like% "developed"])
  dailyloceffects$study = studynm
  dailyloceffects$cropuser = dailyloceffects$pigID %in% longcroppigs

  # Flip coefficients for interpretability
  flipind = dailyloceffects$coef %like% paste0(flippatternbase, collapse="|")
  dailyloceffects[flipind, value:=value * -1]

  alldailyloceffects[[studynm]] = dailyloceffects

  tplot = ggplot(dailyloceffects) + geom_boxplot(aes(x=coef, y=value)) + 
                      geom_point(aes(x=coef, y=value, color=pigID)) + 
                      geom_hline(aes(yintercept=0)) + theme_bw() + 
                      xlab("Coefficient Name") + ylab("Effect size") + 
                      theme(axis.text.x = element_text(angle = 90, hjust = 1))

  ggsave(file.path("../results/effect_plots", 
                paste0("model2maineffects_", studynm, saveappend, ".pdf")), 
                tplot, width=10, height=7)

  ###################################################### 
  ### Step 3: Fit model with seasonal effects        ###
  ######################################################

  cat("Model 3: Seasonal and daily effects", "\n")

  seasonalbestbetas = list()
  seasonaleffects = list()

  for(pig in allpigs){
    
    # Extract a format data
    tdat = dat[pigID == pig]
    cat("Fitting pig", pig, "\n")

    # Set fngrad to 0 if it is NA, otherwise keep it the same
    tdat = remove_nas_from_crop(tdat, cropdistgrad)
   
    if(pig %in% incroppigs){ 
      Xseason = build_seasonal_design_matrix(tdat, stdcols, nonstdcols, 
                            splinecols, seasonalcols, interactions=interactions)
    } else{
      Xseason = build_seasonal_design_matrix(tdat, stdcols, nonstdcols, 
                            splinecols, 
                            seasonalcols[-which(seasonalcols %in% interactions)], 
                            interactions=NULL)
    }
    
    registerDoMC(cores=4)
    fitseason = cv.glmnet(Xseason, y=tdat$z, offset=log(tdat$tau), 
                          family="poisson", alpha=1, nlambda=20, nfolds=5, 
                          parallel=TRUE)
    
    # Save the best fit betas
    bestbetas = fitseason$glmnet.fit$beta[, 
                        which(fitseason$lambda == fitseason$lambda.min), drop=T]
    seasonalbestbetas[[pig]] = bestbetas

    dailyres = list()
    rsums = apply(Xseason[, seasons], 2, sum)
    inseason = names(rsums)[rsums != 0]
    
    for(pattern in plistseas) {
      
      if(length(grep(pattern, names(bestbetas))) > 0){

        dailyres[[pattern]] = list()
        seasbetas = bestbetas[grep(pattern, names(bestbetas))]
        
        for(season in inseason){


          flip = ifelse(pattern %in% flippattern, -1, 1)
          eff = flip*seasbetas[grep(season, names(seasbetas))]

          ind = order(names(eff))
          tod = c("evening", "midday", "morning")

          effdt = data.frame(hour=tod, beta=eff[ind], coef=pattern)
          effdt$season = season
          dailyres[[pattern]][[season]] = effdt
        }

      }
    }
    
    seasonres = do.call(rbind, lapply(dailyres, function(x) do.call(rbind, x)))
    seasonres$pigID = pig
    seasonaleffects[[pig]] = seasonres

  }


  # Plot seasonal + daily effects
  seasoneffects_dt = as.data.table(do.call(rbind, seasonaleffects))
  seasoneffects_dt$season = factor(seasoneffects_dt$season, 
                                    levels=c("winter", "spring", "summer", "fall"))
  seasoneffects_dt$hour = factor(seasoneffects_dt$hour, 
                                        levels=c("morning", "midday", "evening"))

  seasoneffects_dt$study = studynm
  seasoneffects_dt$cropuser = seasoneffects_dt$pigID %in% longcroppigs
  allseasoneffects[[studynm]] = seasoneffects_dt

  seasonmeanpigs = seasoneffects_dt[, list(meanbeta=mean(beta)), 
                                              by=list(coef, season, hour)]
  seasonmeanpigs$season = factor(seasonmeanpigs$season, 
                              levels=c("winter", "spring", "summer", "fall"))
  seasonmeanpigs$hour = factor(seasonmeanpigs$hour, 
                                      levels=c("morning", "midday", "evening"))


  numseason = length(unique(seasonmeanpigs$season))
  tplot = ggplot(data=NULL) + geom_boxplot(data=seasoneffects_dt, 
                                  aes(x=hour, y=beta), alpha=0.3) +
                              geom_point(data=seasoneffects_dt, 
                                  aes(x=hour, y=beta, color=pigID)) +
                              geom_point(data=seasonmeanpigs, aes(x=hour, y=meanbeta), 
                                        size=2, color="red") + 
                              xlab("Time of day") + 
                              ylab("Effect on pig movement") + 
                              theme_bw() + 
                              facet_wrap(~coef + season, scales="free", 
                                            nrow=length(plistseas), 
                                            ncol=numseason)

  ggsave(file.path("../results/effect_plots", 
                paste0("model3seasondailyeffects_", studynm, 
                                saveappend, ".pdf")), 
                tplot, width=5*numseason, height=15)

  # Plot just the seasonal effects
  sbb_dt = do.call(rbind, lapply(names(seasonalbestbetas), function(x) {
                                                      t = melt(seasonalbestbetas[[x]])
                                                      t$coef = rownames(t)
                                                      t$pigID = x
                                                      return(t)}
                                                                  ))
  sbb_dt = as.data.table(sbb_dt)
  sbb_dt$study = studynm
  model3params[[studynm]] = sbb_dt


  # Extract the location columns
  seasonalbb_dt = sbb_dt[((coef %like% "*_loc*") | 
                          (coef %like% "developed")) &
                         !(coef %like% "crw")]

  seasonalbb_dt$season = sapply(strsplit(as.character(seasonalbb_dt$coef), ":"), function(x) x[2])
  seasonalbb_dt$coef1 = sapply(strsplit(as.character(seasonalbb_dt$coef), ":"), function(x) x[1])

  seasonalbb_dt$season[is.na(seasonalbb_dt$season)] = "nonseasonal"
  seasonalbb_dt$season = factor(seasonalbb_dt$season, 
                    levels=c("winter", "spring", "summer", "fall", "nonseasonal"))

  seasonalbb_dt$study = studynm
  seasonalbb_dt$cropuser = seasonalbb_dt$pigID %in% longcroppigs
  seasonalbb_dt = as.data.table(seasonalbb_dt)


  flipind = seasonalbb_dt$coef %like% paste0(flippatternbase, collapse="|")
  seasonalbb_dt[flipind, value:=value * -1]

  allseasonbb[[studynm]] = seasonalbb_dt

  tplot = ggplot(seasonalbb_dt) + geom_boxplot(aes(x=season, y=value)) + 
                          geom_point(aes(x=season, y=value, color=pigID)) + 
                          facet_wrap(~coef1, scales="free") + theme_bw() + 
                          theme(axis.text.x = element_text(angle = 90, hjust = 1))

  ggsave(file.path("../results/effect_plots", 
                paste0("model3seasoneffects_", studynm, saveappend, ".pdf")), 
                tplot, width=10, height=7)

  ###################################################### 
  ### Step 4: Fit model with temp/precip effects     ###
  ######################################################

#   cat("Model 4: Temp/precip and daily effects", "\n")

#   tpbestbetas = list()
#   tpeffects = list()

#   for(pig in longcroppigs) {
    
#     # Extract a format data
#     tdat = dat[pigID == pig]
#     cat("Fitting pig", pig, "\n")

#     # Set fngrad to 0 if it is NA, otherwise keep it the same
#     tdat = remove_nas_from_crop(tdat, cropdistgrad)
    
#     Xtp = build_tempprecip_design_matrix(tdat, stdcols, nonstdcols, 
#                                                  splinecols, seasonalcols,
#                                                  df_hour=df_hour)

#     # Check on correlation between temp a

#     # Check if there is variability in temperature before fitting
#     if(!((sd(Xtp[, "temperature_loc"]) == 0) | (sd(Xtp[, "precipitation_loc"]) == 0))) {

    
#       registerDoMC(cores=4)
#       fittp = cv.glmnet(Xtp, y=tdat$z, offset=log(tdat$tau), 
#                             family="poisson", alpha=1, nlambda=20, nfolds=5, 
#                             parallel=TRUE)

#       bestbetas = fittp$glmnet.fit$beta[, which(fittp$lambda == fittp$lambda.min), drop=T]
#       tpbestbetas[[pig]] = bestbetas
      
#       # Build a monthly temperature and precipitation predictor
#       monthtemp = dat[, list(stemp=scale2(temperature_loc), monthofyear)][, 
#                            list(temperature_loc=mean(stemp)), by=monthofyear]
#       monthprecip = dat[, list(sprecip=scale2(precipitation_loc), monthofyear)][,
#                            list(precipitation_loc=mean(sprecip)), by=monthofyear]
#       monthtp = cbind(monthtemp, monthprecip[, list(precipitation_loc)])
#       monthtp[, ("temperature_loc:precipitation_loc"):=temperature_loc*precipitation_loc]
      
#       splinehour = s(hourofday, bs="cc", k=df_hour)
#       pred_hours = Predict.matrix(smooth.construct2(splinehour, tdat, NULL), 
#                       data.frame(hourofday=sort(unique(tdat$hourofday))))
      
#       dailyres = list()
#       inmonth = unique(tdat$monthofyear)
      
#       tpcols = c("temperature_loc", "precipitation_loc", 
#                                 "temperature_loc:precipitation_loc")
      
#       # Make the basis effects for each time-varying parameter
#       for(pattern in plistseas) {
        
#         dailyres[[pattern]] = list()
#         tpbetas = bestbetas[names(bestbetas) %like% pattern]
        
#         for(month in inmonth){
          
#           basenames = paste0(strsplit(pattern, "*", fixed=T)[[1]], 1:(df_hour - 1))
#           baseeff = pred_hours %*% t(t(tpbetas[basenames]))
          
#           for(tp in tpcols){
            
#             tnames = paste0(strsplit(pattern, "*", fixed=T)[[1]], 1:(df_hour - 1), ":", tp)
#             tpval = as.numeric(monthtp[monthofyear == month, tp, with=F])
#             teff = (pred_hours*tpval) %*% t(t(tpbetas[tnames]))
            
#             baseeff = baseeff + teff
#           }
          
#           flip = ifelse(pattern %in% c("water_grad_*", paste0(cropdistgrad, "_*")), -1, 1)
#           effdt = data.frame(hour=0:23, beta=flip*baseeff, coef=pattern)
#           effdt$monthofyear = month
#           dailyres[[pattern]][[month]] = effdt
#         }
#       }
      
#       monthres = do.call(rbind, lapply(dailyres, function(x) do.call(rbind, x)))
#       monthres$pigID = pig
#       tpeffects[[pig]] = monthres

#     } 

#   } # End pig for loop

#   # Plot the temperature-dependent effects, if there are any pigs that qualify
#   if(length(tpeffects) > 0) {

#     tpeffects_dt = as.data.table(do.call(rbind, tpeffects))
#     tpmeanpig = tpeffects_dt[, list(beta=mean(beta)), by=list(monthofyear, hour, coef)]
#     tpmeanpig$pigID = paste0("meanpig", studynm)
#     tpmeanpig = tpmeanpig[, colnames(tpeffects_dt), with=F]
#     tpeffects_full = rbind(tpeffects_dt, tpmeanpig)
#     tpeffects_full$study = studynm
#     alltpffects[[studynm]] = tpeffects_full

#     numpigs = length(unique(tpeffects_full$pigID))
#     tplot = ggplot(tpeffects_full) + geom_line(aes(x=hour, y=beta, 
#                                       color=as.factor(monthofyear))) + 
#                             xlab("Hour of day") + ylab("Effect on pig movement") + 
#                             theme_bw() + facet_wrap(~pigID + coef, scales="free", 
#                                               ncol=length(plistseas), nrow=numpigs)

#     ggsave(file.path("../results/effect_plots", 
#                 paste0("model4tempprecipeffects_", studynm, saveappend, ".pdf")), 
#                 tplot, width=20, height=20)
#   }

} # End study loop


####################################################################
### Save all of the estimated coefficients for post-hoc analysis ###
####################################################################

append = TRUE # Append to previous results

allresults = list(model1maineffects = allnongam, 
                  model1allparams = model1params,
                  model2dailyeffects = alldailyeffects,
                  model2maineffects = alldailyloceffects,
                  model2allparams = model2params,
                  model3seasonaleffects = allseasoneffects,
                  model3maineffects = allseasonbb,
                  model3allparams = model3params)

for(res in names(allresults)){

  fulldt = as.data.table(do.call(rbind, allresults[[res]]))

  if(append){
    tempdt = fread(file.path("../results/effect_size_data",
                 paste0(res, ".csv")))
    fulldt = rbind(tempdt, fulldt)
  }

  fwrite(fulldt, file.path("../results/effect_size_data",
                 paste0(res, ".csv")), row.names = FALSE)
}


