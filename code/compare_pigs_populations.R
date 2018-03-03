## Computes and saves the data for the CTMC analysis across all pigs and 
## all studies.  This data can then be used to perform the an analysis on
## how resources affect pig movement.
##
## This analysis is broken into 6 steps
## 
## 1. This first step performs study specific cleaning on each data file. This
## is tough to generalize because each study has its own quirks.  clean_studies.R
## contains a sweet of functions that performs this data cleaning for each study.
##
## 2. Once a study is cleaned, step 2 extracts the runs that that meet the
## diff time and length criteria.
##
## 3. The third step then fits a continuous time movement model to each run for
## each pig.  Each run is considered independent for a given pig.
##
## 4. The fourth step formats the necessary rasters for each pig and computes the
## CTMC conversion for each pig, each run, and each imputation. The result is a
## a series of data.tables that, after conversion, can be used to analyze pig movement.
##
## 5. Taking the data.tables from the last step, this step combines the 
## time-varying covariate columns into a single column. Then it selects the 
## necessary formatted columns and concatenates the data.tables across pigs.
##
## 6. The final step incorporates the population-level covariates that primarily
## differ between studies. This step also incorporates time-varying, 
## county-level variables such as drought.
##
## 7. The resulting data.frame is saved. 
##
## Author: Mark Wilber

list.of.packages <- c("data.table", "raster", "ggplot2", "lubridate", "sp",
                        "rgdal", "rgeos", "ctmcmove", "splines", "fda",
                        "yaml", "parallel")
new.packages <- list.of.packages[!(list.of.packages %in% 
																							 installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only=TRUE)

source('pigfxns.R')
source('clean_studies.R')
anal_params = yaml.load_file("analysis_parameters.yml")
sink("log_ctmc.txt") # log file

###########################################################
## Step 1: Clean individual studies for outliers
###########################################################

for(studynm in c("txcamp")){

  #studynm = "tejon"
  cat("Beginning analysis for", studynm, "\n")

  trimdat = eval(parse(text=paste("clean_", studynm, "(plotit=FALSE)", sep=""))) 
  trimdat$datetime = as.POSIXct(trimdat$datetime, tz="GMT")

  ###########################################################
  ## Step 2: Identify the runs in the study that are of sufficient length
  ###########################################################

  maxtime = anal_params$maxtime  # Length between collaring times can't be more than this
  minlength = anal_params$minlength # Run of collaring times must be at least this long

  unqpigs = unique(trimdat$pigID)

  newdat = list()

  # For each pig in the dataset, identify and save the runs with the above 
  # criteria
  for(pignm in unqpigs){
    
    tdat = trimdat[pigID == pignm]
    tdat$run_number = 0 
    runinds = split_runs(tdat$datetime, ctime=maxtime, clength=minlength)
    
    for(i in 1:length(runinds)){
      inds = runinds[[i]]
      tdat$run_number[inds] = i
    }
    
    keepinds = do.call(c, runinds)
    rundat = tdat[keepinds]
    newdat[[pignm]] = rundat
    
  }

  # This is data split up by consecutive runs.
  newdat = do.call(rbind, newdat)

  ###########################################################
  ## Step 3: Fit the continuous-time movement model for each pig and 
  ## each run
  ###########################################################

  unqpigs = unique(newdat$pigID)

  # Loop through and fit a continuous time movement model for each pig on each run
  parallel_pigpaths = function(pignm, newdat, anal_params){

  	cat("Computing continuous path for", pignm, "\n")
    tdat = newdat[pigID == pignm]
  	predpaths = lapply(unique(tdat$run_number), 
  											function(x){ 
  												continuous_path(tdat[run_number == x, 
  														list(x=longitude, y=latitude, datetime=datetime)], 
  														anal_params$timestep, anal_params$method)
  											})
  	return(predpaths)

  }


  allpaths = mclapply(unqpigs, parallel_pigpaths, newdat, anal_params, 
  																								mc.cores=anal_params$cores)
  names(allpaths) = unqpigs

  ###########################################################
  ## Step 4: For each pig, get the raster data and use CTMC to get data.frame
  ###########################################################

  parallel_pigctmc = function(pignm, anal_params, newdat, allpaths, studynm){

  	cat("Beginning CTMC for", pignm, "\n")
  	locvars = anal_params$locvars
  	gradvars = anal_params$gradvars
  	gradxyvar = anal_params$gradxyvar
  	landcovertypes = anal_params$landcovertypes
  	distgrad_types = anal_params$distgrad_types
  	croptypes = anal_params$croptypes
  	fixedtimevars = anal_params$fixedtimevars
  	monthlyvars = anal_params$monthlyvars
  	yearlyvars = anal_params$yearlyvars
  	default_cols = anal_params$default_cols

  	buffer = anal_params$buffer

  	tdat = newdat[pigID == pignm]
    extobj = extent(c(xmin=min(tdat$longitude) - buffer, 
    									xmax=max(tdat$longitude) + buffer, 
    									ymin=min(tdat$latitude) - buffer, 
    									ymax=max(tdat$latitude) + buffer))

    # Get the base raster on which to project all other rasters
    basefiles = Sys.glob(file.path(anal_params$covariate_dir, "croplayer", studynm,
                paste(studynm, anal_params$baserastype, "*.grd", sep="")))
    baseras = crop(raster(basefiles[1]), extobj)
    # plot(baseras)

    loc_stackproj = stack(process_covariates(locvars, studynm, baseras, 
                                   min(tdat$datetime), max(tdat$datetime), 
                                   ext="loc", landcovertypes = landcovertypes,
                                   croptypes=croptypes, 
                                   projectionMethod="bilinear"))
    
    grad_stackproj = stack(process_covariates(gradvars, studynm, baseras,
                                   min(tdat$datetime), max(tdat$datetime), 
                                   ext="grad", landcovertypes = landcovertypes[landcovertypes != "wetland"],
                                   projectionMethod="bilinear"))
    
    gradxy_stackproj = process_covariates(gradxyvar, studynm, baseras,
                                   min(tdat$datetime), max(tdat$datetime),
                                   ext="grad", croptypes=croptypes, 
                                   landcovertypes = landcovertypes, 
                                   distgrad = TRUE, distgrad_types=distgrad_types,
                                   projectionMethod="bilinear")
    
    gradx_stackproj = stack(lapply(gradxy_stackproj, function(m) m$xgrad))
    grady_stackproj = stack(lapply(gradxy_stackproj, function(m) m$ygrad))
    
    # # Ensure projections are the same...lose some edge values this way
    # baseras = loc_stack[[paste('cereals_', year(min(tdat$datetime)), "_loc", sep="")]]
    # loc_stackproj = stack(lapply(loc_stack, 
    # 															function(x) 
    # 																projectRaster(x, baseras, method="ngb")))
    # grad_stackproj = stack(lapply(grad_stack, 
    # 															function(x) 
    # 																projectRaster(x, baseras, method="bilinear")))
    # gradx_stackproj = stack(lapply(gradx_stack, 
    # 															function(x) 
    # 																projectRaster(x, baseras, method="ngb")))
    # grady_stackproj = stack(lapply(grady_stack, 
    # 															function(x) 
    # 																projectRaster(x, baseras, method="ngb")))
    
    # Compute the CTMC for each run and then combine them
    pigpaths = allpaths[[pignm]]
    runglmdat = list()

    for(j in 1:length(pigpaths)){
      
      cat("Run", j, "for", pignm, "\n")
      smallpath = pigpaths[[j]]
      
      # Get reduced extent
      minll = apply(tdat[run_number == j, list(x=longitude, y=latitude)], 2, min)
      maxll = apply(tdat[run_number == j, list(x=longitude, y=latitude)], 2, max)
      smallextobj = extent(c(xmin=minll[1] - buffer, xmax=maxll[1] + buffer, 
      											 ymin=minll[2] - buffer, ymax=maxll[2] + buffer))
      
      tglmdat = fit_ctmcmodel(tdat, pignm, crop(loc_stackproj, smallextobj), 
      												crop(grad_stackproj, smallextobj), xygrad=TRUE, 
                              xgrad.stack=crop(gradx_stackproj, smallextobj), 
                              ygrad.stack = crop(grady_stackproj, smallextobj),
                              predPath=smallpath, 
                              path2ctmcMethod=anal_params$path2ctmcMethod)
      runglmdat[[j]] = tglmdat
       
    }

    pigglmdat = do.call(rbind, runglmdat)
  	cat("Finished CTMC for", pignm, "\n")
    return(as.data.table(pigglmdat))

  }

  glmlist = mclapply(unqpigs, parallel_pigctmc, anal_params, newdat, allpaths, 
                        studynm, mc.cores=anal_params$cores)

  ###########################################################
  ## Step 5: Format the data.frame consistently across pigs and save result
  ###########################################################

  # Format columns and combine all data.tables
  locvars = anal_params$locvars
  gradvars = anal_params$gradvars
  gradxyvar = anal_params$gradxyvar
  landcovertypes = anal_params$landcovertypes
  distgrad_types = anal_params$distgrad_types
  croptypes = anal_params$croptypes
  fixedtimevars = anal_params$fixedtimevars
  monthlyvars = anal_params$monthlyvars
  yearlyvars = anal_params$yearlyvars
  default_cols = anal_params$default_cols

  # A bit ugly, but dropping wetland and canopycover so we don't get duplicates. Will probably want to make this cleaner
  regexcols_locs = unlist(get_regex_columns(locvars, croptypes, landcovertypes, "loc"))
  regexcols_grad = unlist(get_regex_columns(gradvars, c(), landcovertypes[landcovertypes != "wetland"], "grad"))
  regexcols_gradxy = unlist(get_regex_columns(gradxyvar, distgrad_types, landcovertypes[landcovertypes != "canopycover"], "grad"))
  allregexcols = c(regexcols_locs, regexcols_grad, regexcols_gradxy)
  cleancols = sapply(strsplit(allregexcols, ".*", fixed=T), function(x) paste(x[1], x[2], sep="_"))
  timevar_vect = get_timevar_status(allregexcols, croptypes, monthlyvars)

  tallglm = process_glmdata(glmlist, allregexcols, cleancols, timevar_vect)
  tallglm_form = lapply(tallglm, function(x) x[, c(default_cols, cleancols), with=FALSE])
  fullglm = as.data.table(do.call(rbind, tallglm_form))
  fullglm$study = studynm


  ###########################################################
  ## Step 6: Process the population-level covariates
  ###########################################################
  cat("Adding post-processing covariates for", studynm, "\n")
  fullglm = add_covariates(fullglm, anal_params$post_covariates, studynm,
                                    basedir=anal_params$covariate_dir)

  ###########################################################
  ## Step 7: Write the file to disk
  ###########################################################

  fwrite(fullglm, file=paste("../results/glmdata_by_study/", studynm, ".csv", sep=""))
  cat("Completed analysis for", studynm, "\n")

}

# End logging to file
sink()




