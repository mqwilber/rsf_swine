## Script reads in base croplayer for all studies and then computes
## the unique croptypes for each study as defined by NASS

library(raster)
library(data.table)
library(yaml)
source("pigfxns.R")
source("clean_studies.R")
sink("temp.log")

studysum = fread("../data/formatted/study_summary.csv")
croptypes = fread("../data/covariate_data/croplayer/cropgroupings.csv")
croptypes = croptypes[group_name != "nothing"]
anal_params = yaml.load_file("analysis_parameters.yml")

allcrops = list()

for(studynm in studysum$study){

	cat("Working on", studynm, "\n")

	files = Sys.glob(file.path("../data/covariate_data/croplayer/", 
												 studynm, "raw", "*.tif"))

	if(length(files) ==  0) 
		next

	# Load dummy crop rasters
	ras = raster(files[1])

	trimdat = clean_pigs(studynm, plotit=FALSE) 
  trimdat$datetime = as.POSIXct(strptime(trimdat$datetime, 
                                      format="%Y-%m-%d %H:%M:%S", tz="GMT"))

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

    # Only perform split if pig meets criteria
    if(runs(tdat$datetime, ctime=anal_params$maxtime, 
                           clength=anal_params$minlength)){

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
    
  }

  # This is data split up by consecutive runs.
  newdat = do.call(rbind, newdat)

  ###########################################################
  ## Step 3: Fit the continuous-time movement model for each pig and 
  ## each run
  ###########################################################

  unqpigs = unique(newdat$pigID)

  if(length(unqpigs) == 0) next # Skip any study without goodpigs

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

  parallel_pigctmc = function(pignm, anal_params, newdat, allpaths, studynm) {

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
    baseras = crop(ras, extobj)

		loc.stack = stack(baseras, baseras)
		names(loc.stack) = c("crop_loc", "crop_loc2")
		grad.stack = stack(baseras, baseras)
		names(grad.stack) = c("crop_grad", "crop_grad2")


    loc_stackproj = loc.stack
    grad_stackproj = grad.stack
    
    
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
                              crop(grad_stackproj, smallextobj),
                              timestep="30 mins",
                              xygrad=FALSE, 
                              predPath=smallpath, 
                              path2ctmcMethod=anal_params$path2ctmcMethod,
                              directions=anal_params$directions)
      runglmdat[[j]] = tglmdat
       
    }

    pigglmdat = do.call(rbind, runglmdat)
  	cat("Finished CTMC for", pignm, "\n")
    return(as.data.table(pigglmdat))

  }

  glmlist = mclapply(unqpigs, parallel_pigctmc, anal_params, newdat, allpaths, 
                        studynm, mc.cores=anal_params$cores)

  fullglm = do.call(rbind, glmlist)
  incell = fullglm[z == 1]

  # Time spent in crop types by pig
  timeincell = incell[, list(hours=sum(tau) / (60*60)), by=list(pigID, crop_loc)]
  cropvals = croptypes[group_name != "nothing", value]
  timeincrops = timeincell[crop_loc %in% cropvals]
  colnames(timeincrops)[2] = "value"
  timeincrops = merge(timeincrops, croptypes, by="value", all.x=TRUE)
  timeincrops$study = studynm

	allcrops[[studynm]] = timeincrops

}

# Save the crops that are being used for each study
allcrops_dt = do.call(rbind, allcrops)
fwrite(allcrops_dt, "../results/usedcrops_time.csv")
sink()