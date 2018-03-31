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
preprocessed = TRUE


for(studynm in c("tejon", "txcamp", "srel_contact", "tx_tyler_w2", 
                 "tx_tyler_w1", "fl_raoul", "florida", "mo_kurt0",
                 "cali0", "cali2")) {

  #studynm = "tejon"
  cat("Beginning analysis for", studynm, "\n")

  if(!preprocessed) {
    ###########################################################
    ## Step 1: Clean individual studies for outliers
    ###########################################################

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

      # Only perform split if pig meets criteria
      if(runs(tdat$datetime, ctime=anal_params$maxtime, clength=anal_params$minlength)){

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
      basefiles = Sys.glob(file.path(anal_params$covariate_dir, "croplayer", studynm,
                  paste(studynm, anal_params$baserastype, "*.tif", sep="")))
      baseras = crop(raster(basefiles[1]), extobj)

      # Determine the spatial scale of the base raster
      baseras = aggregate(baseras, fact=anal_params$fact)

      loc_stackproj = stack(process_covariates(locvars, studynm, baseras, 
                                     min(tdat$datetime), max(tdat$datetime), 
                                     ext="loc", landcovertypes = landcovertypes,
                                     croptypes=croptypes, 
                                     projectionMethod="bilinear"))
      
      grad_stackproj = stack(process_covariates(gradvars, studynm, baseras,
                                     min(tdat$datetime), max(tdat$datetime), 
                                     ext="grad", landcovertypes = landcovertypes,
                                     nngrad=TRUE, croptypes=croptypes,
                                     projectionMethod="bilinear"))
      
      
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
                                crop(grad_stackproj, smallextobj), xygrad=FALSE, 
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

    ###########################################################
    ## Step 5: Format the data.frame consistently across pigs and save result
    ###########################################################

    # Format columns and combine all data.tables
    locvars = anal_params$locvars
    gradvars = anal_params$gradvars
    # gradxyvar = anal_params$gradxyvar
    landcovertypes = anal_params$landcovertypes
    distgrad_types = anal_params$distgrad_types
    croptypes = anal_params$croptypes
    fixedtimevars = anal_params$fixedtimevars
    monthlyvars = anal_params$monthlyvars
    yearlyvars = anal_params$yearlyvars
    default_cols = anal_params$default_cols

    regexcols_locs = unlist(get_regex_columns(locvars, croptypes, landcovertypes, "loc"))
    regexcols_grad = unlist(get_regex_columns(gradvars, croptypes, landcovertypes, "grad"))
    #regexcols_gradxy = unlist(get_regex_columns(gradxyvar, distgrad_types, landcovertypes[landcovertypes != "canopycover"], "grad"))
    allregexcols = c(regexcols_locs, regexcols_grad)#, regexcols_gradxy)
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
  }
  else {

    cat("Loading preprocessed data for", studynm, "\n")
    fullglm = fread(paste("../results/glmdata_by_study/", studynm, ".csv", sep=""))
  }

  ###########################################################
  ## Step 7: Process pig-specific crop covariates
  ###########################################################

  # For crop use, pigs have individual specific covariates as they tend to use
  # certain portions of crop fields.  This section of the script build these
  # pig-specific crop-use covariates.

  # Make dummy columns for crop distance covariates: location and directional
  croptypes = anal_params$croptypes
  for(ctype in croptypes){

    labloc = paste0(ctype, "dists_loc")
    labgrad = paste0(ctype, "dists_grad")

    fullglm[, (labloc):=NA]
    fullglm[, (labgrad):=NA]
  }

  # Make pig-specific crop covariates for each pig
  for(pignm in unique(fullglm$pigID)) {

    # Get pig-specific data
    cat("Working on pig", pignm, "\n")
    pigind = fullglm$pigID == pignm
    tdat = fullglm[pigind]

    # For each crop type, distance-to-used metric
    for(ctype in croptypes){

      # Identify where the pig has been 
      incrop = tdat[(tdat[['z']] == 1) & (tdat[[paste0(ctype, "_loc")]] == 1)]

      # If it has been in crops, compute distance to nearest gradient to crops
      # Otherwise NA
      if(nrow(incrop) != 0) {

        pts = incrop[, list(x.adj, y.adj)]
        spts = SpatialPoints(pts, 
                      proj4string=crs("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

        # Rasterize in crop points
        basefiles = Sys.glob(file.path(anal_params$covariate_dir, "croplayer", studynm,
                  paste(studynm, anal_params$baserastype, "*.tif", sep="")))

        buffer = 0
        extobj = extent(c(xmin=min(tdat$x.adj) - buffer, 
                        xmax=max(tdat$x.adj) + buffer, 
                        ymin=min(tdat$y.adj) - buffer, 
                        ymax=max(tdat$y.adj) + buffer))
        baseras = crop(raster(basefiles[1]), extobj)

        tras = baseras
        values(tras) = 0
        inras = rasterize(spts, tras)
        values(inras)[!is.na(values(inras))] = 1

        # Convert in-crop raster to polygons using GDAL
        writeRaster(inras, 'temp.tif', format="GTiff", overwrite=F)
        system("gdal_polygonize.py temp.tif -f 'ESRI Shapefile' temp.shp")
        shp = shapefile("temp.shp")
        system("rm temp.*") # Remove after conversion

        # Distance from nearest used crop center.
        centers = gCentroid(shp, byid=T)
        mindists = apply(spDists(rasterToPoints(tras, spatial=T), centers, longlat=T), 1, min)

        # Make distance to nearest raster
        distras = raster()
        extent(distras) = extent(baseras)
        crs(distras) = crs(baseras)
        res(distras) = res(baseras)
        values(distras) = mindists


        # Compute distance gradient 
        grad = rast.grad(distras)

        # Extract gradient values pig location values
        tdatpts = SpatialPoints(tdat[, list(x.current, y.current)], 
                          proj4string=crs("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
        xgrad = extract(grad$rast.grad.x,  tdatpts) 
        ygrad = extract(grad$rast.grad.y, tdatpts)
        gradval = cbind(xgrad, ygrad)

        # Based on how the data is formatted
        #   There are pairs of for that read: west (-1, 0), east (1, 0), north (0, 1), south (0, -1). 
        xstep = rep(c(-1, 1, 0, 0), nrow(tdat) / 4)
        ystep = rep(c(0, 0, 1, -1), nrow(tdat) / 4)
        stepvect = cbind(xstep, ystep)

        # Product of gradient and step direction
        cropdists_grad = (xgrad * xstep) + (ygrad * ystep)
        cropdists_loc = extract(distras, tdatpts)
      }
      else {
        cropdists_grad = array(NA, dim=nrow(tdat))
        cropdists_loc = array(NA, dim=nrow(tdat))
      }

      # Save the pig specific results
      labloc = paste0(ctype, "dists_loc")
      labgrad = paste0(ctype, "dists_grad")
      fullglm[[labloc]][pigind] = cropdists_loc
      fullglm[[labgrad]][pigind] = cropdists_grad

    } # End croptype

  } # End pigID 


  ###########################################################
  ## Step 8: Write the file to disk
  ###########################################################

  fwrite(fullglm, file=paste("../results/glmdata_by_study/", studynm, ".csv", sep=""))
  cat("Completed analysis for", studynm, "\n")

}

# End logging to file
sink()




