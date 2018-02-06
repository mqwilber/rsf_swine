## Functions for analyzing feral swine movement data

library(ctmcmove)
library(geosphere)
library(raster)
library(lubridate)
library(fda)
library(parallel)
library(glmnet)

parallel_pigs = function(pid, alldata, ...){

  args = list(...)
  pigdata = alldata[pigID == pid, ]

}

fit_ctmcmodel = function(pigdata, pigID, loc.stack, grad.stack, timestep="15 mins", 
                                  method="interp", buffer=0.001,
                                  n.mcmc=400, df=NULL, impute=1, 
                                  mc.cores=1, sigma.fixed=NA){
  # Fits a continuous movement path and generates glm data
  #
  # Parameters
  # ----------
  # pigdata : data.frame/data.table
  #     Columns datetime, x (longitude), y (longitude),
  # pigID : str
  #   A unique identifier for a pig.
  # loc.stack : RasterStack
  #   Location rasters
  # grad.stack : RasterStack
  #   Gradient rasters
  # timestep : str
  #   e.g. "15 mins"
  # method: str
  #   Specifies how to compute the continuous path.
  #   Options are:
  #     "interp": Linear interpolation along the path
  #     "bspline_freq": Cubic B-spline without penalization in a frequentist 
  #       framework.
  #     "bspline_bayes": Cubic B-spline with penalization in a Bayesian framework.
  #         Calls ctmcmove::mcmc.fmove
  # impute : int
  #   The number of continuous paths to impute. If "iterp == 1"
  # df : int
  #   Degrees of freedom of the B-spline
  # n.mcmc : int
  #     Number of MCMC samples if Bayesian approach is used
  # buffer : float
  #     Buffer the cropping of the rasters.  Particularly important if using 
  #     a b-spline fit.

  # Get the bounding box around observations and crop raster
  ymin = min(pigdata$y)
  ymax = max(pigdata$y)
  xmin = min(pigdata$x)
  xmax = max(pigdata$x)

  # Add on a dummy raster to loc and grad to prevent column error in ctmc2glm
  dumgrad = unstack(grad.stack)[[1]]
  names(dumgrad) = "dummy_grad"
  grad.stack = stack(list(grad.stack, dumgrad))
  names(dumgrad) = "dummy_loc"
  loc.stack = stack(list(loc.stack, dumgrad))

  # Compute continuous pig path
  predPath = continuous_path(pigdata, timestep, method, n.mcmc=n.mcmc, df=df, 
                                  impute=impute, sigma.fixed=sigma.fixed)

  # Specify where animals can move on the landscape
  # move_field = cras
  # values(move_field) = 1
  # trans = transition(move_field, prod, 4)

  # Make this break more intelligently. 
  imputepath = function(x){

    cat("Working on imputation ", x, " of ", impute, "\n")
    path = list(xy=predPath$path[[x]], t=as.vector(predPath$time))

    # Can be slow...need to be smart about which cells to consider.
    ctmc = path2ctmc(xy=path$xy, t=path$t, rast=grad.stack)

    tglm_data = ctmc2glm(ctmc, loc.stack, grad.stack)
    tglm_data$pigID = pigID
    tglm_data$imputeID = x # Number the imputed distributions
    return(tglm_data)
  }

  if(method == "interp")
      impute = 1

  if(mc.cores > 1)
    full_glmdat = mclapply(1:impute, imputepath, mc.cores=mc.cores)
  else
    full_glmdat = lapply(1:impute, imputepath)


  # full_glmdat = list()
  # for(i in 1:length(predPath$path)){

  #   cat("Working on imputation ", i, " of ", impute, "\n")
  #   path = list(xy=predPath$path[[i]], t=as.vector(predPath$time))

  #   # Can be slow...need to be smart about which cells to consider.
  #   ctmc = path2ctmc(xy=path$xy, t=path$t, rast=cgrad.stack)

  #   tglm_data = ctmc2glm(ctmc, cloc.stack, cgrad.stack)
  #   tglm_data$ID = unique(pigdata$pigID)
  #   full_glmdat[[i]] = tglm_data
  # }

  #return(ctmc)
  return(do.call(rbind, full_glmdat))

}


continuous_path = function(data, timestep, method, impute=1, df=NULL,
                      n.mcmc=400, sigma.fixed=NA){
  # Get predictions for a continuous path from data
  #
  # Parameters
  # ----------
  # data : data.table/data.frame
  #   Columns "x" (e.g. longitude), "y" (e.g. latitude), and "datetime"
  # timestep : str
  #   e.g. "15 mins"
  # method: str
  #   Specifies how to compute the continuous path.
  #   Options are:
  #     "interp": Linear interpolation along the path
  #     "bspline_freq": Cubic B-spline without penalization in a frequentist 
  #       framework.
  #     "bspline_bayes": Cubic B-spline with penalization in a Bayesian framework.
  #         Calls ctmcmove::mcmc.fmove
  # impute : int
  #   The number of continuous paths to impute. If "iterp == 1"
  # df : int
  #   Degrees of freedom of the B-spline
  # n.mcmc : int
  #     Number of MCMC samples if Bayesian approach is used
  #
  # Returns
  # -------
  # : list
  #     $predPaths : A list of impute paths (x, y matrices)
  #     $time : The times that were imputed.

  if(is.null(df))
    df = nrow(data) / 2 # Set the number of knots relative to data length

  mintime = min(data$datetime)
  maxtime = max(data$datetime)
  predTime = seq(mintime, maxtime, by=timestep)
  traw = as.numeric(data$datetime)
  t = traw - mean(traw)
  tpred = as.numeric(predTime) - mean(traw)

  # For unscaling
  mux = mean(data$x); sdx = sd(data$x); muy = mean(data$y); sdy = sd(data$y)

  if(method == "interp"){

    iterpX = approxfun(data$datetime, data$x) # e.g. Longitude
    iterpY = approxfun(data$datetime, data$y) # e.g. Latitude

    predX = iterpX(predTime)
    predY = iterpY(predTime)

    predPaths = list()
    predPaths[[1]] = cbind(x=predX, y=predY)

  } else if(method == "bspline_bayes"){

    xy = as.matrix(data[, list(scale(x), scale(y))])

    # Calculate bspline knots
    bsp = bs(t, df=df)
    knots = attributes(bsp)$knots
    basisfxn = create.bspline.basis(c(min(t), max(t)), breaks=knots, norder=4)

    fit = mcmc.fmove(xy, t, basisfxn, tpred, QQ="CAR", n.mcmc=n.mcmc, 
                  num.paths.save=impute, sigma.fixed = sigma.fixed)

    predPaths = lapply(1:length(fit$pathlist), function(x){
                                        xy = fit$pathlist[[x]]$xy;
                                        x = xy[, 1]*sdx + mux;
                                        y = xy[, 2]*sdy + muy;
                                        return(cbind(x=x, y=y))})

  } else if(method == "bspline_freq"){


    bsp = bs(t, df=df)
    knots = attributes(bsp)$knots

    basisfxn = create.bspline.basis(c(min(t), max(t)), breaks=knots, norder=4)
    Xmat = Matrix(eval.basis(t, basisfxn), sparse=TRUE)

    fitlong = glmnet(Xmat, scale(data$x), alpha=0, lambda=0, standardize = F, intercept=F) # No regularization
    fitlat = glmnet(Xmat, scale(data$y), alpha=0, lambda=0, standardize = F, intercept=F) # No regularization
    Betas1 = as.vector(fitlong$beta)
    Betas2 = as.vector(fitlat$beta)
    residlong = scale(data$x) - predict(fitlong, newx=Xmat) 
    residlat = scale(data$y) - predict(fitlat, newx=Xmat) 
    print("done fitting")

    # Assume same latlong variance for prediction
    sigma2 = mean(c(sum(residlong^2) / (fitlong$nobs - fitlong$df), 
                    sum(residlat^2) / (fitlat$nobs - fitlat$df)))
    print(sigma2)

    # This is brutally slow for large matrices.  Does not scale well...
    Vmat = sigma2*chol2inv(chol(t(Xmat) %*% Xmat))
    L = chol(Vmat) 
    print("done with inversion")

    # Very slow...but converting to the sparse matrix really speeds up imputation
    Xpred = eval.basis(tpred, basisfxn)
    print("done with basis pred")
    #print(class(Xpred))

    predPaths = list()
    for(i in 1:impute){
  
      predX = Xpred %*% (Betas1 + (L %*% rnorm(length(Betas1))))
      predY = Xpred %*% (Betas2 + (L %*% rnorm(length(Betas2))))
      predPaths[[i]] = cbind(x=as.vector(predX*sdx + mux), y=as.vector(predY*sdy + muy))
  
    }
    
  } else{
    stop(paste(method, "is not a known method.\nOptions are 'interp', 'bspline_freq', 
            'bspline_bayes'"))
  }

  return(list(paths=predPaths, time=predTime))
}


getR = function (object, stack.static, stack.grad, normalize.gradients = FALSE, 
                 grad.point.decreasing = TRUE, directions = 4, zero.idx = integer(), 
                 coef){
  # same as ctmcmove::get.rate.matrix, but allows you to use dummy variables. ctmcmove has a bug along these lines.
  if (inherits(stack.static, "Raster")) {
    p.static = nlayers(stack.static)
  }
  else p.static = 0
  p.crw = 0
  if (inherits(stack.grad, "Raster")) {
    p.grad = nlayers(stack.grad)
    stack.gradient = rast.grad(stack.grad)
    if (normalize.gradients) {
      lengths = sqrt(stack.gradient$grad.x^2 + stack.gradient$grad.y^2)
      stack.gradient$grad.x <- stack.gradient$grad.x/lengths
      stack.gradient$grad.y <- stack.gradient$grad.y/lengths
    }
  }
  else {
    p.grad = 0
  }
  p = p.static + p.crw + p.grad
  if (class(stack.static) == "RasterStack") {
    examplerast = stack.static[[1]]
  }
  if (class(stack.static) == "RasterLayer") {
    examplerast = stack.static
  }
  locs = 1:ncell(examplerast)
  nn = ncell(examplerast)
  R = Matrix(0, nrow = nn, ncol = nn, sparse = TRUE)
  adj = adjacent(examplerast, locs, pairs = TRUE, sorted = TRUE, 
                 id = TRUE, directions = directions)
  idx.mot = adj[, 2:3]
  if (p.static > 0) {
    X = data.frame(values(stack.static)[idx.mot[, 1], ])
    for(nm in names(X)){
      X[[nm]] = as.numeric(X[[nm]])
    }
    
  }
  start.cells = idx.mot[, 1]
  adj.cells = idx.mot[, 2]
  xy.cell = xyFromCell(examplerast, start.cells)
  xy.adj = xyFromCell(examplerast, adj.cells)
  v.adj = (xy.adj - xy.cell)/sqrt(apply((xy.cell - xy.adj)^2, 
                                        1, sum))
  if (p.grad > 0) {
    X.grad = v.adj[, 1] * stack.gradient$grad.x[start.cells, 
                                                ] + v.adj[, 2] * stack.gradient$grad.y[start.cells, 
                                                                                       ]
    if (grad.point.decreasing == TRUE) {
      X.grad = -X.grad
    }
    X = cbind(X, X.grad)
  }
  X$tau = 1
  X$crw = 0
  if (missing(coef)) {
    R[idx.mot] = as.vector(predict(object, newdata = X, type = "response"))
  }
  else {
    if (length(coef) != length(coefficients(object))) 
      stop("'coef' vector is not the correct length!")
    R[idx.mot] = exp(X %*% coef)
  }
  if (length(zero.idx) > 0) {
    R[zero.idx, zero.idx] = 0
  }
  R
}


get_direction = function(current, indices){
  # Given a vector of cell indices, determine if a cell is above (1), below (-1), left (-2) or right (2)
  #
  # Parameters
  # ----------
  # current : index of the current cell
  # indices : indices of the four rook-based cells
  #
  # Returns
  # -------
  # Directions of each indices relative to current
  
  diffs = current - indices
  direction = array(NA, dim=length(diffs))
  
  direction[diffs > 1] = 1 # north
  direction[diffs == -1] = 2 # east
  direction[diffs < -1] = -1 # south
  direction[diffs == 1] = -2 # west
  
  return(direction)
}

set_direction_weights = function(past_direction, next_directions){
  # Weight the direction of movement based on a correlated random walk.
  #
  # Following Hank et al. 2015, this is taking the dot product between a unit vectors pointing
  # either up, down, left or right.  Same direction = 1, Opposite direction = -1, Else = 0.
  #
  # Parameters
  # ----------
  # past_direciton : either 1, 2, -1, -2 (see get_direction)
  # next_directions : from get_direction
  
  # Returns
  # -------
  # next_direction : the weights on the next directions (1, -1, 0, 0)
  
  weights = array(0, dim=length(next_directions))
  weights[past_direction == next_directions] = 1
  weights[past_direction == -1*next_directions] = -1
  return(weights)
  
}



sim_traj = function(R, fit, start, steps=100){
  # Simulate the discrete time steps of animal movement based on the CTMC 
  # generator matrix R, and the the GLM model.  
  # 
  # At the moment, the model can't have time dependent covariates.
  #
  # Parameters
  # ----------
  # R : CMTC generators matrix without rates on the diagonals (from get.rate.matrix())
  # fit : Fittted GLM object, with crw as a coefficient.
  # start : index of starting cell.  Should be between 1:nrow(R)
  # steps: Number of steps in the simulation
  #
  # Returns
  # -------
  # list of position indices, interevent times, and step_directions
  
  
  # Build the transition matrix
  zero.idx = which(colSums(R) == 0)
  n.0 = length(zero.idx)
  if (length(zero.idx) > 0) {
    R = R[-zero.idx, -zero.idx]
  }
  n = nrow(R)
  one = matrix(1, n, 1)
  # P = R
  # 
  # # Discrete time transition matrix
  # P = Diagonal(x = 1/as.numeric(P %*% one)) %*% P
  
  # 
  positions = array(NA, dim=steps + 1)
  step_direction = array(NA, dim=steps + 1)
  time = array(NA, dim=steps + 1)
  time[1] = 0
  
  step_direction[1] = 1 # Start moving north
  positions[1] = start
  for(i in 2:(steps + 1)){
    
    
    indices = which(R[positions[i - 1], ] != 0)
    
    # Compute CRW weights
    directions = get_direction(positions[i - 1], indices) # Directions of indices
    weights = set_direction_weights(step_direction[i- 1], directions)
    crw = exp(fit$coefficients['crw']*weights)
    
    rates = R[positions[i - 1], indices]*crw
    time[i] = rexp(1, sum(rates)) # Interevent time
    probs = rates / sum(rates)
    samp = rmultinom(1, 1, probs)
    
    positions[i] = indices[as.logical(samp)]
    step_direction[i] = directions[as.logical(samp)]
    
  }
  
  return(list(positions=positions, steps=step_direction, time=time))
  
}


process_covariates = function(locvars, studynm, extobj, 
          mindate, maxdate, ext,
          cov_path="/Users/mqwilber/Repos/rsf_swine/data/covariate_data",
          timevar=c("temperature", "ndvi")){
	# Compiles lists of rasters to be used as location and gradient covariates
  # in analyses.
  #
  # Parameters
  # ----------
  # locvars : vector of strings
  #   The variables to be used as locations covariates. They must have a matching
  #   folder in cov_path.
  # studynm : str
  #   Name of the study under consideration
  # extobj : Extent object fro raster package
  #   Gives the extent to which to crop the rasters
  # mindate : POSIXct object
  #   Minimum date
  # maxdate : POSIXct object
  #   Maximum date
  # cov_path : Path to the covariate
  #   File path to where the covariate data is stored.
  # timevar : vector of strings
  #   Vector contains the covariates that are time varying.
  # ext : str
  #   To be appended to rasters. Either "loc" (location covariates) or "grad" 
  #   (gradient covariates)
  #
  # Returns
  # -------
  # : list of rasters
  #
  # Notes
  # -----
  # TODO: While all the rasters are cropped in the same way, they are not projected 
  # identically in this function.

  baseras = 
  loclist = list()

  # Extract location covariates
  for(i in 1:length(locvars)){

    loccov = locvars[i]

    if(loccov %in% timevar){ # If there is a time varying covariate

      day(mindate) = 1 
      day(maxdate) = 1
      dates = seq(as.Date(mindate), as.Date(maxdate), by="month")
      months = month(dates)
      years = year(dates)

      # Extract time-dependent rasters
      for(j in 1:length(months)){

        fp = file.path(cov_path, loccov, studynm, paste(studynm, "_", loccov, "_", months[j], "_", years[j], ".grd", sep=""))
        tras = raster(fp)
        loclist[[paste(loccov, months[j], years[j], ext, sep="_")]] = crop(tras, extobj)
      }

    } else {
      # Extract time-independent raster
      tras = raster(file.path(cov_path, loccov, studynm, paste(studynm, "_", loccov, ".grd", sep="")))
      loclist[[paste(loccov, ext, sep="_")]] = crop(tras, extobj)
    }
  }

  return(loclist)

  
}

###### Helper functions ##########


get_distance_vect = function(longlat){
  # Calculate distances between sequential longlat points 
  #
  # Parameters
  # ----------
  # longlat : data.frame with columns longitude and latitude
  #
  # Returns
  # -------
  # : vector of distances
  
  longs = longlat$longitude
  lats = longlat$latitude
  
  # Match distance vectors
  v1 = cbind(longs[-length(longs)], lats[-length(lats)])
  v2 = cbind(longs[-1], lats[-1])
  
  tuples = lapply(1:nrow(v1), function(i) list(v1[i, ], v2[i, ]))
  distvect = lapply(tuples, function(x) distm(x=x[[1]], y=x[[2]], fun=distHaversine))
  # 
  return(do.call(rbind, distvect))
}


runs = function(x, ctime=60, clength=10, ctimemin=0){
  # Given a vector of time differneces X
  # determine if the vector has a run of 
  # clength or greater of diffs < ctime (in minutes). 
  # Parameters
  # ----------
  # x : vector of datetimes
  # ctime : Diffs must be less than this time
  # clength : Length of a run
  
  # Notes
  # -----
  # This will help identify potential trajectories to use for the movement study
  # Returns TRUE or FALSE
  
  deltat = diff(x)
  units(deltat) = "mins"
  inds = (deltat < ctime) & (deltat > ctimemin) # String of booleans
  
  # Count runs
  run_vect = rle(inds)
  res = any((run_vect$values == TRUE) & (run_vect$lengths >= clength))
  return(res)
  
}

