## Functions for analyzing feral swine movement data


list.of.packages <- c("ctmcmove", "geosphere", "raster", "lubridate", "fda", 
                        "parallel", "glmnet", "splines", "fda", "mgcv", "suncalc")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only=TRUE)


fit_ctmcmodel = function(pigdata, pigID, loc.stack, grad.stack, timestep="15 mins", 
                                  method="interp", buffer=0.001,
                                  n.mcmc=400, df=NULL, impute=1, 
                                  mc.cores=1, sigma.fixed=NA,
                                  path2ctmcMethod="ShortestPath",
                                  xygrad=FALSE, xgrad.stack=NULL, 
                                  ygrad.stack=NULL,
                                  grad.point.decreasing=FALSE, 
                                  predPath=NULL, directions=4){
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
  #   e.g. "15 mins",  specifies the equal length time interval to impute.
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
  # path2ctmcMethod: Specifies interpolation method.  Either "ShortestPath", which
  #   uses the shortest graphical path on the raster graph, or
  #   "LinearInterp", which linearly interpolates between observed
  #   locations.  "ShortestPath" is slower, slightly more accurate,
  #   and allows for impassible barriers specified through
  #   "zero.idx".  "LinearInterp" is faster but does not allow for
  #   impassible barriers.
  # xygrad : bool
  #   If TRUE expects, xgrad.stack and ygrad.stack arguments that contain RasterStacks
  #   of pre-computed gradient covariates.  Note that ctmcmove will compute most
  #   gradients for you and this is for covariates such as a gradient pointing 
  #   toward the nearest cropland. 
  # xgrad.stack : Null or RasterStack
  #   RasterStack containing x dimension of covariates
  # ygrad.stack : Null or RasterStack
  #   RasterStack containing y dimension of covariates
  # grad.point.decreasing: Bool
  #   If TRUE, take the negative gradient.

  # Get the bounding box around observations and crop raster
  # ymin = min(pigdata$y)
  # ymax = max(pigdata$y)
  # xmin = min(pigdata$x)
  # xmax = max(pigdata$x)

  # Add on a dummy raster to loc and grad to prevent column error in ctmc2glm
  dumgrad = unstack(grad.stack)[[1]]
  names(dumgrad) = "dummy_grad"
  grad.stack = stack(list(grad.stack, dumgrad))
  names(dumgrad) = "dummy_loc"
  loc.stack = stack(list(loc.stack, dumgrad))

  # Compute continuous pig path

  if(is.null(predPath)){
    predPath = continuous_path(pigdata, timestep, method, n.mcmc=n.mcmc, df=df, 
                                  impute=impute, sigma.fixed=sigma.fixed)
  }

  # Specify where animals can move on the landscape
  # move_field = cras
  # values(move_field) = 1
  # trans = transition(move_field, prod, 4)

  # TODO: Make this break more intelligently. 
  imputepath = function(x){

    cat("Computing imputation", x, "of", impute, "for", pigID, "\n")
    path = list(xy=predPath$path[[x]], t=as.vector(predPath$time))

    # Can be slow...need to be smart about which cells to consider.
    # LinearInterp increases speed, but can't include impassable cells.
    ctmc = path2ctmc(xy=path$xy, t=path$t, rast=grad.stack, 
                method=path2ctmcMethod)

    if(xygrad){
      tglm_data = ctmc2glm_wgrad(ctmc, loc.stack, grad.stack, xgrad.stack, ygrad.stack, 
                        grad.point.decreasing=grad.point.decreasing,
                        directions=directions)
    }
    else{
      tglm_data = ctmc2glm(ctmc, loc.stack, grad.stack, 
                          grad.point.decreasing=grad.point.decreasing,
                          directions=directions)
    }

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
  # sigma.fixed : float
  #   Fix the observation error in the MCMC to sigma.fixed. If NA, this
  #   parameter will be estimated.
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
    cat("done fitting\n")

    # Assume same latlong variance for prediction
    sigma2 = mean(c(sum(residlong^2) / (fitlong$nobs - fitlong$df), 
                    sum(residlat^2) / (fitlat$nobs - fitlat$df)))

    # This is brutally slow for large matrices.  Does not scale well...
    Vmat = sigma2*chol2inv(chol(t(Xmat) %*% Xmat))
    L = chol(Vmat) 
    cat("done with inversion\n")

    # Very slow...but converting to the sparse matrix really speeds up imputation
    Xpred = eval.basis(tpred, basisfxn)
    cat("done with basis pred\n")
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
                 coef, xygrad=TRUE, stack.xgrad=NULL, stack.ygrad=NULL){
  # same as ctmcmove::get.rate.matrix, but allows you to use dummy variables. ctmcmove has a bug along these lines.
  if (inherits(stack.static, "Raster")) {
    p.static = nlayers(stack.static)
  }
  else p.static = 0
  p.crw = 0
  if (inherits(stack.grad, "Raster")) {

    # TODO: Add in additional gradient vectors.
    if(!xygrad){
      p.grad = nlayers(stack.grad)
      stack.gradient = rast.grad(stack.grad)
    } else{
      p.grad = nlayers(stack.grad) + nlayers(stack.xgrad)
      stack.gradient = rast.grad(stack.grad)

      # Add on the additional gradient covariates
      stack.gradient$rast.grad.x = stack(list(stack.gradient$rast.grad.x, stack.xgrad))
      stack.gradient$rast.grad.y = stack(list(stack.gradient$rast.grad.y, stack.ygrad))
      stack.gradient$grad.x = cbind(do.call(cbind, lapply(unstack(stack.xgrad), function(x) values(x))), stack.gradient$grad.x)
      stack.gradient$grad.y = cbind(do.call(cbind, lapply(unstack(stack.ygrad), function(x) values(x))), stack.gradient$grad.y)
      colnames(stack.gradient$grad.x)[1:length(names(stack.xgrad))] = names(stack.xgrad)
      colnames(stack.gradient$grad.y)[1:length(names(stack.ygrad))] = names(stack.ygrad)
    }

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
    # print(X)
    a = as.vector(predict(object, newdata = X, type = "response"))
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


ctmc2glm_wgrad = function (ctmc, stack.static, stack.grad, 
    stack.xgrad, stack.ygrad, crw = TRUE, 
    normalize.gradients = FALSE, 
    grad.point.decreasing = TRUE, include.cell.locations = TRUE, 
    directions = 4, zero.idx = integer()) 
{
    # Same as `ctmc2glm` in the ctmcmove package, but allows for the option of 
    # including a precomputed gradient layer

    p.static = nlayers(stack.static)
    p.crw = 0
    if (crw) {
        p.crw = 1
    }
    if (class(stack.grad) == "RasterLayer" | class(stack.grad) == 
        "RasterStack") {
        p.grad = nlayers(stack.grad) + nlayers(stack.xgrad)
        stack.gradient = rast.grad(stack.grad)

        # Add on the additional gradient covariates
        stack.gradient$rast.grad.x = stack(list(stack.gradient$rast.grad.x, stack.xgrad))
        stack.gradient$rast.grad.y = stack(list(stack.gradient$rast.grad.y, stack.ygrad))
        stack.gradient$grad.x = cbind(do.call(cbind, lapply(unstack(stack.xgrad), function(x) values(x))), stack.gradient$grad.x)
        stack.gradient$grad.y = cbind(do.call(cbind, lapply(unstack(stack.ygrad), function(x) values(x))), stack.gradient$grad.y)
        colnames(stack.gradient$grad.x)[1:length(names(stack.xgrad))] = names(stack.xgrad)
        colnames(stack.gradient$grad.y)[1:length(names(stack.ygrad))] = names(stack.ygrad)

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
    locs = ctmc$ec
    wait.times = ctmc$rt
    notzero.idx = 1:ncell(examplerast)
    if (length(zero.idx) > 0) {
        notzero.idx = notzero.idx[-zero.idx]
    }
    adj = adjacent(examplerast, locs, pairs = TRUE, sorted = TRUE, 
        id = TRUE, directions = directions, target = notzero.idx)
    adj.cells = adj[, 3]
    rr = rle(adj[, 1])
    time.idx = rep(rr$values, times = rr$lengths)
    start.cells = adj[, 2]
    z = rep(0, length(start.cells))
    idx.move = rep(0, length(z))
    diag.move = rep(0, length(locs))
    for (i in 1:(length(locs))) {
        idx.t = which(time.idx == i)
        idx.m = which(adj.cells[idx.t] == locs[i + 1])
        z[idx.t[idx.m]] <- 1
        if (length(idx.m) == 0) {
            diag.move[i] = 1
        }
    }
    tau = rep(wait.times, times = rr$lengths)
    t = rep(ctmc$trans.times, times = rr$lengths)
    if (nlayers(stack.static) > 1) {
        X.static = values(stack.static)[start.cells, ]
    }
    else {
        X.static = matrix(values(stack.static)[start.cells], 
            ncol = 1)
    }
    colnames(X.static) <- names(stack.static)
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
        colnames(X.grad) <- colnames(stack.gradient$grad.x)
    }
    idx.move = which(z == 1)
    idx.move = c(idx.move, length(z))
    v.moves = v.adj[rep(idx.move[1:(length(rr$lengths) - 1)], 
        times = rr$lengths[-1]), ]
    v.moves = rbind(matrix(0, ncol = 2, nrow = rr$lengths[1]), 
        v.moves)
    X.crw = apply(v.moves * v.adj, 1, sum)
    if (crw == FALSE & p.grad > 0) {
        X = cbind(X.static, X.grad)
    }
    if (crw == TRUE & p.grad > 0) {
        X = cbind(X.static, X.grad, X.crw)
        colnames(X)[ncol(X)] = "crw"
    }
    if (crw == FALSE & p.grad == 0) {
        X = cbind(X.static)
    }
    if (crw == TRUE & p.grad == 0) {
        X = cbind(X.static, X.crw)
        colnames(X)[ncol(X)] = "crw"
    }
    if (include.cell.locations) {
        xys = cbind(xy.cell, xy.adj)
        colnames(xys) = c("x.current", "y.current", "x.adj", 
            "y.adj")
        X = cbind(X, xys)
    }
    T = length(wait.times)
    p = ncol(X)
    out = data.frame(z = z, X, tau = tau, t = t)
    T = nrow(out)
    out = out[-((T - 3):T), ]
    out
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
  # generator matrix R, and the GLM model.  
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


process_covariates = function(locvars, studynm, baseras, 
          mindate, maxdate, ext,
          cov_path="/Users/mqwilber/Repos/rsf_swine/data/covariate_data",
          timevar=c("temperature", "ndvi", "precipitation", "snowdepth"),
          croptypes=c("beverage_and_spice", "cereals", "oilseed", "other", 
                      "fruit_and_nuts", "grasses", "leguminous",
                      "root_and_tuber", "sugar", "tobacco", 
                      "vegetables_and_melons"),
          landcovertypes=c("barren_land", "canopycover", "deciduous_forest", 
                           "developed", "evergreen_forest", "mixed_forest", 
                           "wetland"),
          nngrad=FALSE,
          distgrad=FALSE, distgrad_types=numeric(), decay=1, 
          projectionMethod="bilinear"){
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
  # baseras : RasterLayer
  #   Base raster layer to which all other rasters will be projected.
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
  # distgrad : bool
  #   If True, for croptype variables the distance-to-gradient is computed. This
  #   is a global distance measure that simultaneously accounts for 
  #   the distance between all patch types
  # nngrad : bool
  #   If True, extracts the nearest-neighbor gradient raster.  Only applicable
  #   for landcover and croplayer types.
  # distgrad_types : vector
  #   Include croptypes for which you DO want a gradient computed.
  # croptypes : vector
  #   Categories within croplayer locvar.
  # landcovertypes: vector
  #   Categories within landcover locvar.
  # projectionMethod : string
  #   Either "bilinear" or "ngb" as described in projectRaster fxn.
  # decay : float
  #   For distgrad, the level of distance to decay to include when computing 
  #   distance-to-resource gradient.
  #
  # Returns
  # -------
  # : list of rasters
  #   All rasters are projected in the same way


  loclist = list()

  # Extract location covariates
  for(i in 1:length(locvars)){

    loccov = locvars[i]
    cat("Processing", loccov, "\n")
    if(loccov == "croplayer"){
      # This is a special covariate as it is broken into four covariates and is 
      # is on the yearly scale
      minyear = year(mindate)
      maxyear = year(maxdate)

      for(yr in minyear:maxyear){
        for(ctype in croptypes){

          if(!nngrad){

            fp = file.path(cov_path, loccov, studynm, paste(studynm, "_", ctype, 
                                            "_", yr, ".tif", sep=""))
            cat(fp, "\n")
            tras = raster(fp)

            # Compute distance2 gradients for croptypes. 
            if(distgrad){ 

              fp = file.path(cov_path, loccov, studynm, paste(studynm, "_", ctype, 
                                            "_", yr, ".tif", sep=""))
              tras = raster(fp)

              if(ctype %in% distgrad_types){
                cat("Processing gradient for", ctype, "\n")
                loclist[[paste(ctype, yr, ext, sep="_")]] = 
                          get_distance_to_gradient(projectRaster(tras, baseras, 
                              method=projectionMethod), decay=decay)
              }
            } else
              loclist[[paste(ctype, yr, ext, sep="_")]] = 
                            projectRaster(tras, baseras, method=projectionMethod)
          } else{

            fp = file.path(cov_path, loccov, studynm, paste(studynm, "_", ctype, 
                                              "_", yr, ".tif", sep=""))
            cat(fp, "\n")
            tras = raster(fp)
            loclist[[paste(ctype, yr, ext, sep="_")]] = 
                            projectRaster(tras, baseras, method=projectionMethod)
            

          }

        }
      }

    } else if(loccov %in% timevar){ # If there are time varying covariates

      day(mindate) = 1 
      day(maxdate) = 1
      dates = seq(as.Date(mindate), as.Date(maxdate), by="month")
      months = month(dates)
      years = year(dates)

      # Extract time-dependent rasters
      for(j in 1:length(months)){

        fp = file.path(cov_path, loccov, studynm, paste(studynm, "_", loccov, 
                                  "_", months[j], "_", years[j], ".tif", sep=""))
        cat(fp, "\n")
        tras = raster(fp)
        loclist[[paste(loccov, months[j], years[j], ext, sep="_")]] = 
                                              projectRaster(tras, baseras,
                                                         method=projectionMethod)
      }

    } else {
      # Extract time-independent rasters

      if(loccov == "landcover"){

        for(ltype in landcovertypes){

          if(!nngrad){

            if(ltype == "canopycover"){
              fp = file.path(cov_path, loccov, studynm, paste(studynm, "_", 
                                                          ltype, ".tif", sep=""))
            } else{
              # NOTE: Only have wetland and developed rasters here.
              fp = file.path(cov_path, loccov, studynm, paste(studynm, "_", 
                                                          ltype, "_nndistance.tif", sep=""))
            }

            cat(fp, "\n")
            tras = raster(fp)

            if(distgrad){ 

              fp = file.path(cov_path, loccov, studynm, paste(studynm, "_", 
                                                          ltype, ".tif", sep=""))
              tras = raster(fp)
              if(ltype %in% distgrad_types){
                cat("Processing gradient for", ltype, "\n")
                loclist[[paste(ltype, ext, sep="_")]] = 
                          get_distance_to_gradient(projectRaster(tras, baseras, 
                                            method=projectionMethod), decay=decay)
              } 
            } else
              loclist[[paste(ltype, ext, sep="_")]] = projectRaster(tras, baseras, 
                                                          method=projectionMethod)
          } else{

            if(ltype != "canopycover"){

              fp = file.path(cov_path, loccov, studynm, 
                        paste(studynm, "_", ltype, "_nndistance.tif", sep=""))
              cat(fp, "\n")
              tras = raster(fp)
              loclist[[paste(ltype, ext, sep="_")]] = projectRaster(tras, baseras, 
                                                          method=projectionMethod)
            } else{

              fp = file.path(cov_path, loccov, studynm, paste(studynm, "_", 
                                                           ltype, ".tif", sep=""))
              cat(fp, "\n")
              tras = raster(fp)
              loclist[[paste(ltype, ext, sep="_")]] = projectRaster(tras, baseras, 
                                                          method=projectionMethod)
            }
          }
        } # End for loop
      } else {

       
        fp = file.path(cov_path, loccov, studynm, paste(studynm, "_", 
                                                        loccov, ".tif", sep=""))

        # Check it a nearest neighbor distance raster exists. If so, use it.
        # If not use the base file.
        if(loccov == "water"){
          fnnp = file.path(cov_path, loccov, studynm, paste(studynm, "_", 
                                              loccov, "_nndistance.tif", sep=""))
          if(file.exists(fnnp)){
            fp = fnnp
          }
        }

        cat(fp, "\n")
        tras = raster(fp)
        loclist[[paste(loccov, ext, sep="_")]] = projectRaster(tras, baseras, 
                                                        method=projectionMethod)
      } # End else
    } # End else
  } # End for loop

  return(loclist)

  
}

add_covariates = function(fullglm, covariates, studynm, 
        basedir="/Users/mqwilber/Repos/rsf_swine/data/covariate_data"){
  # Adds various population-level covariates onto the data.table after 
  # the initial CTMC processing.
  #
  # Parameters
  # ----------
  # fullglm : data.table 
  #  The processed CTMC data
  # covariates : vector
  #   Strings containing the various covariates to add
  # basedir : string
  #   The location of the covariate data
  #
  # Notes
  # -----
  # There will need to be conditions for each possible covariate

  if("drought" %in% covariates){

    # Match the drought data to the dates in fullglm
    files = Sys.glob(file.path(basedir, "drought", studynm, paste(studynm, "_drought*.csv", sep="")))
    drought = fread(files[1])[, list(DSCI, ValidStart, ValidEnd)]
    drought$ValidStart = as.numeric(as.POSIXct(drought$ValidStart, 
                            origin = '1970-01-01', tz = 'GMT'))
    drought$ValidEnd = as.numeric(as.POSIXct(drought$ValidEnd, 
                            origin = '1970-01-01', tz = 'GMT'))


    get_date = function(t){
      dt = date(as.POSIXct(t, origin = '1970-01-01', tz = 'GMT'))
      return(as.numeric(as.POSIXct(dt, origin = '1970-01-01', tz = 'GMT')))
    }

    fullglm[, c("ValidStart", "ValidEnd"):=list(get_date(t), get_date(t))]
    setkey(drought, ValidStart, ValidEnd)
    setkey(fullglm, ValidStart, ValidEnd)
    fullglm = foverlaps(fullglm, drought, type="within")

    # Remove the Valid* columns
    fullglm = fullglm[, -names(fullglm)[names(fullglm) %like% "Valid"], with=F]

  }

  return(fullglm)
  
}

get_regex_columns = function(lgvars, croptypes, landcovertypes, ext){
  # Get the regex representation of location or gradient columns
  #
  # Parameters 
  # ----------
  # lgvars : vector
  #   Strings of location of gradient variables
  # croptypes : vector
  #   For the "croplayer" variable, strings of the specific croptypes 
  #   that are present
  # landcovertypes : vector
  #   For the "landcover" variable, strings of the specific landcovertypes that 
  #   are present
  # ext : string
  #   Either "grad" or "loc" depending on whether the variables are gradient or 
  #   location variables.
  #
  # Returns
  # -------
  # : vector of strings
  #   The regex representations of the location of gradient variables

  regexnms = list()
  for(lvar in lgvars){

    if(lvar == "croplayer"){
      for(ctype in croptypes){
        regexnms[[ctype]] = paste(ctype, ".*", ext, sep="")
      }
    } else if(lvar == "landcover"){
      for(ltype in landcovertypes){
        regexnms[[ltype]] = paste(ltype, ".*", ext, sep="")
      }
    } else{
      regexnms[[lvar]] = paste(lvar, ".*", ext, sep="")
    }
  }
  return(regexnms)
}

get_timevar_status = function(allregexcols, croptypes, monthlyvars){
  # Get the time-varying status of a column name for each column
  # Returns vector with either "y" or "my"
  #
  # Parameters
  # ----------
  # allregexcols : vector
  #   Vector of regex columns names with ".*" notation.
  #   (i.e. c("temperature.*loc", "temperature.*grad")
  # croptypes : vector
  #   Vector of different croptypes as string
  # monthlyvars : vector
  #   Vector of covariate names that vary by month
  #
  # Returns
  # -------
  # : vector of "my" (for monthly variable) or "y" for yearly or "n" for neither

  timevar_vect = array(NA, length(allregexcols))
  for(i in 1:length(allregexcols)){

    col = allregexcols[i]
    vartype = strsplit(col, ".*", fixed=T)[[1]][1]
    if(vartype %in% croptypes){
      timevar_vect[i] = "y"
    } else if(vartype %in% monthlyvars){
      timevar_vect[i] = "my"
    } else {
      timevar_vect[i] = "n"
    }
  }
  return(timevar_vect)
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
  # : vector of distances in meters
  
  longs = longlat$longitude
  lats = longlat$latitude
  
  # Match distance vectors
  v1 = cbind(longs[-length(longs)], lats[-length(lats)])
  v2 = cbind(longs[-1], lats[-1])
  
  tuples = lapply(1:nrow(v1), function(i) list(v1[i, ], v2[i, ]))
  distvect = lapply(tuples, function(x) distm(x=x[[1]], y=x[[2]], 
                                                    fun=distHaversine))
  # 
  return(do.call(rbind, distvect))
}


runs = function(x, ctime=60, clength=10, ctimemin=0){
  # Given a vector of time differences X
  # determine if the vector has a run of 
  # clength or greater of diffs < ctime (in minutes). 
  #
  # Parameters
  # ----------
  # x : vector of datetimes
  # ctime : Diffs must be less than or equal to this time
  # clength : Length of a run
  #
  # Returns
  # -------
  # : bool
  #     TRUE or FALSE depending if if the vector meets the criteria
  
  # Notes
  # -----
  # This will help identify potential trajectories to use for the movement study
  # Returns TRUE or FALSE
  
  deltat = diff(x)
  units(deltat) = "mins"
  inds = (deltat <= ctime) & (deltat > ctimemin) # String of booleans
  
  # Count runs
  run_vect = rle(inds)
  res = any((run_vect$values == TRUE) & (run_vect$lengths >= clength))
  return(res)
  
}

split_runs = function(x, ctime=60, clength=10, ctimemin=0){
  # Given a vector of times x
  # find the indexes the runs that have a time difference <= ctime and 
  # a length greater than or equal to clength.
  #
  # Parameters
  # ----------
  # x : vector of datetimes
  # ctime : Diffs must be less than this time
  # clength : Length of a run
  #
  # Returns
  # -------
  # : a list of indexes in x that fulfill the criteria

  deltat = diff(x)
  units(deltat) = "mins"
  inds = (deltat <= ctime) & (deltat > ctimemin)

  falseinds = which(inds == FALSE)

  # All data is usable
  if(length(falseinds) == 0)
    return(list(1:length(x)))

  # No data is usable
  if(all(inds == FALSE))
    return(list(numeric()))

  # Extract various runs
  start = 1
  end = length(inds)

  runinds = list()

  for(i in 1:length(falseinds)){
    
    if(i == 1 & length(falseinds) != 1){
      
      if(falseinds[i] == 1)
        run = numeric()
      else
        run = 1:(falseinds[i])
      
      
    } else if(i > 1 & i < length(falseinds)){
      
      run = (falseinds[i - 1] + 1):(falseinds[i])
      
    } else if(i == length(falseinds) & length(falseinds) != 1){
      
      run = (falseinds[i - 1] + 1):(falseinds[i])
      
      if(falseinds[i] != end)
        runfinal = (falseinds[i] + 1):(end + 1)
      else
        runfinal = numeric()
      
    } else if(i == 1 & length(falseinds) == 1){

      run = 2:end

    }
    
    runinds[[i]] = run
    if(i == length(falseinds) & i != 1)
      runinds[[i + 1]] = runfinal
  }


  # Which runinds are of the right length?
  sinds = sapply(runinds, function(x) length(x) >= clength)
  return(runinds[sinds])
  
}

process_glmdata = function(glmdata, regexcol, colnm, format){
  # Sync-columns of glmdata across pigs
  #
  # Parmeters
  # ---------
  # glmdata : list
  #   List of data.tables with different pigIDs
  # regexcol : vector
  #  Regex expression for column names to extract and condense
  # colnm : vector
  #   Same length as regexcol.  The desired name of the regex columns after condensing.
  # format : vector
  #   Same length as regexcol. Either "my" if the variable varies monthly, or "y"
  #   if the variable varies yearly.
  #
  # Returns
  # -------
  # glmdata : list
  #   Returns glmdata with updated columns

  for(i in 1:length(glmdata)){

    gd = as.data.table(glmdata[[i]])

    # Format temperature

    # Format various time-dependent measures
    gd$datetime = as.POSIXct(gd$t, origin = '1970-01-01', tz = 'GMT')
    gd$hour = hour(gd$datetime)
    gd$month = month(gd$datetime)
    gd$year = year(gd$datetime)
    gd$monthyear = paste(gd$month, gd$year, sep="_")

    for(j in 1:length(regexcol)){

      if(format[j] == "my")
        timeref = gd$monthyear
      else
        timeref = gd$year

      gd[, colnm[j]] = condense_cols(gd, regexcol[j], timeref, timetype=format[j])

    }

    glmdata[[i]] = gd

  }

  return(glmdata)

}

condense_cols = function(gd, likename, timeref, timetype="my"){
  # Condense time-varying columns into one column
  #
  # Parameters
  # ----------
  # gd : data.table
  # likename: str, regex expression to match columns
  # timeref : vector with time either years or month_year values
  # timetype : str, "my" for month_year and "y" for year, "n" for not time varying


  temp_cols = names(gd)[names(gd) %like% likename]
  prefix = strsplit(likename, ".*", fixed=TRUE)[[1]][1]

  if(timetype != "n"){

    if(timetype == "my")
      temp_my = lapply(strsplit(sapply(strsplit(temp_cols, prefix), function(x) x[2]), "_"), function(x) paste(x[2], x[3], sep="_"))
    else if(timetype == "y")
      temp_my = lapply(strsplit(sapply(strsplit(temp_cols, prefix), function(x) x[2]), "_"), function(x) as.numeric(x[2]))

    boolmat = do.call(cbind, lapply(temp_my, function(x) x == timeref))
    temp_val = rowSums(boolmat * gd[, temp_cols, with=FALSE]) 

  } else{

    temp_val = gd[, temp_cols, with=FALSE]
  }
 
  return(temp_val)

}


get_distance_to_gradient = function(rast, decay=1){
  # Given a 0 or 1 raster, compute the weighted vector (e.g. gradient) for each 
  # cell that points towards the 1 values.  The 1 values might be crops, water,
  # etc.
  #
  # Parameters
  # ----------
  # raster : RasterLayer object
  #   RasterLayer with 0 and 1 values
  # decay : float
  #   The role of distance decay on influencing the the "pull" of a resource
  #
  # Returns
  # -------
  # : list
  #   xgrad : Raster with xgradient
  #   ygrad : Raster with ygradient
  #
  # Notes
  # -----
  # Need to make this much more efficient OR preprocess

  xras = yras = rast
  longlat_ctype = rasterToPoints(rast, function(x) x == 1)[, c("x", "y"), drop=FALSE]

  # Check if the calculations are going to be outrageous (i.e. > base)
  numcalc = ncell(rast) * nrow(longlat_ctype)
  base = 1.5e7 # ~45 seconds with 1e7, ~16 seconds with 1e6, ~1 minute with 1.5e7
  if(numcalc > base){
    # TODO: What if ncells > base?
    # Draw a random subset of values at which to compute the gradient. 
    # The idea is that this random sample will capture the general influence of 
    # the gradient. Otherwise, this will not be computationally feasible.
    numsamp = floor(base / ncell(rast))
    inds = sample(1:nrow(longlat_ctype), numsamp, replace=FALSE)
    longlat_ctype = longlat_ctype[inds, ] 
  }

  # If there are no crop types return 0 gradients
  if(nrow(longlat_ctype) == 0){

    values(xras) = 0
    values(yras) = 0

  } else{

    longlat_all = rasterToPoints(rast)[, c('x', 'y')]
    grad = t(apply(longlat_all, 1, wgradient, longlat_ctype, decay))
    values(xras) = grad[, 1]
    values(yras) = grad[, 2]
  }

  return(list(xgrad=xras, ygrad=yras))

}

scalar1 = function(x) {x / sqrt(sum(x^2))}

wgradient = function(pt, longlat_ctype, decay){
  # Compute the weighted sum of vectors based on distances.

  dists = as.vector(spDists(longlat_ctype, matrix(pt, nrow=1, ncol=2), longlat=TRUE))
  norm_vects =  t(apply(t(t(longlat_ctype) - pt), 1, scalar1))

  ind = dists != 0 # Drop any 0 distances (i.e. the cell itself)

  if(all(!ind))
    sum_grad = c(0, 0) # If there are no croptypes return 0 gradient.
  else{
    weighted_vects = norm_vects[ind, , drop=FALSE] * exp(-decay*dists[ind]^2) # Squared distance decay?
    sum_grad = colSums(weighted_vects) #scalar1
  }
  return(sum_grad)

}

## Analysis functions

build_design_matrix2 = function(dat, stdcols, nonstdcols, interactions=NULL){
  # Build a design matrix given a data set and columns to standardize and not
  # standardize
  #
  # Parameters 
  # ----------
  # dat : data.table
  #   The dataset
  # stdcols : vector of column names as strings
  #   Columns that should be standardized
  # nonstdcols : vector of column names as strings
  #   Columns that shouldn't be standardized
  # interactions : vector of strings
  #   Interaction strings in the form "crw:crop_loc".  Each columns must also
  #   be in stdcols or nonstdcols
  #
  # Returns
  # -------
  # : matrix
  #   The design matrix

  Xscale = dat[, lapply(.SD, function(x) scale2(x)), .SDcols = stdcols]

  Xfull = as.matrix(cbind(dat[, nonstdcols, with=F], Xscale))

  # Build interactions
  if(!is.null(interactions)){

    for(int in interactions){
      cols = strsplit(int, ":")[[1]]
      Xfull = cbind(Xfull[, cols[1]] * Xfull[, cols[2]], Xfull)
      colnames(Xfull)[1] = int
    }
  }
  return(Xfull)

}

build_daily_design_matrix = function(dat, stdcols, nonstdcols, splinecols,
                                          df_hour=0, interactions=NULL){
  # Build a design matrix for daily factor effect
  #
  # Parameters 
  # ----------
  # dat : data.table
  #   The dataset. Must have a column named hourofday and dayperiod if df_hour=0
  # stdcols : vector of column names as strings
  #   Columns that should be standardized
  # nonstdcols : vector of column names as strings
  #   Columns that shouldn't be standardized
  # splinecols : vector of column names as string
  #   Must be a subset of stdcols and contains the columns that will be
  #   converted into daily basis functions.
  # df_hour : int
  #  Number of basis vectors in cyclic basis function (df_hour - 1). 
  # If df_hour = 0, discrete time groupings are used.
  #
  # Returns
  # -------
  # : matrix
  #   The design matrix with seasonal factor.
  #
  # Notes
  # -----
  # Defaults to cyclic-cubic splines

  # Create hourly dummy variable

  if(df_hour == 0){

    # Build discrete day covariates
    timeofday = c("solarNoon-dusk", "dusk-nadir", "nadir-dawn", "dawn-solarNoon")

    hours_dummies = list()

    for(tod in timeofday){
      hours_dummies[[tod]] = as.numeric(dat$dayperiod == tod)
    }

    hourmat = do.call(cbind, hours_dummies)
    Xtod = hourmat
    todnames = timeofday

  } else{

    # Build continuous hour covariates
    splinehour = s(hourofday, bs="cc", k=df_hour) 
    hourmat = smooth.construct2(splinehour, dat, NULL)$X
    Xtod = hourmat
    todnames = 1:(df_hour - 1)
  }

  # Build non-spline matrix
  Xfull = build_design_matrix2(dat, stdcols, nonstdcols, 
                                      interactions=interactions)

  # Add daily effects onto columns
  for(colname in splinecols){
    hour_cols = Xfull[, colname] * hourmat
    Xtod = cbind(Xtod, hour_cols)
  }

  # Remove splinecols from Xfull
  Xfull_red = Xfull[, -which(colnames(Xfull) %in% splinecols)]
  Xtod_full = cbind(Xfull_red, Xtod)

  # Set column names
  colnames(Xtod_full) = c(colnames(Xfull_red), 
                          paste0("base_hour_", todnames),
                          do.call(c, lapply(splinecols, 
                                        function(x) paste0(x, "_", todnames))))
  return(Xtod_full)
}


build_seasonal_design_matrix = function(dat, stdcols, nonstdcols, 
                                                splinecols, seasonalcols,
                                                df_hour=0, interactions=NULL){
  # Build a design matrix with seasonal effects (winter, spring, summer, fall)
  #
  # Parameters
  # ----------
  # dat : data.table
  #   The dataset. Must have a column named hourofday and monthofyear
  # stdcols : vector of column names as strings
  #   Columns that should be standardized
  # nonstdcols : vector of column names as strings
  #   Columns that shouldn't be standardized
  # splinecols : vector of column names as strings
  #   Must be a subset of stdcols and contains the columns that will be
  #   converted into daily basis functions.
  # seasonalcols : vector of columns names as strings
  #   The columns that will interact with season
  # df_hours : int
  #   The there will be df_hours - 1 basis expansions (i.e. columns) for each
  #   spline. If df_hour == 0, use morning, midday, evening delineation.
  #
  # Returns
  # -------
  # : A design matrix with season effects.

  Xtod_full = build_daily_design_matrix(dat, stdcols, nonstdcols, splinecols, 
                                      df_hour=df_hour, interactions=interactions)

  # Create seasonal dummy variables
  seasons = list("summer" = c(6, 7, 8), "fall" = c(9, 10, 11), 
               "winter"= c(12, 1, 2), "spring"=c(3, 4, 5))
  season_dummies = list()

  for(season in names(seasons)){
    smonths = seasons[[season]]
    season_dummies[[season]] = as.numeric(dat$monthofyear %in% smonths)
  }

  seasonmat = do.call(cbind, season_dummies)
  Xseason = cbind(Xtod_full, seasonmat)

  if(df_hour == 0)
    tod = c("solarNoon-dusk", "dusk-nadir", "nadir-dawn", "dawn-solarNoon")
  else
    tod = 1:(df_hour - 1)

  for(scol in c(seasonalcols, "base_hour")) {

    for(season in names(seasons)){
      if((scol %in% splinecols) | (scol == "base_hour")) {

        varnames = paste0(scol, "_", tod)
        intercols = Xtod_full[, varnames] * seasonmat[, season]
        colnames(intercols) = paste0(varnames, ":", season)
        Xseason = cbind(intercols, Xseason)

      } else{

        intercols = Xtod_full[, scol] * seasonmat[, season]
        Xseason = cbind(intercols, Xseason)
        colnames(Xseason)[1] = paste0(scol, ":", season)

      }
    } # End season loop

    # Remove redundant columns
    if((scol %in% splinecols) | (scol == "base_hour"))
      Xseason = Xseason[, -which(colnames(Xseason) %in% varnames)]
    else
      Xseason = Xseason[, -which(colnames(Xseason) %in% scol)]


  } # End seasonal column loop

  return(Xseason)

}

build_tempprecip_design_matrix = function(dat, stdcols, nonstdcols, 
                                                splinecols, seasonalcols, 
                                                df_hour=5){
  # Build a design matrix with seasonal effects given by temperature and 
  # precipitation interactions
  #
  # Parameters
  # ----------
  # dat : data.table
  #   The dataset. Must have a column named hourofday and monthofyear
  # stdcols : vector of column names as strings
  #   Columns that should be standardized
  # nonstdcols : vector of column names as strings
  #   Columns that shouldn't be standardized
  # splinecols : vector of column names as strings
  #   Must be a subset of stdcols and contains the columns that will be
  #   converted into daily basis functions.
  # seasonalcols : vector of columns names as strings
  #   The columns that will interact with season
  # df_hours : int
  #   The there will be df_hours - 1 basis expansions (i.e. columns) for each
  #   spline.
  #
  # Returns
  # -------
  # : A design matrix with temp and precip effects.

  Xgam_full = build_daily_design_matrix(dat, stdcols, nonstdcols, 
                                    splinecols)

  # Standardize temperature and precip
  tempprecip = dat[, lapply(.SD, function(x) scale2(x)), 
                              .SDcols=c("temperature_loc", "precipitation_loc")]

  tempprecip[, ("temperature_loc:precipitation_loc"):=temperature_loc*precipitation_loc]
  tpmat = as.matrix(tempprecip)
  Xseason = cbind(Xgam_full, tpmat)

  tpcols = c("temperature_loc", "precipitation_loc", 
                        "temperature_loc:precipitation_loc")


  for(scol in c(seasonalcols, "base_hour")) {

    for(tpcol in tpcols){

      if((scol %in% splinecols) | (scol == "base_hour")) {

        varnames = paste0(scol, "_", 1:(df_hour - 1))
        intercols = Xgam_full[, varnames] * tpmat[, tpcol]
        colnames(intercols) = paste0(varnames, ":", tpcol)
        Xseason = cbind(intercols, Xseason)

      } else{

        intercols = Xgam_full[, scol] * tpmat[, tpcol]
        Xseason = cbind(intercols, Xseason)
        colnames(Xseason)[1] = paste0(scol, ":", tpcol)

      }
    }

    # # Remove redundant columns
    # if((scol %in% splinecols) | (scol == "base_hour"))
    #   Xseason = Xseason[, -which(colnames(Xseason) %in% varnames)]
    # else
    #   Xseason = Xseason[, -which(colnames(Xseason) %in% scol)]


  } # End seasonal column loop

  return(Xseason)

}

scale2 = function(x){
  # Allows for vector with no variance to be returned at 0
  if(length(unique(x)) == 1)
    sx = scale(x, scale=FALSE)
  else
    sx = scale(x)
  
  return(sx)
}


remove_nas_from_crop = function(tdat, cropdistgrad){
  # If cropdistgrad columns is NA, set to 0.

  # Set fngrad to 0 if it is NA, otherwise keep it the same
  for(cdist in cropdistgrad){

    if(all(is.na(tdat[[cdist]])))
       fngrad = rep(0, nrow(tdat))
    else
      fngrad = tdat[[cdist]]
    
    tdat[[cdist]] = fngrad
  }

  return(tdat)

}


diffunits = function(x){
  # Helper function to get the differences of a datetime vector in the same units
  # x is a datetime vector
  dt = diff(x)
  units(dt) = "mins"
  return(dt)
}


get_timeofday_groupings = function(sundata){
  # Using 
  #
  # Parameters
  # ----------
  # sundata : data.table with columns date, lat, lon

  dropdata = sundata[!duplicated(sundata$date)]
  suntimes = as.data.table(getSunlightTimes(data=as.data.frame(dropdata), tz="GMT"))
  suntimes = suntimes[, list(date, dusk, nadir, dawn, solarNoon)]

  suntimes$date = as.Date(suntimes$date)

  # Just merge on date as all pigs are in a similar geographical area.
  # mergeddat = merge(sundata, suntimes, by=c("date"))

  return(suntimes)

}


get_speed = function(lon, lat, datetime){
  # Get the average speed (km / h) for every transition. 
  #
  # Parameters
  # ----------
  # lon : longitude
  # lat : latitude
  # datetime : POISXct
  #
  # Returns
  # -------
  # : vector
  #   Average speed for each transition


  mindiffs = as.numeric(diffunits(datetime)) / 60 # to hours
  dists = as.vector(get_distance_vect(as.data.table(data.frame(longitude=lon, 
                        latitude=lat)))) / 1000 # to kilometers

  # Including NA ensures the same length
  return(c(NA, dists / mindiffs))
}

trim_speed = function(data, maxspeed=40){
  # Remove observations that are too fast

  speed = data[order(datetime, pigID), list(speed=get_speed(longitude, latitude, datetime), 
                    dists=c(NA, as.vector(get_distance_vect(
                                            as.data.table(
                                              data.frame(longitude=longitude, 
                                                         latitude=latitude)
                                                          ))) / 1000)), by=pigID]
  fastpig = (speed$speed > maxspeed)
  fastpig[is.na(fastpig)] = FALSE

  return(data[!fastpig])
}











