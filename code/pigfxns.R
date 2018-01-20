## Functions for analyzing feral swine movement data

library(ctmcmove)
library(geosphere)

get_distance_vect = function(longlat){
  # Calculate distances between sequential longlat points 
  #
  # Parameters
  # ----------
  # longlat : data.frame with columns longitude and latitud3e
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
  return(distvect)
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
  # Simulate the discrete time steps of animal movement based on the CTMC generator matrix R, 
  # and the the GLM model.
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
  
  step_direction[i] = 1 # Start moving north
  
  step_direction[i] = 0
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


runs = function(x, ctime=60, clength=10, ctimemin=0){
  # Given a vector of time differneces X
  # determine if the vector has a run of 
  # clength or greater of diffs < ctime (in minutes). 
  # Parameters
  # ----------
  # x : difference vector of times
  # ctime : Diffs must be less than this time
  # clength : Length of a run
  
  # Notes
  # -----
  # This will help identify potential trajectories to use for the movement study
  # Returns TRUE or FALSE
  
  deltat = diff(x)
  inds = (deltat < ctime) & (deltat > ctimemin) # String of booleans
  
  # Count runs
  run_vect = rle(inds)
  res = any((run_vect$values == TRUE) & (run_vect$lengths >= clength))
  return(res)
  
}

