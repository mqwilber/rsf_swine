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

cleanTime = function(){
  # 


}