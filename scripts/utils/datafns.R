# load data from paths
dives.load = function(path, dive_pattern, depth_pattern) {
  
  # identify dive files
  dives = dir(path = path, pattern = dive_pattern, full.names = TRUE)
  depths = dir(path = path, pattern = depth_pattern, full.names = TRUE)
  
  # storage for observed dives
  dives.obs = vector('list', length(dives))
  
  # load dives 
  for(i in 1:length(dives)) {
    
    # read raw data
    d = read.csv(file = dives[i], header = TRUE)
    db = read.csv(file = depths[i], header = TRUE)
    
    # center dive times
    d$times = d$times - d$times[1]
    
    # save
    dive = as.list(d)
    class(dive) = 'dsobs'
    dives.obs[[i]] = list(
      dive = dive,
      depth.bins = db
    )
    
  }
  
  dives.obs
}


# select dives for analysis
fit.ind.fn = function(dives.obs, duration_min, duration_max, holdout, 
                      seed=NULL, holdout_prop, bin_start_max, bin_end_max) {
  
  # pre-select dives by their durations, start and end bins
  all.inds = sapply(1:length(dives.obs), function(ind) {
    d = dives.obs[[ind]]$dive
    dur = max(d$times)
    ifelse(all(dur >= duration_min, 
               dur <= duration_max,
               d$depths[1] <= bin_start_max,
               d$depths[length(d$depths)] <= bin_end_max), ind, NA)
  })
  all.inds = all.inds[!is.na(all.inds)]
  
  # partition dives into training and test sets, if requested
  if(holdout) {
    if(!is.null(seed)) {
      set.seed(seed)
    }
    inds = sample(1:length(all.inds), size = length(all.inds) * holdout_prop)
  } else {
    inds = 1:length(all.inds)
  }
  
  # package result into testing and training sets
  list(
    fit = all.inds[inds],
    validate = all.inds[-inds]
  )
}