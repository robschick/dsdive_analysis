library(dplyr, lib.loc = c('singularity/libs', .libPaths()))

# load data
dive.files = dir(file.path('data', 'tag_endpoints'), pattern = 'ZcTag', 
                 full.names = TRUE)
dive.data = lapply(dive.files, function(f) readRDS(f))

# template bins
depth.bins = read.csv(file.path('data', 'imputed_bins', 'template.csv'))

# sex information for tag
tag.sex = read.csv(file.path('data', 'raw', 'tag_sex.csv'))

# initialize flattened structures
nim_pkg = list(
  data = list(
    depths = NULL,
    times = NULL
  ),
  consts = list(
    endpoint_priors = NULL,
    dive_relations = NULL,
    tag_covariates = NULL
  ),
  inits = list()
)

# merge tag records
for(i in 1:length(dive.data)) {
  
  attach(dive.data[[i]])
  
  # extract sex of tagged whale
  nim_pkg$consts$tag_covariates = c(
    nim_pkg$consts$tag_covariates,
    as.numeric(tag.sex %>% filter(deployid == name) %>% select(sex))
  )
  
  # id's of dives to keep from record
  dive.ids = intersect(which(dive.flags),
                       dive.ranges$dive.id[dive.ranges$type == 'Deep'])
  dive.ids = intersect(dive.ids, 
                       which(dive.ranges$end.ind - dive.ranges$start.ind + 1 > 
                               2))
  
  # id's of endpoints to keep from record
  endpoint.ids = endpoint.inds %>% 
    filter(dive_end %in% dive.ids | dive_start %in% dive.ids) %>% 
    select(endpoint.id) %>% 
    unlist()
  
  # initialize simple storage and lookup for merged id's
  endpoint.inds$merged.id = NA
  
  # merge endpoints
  for(j in 1:nrow(endpoint.inds)) {
    # only process valid endpoints
    if(endpoint.inds$endpoint.id[j] %in% endpoint.ids) {
      
      # construct id for endpoint in merged dataset
      flat_endpoint_id = ifelse(is.null(nim_pkg$consts$endpoint_priors), 0,
                                nrow(nim_pkg$consts$endpoint_priors)) + 1
      
      # store merged endpoint id
      endpoint.inds$merged.id[j] = flat_endpoint_id
      
      # merge endpoint
      nim_pkg$consts$endpoint_priors = rbind(
        nim_pkg$consts$endpoint_priors,
        c(t_lwr = endpoint.inds$t_lwr[j], 
          t_upr = endpoint.inds$t_upr[j])
      )
    }
  }
  
  # merge each dive
  for(j in 1:nrow(dive.ranges)) {
    # only process valid dives
    if(dive.ranges$dive.id[j] %in% dive.ids) {
      
      # construct id for dive in merged dataset
      flat_dive_id = ifelse(is.null(nim_pkg$consts$dive_relations), 0, 
                            nrow(nim_pkg$consts$dive_relations)) + 1
      
      # get index at which first depth will be stored
      flat_depth_ind = length(nim_pkg$data$depths) + 1
      
      # merge observed depths
      nim_pkg$data$depths = c(
        nim_pkg$data$depths,
        depths[dive.ranges$start.ind[j]:dive.ranges$end.ind[j]]
      )
      
      # merge observation times
      nim_pkg$data$times = c(
        nim_pkg$data$times,
        times[dive.ranges$start.ind[j]:dive.ranges$end.ind[j]]
      )
      
      # merge dive record
      nim_pkg$consts$dive_relations = rbind(
        nim_pkg$consts$dive_relations,
        c(T0_endpoint = endpoint.inds$merged.id[dive.ranges$T0_endpoint[j]], 
          T3_endpoint = endpoint.inds$merged.id[dive.ranges$T3_endpoint[j]],
          tag = i,
          depth_first = flat_depth_ind,
          depth_last = length(nim_pkg$data$depths))
      )
    }
  }
  
  detach(dive.data[[i]])
  rm(endpoint.inds)
  
}

rm(dive.files, dive.ids, endpoint.ids, flat_depth_ind, flat_dive_id, i, j,
   flat_endpoint_id)

# locations where there are consecutive deep dives
length(which(
  nim_pkg$consts$dive_relations[-1,'T0_endpoint'] == 
  nim_pkg$consts$dive_relations[1:(nrow(nim_pkg$consts$dive_relations)-1),
                                'T3_endpoint']
))

# extract sizes
nim_pkg$consts$N_tags = length(dive.data)
nim_pkg$consts$N_dives = nrow(nim_pkg$consts$dive_relations)
nim_pkg$consts$N_endpoints = nrow(nim_pkg$consts$endpoint_priors)
nim_pkg$consts$N_bins = nrow(depth.bins)

# export bin widths
nim_pkg$consts$widths = depth.bins$halfwidth * 2

# save output
saveRDS(nim_pkg, 
        file = file.path('data', 'tag_endpoints', 'flattened_endpoints.rds'))
