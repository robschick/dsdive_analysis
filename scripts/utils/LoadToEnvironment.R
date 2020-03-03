LoadToEnvironment <- function(RData, env = new.env()){
  # load an RData file to a(n) (new) environment; return the environment
  load(RData, env)
  return(env) 
}
