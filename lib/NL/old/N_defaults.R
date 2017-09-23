#Default settings
#Hormone cycle simulation
#By Kevin McKee

defaultParams<-function(){
  # Parameters --------------------------------------------------------------
  N<<- 45 #Choose via ModelUnivariate_Sim_SeriesLength.r?
  tSNR<<- 4
  nVars<<-3
  
  # Embedding Settings ------------------------------------------------------
  embedD<<- 4
  theTau<<- 1
  deltaT<<- 1
  
  #Simulation settings: -----------------------------------------------------
  pEventTotal<<-0.0
  pEvent<<- 0.0
  eventType<<-"level"
  eventScale<<-4
  eventNormal<<-F
  
  initCond<<-matrix(0,nVars, 2)
  initCond[,1]<<-rnorm(nVars)
  initCond[,2]<<-rnorm(nVars)
  
  stationary<<-T
  
  
  source("lib/NL/N_matrices.R")
  source("lib/NL/N_model.R")
  # source("N_matrices_RAM.R")
  # source("N_model_RAM.R")
  source("lib/NL/N_constraints.R")
  
  
  indDyn<<- matrix(c(   -.5,   -.3,
                        -.3,   -.15,
                        -.2,   -.2,
                        -.1,   -.1),
               4, 2, byrow=T)
  
  gammas<<-matrix(c(   1,   0.05,   -0.2,     -0.1,
                      -0.3,   1,     0,        0.2,
                      -.05,  -0,     1,       -0.2,
                      -.05,  -0.2,  -0.2,        1  )

   ,4,4, byrow=T)
  
  # gammas.S<<-matrix(c(   1,   0.05,   -0.2,     -0.1,
  #                      -0.3,   1,     0,        0.2,
  #                      -.05,  -0,     1,       -0.2,
  #                      -.05,  -0.2,  -0.2,        1  )
  #                 
  #                 ,4,4, byrow=T)
  # 
  
  indDyn<<-matrix(indDyn[1:nVars,], nVars, 2, byrow=T)
  gammas<<-matrix(gammas[1:nVars,1:nVars], nVars, nVars, byrow=T)




}
defaultParams()




