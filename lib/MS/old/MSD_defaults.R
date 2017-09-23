#Default settings
#Hormone cycle simulation
#By Kevin McKee

defaultParams<-function(){
  # Parameters --------------------------------------------------------------
  N<<- 7*3*3 #Choose via ModelUnivariate_Sim_SeriesLength.r?
  etaX1<<- -0.5
  zetaX1<<- -0.3
  etaX2<<- -0.05
  zetaX2<<- -0.03
  gammaX1<<-  0
  gammaX2<<-  0
  tSNR<<- 8

  # Embedding Settings ------------------------------------------------------
  embedD1<<- 3
  embedD2<<- 4
  theTau<<- 1
  deltaT<<- 1
  
  #Simulation settings: -----------------------------------------------------
  pEvent<<-0.5
  eventScale<<-1
  # eventScaleRatio<<- 1/4
  eventNormal<<-T
  stationary<<-T
  
  
  initCond.X<<- runif(1, 1.5, 4)*sample(c(-1,1),1)
  initCond.dX<<-0
  eventTypeX1<<-"level"
  # eventScaleX1<<-1
  pEventX1<<- 0#0.05
  
  initCond.Y<<- runif(1, 0.5, 1)*sample(c(-1,1),1)
  initCond.dY<<-0
  eventTypeX2<<-"slope"
  # eventScaleX2<<-0.25
  pEventX2<<- 0#0.02
  

  
  # windowN<<- 400 #default, not applicable to series length test
  source("lib/MS/MSD_matrices.R")
  source("lib/MS/MS_model.R")
  source("lib/MS/MS_constraints.R")
  
}
defaultParams()






