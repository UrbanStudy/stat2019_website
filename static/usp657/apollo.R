
rm(list = ls())
setwd("~/qushen26/stat2019_website/static/usp657")
heating.w = read.table('datawide.asc')
names(heating.w) = c("idcase", "depvar", "ic1",  "ic2",  "ic3",  "ic4",  "ic5",  
                     "oc1", "oc2",  "oc3",  "oc4",  "oc5",  "income", "agehed", 
                     "rooms",  "ncoast", "scoast", "mountn", "valley")

library(apollo)
apollo_initialise()

apollo_control=list(
  modelName ="HW1" ,
  modelDescr="USP 657 HW1 Analyzing heating options data",
  indivID="idcase",
  panelData=FALSE
)

database = heating.w

### subset data if necessary
#database = subset(database, SP==1)

apollo_beta=c(asc_1gc  = 0,
              asc_2gr  = 0,
              asc_3ec  = 0,
              asc_4er  = 0,
              asc_5hp  = 0,
              b_ic = 0,
              b_oc = 0)

apollo_fixed = c("asc_1gc")

apollo_inputs = apollo_validateInputs()

apollo_probabilities=function(apollo_beta, apollo_inputs, 
                              functionality="estimate"){
  
  ### Attach inputs and detach after function exit
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))
  
  ### Create list of probabilities P
  P = list()
  
  ### List of utilities: these must use the same names as
  ### in mnl_settings, order is irrelevant.
  V = list()
  V[['gc']] = asc_1gc + b_ic*ic1 + b_oc*oc1
  V[['gr']] = asc_2gr + b_ic*ic2 + b_oc*oc2
  V[['ec']] = asc_3ec + b_ic*ic3 + b_oc*oc3
  V[['er']]=  asc_4er + b_ic*ic4 + b_oc*oc4
  V[['hp']]=  asc_5hp + b_ic*ic5 + b_oc*oc5
  
  ### Define settings for MNL model component
  mnl_settings = list(
    alternatives  = c(gc=1, gr=2, ec=3, er=4, hp=5), 
    avail         = list(gc=1, gr=1, ec=1, er=1, hp=1), 
    choiceVar     = depvar,
    V             = V
  )
  
  ### Compute probabilities using MNL model
  P[['model']] = apollo_mnl(mnl_settings, functionality)
  
  ### Take product across observation for same individual
  #P = apollo_panelProd(P, apollo_inputs, functionality)
  
  ### Prepare and return outputs of function
  P = apollo_prepareProb(P, apollo_inputs, functionality)
  return(P)
}

model = apollo_estimate(apollo_beta, apollo_fixed, 
                        apollo_probabilities, 
                        apollo_inputs)

apollo_modelOutput(model)

# ####################################################### #
#### 7. Postprocessing of results                      ####
# ####################################################### #

### Use the estimated model to make predictions
predictions_base = apollo_prediction(model, 
                                     apollo_probabilities, 
                                     apollo_inputs)

### Now imagine the cost for gc increase by 10% 
### and predict again
database$ic1 = 1.1 * database$ic1
apollo_inputs   = apollo_validateInputs()

predictions_new = apollo_prediction(model, 
                                    apollo_probabilities, 
                                    apollo_inputs)

### Compare predictions
change=(predictions_new-predictions_base)/predictions_base
### Not interested in chosen alternative now, 
### so drop last column
change=change[,-ncol(change)]
### Summary of changes (possible presence of NAs due to
### unavailable alternatives)
summary(change)
