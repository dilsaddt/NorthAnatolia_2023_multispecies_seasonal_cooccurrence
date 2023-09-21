# NorthAnatolia_2023_multispecies_seasonal_cooccurrence
Multispecies occupancy models to investigate seasonal co-occurrence of predator-prey pairs and changes withing these co-occurrences.


In total seven predator-prey pairs were analyzed: 
- gray wolf (_Canis lupus_) and roe deer (_Capreolus capreolus_)
- gray wolf (_Canis lupus_) and red deer (_Cervus elaphus_)
- gray wolf (_Canis lupus_) and wild boar (_Sus scrofa_)
- gray wolf (_Canis lupus_) and European hare (_Lepus europaeus_)

- Eurasian lynx (_Lynx lynx_) and roe deer (_Capreolus capreolus_)
- Eurasian lynx (_Lynx lynx_) and wild boar (_Sus scrofa_)
- Eurasian lynx (_Lynx lynx_) and European hare (_Lepus europaeus_)

Here, as an example only the analysis of gray wolf (_Canis lupus_) and roe deer (_Capreolus capreolus_) is given. 

## Folders

### Static_framework folder

Static multispecies models to check the seasonal co-occurrence.

Seasons (summer and winter) were analyzed separately to account for the seasonal effect with static framework.

cctocl_static_summer_072023.txt = JAGS model text file for the named model type

cctocl_static_summer_072023_github.R = R file to run the JAGS model

cctocl_static_summer_diagnostics_072023_github.R = R file for model diagnostics

same for winter ones.

cctocl_static_prediction_072023_github.R = Predictions based on model outputs, summer and winter together.

#### Static_data folder

cctocl_static_summer_data_github.RData = prepared JAGS data with static detection history of both species in summer season and model variables

cctocl_static_winter_data_github.RData = prepared JAGS data with static detection history of both species in winter season and model variables

### Dynamic_framework folder

Dynamic multispecies models to check the changes in seasonal co-occurrence.

Three different models were applied with season, rural human population density, and elevation covariates.

These models each had their own folders as:
- additive (add) = season + rural human population density + elevaiton
- interaction 1 (int1) = season * rural human population density
- interaction 2 (int2) season * elevation

cctocl_dynamic_<model_type>_072023.txt = JAGS model text file for the named model type

cctocl_dynamic_<model_type>_072023_github_.R = R file to run the named JAGS model

cctocl_dynamic_<model_type>_072023_gof.txt = JAGS model text file for the named model type's goodness-of-fit test

cctocl_dynamic_<model_type>_GOF_072023_github_.R = R file to run the named model's goodness-of-fit test

cctocl_dynamic_<model_type>_diagnostics_github_072023.R = R file for named model's diagnostics

cctocl_dynamic_<model_type>_prediction_github.R = R file for predictions from named model's outputs


cctocl_static_summer_072023.txt = JAGS model text file for the named model type

cctocl_static_summer_072023_github.R = R file to run the JAGS model

cctocl_static_summer_diagnostics_072023_github.R = R file for model diagnostics

same for winter ones.

cctocl_static_prediction_072023_github.R = Predictions based on model outputs, summer and winter together.

#### Dynamic_data folder

cctocl_dynamic_multisp_data_github.RData = prepared JAGS data with dynamic detection history of both species and model variables

