#### LOAD LIBRARIES ####

# Step 0: Install packages using groundhog
install.packages("groundhog")
library("groundhog")

pkg_list <- c(
  "effects", "carData", "DHARMa", "lubridate", "forcats", "stringr",
  "dplyr", "purrr", "readr", "tidyr", "tibble", "ggplot2", "tidyverse",
  "nlme", "piecewiseSEM", "lme4", "lmerTest", "performance")

groundhog.library(pkg_list, "2022-01-01")

# Step 2: Load libraries
library(piecewiseSEM)
library(nlme)
library(tidyverse)
library(DHARMa)
library(effects)
library(lme4)
library(lmerTest)
library(performance)

#
#### LOAD DATA ####

#Stability_data <- read.csv(file.choose())
Stability_data <- read.csv("Data/Stability_metrics.csv")

# LOAD FUNCTIONS
range0_1 <- function(x,na.rm){(x-min(x, na.rm = T)+0.01)/(max(x, na.rm = T)-min(x, na.rm = T)+0.01)}
z_trans <- function(myVar, na.rm){(myVar - mean(myVar, na.rm = T)) / sd(myVar, na.rm = T)}
'%ni%' <- Negate('%in%')

#
#### LINEAR MODELS - SPECIES GROUPS ####

Barnacles1_model<-lmerTest::lmer(Barnacles1_mean ~  
                         Emersion.Rate + I(Emersion.Rate^2) +
                         Emersion.Rate * Tile.Type +
                         #(1|Transect/Station) # Transect explained 0 variance, and was not allowing the computing or r2 values
                         (1|Station)
                       , 
                       Stability_data
); plot(simulateResiduals(Barnacles1_model)); summary(Barnacles1_model); plot(allEffects(Barnacles1_model))
performance::r2(Barnacles1_model)

Barnacles2_model<-lmerTest::lmer(sqrt(Barnacles2_mean) ~  
                         Emersion.Rate + I(Emersion.Rate^2) +
                         Emersion.Rate * Tile.Type +
                         (1|Transect/Station)
                       , 
                       Stability_data
); plot(simulateResiduals(Barnacles2_model)); summary(Barnacles2_model); plot(allEffects(Barnacles2_model))
performance::r2(Barnacles2_model)

Ac_invertebrates_model<-lmerTest::lmer(Ac_invertebrates_mean ~  
                         #Emersion.Rate + I(Emersion.Rate^2) + # not significant, model did not converge
                         Emersion.Rate * Tile.Type +
                         (1|Transect/Station)
                       , 
                       Stability_data
); plot(simulateResiduals(Ac_invertebrates_model)); summary(Ac_invertebrates_model); plot(allEffects(Ac_invertebrates_model))

Eph_macroalgae_model<-lmerTest::lmer(sqrt(Eph_macroalgae_mean) ~  
                         Emersion.Rate + I(Emersion.Rate^2) +
                         Emersion.Rate * Tile.Type +
                           #(1|Transect/Station) # Transect explained 0 variance, and was not allowing the computing or r2 values
                           (1|Transect:Station)
                       , 
                       Stability_data
); plot(simulateResiduals(Eph_macroalgae_model)); summary(Eph_macroalgae_model); plot(allEffects(Eph_macroalgae_model))
performance::r2(Eph_macroalgae_model)

Consumer_model<-lmerTest::lmer(sqrt(Consumer_mean) ~  
                         Emersion.Rate + I(Emersion.Rate^2) +
                         Emersion.Rate * Tile.Type +
                           #(1|Transect/Station) # Transect explained 0 variance, and was not allowing the computing or r2 values
                           (1|Station)
                       , 
                       Stability_data
); plot(simulateResiduals(Consumer_model)); summary(Consumer_model); plot(allEffects(Consumer_model))
performance::r2(Consumer_model)

S_model<-lmerTest::lmer(log(S_mean) ~  
                         Emersion.Rate + I(Emersion.Rate^2) +
                         Emersion.Rate * Tile.Type +
                          #(1|Transect/Station) # Transect explained 0 variance, and was not allowing the computing or r2 values
                          (1|Transect)
                       , 
                       Stability_data
); plot(simulateResiduals(S_model)); summary(S_model); plot(allEffects(S_model))
performance::r2(S_model)

Pielou_model<-lm((Pielou_mean) ~  
                         #Emersion.Rate + I(Emersion.Rate^2) + # not significant, model did not converge
                         Emersion.Rate * Tile.Type #+
                         #(1|Transect/Station) # when including random effects, results do not change qualitatively. Since ggpredict can't compute confidence intervals for this model when including random effects, they had to be removed.
                       , 
                       Stability_data
); plot(simulateResiduals(Pielou_model)); summary(Pielou_model); plot(allEffects(Pielou_model))
performance::r2(Pielou_model)

#
#### LINEAR MODELS - STABILITY DIMENSIONS ####

Pop_stab_model<-lmerTest::lmer(log(Pop_stab) ~  
                         Emersion.Rate + I(Emersion.Rate^2) +
                         Emersion.Rate * Tile.Type + # model did not converge with Emersion.Rate * I(Emersion.Rate^2) * Tile.Type
                         (1|Transect/Station)
                       , 
                       Stability_data
); plot(simulateResiduals(Pop_stab_model)); summary(Pop_stab_model); AIC(Pop_stab_model); plot(allEffects(Pop_stab_model))
performance::r2(Pop_stab_model)

Asynchrony_model<-lmerTest::lmer(Asynchrony_corr_time ~  
                         Emersion.Rate * Tile.Type +
                         (1|Transect/Station)
                       , 
                       Stability_data
); plot(simulateResiduals(Asynchrony_model)); summary(Asynchrony_model); AIC(Asynchrony_model); plot(allEffects(Asynchrony_model))
performance::r2(Asynchrony_model)

Stat_av_model<-lmerTest::lmer(Statistical_averaging ~  
                         Emersion.Rate * Tile.Type +
                           #(1|Transect/Station) # Transect explained 0 variance, and was not allowing the computing or r2 values
                           (1|Transect:Station)
                       , 
                       Stability_data
); plot(simulateResiduals(Stat_av_model)); summary(Stat_av_model); AIC(Stat_av_model); plot(allEffects(Stat_av_model))
performance::r2(Stat_av_model)

Comp_model<-lm((Compositional_similarity) ~  
                         Emersion.Rate * I(Emersion.Rate^2) +
                         Emersion.Rate * Tile.Type #+ # model did not converge with Emersion.Rate * I(Emersion.Rate^2) * Tile.Type
                         #(1|Station) # transect has 0 effect and station does not allow model to converge
                       , 
                       Stability_data
); plot(simulateResiduals(Comp_model)); summary(Comp_model); AIC(Comp_model); plot(allEffects(Comp_model))
performance::r2(Comp_model)

Cov_model<-lmerTest::lmer(log(Cover_stability_alpha_res) ~  
                         Emersion.Rate * I(Emersion.Rate^2) +
                         Emersion.Rate * Tile.Type +
                         (1|Transect/Station)
                       , 
                       Stability_data
); plot(simulateResiduals(Cov_model)); summary(Cov_model); AIC(Cov_model); plot(allEffects(Cov_model))
performance::r2(Cov_model)

#

#### STRUCTURAL EQUATION MODELLING ####

#Trandform metrics
Stability_data_SEM <- Stability_data %>%
  mutate(
    Consumer_mean_log = log(Consumer_mean+0.1),
    Limpet_mean_log = log(Limpet_mean+0.1),
    Barnacles2_mean_sqrt = sqrt(Barnacles2_mean),
    Eph_macroalgae_mean_sqrt = sqrt(Eph_macroalgae_mean),
    S_mean_log = log(S_mean+1),
    Pop_stab_log = log(Pop_stab),
    Cover_stability_alpha_res_log = log(Cover_stability_alpha_res),
    Tile.Type2 = as.numeric(ifelse(Tile.Type=="R",1,0)),
    EmGroup = as.factor(EmGroup),
    Transect=as.factor(Transect),
    Station=as.factor(Station)
  ) %>%
  mutate(
    Emersion.Rate_z = range0_1(z_trans(Emersion.Rate, na.rm=TRUE), na.rm=TRUE),
    Emersion.Rate2_z = range0_1(z_trans(Emersion.Rate^2, na.rm=TRUE), na.rm=TRUE),
    Emersion.Rate_int_z = range0_1(z_trans((Emersion.Rate^2)*Emersion.Rate, na.rm=TRUE), na.rm=TRUE),
    Emersion.Rate_Tile_z = range0_1(z_trans(Emersion.Rate*Tile.Type2, na.rm=TRUE), na.rm=TRUE),
    Consumer_mean_log_z = range0_1(z_trans(Consumer_mean_log, na.rm=TRUE), na.rm=TRUE),
    Limpet_mean_log_z = range0_1(z_trans(Limpet_mean_log, na.rm=TRUE), na.rm=TRUE),
    Barnacles1_mean_z = range0_1(z_trans(Barnacles1_mean, na.rm=TRUE), na.rm=TRUE),
    Barnacles2_mean_sqrt_z = range0_1(z_trans(Barnacles2_mean_sqrt, na.rm=TRUE), na.rm=TRUE),
    Ac_invertebrates_mean_z = range0_1(z_trans(Ac_invertebrates_mean, na.rm=TRUE), na.rm=TRUE),
    Eph_macroalgae_mean_sqrt_z = range0_1(z_trans(Eph_macroalgae_mean_sqrt, na.rm=TRUE), na.rm=TRUE),
    S_mean_log_z = range0_1(z_trans(S_mean_log, na.rm=TRUE), na.rm=TRUE),
    Pielou_mean_z = range0_1(z_trans(Pielou_mean, na.rm=TRUE), na.rm=TRUE),
    Pop_stab_log_z = range0_1(z_trans(Pop_stab_log, na.rm=TRUE), na.rm=TRUE),
    Asynchrony_corr_time_z = range0_1(z_trans(Asynchrony_corr_time, na.rm=TRUE), na.rm=TRUE),
    Statistical_averaging_z = range0_1(z_trans(Statistical_averaging, na.rm=TRUE), na.rm=TRUE),
    Compositional_similarity_z = range0_1(z_trans(Compositional_similarity, na.rm=TRUE), na.rm=TRUE),
    Cover_stability_alpha_res_log_z = range0_1(z_trans(Cover_stability_alpha_res_log, na.rm=TRUE), na.rm=TRUE)
  ) %>% drop_na()

# Fit the SEM model
SEM_model_out<-psem(
    
    lme(Consumer_mean_log_z ~ Emersion.Rate_z + Tile.Type2,  random = ~1|Transect/Station, Stability_data_SEM),
    lme(Barnacles1_mean_z ~ Emersion.Rate_z + Tile.Type2 + Consumer_mean_log_z, random = ~1|Station, Stability_data_SEM), # removed Transect and Station as RE as it had 0 variance
    lme(Barnacles2_mean_sqrt_z ~ Emersion.Rate_z + Tile.Type2 + Consumer_mean_log_z,  random = ~1|Station, Stability_data_SEM), # removed Transect and Station as RE as it had 0 variance
    lme(Eph_macroalgae_mean_sqrt_z ~ Emersion.Rate_z + Tile.Type2 + Consumer_mean_log_z,  random = ~1|Station, Stability_data_SEM), # removed Transect and Station as RE as it had 0 variance
    
    lm(S_mean_log_z ~ Emersion.Rate_z + Tile.Type2 + Consumer_mean_log_z, Stability_data_SEM), # removed Transect and Station as RE as it had 0 variance
    lme(Pielou_mean_z ~ Emersion.Rate_z + Tile.Type2 + Consumer_mean_log_z,  random = ~1|Station, Stability_data_SEM), # removed Transect and Station as RE as it had 0 variance
    
    lm(Pop_stab_log_z ~ Emersion.Rate_z + Tile.Type2 + S_mean_log_z + Pielou_mean_z + Barnacles1_mean_z + Barnacles2_mean_sqrt_z + Consumer_mean_log_z + Eph_macroalgae_mean_sqrt_z, data=Stability_data_SEM),  # removed Transect and Station as RE as it had 0 variance
    lme(Asynchrony_corr_time_z ~ Emersion.Rate_z + Tile.Type2 + S_mean_log_z + Pielou_mean_z + Barnacles1_mean_z + Barnacles2_mean_sqrt_z + Consumer_mean_log_z + Eph_macroalgae_mean_sqrt_z, random = ~1|Transect/Station, Stability_data_SEM),
    lme(Statistical_averaging_z ~ Emersion.Rate_z + Tile.Type2 + S_mean_log_z + Pielou_mean_z + Barnacles1_mean_z + Barnacles2_mean_sqrt_z + Consumer_mean_log_z + Eph_macroalgae_mean_sqrt_z, random = ~1|Transect/Station, data=Stability_data_SEM),
    
    lm(Compositional_similarity_z ~ Emersion.Rate_z + Tile.Type2+ Asynchrony_corr_time_z + Pop_stab_log_z +  Statistical_averaging_z + S_mean_log_z + Pielou_mean_z + Barnacles1_mean_z + Barnacles2_mean_sqrt_z + Consumer_mean_log_z + Eph_macroalgae_mean_sqrt_z, Stability_data_SEM),   # removed Transect and Station as RE as it had 0 variance
    lme(Cover_stability_alpha_res_log_z ~ Emersion.Rate_z + Tile.Type2 + Asynchrony_corr_time_z + Pop_stab_log_z +  Statistical_averaging_z + Compositional_similarity_z,  random = ~1|Transect/Station, Stability_data_SEM),
    
    Consumer_mean_log_z %~~% Cover_stability_alpha_res_log_z,
    
    Barnacles1_mean_z %~~% Barnacles2_mean_sqrt_z,
    Barnacles1_mean_z %~~% Cover_stability_alpha_res_log_z,
    Barnacles1_mean_z %~~% Eph_macroalgae_mean_sqrt_z,
    Barnacles1_mean_z %~~% S_mean_log_z,
    Barnacles1_mean_z %~~% Pielou_mean_z,
    
    Barnacles2_mean_sqrt_z %~~% Cover_stability_alpha_res_log_z,
    Barnacles2_mean_sqrt_z %~~% S_mean_log_z,
    Barnacles2_mean_sqrt_z %~~% Eph_macroalgae_mean_sqrt_z,
    Barnacles2_mean_sqrt_z %~~% Pielou_mean_z,
    
    Eph_macroalgae_mean_sqrt_z %~~% S_mean_log_z,
    Eph_macroalgae_mean_sqrt_z %~~% Pielou_mean_z,
    
    S_mean_log_z %~~% Pielou_mean_z,
    S_mean_log_z %~~% Cover_stability_alpha_res_log_z,
    
    Pop_stab_log_z %~~% Compositional_similarity_z,
    Pop_stab_log_z %~~% Statistical_averaging_z,
    Pop_stab_log_z %~~% Asynchrony_corr_time_z,
    
    Asynchrony_corr_time_z %~~% Compositional_similarity_z,
    Asynchrony_corr_time_z %~~% Statistical_averaging_z,
    
    Compositional_similarity_z %~~% Statistical_averaging_z

)

# SEM Model summary
summary(SEM_model_out)

SEM_one_coefs<-summary(SEM_model_out)$coefficients %>%
  rename(Var.9=9) %>%
  filter(Var.9!='')

SEM_model_out.multigroup <- multigroup(SEM_model_out, group = "EmGroup")

model_H1.multigroup_modelInteractions <- SEM_model_out.multigroup[["anovaInts"]] %>%
  rename("Sig"=6) %>%
  filter(Sig!="") %>%
  ungroup() %>%
  dplyr::select(Response,Predictor,P.Value,Sig) %>%
  mutate(
    Predictor=gsub("\\:EmGroup|\\EmGroup:|_log|_mean","",Predictor),
    Response=gsub("_log|_mean","",Response)
  )

Get_multigroup_coefs <- function(input_dataset) {
  
  extract_dataset<-input_dataset
  
  output<-extract_dataset$group.coefs %>%
    Map(cbind, ., EmGroup=names(.)) %>% # add column indicating which Tile and date each dataframe correspondsto
    do.call(rbind,.)
  
  return(output)
}
SEM_multigroup_coefs<-Get_multigroup_coefs(SEM_model_out.multigroup) %>%
  mutate(Std.Estimate=as.numeric(Std.Estimate)) %>%
  dplyr::select(EmGroup, Response, Predictor, Std.Estimate,P.Value,Var.9)

# Export SEM results
write.csv(SEM_multigroup_coefs, 'SEM_multigroup_coefs.csv')
write.csv(model_H1.multigroup_modelInteractions, 'SEM_multigroup_interactions.csv')
write.csv(SEM_one_coefs, 'SEM_one_coefs.csv')

#