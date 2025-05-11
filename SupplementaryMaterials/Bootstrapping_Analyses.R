library(tidyverse)
library(tidyr)
library(lme4)
library(lmerTest)
library(DHARMa)
library(effects)
library(piecewiseSEM)
library(nlme)

# 1 - Load Datasets ------------

# LOAD FUNCTIONS
range0_1 <- function(x,na.rm){(x-min(x, na.rm = T)+0.01)/(max(x, na.rm = T)-min(x, na.rm = T)+0.01)}
z_trans <- function(myVar, na.rm){(myVar - mean(myVar, na.rm = T)) / sd(myVar, na.rm = T)}
'%ni%' <- Negate('%in%')

#Stability_data <- read.csv(file.choose())
Stability_data <- read.csv("Data/Stability_metrics.csv")

Stability_data_SEM <- Stability_data %>%
  mutate(
    Consumer_mean_log = log(Consumer_mean+0.1),
    Limpet_mean_log = log(Limpet_mean+0.1),
    Barnacles2_mean_sqrt = sqrt(Barnacles2_mean),
    Eph_macroalgae_mean_sqrt = sqrt(Eph_macroalgae_mean),
    S_mean_log = log(S_mean+1),
    Pop_stab_log = log(Pop_stab),
    Cover_stability_alpha_res_log = log(Cover_stability_alpha_res),
    Tile.Type2 = as.numeric(ifelse(Tile.Type=="R",1,0))
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

#
# 2 - Bootstrapping Linear Models --------------

set.seed(123)  # For reproducibility

# Obtain dataset for within station

  resample_data_within <- function(data) {
    
    resampled_data <- data %>%
      mutate(group_id = paste(Transect, Station, sep = "_")) %>%
      nest(data = -group_id) %>%
      sample_n(18) %>%           
      unnest(cols = data)
    
    return(resampled_data)
  }
  
  # Repeat the resampling process 1000 times and store results
  resampled_datasets_within <- replicate(1000, resample_data_within(Stability_data), simplify = FALSE)
  
  # Combine all resamples into one large dataframe with an identifier
  resampled_df_within <- bind_rows(resampled_datasets_within, .id = "resample_id")

# Obtain dataset for across stations

  resample_data_across <- function(data) {
    
    # Step 1: Shuffle data and select one row per (Transect, Station) while alternating S and R
    resampled_data <- data %>%
      group_by(Transect, Station) %>%
      mutate(Random_Number = runif(1)) %>%  # Assign one random number per group
      ungroup() %>%
      arrange(Random_Number) %>%  # Sort rows from low to high
      # Randomly assign "A" or "B" to half of the unique (Transect, Station) groups
      mutate(Group = sample(rep(c("A", "B"), length.out = n_distinct(Transect, Station)))[
        match(paste(Transect, Station), unique(paste(Transect, Station)))
      ]) %>%
      filter((Group == "A" & Tile.Type == "R") | (Group == "B" & Tile.Type == "S"))
    
    return(resampled_data)
  }
  
  # Repeat the resampling process 1000 times and store results
  resampled_datasets_across <- replicate(1000, resample_data_across(Stability_data), simplify = FALSE)
  
  # Combine all resamples into one large dataframe with an identifier
  resampled_df_across <- bind_rows(resampled_datasets_across, .id = "resample_id")

# Compute bootstrapping per dataset

  bootstrap_lms <- function(data_input) {
  
  # Function to fit models on a given resampled dataset and extract estimates
  fit_models <- function(data) {
    # Linear Models (No Random Effects, Since We Are Bootstrapping)
    
    Pop_stab_model <- lm(log(Pop_stab) ~  
                           Emersion.Rate + I(Emersion.Rate^2) +
                           Emersion.Rate * Tile.Type,
                         data = data)
    
    Asynchrony_model <- lm(Asynchrony_corr_time ~  
                             Emersion.Rate * Tile.Type,
                           data = data)
    
    Stat_av_model <- lm(Statistical_averaging ~  
                          Emersion.Rate * Tile.Type,
                        data = data)
    
    Comp_model <- lm(Compositional_similarity ~  
                       Emersion.Rate * I(Emersion.Rate^2) +
                       Emersion.Rate * Tile.Type,
                     data = data)
    
    Cov_model <- lm(log(Cover_stability_alpha_res) ~  
                      Emersion.Rate * I(Emersion.Rate^2) +
                      Emersion.Rate * Tile.Type,
                    data = data)
    
    # Extract summary statistics
    results <- list(
      Pop_stab = coef(summary(Pop_stab_model)),
      Asynchrony = coef(summary(Asynchrony_model)),
      Stat_av = coef(summary(Stat_av_model)),
      Comp = coef(summary(Comp_model)),
      Cov = coef(summary(Cov_model)),
      AIC = c(
        AIC(Pop_stab_model),
        AIC(Asynchrony_model),
        AIC(Stat_av_model),
        AIC(Comp_model),
        AIC(Cov_model)
      )
    )
    
    return(results)
  }
  
  # Apply the function to all resampled datasets
  bootstrapped_results <- lapply(split(data_input, data_input$resample_id), fit_models)
  
  # Convert results to a structured dataframe
  extract_coefs <- function(model_list, model_name) {
    do.call(rbind, lapply(seq_along(model_list), function(i) {
      coef_df <- as.data.frame(model_list[[i]][[model_name]])
      coef_df$resample_id <- names(model_list)[i]
      coef_df$term <- rownames(coef_df)
      coef_df$model <- model_name
      coef_df
    }))
  }
  
  # Extract coefficients and standard errors from all models
  bootstrapped_coefs <- bind_rows(
    extract_coefs(bootstrapped_results, "Pop_stab"),
    extract_coefs(bootstrapped_results, "Asynchrony"),
    extract_coefs(bootstrapped_results, "Stat_av"),
    extract_coefs(bootstrapped_results, "Comp"),
    extract_coefs(bootstrapped_results, "Cov")
  )
  
  # Extract AIC values
  bootstrapped_AICs <- do.call(rbind, lapply(seq_along(bootstrapped_results), function(i) {
    data.frame(
      resample_id = names(bootstrapped_results)[i],
      model = c("Pop_stab", "Asynchrony", "Stat_av", "Comp", "Cov"),
      AIC = bootstrapped_results[[i]]$AIC
    )
  }))
  
  # View Results
  print(head(bootstrapped_coefs))  # Bootstrapped coefficients
  print(head(bootstrapped_AICs))   # Bootstrapped AIC values
  
  # Summarize bootstrap estimates for each coefficient
  summary_bootstrapped_coefs <- bootstrapped_coefs %>%
    group_by(model, term) %>%
    summarise(
      mean_estimate = mean(Estimate),  
      lower_CI = quantile(Estimate, 0.025),  # 2.5% percentile (Lower CI)
      upper_CI = quantile(Estimate, 0.975),  # 97.5% percentile (Upper CI)
      Sig = ifelse(lower_CI > 0 & upper_CI > 0 | lower_CI < 0 & upper_CI < 0, "Yes", "No"),  # Significance check
      .groups = "drop"
    )
  
  # View summarized results
  print(summary_bootstrapped_coefs)
  
  # Merge significance information into bootstrapped_coefs
  bootstrapped_coefs <- bootstrapped_coefs %>%
    left_join(summary_bootstrapped_coefs %>% select(model, term, Sig), by = c("model", "term"))
  
  return(bootstrapped_coefs)
  
  }
  
  boostrapped_lms_within <- bootstrap_lms(resampled_df_within)
  boostrapped_lms_across <- bootstrap_lms(resampled_df_across)
  
  bootstrapped_coefs <- dplyr::bind_rows(
    boostrapped_lms_within %>% dplyr::mutate(source = "within"),
    boostrapped_lms_across %>% dplyr::mutate(source = "across")
  )

bootstrapped_coefs_summary <- bootstrapped_coefs %>%
    group_by(model, term, source, Sig) %>%
    summarise(
      mean_estimate = mean(Estimate),
      lower_CI = quantile(Estimate, 0.025),
      upper_CI = quantile(Estimate, 0.975),
      .groups = "drop"
    ) %>%
    
    filter(term != "(Intercept)") %>%
    
    mutate(
      term_clean = case_when(
        term == "Emersion.Rate" ~ "Emersion rate (linear)",
        term == "I(Emersion.Rate^2)" ~ "Emersion rate (quadratic)",
        term == "Tile.TypeS" ~ "Heterogeneity",
        term == "Emersion.Rate:Tile.TypeS" ~ "Emersion × Heterogeneity",
        term == "Emersion.Rate:I(Emersion.Rate^2)" ~ "Emersion linear × quadratic",
        TRUE ~ term  # fallback in case new levels show up
      ),
      source_clean = case_when(
        source == "across" ~ "Across stations",
        source == "within" ~ "Within stations",
        TRUE ~ source
      ),
      model_clean = case_when(
        model == "Asynchrony" ~ "Species asynchrony",
        model == "Comp" ~ "Compositinal stability",
        model == "Cov" ~ "Cover stability",
        model == "Pop_stab" ~ "Population stability",
        model== "Stat_av" ~ "Statistical averaging"
      )
    )
  
ggplot(bootstrapped_coefs_summary, aes(x = term_clean, y = mean_estimate, color = Sig, shape = source_clean)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point(
      size = 3,
      position = position_dodge(width = 0.75),
      alpha = 0.5
    ) +
    geom_errorbar(
      aes(ymin = lower_CI, ymax = upper_CI),
      width = 0.1,
      position = position_dodge(width = 0.75)
    ) +
    facet_wrap(~model_clean, scales = "free_y", nrow=3) +
    scale_color_manual(values = c("Yes" = "red", "No" = "darkgrey")) +
    labs(
      title = "Bootstrapped SEM Coefficients: Point Estimates with 95% CIs",
      x = "Predictor",
      y = "Mean Coefficient Estimate",
      color = "Significant?",
      shape = "Heterogeneity comparison"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 10, face = "bold"),
      plot.title.position = "plot",
      plot.margin = margin(10, 10, 10, 25)  # Expand left margin
    )
  
#
# 3 - Remove whelks from SEM and compare results ---------

Get_multigroup_coefs <- function(input_dataset) {
    
    extract_dataset<-input_dataset
    
    output<-extract_dataset$group.coefs %>%
      Map(cbind, ., EmGroup=names(.)) %>% # add column indicating which Tile and date each dataframe correspondsto
      do.call(rbind,.)
    
    return(output)
  }
  
# Original multigroup SEM
  
SEM_model_out<-psem(
    
    lme(Consumer_mean_log_z ~ Emersion.Rate_z + Tile.Type2,  random = ~1|Transect/Station, Stability_data_SEM),
    lme(Barnacles1_mean_z ~ Emersion.Rate_z + Tile.Type2 + Consumer_mean_log_z, random = ~1|Transect/Station, Stability_data_SEM),
    lme(Barnacles2_mean_sqrt_z ~ Emersion.Rate_z + Tile.Type2 + Consumer_mean_log_z,  random = ~1|Transect/Station, Stability_data_SEM),
    lme(Eph_macroalgae_mean_sqrt_z ~ Emersion.Rate_z + Tile.Type2 + Consumer_mean_log_z,  random = ~1|Transect/Station, Stability_data_SEM),
    
    lme(S_mean_log_z ~ Emersion.Rate_z + Tile.Type2 + Consumer_mean_log_z,  random = ~1|Transect/Station, Stability_data_SEM),
    lme(Pielou_mean_z ~ Emersion.Rate_z + Tile.Type2 + Consumer_mean_log_z,  random = ~1|Transect/Station, Stability_data_SEM),
    
    lm(Pop_stab_log_z ~ Emersion.Rate_z + Tile.Type2 + S_mean_log_z + Pielou_mean_z + Barnacles1_mean_z + Barnacles2_mean_sqrt_z + Consumer_mean_log_z + Eph_macroalgae_mean_sqrt_z, data=Stability_data_SEM),
    lme(Asynchrony_corr_time_z ~ Emersion.Rate_z + Tile.Type2 + S_mean_log_z + Pielou_mean_z + Barnacles1_mean_z + Barnacles2_mean_sqrt_z + Consumer_mean_log_z + Eph_macroalgae_mean_sqrt_z, random = ~1|Transect/Station, Stability_data_SEM),
    lme(Statistical_averaging_z ~ Emersion.Rate_z + Tile.Type2 + S_mean_log_z + Pielou_mean_z + Barnacles1_mean_z + Barnacles2_mean_sqrt_z + Consumer_mean_log_z + Eph_macroalgae_mean_sqrt_z, random = ~1|Transect/Station, data=Stability_data_SEM),
    
    lm(Compositional_similarity_z ~ Emersion.Rate_z + Tile.Type2+ Asynchrony_corr_time_z + Pop_stab_log_z +  Statistical_averaging_z + S_mean_log_z + Pielou_mean_z + Barnacles1_mean_z + Barnacles2_mean_sqrt_z + Consumer_mean_log_z + Eph_macroalgae_mean_sqrt_z, Stability_data_SEM),
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
  
  summary(SEM_model_out)
  
  SEM_model_out.multigroup <- multigroup(SEM_model_out, group = "EmGroup") # it works under piecewiseSEM Version 2.2.0

  SEM_multigroup_coefs<-Get_multigroup_coefs(SEM_model_out.multigroup) %>%
    mutate(Std.Estimate=as.numeric(Std.Estimate)) %>%
    dplyr::select(EmGroup, Response, Predictor, Std.Estimate,P.Value,Var.9)
  
write.csv(SEM_multigroup_coefs, "Original_SEM.csv")  
  
# Remove whelks

SEM_model_out_nowhelks<-psem(
  
  lme(Limpet_mean_log_z ~ Emersion.Rate_z + Tile.Type2,  random = ~1|Transect/Station, Stability_data_SEM),
  lme(Barnacles1_mean_z ~ Emersion.Rate_z + Tile.Type2 + Limpet_mean_log_z, random = ~1|Transect/Station, Stability_data_SEM),
  lme(Barnacles2_mean_sqrt_z ~ Emersion.Rate_z + Tile.Type2 + Limpet_mean_log_z,  random = ~1|Transect/Station, Stability_data_SEM),
  lme(Eph_macroalgae_mean_sqrt_z ~ Emersion.Rate_z + Tile.Type2 + Limpet_mean_log_z,  random = ~1|Transect/Station, Stability_data_SEM),
  lme(S_mean_log_z ~ Emersion.Rate_z + Tile.Type2 + Limpet_mean_log_z,  random = ~1|Transect/Station, Stability_data_SEM),
  lme(Pielou_mean_z ~ Emersion.Rate_z + Tile.Type2 + Limpet_mean_log_z,  random = ~1|Transect/Station, Stability_data_SEM),
  lm(Pop_stab_log_z ~ Emersion.Rate_z + Tile.Type2 + S_mean_log_z + Pielou_mean_z + Barnacles1_mean_z + Barnacles2_mean_sqrt_z + Limpet_mean_log_z + Eph_macroalgae_mean_sqrt_z, data=Stability_data_SEM),
  lme(Asynchrony_corr_time_z ~ Emersion.Rate_z + Tile.Type2 + S_mean_log_z + Pielou_mean_z + Barnacles1_mean_z + Barnacles2_mean_sqrt_z + Limpet_mean_log_z + Eph_macroalgae_mean_sqrt_z, random = ~1|Transect/Station, Stability_data_SEM),
  lme(Statistical_averaging_z ~ Emersion.Rate_z + Tile.Type2 + S_mean_log_z + Pielou_mean_z + Barnacles1_mean_z + Barnacles2_mean_sqrt_z + Limpet_mean_log_z + Eph_macroalgae_mean_sqrt_z, random = ~1|Transect/Station, data=Stability_data_SEM),
  lm(Compositional_similarity_z ~ Emersion.Rate_z + Tile.Type2+ Asynchrony_corr_time_z + Pop_stab_log_z +  Statistical_averaging_z + S_mean_log_z + Pielou_mean_z + Barnacles1_mean_z + Barnacles2_mean_sqrt_z + Limpet_mean_log_z + Eph_macroalgae_mean_sqrt_z, Stability_data_SEM),
  lme(Cover_stability_alpha_res_log_z ~ Emersion.Rate_z + Tile.Type2 + Asynchrony_corr_time_z + Pop_stab_log_z +  Statistical_averaging_z + Compositional_similarity_z,  random = ~1|Transect/Station, Stability_data_SEM),
  
  Limpet_mean_log_z %~~% Cover_stability_alpha_res_log_z,
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

SEM_model_out.multigroup_nowhelks <- multigroup(SEM_model_out_nowhelks, group = "EmGroup") # it works under piecewiseSEM Version 2.2.0

SEM_multigroup_coefs_nowhelks <- Get_multigroup_coefs(SEM_model_out.multigroup_nowhelks) %>%
  mutate(Std.Estimate=as.numeric(Std.Estimate)) %>%
  dplyr::select(EmGroup, Response, Predictor, Std.Estimate,P.Value,Var.9)


write.csv(SEM_multigroup_coefs_nowhelks, "Original_nowhelks_SEM.csv")

# Compare with Original SEM model

# Combine datasets and add a column indicating their source
combined_data_nowhelks <- bind_rows(
  SEM_multigroup_coefs_nowhelks  %>% dplyr::select(EmGroup, Response, Predictor, "mean_estimate"=Std.Estimate, "Sig"=Var.9) %>% 
    mutate(
      Sig=ifelse(Sig!="", "Yes", "No"),
      lower_CI=mean_estimate,
      upper_CI=mean_estimate,
      dataset = "Original_nowhelks"),
  SEM_multigroup_coefs %>% dplyr::select(EmGroup, Response, Predictor, "mean_estimate"=Std.Estimate, "Sig"=Var.9) %>% 
    mutate(
      Sig=ifelse(Sig!="", "Yes", "No"),
      lower_CI=mean_estimate,
      upper_CI=mean_estimate,
      dataset = "Original_complete")
) %>%
  mutate(
    Predictor = ifelse(Predictor=="Limpet_mean_log_z", "Consumer_mean_log_z", Predictor),
    Response = ifelse(Response=="Limpet_mean_log_z", "Consumer_mean_log_z", Response)
  )

# Create the plot
ggplot(combined_data_nowhelks %>% filter(Response!="")
       , aes(x = Predictor, y = mean_estimate, color = Sig, shape = dataset)) +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_point(size = 3, position = position_dodge(width = 0.75), alpha=0.5) +  # Mean estimates
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI),
                width = 0.1, position = position_dodge(width = 0.75)) +  # 95% CI error bars
  facet_wrap(~ Response + EmGroup, scales = "free_y") +  # Separate by EmGroup & Response
  scale_color_manual(values = c("Yes" = "red", "No" = "darkgrey")) +  # Significant = red, Non-significant = black
  labs(title = "Compare SEM Results With Dataset with only limpets",
       x = "Predictor",
       y = "Mean Standardized Estimate",
       color = "Significant?",
       shape = "Dataset") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for clarity
        strip.text = element_text(size = 10, face = "bold"))

#
# 4 - Simplify Original SEM ---------

SEM_model_out_reduced<-psem(
  
  lme(Consumer_mean_log_z ~ Emersion.Rate_z + Tile.Type2,  random = ~1|Transect/Station, Stability_data_SEM),
  lme(Barnacles1_mean_z ~ Emersion.Rate_z + Tile.Type2 + Consumer_mean_log_z, random = ~1|Transect/Station, Stability_data_SEM),
  lme(Barnacles2_mean_sqrt_z ~ Emersion.Rate_z + Tile.Type2 + Consumer_mean_log_z,  random = ~1|Transect/Station, Stability_data_SEM),
  lme(Eph_macroalgae_mean_sqrt_z ~ Emersion.Rate_z + Tile.Type2 + Consumer_mean_log_z,  random = ~1|Transect/Station, Stability_data_SEM),
  
  lme(S_mean_log_z ~ Emersion.Rate_z + Tile.Type2 + Consumer_mean_log_z,  random = ~1|Transect/Station, Stability_data_SEM),
  
  lm(Pop_stab_log_z ~ Emersion.Rate_z + Tile.Type2 + S_mean_log_z + Barnacles1_mean_z + Barnacles2_mean_sqrt_z + Consumer_mean_log_z + Eph_macroalgae_mean_sqrt_z, data=Stability_data_SEM),
  lme(Asynchrony_corr_time_z ~ Emersion.Rate_z + Tile.Type2 + S_mean_log_z + Barnacles1_mean_z + Barnacles2_mean_sqrt_z + Consumer_mean_log_z + Eph_macroalgae_mean_sqrt_z, random = ~1|Transect/Station, Stability_data_SEM),
  
  lm(Compositional_similarity_z ~ Emersion.Rate_z + Tile.Type2+ Asynchrony_corr_time_z + Pop_stab_log_z +  S_mean_log_z + Barnacles1_mean_z + Barnacles2_mean_sqrt_z + Consumer_mean_log_z + Eph_macroalgae_mean_sqrt_z, Stability_data_SEM),
  lme(Cover_stability_alpha_res_log_z ~ Emersion.Rate_z + Tile.Type2 + Asynchrony_corr_time_z + Pop_stab_log_z + Compositional_similarity_z,  random = ~1|Transect/Station, Stability_data_SEM),
  
  Consumer_mean_log_z %~~% Cover_stability_alpha_res_log_z,
  
  Barnacles1_mean_z %~~% Barnacles2_mean_sqrt_z,
  Barnacles1_mean_z %~~% Cover_stability_alpha_res_log_z,
  Barnacles1_mean_z %~~% Eph_macroalgae_mean_sqrt_z,
  Barnacles1_mean_z %~~% S_mean_log_z,
  
  Barnacles2_mean_sqrt_z %~~% Cover_stability_alpha_res_log_z,
  Barnacles2_mean_sqrt_z %~~% S_mean_log_z,
  Barnacles2_mean_sqrt_z %~~% Eph_macroalgae_mean_sqrt_z,
  
  Eph_macroalgae_mean_sqrt_z %~~% S_mean_log_z,
  
  S_mean_log_z %~~% Cover_stability_alpha_res_log_z,
  
  Pop_stab_log_z %~~% Compositional_similarity_z,
  Pop_stab_log_z %~~% Asynchrony_corr_time_z,
  
  Asynchrony_corr_time_z %~~% Compositional_similarity_z
  
)

SEM_model_out.multigroup_reduced <- multigroup(SEM_model_out_reduced, group = "EmGroup") # it works under piecewiseSEM Version 2.2.0

SEM_multigroup_coefs_reduced<-Get_multigroup_coefs(SEM_model_out.multigroup_reduced) %>%
  mutate(Std.Estimate=as.numeric(Std.Estimate)) %>%
  dplyr::select(EmGroup, Response, Predictor, Std.Estimate,P.Value,Var.9)

write.csv(SEM_multigroup_coefs_reduced, "Original_reduced_SEM.csv")

#
# 5 - Bootstrapping Multigroup SEM -----------------

# Safe wrapper for lme/lm
safe_model <- function(expr) {
  tryCatch(
    expr,
    error = function(e) {
      message("Model fitting failed: ", conditionMessage(e))
      return(NULL)
    }
  )
}

set.seed(123)  # For reproducibility

generate_unique_resamples <- function(data, n_resamples = 1000, min_keep = 46, max_keep = 66) {
  unique_datasets <- list()  # Store unique datasets
  unique_signatures <- character()  # Track dataset signatures for fast lookup
  
  while (length(unique_datasets) < n_resamples) {
    # Step A: Randomly choose how many rows to keep (between 57 and 66)
    num_keep <- sample(min_keep:max_keep, 1)
    
    # Step B: Select those rows to keep
    keep_rows <- sample(nrow(data), num_keep, replace = FALSE)
    
    # Step C: Select the remaining rows randomly with replacement
    replace_rows <- sample(nrow(data), 67 - num_keep, replace = TRUE)
    
    # Create the final dataset
    resampled_df <- bind_rows(data[keep_rows, ], data[replace_rows, ])
    
    # Create a unique signature for this dataset
    dataset_signature <- paste(sort(keep_rows), collapse = "_")  # Use only fixed rows to define uniqueness
    
    # Avoid duplicates
    if (!(dataset_signature %in% unique_signatures)) {
      unique_datasets[[length(unique_datasets) + 1]] <- resampled_df
      unique_signatures <- c(unique_signatures, dataset_signature)
    }
  }
  
  # Convert the list to a single dataframe
  final_resampled_df <- bind_rows(unique_datasets, .id = "resample_id")
  
  return(final_resampled_df)
}
SEM_resampled_df_unique <- generate_unique_resamples(Stability_data_SEM, n_resamples = 5000)

SEM_resampled_df_unique2 <- SEM_resampled_df_unique %>%
  dplyr::select(
    resample_id,
    Transect, Station, Tile.Type2, Emersion.Rate_z, EmGroup,
    Consumer_mean_log_z, Barnacles1_mean_z, Barnacles2_mean_sqrt_z, Eph_macroalgae_mean_sqrt_z,
    S_mean_log_z, Pielou_mean_z, Pop_stab_log_z, Asynchrony_corr_time_z,
    Statistical_averaging_z, Compositional_similarity_z, Cover_stability_alpha_res_log_z
  )

resample_id <- unique(SEM_resampled_df_unique2$resample_id)  # Get unique resample IDs

# Create an empty list to store results
SEM_results_list <- list()

# Bootstrap multigroup SEM loop
for (i in resample_id) {
  
  if (length(SEM_results_list) >= 1000) {
    message("Stopping loop: 1000 valid models collected.")
    break
  }
  
  SEM_filtered_df <- SEM_resampled_df_unique2 %>%
    filter(resample_id == i)
  
  # Fit all models safely
  m1  <- safe_model(lme(Consumer_mean_log_z ~ Emersion.Rate_z + Tile.Type2, random = ~1 | Transect / Station, data = data.frame(SEM_filtered_df)))
  m2  <- safe_model(lme(Barnacles1_mean_z ~ Emersion.Rate_z + Tile.Type2 + Consumer_mean_log_z, random = ~1 | Transect / Station, data = data.frame(SEM_filtered_df)))
  m3  <- safe_model(lme(Barnacles2_mean_sqrt_z ~ Emersion.Rate_z + Tile.Type2 + Consumer_mean_log_z, random = ~1 | Transect / Station, data = data.frame(SEM_filtered_df)))
  m4  <- safe_model(lme(Eph_macroalgae_mean_sqrt_z ~ Emersion.Rate_z + Tile.Type2 + Consumer_mean_log_z, random = ~1 | Transect / Station, data = data.frame(SEM_filtered_df)))
  m5  <- safe_model(lme(S_mean_log_z ~ Emersion.Rate_z + Tile.Type2 + Consumer_mean_log_z, random = ~1 | Transect / Station, data = data.frame(SEM_filtered_df)))
  m6  <- safe_model(lme(Pielou_mean_z ~ Emersion.Rate_z + Tile.Type2 + Consumer_mean_log_z, random = ~1 | Transect / Station, data = data.frame(SEM_filtered_df)))
  m7  <- safe_model(lm(Pop_stab_log_z ~ Emersion.Rate_z + Tile.Type2 + S_mean_log_z + Pielou_mean_z + Barnacles1_mean_z + Barnacles2_mean_sqrt_z + Consumer_mean_log_z + Eph_macroalgae_mean_sqrt_z, data = data.frame(SEM_filtered_df)))
  m8  <- safe_model(lme(Asynchrony_corr_time_z ~ Emersion.Rate_z + Tile.Type2 + S_mean_log_z + Pielou_mean_z + Barnacles1_mean_z + Barnacles2_mean_sqrt_z + Consumer_mean_log_z + Eph_macroalgae_mean_sqrt_z, random = ~1 | Transect / Station, data = data.frame(SEM_filtered_df)))
  m9  <- safe_model(lme(Statistical_averaging_z ~ Emersion.Rate_z + Tile.Type2 + S_mean_log_z + Pielou_mean_z + Barnacles1_mean_z + Barnacles2_mean_sqrt_z + Consumer_mean_log_z + Eph_macroalgae_mean_sqrt_z, random = ~1 | Transect / Station, data = data.frame(SEM_filtered_df)))
  m10 <- safe_model(lm(Compositional_similarity_z ~ Emersion.Rate_z + Tile.Type2 + Asynchrony_corr_time_z + Pop_stab_log_z + Statistical_averaging_z + S_mean_log_z + Pielou_mean_z + Barnacles1_mean_z + Barnacles2_mean_sqrt_z + Consumer_mean_log_z + Eph_macroalgae_mean_sqrt_z, data = data.frame(SEM_filtered_df)))
  m11 <- safe_model(lme(Cover_stability_alpha_res_log_z ~ Emersion.Rate_z + Tile.Type2 + Asynchrony_corr_time_z + Pop_stab_log_z + Statistical_averaging_z + Compositional_similarity_z, random = ~1 | Transect / Station, data = data.frame(SEM_filtered_df)))
  
  # Skip iteration if any model failed
  if (any(sapply(list(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11), is.null))) {
    message("Skipping resample_id ", i, ": one or more models failed to fit.")
    next
  }
  
  # Try to fit the SEM model
  SEM_model_out <- tryCatch({
    psem(
      m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11,
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
      Compositional_similarity_z %~~% Statistical_averaging_z,
      data = SEM_filtered_df
    )
  }, error = function(e) {
    message("Error in psem() for resample_id ", i, ": ", conditionMessage(e))
    return(NULL)
  })
  
  if (is.null(SEM_model_out)) next
  
  # Multigroup analysis
  SEM_model_out_multigroup <- tryCatch({
    multigroup(SEM_model_out, group = "EmGroup")
  }, error = function(e) {
    message("Error in multigroup analysis for resample_id ", i, ": ", conditionMessage(e))
    return(NULL)
  })
  
  if (is.null(SEM_model_out_multigroup)) next
  
  # Coefficient extraction
  Get_multigroup_coefs <- function(input_dataset) {
    input_dataset$group.coefs %>%
      Map(cbind, ., EmGroup = names(.)) %>%
      do.call(rbind, .)
  }
  
  SEM_multigroup_coefs <- tryCatch({
    Get_multigroup_coefs(SEM_model_out_multigroup) %>%
      mutate(Std.Estimate = as.numeric(Std.Estimate)) %>%
      dplyr::select(EmGroup, Response, Predictor, Std.Estimate, P.Value, Var.9) %>%
      mutate(resample_id = i)
  }, error = function(e) {
    message("Error extracting coefficients for resample_id ", i, ": ", conditionMessage(e))
    return(NULL)
  })
  
  if (is.null(SEM_multigroup_coefs)) next
  
  SEM_results_list[[as.character(i)]] <- SEM_multigroup_coefs
}

bootstrapped_SEM_results1 <- SEM_results_list %>%
  map(~ mutate(.x, P.Value = as.numeric(P.Value))) %>%  # Convert P.Value to numeric
  bind_rows()  # Combine results into a single dataframe

summary_SEM_results1 <- bootstrapped_SEM_results1 %>%
  group_by(EmGroup, Response, Predictor) %>%
  summarise(
    mean_estimate = mean(Std.Estimate),  
    lower_CI = quantile(Std.Estimate, 0.025),
    upper_CI = quantile(Std.Estimate, 0.975),
    Sig = ifelse(lower_CI > 0 & upper_CI > 0 | lower_CI < 0 & upper_CI < 0, "Yes", "No"),
    .groups = "drop"
  )

write.csv(summary_SEM_results1, "Bootstrapped_SEM_Original.csv")

ggplot(summary_SEM_results1, aes(x = Predictor, y = mean_estimate, color = Sig)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +  # Mean estimates
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), 
                width = 0.2, position = position_dodge(width = 0.5)) +  # 95% CI error bars
  facet_wrap(~ EmGroup + Response, scales = "free_y") +  # Separate by EmGroup & Response
  scale_color_manual(values = c("Yes" = "red", "No" = "black")) +  # Significant = red, Non-significant = black
  labs(title = "Bootstrapped SEM Results", 
       x = "Predictor", 
       y = "Mean Standardized Estimate", 
       color = "Significant?") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for clarity
        strip.text = element_text(size = 10, face = "bold"))

#
# 6 - Remove unnecessary links in the SEM ------------

fit_multigroup_SEM2 <- function(data){
  
  resample_id <- unique(data$resample_id)  # Get unique resample IDs
  
  # Create an empty list to store results
  SEM_results_list <- list()
  
  for(i in resample_id){
    
    # Stop the loop if the list reaches 100 elements
    if (length(SEM_results_list) >= 1000) {
      message("Stopping loop: 100 valid models collected.")
      break
    }
    
    SEM_filtered_df <- data %>% 
      filter(resample_id == i) %>%
      dplyr::select(
        Transect, Station, Tile.Type2, Emersion.Rate_z, EmGroup,
        Consumer_mean_log_z, Barnacles1_mean_z, Barnacles2_mean_sqrt_z, Eph_macroalgae_mean_sqrt_z,
        S_mean_log_z, Pielou_mean_z, Pop_stab_log_z, Asynchrony_corr_time_z,
        Statistical_averaging_z, Compositional_similarity_z, Cover_stability_alpha_res_log_z
      )
    
    # Try to fit the SEM model and catch errors
    SEM_model_out <- tryCatch({
      psem(
        lme(Consumer_mean_log_z ~ Emersion.Rate_z + Tile.Type2, random = ~1|Transect/Station, SEM_filtered_df),
        lme(Barnacles1_mean_z ~ Emersion.Rate_z + Tile.Type2 + Consumer_mean_log_z, random = ~1|Transect/Station, SEM_filtered_df),
        lme(Barnacles2_mean_sqrt_z ~ Emersion.Rate_z + Tile.Type2 + Consumer_mean_log_z, random = ~1|Transect/Station, SEM_filtered_df),
        lme(Eph_macroalgae_mean_sqrt_z ~ Emersion.Rate_z + Tile.Type2 + Consumer_mean_log_z, random = ~1|Transect/Station, SEM_filtered_df),
        lme(S_mean_log_z ~ Emersion.Rate_z + Tile.Type2 + Consumer_mean_log_z, random = ~1|Transect/Station, SEM_filtered_df),
        lm(Pop_stab_log_z ~ Emersion.Rate_z + Tile.Type2 + S_mean_log_z + Pielou_mean_z + Barnacles1_mean_z + Barnacles2_mean_sqrt_z + Consumer_mean_log_z + Eph_macroalgae_mean_sqrt_z, SEM_filtered_df),
        lme(Asynchrony_corr_time_z ~ Emersion.Rate_z + Tile.Type2 + S_mean_log_z + Pielou_mean_z + Barnacles1_mean_z + Barnacles2_mean_sqrt_z + Consumer_mean_log_z + Eph_macroalgae_mean_sqrt_z, random = ~1|Transect/Station, SEM_filtered_df),
        lm(Compositional_similarity_z ~ Emersion.Rate_z + Tile.Type2+ Asynchrony_corr_time_z + Pop_stab_log_z +  S_mean_log_z + Barnacles1_mean_z + Barnacles2_mean_sqrt_z + Consumer_mean_log_z + Eph_macroalgae_mean_sqrt_z, SEM_filtered_df),
        lme(Cover_stability_alpha_res_log_z ~ Emersion.Rate_z + Tile.Type2 + Asynchrony_corr_time_z + Pop_stab_log_z +  Compositional_similarity_z, random = ~1|Transect/Station, SEM_filtered_df),
        
        Consumer_mean_log_z %~~% Cover_stability_alpha_res_log_z,
        Barnacles1_mean_z %~~% Barnacles2_mean_sqrt_z,
        Barnacles1_mean_z %~~% Cover_stability_alpha_res_log_z,
        Barnacles1_mean_z %~~% Eph_macroalgae_mean_sqrt_z,
        Barnacles1_mean_z %~~% S_mean_log_z,
        Barnacles2_mean_sqrt_z %~~% Cover_stability_alpha_res_log_z,
        Barnacles2_mean_sqrt_z %~~% S_mean_log_z,
        Barnacles2_mean_sqrt_z %~~% Eph_macroalgae_mean_sqrt_z,
        Eph_macroalgae_mean_sqrt_z %~~% S_mean_log_z,
        S_mean_log_z %~~% Cover_stability_alpha_res_log_z,
        Pop_stab_log_z %~~% Compositional_similarity_z,
        Pop_stab_log_z %~~% Asynchrony_corr_time_z,
        Asynchrony_corr_time_z %~~% Compositional_similarity_z
      )
    }, error = function(e) {
      message(paste("Error in resample_id", i, ":", e$message))
      return(NULL)  # Return NULL when an error occurs
    })
    
    # Skip to next iteration if model fitting failed
    if (is.null(SEM_model_out)) next
    
    # Try to run multigroup analysis and catch errors
    SEM_model_out_multigroup <- tryCatch({
      multigroup(SEM_model_out, group = "EmGroup")
    }, error = function(e) {
      message(paste("Error in multigroup analysis for resample_id", i, ":", e$message))
      return(NULL)
    })
    
    # Skip to next iteration if multigroup failed
    if (is.null(SEM_model_out_multigroup)) next
    
    # Function to extract coefficients from multigroup analysis
    Get_multigroup_coefs <- function(input_dataset) {
      extract_dataset <- input_dataset
      output <- extract_dataset$group.coefs %>%
        Map(cbind, ., EmGroup = names(.)) %>%
        do.call(rbind, .)  # Combine into one data frame
      return(output)
    }
    
    # Try to extract coefficients and catch errors
    SEM_multigroup_coefs <- tryCatch({
      Get_multigroup_coefs(SEM_model_out_multigroup) %>%
        mutate(Std.Estimate = as.numeric(Std.Estimate)) %>%
        dplyr::select(EmGroup, Response, Predictor, Std.Estimate, P.Value, Var.9) %>%
        mutate(resample_id = i)  # Add resample ID
    }, error = function(e) {
      message(paste("Error in extracting coefficients for resample_id", i, ":", e$message))
      return(NULL)
    })
    
    # Skip to next iteration if coefficient extraction failed
    if (is.null(SEM_multigroup_coefs)) next
    
    # Store results in list
    SEM_results_list[[as.character(i)]] <- SEM_multigroup_coefs
  }
  
  return(SEM_results_list)  # Return the list instead of binding rows
}

# Run the function
SEM_multigroup_coefs2 <- fit_multigroup_SEM2(SEM_resampled_df_unique)

bootstrapped_SEM_results2 <- SEM_multigroup_coefs2 %>%
  map(~ mutate(.x, P.Value = as.numeric(P.Value))) %>%  # Convert P.Value to numeric
  bind_rows()  # Combine results into a single dataframe

# View results
print(head(bootstrapped_SEM_results))

summary_SEM_results2 <- bootstrapped_SEM_results2 %>%
  group_by(EmGroup, Response, Predictor) %>%
  summarise(
    mean_estimate = mean(Std.Estimate),  
    lower_CI = quantile(Std.Estimate, 0.025),
    upper_CI = quantile(Std.Estimate, 0.975),
    Sig = ifelse(lower_CI > 0 & upper_CI > 0 | lower_CI < 0 & upper_CI < 0, "Yes", "No"),
    .groups = "drop"
  )

ggplot(summary_SEM_results2, aes(x = Predictor, y = mean_estimate, color = Sig)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +  # Mean estimates
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), 
                width = 0.2, position = position_dodge(width = 0.5)) +  # 95% CI error bars
  facet_wrap(~ EmGroup + Response, scales = "free_y") +  # Separate by EmGroup & Response
  scale_color_manual(values = c("Yes" = "red", "No" = "black")) +  # Significant = red, Non-significant = black
  labs(title = "Bootstrapped SEM Results", 
       x = "Predictor", 
       y = "Mean Standardized Estimate", 
       color = "Significant?") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for clarity
        strip.text = element_text(size = 10, face = "bold"))

#
# 7 - Compare with Original SEM -----------------

# Combine datasets and add a column indicating their source
combined_data <- bind_rows(
  summary_SEM_results1 %>% dplyr::select(EmGroup, Response, Predictor, mean_estimate, lower_CI, upper_CI, Sig) %>% mutate(dataset = "Bootstrapped_complete"),
  summary_SEM_results2 %>% dplyr::select(EmGroup, Response, Predictor, mean_estimate, lower_CI, upper_CI, Sig) %>% mutate(dataset = "Bootstrapped_reduced"),
  SEM_multigroup_coefs %>% dplyr::select(EmGroup, Response, Predictor, "mean_estimate"=Std.Estimate, "Sig"=Var.9) %>% 
    mutate(
    Sig=ifelse(Sig!="", "Yes", "No"),
    lower_CI=mean_estimate,
    upper_CI=mean_estimate,
    dataset = "Original_complete"),
  SEM_multigroup_coefs_reduced %>% dplyr::select(EmGroup, Response, Predictor, "mean_estimate"=Std.Estimate, "Sig"=Var.9) %>% 
    mutate(
      Sig=ifelse(Sig!="", "Yes", "No"),
      lower_CI=mean_estimate,
      upper_CI=mean_estimate,
      dataset = "Original_reduced")
)

# Define shape mappings
shape_mapping <- c(
  "Bootstrapped_complete" = 15,  # Filled square
  "Bootstrapped_reduced" = 16,   # Filled circle
  "Original_complete" = 0,       # Plus sign
  "Original_reduced" = 1         # Cross sign
)

# Create the plot
ggplot(combined_data %>% filter(Response!="")
       , aes(x = Predictor, y = mean_estimate, color = Sig, shape = dataset)) +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_point(size = 3, position = position_dodge(width = 0.75), alpha=0.5) +  # Mean estimates
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI),
                width = 0.1, position = position_dodge(width = 0.75)) +  # 95% CI error bars
  facet_wrap(~ Response + EmGroup, scales = "free_y") +  # Separate by EmGroup & Response
  scale_color_manual(values = c("Yes" = "red", "No" = "darkgrey")) +  # Significant = red, Non-significant = black
  scale_shape_manual(values = shape_mapping) +  
  labs(title = "Bootstrapped SEM Results Across Datasets",
       x = "Predictor",
       y = "Mean Standardized Estimate",
       color = "Significant?",
       shape = "Dataset") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for clarity
        strip.text = element_text(size = 10, face = "bold"))
  

# Create a summary graph

combined_data_summary <- combined_data %>%
  separate(dataset, into = c("Method", "Data"), remove=FALSE) %>%
  group_by(Response, Predictor, Data) %>%
  summarise(
    Change_bootstrapped = ifelse(Sig[Method=="Bootstrapped"]=="Yes" & Sig[Method=="Original"]=="No", 1,0),
    Change_original = ifelse(Sig[Method=="Original"]=="Yes" & Sig[Method=="Bootstrapped"]=="No", 1,0)
    ) %>%
  group_by(Data) %>%
  summarise(
    Change_bootstrapped = sum(Change_bootstrapped),
    Change_original = sum (Change_original),
    Total_comparisons = n()
    ) %>%
  pivot_longer(
    cols = Change_bootstrapped:Total_comparisons,
    names_to = "Comparisons",
    values_to = "Number"
  )


combined_data_summary2 <- combined_data %>%
  separate(dataset, into = c("Method", "Data"), remove = FALSE) %>%
  group_by(Response, Predictor, Method) %>%  # Compare within each method now
  summarise(
    Change_complete = ifelse(Sig[Data == "complete"] == "Yes" & Sig[Data == "reduced"] == "No", 1, 0),
    Change_reduced = ifelse(Sig[Data == "reduced"] == "Yes" & Sig[Data == "complete"] == "No", 1, 0),
    .groups = "drop"
  ) %>%
  group_by(Method) %>%  # Summarize across Methods now
  summarise(
    Change_complete = sum(Change_complete),
    Change_reduced = sum(Change_reduced),
    Total_comparisons = n(),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = Change_complete:Total_comparisons,
    names_to = "Comparisons",
    values_to = "Number"
  )

ggplot(combined_data_summary, aes(x = Data, y = Number, fill = Comparisons)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Change_bootstrapped" = "#377EB8",  # Default ggplot Blue
                               "Change_original" = "#E41A1C",    # Default ggplot Red
                               "Total_comparisons" = "lightgrey"),  # Neutral Grey
                    labels = c("Change_bootstrapped" = "Turned significant in bootstrapped model",
                               "Change_original" = "Turned non-significant in bootstrapped model",
                               "Total_comparisons" = "Total number of comparisons assessed")) +
  scale_x_discrete(labels = c("complete" = "Complete Model",
                              "reduced" = "Reduced Model")) +
  labs(title = "Number of tests that changed in Significance due to Bootstrapping",
       x = "Compared Models",
       y = "Number of Changes in Significance",
       fill = "Comparisons") +
  theme_minimal()

ggplot(combined_data_summary2, aes(x = Method, y = Number, fill = Comparisons)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Change_complete" = "#E41A1C",  # Default ggplot Red
                               "Change_reduced" = "#377EB8",   # Default ggplot Blue
                               "Total_comparisons" = "lightgrey"),  # Neutral Grey
                    labels = c("Change_complete" = "Turned non-significant in reduced model",
                               "Change_reduced" = "Turned significant in reduced model",
                               "Total_comparisons" = "Total number of comparisons assessed")) +
  scale_x_discrete(labels = c("Bootstrapped" = "Bootstrapped Model",
                              "Original" = "Original Model")) +
  labs(title = "Number of tests that changed in Significance when Reducing the SEM Model",
       x = "Compared Models",
       y = "Number of Changes in Significance",
       fill = "Comparisons") +
  theme_minimal()



#

# 8 - Obtain bootstrapped CIs for SEM estimates -----------

Bootstrapped_SEM_Original <- read.csv("/Users/jordisola/Library/Mobile Documents/com~apple~CloudDocs/Documents/PhD_thesis/Habitat_stability/Data_analysis/R_outputs/Bootstrapped_SEM_Original.csv")

Obtain_SEM_effects_unique2<-function(input_dataset){
  
  SEM_loop<-input_dataset %>%
    ungroup() %>%
    dplyr::select(EmGroup,Response,Predictor,Std.Estimate) %>%
    mutate(
      Predictor1=NA,
      Predictor2=NA,
      Predictor3=NA,
      Predictor4=NA,
      Predictor5=NA,
      Predictor6=NA,
      Predictor7=NA,
      Predictor8=NA,
      Predictor9=NA,
      Predictor10=NA,
      Std.Estimate1=NA,
      Std.Estimate2=NA,
      Std.Estimate3=NA,
      Std.Estimate4=NA,
      Std.Estimate5=NA,
      Std.Estimate6=NA,
      Std.Estimate7=NA,
      Std.Estimate8=NA,
      Std.Estimate9=NA,
      Std.Estimate10=NA,
      Std.Estimate.Int=NA,
    ) %>%
    relocate(Std.Estimate,.after=Predictor10)
  
  
  SEM_out<-SEM_loop[0,]
  
  Check_predictors<- function(input_dataset1){
    
    input_dataset<-input_dataset1<-SEM_loop %>%filter(EmGroup==a)
    
    data_nrow_start<-1
    data_nrow<-m<-nrow(input_dataset)
    data_nrow2<-q<-nrow(input_dataset)
    
    if(nrow(input_dataset)>0){
      
      for(j in 3:13){
        
        if(nrow(input_dataset)!=data_nrow){
          data_nrow_start<-data_nrow2
          data_nrow2<-q<-nrow(input_dataset)
        }
        
        for(q in data_nrow_start:data_nrow2){  #
          
          for(m in 1:data_nrow){ 
            
            if(!is.na(input_dataset[q,j])){
              
              if(input_dataset[q,j]==input_dataset[m,"Response"]){
                
                t=nrow(input_dataset)
                
                input_dataset[t+1,]<-input_dataset[q,]
                
                p=min(which(is.na(input_dataset[t+1,])))
                
                input_dataset[t+1,p]<-input_dataset[m,"Predictor"]
                
                input_dataset[t+1,p+11]<-input_dataset[m,"Std.Estimate"]
                
                
              }
            }
          }
        }
        
      }
      
    }
    
    return(input_dataset)
    
  }
  
  
  #Responses<-c("Cover_stability_gamma",     "Composition_gamma",         "Asynchrony_corr_spat",      "Turnover_alpha_spatial_av", "Turnover_alpha_time_av",    "Tile_stability_spat"    )
  
  for(a in unique(SEM_loop$EmGroup)){ 
    
    SEM_group<-SEM_loop %>%
      filter(EmGroup==a)%>%
      Check_predictors
    
    SEM_out<-rbind(SEM_out, SEM_group)
    
  } 
  
  #Remove repeated pathways
  
  SEM_out2<-SEM_out %>%
    #remove repeated pathways
    distinct(EmGroup,Response,Predictor,Predictor1,Predictor2,Predictor3,Predictor4,Predictor5,Predictor6,Predictor7,Predictor8,Predictor9,Predictor10,.keep_all=TRUE) %>%
    
    #remove nested pathways
    
    group_by(EmGroup,Response,Predictor) %>%
    mutate(Predictor_n=ifelse(!is.na(Predictor),n(),NA))%>%
    
    group_by(EmGroup,Response,Predictor,Predictor1) %>%
    mutate(Predictor1_n=ifelse(!is.na(Predictor1),n(),NA))%>%
    
    group_by(EmGroup,Response,Predictor,Predictor1,Predictor2) %>%
    mutate(Predictor2_n=ifelse(!is.na(Predictor2),n(),NA))%>%
    
    group_by(EmGroup,Response,Predictor,Predictor1,Predictor2,Predictor3) %>%
    mutate(Predictor3_n=ifelse(!is.na(Predictor3),n(),NA))%>%
    
    group_by(EmGroup,Response,Predictor,Predictor1,Predictor2,Predictor3,Predictor4) %>%
    mutate(Predictor4_n=ifelse(!is.na(Predictor4),n(),NA))%>%
    
    group_by(EmGroup,Response,Predictor,Predictor1,Predictor2,Predictor3,Predictor4,Predictor5) %>%
    mutate(Predictor5_n=ifelse(!is.na(Predictor5),n(),NA))%>%
    
    group_by(EmGroup,Response,Predictor,Predictor1,Predictor2,Predictor3,Predictor4,Predictor5,Predictor6) %>%
    mutate(Predictor6_n=ifelse(!is.na(Predictor6),n(),NA)) %>%
    
    group_by(EmGroup,Response,Predictor,Predictor1,Predictor2,Predictor3,Predictor4,Predictor5,Predictor6,Predictor7) %>%
    mutate(Predictor7_n=ifelse(!is.na(Predictor7),n(),NA)) %>%
    
    group_by(EmGroup,Response,Predictor,Predictor1,Predictor2,Predictor3,Predictor4,Predictor5,Predictor6,Predictor7,Predictor8) %>%
    mutate(Predictor8_n=ifelse(!is.na(Predictor8),n(),NA)) %>%
    
    group_by(EmGroup,Response,Predictor,Predictor1,Predictor2,Predictor3,Predictor4,Predictor5,Predictor6,Predictor7,Predictor8,Predictor9) %>%
    mutate(Predictor9_n=ifelse(!is.na(Predictor9),n(),NA)) %>%
    
    group_by(EmGroup,Response,Predictor,Predictor1,Predictor2,Predictor3,Predictor4,Predictor5,Predictor6,Predictor7,Predictor8,Predictor9,Predictor10) %>%
    mutate(Predictor10_n=ifelse(!is.na(Predictor10),n(),NA)) %>%
    
    mutate(Unique=ifelse(Predictor_n==1 |Predictor1_n==1 | Predictor2_n==1 | Predictor3_n==1 | Predictor4_n==1 | Predictor5_n==1 | Predictor6_n==1
                         | Predictor7_n==1 | Predictor8_n==1 | Predictor9_n==1 | Predictor10_n==1
                         , "Y","N")) %>%
    filter(Unique=="Y") %>%
    
    #remove nested pathways 2
    
    mutate(
      Check_repeated1 = paste0(na.omit(c(EmGroup,Response,Predictor)), collapse=","),
      Check_repeated2 = paste0(na.omit(c(EmGroup,Predictor,Predictor1)), collapse=",")
    ) %>%
    #group_by(EmGroup) %>%
    mutate(
      test_repeated = Check_repeated1 %in% Check_repeated2
    ) %>%
    filter(test_repeated==FALSE)
  
  SEM_out3<-SEM_out2 %>%
    rowwise() %>%
    mutate(
      Last_Predictor=tail(na.omit(c(Predictor,Predictor1,Predictor2,Predictor3,Predictor4,Predictor5,Predictor6,Predictor7,Predictor8,Predictor9,Predictor10)), n=1),
      Path=ifelse(is.na(Predictor1),
                  ifelse(Predictor=="Tile.Type2", "Direct","Covariates"), 
                  ifelse(Last_Predictor=="Tile.Type2","Indirect","Covariates"))
    ) %>%
    rowwise() %>%
    mutate(
      Std.Estimate.Int = prod(na.omit(c(Std.Estimate,Std.Estimate1,Std.Estimate2,Std.Estimate3,Std.Estimate4,Std.Estimate5,Std.Estimate6,Std.Estimate7,Std.Estimate8,Std.Estimate9,Std.Estimate10)))
    ) %>%
    ungroup() %>%
    dplyr::select(EmGroup, Response, Predictor, Last_Predictor, Path, Std.Estimate.Int, Predictor1,Predictor2,Predictor3,Predictor4,Predictor5,Predictor6,Predictor7,Predictor8,Predictor9,Predictor10)
    
  
  return(SEM_out3)
  
}
prepare_sem_effects <- function(data, estimate_col, suffix = NULL) {
  data %>%
    mutate(Std.Estimate = .data[[estimate_col]]) %>%
    Obtain_SEM_effects_unique2() %>%
    rowwise() %>%
    mutate(
      # Safely extract second-last predictor in path
      Second_Last_Predictor = {
        predictors <- na.omit(c(Predictor, Predictor1, Predictor2, Predictor3, Predictor4,
                                Predictor5, Predictor6, Predictor7, Predictor8, Predictor9, Predictor10))
        if (Predictor != Last_Predictor && length(predictors) >= 2) {
          predictors[length(predictors) - 1]
        } else {
          NA_character_
        }
      },
      # Adjust second predictor if Total_cov_sd is the first step
      Second_Predictor = if_else(
        Predictor != Last_Predictor & Predictor == "Total_cov_sd",
        Predictor1,
        Predictor
      )
    ) %>%
    ungroup() %>%
    mutate(
      # Categorize interaction types
      Interaction = case_when(
        Last_Predictor != "Tile.Type2" ~ "Covariates",
        is.na(Second_Last_Predictor) ~ "1",
        Second_Last_Predictor == "Ac_invertebrates_mean" ~ "2",
        Second_Last_Predictor == "Barnacles1_mean" ~ "3",
        Second_Last_Predictor == "S_mean_log" ~ "4",
        Second_Last_Predictor == "Consumer_mean_log" ~ "5",
        TRUE ~ "1"
      ),
      # Pathway classification
      Pathway = case_when(
        Predictor == "Total_cov_mean" ~ "Total_cov_mean",
        Predictor == "Total_cov_sd" ~ "Total_cov_sd",
        TRUE ~ NA_character_
      ),
      # Recode AltGroup from EmGroup ranges
      Tile_alt = recode(EmGroup,
                        "[0.126,0.274)" = "Low",
                        "[0.274,0.383)" = "Mid",
                        "[0.383,0.589]" = "High"),
      # Adjust path label
      Path2 = if_else(Path == "Direct", "Indirect", Path),
      # Ensure factor levels for plotting or ordering
      Path2 = factor(Path2, levels = c("Direct", "Indirect", "Covariates"))
    ) %>%
    # Rename estimate column if suffix is provided
    {
      if (!is.null(suffix)) {
        rename_with(., ~ paste0(., "_", suffix), .cols = c(Std.Estimate.Int))
      } else .
    }
}

Bootstrapped_SEM_Original_Prepared <- Bootstrapped_SEM_Original %>%
  mutate(
    #Response = ifelse(Response=="", "Cover_stability_alpha_res_log_z", Response),
    Sign = ifelse(mean_estimate>0 & lower_CI<0 | mean_estimate<0 & upper_CI>0, "No", "Yes")
  ) %>%
  filter(Sig != "No" & Sign !="No")

# Calculate SEM effects with different estimates
SEM_bootstrapped_effects_Upper <- prepare_sem_effects(Bootstrapped_SEM_Original_Prepared, "upper_CI", "Upper")
SEM_bootstrapped_effects_Lower <- prepare_sem_effects(Bootstrapped_SEM_Original_Prepared, "lower_CI", "Lower")
SEM_bootstrapped_effects_Mean  <- prepare_sem_effects(Bootstrapped_SEM_Original_Prepared, "mean_estimate")

SEM_bootstrapped_effects <- SEM_bootstrapped_effects_Mean %>%
  left_join(SEM_bootstrapped_effects_Lower) %>%
  left_join(SEM_bootstrapped_effects_Upper) %>%
  mutate(Std.Estimate.Int_CI = ifelse(Std.Estimate.Int<0, Std.Estimate.Int_Lower, Std.Estimate.Int_Upper)) %>%
  relocate(Std.Estimate.Int_CI, .after="Std.Estimate.Int")%>% 
  dplyr::select(Tile_alt, Last_Predictor, Predictor, Predictor1, Predictor2, Predictor3, Predictor4, Response, Std.Estimate.Int_CI)
 
SEM_alpha_effects_CI <- SEM_alpha_effects %>%
  dplyr::select(Tile_alt, Last_Predictor, Predictor, Predictor1, Predictor2, Predictor3, Predictor4, Response, Std.Estimate.Int) %>%
  left_join(SEM_bootstrapped_effects) %>%
  mutate(
    Tile_alt = factor(Tile_alt, levels=c("Low", "Mid", "High"))
  )

# Plot: Indirect effects (via Tile.Type2) of each predictor on multiple responses

library(RColorBrewer)
brewer.pal(7, "YlGn")
brewer.pal(5, "BuPu")
brewer.pal(3, "Oranges")

Alpha_cols<-c(brewer.pal(7, "BuPu")[c(2:5)],
              brewer.pal(4, "YlOrRd")[c(1:3)],
              brewer.pal(5, "Greens")[c(2:5)])

response_subset <- c(
  "Compositional_similarity_z",
  "Cover_stability_alpha_res_log_z",
  "Pop_stab_log_z",
  "Asynchrony_corr_time_z"
)


SEM_alpha_effects_CI_stab_mech <- SEM_alpha_effects_CI %>%
  filter(
    Last_Predictor == "Tile.Type2",
    Response %in% response_subset
  ) %>%
  mutate(Response = factor(Response, levels=c(  "Pop_stab_log_z",
                                                "Asynchrony_corr_time_z",
                                                "Compositional_similarity_z",
                                                "Cover_stability_alpha_res_log_z")),
         Predictor=factor(Predictor,levels = c("Barnacles1_mean_z", "Barnacles2_mean_sqrt_z","Eph_macroalgae_mean_sqrt_z","Consumer_mean_log_z","Pielou_mean_z","S_mean_log_z",  
                                      "Pop_stab_log_z","Asynchrony_corr_time_z","Statistical_averaging_z","Compositional_similarity_z","Tile.Type2"))
         )

# Summed indirect effect of each predictor per group
total_effect_by_predictor <- SEM_alpha_effects_CI_stab_mech %>%
  ungroup() %>%
  group_by(Tile_alt, Predictor, Response) %>%
  summarise(Std.Estimate.Int = sum(Std.Estimate.Int), .groups = "drop")

# Add sign, fix missing CI values, and summarize with CI
effect_with_ci_by_predictor <- SEM_alpha_effects_CI_stab_mech %>%
  mutate(
    Sign = ifelse(Std.Estimate.Int > 0, "Positive", "Negative"),
    Std.Estimate.Int_CI = ifelse(is.na(Std.Estimate.Int_CI), Std.Estimate.Int, Std.Estimate.Int_CI)
  ) %>%
  ungroup() %>%
  select(Tile_alt, Predictor, Response, Sign, Std.Estimate.Int, Std.Estimate.Int_CI) %>%
  group_by(Tile_alt, Predictor, Response) %>%
  summarise(
    Std.Estimate.Int = sum(Std.Estimate.Int),
    Std.Estimate.Int_CI = sum(Std.Estimate.Int_CI),
    .groups = "drop"
  )

# Summarize positive/negative effects and compute upper/lower CI for error bars
effect_ci_by_sign <- effect_with_ci_by_predictor %>%
  mutate(Sign = ifelse(Std.Estimate.Int > 0, "Positive", "Negative")) %>%
  group_by(Tile_alt, Response, Sign) %>%
  summarise(
    Std.Estimate.Int = sum(Std.Estimate.Int),
    Std.Estimate.Int_CI = sum(Std.Estimate.Int_CI),
    .groups = "drop"
  ) %>%
  mutate(
    Std.CI = ifelse(Sign == "Positive",
                    Std.Estimate.Int + abs(Std.Estimate.Int - Std.Estimate.Int_CI),
                    Std.Estimate.Int - abs(Std.Estimate.Int - Std.Estimate.Int_CI))
  )

ggplot() +
  # Bar plot of total effects by predictor
  geom_bar(
    data = total_effect_by_predictor,
    mapping = aes(x = Tile_alt, y = Std.Estimate.Int, fill = Predictor),
    stat = "identity",
    position = "stack"
  ) +
  # Error bars showing uncertainty by sign
  geom_errorbar(
    data = effect_ci_by_sign,
    mapping = aes(x = Tile_alt, ymin = Std.CI, ymax = Std.Estimate.Int, group = Sign),
    width = 0.01,
    linewidth = 0.2
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ Response, ncol = 1, labeller = label_wrap_gen(width = 25)) +
  coord_flip() +
  ylim(-0.8,0.8) +
  scale_fill_manual(values = Alpha_cols) +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

# Plot: Indirect effects from each source of variability on cover stability

SEM_alpha_effects_CI_stab_cover <- SEM_alpha_effects_CI %>%
  filter(
    Response == "Cover_stability_alpha_res_log_z"
  )

# Summed indirect effect of each predictor per group
total_effect_by_predictor2 <- SEM_alpha_effects_CI_stab_cover %>%
  ungroup() %>%
  group_by(Tile_alt, Response, Last_Predictor) %>%
  summarise(Std.Estimate.Int=sum(Std.Estimate.Int))

# Add sign, fix missing CI values, and summarize with CI
effect_with_ci_by_predictor2 <- SEM_alpha_effects_CI_stab_cover %>%
  mutate(
    Std.Estimate.Int_CI = ifelse(is.na(Std.Estimate.Int_CI), Std.Estimate.Int, Std.Estimate.Int_CI)
  ) %>%
  ungroup() %>%
  select(Tile_alt, Last_Predictor, Response, Std.Estimate.Int, Std.Estimate.Int_CI) %>%
  group_by(Tile_alt, Response, Last_Predictor) %>%
  summarise(
    Std.Estimate.Int = sum(Std.Estimate.Int),
    Std.Estimate.Int_CI = sum(Std.Estimate.Int_CI),
    .groups = "drop"
  )

# Summarize positive/negative effects and compute upper/lower CI for error bars
effect_ci_by_sign2 <- effect_with_ci_by_predictor2 %>%
  mutate(Sign = ifelse(Std.Estimate.Int > 0, "Positive", "Negative")) %>%
  group_by(Tile_alt, Response, Sign) %>%
  summarise(
    Std.Estimate.Int = sum(Std.Estimate.Int),
    Std.Estimate.Int_CI = sum(Std.Estimate.Int_CI),
    .groups = "drop"
  ) %>%
  mutate(
    Std.CI = ifelse(Sign == "Positive",
                    Std.Estimate.Int + abs(Std.Estimate.Int - Std.Estimate.Int_CI),
                    Std.Estimate.Int - abs(Std.Estimate.Int - Std.Estimate.Int_CI))
  )

ggplot() +
  # Bar plot of total effects by predictor
  geom_bar(
    data = total_effect_by_predictor2,
    mapping = aes(x = Tile_alt, y = Std.Estimate.Int, fill = Last_Predictor),
    stat = "identity",
    position = "stack"
  ) +
  # Error bars showing uncertainty by sign
  geom_errorbar(
    data = effect_ci_by_sign2,
    mapping = aes(x = Tile_alt, ymin = Std.CI, ymax = Std.Estimate.Int, group = Sign),
    width = 0.01,
    linewidth = 0.2
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  ylim(-0.8,0.8) +
  scale_fill_manual(values = Alpha_cols[c(10,9,8)]) +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

#