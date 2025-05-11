library(tidyverse)
library(lme4)
library(DHARMa)
library(lmerTest)
library(piecewiseSEM)
library(ggeffects)
library(effects)

#### LOAD DATA ####

#Stability_data <- read.csv(file.choose())
Stability_data <- read.csv("/Data/Stability_metrics.csv")

# LOAD FUNCTIONS
range0_1 <- function(x,na.rm){(x-min(x, na.rm = T)+0.01)/(max(x, na.rm = T)-min(x, na.rm = T)+0.01)}
z_trans <- function(myVar, na.rm){(myVar - mean(myVar, na.rm = T)) / sd(myVar, na.rm = T)}
'%ni%' <- Negate('%in%')
  
#

# TIME SERIES (Run Heterogeneity_Stability_Data2.R before running this) ####

Tile_data_all <- Raw_community_tax %>%
  mutate(Date2=as.character(as.Date(Date, format = "%d/%m/%Y", origin = "1899-12-30"))) %>%
  left_join(Community_data %>% dplyr::select(Transect, Station, Tile.Type, "Date2"=Date, Emersion.Rate)) %>%
  left_join(Community_data %>% distinct(Season, "Date2"=Date, Emersion.Rate)) %>%
  dplyr::select(Transect, Station, Tile.Type, Emersion.Rate, Date2, Season, "Cover"=Total_cov) %>%
  distinct(Transect, Station, Tile.Type, Emersion.Rate, Date2, Season, Cover) %>%
  mutate(
    Treatment = factor(Tile.Type, 
                       levels=c("R","S"), labels=c("+ Heterogeneity", "- Heterogeneity"))
  ) %>% 
  mutate(Taxa_group="Total_cover")

Tile_data <- Raw_community_tax %>%
  mutate(Date2=as.character(as.Date(Date, format = "%d/%m/%Y", origin = "1899-12-30"))) %>%
  left_join(Community_data %>% dplyr::select(Transect, Station, Tile.Type, Emersion.Rate, "Date2"=Date)) %>%
  left_join(Community_data %>% distinct(Season, "Date2"=Date)) %>%
  filter(Taxa_group!="Bare Substrate" & Taxa_group!="Ac_invertebrates" & Taxa_group!="Other") %>%
  group_by(Transect, Station, Tile.Type, Emersion.Rate, Date2, Season, Taxa_group) %>%
  dplyr::summarise(Cover=sum(Taxa_group_date_sum)) %>%
  bind_rows(Community_data %>% 
              dplyr::select(Transect, Station, Tile.Type, Emersion.Rate, "Date2"=Date, Season, S, Pielou) %>%
              pivot_longer(
                cols = S:Pielou,
                names_to = "Taxa_group",
                values_to = "Cover"
              )
  ) %>%
  bind_rows(Tile_data_all) %>%
  mutate(Taxa_group = ifelse(Taxa_group=="Limpets" | Taxa_group=="Whelks", "Consumers", Taxa_group)) %>%
  mutate(
    Year = year(as.Date(Date2)),
    Taxa_group = factor(Taxa_group, 
                             levels=c("S", "Pielou", "Total_cover", "Barnacles1", "Barnacles2", "Eph_macroalgae", "Consumers"),
                             labels=c("Richness", "Pielou", "Total cover", "A. modestus", "Native barnacles", "Ephemeral macroalgae", "Consumers")),
    Treatment = factor(Tile.Type, 
                       levels=c("R","S"), labels=c("+ Heterogeneity", "- Heterogeneity"))
    )

Time_taxa <- Tile_data %>%
  group_by(Date2, Season, Treatment, Taxa_group) %>%
  dplyr::summarise(
    Cover_mean = mean(Cover),
    Cover_SE = sd(Cover)/sqrt(n())
  )

#Date_model <- lmerTest::lmer(Cover ~ Date2 * Tile.Type + (1|Transect/Station), Tile_data %>% filter(Taxa_group=="Native barnacles")); plot(simulateResiduals(Date_model)); summary(Date_model); plot(allEffects(Date_model))

ggplot() + 
  geom_point(Tile_data, mapping=aes(x=Date2, y=Cover, col=Season, shape=Treatment), alpha=0.4) +
  geom_errorbar(Time_taxa, mapping=aes(x=Date2, ymin=Cover_mean-(Cover_SE*1.96), ymax=Cover_mean+(Cover_SE*1.96)), col="black", width=0.01, alpha=0.5) +
  geom_point(Time_taxa, mapping=aes(x=Date2, y=Cover_mean, shape=Treatment), col="black") +
  #geom_smooth(Time_taxa, mapping=aes(x=Date2, y=Cover_mean, col=Tile.Type), method="lm") +
  facet_wrap(~Taxa_group, scales = "free_y", nrow=2) +
  ylab("") + xlab("Date") +
  theme_bw()


ggplot() + 
  geom_point(Tile_data_all, mapping=aes(x=Date2, y=Cover, col=Season, shape=Treatment), alpha=0.4) +
  geom_errorbar(Time_taxa, mapping=aes(x=Date2, ymin=Cover_mean-(Cover_SE*1.96), ymax=Cover_mean+(Cover_SE*1.96)), col="black", width=0.01, alpha=0.5) +
  geom_point(Time_taxa, mapping=aes(x=Date2, y=Cover_mean, shape=Treatment), col="black") +
  #geom_smooth(Time_taxa, mapping=aes(x=Date2, y=Cover_mean, col=Tile.Type), method="lm") +
  #facet_wrap(~Taxa_group, scales = "free_y") +
  ylab("") + xlab("Date") +
  theme_bw()

#
# SENSITIVITY TEST AFTER REMOVING YEAR PER COMMUNITY COMPONENT ------

Tile_data_Year <- rbind(

# Mean across all years
Tile_data %>%
  group_by(Transect, Station, Tile.Type, Emersion.Rate, Taxa_group) %>%
  summarise(Cover = mean(Cover, na.rm = TRUE), .groups = "drop") %>%
  mutate(Year="All Years"),

# Mean excluding first year
Tile_data %>%
  filter(Date2 %ni% c("2019-10-27", "2020-01-21", "2020-04-12", "2020-07-20")) %>%
  group_by(Transect, Station, Tile.Type, Emersion.Rate, Taxa_group) %>%
  summarise(Cover = mean(Cover, na.rm = TRUE), .groups = "drop") %>%
  mutate(Year="No Year 1"),

# Mean excluding last year
Tile_data %>%
  filter(Date2 %ni% c("2021-10-09", "2022-01-20", "2022-04-01")) %>%
  group_by(Transect, Station, Tile.Type, Emersion.Rate, Taxa_group) %>%
  summarise(Cover = mean(Cover, na.rm = TRUE), .groups = "drop") %>%
  mutate(Year="No Year 2"),

# Mean excluding both first and last years
Tile_data %>%
  filter(Date2 %ni% c("2020-10-17", "2021-01-13", "2021-04-27", "2021-07-13")) %>%
  group_by(Transect, Station, Tile.Type, Emersion.Rate, Taxa_group) %>%
  summarise(Cover = mean(Cover, na.rm = TRUE), .groups = "drop") %>%
  mutate(Year="No Year 3")

) %>%
  mutate(
    Tile.Type = factor(Tile.Type, levels = c("R", "S"), labels = c("+ Heterogeneity", "- Heterogeneity")),
    Taxa_group = factor(Taxa_group, 
                        levels = c("Richness","Pielou","Total cover","A. modestus","Native barnacles","Consumers","Ephemeral macroalgae"),
                        labels = c("Richness","Pielou","Total cover","A. modestus","Native barnacles","Consumers","Ephemeral")
                        )
  ) %>%
  filter(Taxa_group!="Total cover")

Barnacles1_model<-lmerTest::lmer(Cover ~  
                                   Emersion.Rate + I(Emersion.Rate^2) +
                                   Tile.Type *
                                   Year +
                                   (1|Transect/Station)
                                 , 
                                 Tile_data_Year %>% filter(Taxa_group=="A. modestus")
); plot(simulateResiduals(Barnacles1_model)); summary(Barnacles1_model); plot(allEffects(Barnacles1_model))

Barnacles2_model<-lmerTest::lmer(Cover ~  
                                   Emersion.Rate + I(Emersion.Rate^2) +
                                   Emersion.Rate * Tile.Type *
                                   Year + 
                                   (1|Transect/Station)
                                 , 
                                 Tile_data_Year %>% filter(Taxa_group=="Native barnacles")
); plot(simulateResiduals(Barnacles2_model)); summary(Barnacles2_model); plot(allEffects(Barnacles2_model))

Eph_macroalgae_model<-lmerTest::lmer(Cover ~  
                                       Emersion.Rate * I(Emersion.Rate^2) *
                                       Tile.Type *
                                       Year + 
                                       (1|Transect/Station)
                                     , 
                                     Tile_data_Year %>% filter(Taxa_group=="Ephemeral")
); plot(simulateResiduals(Eph_macroalgae_model)); summary(Eph_macroalgae_model); plot(allEffects(Eph_macroalgae_model))

Consumer_model<-lmerTest::lmer(Cover ~  
                                 Emersion.Rate + I(Emersion.Rate^2) +
                                 Emersion.Rate * Tile.Type *
                                 Year + 
                                 (1|Transect/Station)
                               , 
                               Tile_data_Year %>% filter(Taxa_group=="Consumers")
); plot(simulateResiduals(Consumer_model)); summary(Consumer_model); plot(allEffects(Consumer_model))

S_model<-lmerTest::lmer(Cover ~  
                          Emersion.Rate + I(Emersion.Rate^2) +
                          Emersion.Rate * Tile.Type *
                          Year + 
                          (1|Transect/Station)
                        , 
                        Tile_data_Year %>% filter(Taxa_group=="Richness")
); plot(simulateResiduals(S_model)); summary(S_model); plot(allEffects(S_model))

Pielou_model<-lmer(Cover ~  
                     #Emersion.Rate + I(Emersion.Rate^2) + # not significant, model did not converge
                     Emersion.Rate * Tile.Type *
                     Year + 
                     (1|Transect/Station) # when including random effects, results do not change qualitatively. Since ggpredict can't compute confidence intervals for this model when including random effects, they had to be removed.
                   , 
                   Tile_data_Year %>% filter(Taxa_group=="Pielou")
); plot(simulateResiduals(Pielou_model)); summary(Pielou_model); plot(allEffects(Pielou_model))

# Helper to extract predictions and format
predict_on_observed <- function(model, data_subset, taxa_group_label, log_transform = FALSE) {
  # Predict with standard errors
  preds <- predict(model, newdata = data_subset, se.fit = TRUE, re.form = NA)
  
  # Create dataframe with predictions and 95% CIs
  out <- data_subset %>%
    mutate(
      fit = preds$fit,
      se = preds$se.fit,
      predicted = if (log_transform) exp(fit) - 1 else fit,
      conf.low = if (log_transform) exp(fit - 1.96 * se) - 1 else fit - 1.96 * se,
      conf.high = if (log_transform) exp(fit + 1.96 * se) - 1 else fit + 1.96 * se,
      Taxa_group = taxa_group_label
    )
  
  return(out)
}

# Create list of prediction datasets per taxon
predictions_list <- list(
  predict_on_observed(Barnacles1_model, Tile_data_Year %>% filter(Taxa_group == "A. modestus"), "A. modestus", TRUE),
  predict_on_observed(Barnacles2_model, Tile_data_Year %>% filter(Taxa_group == "Native barnacles"), "Native barnacles", TRUE),
  predict_on_observed(Eph_macroalgae_model, Tile_data_Year %>% filter(Taxa_group == "Ephemeral"), "Ephemeral", TRUE),
  predict_on_observed(Consumer_model, Tile_data_Year %>% filter(Taxa_group == "Consumers"), "Consumers"),
  predict_on_observed(S_model, Tile_data_Year %>% filter(Taxa_group == "Richness"), "Richness"),
  predict_on_observed(Pielou_model, Tile_data_Year %>% filter(Taxa_group == "Pielou"), "Pielou")
)

# Combine all predictions
all_preds <- bind_rows(predictions_list)

ggplot(all_preds, aes(x = Emersion.Rate, y = predicted, color = Year)) +
  # Raw data points from Tile_data_Year
  geom_point(
    data = Tile_data_Year,
    aes(x = Emersion.Rate, y = Cover, col = Year),
    alpha = 0.4, shape = 16,
    position = position_jitter(width = 0.05, height = 0)
  ) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high, fill = Year),
    alpha = 0.2, color = NA
  ) +
  # Prediction lines
  geom_line(size = 1) +
  facet_grid(Taxa_group ~ Tile.Type, scales = "free_y") +
  labs(
    x = "Emersion Ratio",
    y = "Observations / Predictions",
    color = "Year group",
    fill = "Year group"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top"
  )

#
# CORRELATE STABILITY PER TILE AFTER EXCLUDING 1ST, 2ND AND 3RD YEARS ####

Get_stability_data <- function(dataset, Selected_year) {
  
  Stability_out<- dataset %>%
    
    mutate(Date_exp = paste(Season, Year)) %>%
    
    mutate(
      Exp_year = ifelse(Date_exp=="Autumn 2019" | Date_exp=="Winter 2020" | Date_exp=="Spring 2020" | Date_exp=="Summer 2020", "Year 1",
                        ifelse(Date_exp=="Autumn 2020" | Date_exp=="Winter 2021" | Date_exp=="Spring 2021" | Date_exp=="Summer 2021", "Year 2",
                               ifelse(Date_exp=="Autumn 2021" | Date_exp=="Winter 2022" | Date_exp=="Spring 2022" | Date_exp=="Summer 2022", "Year 3", NA)))
    ) %>%
    
    filter(Exp_year != Selected_year) %>%
    
    #obtain residuals to calculate detrended alpha cover
    mutate(
      res = residuals(lm(Total_cov ~ Date2))
    ) %>%
    
    # obtain metrics per tile
    group_by(Transect, Station, Tile.Type, Altitude, Emersion.Rate, AltGroup, EmGroup) %>%
    summarise(
      S_mean = mean(S), #mean richness
      Pielou_mean = mean(Pielou), #mean Pielou based on Shannon (Hill_ev) and richness (log(S))
      Cover_stability_alpha = mean(Total_cov)/sd(Total_cov), #alpga stability
      Cover_stability_alpha_res = mean(Total_cov)/sd(res), #alpha stability using residuals
    )

return(Stability_out)

}

Stability_years <- Stability_data %>%
  left_join(Get_stability_data(Community_data, "Year 1") %>% dplyr::select(Transect, Station, Tile.Type, "Stability_minus_year1" = Cover_stability_alpha)) %>%
  left_join(Get_stability_data(Community_data, "Year 2") %>% dplyr::select(Transect, Station, Tile.Type, "Stability_minus_year2" = Cover_stability_alpha)) %>%
  left_join(Get_stability_data(Community_data, "Year 3") %>% dplyr::select(Transect, Station, Tile.Type, "Stability_minus_year3" = Cover_stability_alpha)) %>%
  mutate(Stability_rep = Cover_stability_alpha) %>%
  dplyr::select(Transect, Station, Tile.Type, Stability_minus_year1, Stability_minus_year2, Stability_minus_year3, Stability_rep, Cover_stability_alpha) %>%
  pivot_longer(
    cols = Stability_minus_year1:Stability_minus_year3,
    names_to = 'Stability_type',
    values_to = 'Stability_value'
  ) %>%
  mutate(
    Stability_type = factor(Stability_type, 
                                 levels=c("Stability_minus_year1","Stability_minus_year2","Stability_minus_year3"),
                                 labels=c("Stability minus year 1","Stability minus year 2","Stability minus year 3")),
    Treatment = factor(Tile.Type, 
                       levels=c("R","S"), labels=c("+ Heterogeneity", "- Heterogeneity"))
    )

ggplot(Stability_years) +
  geom_point(aes(x = Cover_stability_alpha, y = Stability_value, col=Stability_type, shape=Treatment), alpha=0.6) +
  geom_line(aes(x = Cover_stability_alpha, y = Stability_rep), col="black", alpha=0.6, linetype="dashed") +
  #geom_smooth(aes(x = Cover_stability_alpha, y = Stability_value, col=Stability_type), method="lm") +
  facet_wrap(~Stability_type, nrow=1, scales = "free_y") +
  ylab("Modified temporal stability") +
  xlab("Temporal stability") +
  ylim(0,30)+
  guides(col="none") +
  theme_bw()

TS_1 <- Stability_years %>% filter(Stability_type=="Stability minus year 1") %>% dplyr::select(Stability_value)
TS_2 <- Stability_years %>% filter(Stability_type=="Stability minus year 2") %>% dplyr::select(Stability_value)
TS_3 <- Stability_years %>% filter(Stability_type=="Stability minus year 3") %>% dplyr::select(Stability_value)
TS_all <- Stability_years %>% filter(Stability_type=="Stability minus year 1") %>% dplyr::select(Cover_stability_alpha)

cor(TS_1, TS_all, method="spearman")
cor(TS_2, TS_all, method="spearman")
cor(TS_3, TS_all, method="spearman")

Stability_years_res <- Stability_data %>%
  left_join(Get_stability_data(Community_data, "Year 1") %>% dplyr::select(Transect, Station, Tile.Type, "Stability_minus_year1" = Cover_stability_alpha_res)) %>%
  left_join(Get_stability_data(Community_data, "Year 2") %>% dplyr::select(Transect, Station, Tile.Type, "Stability_minus_year2" = Cover_stability_alpha_res)) %>%
  left_join(Get_stability_data(Community_data, "Year 3") %>% dplyr::select(Transect, Station, Tile.Type, "Stability_minus_year3" = Cover_stability_alpha_res)) %>%
  mutate(Stability_rep = Cover_stability_alpha_res) %>%
  dplyr::select(Transect, Station, Tile.Type, Stability_minus_year1, Stability_minus_year2, Stability_minus_year3, Stability_rep, Cover_stability_alpha_res) %>%
  pivot_longer(
    cols = Stability_minus_year1:Stability_minus_year3,
    names_to = 'Stability_type',
    values_to = 'Stability_value'
  ) %>%
  mutate(
    Stability_type = factor(Stability_type, 
                                 levels=c("Stability_minus_year1","Stability_minus_year2","Stability_minus_year3"),
                                 labels=c("Stability minus year 1","Stability minus year 2","Stability minus year 3")),
    Treatment = factor(Tile.Type, 
                       levels=c("R","S"), labels=c("+ Heterogeneity", "- Heterogeneity"))
    )

ggplot(Stability_years_res) +
  geom_point(aes(x = Cover_stability_alpha_res, y = Stability_value, col=Stability_type, shape=Treatment), alpha=0.6) +
  geom_line(aes(x = Cover_stability_alpha_res, y = Stability_rep), col="black", alpha=0.6, linetype="dashed") +
  #geom_smooth(aes(x = Cover_stability_alpha_res, y = Stability_value, col=Stability_type), method="lm") +
  facet_wrap(~Stability_type, nrow=1, scales = "free_y") +
  ylab("Modified temporal stability") +
  xlab("Temporal stability") +
  guides(col="none") +
  ylim(0,30)+
  theme_bw()

TS_1_res <- Stability_years_res %>% filter(Stability_type=="Stability minus year 1") %>% dplyr::select(Stability_value)
TS_2_res <- Stability_years_res %>% filter(Stability_type=="Stability minus year 2") %>% dplyr::select(Stability_value)
TS_3_res <- Stability_years_res %>% filter(Stability_type=="Stability minus year 3") %>% dplyr::select(Stability_value)
TS_all_res <- Stability_years_res %>% filter(Stability_type=="Stability minus year 1") %>% dplyr::select(Cover_stability_alpha_res)

cor(TS_1_res, TS_all_res, method="spearman")
cor(TS_2_res, TS_all_res, method="spearman")
cor(TS_3_res, TS_all_res, method="spearman")

#
# CORRELATE TRENDED AND DETRENDED STABILITY DATA ####

Stability_data2 <- Stability_data %>%
  mutate(
    Treatment = factor(Tile.Type, 
                       levels=c("R","S"), labels=c("+ Heterogeneity", "- Heterogeneity"))
  )

ggplot(Stability_data2) + 
  geom_line(aes(x=Cover_stability_alpha, y=Cover_stability_alpha, group=Treatment), col="black", linetype="dashed", alpha=0.5) +
  geom_point(aes(x=Cover_stability_alpha, y=Cover_stability_alpha_res, col=Treatment)) +
  #geom_smooth(aes(x=Cover_stability_alpha, y=Cover_stability_alpha_res, col=Tile.Type), method="lm") +
  ylab("Detrended temporal stability") +
  xlab("Temporal stability") +
  theme_bw()

cor(Stability_data$Cover_stability_alpha, Stability_data$Cover_stability_alpha_res, method="spearman")

#
# SENSITIVITY TEST FOR DETRENDED TEMPORAL STABILTIY LMER ####

Cov_model<-lmerTest::lmer(log(Cover_stability_alpha_res) ~  
                            Emersion.Rate * I(Emersion.Rate^2) +
                            Emersion.Rate * Tile.Type +
                            (1|Transect/Station)
                          , 
                          Stability_data
); plot(simulateResiduals(Cov_model)); summary(Cov_model); AIC(Cov_model); plot(allEffects(Cov_model))

Cov_model<-lmerTest::lmer(log(Cover_stability_alpha) ~  
                            Emersion.Rate * I(Emersion.Rate^2) +
                            Emersion.Rate * Tile.Type +
                            (1|Transect/Station)
                          , 
                          Stability_data
); plot(simulateResiduals(Cov_model)); summary(Cov_model); AIC(Cov_model); plot(allEffects(Cov_model))

#
# CONSUMER LAGGED CORRELATION ####

Raw_community_groups <- Raw_community_tax %>%
  group_by(Transect, Station, Tile.Type, Date, Taxa_group) %>%
  dplyr::summarise(Cover=sum(cover_sum)) %>%
  group_by(Transect, Station, Tile.Type, Taxa_group) %>%
  mutate(Cover_mean = mean(Cover)) %>%
  ungroup() %>%
  mutate(Date2 = z_trans(as.numeric(Date), na.rm=T))
  
Consumers_acf <- Raw_community_groups %>%
  filter(Taxa_group=="Consumers") %>%
  group_by(Transect, Station, Tile.Type) %>%
  summarise(ac = list(acf(Cover, lag.max = 10))) 

Consumers_lag <- Raw_community_groups %>%
  filter(Taxa_group=="Consumers") %>%
  mutate(Cover_res = Cover - Cover_mean) %>%
  mutate(Cover_res_w = (Cover - Cover_mean)/Cover_mean) %>%
  mutate(Tile=paste(Transect, Station)) %>%
  left_join(Community_data %>% dplyr::select(Transect, Station, Tile.Type, EmGroup)) %>%
  mutate(
    EmGroup=factor(EmGroup, 
                   levels=c("[0.1934,0.296]","[0.1394,0.193)", "[0.0664,0.139)"),
                   labels=c("High shore", "Mid shore", "Low shore")
    ),
    Treatment=factor(Tile.Type,
                     levels=c("R","S"),
                     labels=c("+ Heterogeneity", "- Heterogeneity")
    )
  ) %>%
  distinct()


ggplot() +
  #geom_hline(yintercept = 0, linetype="dashed", alpha=0.3) +
  geom_point(Consumers_lag, mapping=aes(x=Cover_mean, y=Cover, col=Treatment)) +
  geom_line(Consumers_lag, mapping=aes(x=Cover_mean, y=Cover_mean), col="black", linetype="dashed") +
  #geom_point(Consumers_lag, mapping=aes(x=Cover_mean, y=Cover, col=Treatment), alpha=0.3) +
  #geom_line(Consumers_lag, mapping=aes(x=Cover_mean, y=Cover_mean), col="black", linetype="dashed") +
  #facet_wrap(~Treatment*EmGroup, scales="free")+
  ylab("Cover") + xlab("Mean cover") +
  theme_bw()

cor(Consumers_lag$Cover_mean, Consumers_lag$Cover)

ggplot() +
  geom_hline(yintercept = 0, linetype="dashed", alpha=0.3) +
  #geom_line(Consumers_lag, mapping=aes(x=Cover_mean, y=Cover_mean, col=Tile.Type), linetype="dashed") +
  #geom_point(Consumers_lag, mapping=aes(x=Cover_mean, y=Cover, col=Tile.Type)) +
  geom_point(Consumers_lag, mapping=aes(x=Cover_mean, y=Cover_res, col=Treatment)) +
  facet_wrap(~Treatment*EmGroup)+
  ylab("Residual (Cover - Mean cover)") + xlab("Mean cover") +
  theme_bw()

ggplot() +
  geom_hline(yintercept = 0, linetype="dashed", alpha=0.3) +
  #geom_line(Consumers_lag, mapping=aes(x=Cover_mean, y=Cover_mean, col=Tile.Type), linetype="dashed") +
  #geom_point(Consumers_lag, mapping=aes(x=Cover_mean, y=Cover, col=Tile.Type)) +
  geom_point(Consumers_lag, mapping=aes(x=Cover_mean, y=Cover_res_w, col=Treatment), alpha=0.3) +
  facet_wrap(~Treatment*EmGroup, scales="free_x")+
  ylab("Residual ([Cover - Mean cover]/Mean cover)") + xlab("Mean cover") +
  theme_bw() + theme(legend.position = "none")


# %>%
  arrange(Transect, Station, Tile.Type, Date) %>% 
  group_by(Transect, Station, Tile.Type) %>% 
  mutate(
    lag_Cover1 = lag(Cover, 1L),
    lag_Cover2 = lag(Cover, 2L),
    lag_Cover3 = lag(Cover, 3L),
    lag_Cover4 = lag(Cover, 4L),
    lag_Cover5 = lag(Cover, 5L),
    lag_Cover6 = lag(Cover, 6L),
    lag_Cover7 = lag(Cover, 7L),
    lag_Cover8 = lag(Cover, 8L),
    lag_Cover9 = lag(Cover, 9L),
    lag_Cover10 = lag(Cover, 10L)
    )

Consumers_lag_plot <- Consumers_lag %>%
  pivot_longer(
    cols=lag_Cover1:lag_Cover10,
    names_to="Lag",
    values_to="Lag_Cover"
  ) %>%
  mutate(Lag=factor(Lag,
                    levels=c("lag_Cover1","lag_Cover2","lag_Cover3","lag_Cover4","lag_Cover5","lag_Cover6","lag_Cover7","lag_Cover8","lag_Cover9","lag_Cover10")
                    ))


ggplot(Consumers_lag_plot, aes(x=Cover, y=Lag_Cover, col=Lag, shape=Tile.Type)) +
  geom_line(Consumers_lag_plot, mapping=aes(x=Cover, y=Cover), col="black", linetype="dashed")+
  geom_point()+
  facet_wrap(~Lag) +
  theme_bw()

Get_cor_values <- function(dataset){
  
  Consumer_out <- data.frame(
    "Time_lag"=c(1:8)
    ,
    "Spearman_values"=c(
      cor(dataset$Cover, dataset$lag_Cover1, method="spearman", use="pairwise.complete.obs"),
      cor(dataset$Cover, dataset$lag_Cover2, method="spearman", use="pairwise.complete.obs"),
      cor(dataset$Cover, dataset$lag_Cover3, method="spearman", use="pairwise.complete.obs"),
      cor(dataset$Cover, dataset$lag_Cover4, method="spearman", use="pairwise.complete.obs"),
      cor(dataset$Cover, dataset$lag_Cover5, method="spearman", use="pairwise.complete.obs"),
      cor(dataset$Cover, dataset$lag_Cover6, method="spearman", use="pairwise.complete.obs"),
      cor(dataset$Cover, dataset$lag_Cover7, method="spearman", use="pairwise.complete.obs"),
      cor(dataset$Cover, dataset$lag_Cover8, method="spearman", use="pairwise.complete.obs")
      )
  )
  
  return(Consumer_out)
  
}

Consumer_correlations <- Consumers_lag %>%
  mutate(Tile=factor(paste(Transect, Station, Tile.Type))) %>%
  split(.$Tile) %>% # separate dataset into multiple dataset - each dataset per tile
  lapply(function(x) {Get_cor_values(x)}) %>%
  Map(cbind, ., Tile=names(.)) %>% # add column indicating which Tile and date each dataframe corresponds to
  map_dfr(., as.list ) %>% # join dataframes within list into one single dataframe
  separate(Tile, into=c("Transect", "Station", "Tile.Type"), remove=FALSE) %>%
  mutate(
    Transect=as.integer(Transect),
    Station=as.integer(Station)
  ) %>%
  left_join(Community_data %>% dplyr::select(Transect, Station, Tile.Type, EmGroup)) %>%
  distinct() %>%
  mutate(
    EmGroup=factor(EmGroup, 
                   levels=c("[0.1934,0.296]","[0.1394,0.193)", "[0.0664,0.139)"),
                   labels=c("High shore", "Mid shore", "Low shore")
                   ),
    Treatment=factor(Tile.Type,
                     levels=c("R","S"),
                     labels=c("+ Heterogeneity", "- Heterogeneity")
                     )
  )

Consumer_correlations_mean <- Consumer_correlations %>%
  group_by(EmGroup, Treatment, Time_lag) %>%
  dplyr::summarise(
    Spearman_value = mean(Spearman_values, na.rm=T),
    Spearman_Se = sd(Spearman_values, na.rm=T)/sqrt(sum(!is.na(Spearman_values)))
  )

ggplot() +
  geom_hline(yintercept=0, linetype="dotted", alpha=0.3) +
  geom_point(Consumer_correlations, mapping=aes(x=Time_lag, y=Spearman_values, col=Treatment, group=Tile), alpha=0.3) +
  geom_point(Consumer_correlations_mean, mapping=aes(x=Time_lag, y=Spearman_value, shape=Treatment), col="black") +
  geom_errorbar(Consumer_correlations_mean, mapping=aes(x=Time_lag, ymin=Spearman_value-Spearman_Se, ymax=Spearman_value+Spearman_Se), col="black", width=0.1) +
  facet_wrap(~Treatment*EmGroup)+
  ylab("Spearman correlation") +
  xlab("Time Lag") +
  theme_bw() +
  theme(legend.position = "none")

High_Het <- Consumer_correlations %>% filter(EmGroup=="High shore" & Treatment=="+ Heterogeneity")
Mid_Het <- Consumer_correlations %>% filter(EmGroup=="Mid shore" & Treatment=="+ Heterogeneity")
Low_Het <- Consumer_correlations %>% filter(EmGroup=="Low shore" & Treatment=="+ Heterogeneity")

High_NoHet <- Consumer_correlations %>% filter(EmGroup=="High shore" & Treatment=="- Heterogeneity")
Mid_NoHet <- Consumer_correlations %>% filter(EmGroup=="Mid shore" & Treatment=="- Heterogeneity")
Low_NoHet <- Consumer_correlations %>% filter(EmGroup=="Low shore" & Treatment=="- Heterogeneity")

cor(High_Het$Time_lag, High_Het$Spearman_values, method="spearman", use="pairwise.complete.obs")
cor(High_NoHet$Time_lag, High_NoHet$Spearman_values, method="spearman", use="pairwise.complete.obs")

cor(Mid_Het$Time_lag, Mid_Het$Spearman_values, method="spearman", use="pairwise.complete.obs")
cor(Mid_NoHet$Time_lag, Mid_NoHet$Spearman_values, method="spearman", use="pairwise.complete.obs")

cor(Low_Het$Time_lag, Low_Het$Spearman_values, method="spearman", use="pairwise.complete.obs")
cor(Low_NoHet$Time_lag, Low_NoHet$Spearman_values, method="spearman", use="pairwise.complete.obs")


#
# CONSUMER LAGGED CORRELATION - ABUNDANCE ####

Consumer_Abundance<-read.csv("/Users/jordisola/Library/Mobile Documents/com~apple~CloudDocs/Documents/PhD_thesis/Habitat_stability/Data_analysis/CommunityAnalysis/Species_psd_raw_separated.csv") %>%
  dplyr::select(Transect, Station, Tile.Type, Date, Taxa) %>%
  filter(Taxa=="PV" | Taxa=="PV+H" | Taxa=="PD" | Taxa=="PSP" | Taxa=="NL" | Taxa=="OA") %>%
  group_by(Transect,Station,Tile.Type,Date) %>%
  summarise(Abundance = n()) %>%
  group_by(Transect,Station,Tile.Type) %>%
  mutate(Abundance_mean = mean(Abundance)) %>%
  left_join(Stability_data %>% ungroup() %>% dplyr::select(Transect,Station,Tile.Type,"Cover_mean"=Consumer_mean)) %>%
  ungroup() %>%
  mutate(
    Abundance_01 = range0_1(Abundance, na.rm=T),
    Abundance_mean_01 = range0_1(Abundance_mean, na.rm=T),
    Cover_mean_01 = range0_1(Cover_mean, na.rm=T)
  )

ggplot() +
  geom_point(Consumer_Abundance, mapping=aes(x=Abundance_mean, y=Abundance, col=Tile.Type)) +
  geom_line(Consumer_Abundance, mapping=aes(x=Abundance_mean, y=Abundance_mean), col="black", linetype="dashed") +
  ylab("Abundance") + xlab("Mean abundance") +
  theme_bw()

cor(Consumer_Abundance$Abundance_mean_01, Consumer_Abundance$Abundance_01)

#