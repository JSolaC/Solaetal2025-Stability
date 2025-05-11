# LOAD LIBRARIES
library(tidyverse)
library(Hmisc)
library(labdsv)
library(betapart)
library(reshape2)

#
#### LOAD FUNCTIONS ####
range0_1 <- function(x,na.rm){(x-min(x, na.rm = T)+0.01)/(max(x, na.rm = T)-min(x, na.rm = T)+0.01)}
z_trans <- function(myVar, na.rm){(myVar - mean(myVar, na.rm = T)) / sd(myVar, na.rm = T)}
'%ni%' <- Negate('%in%')
sd_true<-function(x){sd(x, na.rm=TRUE)}
sum_true<-function(x){sum(x, na.rm=TRUE)}
Dissimilarity_to_similarity<-function(x){ ( x*(-1) ) + 1}

#
#### LOAD DATASETS ####

# SPECIES-LEVEL DATA (per taxa, tile and time point)
Raw_community_tax<-read.csv("Data/Datasets/Raw_community_tax.csv") %>%
  
  mutate(
    Date = as.Date(Date, format = "%d/%m/%Y"),
    Taxa_group = ifelse(Taxa=="AM" ,"Barnacles1",
                         ifelse(Taxa=="SB" | Taxa=="BC" | Taxa=="PER"| Taxa=="CHM", "Barnacles2",
                                ifelse(Taxa=="PSP" | Taxa=="PV" | Taxa=="PD", "Limpets",
                                ifelse(Taxa=="NL", "Whelks",
                                       ifelse(Taxa=="ME" | Taxa=="SA" | Taxa=="SPT", "Ac_invertebrates",
                                              ifelse(Taxa=="US" | Taxa=="US1" | Taxa=="US2", "Eph_macroalgae", 
                                                     ifelse(Taxa!="BS","Other", "Bare Substrate")))))))
  ) %>%
  
  # create variables to later sort data into groups 
  group_by(Transect, Station, Tile.Type, Taxa) %>%
  mutate(
    Tile=paste(Transect, Station, Tile.Type),
    Taxa_Tile=paste(Taxa, Transect, Station, Tile.Type)
  ) %>%
  
  #get values per date - group taxa
  group_by(Transect, Station, Tile.Type, Date, Taxa_group) %>%
  mutate(Taxa_group_date_sum=sum(cover_sum)) %>%
  
  # get values per data - overall cover
  group_by(Transect, Station, Tile.Type, Date) %>%
  mutate(Total_cov = sum(cover_sum)) 

# COMMUNITY-LEVEL DATA (per tile and time point)
Community_data<-read.csv("Data/Datasets/Community_data.csv") %>%
  
  dplyr::select(-Emersion.Rate) %>%
  
  left_join(read.csv("Data/Datasets/ShoreHeight.csv") %>%
              dplyr::select(Altitude, Emersion.Rate), by="Altitude") %>%
  
  mutate(
    Pielou = Hill_ev/log(S),
    Total_cov = as.numeric(Total_cov),
    Date2 = z_trans(as.numeric(as.Date(Date, format = "%Y-%m-%d")), na.rm=T),
    Tile.Type2 = ifelse(Tile.Type=="R",1,0),
    EmGroup = cut2(Emersion.Rate, g=3)
  )

#
#### OBTAIN DATA - STABILITY METRICS ####

# PREPARE DATA: Obtain main stats descriptors per species
Species_metrics <- Raw_community_tax  %>%
  
  filter(Taxa!="BS") %>%
  
  # Per taxa, obtain the mean value and SD per species
  group_by(Transect, Station, Tile.Type, Taxa) %>%
  mutate(
    Species_mean = mean(cover_sum), # mean cover over time
    Species_sd = sd(cover_sum), # variability in cover over time
    n_timepoints2 = n() # how many time points?
  )%>% 
  distinct(Transect, Station, Tile.Type, Taxa, .keep_all=TRUE) %>% # this is done so as to a) summmarise and b) keep all columns
  
  # Over time, compute the mean and sd across species - CHECK THIS WHOLE SECTION!
  group_by(Transect, Station, Tile.Type) %>%
  mutate(
    Species_mean_sum = sum(Species_mean), # needed to compute Population Stability (as per Species Stability) and for Species Asynchrony
    Species_sd_sum = sum(Species_sd[n_timepoints2>2]), # needed to compute Species Stability, Population Stability and Statistical Averaging
    Species_sd2_sum = sqrt(sum((Species_sd[n_timepoints2>2]^2))), # needed to compute Statistical Averaging
    Total_cov_sd=sd(Total_cov),
    Total_cov_mean=mean(Total_cov) 
    ) %>% # needed to compute Population Stability 
  
  group_by(Transect, Station, Tile.Type, Taxa_group) %>%
  mutate(
    Taxa_group_mean_time = mean(Taxa_group_date_sum)
    ) %>%
  
  # Per species, but using metrics calculated across species as well
  ungroup() %>%
  mutate(
    Species_weighted_cover = (Species_mean/Species_mean_sum), # needed to compute species asynchrony in Dataset 1
  ) %>%
  
  dplyr::select(Transect, Station, Tile.Type, Tile, Taxa, Taxa_group,
                Total_cov_mean, Total_cov_sd, Taxa_group_mean_time,
                n_timepoints2, Species_mean, Species_mean_sum, Species_sd_sum, Species_sd2_sum, Species_weighted_cover)

# DATASET 1: Species asynchrony following CRAVEN ET AL 2018 and Blüthgen ET AL 2016
Time_asynchrony <- Raw_community_tax %>%
  
  filter(Taxa!="BS") %>%
  
  # To allow for reliable corr values, spp. with less than 3 observations cannot be accounted for
  group_by(Transect, Station, Tile.Type, Taxa) %>%
  filter(length(unique(Date))>2) %>% 
  ungroup() %>%
  
  # Obtain a list with tile dataset per entry (incl. all dates per tile within each dataset)
  dplyr::select(Tile, Date, Taxa, cover_sum) %>% # select the variables that should be included within each dataset
  mutate(
    Tile = factor(Tile),
    Date=as.character(Date)
    ) %>% # for the split function to work - Tile needs to be of the type 'factor'
  split(.$Tile) %>% # separate dataset into multiple dataset - each dataset per tile
  
  # calculate correlation coefficients per tile - accounting for all dates within each tile
  lapply(function(x) {as.data.frame(x)[,-1]}) %>% #remove column re. tile
  lapply(function(x) {matrify(x)}) %>% # obtain matrix
  lapply(function(x) {replace(x,x==0, NA)}) %>% # replace 0 with NA, since that is data that I do not have (cover = 0 means that taxa were not found there)
  lapply(function(x) {cor(x, use="pairwise.complete.obs")}) %>% # obtain correlation coefficients, while omitting NAs. 
  # Warning messages refer to correlations when only one value is repeated over time for one of the two taxa. Since they are unvariable and they would not contribute to asynchrony ('value = 0'), then it is ok since it won't compute for corr cums but can still account for them bc I count number of rows to quantify number of corrs / number of taxa correlated in each instance
  lapply(function(x) {ifelse(x==diag(x),NA,x)}) %>% # remove elements in the diagonal (autocorrelations)
  lapply(function(x) {as.data.frame(x)}) %>% # convert to dataframe
  lapply(function(x) {transform(x, sd_cor=apply(x, 1, sd_true))}) %>% #calculate sd per row
  # sum correlation values per row - I did subtract sd_cor
  lapply(function(x) {transform(x, sum_cor=apply(x, 1, sum_true)-sd_cor)}) %>% # sum correlation values per row
  lapply(function(x) {transform(x, n_cor=nrow(x))}) %>% # how many correlations are there?
  lapply(function(x) {transform(x, Taxa=rownames(x))}) %>% # introduce the names of the rows: taxa
  Map(cbind, ., Tile=names(.)) %>% # add column indicating which Tile and date each dataframe corresponds to
  map_dfr(., as.list ) %>% # join dataframes within list into one single dataframe
  dplyr::select(Tile,Taxa,n_cor,sum_cor,sd_cor) %>% # select variables needed
  separate(Tile, into=c("Transect", "Station", "Tile.Type")) %>% # Obtain Plot ID (Tile ID) variables from 'summary' variable 'Tile'
  filter(!is.na(sum_cor)) %>% # There is a case (CHM in tile 5 3 S) that returns NA since it has the same value 0.02 for all three dates when it was found. It needs to be removed or else Asynchrony turns NA at the end for that tile.
  
  left_join( # Join with species_metrics to be able to weight individual asynchrony metrics per species following weighted cover per species (following Blüthgen et al 2016)
    Species_metrics %>% 
      select(Transect, Station, Tile.Type, Taxa, Species_weighted_cover) %>% 
      distinct() %>%
      mutate_at(c(1,2), as.character)
  ) %>%
  
  # to avoid introducing NA, Inf or accounting for 'rare' species, only cases where at least three diff cover values per species are found
  group_by(Transect, Station, Tile.Type, Taxa) %>%
  # as theory predicts that asynchrony happens due to a) competition or b) asynchronous response to env conditions, we remove consumers from the dataset to avoid confounding competition with consumer effects
  # also, remove bare substrate again to make sure it is not in the dataset
  filter(Taxa!="NL" & Taxa!="PSP" & Taxa!="PV" & Taxa!="PD" & Taxa!="BS") %>%
  mutate(Asynchrony_species =  (sum(sum_cor)/n_cor) * Species_weighted_cover) %>%
  group_by(Transect,Station,Tile.Type) %>%
  summarise(Asynchrony_corr_time = sum(Asynchrony_species) * (-1))

# DATASET 2: Compute population stability metrics and statistical averaging
Pop_metrics <- Species_metrics %>%

#make sure to use stats summary for every metric - otherwise number of rows is not reduced!
group_by(Transect,Station,Tile.Type) %>%
  mutate(
    Species_stability = Species_mean_sum / Species_sd_sum,
    Pop_stab = Total_cov_mean/Species_sd_sum, # R2=0.98 in comparison with 'Species_stability'
    Statistical_averaging = Species_sd_sum/Species_sd2_sum,
    Asynchrony_ratio_time = (Total_cov_sd^2)/(Species_sd_sum^2),
    Whelk_mean=sum((Species_mean[Taxa_group=="Whelks"])),
    Limpet_mean=sum((Species_mean[Taxa_group=="Limpets"])),
    Consumer_mean=sum((Species_mean[Taxa_group=="Limpets"|Taxa_group=="Whelks"])),
    Barnacles1_mean=sum((Species_mean[Taxa_group=="Barnacles1"])),
    Barnacles2_mean=sum((Species_mean[Taxa_group=="Barnacles2"])),
    Ac_invertebrates_mean=sum((Species_mean[Taxa_group=="Ac_invertebrates"])),
    Eph_macroalgae_mean=sum((Species_mean[Taxa_group=="Eph_macroalgae"]))
  ) %>%
  dplyr::distinct(Tile, .keep_all=TRUE) %>%
  select(
    Transect, Station, Tile.Type,
    Consumer_mean, Limpet_mean, Whelk_mean, Ac_invertebrates_mean, Eph_macroalgae_mean, Barnacles1_mean, Barnacles2_mean,
    Species_stability, Pop_stab, Statistical_averaging, Asynchrony_ratio_time, n_timepoints2
  )

# DATASET 3: Compute mean composition stability across consecutive time points within a tile
Composition_out <- Raw_community_tax %>%
  
  # select variables and filter out bare substrate from dataset
  ungroup() %>%
  dplyr::select("Tile.Date", "Taxa", "cover_sum") %>%
  filter(Taxa!="BS") %>% # make sure you remove bare substrate
  
  # turn dataframe into a matrix
  as.data.frame() %>%
  matrify() %>%
  
  # calculate beta diversity metrics
  beta.pair.abund(.,index.family="bray") %>%
  pluck("beta.bray") %>% #select metric
  
  # calculate similarity from dissimilarity (as per Hildebrand et al. 2016)
  as.matrix() %>%
  as.data.frame() %>%
  mutate(across(.cols=everything(), Dissimilarity_to_similarity)) %>%
  
  #Obtain beta diversity metrics - comapring 
  list() %>% 
  lapply(function(x) {as.matrix(x)}) %>%
  lapply(function(x) {melt(x)}) %>%
  do.call(cbind,.) %>%
  rename("Var1" = 1,"Var2" = 2,"Bray" = 3) %>%
  as.data.frame() %>%
  filter(Var1 != Var2) %>%
  separate(Var1, into=c("Transect", "Station", "Tile.Type", "Date"), sep = " ") %>%
  separate(Var2, into=c("Transect2", "Station2", "Tile.Type2", "Date2"), sep = " ") %>%
  mutate(Date = as.Date(Date, format = "%d/%m/%Y")) %>%
  mutate(Date2 = as.Date(Date2, format = "%d/%m/%Y")) %>%
  
  #select tiles for time comparisons. Requisites: same tile over time, but only comparing among subsequent dates (e.g. dates 1 to 2, 2 to 3, etc.)
  #select subsequent dates separately to make sure this step is solid (crucial for this time comparison)
  #filter(Transect==Transect2 & Station==Station2 & Tile.Type==Tile.Type2) %>%
  filter(
    
    Transect==Transect2 & Station==Station2 & Tile.Type==Tile.Type2 &
    
      (
      (Date == "2019-10-27" & Date2 == "2020-01-21") |
      (Date == "2020-01-21" & Date2 == "2020-04-12") |
      (Date == "2020-04-12" & Date2 == "2020-07-20") |
      (Date == "2020-07-20" & Date2 == "2020-10-17") |
      (Date == "2020-10-17" & Date2 == "2021-01-13") |
      (Date == "2021-01-13" & Date2 == "2021-04-27") |
      (Date == "2021-04-27" & Date2 == "2021-07-13") |
      (Date == "2021-07-13" & Date2 == "2021-10-09") |
      (Date == "2021-10-09" & Date2 == "2022-01-20") |
      (Date == "2022-01-20" & Date2 == "2022-04-01")
      )
    
  ) %>%
  group_by(Transect, Station, Tile.Type) %>%
  summarise(
    Compositional_similarity = mean(Bray),
    n_composition = n() # safety check to see if the number of time points included influenced mean compositional stability (R2=0.08)
  )

# FINAL DATASET: Calculate community metrics and join all datasets
Stability_metrics<- Community_data %>%
  
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
  ) %>%
  
  
  left_join(Time_asynchrony %>%
              mutate_at(c(1,2), as.integer)) %>%
  left_join(Pop_metrics %>%
              mutate_at(c(1,2), as.integer)) %>%
  left_join(Composition_out %>%
              mutate_at(c(1,2), as.integer))

#
#### OUTPUT DATA ####

write.csv(Stability_metrics, "Data/Stability_metrics.csv")

#