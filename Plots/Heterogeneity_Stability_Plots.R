### LOAD LIBRARIES
library(tidyverse)
library(ggeffects)
library(ggforce)

### LOAD FUNCTIONS 

Plot_model<-function(input_model, input_dataset, Factor1, Factor2, Factor1_dataset,  Response_variable, y_title, x_title) {

mydf <- ggpredict(input_model, terms = c(Factor1, Factor2))
mydf2 <- input_dataset[,c(Factor1_dataset, Response_variable, Factor2)]
names(mydf2)<-c('x', 'predicted', 'group')
plot<-ggplot() +
  geom_vline(xintercept=0.274, linetype="dashed")+
  geom_vline(xintercept=0.383, linetype="dashed")+
  geom_ribbon(data=mydf, aes(x = x, ymax=conf.high, ymin=conf.low, fill=as.factor(group)), alpha=.2)+
  geom_point(data=mydf2, aes(x=x, y=predicted, colour=as.factor(group)), alpha=.5)+
  geom_line(data=mydf, aes(x = x, y = predicted, colour =as.factor(group))) +
  ylab(y_title)+
  xlab(x_title)+
  theme_bw() +
  theme(legend.position = "none")

return(plot)

}

#
#### SPECIES GROUPS ####

#4x3.5
Plot_model(Barnacles1_model, Stability_data, "Emersion.Rate [all]", "Tile.Type", "Emersion.Rate", "Barnacles1_mean", "Barnacles1_mean", "Emersion.Rate")
Plot_model(Barnacles2_model, Stability_data, "Emersion.Rate [all]", "Tile.Type", "Emersion.Rate", "Barnacles2_mean", "Barnacles2_mean", "Emersion.Rate")
#Plot_model(Ac_invertebrates_model, Stability_data, "Emersion.Rate [all]", "Tile.Type", "Emersion.Rate", "Ac_invertebrates_mean", "Ac_invertebrates_mean", "Emersion.Rate")
Plot_model(Consumer_model, Stability_data, "Emersion.Rate [all]", "Tile.Type", "Emersion.Rate", "Consumer_mean", "Consumer_mean", "Emersion.Rate")
Plot_model(Eph_macroalgae_model, Stability_data, "Emersion.Rate [all]", "Tile.Type", "Emersion.Rate", "Eph_macroalgae_mean", "Eph_macroalgae_mean", "Emersion.Rate")

# 4 x4 
Plot_model(S_model, Stability_data, "Emersion.Rate [all]", "Tile.Type", "Emersion.Rate", "S_mean", "S_mean", "Emersion.Rate")
Plot_model(Pielou_model, Stability_data, "Emersion.Rate [all]", "Tile.Type", "Emersion.Rate", "Pielou_mean", "Pielou_mean", "Emersion.Rate")

#
#### STABILITY DIMENSIONS ####

Plot_model(Pop_stab_model, Stability_data, "Emersion.Rate [all]", "Tile.Type", "Emersion.Rate", "Pop_stab", "Pop_stab", "Emersion.Rate")
Plot_model(Asynchrony_model, Stability_data, "Emersion.Rate [all]", "Tile.Type", "Emersion.Rate", "Asynchrony_corr_time", "Asynchrony_corr_time", "Emersion.Rate")
Plot_model(Stat_av_model, Stability_data, "Emersion.Rate [all]", "Tile.Type", "Emersion.Rate", "Statistical_averaging", "Statistical_averaging", "Emersion.Rate")
Plot_model(Comp_model, Stability_data, "Emersion.Rate [all]", "Tile.Type", "Emersion.Rate", "Compositional_similarity", "Compositional_similarity", "Emersion.Rate")
Plot_model(Cov_model, Stability_data, "Emersion.Rate [all]", "Tile.Type", "Emersion.Rate", "Cover_stability_alpha_res", "Cover_stability_alpha_res", "Emersion.Rate")

#

#### SEM PATHWAYS ####

SEM_coefs <- read.csv('/Users/jordisola/Library/Mobile Documents/com~apple~CloudDocs/Documents/PhD_thesis/Habitat_stability/Data_analysis/CommunityAnalysis/Publication/Models/SEM_multigroup_coefs.csv') %>%
  filter(Var.9!='')

SEM_one_coefs <- read.csv('/Users/jordisola/Library/Mobile Documents/com~apple~CloudDocs/Documents/PhD_thesis/Habitat_stability/Data_analysis/CommunityAnalysis/Publication/Models/SEM_one_coefs.csv') %>%
  filter(Var.9!='') %>%
  mutate(EmGroup="A")

Tile_coefs <- SEM_coefs %>%
  filter(Predictor!="Emersion.Rate_z")

Obtain_SEM_effects_unique<-function(input_dataset){
  
  SEM_loop<-input_dataset %>%
    ungroup() %>%
    dplyr::select(Tile_alt,Response,Predictor,Std.Estimate) %>%
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
    
    input_dataset<-input_dataset1<-SEM_loop %>%filter(Tile_alt==a)
    
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
  
  
  for(a in unique(SEM_loop$Tile_alt)){ 
    
    SEM_group<-SEM_loop %>%
      filter(Tile_alt==a)%>%
      Check_predictors
    
    SEM_out<-rbind(SEM_out, SEM_group)
    
  } 
  
  #Remove repeated pathways
  
  SEM_out2<-SEM_out %>%
    #remove repeated pathways
    distinct(Tile_alt,Response,Predictor,Predictor1,Predictor2,Predictor3,Predictor4,Predictor5,Predictor6,Predictor7,Predictor8,Predictor9,Predictor10,.keep_all=TRUE) %>%
    
    #remove nested pathways
    
    group_by(Tile_alt,Response,Predictor) %>%
    mutate(Predictor_n=ifelse(!is.na(Predictor),n(),NA))%>%
    
    group_by(Tile_alt,Response,Predictor,Predictor1) %>%
    mutate(Predictor1_n=ifelse(!is.na(Predictor1),n(),NA))%>%
    
    group_by(Tile_alt,Response,Predictor,Predictor1,Predictor2) %>%
    mutate(Predictor2_n=ifelse(!is.na(Predictor2),n(),NA))%>%
    
    group_by(Tile_alt,Response,Predictor,Predictor1,Predictor2,Predictor3) %>%
    mutate(Predictor3_n=ifelse(!is.na(Predictor3),n(),NA))%>%
    
    group_by(Tile_alt,Response,Predictor,Predictor1,Predictor2,Predictor3,Predictor4) %>%
    mutate(Predictor4_n=ifelse(!is.na(Predictor4),n(),NA))%>%
    
    group_by(Tile_alt,Response,Predictor,Predictor1,Predictor2,Predictor3,Predictor4,Predictor5) %>%
    mutate(Predictor5_n=ifelse(!is.na(Predictor5),n(),NA))%>%
    
    group_by(Tile_alt,Response,Predictor,Predictor1,Predictor2,Predictor3,Predictor4,Predictor5,Predictor6) %>%
    mutate(Predictor6_n=ifelse(!is.na(Predictor6),n(),NA)) %>%
    
    group_by(Tile_alt,Response,Predictor,Predictor1,Predictor2,Predictor3,Predictor4,Predictor5,Predictor6,Predictor7) %>%
    mutate(Predictor7_n=ifelse(!is.na(Predictor7),n(),NA)) %>%
    
    group_by(Tile_alt,Response,Predictor,Predictor1,Predictor2,Predictor3,Predictor4,Predictor5,Predictor6,Predictor7,Predictor8) %>%
    mutate(Predictor8_n=ifelse(!is.na(Predictor8),n(),NA)) %>%
    
    group_by(Tile_alt,Response,Predictor,Predictor1,Predictor2,Predictor3,Predictor4,Predictor5,Predictor6,Predictor7,Predictor8,Predictor9) %>%
    mutate(Predictor9_n=ifelse(!is.na(Predictor9),n(),NA)) %>%
    
    group_by(Tile_alt,Response,Predictor,Predictor1,Predictor2,Predictor3,Predictor4,Predictor5,Predictor6,Predictor7,Predictor8,Predictor9,Predictor10) %>%
    mutate(Predictor10_n=ifelse(!is.na(Predictor10),n(),NA)) %>%
    
    mutate(Unique=ifelse(Predictor_n==1 |Predictor1_n==1 | Predictor2_n==1 | Predictor3_n==1 | Predictor4_n==1 | Predictor5_n==1 | Predictor6_n==1
                         | Predictor7_n==1 | Predictor8_n==1 | Predictor9_n==1 | Predictor10_n==1
                         , "Y","N")) %>%
    filter(Unique=="Y") %>%
    
    #remove nested pathways 2
    
    mutate(
      Check_repeated1 = paste0(na.omit(c(Tile_alt,Response,Predictor)), collapse=","),
      Check_repeated2 = paste0(na.omit(c(Tile_alt,Predictor,Predictor1)), collapse=",")
    ) %>%
    #group_by(Tile_alt) %>%
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
    )
  
  return(SEM_out3)
  
}

SEM_alpha_effects<-SEM_coefs %>%
  mutate(Tile_alt=EmGroup) %>%
  Obtain_SEM_effects_unique %>%
  rowwise() %>%
  mutate(
    Second_Last_Predictor=ifelse(Predictor!=Last_Predictor,tail(na.omit(c(Predictor,Predictor1,Predictor2,Predictor3,Predictor4,Predictor5,Predictor6,Predictor7,Predictor8,Predictor9,Predictor10)), n=2),NA),
    Second_Predictor=ifelse(Predictor!=Last_Predictor & Predictor=="Total_cov_sd",Predictor1,Predictor)
  ) %>% 
  mutate(
    Interaction = ifelse(Last_Predictor!="Tile.Type2", "Covariates",
                         ifelse(is.na(Second_Last_Predictor),"1",
                                ifelse(Second_Last_Predictor=="Ac_invertebrates_mean", "2",
                                       ifelse(Second_Last_Predictor=="Barnacles1_mean", "3",
                                              ifelse(Second_Last_Predictor=="S_mean_log","4",
                                                     ifelse(Second_Last_Predictor=="Consumer_mean_log", "5",
                                                            "1")))))),
    Pathway=ifelse(Predictor=="Total_cov_mean", "Total_cov_mean",
                   ifelse(Predictor=="Total_cov_sd", "Total_cov_sd",NA)),
    AltGroup=ifelse(Tile_alt=="[0.693,0.764)", "Low",
                    ifelse(Tile_alt=="[0.764,0.819)", "Mid",
                           ifelse(Tile_alt=="[0.819,0.923]", "High",NA))),
    Tile_alt=ifelse(Tile_alt=="[0.0664,0.139)", "Low",
                    ifelse(Tile_alt=="[0.1394,0.193)", "Mid",
                           ifelse(Tile_alt=="[0.1934,0.296]", "High",NA))),
    Path2=ifelse(Path=="Direct", "Indirect",Path)
  ) %>%
  mutate(
    Path2=factor(Path2,c("Direct","Indirect","Covariates")),
    AltGroup=factor(AltGroup,c("Low","Mid","High")),
    Tile_alt=factor(Tile_alt,c("Low","Mid","High"))
  )

library(RColorBrewer)
brewer.pal(7, "YlGn")
brewer.pal(5, "BuPu")
brewer.pal(3, "Oranges")

Alpha_cols<-c(brewer.pal(7, "BuPu")[c(2:5)],
              brewer.pal(4, "YlOrRd")[c(1:3)],
              brewer.pal(5, "Greens")[c(2:5)])

ggplot(SEM_alpha_effects  %>%
         filter(Last_Predictor=="Tile.Type2") %>%
         filter((Response=="Compositional_similarity_z" | Response=="Cover_stability_alpha_res_log_z" | Response=="Pop_stab_log_z" | Response=="Asynchrony_corr_time_z")) %>% 
         mutate(Response=factor(Response, c("Pop_stab_log_z", "Asynchrony_corr_time_z", "Statistical_averaging_z", "Compositional_similarity_z", "Cover_stability_alpha_res_log_z"))) %>%
         group_by(Tile_alt, Predictor, Response) %>%
         summarise(Std.Estimate.Int=sum(Std.Estimate.Int)) %>%
         mutate(Predictor=factor(Predictor,c("Ac_invertebrates_mean_z","Barnacles1_mean_z", "Barnacles2_mean_sqrt_z","Eph_macroalgae_mean_sqrt_z","Consumer_mean_log_z","Pielou_mean_z","S_mean_log_z",  
                                             "Pop_stab_log_z","Asynchrony_corr_time_z","Statistical_averaging_z","Compositional_similarity_z","Tile.Type2")))
       , aes(x=Tile_alt, y=Std.Estimate.Int, fill=Predictor)) +
  geom_bar(position="stack", stat="identity")+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  facet_wrap(~Response,ncol=1)+
  theme_bw() +
  coord_flip() +
  scale_y_continuous(position = "right", lim=c(-0.45,0.6)) +
  scale_fill_manual(values=Alpha_cols) +
  theme(strip.background = element_blank()
        , strip.text.x = element_blank()
  )

ggplot(SEM_alpha_effects  %>%
         #filter(Last_Predictor=="Tile.Type2") %>%
         filter((Response=="Cover_stability_alpha_res_log_z")) %>% 
         group_by(Tile_alt, Response, Last_Predictor) %>%
         summarise(Std.Estimate.Int=sum(Std.Estimate.Int))
       , aes(x=Tile_alt, y=Std.Estimate.Int, fill=Last_Predictor)) +
  geom_bar(position="stack", stat="identity")+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  facet_wrap(~Response,ncol=1)+
  theme_bw() +
  coord_flip() +
  scale_y_continuous(position = "right", lim=c(-0.45,0.6)) +
  scale_fill_manual(values=Alpha_cols[c(10,9,8)]) +
  theme(strip.background = element_blank()
        , strip.text.x = element_blank()
  )


ggplot(SEM_alpha_effects  %>%
         filter((Response=="Compositional_similarity_z" |Response=="Cover_stability_alpha_res_log_z")) %>% 
         group_by(Tile_alt, Last_Predictor, Response) %>%
         summarise(Std.Estimate.Int=sum(Std.Estimate.Int)) #%>%
       # mutate(Last_Predictor=factor(Last_Predictor, c("Emersion.Rate", "Emersion.Rate2", "Tile.Type2")))
       , aes(x=Tile_alt, y=Std.Estimate.Int, fill=Last_Predictor)) +
  geom_bar(position="stack", stat="identity")+
  facet_wrap(~Response,ncol=1)+
  #ylim(-0.75, 2.75)+
  theme_bw() +
  coord_flip() +
  scale_y_continuous(position = "right") +
  scale_fill_manual(values=Alpha_cols[c(10,9,8)]) +
  theme(strip.background = element_blank()
        #, strip.text.x = element_blank()
  )

#