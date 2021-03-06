#### a hierarchical GAM analysis of estimated population trajectories
### trajectories from the Rosenberg et al. 2019 paper for all non-BBS species
#### and from the CWS GAMYE analyis of the 2019 BBS data for all BBS species
library(tidyverse)
library(ggrepel)
library(ggforce)
library(rjags)
library(jagsUI)
library(mgcv)
source("utility_functions.R")

ind = "index"
lci = "lci"
uci = "uci"
year = "year"
base.yr = 1970
popsource = "Pop.source"

set.seed(2019)


spslist = read.csv("data/Rosenberg et al species list.csv")

## merge with generation time file
# loading generation times from Bird et al 2020 supplemental ---------------------------

gens = read.csv("data/cobi13486-sup-0004-tables4.csv")

gens <- gens %>% select(Scientific_name,
                        GenLength)

spslist <- left_join(spslist,gens,by = c("sci_name" = "Scientific_name"))

# reconciling scientific names --------------------------------------------
sps_nomatch <- spslist[which(is.na(spslist$GenLength)),"species"]
fullgensnames = read.csv("data/cobi13486-sup-0001-tables1.csv")
sps_altmatch <- sps_nomatch[sps_nomatch %in% unique(fullgensnames$Common_name)]
for(s1 in sps_altmatch){
  scin1 = fullgensnames[which(fullgensnames$Common_name == s1),"Scientific_name"]
  spslist[which(spslist$species == s1),"GenLength"] <- gens[which(gens$Scientific_name == scin1),"GenLength"]
}


sp_still_miss <- spslist[which(is.na(spslist$GenLength)),"species"]


spslist[which(spslist$species == "American Three-toed Woodpecker"),"GenLength"] <- gens[which(gens$Scientific_name == "Picoides tridactylus"),"GenLength"]
spslist[which(spslist$species == "Black-necked Stilt"),"GenLength"] <- gens[which(gens$Scientific_name == "Himantopus himantopus"),"GenLength"]
spslist[which(spslist$species == "Black Oystercatcher"),"GenLength"] <- gens[which(gens$Scientific_name == "Haematopus ater"),"GenLength"]
spslist[which(spslist$species == "Green Heron"),"GenLength"] <- gens[which(gens$Scientific_name == "Butorides striata"),"GenLength"]
spslist[which(spslist$species == "Hoary Redpoll"),"GenLength"] <- gens[which(gens$Scientific_name == "Acanthis flammea"),"GenLength"]
spslist[which(spslist$species == "Woodhouse's Scrub-Jay"),"GenLength"] <- gens[which(gens$Scientific_name == "Aphelocoma coerulescens"),"GenLength"]

spslist[which(is.na(spslist$GenLength)),"species"]

write.csv(spslist[,c("species","sci_name","GenLength")],"Generation_lengths_Rosenberg_merge_Birdetal.csv",row.names = FALSE)

spsBBS = spslist[grepl(pattern = "BBS",spslist$Tr_source),"species"]
spsCBC = spslist[grepl(pattern = "CBC",spslist$Tr_source),"species"]
spsShorebird = spslist[grepl(pattern = "Mig7416",spslist$Tr_source),"species"]

indicesRosen = read.csv("data/Rosenberg et al annual indices of abundance.csv",
                      stringsAsFactors = F)
yr_span <- unique(indicesRosen[,c("species","firstyear","lastyear")])

# ## read in the full suite of annual indices from the GAMYE model - CWS 2019 BBS estimates
# ## select continental long-term indices and export to smaller csv file to include in Git repo
# ## does not need to be run because full csv file ~ 1GB in size
# indicesgamye = read.csv("data/arch/2019All BBS indices continent and national.csv",
#                       stringsAsFactors = F)
# indicesgamye[which(indicesgamye$species == "Dark-eyed Junco (all forms)"),"species"] <- "Dark-eyed Junco"
# indicesgamye[which(indicesgamye$species == "Northern Flicker (all forms)"),"species"] <- "Northern Flicker"
# indicesgamye[which(indicesgamye$species == "Red-tailed Hawk (all forms)"),"species"] <- "Red-tailed Hawk"
# indicesgamye[which(indicesgamye$species == "Yellow-rumped Warbler (all forms)"),"species"] <- "Yellow-rumped Warbler"
# ## select the continental estimates for the BBS species
# indicesraw <- indicesgamye %>%
#   filter(Region == "Continental",
#          Trend_Time == "Long-term",
#          species %in% spsBBS)
# 
# write.csv(indicesraw,"data/2019_BBS_continental_indices_gamye.csv",row.names = FALSE)
# rm(list = "indicesgamye")


indicesraw <- read.csv("data/2019_BBS_continental_indices_gamye.csv")

indicesraw <- left_join(indicesraw,yr_span,by = "species")

sp_wnew_bbs <- unique(indicesraw$species)
sp_missing_new_bbs <- spsBBS[-which(spsBBS %in% sp_wnew_bbs)]
sp_missing_new_bbs
#[1] "Common Ground-Dove"     "Eurasian Collared-Dove" "Pacific Golden-Plover" not included in model output - model failures


indicesraw <- indicesraw %>% 
  group_by(species) %>% 
  filter(Year >= firstyear,
         species %in% spsBBS)

indicesraw$lci.raw <- indicesraw$Index_q_0.025
indicesraw$uci.raw <- indicesraw$Index_q_0.975
indicesraw$index.raw <- indicesraw$Index
indicesraw$year <- indicesraw$Year

indicesraw = indicesraw[order(indicesraw$species,indicesraw$year),]

indicesBBS <- indicesraw %>% select(species,year,index.raw,lci.raw,uci.raw,firstyear,lastyear)

sp_wnew_bbs <- unique(indicesBBS$species)
# merge bbs indices with remaining Rosenberg indices ---------------------
spRosen <- unique(indicesRosen$species)
spRosen <- spRosen[-which(spRosen %in% sp_wnew_bbs)]
## drop the species with new BBS from the Rosenberg indices
indicesRosen <- filter(indicesRosen,
                       species %in% spRosen)



# Inserting updated CBC indices -------------------------------------------


indicesraw <- read.csv("data/cbc_indices_Continental_1966-2019.csv")

indicesraw[which(grepl(indicesraw$scientific_name_cbc,pattern = "Branta hutchinsii canadensis")),"english_name"] <- "Canada Goose"
tmp = indicesraw[which(grepl(indicesraw$scientific_name_cbc,pattern = "Aechmophorus")),]
indicesraw[which(grepl(indicesraw$scientific_name_cbc,pattern = "Aechmophorus")),"english_name"] <- "Clark's Grebe"
tmp$english_name <- "Western Grebe"

indicesraw <- bind_rows(indicesraw,tmp)


indicesraw <- left_join(indicesraw,yr_span,by = c("english_name" = "species"))

sp_wnew_cbc <- unique(indicesraw$english_name)
sp_missing_new_cbc <- spsCBC[-which(spsCBC %in% sp_wnew_cbc)]
sp_missing_new_cbc


indicesraw <- indicesraw %>% 
  filter(count_year >= firstyear)

indicesraw$lci.raw <- indicesraw$estimate_lcl 
indicesraw$uci.raw <- indicesraw$estimate_ucl 
indicesraw$index.raw <- indicesraw$estimate_median 
indicesraw$year <- indicesraw$count_year
indicesraw$species <- indicesraw$english_name

indicesraw = indicesraw[order(indicesraw$species,indicesraw$year),]

indicesCBC <- indicesraw %>% select(species,year,index.raw,lci.raw,uci.raw,firstyear,lastyear) %>% 
  filter(species %in% spsCBC)

sp_wnew_cbc <- unique(indicesCBC$species)


# merge CBC indices with remaining Rosenberg indices ---------------------
spRosen <- unique(indicesRosen$species)
spRosen <- spRosen[-which(spRosen %in% spsCBC)]

indicesRosen <- filter(indicesRosen,
                       species %in% spRosen)




# Inserting new shorebird migration survey trajectories -------------------


indicesraw <- read.csv("data/Shorebird_indices_Continental_gamma_t.csv")

#indicesraw <- read.csv("data/Shorebird_indices_Continental_1980-2019.csv")

yr_span[which(yr_span$species %in% spsShorebird),"firstyear"] <- 1980

indicesraw <- left_join(indicesraw,yr_span,by = c("species" = "species"))

sp_wnew_Shorebird <- unique(indicesraw$species)
sp_missing_new_Shorebird <- spsShorebird[-which(spsShorebird %in% sp_wnew_Shorebird)]
sp_missing_new_Shorebird

indicesraw <- indicesraw %>% 
  filter(year >= firstyear,
         parm == "N")

indicesraw$lci.raw <- indicesraw$lci 
indicesraw$uci.raw <- indicesraw$uci 
indicesraw$index.raw <- indicesraw$median 

indicesraw = indicesraw[order(indicesraw$species,indicesraw$year),]

indicesShorebird <- indicesraw %>% select(species,year,index.raw,lci.raw,uci.raw,firstyear,lastyear) %>% 
  filter(species %in% spsShorebird)

sp_wnew_Shorebird <- unique(indicesShorebird$species)


# merge bbs indices with remaining Rosenberg indices ---------------------
spRosen <- unique(indicesRosen$species)
spRosen <- spRosen[-which(spRosen %in% spsShorebird)]

indicesRosen <- filter(indicesRosen,
                       species %in% spRosen)









indicesAll <- bind_rows(indicesBBS,indicesCBC,indicesShorebird,indicesRosen)








sp_all <- unique(indicesAll$species)


# GAM smoothing of indices ------------------------------------------------

indicesAll <- indicesAll %>% 
  mutate(lind = log(index.raw),
         lsd = (log(uci.raw)-log(lci.raw))/(1.96*2))




#nyears_recent = 15

# Species loop ------------------------------------------------------------


jj = 0
for(ss in sp_all[1:length(sp_all)]){
  jj = jj+1
  
  nyears_recent <- min(max(ceiling(spslist[which(spslist$species == ss),"GenLength"]*3),10),25)

  
  
  wss = which(indicesAll$species == ss)
  tmp = indicesAll[wss,]
  
  fyr = max(c(unique(tmp$firstyear),min(tmp$year)))
  lyr = max(c(unique(tmp$lastyear),max(tmp$year)))
  
  torep = which(indicesAll$species == ss & (indicesAll$year >= fyr & indicesAll$year <= lyr))
  
  if(any(!is.na(tmp$lind))){
    wprec = T
  }else{
    wprec = F
  }
  
  tmpd = tmp[which(!is.na(tmp$lind)),]
  
 
  nknots = min(c(13,max(floor(nrow(tmpd)/3),3)))
  if(ss %in% c("Cackling Goose","Greater White-fronted Goose","Trumpeter Swan")){nknots = 4} 
  
  
  
  form = as.formula(paste("lind","~",
                          "s(year,k =",nknots,")"))
  
  
  
  
  ncounts = nrow(tmpd)
  ###### building the GAM basis function
  # gam basis functions created using function jagam from mgcv package          
  
  yminy = min(tmpd$year,na.rm = T)
  ymaxy = max(tmpd$year,na.rm = T)
  yearvec = tmpd$year-(yminy-1)
  yrs = seq(yminy,ymaxy,by = 1)
  nyears = length(yrs)
  ymin = min(yearvec)
  ymax = max(yearvec)
  lindex = tmpd$lind
  preci = 1/(tmpd$lsd^2)  
  
  if(any(is.na(preci))){
    w_sub <- which(is.na(preci))
    preci[w_sub] = 1/((lindex[w_sub]*0.15)^2) # if no uncertainty values - replace with 15% CV
  }
  
  preddat = data.frame(lind = 1,
                       year = yrs)
  
  
  nyears_recent <- min(floor(nyears*0.75),nyears_recent)

    
    gamprep = jagam(formula = form,
                    data = tmpd,
                    file = "tempgam.txt",
                    centred = T)
    
    gamprep.pred = jagam(formula = form,
                         data = preddat,
                         file = "tempgampred.txt",
                         centred = T)
    
    
    
    dat = list(X = gamprep$jags.data$X,
               S1 = gamprep$jags.data$S1,
               ncounts = nrow(tmpd),
               lindex = lindex,
               nknots = nknots,
               preci = preci,
               X.pred = gamprep.pred$jags.data$X,
               nyears = nyears,
               zero = gamprep$jags.data$zero,
               Ys = as.integer(ymin),
               Ye = as.integer(ymax),
               Yb = as.integer(ymax-nyears_recent))
    
    
    mgo = jagsUI(data = dat,
                 model.file = paste0("models/JAGS_model_GAM_Asymetric.R"),
                 n.chains = 3,
                parameters.to.save = c("ind.pred",
                         "C1",
                         "C2",
                         "C3",
                         "T1", # early short-term change
                         "T2", # recent short-term change
                         "T3", # long-term change
                         "Tdif",
                         "Tdif_neg"
                         #"rho"
                         #"mu",
                         #"b"
                       ),
                       n.iter = 20000,
                       n.burnin = 10000,
                       n.thin = 10,
                parallel = T)
    
    
    #mgosum = data.frame(mgo$summary)
    mgosum2 = summary(mgo$samples,quantiles = c(0.025,0.05,0.25,0.5,0.75,0.95,0.975))
    mgosum = as.data.frame(cbind(mgosum2$statistics,mgosum2$quantiles))
    #names(mgosum) <- gsub(names(mgosum),pattern = "%",replacement = ".",fixed = TRUE)
    names(mgosum) <- c("mean","sd","NaiveSE","T_S_SE",
                       "X2_5_","X5_","X25_","X50_","X75_","X95_","X97_5_")
    
    mgosum$parameter <- row.names(mgosum)
    
    mgosum$species = ss
    mgosum$nyears_recent <- nyears_recent
    mgosum$firstyear <- yminy
    mgosum$lastyear <- ymaxy
    mgosum$Three_gen_time_used <- nyears_recent
    mgosum$start_recent_trend <- ymaxy-nyears_recent
    
    if(jj == 1){
      out = mgosum
    }else{
      out <- bind_rows(out,mgosum)
    }
    
    print(round(jj/nrow(spslist),2))

}

save(list = c("out"),file = "output/temp_out_all_species.RData")





# END of the species loop -------------------------------------------------



### New figures to select declining species, 
# show trend-windows for the early and late. 
# Colour the trajectories based on the probability of accelerating decline.

load("output/temp_out_all_species.RData")
#names(out)[which(grepl(names(out),pattern = ".",fixed = T))] <- paste0(gsub(names(out)[which(grepl(names(out),pattern = ".",fixed = T))] ,pattern = ".",replacement = "_", fixed = TRUE))

mu_ <- filter(out,parameter %in% c(paste0("ind.pred[",1:53,"]")))
mu_ <- mutate(mu_,yr = jags_dim(dat = mu_, cl = "parameter",var = "ind.pred"),
              est = "mu")
# fyr_sp <- unique(indicesAll[,c("species","firstyear")])
# for(i in 1:nrow(fyr_sp)){
#   sp = as.character(fyr_sp[i,"species"])
#   fy = (fyr_sp[i,"firstyear"])
#   Tfy = min(indicesAll[which(indicesAll$species == sp & !is.na(indicesAll$index.raw)),"year"])
#   if(Tfy != fy){
#     fyr_sp[i,"firstyear"] <- Tfy
#   }
#   chp <- unique(mu_[which(mu_$species == sp),"nyears_recent"])
#   fyr_sp[i,"Three_gen_time_used"] <- chp
#   fyr_sp[i,"lastyear"] <- max(indicesAll[which(indicesAll$species == sp & !is.na(indicesAll$index.raw)),"year"])
#   fyr_sp[i,"start_recent_trend"] <- fyr_sp[i,"lastyear"]-chp
#     
# }
# mu_ <- left_join(mu_,fyr_sp,by = "species")

fyr_sp <- unique(mu_[,c("species","firstyear","lastyear","start_recent_trend","Three_gen_time_used")])
mu_$year <- mu_$yr+(mu_$firstyear-1)






Tdif <- filter(out,parameter %in% c(paste0("Tdif")))
Tdif_neg <- filter(out,parameter %in% c(paste0("Tdif_neg")))

T1 <- filter(out,parameter %in% c(paste0("T1")))
T2 <- filter(out,parameter %in% c(paste0("T2")))
T3 <- filter(out,parameter %in% c(paste0("T3")))


prob_annot <- select(Tdif_neg,mean,species)
names(prob_annot)[1] <- "mean_p"


# Using 90% CIs -----------------------------------------------------------

bdif <- select(Tdif,X5_,mean,X95_,species)
bdif <- left_join(bdif,prob_annot,by = "species")
bdif$Difference_late_minus_early_trend <- round(bdif$mean,1)
bdif$UCI_90_Difference_late_minus_early_trend <- round(bdif$X95_,1)
bdif$LCI_90_Difference_late_minus_early_trend <- round(bdif$X5_,1)
bdif$prob_decreasing_trend <- round(bdif$mean_p,2)

T1$early_trend <- round(T1$mean,1)
T1$LCI_early_trend <- round(T1$X5_,1)
T1$UCI_early_trend <- round(T1$X95_,1)

T2$late_trend <- round(T2$mean,1)
T2$LCI_late_trend <- round(T2$X5_,1)
T2$UCI_late_trend <- round(T2$X95_,1)

T3$long_term_trend <- round(T3$mean,1)
T3$LCI_long_term_trend <- round(T3$X5_,1)
T3$UCI_long_term_trend <- round(T3$X95_,1)



bdif <- left_join(bdif,T3[,c("species","long_term_trend","LCI_long_term_trend","UCI_long_term_trend")])
bdif <- left_join(bdif,T1[,c("species","early_trend","LCI_early_trend","UCI_early_trend")])
bdif <- left_join(bdif,T2[,c("species","late_trend","LCI_late_trend","UCI_late_trend")])



bdif <- left_join(bdif,fyr_sp)

sp_by_survey <- unique(spslist[,c("species","Tr_source")])

sp_by_survey[which(sp_by_survey$species %in% sp_wnew_bbs),"Tr_source"] <- "BBS_GAM"
sp_by_survey[which(sp_by_survey$species %in% sp_wnew_cbc),"Tr_source"] <- "CBC_2019"
sp_by_survey[which(sp_by_survey$species %in% sp_wnew_Shorebird),"Tr_source"] <- "Shorebird_GAM"

# spsBBS = spslist[grepl(pattern = "BBS",spslist$Tr_source),"species"]
# spsCBC = spslist[grepl(pattern = "CBC",spslist$Tr_source),"species"]
# spsShorebird = spslist[grepl(pattern = "Mig7416",spslist$Tr_source),"species"]

bdif <- left_join(bdif,sp_by_survey,by = "species")

write.csv(bdif[,c("species","firstyear","Three_gen_time_used","lastyear","start_recent_trend",
                  "Difference_late_minus_early_trend","LCI_90_Difference_late_minus_early_trend","UCI_90_Difference_late_minus_early_trend",
                  "prob_decreasing_trend",
                  "long_term_trend","LCI_long_term_trend","UCI_long_term_trend",
                  "early_trend","LCI_early_trend","UCI_early_trend",
                  "late_trend","LCI_late_trend","UCI_late_trend",
                  "Tr_source")],"output/Differences_in_Trends_GAM_all_species.csv",row.names = F)


bdif_plot <- bdif %>% filter(Tr_source %in% c("BBS_GAM","CBC_2019","Shorebird_GAM"))
distr_change <- ggplot(data = bdif_plot,aes(x = prob_decreasing_trend))+
  geom_histogram(bins = 10)+
  facet_wrap(~Tr_source,scales = "free_y")

png_mag = 5

png(paste0("output/images/probability_change_among_surveys_demo.png"),
    res = 100*png_mag,
    width = 480*png_mag,
    height = 480*png_mag)
print(distr_change)
dev.off()



sp_annot <- fyr_sp %>% mutate(year1 = firstyear,#max(firstyear,floor(start_recent_trend-(Three_gen_time_used))),
                              year2 = lastyear,#floor(start_recent_trend+(Three_gen_time_used)),
                              year3 = floor(firstyear+((start_recent_trend-firstyear)/2)),
                              year4 = start_recent_trend)
sp_annot <- as.data.frame(sp_annot)
for(sp in sp_annot$species){
  i = which(sp_annot$species == sp)
  tmp = indicesAll[which(indicesAll$species == sp),]
  tmpmu = mu_[which(mu_$species == sp),]
  y1 = sp_annot[i,"year1"]
  y2 = sp_annot[i,"year2"]
  y3 = sp_annot[i,"year3"]
  y4 = sp_annot[i,"year4"]
  
  yup = max(c(quantile(tmp$uci.raw,0.95,na.rm = T),max(tmp$index.raw,na.rm = T)*1.1),na.rm = T)
  yr = min(c(tmp$lci.raw,tmp$index.raw,yup*0.1),na.rm = T)
  
  sp_annot[i,"index1"] <- tmpmu[which(tmpmu$year == y1),"mean"]
  sp_annot[i,"index2"] <- tmpmu[which(tmpmu$year == y2),"mean"]
  sp_annot[i,"index4"] <- tmpmu[which(tmpmu$year == y4),"mean"]
  sp_annot[i,"index3"] <- yr
  try(sp_annot[i,"index5"] <- tmpmu[which(tmpmu$year == y4),"mean"])
  
  if(is.na(sp_annot[i,"index1"])){sp_annot[i,"index1"]<- yup*0.7} 
  if(is.na(sp_annot[i,"index2"])){sp_annot[i,"index2"]<- sp_annot[i,"index1"]} 
  if(is.na(sp_annot[i,"index3"])){sp_annot[i,"index3"]<- sp_annot[i,"index1"]} 
  
}

sp_annot <- sp_annot %>% select(species,
                                year1,
                                year2,
                                year3,
                                year4,
                                index1,
                                index2,
                                index3,
                                index4,
                                index5) %>% 
  left_join(.,bdif,by = "species") %>% 
  mutate(lab1 = paste0(early_trend,"% ",year1,"-",start_recent_trend),
         lab2 = paste0(late_trend,"% ",start_recent_trend,"-",lastyear),
         lab4 = paste0("p_down = ",round(mean_p,2)),
         lab3 = paste0(long_term_trend,"%",firstyear,"-",lastyear))

plot_ts1 <- NULL
plot_ts2 <- NULL
for(sp in sp_annot$species){
  if(sp == "Brant"){next}
  i = which(sp_annot$species == sp)
  ys <- c((sp_annot[i,"firstyear"]),
          sp_annot[i,"start_recent_trend"])
  tmp <- mu_[which(mu_$species == sp & mu_$year %in% ys),]
  tmp$period = "early"
  plot_ts1 <- bind_rows(plot_ts1,tmp)
  ys <- c(sp_annot[i,"start_recent_trend"],
          sp_annot[i,"lastyear"])
  tmp <- mu_[which(mu_$species == sp & mu_$year %in% ys),]
  tmp$period = "late"
  plot_ts2 <- bind_rows(plot_ts2,tmp)
  
  
}

levels = c("early","late","long_term","p_downturn","long_term_decline","p_downturn_high")
labs_rep = NULL
for(j in 1:4){
  nm <- c("species","year","index","lab","long_term_trend","mean_p")
  tmp <- sp_annot[,c(nm[1],paste0(nm[2:4],j),nm[5],nm[6])]
  names(tmp) <- nm
  tmp$group <- factor(levels[j],levels = levels,ordered = TRUE)
  labs_rep <- bind_rows(labs_rep,tmp)
}
labs_rep$group[which(labs_rep$long_term_trend < -0.5 &
                       labs_rep$group == "long_term")] <- "long_term_decline"
labs_rep$group[which(labs_rep$long_term_trend < -0.5 &
                       labs_rep$mean_p > 0.75 &
                       labs_rep$group == "p_downturn")] <- "p_downturn_high"

cls = scales::viridis_pal(alpha = 0.8,begin = 0,end = 0.8,option = "magma",direction = -1)(5)

plot(1:5,col = cls,pch = 20,cex = 4)

cls <- c(cls[c(1,3)],grey(0.75),grey(0.75),cls[c(2,5)])
names(cls) <- levels


# i90 <- filter(indicesAll,year == 2000)
# bdiflab <- left_join(i90,bdif,by = "species")
# #bdiflab$lab = paste(bdiflab$Difference_late_minus_early_trend,":",bdiflab$LCI_90_Difference_late_minus_early_trend,"-",bdiflab$UCI_90_Difference_late_minus_early_trend,"p =",bdiflab$prob_decreasing_trend)
# bdiflab$lab = paste0(bdiflab$long_term_trend,"%","pAccelDecl = ",bdiflab$prob_decreasing_trend," dif=",bdiflab$Difference_late_minus_early_trend)
# bdiflab$index.raw <- bdiflab$index.raw*1.1
txt_s <- 3#font size in plot below
npag <- ceiling(nrow(bdif)/9)
clm = scales::viridis_pal(alpha = 1,direction = 1,end = 0.8,begin = 0.2)(2)
clfil <- rgb(255, 255, 255, max = 255, alpha = floor(255*0.35)) #mostly transparent white

pdf("output/changepoint_graphs_gam_all_species.pdf",
    height = 8.5,width = 11)
for(ij in 1:npag){
  trajs <- ggplot(data = mu_,aes(x = year,y = mean))+
    geom_pointrange(data = indicesAll,inherit.aes = FALSE,
                    aes(x = year,y = index.raw,ymax = uci.raw,ymin = lci.raw),
                    alpha = 0.1, size = 0.3)+
    geom_ribbon(aes(ymin = X5_,ymax = X95_),alpha = 0.2,fill = clm[2])+
    geom_line(colour = clm[2])+
    geom_point(data = sp_annot,inherit.aes = FALSE,
               aes(x = year4,y = index5),colour = "red")+
    geom_line(data = plot_ts1,inherit.aes = FALSE,
               aes(x = year,y = mean),colour = cls[1])+
    geom_line(data = plot_ts2,inherit.aes = FALSE,
              aes(x = year,y = mean),colour = cls[2])+
    geom_label_repel(data = labs_rep,inherit.aes = FALSE,
                     aes(x = year, y = index, colour = group,label = lab),
                     min.segment.length = 0.1, #vjust = 1,
                     size = txt_s,fill = clfil,
                     )+
    scale_colour_manual(values = cls)+
    scale_y_continuous(trans = "log",labels = scales::comma)+
    theme_classic()+
    theme(legend.position = "none")+
    facet_wrap_paginate(facets = ~species,scales = "free_y",nrow = 3,ncol = 3,page = ij)
  print(trajs)
}
dev.off()



# plotting graphs for selected species ------------------------------------

sp_sel <- c("Barn Swallow",
            "Bank Swallow",
            "Chimney Swift",
            "Common Grackle",
            "Greater Yellowlegs",
            "Lesser Yellowlegs",
            "Red Knot",
            "Sanderling",
            "Canada Warbler",
            "Loggerhead Shrike",
            "Clark's Grebe",
            "Snowy Owl")
mu_sel <- mu_ %>% filter(species %in% sp_sel)
indicesAll_sel <- indicesAll %>% filter(species %in% sp_sel)
sp_annot_sel <- sp_annot %>% filter(species %in% sp_sel)
labs_rep_sel <- labs_rep %>% filter(species %in% sp_sel)
plot_ts1_sel <- plot_ts1 %>% filter(species %in% sp_sel)
plot_ts2_sel <- plot_ts2 %>% filter(species %in% sp_sel)


pdf("output/changepoint_graphs_gam_selected_species.pdf",
    height = 8.5,width = 11)
for(ij in 1:3){
  trajs <- ggplot(data = mu_sel,aes(x = year,y = mean))+
    geom_pointrange(data = indicesAll_sel,inherit.aes = FALSE,
                    aes(x = year,y = index.raw,ymax = uci.raw,ymin = lci.raw),
                    alpha = 0.1, size = 0.3)+
    geom_ribbon(aes(ymin = X5_,ymax = X95_),alpha = 0.2,fill = clm[2])+
    geom_line(colour = clm[2])+
    geom_point(data = sp_annot_sel,inherit.aes = FALSE,
               aes(x = year4,y = index5),colour = "red")+
    geom_line(data = plot_ts1_sel,inherit.aes = FALSE,
              aes(x = year,y = mean),colour = cls[1])+
    geom_line(data = plot_ts2_sel,inherit.aes = FALSE,
              aes(x = year,y = mean),colour = cls[2])+
    geom_label_repel(data = labs_rep_sel,inherit.aes = FALSE,
                     aes(x = year, y = index, colour = group,label = lab),
                     min.segment.length = 0.1, #vjust = 1,
                     size = txt_s,fill = clfil,
    )+
    scale_colour_manual(values = cls)+
    scale_y_continuous(trans = "log",labels = scales::comma)+
    #scale_y_continuous(labels = scales::comma)+
    theme_classic()+
    theme(legend.position = "none")+
    facet_wrap_paginate(facets = ~species,scales = "free_y",nrow = 2,ncol = 2,page = ij)
  print(trajs)
}
dev.off()



# single species graphs ---------------------------------------------------
sp_sel1 <- "Barn Swallow"

for(sp_sel1 in sp_sel){

mu_sel <- mu_ %>% filter(species %in% sp_sel1)
indicesAll_sel <- indicesAll %>% filter(species %in% sp_sel1)
sp_annot_sel <- sp_annot %>% filter(species %in% sp_sel1)
labs_rep_sel <- labs_rep %>% filter(species %in% sp_sel1)
plot_ts1_sel <- plot_ts1 %>% filter(species %in% sp_sel1)
plot_ts2_sel <- plot_ts2 %>% filter(species %in% sp_sel1)

yup = max(c(max(indicesAll_sel$index.raw,na.rm = TRUE)*1.1,max(indicesAll_sel$uci.raw,na.rm = TRUE)))
ydn = min(c(indicesAll_sel$lci.raw,yup*0.1),na.rm = T)


png_mag = 5

  trajs <- ggplot(data = mu_sel,aes(x = year,y = mean))+
    geom_pointrange(data = indicesAll_sel,inherit.aes = FALSE,
                    aes(x = year,y = index.raw,ymax = uci.raw,ymin = lci.raw),
                    alpha = 0.4, size = 0.3)+
    ylab("Estimated population trajectory")+
    xlab("")+
    labs(title = sp_sel1)+
    scale_y_continuous(trans = "log",labels = scales::comma)+
    scale_x_continuous(limits = c(1970,2020),breaks = seq(1970,2020,by = 10))+
    #scale_y_continuous(labels = scales::comma)+
    theme_classic()+
    coord_cartesian(ylim = c(ydn,yup))+
    theme(legend.position = "none")

  stp = 1
  png(paste0("output/images/",sp_sel1,"_",stp,"_example.png"),
      res = 100*png_mag,
      width = 480*png_mag,
      height = 480*png_mag)
  print(trajs)
  dev.off()
  
  trajs <- ggplot(data = mu_sel,aes(x = year,y = mean))+
    geom_pointrange(data = indicesAll_sel,inherit.aes = FALSE,
                    aes(x = year,y = index.raw,ymax = uci.raw,ymin = lci.raw),
                    alpha = 0.2, size = 0.3)+
    ylab("Estimated population trajectory")+
    xlab("")+
    labs(title = sp_sel1)+
    scale_y_continuous(trans = "log",labels = scales::comma)+
    scale_x_continuous(limits = c(1970,2020),breaks = seq(1970,2020,by = 10))+
    #scale_y_continuous(labels = scales::comma)+
    theme_classic()+
    coord_cartesian(ylim = c(ydn,yup))+
    theme(legend.position = "none")
  
  
  trajs <- trajs +
    geom_ribbon(aes(ymin = X5_,ymax = X95_),alpha = 0.2,fill = clm[2])+
    geom_line(colour = clm[2])
  
  stp = stp+1
  png(paste0("output/images/",sp_sel1,"_",stp,"_example.png"),
      res = 100*png_mag,
      width = 480*png_mag,
      height = 480*png_mag)
  print(trajs)
  dev.off()
  
  trajs <- trajs+
    geom_point(data = sp_annot_sel,inherit.aes = FALSE,
               aes(x = year4,y = index5),colour = "red") +
    geom_line(data = plot_ts2_sel,inherit.aes = FALSE,
              aes(x = year,y = mean),colour = cls[2])
  
  stp = stp+1
  png(paste0("output/images/",sp_sel1,"_",stp,"_example.png"),
      res = 100*png_mag,
      width = 480*png_mag,
      height = 480*png_mag)
  print(trajs)
  dev.off()
  
  
  trajs <- trajs +
  geom_line(data = plot_ts1_sel,inherit.aes = FALSE,
            aes(x = year,y = mean),colour = cls[1])
  
  stp = stp+1
  png(paste0("output/images/",sp_sel1,"_",stp,"_example.png"),
      res = 100*png_mag,
      width = 480*png_mag,
      height = 480*png_mag)
  print(trajs)
  dev.off()
  

  trajs <- trajs +
    geom_label_repel(data = labs_rep_sel,inherit.aes = FALSE,
                     aes(x = year, y = index, colour = group,label = lab),
                     min.segment.length = 0.1, #vjust = 1,
                     size = txt_s,fill = clfil)+
    scale_colour_manual(values = cls)
  
  
  stp = stp+1
  png(paste0("output/images/",sp_sel1,"_",stp,"_example.png"),
      res = 100*png_mag,
      width = 480*png_mag,
      height = 480*png_mag)
  print(trajs)
  dev.off()
  
  

}
