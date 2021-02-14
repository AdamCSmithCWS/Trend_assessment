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



spsBBS = spslist[grepl(pattern = "BBS",spslist$Tr_source),"species"]

indicesRosen = read.csv("data/Rosenberg et al annual indices of abundance.csv",
                      stringsAsFactors = F)
yr_span <- unique(indicesRosen[,c("species","firstyear","lastyear")])

# ## read in the full suite of annual indices from the GAMYE model - CWS 2019 BBS estimates
# ## select continental long-term indices and export to smaller csv file to include in Git repo
# ## does not need to be run because full csv file ~ 1GB in size
# indicesgamye = read.csv("data/arch/All 2019 BBS indices.csv",
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
  filter(Year >= firstyear & Year <= lastyear)

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

indicesRosen <- filter(indicesRosen,
                       species %in% spRosen)


indicesAll <- bind_rows(indicesBBS,indicesRosen)

sp_all <- unique(indicesAll$species)


# GAM smoothing of indices ------------------------------------------------

indicesAll <- indicesAll %>% 
  mutate(lind = log(index.raw),
         lsd = (log(uci.raw)-log(lci.raw))/(1.96*2))




#nyears_recent = 15

jj = 0
for(ss in sp_all[1:length(sp_all)]){
  jj = jj+1
  
  nyears_recent <- min(max(ceiling(spslist[which(spslist$species == ss),"GenLength"]*3),10),25)

  
  
  wss = which(indicesAll$species == ss)
  tmp = indicesAll[wss,]
  
  fyr = unique(tmp$firstyear)
  lyr = unique(tmp$lastyear)
  
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
                 model.file = paste0("models/JAGS_model_GAM_symetric.R"),
                 n.chains = 3,
                parameters.to.save = c("ind.pred",
                         "C1",
                         "C2",
                         "C3",
                         "T1",
                         "T2",
                         "T3",
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
    
    
    mgosum = data.frame(mgo$summary)
    
    
    mgosum$parameter <- row.names(mgosum)
    
    mgosum$species = ss
    mgosum$nyears_recent <- nyears_recent
    
    if(jj == 1){
      out = mgosum
    }else{
      out <- bind_rows(out,mgosum)
    }
    
    print(round(jj/nrow(spslist),2))

}

save(list = c("out"),file = "output/temp_out_all_species.RData")


### New figures to select declining species, 
# show trend-windows for the early and late. 
# Colour the trajectories based on the probability of accelerating decline.


names(out)[which(grepl(names(out),pattern = ".",fixed = T))] <- paste0(gsub(names(out)[which(grepl(names(out),pattern = ".",fixed = T))] ,pattern = ".",replacement = "_", fixed = TRUE))

mu_ <- filter(out,parameter %in% c(paste0("ind.pred[",1:53,"]")))
mu_ <- mutate(mu_,yr = jags_dim(dat = mu_, cl = "parameter",var = "ind.pred"),
              est = "mu")
fyr_sp <- unique(indicesAll[,c("species","firstyear")])
mu_ <- left_join(mu_,fyr_sp,by = "species")
mu_$year <- mu_$yr+(mu_$firstyear-1)






Tdif <- filter(out,parameter %in% c(paste0("Tdif")))
Tdif_neg <- filter(out,parameter %in% c(paste0("Tdif_neg")))

T1 <- filter(out,parameter %in% c(paste0("T1")))
T2 <- filter(out,parameter %in% c(paste0("T2")))
T3 <- filter(out,parameter %in% c(paste0("T3")))


prob_annot <- select(Tdif_neg,mean,species)
names(prob_annot)[1] <- "mean_p"
bdif <- select(Tdif,X2_5_,mean,X97_5_,species)
bdif <- left_join(bdif,prob_annot,by = "species")
bdif$Difference_late_minus_early_trend <- round(bdif$mean,1)
bdif$UCI_90_Difference_late_minus_early_trend <- round(bdif$X97_5_,1)
bdif$LCI_90_Difference_late_minus_early_trend <- round(bdif$X2_5_,1)
bdif$prob_decreasing_trend <- round(bdif$mean_p,2)

T1$early_trend <- round(T1$mean,1)
T1$LCI_early_trend <- round(T1$X2_5_,1)
T1$UCI_early_trend <- round(T1$X97_5_,1)

T2$late_trend <- round(T2$mean,1)
T2$LCI_late_trend <- round(T2$X2_5_,1)
T2$UCI_late_trend <- round(T2$X97_5_,1)

T3$long_term_trend <- round(T3$mean,1)
T3$LCI_long_term_trend <- round(T3$X2_5_,1)
T3$UCI_long_term_trend <- round(T3$X97_5_,1)



bdif <- left_join(bdif,T3[,c("species","long_term_trend","LCI_long_term_trend","UCI_long_term_trend")])
bdif <- left_join(bdif,T1[,c("species","early_trend","LCI_early_trend","UCI_early_trend")])
bdif <- left_join(bdif,T2[,c("species","late_trend","LCI_late_trend","UCI_late_trend")])


write.csv(bdif[,c("species","Difference_late_minus_early_trend","LCI_90_Difference_late_minus_early_trend","UCI_90_Difference_late_minus_early_trend",
                  "prob_decreasing_trend",
                  "long_term_trend","LCI_long_term_trend","UCI_long_term_trend",
                  "early_trend","LCI_early_trend","UCI_early_trend",
                  "late_trend","LCI_late_trend","UCI_late_trend")],"output/Differences_in_Trends_GAM_altBBS.csv",row.names = F)

i90 <- filter(indicesAll,year == 2000)
bdiflab <- left_join(i90,bdif,by = "species")
#bdiflab$lab = paste(bdiflab$Difference_late_minus_early_trend,":",bdiflab$LCI_90_Difference_late_minus_early_trend,"-",bdiflab$UCI_90_Difference_late_minus_early_trend,"p =",bdiflab$prob_decreasing_trend)
bdiflab$lab = paste0(bdiflab$long_term_trend,"%","pAccelDecl = ",bdiflab$prob_decreasing_trend," dif=",bdiflab$Difference_late_minus_early_trend)
bdiflab$index.raw <- bdiflab$index.raw*1.1

npag <- ceiling(nrow(bdif)/9)
pdf("output/changepoint_graphs_gam_altBBS.pdf",
    height = 8.5,width = 11)
for(ij in 1:npag){
  trajs <- ggplot(data = mu_,aes(x = year,y = mean))+
    geom_pointrange(data = indicesAll,inherit.aes = FALSE,aes(x = year,y = index.raw,ymax = (uci.raw),ymin = (lci.raw)),alpha = 0.2, size = 0.5)+
    geom_ribbon(aes(ymin = X2_5_,ymax = X97_5_,fill = est),alpha = 0.2)+
    geom_line(aes(colour = est))+
    scale_color_viridis_d(aesthetics = c("colour","fill"),direction = 1,end = 0.8,begin = 0.2)+
    geom_text(data = bdiflab,inherit.aes = FALSE,aes(x = year,y = index.raw,label = lab))+
    scale_y_continuous(trans = "log",labels = scales::comma)+
    facet_wrap_paginate(facets = ~species,scales = "free_y",nrow = 3,ncol = 3,page = ij)
  print(trajs)
}
dev.off()

