library(dlnm)
library(mgcv)
library(dplyr)
library(data.table)
library(openxlsx)
library(readxl)
library(forecast)
library(splines)

#### 0.read data ####
data = fread("../Data_clean/Hospitalization/data0325.csv")
cc = table(data$county) %>% as.data.frame()
cc = cc$Var1[which(cc$Freq>700)] %>% as.character() %>% as.numeric()
data = data %>% subset(county %in% cc)

dis = names(data)[3:11]

## ..single event
data$drought = ifelse(data$event == 1 & data$DF!=1, 1, 0)
data$flood = ifelse(data$event == 2 & data$DF!=1, 1, 0)



#### 1.COV ####
edu = fread("../Data_clean/Covs/Education.csv")
fpl = fread("../Data_clean/Covs/FPL_2020.csv")
gdp = fread("../Data_clean/Covs/GDP_2005_2024_prediction.csv") %>% subset(Year %in% 2016:2024)
bed = fread("../Data_clean/Covs/Hospital_Bed.csv")
ndvi = fread("../Data_clean/Covs/NDVI_2015_2024.csv")
aging = fread("../Data_clean/Covs/Population_AgingRate.csv")
dens = fread("../Data_clean/Covs/Population_Density.csv")
sexr = fread("../Data_clean/Covs/Population_SexRatio.csv")
bed = bed %>% left_join(dens %>% dplyr::select(Year,County,Population), by = c("County","Year"))
bed$Bed.Num = bed$Bed.Num / bed$Population
gdp = gdp %>% left_join(dens %>% dplyr::select(Year,County,Population), by = c("County","Year"))
gdp$gdp_interp = gdp$gdp_interp / gdp$Population
fpl$FPL = ifelse(fpl$FPL %in% c("","0","No data","10–20"),"1",fpl$FPL)
fpl$FPL = ifelse(fpl$FPL %in% c("20–30"),"2",fpl$FPL)
fpl$FPL = ifelse(fpl$FPL %in% c("30–50"),"3",fpl$FPL)
fpl$FPL = ifelse(fpl$FPL %in% c("50–100"),"4",fpl$FPL)
fpl$FPL = ifelse(fpl$FPL %in% c("100–200"),"5",fpl$FPL)
fpl$FPL = ifelse(fpl$FPL %in% c("≥200"),"6",fpl$FPL)
fpl$FPL = as.numeric(fpl$FPL)

#### 2. DLNM ####
county = unique(data$county)
data.slice = lapply(county,function(x){
  data[data$county==x,]})
names(data.slice) = county

par(mfrow=c(3,3))
rr = list()
for (d in 1:length(dis)) {
  
  # coefall = matrix(NA,length(data.slice),1,
  #                 dimnames = list(county,paste("b",1,sep = "")))
  vcovall = vector("list",length(data.slice))
  names(vcovall) = county
  
  coefall = matrix(NA,length(data.slice),3,
                   dimnames = list(county,paste("b",1:3,sep = "")))
  
  
  for (i in seq(length(data.slice))) {
    print(paste(dis[d],i,county[i],sep="-"))
    sub = data.slice[[i]]
    sub = as.data.frame(sub)
    sub = sub[order(sub$date),]
    sub$time = c(1:nrow(sub))
    
    sub = sub %>% subset(select =  c("date","county",dis[d],"windows","drought","flood","tmlag021","rh","holiday","time"))
    names(sub)[3] = "case"
    zero = sum(sub$case==0)
    zerop = 100*zero / nrow(sub)
    
    
    if(zerop<20){
      imput = ifelse(sub$case==0, NA, sub$case)
      imput = na.interp(imput)
      sub$case = round(as.numeric(imput),0)
      sub$dow = lubridate::wday(sub$date)
      # dfc = crossbasis(sub$DF,
      #                 lag=30,
      #                 argvar=list(fun = "integer"),
      #                 arglag=list(knots = logknots(30,3))
      #                 )
      dfc = crossbasis(sub$windows,
                       lag=60,
                       argvar=list(fun = "integer"),
                       arglag=list(fun = "ns", df=3)
      )
      
      cb.drought = crossbasis(sub$drought,
                              lag=60,
                              argvar=list(fun = "integer"),
                              arglag=list(fun = "ns", df=3))
      cb.flood   = crossbasis(sub$flood,
                              lag=60,
                              argvar=list(fun = "integer"),
                              arglag=list(fun = "ns", df=3))
      dflong = round(nrow(sub) * (7/365),0)
      
      model =  gam(case~ dfc + cb.drought+cb.flood+
                     ns(tmlag021,3) + ns(rh,3) +
                     factor(holiday)+  factor(dow) + ns(time,df=dflong),
                   family=quasipoisson,
                   sub,
                   na.action="na.exclude")
      
      model_all =  crossreduce(dfc,model, cen=0)
      ### ..coef and vcov
      # coefall[i,] <- coef(model_all)
      # vcovall[[i]] <- vcov(model_all)
      coefall[i,] <- coef(model_all)
      vcovall[[i]] <- vcov(model_all)
    }
    
    if(zerop>=20){
      coefall[i,] <- NA
      vcovall[[i]] <- NA
    }
    
  }
  
  nouse = c()
  for (i in 1:length(vcovall)) {
    print(i)
    wd = dim(vcovall[[i]])[[1]]
    if(is.null(wd)) wd=0
    if(wd!=3) {nouse = c(i,nouse)}
  }
  
  coefall[nouse,] = NA
  for (i in nouse) {
    vcovall[[i]] = NA
  }
  
  ## ..meta
  coefall.2 = coefall[c(which(!is.na(coefall[,1]))),]
  rn = row.names(coefall.2)
  vcovall2 = list()
  kk = 1
  for (kk in 1:length(rn)) {
    vcovall2[[kk]] = vcovall[[rn[kk]]]
  }
  
  
  ## ..cov
  rn = as.numeric(rn)
  edu.1 = edu[which(edu$区划码 %in% rn)]; names(edu.1) = c("code","edu")
  fpl.1 = fpl[which(fpl$county %in% rn)]; names(fpl.1) = c("code","fpl")
  gdp.1 = gdp[which(gdp$County %in% rn & gdp$Year %in% 2016:2024)]
  gdp.1 = gdp.1 %>% group_by(County) %>% summarise(gdp = mean(gdp_interp, na.rm = T))
  names(gdp.1)[1] = "code"
  bed.1 = bed[which(bed$County %in% rn & bed$Year %in% 2016:2024)]
  bed.1 = bed.1 %>% group_by(County) %>% summarise(bed = mean(Bed.Num, na.rm = T))
  ndvi.1 = ndvi[which(ndvi$county.code %in% rn & ndvi$year %in% 2016:2024)]
  ndvi.1 = ndvi.1 %>% group_by(county.code) %>% summarise(NDVI = mean(NDVI, na.rm = T))
  aging.1 = aging[which(aging$County %in% rn & aging$Year %in% 2016:2024)]
  aging.1 = aging.1 %>% group_by(County) %>% summarise(Agingrate = mean(AgingRate, na.rm = T))
  dens.1 = dens[which(dens$County %in% rn & dens$Year %in% 2016:2024)]
  dens.1 = dens.1 %>% group_by(County) %>% summarise(Pop.Density = mean(Pop_Density, na.rm = T))
  sexr.1 = sexr[which(sexr$County %in% rn & sexr$Year %in% 2016:2024)]
  sexr.1 = sexr.1 %>% group_by(County) %>% summarise(SexRatio = mean(SexRatio, na.rm = T))
  
  names(edu.1)[1] <- names(fpl.1)[1] <- names(gdp.1)[1] <- names(ndvi.1)[1] <- "County"
  cov.data = data.frame(County = rn, order = 1:length(rn))
  cov.data = cov.data %>%
    merge(edu.1,by="County",all.x = T) %>% 
    merge(fpl.1,by="County",all.x = T) %>% 
    merge(gdp.1,by="County",all.x = T) %>% 
    merge(bed.1,by="County",all.x = T) %>% 
    merge(ndvi.1,by="County",all.x = T) %>% 
    merge(aging.1,by="County",all.x = T) %>% 
    merge(dens.1,by="County",all.x = T) %>% 
    merge(sexr.1,by="County",all.x = T)
  cov.data = cov.data %>% arrange(order)
  cov.data$edu[which(is.na(cov.data$edu))] = median(cov.data$edu,na.rm = T)
  cov.data = cov.data %>% mutate(
    edu.2 = scale(edu),
    fpl.2 = scale(fpl),
    gdp.2 = scale(gdp),
    bed.2 = scale(bed),
    NDVI.2 = scale(NDVI),
    Agingrate.2 = scale(Agingrate),
    Pop.Density.2 = scale(Pop.Density),
    SexRatio.2 = scale(SexRatio),
    
  
  )

  metaResult = mixmeta::mixmeta(coefall.2 ~ gdp.2+NDVI.2+Agingrate.2+fpl.2+edu.2,
                                vcovall2,data = cov.data,random = ~1|County,method = "reml")
  ## .. meta end
  
  tperc = rep(c(0,1,2,3),500)
  cb =  crossbasis(tperc,
                   lag=60,
                   argvar = list(fun="integer"),
                   arglag=list(fun = "ns", df=3)
  )
  
  
  
  cbvar = do.call("onebasis",c(list(x = 0:3), attr(cb,"argvar")))
  pre = crosspred(cbvar, coef = coef(metaResult)[1:3], vcov = vcov(metaResult)[1:3,1:3], 
                  at=0:3, cen = 0,
                  model.link = "log", cumul = T)
  
  rr.data = data.frame(
    Disease = dis[d],
    Event = c("Dry","Interval","Wet"),
    RR = pre[["allRRfit"]][-1],
    RRL = pre[["allRRlow"]][-1],
    RRH = pre[["allRRhigh"]][-1]
  )
  rr[[d]] = rr.data
  
  plot(pre, main = paste0(dis[d]), col = "red")
}

rr = do.call(rbind,rr)
write.xlsx(rr,"../Text_and_figures/2_DLNM/All.xlsx")


