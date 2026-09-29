
load("Work/Data_clean/RData/SWAP_allCell.RData")


search_day <- function(type, data, start_index){
  
  count <- 1
  i <- start_index
  q5 = quantile(data,0.05,na.rm=T)
  q95 = quantile(data,0.95,na.rm=T)
  
  ## abs values
  q5 = -0.5
  q95 = 0.5
  
  repeat{
    
    i <- i + 1
    
    if(i > length(data)){
      return(list(success=FALSE, stay=i))
    }
    
    v <- data[i]
    
    if(type == "drought"){
      if(v <= q5){ #-0.5
        count <- count + 1
      } else{
        return(list(success=FALSE, stay=i))
      }
    }
    
    if(type == "flood"){
      if(v >= q95){ #0.5
        count <- count + 1
      } else{
        return(list(success=FALSE, stay=i))
      }
    }
    
    if(count >= 10){
      return(list(success=TRUE, stay=i))
    }
    
  }
}
cal_single_event <- function(data, type){
  
  events <- list()
  q5 = quantile(data,0.05,na.rm=T)
  q95 = quantile(data,0.95,na.rm=T)
  
  ## abs values
  q5 = -0.5
  q95 = 0.5
  i <- 1
  mode <- ifelse(type=="gan","search drought start","search flood start")
  
  while(i <= length(data)){
    
    v <- data[i]
    
    if(mode=="search drought start"){
      
      if(v <= q5){ #-0.5
        
        res <- search_day("drought",data,i)
        
        if(res$success){
          
          start <- res$stay - 9
          mode <- "search drought end"
          i <- res$stay + 1
          
        } else{
          i <- res$stay
        }
        
      } else{
        i <- i + 1
      }
      
    }
    
    else if(mode=="search drought end"){
      
      if(v > q5){ #-0.5
        
        end <- i-1
        
        events[[length(events)+1]] <- list(
          start = start,
          end = end,
          continue = end-start+1,
          s = sum(data[start:end],na.rm=TRUE)
        )
        
        mode <- "search drought start"
        
      }
      
      i <- i+1
      
    }
    
    else if(mode=="search flood start"){
      
      if(v >= q95){ #0.5
        
        res <- search_day("flood",data,i)
        
        if(res$success){
          
          start <- res$stay - 9
          mode <- "search flood end"
          i <- res$stay + 1
          
        } else{
          i <- res$stay
        }
        
      } else{
        i <- i+1
      }
      
    }
    
    else if(mode=="search flood end"){
      
      if(v < q95){ #0.5
        
        end <- i-1
        
        events[[length(events)+1]] <- list(
          start=start,
          end=end,
          continue=end-start+1,
          s=sum(data[start:end],na.rm=TRUE)
        )
        
        mode <- "search flood start"
        
      }
      
      i <- i+1
      
    }
    
  }
  
  return(events)
}

merge_events <- function(events,data){
  
  repeat{
    
    merged <- FALSE
    
    for(i in 2:length(events)){
      
      last <- events[[i-1]]
      now  <- events[[i]]
      
      if(now$start - last$end -1 <=2){
        
        interval <- (last$end+1):(now$start-1)
        
        if(length(interval)==0){
          s_interval <- 0
        }else{
          s_interval <- sum(abs(1-abs(data[interval])))
        }
        
        if(s_interval / last$s <=0.2){
          
          new_event <- list(
            start=last$start,
            end=now$end,
            continue=now$end-last$start+1,
            interval=interval,
            s=last$s+now$s
          )
          
          events[[i-1]] <- new_event
          events[[i]] <- NULL
          
          merged <- TRUE
          break
          
        }
        
      }
      
    }
    
    if(!merged) break
    
  }
  
  return(events)
}
get_gan_to_shi <- function(gan,shi,data){
  
  res <- list()
  
  for(g in gan){
    
    for(s in shi){
      
      interval <- s$start - g$end -1
      
      if(interval<=5 & interval>=-5){
        
        intensity <- abs(
          (
            sum(data[s$start:(s$start+interval)]) -
              sum(data[(g$end-interval):g$end])
          )/interval
        )
        
        res[[length(res)+1]] <- data.frame(
          start=g$start,
          end=s$end,
          continue=s$end-g$start+1,
          breakpoint=g$end,
          interval=interval,
          intensity=intensity
        )
        
        break
      }
      
    }
    
  }
  
  return(res)
}
get_shi_to_gan <- function(gan_events, shi_events, data){
  
  results <- list()
  
  for(shi in shi_events){
    
    for(gan in gan_events){
      
      interval <- gan$start - shi$end - 1
      
      if(interval <= 5 & interval >= -5){
        
        intensity <- abs(
          (
            sum(data[gan$start:(gan$start+interval)], na.rm=TRUE) -
              sum(data[(shi$end-interval):shi$end], na.rm=TRUE)
          ) / interval
        )
        
        results[[length(results)+1]] <- data.frame(
          start = shi$start,
          end = gan$end,
          continue = gan$end - shi$start + 1,
          breakpoint = shi$end,
          interval=interval,
          intensity = intensity
        )
        
        break
      }
      
    }
    
  }
  
  return(results)
}


detect_fd <- function(swap){
  
  drought <- cal_single_event(swap,"gan")
  flood   <- cal_single_event(swap,"shi")
  
  if(length(drought)==0 | length(flood)==0){
    return(NULL)
  }
  
  drought <- merge_events(drought,swap)
  flood   <- merge_events(flood,swap)
  
  events <- get_gan_to_shi(drought,flood,swap)
  
  return(events)
}


### ..test
swap_point =  lstsave[[1]]
plot(swap_point, type="l")
abline(h=c(-0.5,0.5), col="red")
drought_events <- cal_single_event(swap_point, "gan")
flood_events <- cal_single_event(swap_point, "shi")
drought_events2 <- merge_events(drought_events, swap_point)
flood_events2   <- merge_events(flood_events, swap_point)
d2f <- get_gan_to_shi(drought_events2,
                            flood_events2,
                            swap_point)
f2d <- get_shi_to_gan(drought_events2,
                            flood_events2,
                            swap_point)



#### loop for each cell ####
load("Posted/Datasets/Derived/all_prec.RData")
is_notna = which(!is.na(lst[[1]]))
datetime = seq.Date(as.IDate("2000-01-01"),as.IDate("2025-12-31"),by="day")
datetime = datetime[!grepl("02-29",datetime)]
#' @lstsave 储存所有格点所有天的SWAP
#' @is_notna 所有cell的序号 

library(dplyr)
library(data.table)
myfun = function(i){
  
  swap_point = lstsave[[i]]
  drought_events <- cal_single_event(swap_point, "gan")
  flood_events <- cal_single_event(swap_point, "shi")
  drought_events2 <- merge_events(drought_events, swap_point)
  flood_events2   <- merge_events(flood_events, swap_point)
  d2f <- get_gan_to_shi(drought_events2,
                        flood_events2,
                        swap_point)
  f2d <- get_shi_to_gan(drought_events2,
                        flood_events2,
                        swap_point)
  save.d2f = data.frame(
    start = 0,
    end = 0,
    continue = 0,
    breakpoint = 0,
    interval = 0,
    intensity = 0,
    type = "d2f"
  )
  save.f2d = data.frame(
    start = 0,
    end = 0,
    continue = 0,
    breakpoint = 0,
    interval = 0,
    intensity = 0,
    type = "f2d"
  )
  if(length(d2f)!=0){
    save.d2f = do.call(rbind,d2f) %>% mutate(type = "d2f")
    # save.d2f$start = datetime[save.d2f$start]
    # save.d2f$end = datetime[save.d2f$end]
  }
  if(length(f2d)!=0){
    save.f2d = do.call(rbind,f2d) %>% mutate(type = "f2d")
    # save.f2d$start = datetime[save.f2d$start]
    # save.f2d$end = datetime[save.f2d$end]
  }
  
  
  saved = rbind(save.d2f,save.f2d)
  fwrite(saved,paste0("X:/Drought_Flood/Work/Data_clean/DAFF_abs/events/cell_",is_notna[i],".csv"))
 
}
library(parallel)
detectCores() #查看有多少核心
cl.cores <- 10  #使用的核心数
cl <- makeCluster(cl.cores)  #启用核心

clusterExport(cl, c("is_notna","datetime","lstsave","myfun","cal_single_event","merge_events",
                    "search_day","get_gan_to_shi","get_shi_to_gan"))#全局变量
#设置全局package
clusterEvalQ(cl, c(library(dplyr),
                   library(data.table))
)

x <- 1:length(is_notna)
system.time({
  parLapply(cl,x, myfun)
})
stopCluster(cl)



currf = list.files("X:/Drought_Flood/Work/Data_clean/DAFF_abs/events/")
currf = gsub("cell_","",currf)
currf = gsub(".csv","",currf)
currf = as.numeric(currf)
notyet = is_notna[is_notna %in% currf == F]
ii = which(is_notna %in% notyet)
# 
# for (i in 1:length(is_notna)) {
#   cat(sprintf("\rProcessed %d", i))
#   myfun(i)
# }

for (i in ii) {
  cat(sprintf("\rProcessed %d", i))
  myfun(i)
}



