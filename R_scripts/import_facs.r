
library(flowCore)
library(data.table)


autoLgclmod <- function(x, channels, m = 4.5, q = 0.05) {
  
  if (missing(channels)) 
    stop("Please specify the channels to be logicle transformed")
  indx <- channels %in% colnames(x)
  if (!all(indx)) 
    stop(paste("Channels", channels[!indx], "were not found in the FCS file.\n ", 
               sep = " "))
  
  trans <- lapply(channels, function(p) {
    data <- x[, p][!(is.na(x[, p]))]
    w <- 0
    t <- max(data)
    ndata <- data[data < 0]
    ## use 1.5 * IQR to filter outliers in negative values
    nThres <- quantile(ndata, 0.25) - 1.5 * IQR(ndata)
    ndata <- ndata[ndata >= nThres]
    transId <- paste(p, "autolgclTransform", sep = "_")
    
    if (length(ndata)) {
      r <- .Machine$double.eps + quantile(ndata, q)
      ## Check to avoid failure of negative w
      if (10^m * abs(r) <= t) {
        w <- 0  
      } else {
        w <- (m - log10(t/abs(r)))/2
        if(is.nan(w) || w>2) {
          warning(paste0("autoLgcl failed for channel: ", p, "; using default logicle transformation!"))
          w <- 0.1
          t <- 4000 
          m <- 4.5 
        }
      }
    }
    logicleTransform(transformationId = transId, 
                     w = w, t = t, m = m, a = 0)
  })
  transformList(channels, trans)
}


scaleData <- function(x, range = c(0, 4.5)) {
  #(x - min(x))/(max(x) - min(x)) * (range[2] - range[1]) + range[1]
  (x - 0)/(max(x,na.rm = TRUE) - 0) * (range[2] - range[1]) + range[1]
}





mergeAutoLgcl <- function(files,
                           pick = "all",
                           fileIDs = c(),
                           stainStatus = c(),
                           flowGroup,
                           sampleIDs=c()
                           ){
  q = 0.05
  if(pick == "all"){
    allEvents = TRUE
  }else{allEvents = FALSE}
  
  flowSets <- list()
  i=1
  for(sp in unique(sampleIDs)){
    print(sp)
    
    flowSets[[i]] <- read.flowSet(files[sampleIDs == sp],transformation = NULL)
    
    i=i+1
  }
  
  transformedFs <- list()
  j=1
  for(sp in unique(sampleIDs)){
    print(j)
    data <- data.frame()
    for( i in c(1:length(files[sampleIDs == sp]))){
      data <- rbind(data,flowSets[[j]][[i]]@exprs)
    }
    
    marker_id <- seq_along(colnames(data))
    
    size_channels <- grep("FSC|SSC", colnames(data), ignore.case = TRUE)
    transMarker_id <- setdiff(marker_id, size_channels)
    
    
    trans <- autoLgclmod(data, channels = colnames(data)[transMarker_id], q = q)
    
    transformedFs[[j]] <- transform(flowSets[[j]],trans)
    j=j+1
    
    
    
    
  }
  
  dataTransList <- list()
  j=1
  for(sp in unique(sampleIDs)){
    
    data <- data.frame()
    for( i in c(1:length(files[sampleIDs == sp]))){
      fsdata <- as.data.frame(transformedFs[[j]][[i]]@exprs)
      fsdata$stainInfo <- rep(stainStatus[sampleIDs == sp][i],nrow(fsdata))
      fsdata$fileID <- fileIDs[sampleIDs == sp][i]
      fsdata$sampleID <- sp
      fsdata$groupID <- flowGroup[sampleIDs == sp][[1]]
      data <- rbind(data,fsdata)
      
    }
    
    dataTransList[[j]] <- data
    
    j=j+1
  }
  
  rm(transformedFs)
  
  dataTrans <- data.table::rbindlist(dataTransList,fill=TRUE)
  dataTrans <- as.data.frame(dataTrans)
  
  #exclude_channels <- grep("Time|Event", colnames(dataTrans), ignore.case = TRUE)
  marker_id <- seq_along(colnames(dataTrans))
  
  size_channels <- grep("FSC|SSC", colnames(dataTrans), ignore.case = TRUE)
  #transMarker_id <- setdiff(marker_id, size_channels)
  
  ##apply linear transformation to scatter channels
  if(length(size_channels) > 0){
    if(any(size_channels %in% marker_id)){
      used_size_channel <- size_channels[size_channels %in% marker_id]
      #used_size_channel_id <- match(used_size_channel, marker_id)
      dataTrans[ ,used_size_channel] <- apply(dataTrans[ , used_size_channel, drop=FALSE], 2, 
                                              function(x) scaleData(x, range=c(0, 4.5)))
    }
  }
  
  dataTrans <- dataTrans[,!(colnames(dataTrans) %in% c("Time"))]
  
  if(isFALSE(allEvents)){
    data.picked <- data.frame()
    for(i in 1:(length(fileIDs)/2)){
      dataSample <- dataTrans[dataTrans$sampleID %in% c(fileIDs[i],fileIDs[i+length(fileIDs)/2]),]
      
      dataS <- dataSample[dataSample$stainInfo == "stained",]
      dataU <- dataSample[dataSample$stainInfo == "unstained",]
      
      
      dataS <- dataS[sample(nrow(dataS),pick),]
      
      dataU <- dataU[sample(nrow(dataU),pick),]
      data.picked <- rbind(data.picked,dataS,dataU)
    }
    dataTrans <- data.picked
  }
  
  dataTrans <- split(dataTrans,dataTrans$groupID)
  return(dataTrans)
}




