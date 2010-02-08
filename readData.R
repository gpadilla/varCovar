readData <- function(filePath="28of28/",pattern=".txt$"){

  fn <- list.files(filePath, pattern=pattern)

  noOfFiles <- length(fn)



  for (k in 1:noOfFiles){

    tmpDat <- read.table(paste(filePath,fn[k],sep=""), skip=1)

    if(k ==1){
      ## data matrix
      datAlpha <- matrix(ncol=length(tmpDat[,3]), nrow=noOfFiles)
      datNu <- matrix(ncol=length(tmpDat[,1]), nrow=noOfFiles)

      datAlpha[k,] <- tmpDat[,3]
      datNu[k,] <-  tmpDat[,1]

      ## Unsicherheitsmatrix
      unsAlpha <- matrix(ncol=length(tmpDat[,4]),nrow=noOfFiles)
      unsNu    <-   matrix(ncol=length(tmpDat[,2]), nrow=noOfFiles)
      unsAlpha[k,] <- tmpDat[,4]
      unsNu[k,] <- tmpDat[,2]

    }else{
      datAlpha[k,] <- tmpDat[,3]
      datNu[k,]    <- tmpDat[,1]
      unsAlpha[k,] <- tmpDat[,4]
      unsNu[k,]        <- tmpDat[,2]
    }

  }

  return( list(datAlpha=datAlpha,
               datNu=datNu,
               unsAlpha=unsAlpha,
               unsNu=unsNu)
         )
}

centerCurve <- function(resList){

  noOfQuadradicFitPoints <- 100



  noOfScans <- length(resList$datAlpha[,1])
  noOfPoints <- length(resList$datAlpha[1,])

  maxIdxAlphaVec <- rep(NA,noOfScans)
  maxNuVec <- rep(NA,noOfScans)
  ## sammeln der Informat. zur Zentrierung
  for(i in 1:noOfScans){

    maxIdx <- which.max(resList$datAlpha[i,])

    ax <- (maxIdx - noOfQuadradicFitPoints):(maxIdx + noOfQuadradicFitPoints)
    ay <- resList$datAlpha[i,ax]

    g2 <- lm(ay ~ ax +I(ax^2))
    ## der mittelwert der 2 nullstellen ist das exakte maximum
    maxIdxAlphaVec[i] <- round(mean(Re(polyroot(c(g2$coef[1],g2$coef[2],g2$coef[3])))),0)
    maxNuVec[i] <- resList$datNu[i,maxIdxAlphaVec[i]]
  }
  centerPos <- round(mean(maxIdxAlphaVec),0)
  shiftVec <-  centerPos - maxIdxAlphaVec

  ## resIdx ist so lang, dass alle zugeschnittenen alphas platz haben

  reducedVecLength <- noOfPoints+min(shiftVec) - max(shiftVec)

  resIdx <- seq(1,reducedVecLength)

  datAlpha <- matrix(ncol=reducedVecLength, nrow=noOfScans)
  datNu    <- matrix(ncol=reducedVecLength, nrow=noOfScans)
  unsAlpha <- matrix(ncol=reducedVecLength, nrow=noOfScans)
  unsNu    <- matrix(ncol=reducedVecLength, nrow=noOfScans)

  ## zentrieren der Einzelkurven
  for(i in 1:noOfScans){

    if(shiftVec[i] > 0){
      l <- shiftVec[i] + 1
      m <- shiftVec[i] + reducedVecLength
      print(c(l,m))
      datAlpha[i,] <- resList$datAlpha[i, l:m]
    }else{
      #
      l <- ((noOfPoints+shiftVec[i]) -reducedVecLength) +1
      m <- (noOfPoints+shiftVec[i])

      datAlpha[i,] <-resList$datAlpha[i, l:m]
    }

  }


  return( list(datAlpha=datAlpha,
               datNu=datNu,
               unsAlpha=unsAlpha,
               unsNu=unsNu,
               shiftVec=shiftVec)
         )
}



## ==============================
calVarCovar <- function(datMat, mvMat){

noOfScans <-length(datMat[,1])
lengthOfScan <-length(datMat[1,])

resMat<- matrix(ncol=lengthOfScan ,nrow=lengthOfScan )

for(i in 1:lengthOfScan ){
    for(k in 1:lengthOfScan ){

      resMat[i,k] <- sum( (datMat[,i] - mvMat[i])*(datMat[,k] - mvMat[k]))/noOfScans
    }
  }
return(resMat)
}

calSdMv <- function(datMat, unsMat){

  noOfScans <-length(datMat[,1])
  lengthOfScan <-length(datMat[1,])


  mvMat <- NULL
  sdMat <- NULL
  mvUnsMat <- NULL
  for(l in 1:lengthOfScan){
    mvMat[l] <- mean(datMat[,l])
    mvUnsMat[l]<- mean(unsMat[,l])
    sdMat[l] <- sd(datMat[,l])
  }

return(list(mvMat = mvMat,
            sdMat= sdMat,
            mvUnsMat=mvUnsMat
            ))
}

reduceBy <- 10

resList <- readData()
idx <- seq(1,length(resList$datMat[1,]), reduceBy)

datMat <-  resList$datMat[,idx]
unsMat <-  resList$unsMat[,idx]

msd <- calSdMv(datMat, unsMat)
mvMat <- msd$mvMat
sdMat <- msd$sdMat
unsMvMat <- msd$mvUnsMat

ua <- ((sdMat)^2)^.5 # +  (unsMvMat)^2)^.5

uaMat <- outer(ua,ua,function(x,y){return(x*y)})

res <- calVarCovar(datMat, mvMat)/uaMat

colS <-  heat.colors(10)

image(idx,idx,res,col=colS,
      xlab=expression(paste( nu,", ",alpha)),
      ylab=expression(paste( nu,", ",alpha)))
##contour(idx,idx,res, add=TRUE)
legend(.0,idx[length(idx)],
       round(seq(min(res), max(res),(max(res) -min(res))/(reduceBy -1)),3),
       col=colS,
       lty=rep(1,reduceBy),
       lwd=10)

