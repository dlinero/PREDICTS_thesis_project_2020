

PlotErrBar_interactions_modi<-function(model
                                  , resp = ""
                                  , Effect1
                                  , Effect2
                                  , off= -0.4
                                  , off_increment = 0.2
                                  , pointtype = c(16,17,15)
                                  , legend = TRUE
                                  , leg.pos = "topright"
                                  , blackwhite = TRUE
                                  , ylims = c(-5,3)
                                  , srttxt =45
){ 
  
  # I've written this assuming that Effect 1 is the land-use, and Effect 2 is the interacting factor.
  # Personally, I like colours and points to refer to the second variable, not land-use, so that's how I've written the function.
  
  require(roquefort)
  library(RColorBrewer)
  
  exp1levs<-levels(model@frame[,Effect1])
  exp2levs<-levels(model@frame[,Effect2])
  
  ## set up pch and colours for plotting
  if(missing(pointtype)){
    pointtype <- rep(16,length(exp2levs))
  }
  
  if(blackwhite == "FALSE"){
    if(length(exp2levs)>2){
      cols<-brewer.pal(n = length(exp2levs),name = "Set1")
    }else{
      cols<-brewer.pal(n = 3,name = "Set1")[1:2]
    }
  }else{
    cols = rep("#000000", length(exp2levs))
  }
  
  #Plot the reference levels with offset
  PlotErrBarADP(model = model,data = model@frame, 
                responseVar = resp, logLink = "n", 
                catEffects = Effect1
                , offset = off
                , forPaper = TRUE
                , errbar.cols = cols[1]
                , ylim = ylims
                , pt.pch = pointtype[1])
  
  
  # extract variance covariance matrix
  var.cov<-vcov(model)
  
  # start off counter for point placement
  o<-off+off_increment
  for(i in 2:length(exp2levs)){ # for each level of your interacting factor (except the reference level)
    
    ref<-fixef(model)[paste(Effect2,exp2levs[i], sep = "")] # get the main coefficient for Effect 2
    
    # plot ref level rescaled to zero
    errbar(1+o,ref-ref,ref-ref,ref-ref,
           add=TRUE,col=cols[i],errbar.col=cols[i],pch=pointtype[i])
    
    for(s in 2:length(exp1levs)){# for each level of habuse (except for the reference level)
      
      # get the coef by adding the main land use effect, the main trait effect, and the interaction
      coef<-fixef(model)[paste(Effect1,exp1levs[s], sep = "")]+
        fixef(model)[paste(Effect2,exp2levs[i], sep = "")]+
        fixef(model)[paste(Effect1,exp1levs[s],":",Effect2,exp2levs[i], sep = "")]
      
      cov<-var.cov[which(row.names(var.cov)==paste(Effect1,exp1levs[s], sep = "")),
                   which(names(var.cov[1,])==paste(Effect1,exp1levs[s],":",Effect2,exp2levs[i], sep = ""))]
      
      
      # extract standard errors
      se1<-se.fixef(model)[which(names(fixef(model))==paste(Effect1,exp1levs[s], sep = ""))]
      se2<-se.fixef(model)[which(names(fixef(model))==paste(Effect1,exp1levs[s],":",Effect2,exp2levs[i], sep = ""))]
      
      ## create standard error for plotting
      se.plot<-sqrt(se1^2+se2^2+2*cov)
      
      errbar(s+o,coef-ref,coef+1.96*se.plot-ref,coef-1.96*se.plot-ref,
             add=TRUE,col=cols[i],errbar.col=cols[i],pch = pointtype[i])
    }
    # increase counter for point placement
    o<-o+off_increment
    
  }
  
  if(legend == TRUE){
    legend(leg.pos, legend = exp2levs, col = cols, title = Effect2, lty = 1, pch = pointtype)
  }

}
