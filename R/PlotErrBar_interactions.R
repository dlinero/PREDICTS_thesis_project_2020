PlotErrBarADP<-function(model,data,responseVar,seMultiplier=1.96,outDir=NULL,
                     logLink="n",catEffects=NULL,
                     contEffects=list(),contEffectsLabels=NULL,
                     otherCatEffects=list(),
                     forPaper=FALSE,align=FALSE,secdAge=FALSE,
                     xtext.srt=0,ylim=NA,order=NULL,rescale=NULL,
                     errbar.cols=NULL,pt.pch=NULL,
                     params=list(),add=FALSE,offset=0){
  
  if (!is.null(contEffectsLabels)){
    if (length(contEffects)!=length(contEffectsLabels)){
      stop("Labels for continuous effects must be of the same number as the effects")
    }
  }
  
  # For certain recognized factors, define the default ordering and colour
  predefined.factors<-list()
  
  if (secdAge){
    predefined.factors$UI<-data.frame(UI=c("Primary Vegetation Minimal use",
                                           "Primary Vegetation Light use",
                                           "Primary Vegetation Intense use",
                                           "Mature Secondary Vegetation Minimal use",
                                           "Mature Secondary Vegetation Light use",
                                           "Intermediate Secondary Vegetation Minimal use",
                                           "Intermediate Secondary Vegetation Light use",
                                           "Young Secondary Vegetation Minimal use",
                                           "Young Secondary Vegetation Light use",
                                           "Plantation forest Minimal use",
                                           "Plantation forest Light use",
                                           "Plantation forest Intense use",
                                           "Cropland Minimal use",
                                           "Cropland Light use",
                                           "Cropland Intense use",
                                           "Pasture Minimal use",
                                           "Pasture Light use",
                                           "Pasture Intense use",
                                           "Urban Minimal use",
                                           "Urban Light use",
                                           "Urban Intense use"),
                                      col=c(rep("#66A61E",3),rep("#147659",2),
                                            rep("#1B9E77",2),rep("#8ecfbc",2),
                                            rep("#7570B3",3),rep("#E6AB02",3),
                                            rep("#D95F02",3),rep("#E7298A",3)))
    predefined.factors$LandUse<-data.frame(LandUse=c("Primary Vegetation",
                                                     "Mature Secondary Vegetation",
                                                     "Intermediate Secondary Vegetation",
                                                     "Young Secondary Vegetation",
                                                     "Plantation forest",
                                                     "Cropland",
                                                     "Pasture",
                                                     "Urban"),
                                           col=c("#66A61E","#147659","#1B9E77","#8ecfbc",
                                                 "#7570B3","#E6AB02","#D95F02","#E7298A"))
  } else {
    predefined.factors$UI<-data.frame(UI=c("Primary Vegetation Minimal use",
                                           "Primary Vegetation Light use",
                                           "Primary Vegetation Intense use",
                                           "Secondary Vegetation Minimal use",
                                           "Secondary Vegetation Light use",
                                           "Secondary Vegetation Intense use",
                                           "Plantation forest Minimal use",
                                           "Plantation forest Light use",
                                           "Plantation forest Intense use",
                                           "Cropland Minimal use",
                                           "Cropland Light use",
                                           "Cropland Intense use",
                                           "Pasture Minimal use",
                                           "Pasture Light use",
                                           "Pasture Intense use",
                                           "Urban Minimal use",
                                           "Urban Light use",
                                           "Urban Intense use"),
                                      col=c(rep("#66A61E",3),rep("#1B9E77",3),
                                            rep("#7570B3",3),rep("#E6AB02",3),
                                            rep("#D95F02",3),rep("#E7298A",3)))
    predefined.factors$LandUse<-data.frame(LandUse=c("Primary vegetation",
                                                     "Secondary vegetaion",
                                                     "Secondary Vegetation",
                                                     "Null Secondary Vegetation",
                                                     "Plantation forest",
                                                     "1st generation",
                                                     "2nd generation",
                                                     "Pasture",
                                                     "Urban"),
                                           col=c("#66A61E","#1B9E77","#1B9E77","#1B9E77","#7570B3",
                                                 "#E6AB02","#E6AB02","#D95F02","#E7298A"))
  }
  predefined.factors$UseIntensity<-data.frame(UseIntensity=c("Minimal use","Light use",
                                                             "Intense use"),
                                              col=c("#336600","#006666","#CC0099"))
  
  labels<-character(0)
  coef.labels<-character(0)
  # Get the names of the labels and coefficient names for each factor
  # If align is specified, then use the full set of labels defined above
  for (e in catEffects){
    if (align & (e %in% names(predefined.factors))){
      eval(substitute(names<-paste(predefined.factors$x$x),list(x=e)))
      labels<-c(labels,names)
      coef.labels<-c(coef.labels,paste(e,names,sep=""))
    } else {
      eval(substitute(names<-levels(model@frame[,e])[order(match(tolower(
        levels(model@frame[,e])),tolower(predefined.factors$x$x)))],list(x=e)))
      labels<-c(labels,names)
      coef.labels<-c(coef.labels,paste(e,names,sep=""))
    }
  }
  
  # Get coefficient and standard error estimates from the model
  o<-match(tolower(coef.labels),tolower(names(fixef(model))))
  y<-fixef(model)[o]
  yplus<-y+se.fixef(model)[o]*seMultiplier
  yminus<-y-se.fixef(model)[o]*seMultiplier
  
  # For each categorical effect, get the reference level
  for (e in catEffects){
    ref.name<-paste(e,levels(model@frame[,e])[1],sep="")
    o<-match(tolower(ref.name),tolower(coef.labels))
    y[o]<-0
    yplus[o]<-0
    yminus[o]<-0
    
  }
  
  if (!is.null(order)){
    labels<-labels[order]
    coef.labels<-coef.labels[order]
    y<-y[order]
    yplus<-yplus[order]
    yminus<-yminus[order]
  }
  
  if(!is.null(rescale)){
    y<-y+rescale
    yplus<-yplus+rescale
    yminus<-yminus+rescale
  }
  
  if (logLink=="e"){
    y<-(exp(y)*100)-100
    yplus<-(exp(yplus)*100)-100
    yminus<-(exp(yminus)*100)-100
  } else if (logLink=="10") {
    y<-(10^(y)*100)-100
    yplus<-(10^(yplus)*100)-100
    yminus<-(10^(yminus)*100)-100
  } else if (logLink=="b"){
    y<-((1/(1+exp(-(y))))*200)-100
    yplus<-((1/(1+exp(-(yplus))))*200)-100
    yminus<-((1/(1+exp(-(yminus))))*200)-100
  } else if (logLink=="n"){
    
  } else {
    stop("Error: the specified log link is not supported")
  }
  
  
  if (length(contEffects)>0){
    for (i in 1:length(contEffects)){
      # Add this continuous effect to the labels
      if (is.null(contEffectsLabels)){
        labels<-c(labels,"",names(contEffects)[i],"")
      } else {
        labels<-c(labels,"",contEffectsLabels[i],"")
      }
      
      
      # Make a data frame with the minimum and maximum values onto which to project the model
      if (contEffects[i]==0){
        
      } else if (contEffects[i]==1){
        eval(substitute(newdat<-data.frame(c(min(data$x),median(data$x),max(data$x))),list(x=names(contEffects)[i])))
      } else if (contEffects[i]==2){
        eval(substitute(newdat<-data.frame(c(max(data$x),median(data$x),min(data$x))),list(x=names(contEffects)[i])))
      } else {
        stop("Continuous effect order must be either 0 (plot nothing), 1 (ascending) or 2 (descending)")
      }
      
      
      if (contEffects[i]==0){
        y<-c(y,rep(NA,3))
        yplus<-c(yplus,rep(NA,3))
        yminus<-c(yminus,rep(NA,3))
      } else {
        
        
        names(newdat)<-names(contEffects)[i]
        
        # Loop over other continuous effects and add median values
        for (e in names(contEffects)){
          if (e != names(contEffects)[i]){
            newdat[,e]<-median(data[,e])
          }
        }
        
        for (e in 1:length(otherCatEffects)){
          newdat[,names(otherCatEffects)[e]]<-factor(
            otherCatEffects[e],levels=levels(model@frame[,names(otherCatEffects)[e]]))
        }
        
        # Loop over factors and add reference values
        for (e in catEffects){
          newdat[,e]<-factor(levels(model@frame[,e])[1],levels=levels(model@frame[,e]))
        }
        
        newdat[,names(model@frame)[1]] <- 0
        
        mm<-model.matrix(terms(model),newdat)
        pvar1 <- diag(mm %*% base::tcrossprod(as.matrix(vcov(model)),mm))
        newdat$y<-mm %*% fixef(model)
        newdat$yplus<-newdat$y+sqrt(pvar1)*seMultiplier
        newdat$yminus<-newdat$y-sqrt(pvar1)*seMultiplier
        
        if (logLink=="e"){
          temp.y<-c(0,(((exp(newdat$y[2])/exp(newdat$y[1]))*100)-100),(((exp(newdat$y[3])/exp(newdat$y[1]))*100)-100))
          temp.yplus<-c(0,(((exp(newdat$yplus[2])/exp(newdat$y[1]))*100)-100),(((exp(newdat$yplus[3])/exp(newdat$y[1]))*100)-100))
          temp.yminus<-c(0,(((exp(newdat$yminus[2])/exp(newdat$y[1]))*100)-100),(((exp(newdat$yminus[3])/exp(newdat$y[1]))*100)-100))
        } else if (logLink=="10"){
          temp.y<-c(0,(((10^(newdat$y[2])/10^(newdat$y[1]))*100)-100),(((10^(newdat$y[3])/10^(newdat$y[1]))*100)-100))
          temp.yplus<-c(0,(((10^(newdat$yplus[2])/10^(newdat$y[1]))*100)-100),(((10^(newdat$yplus[3])/10^(newdat$y[1]))*100)-100))
          temp.yminus<-c(0,(((10^(newdat$yminus[2])/10^(newdat$y[1]))*100)-100),(((10^(newdat$yminus[3])/10^(newdat$y[1]))*100)-100))
        } else if (logLink=="b"){
          y<-((1/(1+exp(-(y))))*200)-100
          yplus<-((1/(1+exp(-(yplus))))*200)-100
          yminus<-((1/(1+exp(-(yminus))))*200)-100
        } else if (logLink=="n") {
          temp.y<-c(0,(newdat$y[2]-newdat$y[1]),(newdat$y[3]-newdat$y[1]))
          temp.yplus<-c(0,(newdat$yplus[2]-newdat$y[1]),(newdat$yplus[3]-newdat$y[1]))
          temp.yminus<-c(0,(newdat$yminus[2]-newdat$y[1]),(newdat$yminus[3]-newdat$y[1]))
        } else {
          stop("Error: the specified log link is not supported")
        }
        
        
        y<-c(y,temp.y)
        yplus<-c(yplus,temp.yplus)
        yminus<-c(yminus,temp.yminus)
        
        
      }
      
    }
    
  }
  
  predRange<-max(yplus,na.rm=T)-min(yminus,na.rm=T)
  if (all(is.na(ylim))){
    plotLims<-c(min(yminus,na.rm=T)-0.35*predRange,max(yplus,na.rm=T)+0.05*predRange)
  } else {
    plotLims<-ylim
    predRange<-ylim[2]-ylim[1]
  }
  
  
  if (!is.null(outDir)){
    png(paste(outDir,"/LU effects_",responseVar,".png",sep=""),width=22.86/2.54,
        height=12.57/2.54,units="in",res=1200)
  }
  
  if(forPaper){
    par(mar=c(0.2,3.5,0.2,0.2))
    par(cex.lab=1)
    par(cex.axis=0.7)
    cex.pt<-0.5
    par(mgp=c(2,1,0))
    cex.txt<-0.5
  } else {
    par(mar=c(1,5,1,1))
    par(cex.lab=1.8)
    par(cex.axis=1.5)
    cex.pt<-1
    par(mgp=c(3.5,1,0))
    cex.txt<-0.8
  }    
  par(las=1)
  
  for (p in names(params)){
    eval(parse(text=gsub("x",p,"par(x=params$x)")))
  }
  
  if (responseVar != ""){
    if (logLink=="n"){
      ylabel=paste("Relative",responseVar)
    } else {
      ylabel=paste(responseVar,"change (%)")
    }
    
  } else {
    ylabel<-""
  }
  
  allLevels<-unlist(lapply(predefined.factors,function(x) return(x[,1])))
  allCols<-unlist(lapply(predefined.factors,function(x) return(x[,2])))
  
  plot.cols<-paste(allCols[match(tolower(labels),tolower(allLevels))])
  plot.cols[plot.cols=="NA"]<-"black"
  
  if (!is.null(errbar.cols)){
    plot.cols<-errbar.cols
  }
  
  if(!add){
    plot(-1,-1,xlim=c(0.5,length(labels)+0.5),ylim=plotLims,xlab=NA,ylab=
           ylabel,xaxt="n",bty="n")
  }
  
  
  
  for (i in 1:length(catEffects)){
    if ((i!=length(catEffects)) | (length(contEffects)>0)){
      abline(v=max(grep(catEffects[i],coef.labels))+0.5,lwd=1,
             col="dark grey")
    }
    
  }
  if (any(contEffects>0)){
    if (length(contEffects)>0){
      for (i in 1:length(contEffects)){
        if (i != length(contEffects)){
          if (is.null(contEffectsLabels)){
            abline(v=which(labels==names(contEffects)[i])+1.5,lwd=1,
                   col="dark grey")
          } else {
            abline(v=which(labels==contEffectsLabels[i])+1.5,lwd=1,
                   col="dark grey")
          }
          
        }
      }
      
    }
    
  }
  
  # if (!add){
  #   if (TRUE %in% grepl("primary",tolower(labels))){
  #     if(!forPaper){
  #       rect(min(grep("primary",tolower(labels)))-0.5,plotLims[1],
  #            max(grep("primary",tolower(labels)))+0.5,plotLims[2],
  #            col="#66A61E33",border=NA)
  #     }
  #     text(mean(grep("primary",tolower(labels))),plotLims[1]+0.05*predRange,
  #          "Primary",col="black",cex=cex.txt,srt=xtext.srt,xpd=TRUE)
  #   }
  #   if (secdAge){
  #     if (TRUE %in% grepl("mature secondary",tolower(labels))){
  #       if (!forPaper){
  #         rect(min(grep("mature secondary",tolower(labels)))-0.5,plotLims[1],
  #              max(grep("mature secondary",tolower(labels)))+0.5,plotLims[2],
  #              col="#14765933",border=NA)
  #       }
  #       text(mean(grep("mature secondary",tolower(labels))),plotLims[1]+0.05*
  #              predRange,"MSV",col="black",cex=cex.txt,srt=xtext.srt,xpd=TRUE)
  #     }
  #     if (TRUE %in% grepl("intermediate secondary",tolower(labels))){
  #       if (!forPaper){
  #         rect(min(grep("intermediate secondary",tolower(labels)))-0.5,plotLims[1],
  #              max(grep("intermediate secondary",tolower(labels)))+0.5,plotLims[2],
  #              col="#1B9E7733",border=NA)
  #       }
  #       text(mean(grep("intermediate secondary",tolower(labels))),plotLims[1]+0.05*
  #              predRange,"ISV",col="black",cex=cex.txt,srt=xtext.srt,xpd=TRUE)
  #     }
  #     if (TRUE %in% grepl("young secondary",tolower(labels))){
  #       if (!forPaper){
  #         rect(min(grep("young secondary",tolower(labels)))-0.5,plotLims[1],
  #              max(grep("young secondary",tolower(labels)))+0.5,plotLims[2],
  #              col="#8ecfbc33",border=NA)
  #       }
  #       text(mean(grep("young secondary",tolower(labels))),plotLims[1]+0.05*
  #              predRange,"YSV",col="black",cex=cex.txt,srt=xtext.srt,xpd=TRUE)
  #     }
  #     
  #   } else {
  #     if (TRUE %in% grepl("secondary",tolower(labels))){
  #       if (!forPaper){
  #         rect(min(grep("secondary",tolower(labels)))-0.5,plotLims[1],
  #              max(grep("secondary",tolower(labels)))+0.5,plotLims[2],
  #              col="#1B9E7733",border=NA)
  #       }
  #       text(mean(grep("secondary",tolower(labels))),plotLims[1]+0.05*predRange,
  #            "Secondary",col="black",cex=cex.txt,srt=xtext.srt,xpd=TRUE)
  #     }
  #     
  #   }
  #   if (TRUE %in% grepl("plantation",tolower(labels))){
  #     if (!forPaper){
  #       rect(min(grep("plantation",tolower(labels)))-0.5,plotLims[1],
  #            max(grep("plantation",tolower(labels)))+0.5,plotLims[2],
  #            col="#7570B333",border=NA)
  #     }
  #     text(mean(grep("plantation",tolower(labels))),plotLims[1]+0.05*predRange,
  #          "Plantation",col="black",cex=cex.txt,srt=xtext.srt,xpd=TRUE)
  #   }
  #   if (TRUE %in% grepl("cropland",tolower(labels))){
  #     if (!forPaper){
  #       rect(min(grep("cropland",tolower(labels)))-0.5,plotLims[1],
  #            max(grep("cropland",tolower(labels)))+0.5,plotLims[2],
  #            col="#E6AB0233",border=NA)
  #     }
  #     text(mean(grep("cropland",tolower(labels))),plotLims[1]+0.05*predRange,
  #          "Cropland",col="black",cex=cex.txt,srt=xtext.srt,xpd=TRUE)
  #   }
  #   if (TRUE %in% grepl("pasture",tolower(labels))){
  #     if (!forPaper){
  #       rect(min(grep("pasture",tolower(labels)))-0.5,plotLims[1],
  #            max(grep("pasture",tolower(labels)))+0.5,plotLims[2],
  #            col="#D95F0233",border=NA)
  #     }
  #     text(mean(grep("pasture",tolower(labels))),plotLims[1]+0.05*predRange,
  #          "Pasture",col="black",cex=cex.txt,srt=xtext.srt,xpd=TRUE)
  #   }
  #   if (TRUE %in% grepl("urban",tolower(labels))){
  #     if (!forPaper){
  #       rect(min(grep("urban",tolower(labels)))-0.5,plotLims[1],
  #            max(grep("urban",tolower(labels)))+0.5,plotLims[2],
  #            col="#E7298A33",border=NA)
  #     }
  #     text(mean(grep("urban",tolower(labels))),plotLims[1]+0.05*predRange,
  #          "Urban",col="black",cex=cex.txt,srt=xtext.srt,xpd=TRUE)
  #   }
  #   if (any(regexpr("Minimal",labels)==1)){
  #     if (!forPaper){
  #       rect(min(which(regexpr("Minimal",labels)==1))-0.5,plotLims[1],
  #            max(which(regexpr("Minimal",labels)==1))+0.5,plotLims[2],
  #            col="#33660033",border=NA)
  #     }
  #     text(mean(grep("Minimal",labels)),plotLims[1]+0.05*predRange,
  #          "Minimal",col="black",cex=cex.txt,srt=xtext.srt,xpd=TRUE)
  #   }
  #   if (any(regexpr("Light",labels)==1)){
  #     if (!forPaper){
  #       rect(min(which(regexpr("Light",labels)==1))-0.5,plotLims[1],
  #            max(which(regexpr("Light",labels)==1))+0.5,plotLims[2],
  #            col="#00666633",border=NA)
  #     }
  #     text(mean(grep("Light",labels)),plotLims[1]+0.05*predRange,
  #          "Light",col="black",cex=cex.txt,srt=xtext.srt,xpd=TRUE)
  #   }
  #   if (any(regexpr("Intense",labels)==1)){
  #     if (!forPaper){
  #       rect(min(which(regexpr("Intense",labels)==1))-0.5,plotLims[1],
  #            max(which(regexpr("Intense",labels)==1))+0.5,plotLims[2],
  #            col="#CC009933",border=NA)
  #     }
  #     text(mean(grep("Intense",labels)),plotLims[1]+0.05*predRange,
  #          "Intense",col="black",cex=cex.txt,srt=xtext.srt,xpd=TRUE)
  #   }
  #   
  #   if (length(which(plot.cols=="black"))>0){
  #     text(which(plot.cols=="black"),plotLims[1]+0.05*predRange,
  #          labels[which(plot.cols=="black")],cex=cex.txt,srt=xtext.srt,xpd=TRUE)
  #     
  #   }
    
    if(logLink=="n"){
      abline(h=0,col="dark grey")
    } else {
      abline(h=0,col="dark grey")
    }
    
  #}
  
  #   if (("UseIntensity" %in% names(model@frame)) | ("UI" %in% names(model@frame))){
  #     intensities<-character(length(labels))
  #     intensities[which(grepl("Minimal",labels))]<-"Min"
  #     intensities[which(grepl("Light",labels))]<-"Lt"
  #     intensities[which(grepl("Intense",labels))]<-"Int"
  #     
  #     if (secdAge) intensities[(intensities=="Lt") & (grepl("Secondary",labels))]<-"L-I"
  #     
  #     if ((logLink=="10") | (logLink=="e")){
  #       text(1:length(intensities),yminus-(predRange/11),intensities,col=
  #              "black",cex=cex.txt,srt=45)
  #     } else {
  #       text(1:length(intensities),yminus-(predRange/11),intensities,col=
  #              "black",cex=cex.txt,srt=45)
  #     }
  #     
  #   }
  #   
  if (length(contEffects>0)){
    if ((logLink=="10") | (logLink=="e")){
      text((which(labels=="")[1]:(which(labels=="")[1]+
                                    length(contEffects)*3))+offset,
           yminus[which(labels=="")[1]:(which(labels=="")[1]+
                                          length(contEffects)*3)]-(predRange/11),
           rep(c("Lo","Med","Hi"),length(contEffects)),
           col="black",cex=cex.txt,srt=45)
    } else {
      text((which(labels=="")[1]:(which(labels=="")[1]+
                                    length(contEffects)*3))+offset,
           yminus[which(labels=="")[1]:(which(labels=="")[1]+
                                          length(contEffects)*3)]-(predRange/11),
           rep(c("Lo","Med","Hi"),length(contEffects)),
           col="black",cex=cex.txt,srt=45)
    }
    
  }
  
  if(is.null(pt.pch)){
    pt.pch<-16
  }
  
  errbar((1:length(labels))+offset,y,yplus,yminus,col=plot.cols,errbar.col=plot.cols,add=T,pch=pt.pch)
  
  if (!is.null(outDir)){
    dev.off()
  }
}

errbar<-function (x, y, yplus, yminus, cap = 0.015, main = NULL, sub = NULL, 
          xlab = as.character(substitute(x)), ylab = if (is.factor(x) || 
                                                           is.character(x)) "" else as.character(substitute(y)), 
          add = FALSE, lty = 1, type = "p", ylim = NULL, lwd = 1, pch = 16, pointsize = 1.5,
          errbar.col = par("fg"), Type = rep(1, length(y)), ...) 
{
  if (is.null(ylim)) 
    ylim <- range(y[Type == 1], yplus[Type == 1], yminus[Type == 
                                                           1], na.rm = TRUE)
  if (is.factor(x) || is.character(x)) {
    x <- as.character(x)
    n <- length(x)
    t1 <- Type == 1
    t2 <- Type == 2
    n1 <- sum(t1)
    n2 <- sum(t2)
    omai <- par("mai")
    mai <- omai
    mai[2] <- max(strwidth(x, "inches")) + 0.25
    par(mai = mai)
    on.exit(par(mai = omai))
    plot(NA, NA, xlab = ylab, ylab = "", xlim = ylim, ylim = c(1, 
                                                               n + 1), axes = FALSE, ...)
    axis(1)
    w <- if (any(t2)) 
      n1 + (1:n2) + 1
    else numeric(0)
    axis(2, at = c(seq.int(length.out = n1), w), labels = c(x[t1], 
                                                            x[t2]), las = 1, adj = 1)
    points(y[t1], seq.int(length.out = n1), pch = pch, type = type, cex = pointsize)
    segments(yplus[t1], seq.int(length.out = n1), yminus[t1], 
             seq.int(length.out = n1), lwd = lwd, lty = lty, col = errbar.col)
    if (any(Type == 2)) {
      abline(h = n1 + 1, lty = 2, ...)
      offset <- mean(y[t1]) - mean(y[t2])
      if (min(yminus[t2]) < 0 & max(yplus[t2]) > 0) 
        lines(c(0, 0) + offset, c(n1 + 1, par("usr")[4]), 
              lty = 2, ...)
      points(y[t2] + offset, w, pch = pch, type = type, cex = pointsize)
      segments(yminus[t2] + offset, w, yplus[t2] + offset, 
               w, lwd = lwd, lty = lty, col = errbar.col)
      at <- pretty(range(y[t2], yplus[t2], yminus[t2]))
      axis(side = 3, at = at + offset, labels = format(round(at, 
                                                             6)))
    }
    return(invisible())
  }
  if (add) 
    points(x, y, pch = pch, type = type, cex = pointsize)
  else plot(x, y, ylim = ylim, xlab = xlab, ylab = ylab, pch = pch, 
            type = type, ...)
  xcoord <- par()$usr[1:2]
  smidge <- cap * (xcoord[2] - xcoord[1])/2
  segments(x, yminus, x, yplus, lty = lty, lwd = lwd, col = errbar.col)
  if (par()$xlog) {
    xstart <- x * 10^(-smidge)
    xend <- x * 10^(smidge)
  }
  else {
    xstart <- x - smidge
    xend <- x + smidge
  }
  segments(xstart, yminus, xend, yminus, lwd = lwd, lty = lty, 
           col = errbar.col)
  segments(xstart, yplus, xend, yplus, lwd = lwd, lty = lty, 
           col = errbar.col)
  return(invisible())
}


PlotErrBar_interactions<-function(model
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
  
  for(i in 1:length(exp1levs)){
    xtxt = (i + off + (off_increment * length(exp2levs))/2)
    text(exp1levs[i], x = xtxt, y = min(ylims), srt = srttxt)
    
  }
}
