
compareFitsPlotOne=function(dataSpline,dataLog,MIC,DIA,xcens,ycens,MICBrkpt,flipGraph,logConvert){
      
  # save.image("compare.Rdata")
  # load('compare.Rdata')
  
  logData=data.frame(grid=xgrid,fit=dataLog$fit,dens=dataLog$dens/sum(dataLog$dens),group="Logistic")
  splineData=data.frame(grid=xgrid,fit=dataSpline$fit,dens=dataSpline$dens/sum(dataSpline$dens),group="Spline")
  
  xobs=MIC
  yobs=DIA
  xobs1=xobs
  DIA1=yobs
  xobs[xcens==1 & xobs==max(xobs)]=max(xobs)+1
  xobs[xcens==-1 & xobs==min(xobs)]=min(xobs)-1
  yobs[ycens==1 & yobs==max(yobs)]=max(yobs)+1
  yobs[ycens==-1 & yobs==min(yobs)]=min(yobs)-1
  a1=data.frame(xobs,yobs)
  a1=count(a1,c('xobs','yobs'))
  
  if(flipGraph=='Yes' && logConvert==TRUE){
    pltRel=ggplot(data=a1,aes(x=xobs,y=yobs,label=freq))+geom_text(size=3.2,color='black')+
      geom_line(data=rbind(logData,splineData),aes(x=grid,y=fit,group=group,color=group),inherit.aes = FALSE)+
      # geom_ribbon(data=gxDat,aes(x=xgrid,ymin=lower, ymax=upper),fill='cyan4',alpha=.5,inherit.aes = FALSE)+
      geom_vline(xintercept=MICBrkpt,lty=2,alpha=.5)+
      scale_x_continuous(breaks = seq(min(xobs1)-1,max(xobs1)+1,by=1),
                         labels = c(paste("<",min(xobs1),sep=''),seq(min(xobs1),max(xobs1),by=1), paste(">",max(xobs1),sep='')),
                         limits = c(min(xobs1)-1,max(xobs1)+1))+
      scale_y_continuous(breaks = seq(min(DIA1)-1,max(DIA1)+1,by=1),
                         labels = c(paste("<",min(DIA1),sep=''),seq(min(DIA1),max(DIA1),by=1), paste(">",max(DIA1),sep='')),
                         limits = c(min(DIA1)-1,max(DIA1)+1))+
      labs(title='MIC/DIA Relationship',y='DIA (mm)',x="MIC (Dilution Test in ug/mL)",color='')+
      theme(legend.position='none')
    
    plotDens=ggplot(rbind(logData,splineData),aes(x=grid,y=dens,group=group,colour=group))+geom_line()+
      # geom_ribbon(aes(ymin=lower, ymax=upper),alpha=.5,fill='cyan4')+
      geom_vline(xintercept=MICBrkpt,lty=2,alpha=.5)+
      scale_x_continuous(breaks = seq(round(min(xgrid)),round(max(xgrid)),by=1),limits=c(min(xgrid),max(xgrid)))+
      labs(title='MIC Density',y='',x='',color='')
    
  }
  if(flipGraph=='Yes' && logConvert==FALSE){
    
    MICBrkpt=2^MICBrkpt
    a1$xobs=2^a1$xobs
    MIC2=MIC
    MICTemp=c(min(MIC2)-1,min(MIC):max(MIC),max(MIC2)+1)
    MICTemp=2^MICTemp
    x=2^(min(MIC2):max(MIC2))
    logData$grid=2^logData$grid
    splineData$grid=2^splineData$grid
        
    pltRel=ggplot(data=a1,aes(x=xobs,y=yobs,label=freq))+geom_text(size=3.2,color='black')+
      geom_line(data=rbind(logData,splineData),aes(x=grid,y=fit,group=group,color=group),inherit.aes = FALSE)+
      # geom_ribbon(data=gxDat,aes(x=xgrid,ymin=lower, ymax=upper),fill='cyan4',alpha=.5,inherit.aes = FALSE)+
      geom_vline(xintercept=MICBrkpt,lty=2,alpha=.5)+
      scale_x_continuous(trans=log2_trans(),
                         limits=c(min(MICTemp),max(MICTemp)),
                         breaks=MICTemp,
                         labels=c(paste("<",min(x),sep=''),sort(unique(x)), paste(">",max(x),sep='')))+
      scale_y_continuous(breaks = seq(min(DIA1)-1,max(DIA1)+1,by=1),
                         labels = c(paste("<",min(DIA1),sep=''),seq(min(DIA1),max(DIA1),by=1), paste(">",max(DIA1),sep='')),
                         limits = c(min(DIA1)-1,max(DIA1)+1))+
      labs(title='MIC/DIA Relationship',y='DIA (mm)',x="MIC (Dilution Test in ug/mL)",color='')+
      theme(legend.position='none')
    
    
    MICTemp=(min(MIC)-3):(max(MIC)+3)
    MICTemp=2^MICTemp
    
    plotDens=ggplot(rbind(logData,splineData),aes(x=grid,y=dens,group=group,colour=group))+geom_line()+
      # geom_ribbon(aes(ymin=lower, ymax=upper),alpha=.5,fill='cyan4')+
      geom_vline(xintercept=MICBrkpt,lty=2,alpha=.5)+
      scale_x_continuous(trans=log2_trans(),
                         limits=c(min(2^xgrid),max(2^xgrid)),
                         breaks=MICTemp)+
      labs(title='MIC Density',y='',x='',color='')
  }
  if(flipGraph=='No' &&  logConvert==TRUE){

    pltRel=ggplot(data=a1,aes(x=yobs,y=xobs,label=freq))+geom_text(size=3.2,color='black')+
      geom_line(data=rbind(logData,splineData),aes(y=grid,x=fit,group=group,color=group),inherit.aes = FALSE)+
      geom_vline(xintercept=MICBrkpt,lty=2,alpha=.5)+
      scale_y_continuous(breaks = seq(min(xobs1)-1,max(xobs1)+1,by=1),
                         labels = c(paste("<",min(xobs1),sep=''),seq(min(xobs1),max(xobs1),by=1), paste(">",max(xobs1),sep='')),
                         limits = c(min(xobs1)-1,max(xobs1)+1))+
      scale_x_continuous(breaks = seq(min(DIA1)-1,max(DIA1)+1,by=1),
                         labels = c(paste("<",min(DIA1),sep=''),seq(min(DIA1),max(DIA1),by=1), paste(">",max(DIA1),sep='')),
                         limits = c(min(DIA1)-1,max(DIA1)+1))+
      labs(title='MIC/DIA Relationship',x='DIA (mm)',y="MIC (Dilution Test in ug/mL)",color='')+
      theme(legend.position='none')
    
    
    plotDens=ggplot(rbind(logData,splineData),aes(x=grid,y=dens,group=group,colour=group))+geom_line()+
      geom_vline(xintercept=MICBrkpt,lty=2,alpha=.5)+
      scale_x_continuous(breaks = seq(round(min(xgrid)),round(max(xgrid)),by=1),limits=c(min(xgrid),max(xgrid)))+
      labs(title='MIC Density',y='',x='',color='')
  }
  if(flipGraph=='No' &&  logConvert==FALSE){
    
    MICBrkpt=2^MICBrkpt
    a1$xobs=2^a1$xobs
    MIC2=MIC
    MICTemp=c(min(MIC2)-1,min(MIC):max(MIC),max(MIC2)+1)
    MICTemp=2^MICTemp
    x=2^(min(MIC2):max(MIC2))
    logData$grid=2^logData$grid
    splineData$grid=2^splineData$grid
        
    pltRel=ggplot(data=a1,aes(x=yobs,y=xobs,label=freq))+geom_text(size=3.2,color='black')+
      geom_line(data=rbind(logData,splineData),aes(y=grid,x=fit,group=group,color=group),inherit.aes = FALSE)+
      geom_vline(xintercept=MICBrkpt,lty=2,alpha=.5)+
      scale_y_continuous(trans=log2_trans(),
                         limits=c(min(MICTemp),max(MICTemp)),
                         breaks=MICTemp,
                         labels=c(paste("<",min(x),sep=''),sort(unique(x)), paste(">",max(x),sep='')))+
      scale_x_continuous(breaks = seq(min(DIA1)-1,max(DIA1)+1,by=1),
                         labels = c(paste("<",min(DIA1),sep=''),seq(min(DIA1),max(DIA1),by=1), paste(">",max(DIA1),sep='')),
                         limits = c(min(DIA1)-1,max(DIA1)+1))+
      labs(title='MIC/DIA Relationship',x='DIA (mm)',y="MIC (Dilution Test in ug/mL)",color='')+
      theme(legend.position='none')
    
    MICTemp=(min(MIC)-3):(max(MIC)+3)
    MICTemp=2^MICTemp
    
    plotDens=ggplot(rbind(logData,splineData),aes(x=grid,y=dens,group=group,colour=group))+geom_line()+
      geom_vline(xintercept=MICBrkpt,lty=2,alpha=.5)+
      scale_x_continuous(trans=log2_trans(),
                         limits=c(min(2^xgrid),max(2^xgrid)),
                         breaks=MICTemp)+
      labs(title='MIC Density',y='',x='',color='')
  }    
  
  fit <- gridExtra::arrangeGrob(pltRel, plotDens,ncol=1, heights=c(7,5))
  
  return(fit)
}


probDIAClassOne=function(data,MICBrkpt,DIA1,DIA2){
    
  MICBrkpt=MICBrkpt-.5
  MTrue=MICBrkpt-.5
  
  # cat(file=stderr(),MICBrkpt,MTrue,DIA1,"\n")
  # write.csv(data,"data.csv") 
  
  f=data$fit
  density=data$dens
  
  e=sqrt(.5^2+.5^2)
  t=sqrt(1.5^2+1.5^2)
  
  grid=xgrid[xgrid>=min(MIC)-1 & xgrid<=max(MIC)+1]
  weights=density[xgrid>=min(MIC)-1 & xgrid<=max(MIC)+1]
  f=f[xgrid>=min(MIC)-1 & xgrid<=max(MIC)+1]
  weights=weights/sum(weights)
  
  DIABrkpt1=DIA1
  if(is.na(DIA2)==TRUE){
    DIABrkpt2=0
  }else{
    DIABrkpt2=DIA2 
  }
  
  probDIAC=rep(NA,length(grid))
  probDIAIC=rep(NA,length(grid))
  probDIAC1=rep(NA,length(grid))
  probDIAIC1=rep(NA,length(grid))
  probMIC=rep(NA,length(grid))
  
  ### DIA Probs
  for(i in 1:length(grid)){
    d=f[i]
    m=grid[i]
    
    ### MIC
    #Susceptible
    if (m<=MTrue)  probMIC[i]=pnorm(MICBrkpt,m,e)
    #resistant
    if(m>MTrue) probMIC[i]=1-pnorm(MICBrkpt,m,e)
    
    ###DIA
    #Susceptible
    if(m<=MTrue){
      probDIAC[i]=1-pnorm(DIABrkpt1-.5,d,t)
      probDIAC1[i]=1-pnorm(DIABrkpt2-.5,d,t)
    #Resistant
    }else{
      probDIAC[i]=pnorm(DIABrkpt1-.5,d,t)
      probDIAC1[i]=pnorm(DIABrkpt2-.5,d,t)
    }
    probDIAIC[i]=1-probDIAC[i]
    probDIAIC1[i]=1-probDIAC1[i]
  }

  ### Print DIA
  if(is.na(DIA2)==FALSE){
    cat('DIA Breakpoint \n')
    cat('Set 1: ',DIABrkpt1,'\n',sep='')
    cat('Set 2: ',DIABrkpt2,'\n',sep='')
	  CorrectM=sum(probMIC)/length(grid)
    Correct1=sum(probDIAC)/length(grid)
    Correct2=sum(probDIAC1)/length(grid)
    cat('\nProbability Correct Classification \n')
    temp=data.frame(CorrectM,Correct1,Correct2)
    temp[,1:3]=round(temp[,1:3],digits=3)
    temp[,1:3]=as.character(temp[,1:3])
    names(temp)=c('MIC Correct','DIA 1 Correct','DIA 2 Correct')
    name.width <- max(sapply(names(temp), nchar))
    names(temp) <- format(names(temp), width = name.width, justify = "centre")
    print(format(temp, width = name.width, justify = "centre"),row.names=FALSE,quote=FALSE)
    
	  CorrectM=sum(probMIC*weights)
    Correct1=sum(probDIAC*weights)
    Correct2=sum(probDIAC1*weights)
    cat('\nProbability Correct Classification Weighted by Isolate Distribution \n')
    temp=data.frame(CorrectM,Correct1,Correct2)
    temp[,1:3]=round(temp[,1:3],digits=4)
    temp[,1:3]=as.character(temp[,1:3])
    names(temp)=c('MIC Correct','DIA 1 Correct','DIA 2 Correct')
    name.width <- max(sapply(names(temp), nchar))
    names(temp) <- format(names(temp), width = name.width, justify = "centre")
    print(format(temp, width = name.width, justify = "centre"),row.names=FALSE,quote=FALSE)
    
    a1=data.frame(grid,probDIAC,probDIAC1,probMIC,weights)
  }else{
    cat('DIA Breakpoint: ',DIABrkpt1,'\n',sep='')
    Correct1=sum(probDIAC)/length(grid)
	  CorrectM=sum(probMIC)/length(grid)
    cat('\nProbability Correct Classification \n')
    temp=data.frame(CorrectM,Correct1)
    temp[,1:2]=round(temp[,1:2],digits=3)
    temp[,1:2]=as.character(temp[,1:2])
    names(temp)=c('MIC Correct','DIA Correct')
    name.width <- max(sapply(names(temp), nchar))
    names(temp) <- format(names(temp), width = name.width, justify = "centre")
    print(format(temp, width = name.width, justify = "centre"),row.names=FALSE,quote=FALSE)
    
    Correct1=sum(probDIAC*weights)
	  CorrectM=sum(probMIC*weights)
    cat('\nProbability Correct Classification Weighted by Isolate Distribution \n')
    temp=data.frame(CorrectM,Correct1)
    temp[,1:2]=round(temp[,1:2],digits=4)
    temp[,1:2]=as.character(temp[,1:2])
    names(temp)=c('MIC Correct','DIA Correct')
    name.width <- max(sapply(names(temp), nchar))
    names(temp) <- format(names(temp), width = name.width, justify = "centre")
    print(format(temp, width = name.width, justify = "centre"),row.names=FALSE,quote=FALSE)
    
    a1=data.frame(grid,probDIAC,probDIAC1,probMIC,weights)
  }
  
  return(a1)
  
}

plotUnderlyingDistibutionOne=function(a1,logConvert,MICBrkpt){
  
  MICBrkpt=MICBrkpt-1
  
  if(logConvert==TRUE){
    plt=ggplot(a1,aes(x=grid,y=weights))+geom_area(alpha=.3,fill='skyblue2',color='black')+
      labs(x='',y='',title='MIC Isolate Distribution')+
      geom_vline(xintercept=MICBrkpt,lty=2,alpha=1)+
      scale_x_continuous(breaks = seq(round(min(xgrid)),round(max(xgrid)),by=1),limits=c(min(xgrid),max(xgrid)))
  }else{
    
    MICBrkpt=2^MICBrkpt
    a1$grid=2^a1$grid
    MICTemp=(min(MIC)-3):(max(MIC)+3)
    MICTemp=2^MICTemp
    
    plt=ggplot(a1,aes(x=grid,y=weights))+geom_area(alpha=.3,fill='skyblue2',color='black')+
      labs(x='',y='',title='MIC Isolate Distribution')+
      geom_vline(xintercept=MICBrkpt,lty=2,alpha=1)+
      scale_x_continuous(trans=log2_trans(),
                         limits=c(min(2^xgrid),max(2^xgrid)),
                         breaks=MICTemp)
  }
  
  return(plt)
  
}

plotProbDIAClassOne=function(a1,MICBrkpt,DIA1,DIA2,logConvert){
  
  MICBrkpt=MICBrkpt-1
  xgrid=a1$grid
  probDIAC=a1$probDIAC
  if(is.na(DIA2)==FALSE)
    probDIAC1=a1$probDIAC1
  probMIC=a1$probMIC
  
  
  grid1=xgrid[xgrid<MICBrkpt-.02]
  grid2=xgrid[xgrid>MICBrkpt+.02]
  
  probMIC1=probMIC[xgrid<MICBrkpt-.02]
  probMIC2=probMIC[xgrid>MICBrkpt+.02]
  
  ProbCorrectDIA1=probDIAC[xgrid<MICBrkpt-.02]
  ProbCorrectDIA2=probDIAC[xgrid>MICBrkpt+.02]
  
  if(is.na(DIA2)==FALSE){
    ProbCorrectDIA11=probDIAC1[xgrid<MICBrkpt-.02]
    ProbCorrectDIA21=probDIAC1[xgrid>MICBrkpt+.02]
  }
  
  
  if(is.na(DIA2)==FALSE){
    a1=data.frame(grid=c(grid1,grid1,grid1),probC=c(probMIC1,ProbCorrectDIA1,ProbCorrectDIA11),
                  group=c(rep(0,length(grid1)),rep(1,length(grid1)),rep(2,length(grid1))))
    names(a1)=c('grid','Correct','Breakpoints')
    a2=data.frame(grid=c(grid2,grid2,grid2),probC=c(probMIC2,ProbCorrectDIA2,ProbCorrectDIA21),
                  group=c(rep(0,length(grid2)),rep(1,length(grid2)),rep(2,length(grid2))))
    names(a2)=c('grid','Correct','Breakpoints')
    a1=melt(a1,id=c(1,3))
    a1[a1$Breakpoints==0,2]=paste('MIC ')
    a1[a1$Breakpoints==1,2]=paste('DIA 1: ',DIA1,sep='')
    a1[a1$Breakpoints==2,2]=paste('DIA 2: ',DIA2,sep='')
    a2=melt(a2,id=c(1,3))
    a2[a2$Breakpoints==0,2]=paste('MIC ')
    a2[a2$Breakpoints==1,2]=paste('DIA 1: ',DIA1,sep='')
    a2[a2$Breakpoints==2,2]=paste('DIA 2: ',DIA2,sep='')

    a1$Breakpoints=factor(a1$Breakpoints,levels=c("MIC ",paste('DIA 1: ',DIA1,sep=''),
          paste('DIA 2: ',DIA2,sep='')))
    a2$Breakpoints=factor(a2$Breakpoints,levels=c("MIC ",paste('DIA 1: ',DIA1,sep=''),
          paste('DIA 2: ',DIA2,sep='')))

    
    if(logConvert==TRUE){
      fit=ggplot(a1,aes(x=grid,y=value,color=Breakpoints))+geom_line()+
        geom_line(data=a2,aes(x=grid,y=value,color=Breakpoints))+
        geom_vline(xintercept=MICBrkpt,lty=2,alpha=.5)+
        labs(x='MIC',y='Probability',title='Probability of Correct Classification (non-weighted)')+
        ylim(0,1)+
        scale_x_continuous(breaks = seq(round(min(xgrid)),round(max(xgrid)),by=1),limits=c(min(MIC)-1,max(MIC)+1))
    }
    if(logConvert==FALSE){
      
      MICBrkpt=2^MICBrkpt
      a1$grid=2^a1$grid
      a2$grid=2^a2$grid
      MICTemp=(min(MIC)-1):(max(MIC)+1)
      MICTemp=2^MICTemp

      fit=ggplot(a1,aes(x=grid,y=value,color=Breakpoints))+geom_line()+
        geom_line(data=a2,aes(x=grid,y=value,color=Breakpoints))+
        geom_vline(xintercept=MICBrkpt,lty=2,alpha=.5)+
        labs(x='MIC',y='Probability',title='Probability of Correct Classification (non-weighted)')+
        ylim(0,1)+
        scale_x_continuous(trans=log2_trans(),
                           limits=c(min(MICTemp),max(MICTemp)),
                           breaks=MICTemp)
    }
  }
  if(is.na(DIA2)==TRUE){
    a1=data.frame(grid=c(grid1,grid1),probC=c(probMIC1,ProbCorrectDIA1),
                  group=c(rep(0,length(grid1)),rep(1,length(grid1))))
    names(a1)=c('grid','Correct','Breakpoints')
    a2=data.frame(grid=c(grid2,grid2),probC=c(probMIC2,ProbCorrectDIA2),
                  group=c(rep(0,length(grid2)),rep(1,length(grid2))))
    names(a2)=c('grid','Correct','Breakpoints')
    a1=melt(a1,id=c(1,3))
    a1[a1$Breakpoints==1,2]=paste('DIA: ',DIA1,sep='')
    a1[a1$Breakpoints==0,2]=paste('MIC ')
    a2=melt(a2,id=c(1,3))
    a2[a2$Breakpoints==1,2]=paste('DIA: ',DIA1,sep='')
    a2[a2$Breakpoints==0,2]=paste('MIC ')

    a1$Breakpoints=factor(a1$Breakpoints,levels=c("MIC ",paste('DIA: ',DIA1,sep='')))
    a2$Breakpoints=factor(a2$Breakpoints,levels=c("MIC ",paste('DIA: ',DIA1,sep='')))
    
    if(logConvert==TRUE){
      fit=ggplot(a1,aes(x=grid,y=value,color=Breakpoints))+geom_line()+
        geom_line(data=a2,aes(x=grid,y=value,color=Breakpoints))+
        geom_vline(xintercept=MICBrkpt,lty=2,alpha=.5)+
        labs(x='MIC',y='Probability',title='Probability of Correct Classification (non-weighted)')+
        ylim(0,1)+
        scale_x_continuous(breaks = seq(round(min(xgrid)),round(max(xgrid)),by=1),limits=c(min(MIC)-1,max(MIC)+1))
    }
    if(logConvert==FALSE){
      MICBrkpt=2^MICBrkpt
      a1$grid=2^a1$grid
      a2$grid=2^a2$grid
      MICTemp=(min(MIC)-1):(max(MIC)+1)
      MICTemp=2^MICTemp
      
      fit=ggplot(a1,aes(x=grid,y=value,color=Breakpoints))+geom_line()+
        geom_line(data=a2,aes(x=grid,y=value,color=Breakpoints))+
        geom_vline(xintercept=MICBrkpt,lty=2,alpha=.5)+
        labs(x='MIC',y='Probability',title='Probability of Correct Classification (non-weighted)')+
        ylim(0,1)+
        scale_x_continuous(trans=log2_trans(),
                           limits=c(min(MICTemp),max(MICTemp)),
                           breaks=MICTemp)
    }
  }
  
  
  return(fit)
  
  
}

findDIACOne=function(yobs,gridx,weights,fit,MICBrkpt,xsig,ysig){
  
  DIABrkpt=0; index=0
  lgrid=length(gridx)
  MICBrkpt=MICBrkpt-.5
  MTrue=MICBrkpt-.5
  # cat(file=stderr(),"C Call :",MICBrkpt,MTrue,"\n")
  minDIA=min(yobs)
  maxDIA=max(yobs)
  storage.mode(gridx) <- "double"
  storage.mode(weights) <- "double"
  storage.mode(fit) <- "double"
  storage.mode(MICBrkpt) <- "double"
  storage.mode(MTrue) <- "double"
  storage.mode(xsig) <- "double"
  storage.mode(ysig) <- "double"
  storage.mode(minDIA) <- "double"
  storage.mode(maxDIA) <- "double"
  storage.mode(lgrid) <- "integer"
  storage.mode(DIABrkpt) <- "double"
  storage.mode(index) <- "double"
  temp=.C("findDIATrueOne",gridx,weights,fit,MICBrkpt,MTrue,
          xsig,ysig,minDIA,maxDIA,lgrid,DIABrkpt,index)
  DIABrkpt=temp[[11]]
  index=temp[[12]]
  #   print(c(D1,D2,index))
  
  return(list(DIABrkpt=DIABrkpt,index=index))
  
}

