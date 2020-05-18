

plotSpline=function(MIC,DIA,xcens,ycens,M1,M2,D1Spline,D2Spline,
                    medianDensity,medFit,lowerFit,upperFit,flipGraph,logConvert){
    
  spline_fit_data=data.frame(grid=xgrid,Median=medFit,Lower=lowerFit,Upper=upperFit)

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
  
  M1=M1+.5
  M2=M2-.5
  
  if(flipGraph=='Yes' && logConvert==TRUE){
    fit=ggplot(a1,aes(xobs,yobs))+geom_text(aes(label=freq),size=5)+
      geom_ribbon(data=spline_fit_data,aes(x=grid,ymin=Lower,ymax=Upper),
                  fill="skyblue2",alpha=.7, inherit.aes = FALSE)+
      geom_line(data=spline_fit_data,aes(x=grid,y=Median),color='blue',inherit.aes = FALSE)+
      geom_vline(xintercept=M1,lty=2,alpha=.5)+
      geom_vline(xintercept=M2,lty=2,alpha=.5)+
      geom_hline(yintercept=D1Spline+.5,lty=2,alpha=.5)+
      geom_hline(yintercept=D2Spline-.5,lty=2,alpha=.5)+
      labs(x='MIC (Dilution Test in log(ug/mL))',y='DIA (Diffusion Test in mm)')+
      scale_x_continuous(breaks = seq(min(xobs1)-1,max(xobs1)+1,by=1),
          labels = c(paste("<",min(xobs1),sep=''),seq(min(xobs1),max(xobs1),by=1), paste(">",max(xobs1),sep='')),
          limits = c(min(xobs1)-1,max(xobs1)+1))+
      scale_y_continuous(breaks = seq(min(DIA1)-1,max(DIA1)+1,by=1),
             labels = c(paste("<",min(DIA1),sep=''),seq(min(DIA1),max(DIA1),by=1), paste(">",max(DIA1),sep='')),
             limits = c(min(DIA1)-1,max(DIA1)+1))+
      theme(
           legend.position='none',
           legend.text = element_text(size = 15))
  }
  if(flipGraph=='Yes' && logConvert==FALSE){
    M1=2^M1
    M2=2^M2
    a1$xobs=2^a1$xobs
    MIC2=MIC
    MICTemp=c(min(MIC2)-1,min(MIC2):max(MIC2),max(MIC2)+1)
    MICTemp=2^MICTemp
    x=2^(min(MIC2):max(MIC2))
    spline_fit_data$grid=2^spline_fit_data$grid
    fit=ggplot(a1,aes(xobs,yobs))+geom_text(aes(label=freq),size=5)+
      geom_ribbon(data=spline_fit_data,aes(x=grid,ymin=Lower,ymax=Upper),
                  fill="skyblue2",alpha=.7, inherit.aes = FALSE)+
      geom_line(data=spline_fit_data,aes(x=grid,y=Median),color='blue',inherit.aes = FALSE)+
      geom_vline(xintercept=M1,lty=2,alpha=.5)+
      geom_vline(xintercept=M2,lty=2,alpha=.5)+
      geom_hline(yintercept=D1Spline+.5,lty=2,alpha=.5)+
      geom_hline(yintercept=D2Spline-.5,lty=2,alpha=.5)+
      labs(x='MIC (Dilution Test in ug/mL)',y='DIA (Diffusion Test in mm)')+
      scale_x_continuous(trans=log2_trans(),
           limits=c(min(MICTemp),max(MICTemp)),
           breaks=MICTemp,
           labels=c(paste("<",min(x),sep=''),sort(unique(x)), paste(">",max(x),sep='')))+
      scale_y_continuous(breaks = seq(min(DIA1)-1,max(DIA1)+1,by=1),
           labels = c(paste("<",min(DIA1),sep=''),seq(min(DIA1),max(DIA1),by=1), paste(">",max(DIA1),sep='')),
           limits = c(min(DIA1)-1,max(DIA1)+1))+
      theme(
        legend.position='none',
        legend.text = element_text(size = 15))    
    
  }
  if(flipGraph=='No' && logConvert==TRUE){
    fit=ggplot(a1,aes(yobs,xobs))+geom_text(aes(label=freq),size=5)+
      geom_ribbon(data=spline_fit_data,aes(y=grid,xmin=Lower,xmax=Upper),
                  fill="skyblue2",alpha=.7, inherit.aes = FALSE)+
      geom_line(data=spline_fit_data,aes(y=grid,x=Median),color='blue',inherit.aes = FALSE)+
      geom_hline(yintercept=M1,lty=2,alpha=.5)+
      geom_hline(yintercept=M2,lty=2,alpha=.5)+
      geom_vline(xintercept=D1Spline+.5,lty=2,alpha=.5)+
      geom_vline(xintercept=D2Spline-.5,lty=2,alpha=.5)+
      labs(y='MIC (Dilution Test in log(ug/mL))',x='DIA (Diffusion Test in mm)')+
      scale_y_continuous(breaks = seq(min(xobs1)-1,max(xobs1)+1,by=1),
          labels = c(paste("<",min(xobs1),sep=''),seq(min(xobs1),max(xobs1),by=1), paste(">",max(xobs1),sep='')),
          limits = c(min(xobs1)-1,max(xobs1)+1))+
      scale_x_continuous(breaks = seq(min(DIA1)-1,max(DIA1)+1,by=1),
           labels = c(paste("<",min(DIA1),sep=''),seq(min(DIA1),max(DIA1),by=1), paste(">",max(DIA1),sep='')),
           limits = c(min(DIA1)-1,max(DIA1)+1))+
      theme(
          legend.position='none',
          legend.text = element_text(size = 15))
  }
  if(flipGraph=='No' && logConvert==FALSE){
    M1=2^M1
    M2=2^M2
    a1$xobs=2^a1$xobs
    MIC2=MIC
    MICTemp=c(min(MIC2)-1,min(MIC2):max(MIC2),max(MIC2)+1)
    MICTemp=2^MICTemp
    x=2^(min(MIC2):max(MIC2))
    spline_fit_data$grid=2^spline_fit_data$grid
    fit=ggplot(a1,aes(yobs,xobs))+geom_text(aes(label=freq),size=5)+
      geom_ribbon(data=spline_fit_data,aes(y=grid,xmin=Lower,xmax=Upper),
                  fill="skyblue2",alpha=.7, inherit.aes = FALSE)+
      geom_line(data=spline_fit_data,aes(y=grid,x=Median),color='blue',inherit.aes = FALSE)+
      geom_hline(yintercept=M1,lty=2,alpha=.5)+
      geom_hline(yintercept=M2,lty=2,alpha=.5)+
      geom_vline(xintercept=D1Spline+.5,lty=2,alpha=.5)+
      geom_vline(xintercept=D2Spline-.5,lty=2,alpha=.5)+
      labs(y='MIC (Dilution Test in ug/mL)',x='DIA (Diffusion Test in mm)')+
      scale_y_continuous(trans=log2_trans(),
             limits=c(min(MICTemp),max(MICTemp)),
             breaks=MICTemp,
             labels=c(paste("<",min(x),sep=''),sort(unique(x)), paste(">",max(x),sep='')))+
      scale_x_continuous(breaks = seq(min(DIA1)-1,max(DIA1)+1,by=1),
             labels = c(paste("<",min(DIA1),sep=''),seq(min(DIA1),max(DIA1),by=1), paste(">",max(DIA1),sep='')),
             limits = c(min(DIA1)-1,max(DIA1)+1))+
      theme(
        legend.position='none',
        legend.text = element_text(size = 15))
  }
  
  
  return(fit)
}



plotSplineOne=function(MIC,DIA,xcens,ycens,MICBrkpt,DIABrkpt,
                    medianDensity,medFit,lowerFit,upperFit,flipGraph,logConvert){
  
  spline_fit_data=data.frame(grid=xgrid,Median=medFit,Lower=lowerFit,Upper=upperFit)
  
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
    fit=ggplot(a1,aes(xobs,yobs))+geom_text(aes(label=freq),size=5)+
      geom_ribbon(data=spline_fit_data,aes(x=grid,ymin=Lower,ymax=Upper),
                  fill="skyblue2",alpha=.7, inherit.aes = FALSE)+
      geom_line(data=spline_fit_data,aes(x=grid,y=Median),color='blue',inherit.aes = FALSE)+
      geom_vline(xintercept=MICBrkpt,lty=2,alpha=.5)+
      geom_hline(yintercept=DIABrkpt,lty=2,alpha=.5)+
      labs(x='MIC (Dilution Test in log(ug/mL))',y='DIA (Diffusion Test in mm)')+
      scale_x_continuous(breaks = seq(min(xobs1)-1,max(xobs1)+1,by=1),
                         labels = c(paste("<",min(xobs1),sep=''),seq(min(xobs1),max(xobs1),by=1), paste(">",max(xobs1),sep='')),
                         limits = c(min(xobs1)-1,max(xobs1)+1))+
      scale_y_continuous(breaks = seq(min(DIA1)-1,max(DIA1)+1,by=1),
                         labels = c(paste("<",min(DIA1),sep=''),seq(min(DIA1),max(DIA1),by=1), paste(">",max(DIA1),sep='')),
                         limits = c(min(DIA1)-1,max(DIA1)+1))+
      theme(
        legend.position='none',
        legend.text = element_text(size = 15))
  }
  if(flipGraph=='Yes' && logConvert==FALSE){
    MICBrkpt=2^MICBrkpt
    a1$xobs=2^a1$xobs
    MIC2=MIC
    MICTemp=c(min(MIC2)-1,min(MIC):max(MIC),max(MIC2)+1)
    MICTemp=2^MICTemp
    x=2^(min(MIC2):max(MIC2))
    spline_fit_data$grid=2^spline_fit_data$grid
    fit=ggplot(a1,aes(xobs,yobs))+geom_text(aes(label=freq),size=5)+
      geom_ribbon(data=spline_fit_data,aes(x=grid,ymin=Lower,ymax=Upper),
                  fill="skyblue2",alpha=.7, inherit.aes = FALSE)+
      geom_line(data=spline_fit_data,aes(x=grid,y=Median),color='blue',inherit.aes = FALSE)+
      geom_vline(xintercept=MICBrkpt,lty=2,alpha=.5)+
      geom_hline(yintercept=DIABrkpt,lty=2,alpha=.5)+
      labs(x='MIC (Dilution Test in ug/mL)',y='DIA (Diffusion Test in mm)')+
      scale_x_continuous(trans=log2_trans(),
                         limits=c(min(MICTemp),max(MICTemp)),
                         breaks=MICTemp,
                         labels=c(paste("<",min(x),sep=''),sort(unique(x)), paste(">",max(x),sep='')))+
      scale_y_continuous(breaks = seq(min(DIA1)-1,max(DIA1)+1,by=1),
                         labels = c(paste("<",min(DIA1),sep=''),seq(min(DIA1),max(DIA1),by=1), paste(">",max(DIA1),sep='')),
                         limits = c(min(DIA1)-1,max(DIA1)+1))+
      theme(
        legend.position='none',
        legend.text = element_text(size = 15))    
    
  }
  if(flipGraph=='No' && logConvert==TRUE){
    fit=ggplot(a1,aes(yobs,xobs))+geom_text(aes(label=freq),size=5)+
      geom_ribbon(data=spline_fit_data,aes(y=grid,xmin=Lower,xmax=Upper),
                  fill="skyblue2",alpha=.7, inherit.aes = FALSE)+
      geom_line(data=spline_fit_data,aes(y=grid,x=Median),color='blue',inherit.aes = FALSE)+
      geom_hline(yintercept=MICBrkpt,lty=2,alpha=.5)+
      geom_vline(xintercept=DIABrkpt,lty=2,alpha=.5)+
      labs(y='MIC (Dilution Test in log(ug/mL))',x='DIA (Diffusion Test in mm)')+
      scale_y_continuous(breaks = seq(min(xobs1)-1,max(xobs1)+1,by=1),
                         labels = c(paste("<",min(xobs1),sep=''),seq(min(xobs1),max(xobs1),by=1), paste(">",max(xobs1),sep='')),
                         limits = c(min(xobs1)-1,max(xobs1)+1))+
      scale_x_continuous(breaks = seq(min(DIA1)-1,max(DIA1)+1,by=1),
                         labels = c(paste("<",min(DIA1),sep=''),seq(min(DIA1),max(DIA1),by=1), paste(">",max(DIA1),sep='')),
                         limits = c(min(DIA1)-1,max(DIA1)+1))+
      theme(
        legend.position='none',
        legend.text = element_text(size = 15))
  }
  if(flipGraph=='No' && logConvert==FALSE){
    MICBrkpt=2^MICBrkpt
    a1$xobs=2^a1$xobs
    MIC2=MIC
    MICTemp=c(min(MIC2)-1,min(MIC):max(MIC),max(MIC2)+1)
    MICTemp=2^MICTemp
    x=2^(min(MIC2):max(MIC2))
    spline_fit_data$grid=2^spline_fit_data$grid
    fit=ggplot(a1,aes(yobs,xobs))+geom_text(aes(label=freq),size=5)+
      geom_ribbon(data=spline_fit_data,aes(y=grid,xmin=Lower,xmax=Upper),
                  fill="skyblue2",alpha=.7, inherit.aes = FALSE)+
      geom_line(data=spline_fit_data,aes(y=grid,x=Median),color='blue',inherit.aes = FALSE)+
      geom_hline(yintercept=MICBrkpt,lty=2,alpha=.5)+
      geom_vline(xintercept=DIABrkpt,lty=2,alpha=.5)+
      labs(y='MIC (Dilution Test in ug/mL)',x='DIA (Diffusion Test in mm)')+
      scale_y_continuous(trans=log2_trans(),
                         limits=c(min(MICTemp),max(MICTemp)),
                         breaks=MICTemp,
                         labels=c(paste("<",min(x),sep=''),sort(unique(x)), paste(">",max(x),sep='')))+
      scale_x_continuous(breaks = seq(min(DIA1)-1,max(DIA1)+1,by=1),
                         labels = c(paste("<",min(DIA1),sep=''),seq(min(DIA1),max(DIA1),by=1), paste(">",max(DIA1),sep='')),
                         limits = c(min(DIA1)-1,max(DIA1)+1))+
      theme(
        legend.position='none',
        legend.text = element_text(size = 15))
  }
  
  
  return(fit)
}



