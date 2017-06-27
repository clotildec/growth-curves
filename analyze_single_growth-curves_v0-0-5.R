# analysis growth curves #

# INPUT:
# data frames from code "make-data-frames", & associated table with cell cycle informations
# OUTPUT:
# the velocity analyzed on each individual track after having filtered outliers and smoothed the raw curves
# the first part of this code shows the two methods that were tried (method 2 was kept for the paper) and how to play around with the tuning parameters
# the second part of this code is the actual workflow for the analysis that was done for the paper

# import dataframes ####
df=read.csv("E:/Analysis/00_GROWTH-CURVES/00_SINGLE CURVES/00-DATA/merge_data_frames/151016_151028_rosco-ctrl.csv",stringsAsFactors = FALSE)
df=df[,c(1,6,7,9,10:12,17:21)] # we import columns labelled : c("ID" , "Frame" ,"Cell_Intensisty", "CellIdInAutoTracking", "sum_im.median","sum_mask.median","median" , "index_c", "event" ,"cond","date" ,"index"        
colnames(df)[3]="Volume1" # we rename "Cell_Intensisty" which is the name of volume columns in the MATLAB software for FXm measurement with "Volume1" 
table=read.csv("E:/Analysis/00_GROWTH-CURVES/00_SINGLE CURVES/00-DATA/merge_data_frames/151016_151028_rosco-ctrl_table.csv",stringsAsFactors = FALSE)

dirout="E:/Analysis/00_GROWTH-CURVES/00_SINGLE CURVES/170320_analysis_HeLa-rosco_151016-151028/151016_only/"
title1="151016-rosco-df2-v2-QRfiltering"



############### 1) Small prompt to check the trajectories one by one if needed ###########################################################
d<-df
# X11()
table$keep=NA
{
for (i in unique(d$index_c)){
  temp=d[which(d$index_c==i),]
  if(nrow(temp)>40){
    y=c(temp$Volume1,temp$av1,temp$av)
    x=rep(temp$tnormi,3)
    z=c(rep(1,nrow(temp)),rep(2,nrow(temp)),rep(3,nrow(temp)))
    plot.new()
    p=ggplot(data.frame(x,y,z),aes(x=x,y=y,color=factor(z)))+
      geom_line()+
      ggtitle(i)+
      scale_x_continuous(limits=c(1,24))+
      scale_y_continuous(limits=c(2000,4000))+
      theme_pub()
    print(p)
    a=readline("good ?")
    table$keep[which(table$index_c==i)]=a
  }
}
}
     


############### 2) Different methods to filter and smooth the curves and finally extract growth speed #########################
# (code for optimization and vizualization of each of the steps used for the smoothing, the clean code for the final analysis workflow is after

d<-df
#colours=c("#0000FF", "#050BF3", "#0B17E7", "#1122DC", "#172ED0", "#1C39C5", "#2245B9", "#2851AD", "#2E5CA2", "#336896", "#39738B", "#3F7F7F", "#458B73", "#4B9668", "#50A25C", "#56AD51", "#5CB945", "#62C539", "#67D02E", "#6DDC22", "#73E717", "#79F30B", "#7FFF00", "#84F901", "#8AF303", "#90ED04", "#96E706", "#9CE108", "#A1DC09", "#A7D60B", "#ADD00D", "#B3CA0E", "#B9C410", "#BFBF12", "#C4B913", "#CAB315", "#D0AD16", "#D6A718", "#DCA11A", "#E19C1B", "#E7961D", "#ED901F", "#F38A20", "#F98422", "#FF7F24")
colours=c("#BDFF00", "#D3EA00", "#E9D500", "#FFC100", "#FFA222", "#FF8444", "#FF6666", "#CA8099", "#969ACC", "#62B4FF", "#558DDA", "#4967B6","#9F0202", "#A10426", "#A4074A", "#A60A6E", "#A90D93", "#BB2970", "#CD454E", "#DF612B", "#F27E09", "#E59D1B", "#D8BD2E", "#CCDD41",
          "#BDFF00", "#D3EA00", "#E9D500", "#FFC100", "#FFA222", "#FF8444", "#FF6666", "#CA8099", "#969ACC", "#62B4FF", "#558DDA", "#4967B6","#9F0202", "#A10426", "#A4074A", "#A60A6E", "#A90D93", "#BB2970", "#CD454E", "#DF612B", "#F27E09", "#E59D1B", "#D8BD2E", "#CCDD41")


####### mtd 1 ################
# V0 = raw curve - 1) V1 = smoothed V0 - 2) V2 = V0 - outliers from V1 Filtering based on quartile range estimation on sliding windows for outlier removal 3) V3 = smoothed V2  (output = pdf file with the different steps for each curve) ####
pdf("E:/Analysis/00_GROWTH-CURVES/00_SINGLE CURVES/170424_new-fit_results/all-complete-tracks-part2.pdf",onefile = TRUE,width=2.5*3,height=1.9*2)

# parameters to tune : ####
dtvelocity=4 # will be multiplide by 2 +1
dtav2=3 # will be multiplide by 2 +1
QRFW = 0 # Quartile rangge filtering  (distance from the quartile*the IQR)
dtf = 11 # dt window quartile range filtering
dtav1=1 # multiplide by 2+1 dt of smoothing for first step

blank <- grid.rect(gp=gpar(col="white"))
raw=list()
dtplotav=list()
dtplotplotav=list()
dtrsquare=list()
dtrsquareav=list()
av5=list()
av7=list()
l=0


for(n in unique(d$cond)){
  if(n=="ctrl"){
    limy=c(1000,4000)
  } else if (n=="rosco"){
    limy=c(2000,5000)
  }
  for(nn in unique(d$date[which(d$cond==n)])){
    for (i in unique(d$index_c[which(d$cond==n & d$date==nn)])){
      l=l+1
        temp=d[which(d$index_c==i & d$cond==n & d$date==nn),]
        
        # step 1 = averaging (V1) ####
        temp$V1=NA
        dt=dtav1
        for (k in c((1+dt):(nrow(temp)-dt))){
          temp$V1[k]=mean(temp$Volume1[c((k-dt):(k+dt))],na.rm = TRUE)
        }
        # step 2 = outliers removal (V2) from (V1), the values kept are V0 ####
        # using qurtile filtering
        dt=((dtf-1)/2)
        seq=which(is.na(temp$V1)==FALSE)
        temp$V2=NA
        for (j in c((1+dt):(length(seq)-dt))){
          hist=temp$V1[c((seq[(j-dt)]):(seq[(j+dt)]))]
          if(length(hist[which(is.na(hist)==FALSE)])>2){
            hist2=which(hist<(quantile(hist,na.rm=TRUE)[[4]]+QRFW*IQR(hist,na.rm=TRUE)) & hist>(quantile(hist,na.rm=TRUE)[[2]]-QRFW*IQR(hist,na.rm=TRUE)))
            temp$V2[(c((seq[(j-dt)]):(seq[(j+dt)]))[hist2])]=temp$Volume1[(c((seq[(j-dt)]):(seq[(j+dt)]))[hist2])]
          } else {}
        }

        # step 3 = averaging of V0 - V2 ####
        
        # average=7
        temp=cbind(temp,"av7"=rep(NA,nrow(temp)))
        dt=dtav2
        for (k in c((1+dt):(nrow(temp)-dt))){
          temp[k,ncol(temp)]=mean(temp$V2[c((k-dt):(k+dt))],na.rm = TRUE)
        }
          dt=dtvelocity
          k=dt
          m=data.frame(a=rep(NA,nrow(temp)),b=rep(NA,nrow(temp)),c=rep(NA,nrow(temp)),d=rep(NA,nrow(temp)),dvdt=rep(NA,nrow(temp)),R2=rep(NA,nrow(temp)),dvdt2=rep(NA,nrow(temp)))
          colnames(m)=c(paste0("x1-av-",k),paste0("y1-av-",k),paste0("x2-av-",k),paste0("y2-av-",k),paste0("advdt-av-",k),paste0("R2-av-",k),paste0("bdvdt-av-",k))
          temp<-cbind(temp,m)
          for (j in c((1+dt):(nrow(temp)-dt))){
            x=temp$tnormi[c((j-dt):(j+dt))]
            y=temp$av7[c((j-dt):(j+dt))]
            temp[j,(ncol(temp)-2)]=lmrob(y~x)$coeff[[2]]
            temp[j,(ncol(temp))]=lmrob(y~x)$coeff[[1]]
            # temp[j,(ncol(temp)-1)]=summary(lmrob(y~x))$r.squared
            temp[j,(ncol(temp)-1)]=summary(lmrob(y~x))$sigma
            # temp[j,(ncol(temp)-6)]=temp$tnormi[(j-dt)] #x1
            temp[j,(ncol(temp)-6)]=temp$tnormi[j]-0.5 #x1
            temp[j,(ncol(temp)-5)]=temp[j,(ncol(temp)-2)]*temp[j,(ncol(temp)-6)]+temp[j,(ncol(temp))]
            # temp[j,(ncol(temp)-4)]=temp$tnormi[(j+dt)] #x1
            temp[j,(ncol(temp)-4)]=temp$tnormi[j]+0.5
            temp[j,(ncol(temp)-3)]=temp[j,(ncol(temp)-2)]*temp[j,(ncol(temp)-4)]+temp[j,(ncol(temp))]
          }

        # dtplotav[[l]]=list()
        # dtplotplotav[[l]]=list()
        # dtrsquareav[[l]]=list()

        p=dtvelocity

          col=c(which(colnames(temp)=="Volume1"),which(colnames(temp)=="V2"),which(colnames(temp)=="av7"),which(colnames(temp)=="tnormi"),which(colnames(temp)==paste0("x1-av-",p)):(which(colnames(temp)==paste0("x1-av-",p))+5))
          temptemp2=temp[,col]
          colnames(temptemp2)=c("Volume1","av1","av7","tnormi","x","y","xend","yend","dvdt","r2")
          dtplotav[[l]]=ggplot(temptemp2,aes(x=tnormi,y=dvdt))+
            geom_line(size=0.3)+
            # geom_segment(aes(x=x,y=y,xend=xend,yend=yend),size=0.2)+
            theme_pub()+
            labs(x="v (smoothed w=7)",y="dvdt",title=paste0(n,"-cell ",i,"-window fit:",(p*2+1)))+
            theme(legend.position="none")+
            scale_y_continuous(limits=c(-100,300))#+
          # scale_x_continuous(limits=c(1000,3500))
          dtplotplotav[[l]]=ggplot(temptemp2,aes(x=tnormi,y=av1))+
            geom_point(size=0.1,aes(x=tnormi,y=Volume1),colour=table_aes_ggplot$color[3])+
            geom_point(size=0.1)+
            geom_segment(aes(x=x,y=y,xend=xend,yend=yend,colour=factor(tnormi)),size=0.2)+
            theme_pub()+
            labs(x="time",y="volume (smoothed w=7)",title=paste0(n,"-cell ",i,"-window fit:",(p*2+1)))+
            theme(legend.position="none")#+
          # scale_x_continuous(limits=c(0,20))+
          # scale_y_continuous(limits=c(1000,3500))
          dtrsquareav[[l]]=ggplot(temptemp2,aes(x=tnormi,y=r2))+
            geom_line(size=0.3)+
            theme_pub()+
            labs(x="tnormi",y="residuals std",title=paste0(n,"-cell ",i,"-window fit:",(p*2+1)))+
            theme(legend.position="none")+
            scale_y_continuous(limits=c(0,50))
        
        
        # raw[[l]]=ggplot(data=data.frame(x=rep(temp$tnormi,2),y=c(temp$Volume1,temp$av1),z=c(rep("raw",nrow(temp)),rep("iqr",nrow(temp)))),
        #                 aes(x=x,y=y,colour=factor(z)))+
        #   geom_line(size=0.3)+
        #   scale_colour_manual(values=table_aes_ggplot$color[c(3,2)])+
        #   labs(x="time",y="volume",title=i)+
        #   theme_pub()+
        #   scale_y_continuous(limits=c(800,3500))+
        #   scale_x_continuous(limits=c(0,20))#+
        # # theme(legend.position = "bottom")
        raw[[l]]=ggplot(temp,aes(x=tnormi,y=Volume1))+
          # geom_point(aes(x=tnormi,y=av1),size=0.3,colour=table_aes_ggplot$color[2])+
          geom_line(size=0.3,colour=table_aes_ggplot$color[1])+
          # scale_colour_manual(values=table_aes_ggplot$color[c(3,2)])+
          labs(x="time",y="volume",title="raw + QR filtering")+
          theme_pub()+
          scale_y_continuous(limits=limy)+
          scale_x_continuous(limits=c(0,20))#+
        # theme(legend.position = "bottom")
        # 
        av7[[l]]=ggplot(temptemp2,aes(x=tnormi,y=av7))+
          geom_point(aes(x=tnormi,y=Volume1),size=0.3,colour=table_aes_ggplot$color[3])+
          geom_point(aes(x=tnormi,y=av1),size=0.3,colour=table_aes_ggplot$color[1])+
          geom_line(size=0.3)+
          scale_colour_manual(values=table_aes_ggplot$color[c(1,2)])+
          labs(x="time",y="volume",title="RAW + QRF on smoothed +AV dt=7")+
          theme_pub()+
          scale_y_continuous(limits=limy)+
          scale_x_continuous(limits=c(0,20))#+
        # theme(legend.position = "bottom")
       
        
        grid.arrange(arrangeGrob(grobs=list(raw[[l]],av7[[l]],blank),ncol=3),
                     # arrangeGrob(grobs=dtplot[[l]],ncol=1),
                     # arrangeGrob(grobs=dtplotplot[[l]],ncol=1),
                     # arrangeGrob(grobs=dtrsquare[[l]],ncol=1),
                     # arrangeGrob(grobs=dtplotav[[l]],ncol=1),
                     # arrangeGrob(grobs=dtplotplotav[[l]],ncol=1),
                     # arrangeGrob(grobs=dtrsquareav[[l]],ncol=1),
                     arrangeGrob(grobs=list(dtplotav[[l]],dtplotplotav[[l]],dtrsquareav[[l]]),ncol=3),
                     top=paste0("1) smoothing, dt=",dtav1*2+1,"; 2) QRFiltering:",QRFW,"*IQR range, dt=",dtf, 
                                "\n3) smoothing on (raw-2): dt=",dtav2*2+1,"; 4) robust fit, dt=",dtvelocity*2+1),
                     nrow=2)
        # } else if(l==35){
        #   dev.off()
      # } else {}
      }
      
      
      # ddf=rbind(ddf,temp=
    
  }
}
dev.off()





####### method 2 ##################################
# V0 = raw curve - 1) V1 = averaged V0 - 2) V2 = V0 - outliers from V1 estimated by measuring the distance from the mean for each point # 3) V3 = averaged V2
# (also described in SI of the paper)
pdf("v16_rel-CV_deviation-from-local-hist-median_v3.pdf",onefile = TRUE,width=2.3*2,height=1.84*3)

# parameters to tune : ####
dtvelocity=4 # will be multiplide by 2 +1
dtav2=3 # will be multiplide by 2 +1
SD = rev(c(1,1.5,2)) # Quartile rangge filtering  (distance from the quartile*the IQR)
dtf = 11 # dt window quartile range filtering
dtav1=c(2,4,5) # multiplide by 2+1 dt of smoothing for first step
title="abs((xi/local median -1))<n"

# blank <- grid.rect(gp=gpar(col="white"))
l=0
# relative values for SD, but delta normalized by local median
for(n in unique(d$cond)){
  if(n=="ctrl"){
    limy=c(1000,4000)
  } else if (n=="rosco"){
    limy=c(2000,5000)
  }
  for(nn in unique(d$date[which(d$cond==n)])){
    for (i in unique(d$index_c[which(d$cond==n & d$date==nn & d$index_c<1000)])){
      l=l+1
      temp=d[which(d$index_c==i & d$cond==n & d$date==nn),]
      temp$tnormi=c(0:(nrow(temp)-1))*10/60
      
      V1=list()
      delta=list()
      v2=list()
      V1plot=list()
      deltaplot=list()
      V2plot=list()
      
      for (p in c(1:length(dtav1))){
        V1[[p]]=rep(NA,nrow(temp))
        delta[[p]]=rep(NA,nrow(temp))
        #step 1 = averaging (V1) ####
        dt=dtav1[p]
        for (k in c((1+dt):(nrow(temp)-dt))){
          V1[[p]][k]=mean(temp$Volume1[c((k-dt):(k+dt))],na.rm = TRUE)
        }
        
        # step 2 = calculate distance to average ####
        dt=2
        for (k in c((1+dt):(nrow(temp)-dt))){
          delta[[p]][k]=(temp$Volume1[k]-V1[[p]][k])/median(V1[[p]][c((k-dt):(k+dt))],na.rm=TRUE)
        }
        seq=c(c(1:dt),c((nrow(temp)-dt+1):nrow(temp)))
        for (k in seq[1:(length(seq)/2)]){
          delta[[p]][k]=(temp$Volume1[k]-V1[[p]][k])/median(V1[[p]][c(k:(k+dt*2))],na.rm=TRUE)
        }
        for (k in seq[(length(seq)/2+1):length(seq)]){
          delta[[p]][k]=(temp$Volume1[k]-V1[[p]][k])/median(V1[[p]][c((k-dt*2):k)],na.rm=TRUE)
        }
        
      
        # step 2 part 2 = outliers removal (V2) from (V1), the values kept are V0 ####
        V2[[p]]=list()
        V2res=NULL
        V2lab=NULL
        mu=mean(delta[[p]],na.rm = TRUE)
        sd=sd(delta[[p]],na.rm = TRUE)
        for (pp in c(1:length(SD))){
          dt=((dtf-1)/2)
          filter=SD[pp]
          
          V2[[p]][[pp]]=rep(NA,nrow(temp))
          val=which(abs(delta[[p]])<mu+filter*sd)
          V2[[p]][[pp]][val]=temp$Volume1[val]
          V2res=c(V2res,V2[[p]][[pp]])
          V2lab=c(V2lab,rep(pp,length(V2[[p]][[pp]])))
        }

        if(l<35){ 
        V1plot[[p]]=ggplot(data.frame(x=temp$tnormi,y=V1[[p]]),aes(x=x,y=y))+
          geom_line(size=0.3,colour=table_aes[[2]]$color[p])+
          labs(x="time (hrs)",y=paste0("volume, dt smooth=",dtav1[p]*2+1))+
          scale_y_continuous(limits=limy)+
          scale_x_continuous(limits=c(0,24))+
          theme_pub()
        deltaplot[[p]]=ggplot(data.frame(x=temp$tnormi,y=delta[[p]]),aes(x=y))+
          geom_histogram(aes(y=..count../sum(..count..)))+
          labs(y=paste0("delta, dt smooth=",dtav1[p]*2+1))+
          # scale_y_continuous(limits=limy)+
          scale_x_continuous(limits=c(-0.06,0.06))+
          scale_y_continuous(limits=c(0,0.3))+
          theme_pub()
        for (pp in c(1:length(SD))){
          deltaplot[[p]]=deltaplot[[p]]+
            geom_vline(xintercept=c((mu-SD[pp]*sd),(mu+SD[pp]*sd)),colour=table_aes[[1]]$color[(pp+1)],size=0.3)
          }
        V2plot[[p]]=ggplot(data.frame(x=rep(temp$tnormi,(length(SD)+1)),y=c(temp$Volume1,V2res),z=c(rep(0,nrow(temp)),V2lab)),
                           aes(x=x,y=y,colour=factor(z)))+
          geom_point(size=0.4)+
          scale_colour_manual(values=table_aes[[1]]$color[c(1:(length(SD)+1))],labels=c("raw",as.character(SD)),name="mean+SD*:")+
          theme_pub()+
          scale_y_continuous(limits=limy)+
          scale_x_continuous(limits=c(0,24))+
          labs(x="time(hrs)",y="Volume")+
          theme(legend.position="none")
        } else {}
      }
      
      if(l<35){ 
      raw=ggplot(temp,aes(x=tnormi,y=Volume1))+
        geom_point(size=0.5,colour=table_aes[[1]]$color[1])+
        theme_pub()+
        scale_y_continuous(limits=limy)+
        scale_x_continuous(limits=c(0,24))+
        labs(x="time",y="Volume RAW")
      
      
 
      grid.arrange(arrangeGrob(grobs=list(blank,raw,blank),ncol=3),
                   arrangeGrob(grobs=list(V1plot[[1]],deltaplot[[1]],V2plot[[1]]),ncol=3),
                   arrangeGrob(grobs=list(V1plot[[2]],deltaplot[[2]],V2plot[[2]]),ncol=3),
                   arrangeGrob(grobs=list(V1plot[[3]],deltaplot[[3]],V2plot[[3]]),ncol=3),
                   nrow=4,
                   top=title
                   # top=paste0("1) smoothing, dt=",dtav1*2+1,"; 2) QRFiltering:",QRFW,"*IQR range, dt=",dtf, 
                   #            "\n3) smoothing on (raw-2): dt=",dtav2*2+1,"; 4) robust fit, dt=",dtvelocity*2+1),
                  )
      } else {
        dev.off()
      }
    }
  }
}

# abs values for SD filtering
SD=rev(c(0.025,0.035,0.45))
l=0
for(n in unique(d$cond)){
  if(n=="ctrl"){
    limy=c(1000,4000)
  } else if (n=="rosco"){
    limy=c(2000,5000)
  }
  for(nn in unique(d$date[which(d$cond==n)])){
    for (i in unique(d$index_c[which(d$cond==n & d$date==nn & d$index_c<1000)])){
      l=l+1
      temp=d[which(d$index_c==i & d$cond==n & d$date==nn),]
      temp$tnormi=c(0:(nrow(temp)-1))*10/60
      
      V1=list()
      delta=list()
      v2=list()
      V1plot=list()
      deltaplot=list()
      V2plot=list()
      
      for (p in c(1:length(dtav1))){
        V1[[p]]=rep(NA,nrow(temp))
        delta[[p]]=rep(NA,nrow(temp))
        #step 1 = averaging (V1) ####
        dt=dtav1[p]
        for (k in c((1+dt):(nrow(temp)-dt))){
          V1[[p]][k]=median(temp$Volume1[c((k-dt):(k+dt))],na.rm = TRUE)
        }
        
        # step 2 = calculate distance to average ####
        # dt=2
        # for (k in c((1+dt):(nrow(temp)-dt))){
        #   delta[[p]][k]=(temp$Volume1[k]-V1[[p]][k])/median(V1[[p]][c((k-dt):(k+dt))],na.rm=TRUE)
        # }
        # seq=c(c(1:dt),c((nrow(temp)-dt+1):nrow(temp)))
        # for (k in seq[1:(length(seq)/2)]){
        #   delta[[p]][k]=(temp$Volume1[k]-V1[[p]][k])/median(V1[[p]][c(k:(k+dt*2))],na.rm=TRUE)
        # }
        # for (k in seq[(length(seq)/2+1):length(seq)]){
        #   delta[[p]][k]=(temp$Volume1[k]-V1[[p]][k])/median(V1[[p]][c((k-dt*2):k)],na.rm=TRUE)
        # }
        delat[[p]]=(temp$Volume1-V1[[p]])/mean(V1[[p]],na.rm=TRUE)
        
        
        # step 2 part 2 = outliers removal (V2) from (V1), the values kept are V0 ####
        V2[[p]]=list()
        V2res=NULL
        V2lab=NULL
        # mu=mean(delta[[p]],na.rm = TRUE)
        # sd=sd(delta[[p]],na.rm = TRUE)
        mu=0
        for (pp in c(1:length(SD))){
          dt=((dtf-1)/2)
          filter=SD[pp]
          
          V2[[p]][[pp]]=rep(NA,nrow(temp))
          val=which(abs(delta[[p]])<mu+filter)
          V2[[p]][[pp]][val]=temp$Volume1[val]
          V2res=c(V2res,V2[[p]][[pp]])
          V2lab=c(V2lab,rep(pp,length(V2[[p]][[pp]])))
        }
        
        if(l<35){ 
          V1plot[[p]]=ggplot(data.frame(x=temp$tnormi,y=V1[[p]]),aes(x=x,y=y))+
            geom_line(size=0.3,colour=table_aes[[2]]$color[p])+
            labs(x="time (hrs)",y=paste0("volume, dt smooth=",dtav1[p]*2+1))+
            scale_y_continuous(limits=limy)+
            scale_x_continuous(limits=c(0,24))+
            theme_pub()
          deltaplot[[p]]=ggplot(data.frame(x=temp$tnormi,y=delta[[p]]),aes(x=y))+
            geom_histogram(aes(y=..count../sum(..count..)))+
            labs(y=paste0("delta, dt smooth=",dtav1[p]*2+1))+
            # scale_y_continuous(limits=limy)+
            scale_x_continuous(limits=c(-0.06,0.06))+
            scale_y_continuous(limits=c(0,0.3))+
            theme_pub()
          for (pp in c(1:length(SD))){
            deltaplot[[p]]=deltaplot[[p]]+
              geom_vline(xintercept=c((mu-SD[pp]),(mu+SD[pp])),colour=table_aes[[1]]$color[(pp+1)],size=0.3)
          }
          V2plot[[p]]=ggplot(data.frame(x=rep(temp$tnormi,(length(SD)+1)),y=c(temp$Volume1,V2res),z=c(rep(0,nrow(temp)),V2lab)),
                             aes(x=x,y=y,colour=factor(z)))+
            geom_point(size=0.4)+
            scale_colour_manual(values=table_aes[[1]]$color[c(1:(length(SD)+1))],labels=c("raw",as.character(SD)),name="mean+SD*:")+
            theme_pub()+
            scale_y_continuous(limits=limy)+
            scale_x_continuous(limits=c(0,24))+
            labs(x="time(hrs)",y="Volume")+
            theme(legend.position="none")
        } else {}
      }
      
      if(l<35){ 
        raw=ggplot(temp,aes(x=tnormi,y=Volume1))+
          geom_point(size=0.5,colour=table_aes[[1]]$color[1])+
          theme_pub()+
          scale_y_continuous(limits=limy)+
          scale_x_continuous(limits=c(0,24))+
          labs(x="time",y="Volume RAW")
        
        
        
        grid.arrange(arrangeGrob(grobs=list(blank,raw,blank),ncol=3),
                     arrangeGrob(grobs=list(V1plot[[1]],deltaplot[[1]],V2plot[[1]]),ncol=3),
                     arrangeGrob(grobs=list(V1plot[[2]],deltaplot[[2]],V2plot[[2]]),ncol=3),
                     arrangeGrob(grobs=list(V1plot[[3]],deltaplot[[3]],V2plot[[3]]),ncol=3),
                     nrow=4,
                     top=title
                     # top=paste0("1) smoothing, dt=",dtav1*2+1,"; 2) QRFiltering:",QRFW,"*IQR range, dt=",dtf, 
                     #            "\n3) smoothing on (raw-2): dt=",dtav2*2+1,"; 4) robust fit, dt=",dtvelocity*2+1),
        )
      } else {
        dev.off()
      }
    }
  }
}

# position of value from local histogramm
SD=rev(c(0.035,0.045,0.055))
l=0
for(n in unique(d$cond)){
  if(n=="ctrl"){
    limy=c(1000,4000)
  } else if (n=="rosco"){
    limy=c(2000,5000)
  }
  for(nn in unique(d$date[which(d$cond==n)])){
    for (i in unique(d$index_c[which(d$cond==n & d$date==nn & d$index_c<1000)])){
      l=l+1
      temp=d[which(d$index_c==i & d$cond==n & d$date==nn),]
      temp$tnormi=c(0:(nrow(temp)-1))*10/60
      
      V1=list()
      SD1=list()
      delta=list()
      v2=list()
      V1plot=list()
      deltaplot=list()
      V2plot=list()
      
      for (p in c(1:length(dtav1))){
        V1[[p]]=rep(NA,nrow(temp))
        delta[[p]]=rep(NA,nrow(temp))
        SD1[[p]]=rep(NA,nrow(temp))
        #step 1 = averaging (V1) ####
        dt=dtav1[p]
        for (k in c((1+dt):(nrow(temp)-dt))){
          V1[[p]][k]=median(temp$Volume1[c((k-dt):(k+dt))],na.rm = TRUE)
          SD1[[p]][k]=sd(temp$Volume1[c((k-dt):(k+dt))],na.rm = TRUE)
        }
        
        # step 2 = calculate distance to average ####
        # dt=2
        # for (k in c((1+dt):(nrow(temp)-dt))){
        #   delta[[p]][k]=(temp$Volume1[k]-V1[[p]][k])/median(V1[[p]][c((k-dt):(k+dt))],na.rm=TRUE)
        # }
        # seq=c(c(1:dt),c((nrow(temp)-dt+1):nrow(temp)))
        # for (k in seq[1:(length(seq)/2)]){
        #   delta[[p]][k]=(temp$Volume1[k]-V1[[p]][k])/median(V1[[p]][c(k:(k+dt*2))],na.rm=TRUE)
        # }
        # for (k in seq[(length(seq)/2+1):length(seq)]){
        #   delta[[p]][k]=(temp$Volume1[k]-V1[[p]][k])/median(V1[[p]][c((k-dt*2):k)],na.rm=TRUE)
        # }
        # delta[[p]]=(temp$Volume1-V1[[p]])/mean(V1[[p]],na.rm=TRUE)
        
        
        # step 2 part 2 = outliers removal (V2) from (V1 +/- n*SD), the values kept are V0 ####
        V2[[p]]=list()
        V2res=NULL
        V2lab=NULL
        # abszolute:
        # mu=mean(delta[[p]],na.rm = TRUE)
        # sd=sd(delta[[p]],na.rm = TRUE)
        
        # relative:
        mu=0
        for (pp in c(1:length(SD))){
          dt=((dtf-1)/2)
          filter=SD[pp]
          
          V2[[p]][[pp]]=rep(NA,nrow(temp))
          val=which(abs(temp$Volume1/V1[[p]]-1)<mu+filter)
          V2[[p]][[pp]][val]=temp$Volume1[val]
          V2res=c(V2res,V2[[p]][[pp]])
          V2lab=c(V2lab,rep(pp,length(V2[[p]][[pp]])))
        }
        
        # if(l<35){ 
          # V1plot[[p]]=ggplot(data.frame(x=temp$tnormi,y=V1[[p]]),aes(x=x,y=y))+
          #   geom_line(size=0.3,colour=table_aes[[2]]$color[p])+
          #   labs(x="time (hrs)",y=paste0("volume, dt smooth=",dtav1[p]*2+1))+
          #   scale_y_continuous(limits=limy)+
          #   scale_x_continuous(limits=c(0,24))+
          #   theme_pub()
          deltaplot[[p]]=ggplot(data.frame(x=temp$tnormi,y=(temp$Volume1/V1[[p]]-1)),aes(x=y))+
            geom_histogram(aes(y=..count../sum(..count..)))+
            labs(y=paste0("(xi/local mean -1), dt hist=",dtav1[p]*2+1))+
            # scale_y_continuous(limits=limy)+
            scale_x_continuous(limits=c(-0.06,0.06))+
            scale_y_continuous(limits=c(0,0.3))+
            theme_pub()
          for (pp in c(1:length(SD))){
            deltaplot[[p]]=deltaplot[[p]]+
              geom_vline(xintercept=c((mu-SD[pp]),(mu+SD[pp])),colour=table_aes[[1]]$color[(pp+1)],size=0.3)
          }
          V2plot[[p]]=ggplot(data.frame(x=rep(temp$tnormi,(length(SD)+1)),y=c(temp$Volume1,V2res),z=c(rep(0,nrow(temp)),V2lab)),
                             aes(x=x,y=y,colour=factor(z)))+
            geom_point(size=0.4)+
            scale_colour_manual(values=table_aes[[1]]$color[c(1:(length(SD)+1))],labels=c("raw",as.character(SD)),name="mean+SD*:")+
            theme_pub()+
            scale_y_continuous(limits=limy)+
            scale_x_continuous(limits=c(0,24))+
            labs(x="time(hrs)",y="Volume")+
            theme(legend.position="none")
        # } else {}
      }
      
      # if(l<35){ 
        raw=ggplot(temp,aes(x=tnormi,y=Volume1))+
          geom_point(size=0.5,colour=table_aes[[1]]$color[1])+
          theme_pub()+
          scale_y_continuous(limits=limy)+
          scale_x_continuous(limits=c(0,24))+
          labs(x="time",y="Volume RAW")
        
        
        
        grid.arrange(
          # arrangeGrob(grobs=list(blank,raw),ncol=2),
                     arrangeGrob(grobs=list(deltaplot[[1]],V2plot[[1]]),ncol=2),
                     arrangeGrob(grobs=list(deltaplot[[2]],V2plot[[2]]),ncol=2),
                     arrangeGrob(grobs=list(deltaplot[[3]],V2plot[[3]]),ncol=2),
                     nrow=4,
                     top=title
                     # top=paste0("1) smoothing, dt=",dtav1*2+1,"; 2) QRFiltering:",QRFW,"*IQR range, dt=",dtf, 
                     #            "\n3) smoothing on (raw-2): dt=",dtav2*2+1,"; 4) robust fit, dt=",dtvelocity*2+1),
        )
      # } else {
        # dev.off()
      # }
    }
  }
}

# position of value from local histogramm, n* local CV
title="abs((xi/local median -1))<n*local cv"

SD=rev(c(0.0004,0.0005,0.0008))
title=paste0("abs((xi/local median -1))<n*local cv\n",as.character(SD))
l=0
for(n in unique(d$cond)){
  if(n=="ctrl"){
    limy=c(1000,4000)
  } else if (n=="rosco"){
    limy=c(2000,5000)
  }
  for(nn in unique(d$date[which(d$cond==n)])){
    for (i in unique(d$index_c[which(d$cond==n & d$date==nn & d$index_c<1000)])){
      l=l+1
      temp=d[which(d$index_c==i & d$cond==n & d$date==nn),]
      temp$tnormi=c(0:(nrow(temp)-1))*10/60
      
      V1=list()
      SD1=list()
      delta=list()
      v2=list()
      V1plot=list()
      deltaplot=list()
      V2plot=list()
      
      for (p in c(1:length(dtav1))){
        V1[[p]]=rep(NA,nrow(temp))
        delta[[p]]=rep(NA,nrow(temp))
        SD1[[p]]=rep(NA,nrow(temp))
        #step 1 = averaging (V1) ####
        dt=dtav1[p]
        for (k in c((1+dt):(nrow(temp)-dt))){
          V1[[p]][k]=median(temp$Volume1[c((k-dt):(k+dt))],na.rm = TRUE)
          SD1[[p]][k]=sd(temp$Volume1[c((k-dt):(k+dt))],na.rm = TRUE)
        }
        
        # step 2 = calculate distance to average ####
        # dt=2
        # for (k in c((1+dt):(nrow(temp)-dt))){
        #   delta[[p]][k]=(temp$Volume1[k]-V1[[p]][k])/median(V1[[p]][c((k-dt):(k+dt))],na.rm=TRUE)
        # }
        # seq=c(c(1:dt),c((nrow(temp)-dt+1):nrow(temp)))
        # for (k in seq[1:(length(seq)/2)]){
        #   delta[[p]][k]=(temp$Volume1[k]-V1[[p]][k])/median(V1[[p]][c(k:(k+dt*2))],na.rm=TRUE)
        # }
        # for (k in seq[(length(seq)/2+1):length(seq)]){
        #   delta[[p]][k]=(temp$Volume1[k]-V1[[p]][k])/median(V1[[p]][c((k-dt*2):k)],na.rm=TRUE)
        # }
        # delta[[p]]=(temp$Volume1-V1[[p]])/mean(V1[[p]],na.rm=TRUE)
        
        
        # step 2 part 2 = outliers removal (V2) from (V1 +/- n*SD), the values kept are V0 ####
        V2[[p]]=list()
        V2res=NULL
        V2lab=NULL
        # abszolute:
        # mu=mean(delta[[p]],na.rm = TRUE)
        # sd=sd(delta[[p]],na.rm = TRUE)
        
        # relative:
        mu=0
        for (pp in c(1:length(SD))){
          dt=((dtf-1)/2)
          filter=SD[pp]
          
          V2[[p]][[pp]]=rep(NA,nrow(temp))
          val=which(abs(temp$Volume1/V1[[p]]-1)<V1[[p]]/SD1[[p]]*filter)
          V2[[p]][[pp]][val]=temp$Volume1[val]
          V2res=c(V2res,V2[[p]][[pp]])
          V2lab=c(V2lab,rep(pp,length(V2[[p]][[pp]])))
        }
        
        # if(l<35){ 
        # V1plot[[p]]=ggplot(data.frame(x=temp$tnormi,y=V1[[p]]),aes(x=x,y=y))+
        #   geom_line(size=0.3,colour=table_aes[[2]]$color[p])+
        #   labs(x="time (hrs)",y=paste0("volume, dt smooth=",dtav1[p]*2+1))+
        #   scale_y_continuous(limits=limy)+
        #   scale_x_continuous(limits=c(0,24))+
        #   theme_pub()
        deltaplot[[p]]=ggplot(data.frame(x=temp$tnormi,y=(temp$Volume1/V1[[p]]-1)),aes(x=y))+
          geom_histogram(aes(y=..count../sum(..count..)))+
          labs(y=paste0("dt hist=",dtav1[p]*2+1))+
          # scale_y_continuous(limits=limy)+
          scale_x_continuous(limits=c(-0.1,0.1))+
          scale_y_continuous(limits=c(0,0.3))+
          theme_pub()
        # for (pp in c(1:length(SD))){
        #   deltaplot[[p]]=deltaplot[[p]]+
        #     geom_vline(xintercept=c((mu-SD[pp]),(mu+SD[pp])),colour=table_aes[[1]]$color[(pp+1)],size=0.3)
        # }
        V2plot[[p]]=ggplot(data.frame(x=rep(temp$tnormi,(length(SD)+1)),y=c(temp$Volume1,V2res),z=c(rep(0,nrow(temp)),V2lab)),
                           aes(x=x,y=y,colour=factor(z)))+
          geom_point(size=0.4)+
          scale_colour_manual(values=table_aes[[1]]$color[c(1:(length(SD)+1))],labels=c("raw",as.character(SD)),name="mean+SD*:")+
          theme_pub()+
          scale_y_continuous(limits=limy)+
          scale_x_continuous(limits=c(0,24))+
          labs(x="time(hrs)",y="Volume")+
          theme(legend.position="none")
        # } else {}
      }
      
      # if(l<35){ 
      raw=ggplot(temp,aes(x=tnormi,y=Volume1))+
        geom_point(size=0.5,colour=table_aes[[1]]$color[1])+
        theme_pub()+
        scale_y_continuous(limits=limy)+
        scale_x_continuous(limits=c(0,24))+
        labs(x="time",y="Volume RAW")
      
      
      
      grid.arrange(
        # arrangeGrob(grobs=list(blank,raw),ncol=2),
        arrangeGrob(grobs=list(deltaplot[[1]],V2plot[[1]]),ncol=2),
        arrangeGrob(grobs=list(deltaplot[[2]],V2plot[[2]]),ncol=2),
        arrangeGrob(grobs=list(deltaplot[[3]],V2plot[[3]]),ncol=2),
        nrow=3,
        top=title
        # top=paste0("1) smoothing, dt=",dtav1*2+1,"; 2) QRFiltering:",QRFW,"*IQR range, dt=",dtf, 
        #            "\n3) smoothing on (raw-2): dt=",dtav2*2+1,"; 4) robust fit, dt=",dtvelocity*2+1),
      )
      # } else {
      # dev.off()
      # }
    }
  }
}
dev.off()

# vizualize better
dtav1=5
SD=rev(c(0.004,0.005,0.008))
SD=0.005
# pdf("v16_rel-CV_deviation-from-local-hist-median_visualize-alluncomplete.pdf",onefile = TRUE,width=2.3*1,height=2.2*1)

# check one by one the curves from method 2 ####
df=read.csv("E:/Analysis/00_GROWTH-CURVES/00_SINGLE CURVES/170515/1-make-data-frame/mtd2_dtav1-5_SD-filtering0p004_dtav3_dtvelocity4_complete-uncomplete.csv",stringsAsFactors = FALSE)
df=df[,c(1,6,7,9,10:12,17:21)]
colnames(df)[3]="Volume1"
d<-df
table=read.csv("E:/Analysis/00_GROWTH-CURVES/00_SINGLE CURVES/00-DATA/merge_data_frames_v2/151016_151028_rosco-ctrl_table.csv")
df=read.csv("E:/Analysis/00_GROWTH-CURVES/00_SINGLE CURVES/00-DATA/merge_data_frames_v2/151016_151028_rosco-ctrl.csv")
rm=NULL
keep=NULL
redo=NULL
title=paste0("abs((xi/local median -1))<n*local cv\ndt=",dtav1*2+1)
l=0
{
for(n in unique(d$cond)){
  if(n=="ctrl"){
    limy=c(1000,4000)
  } else if (n=="rosco"){
    limy=c(2000,5000)
  }
  for(nn in unique(d$date[which(d$cond==n)])){
    for (i in unique(d$index_c[which(d$cond==n & d$date==nn & d$index_c<1000)])){
      
      temp=d[which(d$index_c==i & d$cond==n & d$date==nn),]
      if(nrow(temp)>30){
        t=table[which(table$cond==n & table$date==nn & table$index_c==i),]
        
      l=l+1
      temp$tnormi=c(0:(nrow(temp)-1))*10/60
      
      V1=list()
      SD1=list()
      delta=list()
      V2=list()
      V1plot=list()
      deltaplot=list()
      V2plot=list()
      
      for (p in c(1:length(dtav1))){
        V1[[p]]=rep(NA,nrow(temp))
        delta[[p]]=rep(NA,nrow(temp))
        SD1[[p]]=rep(NA,nrow(temp))
        #step 1 = averaging (V1) ####
        dt=dtav1[p]
        for (k in c((1+dt):(nrow(temp)-dt))){
          V1[[p]][k]=median(temp$Volume1[c((k-dt):(k+dt))],na.rm = TRUE)
          SD1[[p]][k]=sd(temp$Volume1[c((k-dt):(k+dt))],na.rm = TRUE)
        }
        seq=c(c(1:dt),c((nrow(temp)-dt+1):nrow(temp)))
        for (k in seq[1:(length(seq)/2)]){
          V1[[p]][k]=median(temp$Volume1[c(k:(k+dt*2))],na.rm = TRUE)
          SD1[[p]][k]=sd(temp$Volume1[c(k:(k+dt*2))],na.rm = TRUE)
        }
        for (k in seq[(length(seq)/2+1):length(seq)]){
          V1[[p]][k]=median(temp$Volume1[c((k-dt*2):k)],na.rm = TRUE)
          SD1[[p]][k]=sd(temp$Volume1[c((k-dt*2):k)],na.rm = TRUE)
        }
        
        # step 2 part 2 = outliers removal (V2) from (V1 +/- n*SD), the values kept are V0 ####
        V2[[p]]=list()
        V2res=NULL
        V2lab=NULL
        # abszolute:
        # mu=mean(delta[[p]],na.rm = TRUE)
        # sd=sd(delta[[p]],na.rm = TRUE)
        
        # relative:
        mu=0
        plots=list()
        for (pp in c(1:length(SD))){
          # dt=((dtf-1)/2)
          filter=SD[pp]
          
          V2[[p]][[pp]]=rep(NA,nrow(temp))
          val=which(abs(temp$Volume1/V1[[p]]-1)<V1[[p]]/SD1[[p]]*filter)
          V2[[p]][[pp]][val]=temp$Volume1[val]
          V2res=c(V2res,V2[[p]][[pp]])
          V2lab=c(V2lab,rep(pp,length(V2[[p]][[pp]])))
          plots[[pp]]=ggplot(data.frame(y=c(temp$Volume1,V2[[p]][[pp]]),x=rep(temp$tnormi,2),z=c(rep(1,nrow(temp)),rep(2,nrow(temp)))),
                             aes(x=x,y=y,colour=factor(z)))+
            geom_point(size=0.3)+
            scale_colour_manual(values=table_aes[[1]]$color[c(3,1)])+
            labs(x="time",y=paste0("filter:",filter,"* local CV"),title=paste0(temp$date[1],"-",temp$cond[1],"-",temp$index_c[1]))+
            scale_y_continuous(limits=limy)+
            scale_x_continuous(limits=c(0,24))+
            theme_pub()+
            theme(legend.position = "none")
          print(plots[[pp]])
        }
        a<-readline(prompt = "c=keep,n=remove,r=check")
        
        if(a=="r"){
          p1=ggplot(data.frame(y=c(temp$Volume1,V2[[p]][[pp]]),x=rep(temp$tnormi,2),z=c(rep(1,nrow(temp)),rep(2,nrow(temp)))),
                     aes(x=x,y=y,colour=factor(z)))+
            geom_point(size=0.3)+
            scale_colour_manual(values=table_aes[[1]]$color[c(3,1)])+
            labs(x="time",y=paste0("filter:",filter,"* local CV"),title=paste0(temp$date[1],"-",temp$cond[1],"-",temp$index_c[1]))+
            # scale_y_continuous(limits=limy)+
            scale_x_continuous(limits=c(0,24))+
            theme_pub()+
            theme(legend.position = "none")
          print(p1)
          a=readline(prompt = "c=keep,n=remove")
        } 
        if(a=="c"){
          keep=rbind(keep,t)
        } else if(a=="n"){
          rm=rbind(rm,t)
        } else {
          redo=rbind(redo,t)
        }
      }
      }
    }
  }
}
}
# dev.off()

# write.csv(rm,"rm_complete_table.csv",row.names = FALSE)
# write.csv(keep,"keep_complete_table.csv",row.names=FALSE)

d=d[-which(d$index %in% rm$index),]
# table<-keep


####################################################################################################################################
######################## WORKFLOW for the analysis of the paper ####################################################################

# 5.1) bins by initial size ####
# 5.1.1) import dataframes ####
setwd("E:/Analysis/00_GROWTH-CURVES/00_SINGLE CURVES/170515/1-make-data-frame/outliers/")
table=rbind(read.csv("keep_complete_table.csv",stringsAsFactors = FALSE),read.csv("keep_uncomplete_table.csv",stringsAsFactors = FALSE))
df=read.csv("E:/Analysis/00_GROWTH-CURVES/00_SINGLE CURVES/00-DATA/merge_data_frames_v2/151016_151028_rosco-ctrl.csv",stringsAsFactors = FALSE)

df=df[,c(1,6,7,9,10:12,17:21)]
colnames(df)[3]="Volume1"
d<-df
d<-d[which(d$index %in% table$index),]
# write.csv(d,"df_all-complete-uncomplete-minus-weirdcurves_raw.csv",row.names = FALSE)

# 5.1.1 recalculate growth speed version 1 of smoothing ####
# d=d[which(d$phase==1),]
dtvelocity=4 # will be multiplide by 2 +1
dtav2=3 # will be multiplide by 2 +1
QRFW = 0 # Quartile rangge filtering  (distance from the quartile*the IQR)
dtf = 11 # dt window quartile range filtering
dtav1=1 # multiplide by 2+1 dt of smoothing for first step
d<-df_all_raw
pdf(paste0("mtd1-dtav1-",dtav1,"_QRF",sub(pattern = "[.]",QRFW,replacement = "p"),"_dtav",dtav2,"_dtvelocity",dtvelocity,"_complete.pdf"),onefile=TRUE,height=1.84,width=4.4)
df=NULL

for (i in unique(table$index[which(table$index<1000)])){
  temp=d[which(d$index==i),]
  temp$tnormi=c(0:(nrow(temp)-1))*10/60
  if(temp$cond[1]=="ctrl"){
    limy=c(1000,4000)
  } else {
    limy=c(2000,5000)
  }
    if(nrow(temp)>(dtvelocity*2)){
      # step 1
      temp$V1=NA
      dtav=dtav1
      for(j in c((1+dtav):(nrow(temp)-dtav))){
        temp$V1[j]=mean(temp$Volume1[c((j-dtav):(j+dtav))],na.rm = TRUE)
      }
      # step 2 = outliers removal (V2) from (V1), the values kept are V0 ####
      # qurtile filtering
      dt=((dtf-1)/2)
      seq=which(is.na(temp$V1)==FALSE)
      temp$V2=NA
      for (j in c((1+dt):(length(seq)-dt))){
        hist=temp$V1[c((seq[(j-dt)]):(seq[(j+dt)]))]
        if(length(hist[which(is.na(hist)==FALSE)])>2){
          hist2=which(hist<(quantile(hist,na.rm=TRUE)[[4]]+QRFW*IQR(hist,na.rm=TRUE)) & hist>(quantile(hist,na.rm=TRUE)[[2]]-QRFW*IQR(hist,na.rm=TRUE)))
          temp$V2[(c((seq[(j-dt)]):(seq[(j+dt)]))[hist2])]=temp$Volume1[(c((seq[(j-dt)]):(seq[(j+dt)]))[hist2])]
        } else {}
      }
      
      # step 3 = averaging on V0 - points removed on V2
      temp$av=NA
      dtav=dtav2
      for(j in c((1+dtav):(nrow(temp)-dtav))){
        temp$av[j]=mean(temp$V2[c((j-dtav):(j+dtav))],na.rm = TRUE)
      }
      
      k=dtvelocity
      m=data.frame(a=rep(NA,nrow(temp)),b=rep(NA,nrow(temp)),c=rep(NA,nrow(temp)),d=rep(NA,nrow(temp)),dvdt=rep(NA,nrow(temp)),R2=rep(NA,nrow(temp)),dvdt2=rep(NA,nrow(temp)))
      colnames(m)=c(paste0("x1-av-",k),paste0("y1-av-",k),paste0("x2-av-",k),paste0("y2-av-",k),paste0("advdt-av-",k),paste0("R2-av-",k),paste0("bdvdt-av-",k))
      temp<-cbind(temp,m)
      for (j in c((1+dt):(nrow(temp)-dt))){
        x=temp$tnormi[c((j-dt):(j+dt))]
        y=temp$av[c((j-dt):(j+dt))]
        temp[j,(ncol(temp)-2)]=lmrob(y~x)$coeff[[2]]
        temp[j,(ncol(temp))]=lmrob(y~x)$coeff[[1]]
        # temp[j,(ncol(temp)-1)]=summary(lmrob(y~x))$r.squared
        temp[j,(ncol(temp)-1)]=summary(lmrob(y~x))$sigma
        # temp[j,(ncol(temp)-6)]=temp$tnormi[(j-dt)] #x1
        temp[j,(ncol(temp)-6)]=temp$tnormi[j]-0.5 #x1
        temp[j,(ncol(temp)-5)]=temp[j,(ncol(temp)-2)]*temp[j,(ncol(temp)-6)]+temp[j,(ncol(temp))]
        # temp[j,(ncol(temp)-4)]=temp$tnormi[(j+dt)] #x1
        temp[j,(ncol(temp)-4)]=temp$tnormi[j]+0.5
        temp[j,(ncol(temp)-3)]=temp[j,(ncol(temp)-2)]*temp[j,(ncol(temp)-4)]+temp[j,(ncol(temp))]
      }
      
      df=rbind(df,temp)
        
      p1=ggplot(data.frame(x=rep(temp$tnormi,2),y=c(temp$Volume1,temp$V2),z=c(rep(1,nrow(temp)),rep(2,nrow(temp)))),
                aes(x=x,y=y,colour=factor(z)))+
        geom_point(size=0.3)+
        scale_colour_manual(values=table_aes[[1]]$color[c(3,1)])+
        labs(x="time",y=paste0("volume"),title=paste0(temp$date[1],"-",temp$cond[1],"-",temp$index_c[1]))+
        scale_y_continuous(limits=limy)+
        scale_x_continuous(limits=c(0,24))+
        theme_pub()+
        theme(legend.position = "none")
      p2=ggplot(temp,aes(x=tnormi,y=Volume1))+
        geom_point(size=0.3)+
        geom_segment(aes(x=temp$`x1-av-4`,y=temp$`y1-av-4`,xend=temp$`x2-av-4`,yend=temp$`y2-av-4`,colour=factor(tnormi)),size=0.1)+
        labs(x="time",y=paste0("volume"),title=paste0(temp$date[1],"-",temp$cond[1],"-",temp$index_c[1]))+
        # scale_y_continuous(limits=limy)+
        # scale_x_continuous(limits=c(0,24))+
        theme_pub()+
        theme(legend.position = "none")
      grid.arrange(arrangeGrob(grobs=list(p1,p2),nrow=1))
    }
    }

dev.off()

# write.csv(table,paste0("table_all",sub(pattern = "[.]",QRFW,replacement = "p"),"_dtav",dtav2,"_dtvelocity",dtvelocity,".csv"),row.names = FALSE)
write.csv(df,paste0("mtd1_","qr", gsub(pattern = "[.]",QRFW,replacement = "p"),"_dtav",dtav2,"_dtvelocity",dtvelocity,"_complete.csv"),row.names=FALSE)

# 5.1.1 recalculate growth speed, version 2 of smoothing ####
dtvelocity=4
dtav2=3
dtav1=5
SD=0.004
df=NULL
d<-df_all_raw
pdf(paste0("mtd2-dtav1-",dtav1,"_SD-filtering",sub(pattern = "[.]",SD,replacement = "p"),"_dtav",dtav2,"_dtvelocity",dtvelocity,"_complete.pdf"),onefile=TRUE,height=1.84,width=4.4)

for (i in unique(table$index[which(table$index<1000)])){
  temp=d[which(d$index==i),]
  temp$tnormi=c(0:(nrow(temp)-1))*10/60
  if(temp$cond[1]=="ctrl"){
    limy=c(1000,4000)
  } else {
    limy=c(2000,5000)
  }
  
  # step1 = smoothing version 2:
  dt=dtav1
  temp$V1=NA
  temp$SD1=NA
  temp$V2=NA
  for (k in c((1+dt):(nrow(temp)-dt))){
    temp$V1[k]=median(temp$Volume1[c((k-dt):(k+dt))],na.rm = TRUE)
    temp$SD1[k]=sd(temp$Volume1[c((k-dt):(k+dt))],na.rm = TRUE)
  }
  seq=c(c(1:dt),c((nrow(temp)-dt+1):nrow(temp)))
  for (k in seq[1:(length(seq)/2)]){
    temp$V1[k]=median(temp$Volume1[c(k:(k+dt*2))],na.rm = TRUE)
    temp$SD1[k]=sd(temp$Volume1[c(k:(k+dt*2))],na.rm = TRUE)
  }
  for (k in seq[(length(seq)/2+1):length(seq)]){
    temp$V1[k]=median(temp$Volume1[c((k-dt*2):k)],na.rm = TRUE)
    temp$SD1[k]=sd(temp$Volume1[c((k-dt*2):k)],na.rm = TRUE)
  }
  
  # step 2 outliers removal (V2) from (V1 +/- n*SD), the values kept are V0 ####
  filter=SD
  val=which(abs(temp$Volume1/temp$V1-1)<temp$V1/temp$SD1*filter)
  temp$V2[val]=temp$Volume1[val]
  
  # step 3 smoothing for the velocity calculation
  temp$av=NA
  dtav=dtav2
  for(j in c((1+dtav):(nrow(temp)-dtav))){
    temp$av[j]=mean(temp$V2[c((j-dtav):(j+dtav))],na.rm = TRUE)
  }
  
  k=dtvelocity
  m=data.frame(a=rep(NA,nrow(temp)),b=rep(NA,nrow(temp)),c=rep(NA,nrow(temp)),d=rep(NA,nrow(temp)),dvdt=rep(NA,nrow(temp)),R2=rep(NA,nrow(temp)),dvdt2=rep(NA,nrow(temp)))
  colnames(m)=c(paste0("x1-av-",k),paste0("y1-av-",k),paste0("x2-av-",k),paste0("y2-av-",k),paste0("advdt-av-",k),paste0("R2-av-",k),paste0("bdvdt-av-",k))
  temp<-cbind(temp,m)
  for (j in c((1+dt):(nrow(temp)-dt))){
    x=temp$tnormi[c((j-dt):(j+dt))]
    y=temp$av[c((j-dt):(j+dt))]
    temp[j,(ncol(temp)-2)]=lmrob(y~x)$coeff[[2]]
    temp[j,(ncol(temp))]=lmrob(y~x)$coeff[[1]]
    # temp[j,(ncol(temp)-1)]=summary(lmrob(y~x))$r.squared
    temp[j,(ncol(temp)-1)]=summary(lmrob(y~x))$sigma
    # temp[j,(ncol(temp)-6)]=temp$tnormi[(j-dt)] #x1
    temp[j,(ncol(temp)-6)]=temp$tnormi[j]-0.5 #x1
    temp[j,(ncol(temp)-5)]=temp[j,(ncol(temp)-2)]*temp[j,(ncol(temp)-6)]+temp[j,(ncol(temp))]
    # temp[j,(ncol(temp)-4)]=temp$tnormi[(j+dt)] #x1
    temp[j,(ncol(temp)-4)]=temp$tnormi[j]+0.5
    temp[j,(ncol(temp)-3)]=temp[j,(ncol(temp)-2)]*temp[j,(ncol(temp)-4)]+temp[j,(ncol(temp))]
  }
  
  df=rbind(df,temp)
  
  p1=ggplot(data.frame(x=rep(temp$tnormi,2),y=c(temp$Volume1,temp$V2),z=c(rep(1,nrow(temp)),rep(2,nrow(temp)))),
            aes(x=x,y=y,colour=factor(z)))+
    geom_point(size=0.3)+
    scale_colour_manual(values=table_aes[[1]]$color[c(3,1)])+
    labs(x="time",y=paste0("volume filtered"),title=paste0(temp$date[1],"-",temp$cond[1],"-",temp$index_c[1]))+
    scale_y_continuous(limits=limy)+
    scale_x_continuous(limits=c(0,24))+
    theme_pub()+
    theme(legend.position = "none")
  p2=ggplot(temp,aes(x=tnormi,y=Volume1))+
    geom_point(size=0.3)+
    geom_segment(aes(x=temp$`x1-av-4`,y=temp$`y1-av-4`,xend=temp$`x2-av-4`,yend=temp$`y2-av-4`,colour=factor(tnormi)),size=0.1)+
    labs(x="time",y=paste0("volume filterd + smoothed"),title=paste0(temp$date[1],"-",temp$cond[1],"-",temp$index_c[1]))+
    # scale_y_continuous(limits=limy)+
    # scale_x_continuous(limits=c(0,24))+
    theme_pub()+
    theme(legend.position = "none")
  grid.arrange(arrangeGrob(grobs=list(p1,p2),nrow=1))
}
dev.off()
write.csv(df,paste0("mtd2_","dtav1-",dtav1,"_SD-filtering",sub(pattern = "[.]",SD,replacement = "p"),"_dtav",dtav2,"_dtvelocity",dtvelocity,"_complete.csv"),row.names = FALSE)



# 5.1.1.2 make dataframe and table with the proper G1 inforamtions ####
setwd("E:/Analysis/00_GROWTH-CURVES/00_SINGLE CURVES/170515/1-make-data-frame/")
d=rbind(read.csv("mtd2_dtav1-5_SD-filtering0p004_dtav3_dtvelocity4_uncomplete.csv",stringsAsFactors = FALSE),
        read.csv("mtd2_dtav1-5_SD-filtering0p004_dtav3_dtvelocity4_complete.csv",stringsAsFactors = FALSE))
# write.csv(d,"mtd2_dtav1-5_SD-filtering0p004_dtav3_dtvelocity4_complete-uncomplete.csv",row.names = FALSE)
table=read.csv("outliers/all_kept_table.csv",stringsAsFactors = FALSE)
# write.csv(table,"table_all_info.csv",row.names = FALSE)
res=NULL
table$vbirth=NA
table$vm=NA
table$vg1=NA
d$phase=NA
d$tnormg1=NA
weird=NA
for (i in table$index){
  t=table[which(table$index==i),]
  if(is.na(t$tbirth)==FALSE){
    t$vbirth=d$av[which(d$index==i & is.na(d$av)==FALSE)[1]]
  }
  if(is.na(t$tg1)==FALSE && t$tg1%in% d$Frame[which(d$index==i)]){
    t$vg1=d$av[which(d$index==i & d$Frame==t$tg1)]
    d$tnormg1[which(d$index==i)]=(d$Frame[which(d$index==i)]-t$tg1)*10/60
    if(length(which(d$index==i & d$Frame<=t$tg1))>0){
      d$phase[which(d$index==i & d$Frame<=t$tg1)]=1
    }
    if(length(which(d$index==i & d$Frame>t$tg1))>0){
      d$phase[which(d$index==i & d$Frame>t$tg1)]=2
    }
  }
  if(is.na(t$tm)==FALSE && t$tm %in% d$Frame[which(d$index==i)]){
    t$vm=d$av[max(which(d$index==i & d$Frame<=t$tm & is.na(d$av)==FALSE))]
    if(max(d$Frame[which(d$index==i)])>t$tm){
      weird=c(weird,i)
      print(length(which(d$index==i & d$Frame>t$tm)))
      # d=d[-which(d$index==i && d$Frame>t$tm),]
    } else {}
  }
  res=rbind(res,t)
}
# res[which(res$index%in% weird),which(colnames(res)%in%c("tg1","tm"))]=NA
write.csv(res,"E:/Analysis/00_GROWTH-CURVES/00_SINGLE CURVES/170515/1-make-data-frame/table_all_info.csv",row.names = FALSE)
table<-res
write.csv(d,"E:/Analysis/00_GROWTH-CURVES/00_SINGLE CURVES/170515/1-make-data-frame/mtd2_dtav1-5_SD-filtering0p004_dtav3_dtvelocity4_complete-uncomplete.csv",row.names = FALSE)


# 5.1.1.3 vizualize all single plots ####
p=list()
l=0
for (i in table$index[which(table$index_c<1000 & table$cond=="ctrl")]){
  temp=d[which(d$index==i),]
  if(temp$cond[1]=="ctrl"){
    limy=c(1000,4000)
  } else {
    limy=c(2000,5000)
  }
  
  l=l+1
  p[[l]]=ggplot(temp,aes(x=tnormi,y=V2))+
    geom_point(size=0.3)+
    labs(x="time",y="volume (mtd2)",title=paste0(temp$cond[1],"-",temp$date[1],"-",i))+
    scale_y_continuous(limits=limy)+
    scale_x_continuous(limits=c(0,24))+
    theme_pub()
  if(is.na(table$tg1[which(table$index==i)])==FALSE){
    p[[l]]=p[[l]]+
      geom_vline(xintercept = temp$tnormi[which(temp$Frame==table$tg1[which(table$index==i)])],size=0.3)
  }
}
ncol=2
nrow=ceiling(length(p)/ncol)
pdf("E:/Analysis/00_GROWTH-CURVES/00_SINGLE CURVES/170515/1-make-data-frame/all_ctrl_plots_complete_smooth-mtd2.pdf",height=1.84*nrow,width=2.2*ncol)
grid.arrange(arrangeGrob(grobs=p,nrow=nrow))
dev.off()


# 5.1.2 bins ####

setwd("E:/Analysis/00_GROWTH-CURVES/00_SINGLE CURVES/170515/hyp4_growth speed vs vbirth/")
df=read.csv("E:/Analysis/00_GROWTH-CURVES/00_SINGLE CURVES/170515/1-make-data-frame/mtd2_dtav1-5_SD-filtering0p004_dtav3_dtvelocity4_complete-uncomplete.csv",stringsAsFactors = FALSE)
df=df[-which(df$index==1022),]
df=df[which(df$index %in%table$index[which(is.na(table$vbirth)==FALSE)]),]

# first method absolute bins ####
# bins initial size:
hist(df$V2[which(df$tnormi==df$tnormi[4])])
bins=c(1500,1900)
binnames=c("<1500","1500-1900",">1900")
bintitle="volume initial"
laby="dvdt"
# laby=expression(paste("<",delta,"v/",delta,"t >"," (µm"^3,".hrs"^-1,")"))

# bins inital dvdti:
hist(df$dvdti[which(df$tnormi==df$tnormi[5])])
bins=c(-50,50,150)
binnames=c("<-50","50-50","50-150",">150")
bintitle="dvdt initial"

# bins by dt G1:
bins=c(4,6,8)
binnames=c("<4","4-6","6-8",">8")
bintitle="DT G1"
df$dtg1=-df$tnormg1

sort_bin="V2"
df$z=NA
C=which(colnames(df)==sort_bin)
for (j in unique(df$index)){
  val=df[which(df$index==j & is.na(df[,C])==FALSE)[1],C]
  b=1
  while((val>bins[b])==TRUE & b<=length(bins)){
    b=b+1
  }
  df$z[which(df$index==j)]=b
}
# df$dvdt_initial=df$z
# df$bin_dtg1=df$z
df$bin_vinitial=df$z


# bins all timepoints together by total size
hist(df$V2)
bins=seq(1300,4000,200)
binnames=NULL
sort_bin="V2"
df$z=NA
df=df[-which(is.na(df$V2)==TRUE),]
C=which(colnames(df)==sort_bin)
for (i in c(1:(length(bins)+1))){
  if(i==(length(bins)+1)){
    list=which(df[,C]>=bins[length(bins)])
    binnames=c(binnames,paste0(">",bins[(i-1)]))
  } else if(i==1){
    list=which(df[,C]<=bins[1])
    binnames=c(binnames,paste0("<",bins[i]))
  } else {
    list=which(df[,C]<=bins[i] & df[,C]>bins[(i-1)])
    binnames=c(binnames,paste0(bins[(i-1)],"-",bins[i]))
  }
  df$z[list]=i
}
binnames1<-binnames
df$bins_inst_size=df$z


temp=data.frame(y=df$advdt.av.4,z=df$z,z2=df$tnormi,z3=df$tnormg1,av1=df$V2,index=df$index)
temp=temp[which(temp$z2<5),]
temp=temp[which(temp$z3>-5),]


temp=data.frame(y=df$advdt.av.4,z=df$bins_inst_size,z3=df$bin_vinitial,tnormi=df$tnormi,tnormg1=df$tnormg1,cond=df$cond)
# write.csv(temp,"df_bin_vbirth_vinst.csv",row.names = FALSE)

# calculate x y coordinates of the boxplot ####
d<-temp
for (i in unique(d$cond)){
  for (j in sort(unique(d$z3))){
    temp=d[which(d$cond==i & d$z3==j),]
    for(k in sort(unique(temp$y))){
      mean=c(mean,mean(temp$y[which(temp$z==k)],na.rm=TRUE))
      sd=c(sd,sd(temp$y[which(temp$z==k)],na.rm=TRUE))
      n=c(n,length(temp$y[which(temp$z==k)]))
      x=c(x,k)
    }
  }
}


# new method fo bins choosing ####
# bins initial size:
cond=c("ctrl","rosco")
table=read.csv("E:/Analysis/00_GROWTH-CURVES/00_SINGLE CURVES/170515/1-make-data-frame/table_all_info.csv",stringsAsFactors = FALSE)
table=table[which(is.na(table$vbirth)==FALSE),]
table=table[-which(table$index==1022),]
bins=list()
for (i in c(1:length(cond))){
  x=table$vbirth[which(table$cond==cond[i])]
  m=mean(table$vbirth[which(table$cond==cond[i])],na.rm = TRUE)
  # bins[[i]]=c((m-0.3*m),(m+0.3*m))
  bins[[i]]=c(quantile(x,0.2)[[1]],quantile(x,0.8)[[1]])
}
d=list()
for (i in c(1:length(cond))){
  temp=df[which(df$index %in% table$index[which(table$cond==cond[i] & table$vbirth<=bins[[i]][1])]),]
  temp$bin_vinitial=1
  temp2=df[which(df$index %in% table$index[which(table$cond==cond[i] & table$vbirth>=bins[[i]][2])]),]
  temp2$bin_vinitial=2
  temp3=df[which(df$index %in% table$index[which(table$cond==cond[i] & !(table$index %in% unique(c(temp$index,temp2$index))))]),]
  temp3$bin_vinitial=3
  d[[i]]=rbind(temp,temp2,temp3)
}

# for grouped measurement:
x=table$vbirth
bins[[3]]=c(quantile(x,0.2)[[1]],quantile(x,0.8)[[1]])

cond[3]="all"
i=3
temp=df[which(df$index %in% table$index[which(table$vbirth<=bins[[i]][1])]),]
temp$bin_vinitial=1
temp2=df[which(df$index %in% table$index[which(table$vbirth>=bins[[i]][2])]),]
temp2$bin_vinitial=2
temp3=df[which(df$index %in% table$index[which(!(table$index %in% unique(c(temp$index,temp2$index))))]),]
temp3$bin_vinitial=3
d[[i]]=rbind(temp,temp2,temp3)

# calculate the slope for each subgroup #####
d=df[which(df$cond=="ctrl"),]
table_aes_ggplot=table_aes[[1]]
table_aes_ggplot$size=0.3
# table_aes_ggplot$color=rev(c('#939393','#d3d3d3',"#d7d7d7"))


for (limindex in c(0,3,5,8)){
  results=NULL
# plot1=plot_sgp_bm_by_cond(d[which(df$phase==1),],"advdt.av.4","av","dvdt","volume","dv dt vs v bin by vi",sort_cond = "bin_vinitial",eb = FALSE,sp = FALSE,"index",limindex = limindex)$plot1_pub+
#   scale_y_continuous(limits=c(0,200))+
#   scale_x_continuous(limits=c(1200,3500))
  plot1=plot_sgp_bm_by_cond(df[which(df$cond=="ctrl" & df$phase==1),],"advdt.av.4","av","dvdt","volume","dv dt vs v bin by vi",sort_cond = "bin_vinitial",eb = FALSE,sp = FALSE,"index",limindex = limindex)$plot1_pub+
    scale_y_continuous(limits=c(0,200))+
    scale_x_continuous(limits=c(1200,3500))
  plot2=plot_sgp_bm_by_cond(df[which(df$cond=="rosco" & df$phase==1),],"advdt.av.4","av","dvdt","volume","dv dt vs v bin by vi",sort_cond = "bin_vinitial",eb = FALSE,sp = FALSE,"index",limindex = limindex)$plot1_pub+
    scale_y_continuous(limits=c(0,200))+
    scale_x_continuous(limits=c(1200,3500))
pdf(paste0("G1_dvdt-vs-vi_sep_limindex-",limindex,".pdf"),height=1.84,width=4*2)
grid.arrange(plot1,plot2,nrow=1)
# print(plot1)
dev.off()

m=plot_sgp_bm_by_cond(d[which(df$cond=="ctrl" & df$phase==1),],"advdt.av.4","av","dvdt","volume","dv dt vs v bin by vi",sort_cond = "bin_vinitial",eb = FALSE,sp = FALSE,"index",limindex = limindex)$table_xlsx
m$cond="ctrl"
results=rbind(results,m)

m=plot_sgp_bm_by_cond(d[which(df$cond=="rosco" & df$phase==1),],"advdt.av.4","av","dvdt","volume","dv dt vs v bin by vi",sort_cond = "bin_vinitial",eb = FALSE,sp = FALSE,"index",limindex = limindex)$table_xlsx
m$cond="rosco"
results=rbind(results,m)

results$limindex=limindex
write.xlsx(x = results,paste0("G1_dvdt-vs-vi_separated_limindex-",limindex,".xlsx"),sheetName = "1")
}

labels=c("<25%",">75%","25-75%")
labels=c('<20%',">80%","20-80")
for (limindex in c(0,2,4)){
  plots=list()
  results=NULL
  for (i in c(1:length(cond))){
    plots[[i]]=plot_sgp_bm_by_cond(d[[i]][which(d[[i]]$phase==1),],"advdt.av.4","av","dvdt","volume","dv dt vs v bin by vi",sort_cond = "bin_vinitial",eb = FALSE,sp = FALSE,"index",limindex = limindex,labs=labels)$plot1_pub+
      scale_y_continuous(limits=c(0,200))+
      scale_x_continuous(limits=c(1200,3500))+
      ggtitle(paste0("- bins with <",limindex, "\ndifferent cells"))
    m=plot_sgp_bm_by_cond(d[[i]][which(d[[i]]$phase==1),],"advdt.av.4","av","dvdt","volume","dv dt vs v bin by vi",sort_cond = "bin_vinitial",eb = FALSE,sp = FALSE,"index",limindex = limindex,labs=labels)$table_xlsx
    m$cond=cond[i]
    results=rbind(results,m)
    
  }
  pdf(paste0("G1_dvdt-vs-vi_sep_limindex-",limindex,".pdf"),height=2,width=2.6*2)
  grid.arrange(plots[[1]],plots[[2]],nrow=1)
  dev.off()
  results$limindex=limindex
  write.xlsx(x = results,paste0("G1_dvdt-vs-vi_separated_limindex-",limindex,".xlsx"),sheetName = "1")
}

l=0
i=3
results=NULL
for(limindex in c(0,5,8)){
  l=l+1
  plots[[l]]=plot_sgp_bm_by_cond(d[[i]][which(d[[i]]$phase==1),],"advdt.av.4","av","dvdt","volume","dv dt vs v bin by vi",sort_cond = "bin_vinitial",eb = FALSE,sp = FALSE,"index",limindex = limindex,labs=labels)$plot1_pub+
    scale_y_continuous(limits=c(0,200))+
    scale_x_continuous(limits=c(1200,3500))+
    ggtitle(paste0("- bins with <",limindex, "\ndifferent cells"))
  m=plot_sgp_bm_by_cond(d[[i]][which(d[[i]]$phase==1),],"advdt.av.4","av","dvdt","volume","dv dt vs v bin by vi",sort_cond = "bin_vinitial",eb = FALSE,sp = FALSE,"index",limindex = limindex,labs=labels)$table_xlsx
  m$cond=cond[i]
  m$limindex=limindex
  results=rbind(results,m)
}
  pdf(paste0("G1_dvdt-vs-vi_together_v2",".pdf"),height=1.84,width=2.6*length(plots))
  grid.arrange(arrangeGrob(grobs=plots,nrow=1))
  dev.off()
  write.xlsx(x = results,paste0("G1_dvdt-vs-vi_together_v2",".xlsx"),sheetName = "1")

# alternatively, to plot everything by yourself: ####
  i=3
  lim=0
  a=3
  limindex=5
da=plot_sgp_bm_by_cond(d[[i]][which(d[[i]]$phase==1),],"advdt.av.4","av","dvdt","volume","dv dt vs v bin by vi",sort_cond = "bin_vinitial",eb = FALSE,sp = FALSE,"index",limindex = limindex,labs=labels)$mat
table_xlsx=plot_sgp_bm_by_cond(d[[i]][which(d[[i]]$phase==1),],"advdt.av.4","av","dvdt","volume","dv dt vs v bin by vi",sort_cond = "bin_vinitial",eb = FALSE,sp = FALSE,"index",limindex = limindex,labs=labels)$table_xlsx

temp=NULL
for (i in c(1:(length(da)-1))){
  temp=rbind(temp,da[[i]])
}
# temp=temp[-which(temp$`n index`<limindex),]
plot=ggplot(temp[-which(temp$`n index`<limindex),],aes(x=x,y=y,colour=factor(zz)))+
  geom_point(size=1.5)+
  scale_colour_manual(values=table_aes_ggplot$color[c(1:length(unique(temp$zz)))],labels=labels)+
  theme_pub()
temp=temp[-which(temp$`n index`<limindex),]
sink("summary_lmrob_dvdt-vs-v-bin-vi_limindex5.txt")
for (j in c(1:length(unique(temp$zz)))){
  i=unique(temp$zz)[j]
  x=temp$x[which(temp$zz==i)]
  y=temp$y[which(temp$zz==i)]
  plot=plot+
    geom_abline(slope=lmrob(y~x)$coeff[[2]],intercept = lmrob(y~x)$coeff[[1]],colour=table_aes_ggplot$color[j],size=0.3)
  print(i)
  print(summary(lmrob(y~x)))
}
sink()
setEPS(width=4,height=1.6,family="ArialMT")
postscript(file=paste0("G1_dvdt-vs-vi_together_v1","20-80quartiles","_limindex",limindex,".eps"),height=1.6,width=4,family="ArialMT")
print(plot)
dev.off()

# histogram of volume birth distribution
dens=density(table$vbirth)
df <- data.frame(x=dens$x, y=dens$y)
# probs <- c(0.1, 0.25, 0.5, 0.75, 0.9)
probs <-c(0.2,0.8)
quantiles <- quantile(table$vbirth, prob=probs)
df$quant <- factor(findInterval(df$x,quantiles))
hist=ggplot(df, aes(x,y)) + 
  geom_line() + 
  geom_ribbon(aes(ymin=0, ymax=y, fill=quant)) + 
  scale_x_continuous(breaks=quantiles) + 
  scale_fill_brewer(guide="none")+
  theme_pub()+
  theme(aspect.ratio=0.4)
setEPS(width=4,height=(1.6/2),family="ArialMT")
postscript(file=paste0("G1_dvdt-vs-vi_together_v1","20-80quartiles","_histogram.eps"),height=(1.6/2),width=4,family="ArialMT")
print(hist)
dev.off()

hist=ggplot(table,aes(x=vbirth))+
  geom_histogram(aes(y=..density..), binwidth=1, colour="black", fill="white") + 
  geom_density(fill=NA, colour="royalblue",aes(y = ..count..))

# plots ####
# lineplot
p<-ggplot(temp,aes(y=y,x=av1,group=interaction(factor(index),factor(z)),colour=factor(z)))+
  geom_line(size=0.1)+theme_pub()+
  scale_colour_manual(values=table_aes_ggplot$color[c(1:length(unique(temp$z)))],name=bintitle,labels=binnames[sort(unique(temp$z))])+
  coord_cartesian(ylim=c(-200,400))+
  guides(index=FALSE)


# boxplot
p<-ggplot(temp,aes(factor(z),y,group=interaction(factor(z),factor(z3)),fill=factor(z),alpha=factor(z3)))+
  geom_boxplot(outlier.size=0.1,size=0.1)+
  scale_fill_manual(values=table_aes[[1]]$color[1:(length(bins)+1)])+
  scale_x_discrete(labels=binnames)+
  coord_cartesian(ylim=c(-200,400))+
  labs(x="initial volume",
       y=laby,
       title="All G1 phase, rosco and ctrl together, bin by inital-dvdt")+
  theme_pub()+
  theme(legend.position="none",
        aspect.ratio=0.3)+
  geom_hline(yintercept = 100,size=0.3,linetype=3)+
  geom_hline(yintercept = 50,size=0.3,linetype=3)

p<-ggplot(temp,aes(factor(z),y,group=interaction(factor(z),factor(z3)),fill=factor(z3)))+
  geom_boxplot(outlier.size=0.1,size=0.1)+
  scale_fill_manual(values=table_aes[[1]]$color[1:(length(bins)+1)],labels=binnames,name=bintitle)+
  scale_x_discrete(labels=binnames1)+
  coord_cartesian(ylim=c(-200,400))+
  labs(x="volume",
       y=laby,
       title="All G1 phase, rosco and ctrl together, bin by dt g1")+
  theme_pub()+
  theme(
        aspect.ratio=0.3)+
  geom_hline(yintercept = 100,size=0.3,linetype=3)+
  geom_hline(yintercept = 50,size=0.3,linetype=3)

p<-ggplot(temp,aes(factor(z),y,fill=factor(z3),alpha=factor(cond)))+
  geom_boxplot(outlier.size=0.1,size=0.1)+
  scale_fill_manual(values=table_aes[[1]]$color[1:(length(bins)+1)])+
  scale_x_discrete(labels=binnames1)+
  coord_cartesian(ylim=c(-200,400))+
  labs(x="initial volume",
       y=laby,
       title="All G1 phase, rosco and ctrl together, bin by inital-vi")+
  theme_pub()+
  theme(aspect.ratio = 0.3)#+
  # theme(legend.position="none",
        # aspect.ratio=0.3)+
  # geom_hline(yintercept = 100,size=0.3,linetype=3)+
  # geom_hline(yintercept = 50,size=0.3,linetype=3)


pdf(paste0("dvdt-vs-v_grouped by_vbirth.pdf"),height=2*1,width=9)
# grid.arrange(p1,p2,p3,nrow=3)
p
dev.off()









########### modif functions for plots #######################
plot_sgp_bm_by_cond<-function(df,paramy,paramx,title_paramy,title_paramx,title_plot,sort_cond,eb,sp,index,limindex,labs){
  K=which(colnames(df)==paramy)
  L=which(colnames(df)==paramx)
  M=which(colnames(df)==sort_cond)
  N=which(colnames(df)==index)
  x=df[,L][which(is.na(df[,K])==F & is.na(df[,L])==F)]
  y=df[,K][which(is.na(df[,K])==F & is.na(df[,L])==F)]
  z=df[,M][which(is.na(df[,K])==F & is.na(df[,L])==F)]
  index=df[,N][which(is.na(df[,K])==F & is.na(df[,L])==F)]
  d = data.frame(x,y,z,index)
  if(is.numeric(d$z)==FALSE){
    d$z=as.character(d$z)
  }
  d= d [with(d,order(z)),]
  nobs=nrow(d)
  zlist=unique(d$z)
  
  # set aesthetics for the plot ####
  barheight=abs(max(d$y)-min(d$y))/40
  barwidth=abs(max(d$x)-min(d$x))/40
  temp=c(min(d$y),min(d$x),max(d$y),max(d$x))
  coord=c(round(min(temp),-2)-100,round(max(temp),-2)+100)
  
  # calculate the binned values for each condition
  dsave<-d
  d$n=NA
  d$sdx=NA
  d$semx=NA
  d$sdy=NA
  d$semy=NA
  d$zz=NA
  d$`n index`=NA
  da=list()
  binw_list=list()
  m=0
  shapes=NULL
  colors=NULL
  sizes=NULL
  fill=NULL
  linetype=NULL
  labels=NULL
  limindex=limindex
  lim=1
  
  binn=c(round(min(d$x,na.rm=TRUE),-2),100)
  for (i in c(1:(length(zlist)+1))){
    if(i<length(zlist)+1){
      xtest=dsave$x[which(dsave$z==zlist[i])]
      ytest=dsave$y[which(dsave$z==zlist[i])]
      index=dsave$index[which(dsave$z==zlist[i])]
      da[[i]]=calculate_mean_bin_modif(xtest,ytest,binn,mtd=4,index)$da
      binw_list[[i]]=calculate_mean_bin_modif(xtest,ytest,binn,mtd=4,index)$binw
      # da[[i]]=da[[i]][which(is.na(da[[i]]$x)==FALSE),]
      da[[i]]=da[[i]][which(da[[i]]$`n index`>limindex),]
      da[[i]]$index=NA
      da[[i]]$z=zlist[i]
      da[[i]]$zz=m
      d$zz[which(d$z==zlist[i])]=m+1
      m=m+2
      # colors=c(colors,color_stat,table_aes_ggplot$color[i])
      colors=c(colors,table_aes_ggplot$color_stat[i],table_aes_ggplot$color[i])
      shapes=c(shapes,table_aes_ggplot$shape_bin[i],table_aes_ggplot$shape_scp[i])
      sizes=c(sizes,0.5,table_aes_ggplot$size[i])
      fill=c(fill,table_aes_ggplot$fill_bin[i],table_aes_ggplot$color[i])
      labels=c(labels,"median bins",as.character(zlist[i]))
      linetype=c(linetype,table_aes_ggplot$linetype[i])
      d<-rbind(d,da[[i]])
    } else {
      xtest=dsave$x
      ytest=dsave$y
      index=dsave$index
      da[[i]]=calculate_mean_bin_modif(xtest,ytest,binn,mtd=4,index)$da
      da[[i]]=da[[i]][which(da[[i]]$`n index`>limindex),]
      da[[i]]$z="all"
      da[[i]]$zz=m
      da[[i]]$index=NA
      # for all together, we don't want to plot them on top of the others so we don't bind it with the dataframe d
    }
  }
  # mat with results for table ####
  mat=matrix(ncol=8,nrow=(length(zlist)+1))
  mat=as.data.frame(mat)
  for (i in c(1:(length(zlist)+1))){
    if(i<(length(zlist)+1)){
      mat[i,1]=signif(mean(dsave$y[which(dsave$z==zlist[i])]),a)
      mat[i,2]=signif(sd(dsave$y[which(dsave$z==zlist[i])]),a)
      mat[i,3]=signif(cor(dsave$x[which(dsave$z==zlist[i])],dsave$y[which(dsave$z==zlist[i])]),a)
      mat[i,4]=signif(cor.test(dsave$x[which(dsave$z==zlist[i])],dsave$y[which(dsave$z==zlist[i])])$p.value,a)
      mat[i,7]=length(dsave$y[which(dsave$z==zlist[i])])
      if(nrow(da[[i]])>3){
        testb=lm(da[[i]]$y~da[[i]]$x,weights=da[[i]]$n)
        testb$df.residual = with(testb,sum(weights)-length(coefficients))
        mat[i,8]=signif(summary(testb)$coefficients[2,4],a) # pval a
        if(is.na(mat[i,8])==FALSE && (mat[i,8]<0.05)==TRUE){
          mat[i,5]=signif(coef(testb)[2],a)
          mat[i,6]=signif(coef(testb)[1],a)
          rm(testb)
        } else {}
      }else{}
    } else {
      mat[i,1]=signif(mean(dsave$y),a)
      mat[i,2]=signif(sd(dsave$y),a)
      mat[i,3]=signif(cor(dsave$x,dsave$y),a)
      mat[i,4]=signif(cor.test(dsave$x,dsave$y)$p.value,a)
      mat[i,7]=length(dsave$y)
      if(nrow(da[[i]])>3){
        testb=lm(da[[i]]$y~da[[i]]$x,weights=da[[i]]$n)
        testb$df.residual = with(testb,sum(weights)-length(coefficients))
        mat[i,8]=signif(summary(testb)$coefficients[2,4],a)
        if(is.na(mat[i,8])==FALSE && (mat[i,8]<0.05)==TRUE){
          mat[i,5]=signif(coef(testb)[2],a)
          mat[i,6]=signif(coef(testb)[1],a)
        } else {}
      }else{}
    }
  }
  
  colnames(mat)=c("meany","sdy","Pearson","pval","alpha","beta","n","pval a")
  rownames(mat)=c(as.character(zlist),"all")
  
  da_list=da
  name_plot=title_plot
  m=table_stat(M=M,d = d,df = df,zlist = zlist, name_plot,mtd = "bm",binw = binw_list,da = da_list,lim = lim)
  m$n=m$n-m$`bin number<lim`
  m$mtd="lm_medbin_weight-n"
  dd<-d
  
  # plot general characteristics,
  if(sp==FALSE){
    d=d[which(d$zz %% 2 == 0),]
  } else {}
  
  plot1 <- ggplot(d,aes(x=x,y=y, color=factor(zz),shape=factor(zz),size=factor(zz),fill=factor(zz)))#,linetype=factor(zz)))
  
  zlist2=unique(d$z) ##so if there were not enough bins, it was removed
  
  if(length(labs)>0){
    labels<-labs
  }
  
  if(sp==TRUE){
    plot1<-plot1+
      geom_point() +
      labs(y=title_paramy, 
           x=title_paramx)+
      scale_colour_manual(values = colors, name=sort_cond,labels=labels) +
      scale_shape_manual(values = shapes,guide=FALSE) +
      scale_size_manual(values = sizes,guide=FALSE) +
      scale_fill_manual(values = fill, guide=FALSE) 
    # scale_linetype_manual(values=c(1,1,3,3))
  } else {
    colors=table_aes_ggplot$color_stat[1:length(zlist2)]
    shapes=table_aes_ggplot$shape_bin[1:length(zlist2)]
    fill=table_aes_ggplot$fill_bin[1:length(zlist2)]
    sizes=rep(3,length(zlist2))
    plot1<-plot1+
      geom_point(size=0.5) +
      labs(y=title_paramy, 
           x=title_paramx)+
      scale_colour_manual(values = colors, name=sort_cond,labels=labels) +
      scale_shape_manual(values = shapes,guide=FALSE) +
      scale_size_manual(values = sizes,guide=FALSE) +
      scale_fill_manual(values=fill, guide=FALSE)
  }
  if(eb==TRUE){
    for (i in c(1:length(zlist2))){
      dtemp=d[which(d$z==zlist2[i]),]
      plot1<-plot1+
        geom_errorbar(data=dtemp,aes(ymin=y-sdy, ymax=y+sdy), width=barwidth,colour=table_aes_ggplot$color_stat[i],linetype=1,size=0.25) +
        geom_errorbarh(data=dtemp,aes(xmin=x-sdx, xmax=x+sdx), height=barheight,colour=table_aes_ggplot$color_stat[i],linetype=1,size=0.25)
    }
  } else {}
  plot1_pub<-plot1+
    scale_x_continuous(expand = c(0.1, 0))+
    scale_y_continuous(expand=c(0.1,0))+
    theme_pub()
  plot1_png<-plot1+
    scale_x_continuous(expand = c(0.1, 0))+
    scale_y_continuous(expand=c(0.1,0))+
    theme_png()
  
  for (i in c(1:length(zlist))){
    if(is.na(mat[i,8])==FALSE && (mat[i,8]<0.05)==TRUE){ # test on pval alpha lmrob
      # if ((mat[i,4]<0.05)==TRUE){ # test on pval pearson
      plot1_pub<-plot1_pub+
        geom_abline(intercept = mat[i,6], slope = mat[i,5],size=0.3,linetype =1,colour=table_aes_ggplot$color_stat[i])
      plot1_png<-plot1_png+
        geom_abline(intercept = mat[i,6], slope = mat[i,5],size=0.45,linetype =1,colour=table_aes_ggplot$color_stat[i])
    }
  }
  
  
  return(list("plot1_pub"=plot1_pub,"plot1_png"=plot1_png,"limcoord"=temp,"table"=mat,"table_xlsx"=m,"dsave"=dd, "mat"=da))
} ### ok median bins for each condition

calculate_mean_bin_modif<-function(xtest,ytest,binn,mtd,index){
  if(mtd==1){
    binw=(max(xtest)-min(xtest))/binn # for sgp
  } else if (mtd==2) {
    binw=ceiling((max(xtest)-min(xtest))/binn) # for other plots
  } else if (mtd==3) {
    mm=mean(xtest)
    sd=sd(xtest)
    max=max(xtest[which(xtest<(3*sd+mm))])
    min=min(xtest[which(xtest>(3*sd-mm))])
    binw=(max-min)/binn
  } else if (mtd==4){
    bin=binn[1]+binn[2]
    binstart=binn[1]
    binw=binn[2]
    k=1
    while((binstart+binw)<max(xtest,na.rm=TRUE)){
      binstart=binstart+binw
      k=k+1
    }
    binn<-k
  }
  mat=matrix(nrow=binn,ncol=9)
  m=0
  dd=data.frame(xtest,ytest,index)
  colnames(dd)=c("x","y","index")
  if(mtd==3){
    bin=min(xtest[which(xtest>(3*sd-mm))])+binw
    maxmax=max(xtest[which(xtest<(3*sd+mm))])
  } else if(mtd==4){
    maxmax=binstart+binw
  } else {
    bin=min(dd$x)+binw ## uper limit of the interval considered
    maxmax=max(dd$x)
  }
  
  while(bin<=maxmax){
    m=m+1
    tempx=dd$x[which(dd$x<bin & dd$x>=(bin-binw))] # intervals are [;[
    tempy=dd$y[which(dd$x<bin & dd$x>=(bin-binw))]
    # print(tempy)
    # print(tempx)
    mat[m,1]=mean(tempx)
    mat[m,2]=mean(tempy)
    mat[m,4]=length(tempx)
    if (length(tempx)<lim){ # remove binnings that have a too low number of values
      mat[m,1]=NA
      mat[m,2]=NA
    } else {
      mat[m,5]=sd(tempx)
      mat[m,6]=sd(tempx)/sqrt(length(tempx))
      mat[m,7]=sd(tempy)
      mat[m,8]=sd(tempy)/sqrt(length(tempy))
      mat[m,9]=length(unique(dd$index[which(dd$x<bin & dd$x>=(bin-binw))]))
    }
    bin=bin+binw
  }
  da = as.data.frame(mat)
  colnames(da)=c("x","y","z","n","sdx","semx","sdy","semy","n index")
  return(list("da"=da,"binw"=binw))
}
