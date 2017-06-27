# INPUT = folder (dirin) with all the ".csv" files imported from the MATLAB FXM software
# OUTPUT = a data frame "DF1" with all the complete trahectories and "df2" with all the incomplete trajectories 
# + table1 and table2 (coresponding respectively to df1 and df2) with dditional informations about cell cycle progression/measurement at key timepoints (G1/S transition, birth and mitosis)
# the code alternates automated steps and steps where you need to vizualize single curves and mark the key points (a window pops up and guides you through the steps)
# complete trajectories are numbered from 1 to 999 ("index_c"), uncompomplete trahectories are  numbered from 1000 and on

# import single-trajectories from raw-automated tracking in Matlab ####
dirin="E:/Analysis/00_GROWTH-CURVES/00_SINGLE CURVES/00-DATA/151028-traj_HeLa-rosco_170201/tracks/"
dirout="E:/Analysis/00_GROWTH-CURVES/00_SINGLE CURVES/00-DATA/151028-traj_HeLa-rosco_170201/results_df/"
dir.create(dirout)
list=list.files(dirin)

#################### 1) MAKE the DATAFRAMES with Complete/Uncomplete trajectories ####
# at the end of this step, df1 has complete trajetories (including the overshoot),
# df2 has the uncomplet trajectories (without the overhdoot) + the indication of whether they should be align to "s" start","m" mitosis or "p" nothing in particular

df1=NULL
df2=NULL
# index = number of the cell ("_C" for complete trajectories and "_u" for uncomplete trajectories)
# index_c=0
# index_u=1000
# index_c=read.csv(paste0(dirout,"store_variables/",list.files(paste0(dirout,"store_variables/"))[which(grepl("index_c",list.files(paste0(dirout,"store_variables/")))==TRUE)]))[1,2]
# index_u=read.csv(paste0(dirout,"store_variables/",list.files(paste0(dirout,"store_variables/"))[which(grepl("index_u",list.files(paste0(dirout,"store_variables/")))==TRUE)]))[1,2]

list_analysis=15
df1=NULL
df2=NULL
list_redo=NULL
list_check_movie=NULL
list_check_index=NULL
x11()
{
  # for(i in c(2:length(list))){
  for (i in list_analysis){
    d<-read.csv(paste0(dirin,list[i]))#,stringsAsFactors = FALSE)
    colnames(d)[c(10:13)]=c("sum_im-median","sum_mask-median","median","mask")
      for (j in unique(d$CellIdInAutoTracking)){
      # for (j in unique(d$CellIdInAutoTracking)[c(13:length(unique(d$CellIdInAutoTracking)))]){
      temp=d[which(d$CellIdInAutoTracking==j),]
      if(nrow(temp)>40){
        y=temp$Cell_Intensisty
        x=temp$Frame
        plot.new()
        plot(x,y)
        # colnames(temp)[c(10:13)]=c("sum_im-median","sum_mask-median","median","mask")
        a=readline("keep, complete cellcycle = c, v=check family, check raw image=y, reject = n, made a mistake with the previous plot = p")
        if(a=="p") {
          list_redo=rbind(list_redo,d[which(d$CellIdInAutoTracking==(tempj)),])
          a=readline(" now, for this current plot: keep, complete cc = c, v=check family, check raw image=y, reject = n")
        }
        if(a=="c"){
          temp$index_c=index_c+1
          index_c=index_c+1
          df1=rbind(df1,temp)
        } else if (a=="y"){
          list_check_movie=rbind(list_check_movie,list[i])
          list_check_index=rbind(list_check_index,j)
        } else if (a=="v"){
          lineage=unique(temp$lineage)
          number=unique(temp$number_in_lineage)
          if(number %% 2 == 0){
            mother=number/2
          } else {
            mother=(number-1)/2
          }
          if(mother>1 && (mother %%2 ==0)){
            grandm=mother/2
          } else {
            grandm=(mother-1)/2
          }
          daughters=c(number*2,number*2+1)
          family=c(grandm,mother,daughters,(daughters[1]*2),(daughters[1]*2+1),(daughters[2]*2),(daughters[2]*2+1))
          temp2<-temp
          for (f in family){
            if(f %in% d$number_in_lineage[which(d$lineage==lineage)]){
              if(f> number){
                temp2=rbind(temp2,d[which(d$lineage==lineage & d$number_in_lineage==f)[c(1:min(c(20,length(d$number_in_lineage[which(d$lineage==lineage & d$number_in_lineage==f)]))))],])
              } else {
                temp2=rbind(temp2,d[which(d$lineage==lineage & d$number_in_lineage==f),])
              }
            }
          }
          temp2$number_in_lineage[which(temp2$number_in_lineage==number)]=0
          y=temp2$Cell_Intensisty
          x=temp2$Frame
          z=temp2$number_in_lineage
          plot.new()
          plot(x,y,col=factor(z))
          rm(daughters,mother,grandm,family)
          a2=readline("c=keep and cut, n=reject")
          if(a2=="c"){
            out <- sapply(list(x,y),"[",identify(x,y))
            if(nrow(out) %%2 ==0){
              for (f in seq(1,nrow(out),2)){
                temp2=temp[c(which(temp$Frame==out[f,1]):which(temp$Frame==out[(f+1),1])),]
                temp2$index_c=index_u+1
                index_u=index_u+1
                temp2$event=readline("s=start, m=mitosis, p=just a portion, c=complete")
                df2=rbind(df2,temp2)
              }
            } else {
              list_redo=rbind(list_redo,d[which(d$CellIdInAutoTracking==(j)),])
            }
          }
        }
        tempj=j
      }
    }
    
  }
}


# correct if there was a "p" ####
x11()
{
for (j in unique(list_redo$CellIdInAutoTracking)){
  # for (j in unique(d$CellIdInAutoTracking)[c(13:length(unique(d$CellIdInAutoTracking)))]){
  if((j %in% unique(df1$CellIdInAutoTracking)==TRUE)){
    df1<-df1[-which(df1$CellIdInAutoTracking==j),]
  } else if((j %in% unique(df2$CellIdInAutoTracking)==TRUE)){
    df2<-df2[-which(df2$CellIdInAutoTracking==j),]
  }
  temp=d[which(d$CellIdInAutoTracking==j),]
  if(nrow(temp)>40){
    y=temp$Cell_Intensisty
    x=temp$Frame
    plot.new()
    plot(x,y)
    # colnames(temp)[c(10:13)]=c("sum_im-median","sum_mask-median","median","mask")
    a=readline("keep, complete cellcycle = c, v=check family, check raw image=y, reject = n, made a mistake with the previous plot = p")
    if(a=="p") {
      list_redo=rbind(list_redo,d[which(d$CellIdInAutoTracking==(j-1)),])
      a=readline(" now, for this current plot: keep, complete cc = c, v=check family, check raw image=y, reject = n")
    }
    if(a=="c"){
      temp$index_c=index_c+1
      index_c=index_c+1
      df1=rbind(df1,temp)
    } else if (a=="y"){
      list_check_movie=rbind(list_check_movie,list[i])
      list_check_index=rbind(list_check_index,j)
    } else if (a=="v"){
      lineage=unique(temp$lineage)
      number=unique(temp$number_in_lineage)
      if(number %% 2 == 0){
        mother=number/2
      } else {
        mother=(number-1)/2
      }
      if(mother>1 && (mother %%2 ==0)){
        grandm=mother/2
      } else {
        grandm=(mother-1)/2
      }
      daughters=c(number*2,number*2+1)
      family=c(grandm,mother,daughters,(daughters[1]*2),(daughters[1]*2+1),(daughters[2]*2),(daughters[2]*2+1))
      temp2<-temp
      for (f in family){
        if(f %in% d$number_in_lineage[which(d$lineage==lineage)]){
          if(f> number){
            temp2=rbind(temp2,d[which(d$lineage==lineage & d$number_in_lineage==f)[c(1:min(c(20,length(d$number_in_lineage[which(d$lineage==lineage & d$number_in_lineage==f)]))))],])
          } else {
            temp2=rbind(temp2,d[which(d$lineage==lineage & d$number_in_lineage==f),])
          }
        }
      }
      temp2$number_in_lineage[which(temp2$number_in_lineage==number)]=0
      y=temp2$Cell_Intensisty
      x=temp2$Frame
      z=temp2$number_in_lineage
      plot.new()
      plot(x,y,col=factor(z))
      rm(daughters,mother,grandm,family)
      a2=readline("c=keep and cut, n=reject")
      if(a2=="c"){
        out <- sapply(list(x,y),"[",identify(x,y))
        if(nrow(out) %%2 ==0){
          for (f in seq(1,nrow(out),2)){
            temp2=temp[c(which(temp$Frame==out[f,1]):which(temp$Frame==out[(f+1),1])),]
            temp2$index_c=index_u+1
            index_u=index_u+1
            temp2$event=readline("s=start, m=mitosis, p=just a portion, c=complete")
            df2=rbind(df2,temp2)
          }
        } else {
          list_redo=rbind(list_redo,d[which(d$CellIdInAutoTracking==(j-1)),])
        }
      }
    }
  }
}
}


# Save at any point when you stopped in the process ####
dirtemp=paste0(dirout,str_split(list[i],".csv")[[1]][1],"/")
dir.create(dirtemp)

temp<-df2[which(df2$event=="c"),]
# index_c=max(df1$index_c)
if(nrow(temp)>0){
  for (j in unique(temp$index_c)){
    temp$index_c[which(temp$index_c==j)]=index_c+1
    index_c=index_c+1
  }
  temp<-temp[,-which(colnames(temp)=="event")]
  df1=rbind(df1,temp)
  df2<-df2[-which(df2$event=="c"),]  
} else {
}

list_save=list(df1,df2,index_c,index_u,list_check_index,list_check_movie)
listnames=c("df1","df2","index_c","index_u","list_check_index","list_check_movie")
for (j in c(1:length(list_save))){
write.csv(list_save[[j]],paste0(dirtemp,str_split(list[i],".csv")[[1]][1],"_",listnames[j],".csv"))  
}


# small check ####
x11()
{
for (i in unique(df2$index_c)){
  plot.new()
  plot(y=df2$Cell_Intensisty[which(df2$index_c==i)],x=df2$Frame[which(df2$index_c==i)])
  a=readline("c")
}
}



# build dataframe with everyone ####
list=list.dirs(dirout,full.names = FALSE)
df1=NULL
df2=NULL
checkd1=NULL
checkd2=NULL
for (i in c(2:(length(list)-1))){
# for (i in listfiles){
  if(TRUE %in% grepl("df1",list.files(paste0(dirout,list[i],"/")))){
  d1=read.csv(paste0(dirout,list[i],"/",list.files(paste0(dirout,list[i],"/"))[which(grepl("df1",list.files(paste0(dirout,list[i],"/")))==TRUE)]))
  } else {
    d1=NULL
  }
  d2=read.csv(paste0(dirout,list[i],"/",list.files(paste0(dirout,list[i],"/"))[which(grepl("df2",list.files(paste0(dirout,list[i],"/")))==TRUE)]))
  if(is.data.frame(d1)==TRUE && nrow(d1)>0){
    df1=rbind(df1,d1)
    checkd1=c(checkd1,min(d1$index_c),max(d1$index_c))
  }
  if(nrow(d2)>0){
    df2=rbind(df2,d2)
    checkd2=c(checkd2,min(d2$index_c),max(d2$index_c))
  }
}
write.csv(df1,paste0(dirout,"df1_complete-tracks.csv"),row.names = FALSE)
write.csv(df2,paste0(dirout,"df2_uncomplete-tracks.csv"),row.names = FALSE)





#################### 2) ADD the FUCCI information ####
#for complete trajectories : ####
d=read.csv("E:/Analysis/00_GROWTH-CURVES/00_SINGLE CURVES/00-DATA/151028-traj_HeLa-control_170102/results_df/df1_complete-tracks.csv",stringsAsFactors = FALSE)
mat=matrix(nrow=0,ncol=5)
mat=as.data.frame(mat)
colnames(mat)=c("ID","index_c","tbirth","tg1","tm")
info_stage<-mat
rm(mat)
df=NULL
list_redo=NULL
list_check=NULL
x11()
{
for (i in unique(d$index_c)){#[c(7:length(unique(d$index_c)))]
  temp=d[which(d$index_c==i),]
  plot.new()
  plot(y=temp$sum_im.median,x=temp$Frame)
  a=readline("check fucci signla: c=keep and cut, n=exclude, v=NA, y=check the raw movie, p= made a mistake with the previous plot")
  if(a=="p"){
    list_redo=rbind(list_redo,d[which(d$index_c==unique(d$index_c)[(i-1)])])
    a=readline(" now for this plot check fucci signla: c=keep and cut, n=exclude, v=NA, p= made a mistake with the previous plot")
  }
  
  if(a=="c"){
    print("enter G1/S transition")
    out <- sapply(list(temp$Frame,temp$sum_im.median),"[",identify(temp$Frame,temp$sum_im.median))
    if(length(out)==2){
      tg1=out[1]
      plot.new()
      plot(y=temp$Cell_Intensisty,x=temp$Frame)
      abline(v=tg1)
      a2=readline("plot ok? c=keep, n=exclude,p=redo")
    } else {
      list_redo=rbind(list_redo,temp)
    }
  } else if (a=="v"){
    tg1=NA
    plot.new()
    plot(y=temp$Cell_Intensisty,x=temp$Frame)
    a2=readline("plot ok? c=keep, n=exclude,p=redo,m=do some cuts")
  } else if(a=="y"){
    list_check=rbind(list_check,temp)
  }
  if(exists("a2")==TRUE && a2=="c"){
    info_stage[nrow(info_stage)+1,]=NA
    info_stage[nrow(info_stage),c(3,5)]=c(temp$Frame[1],temp$Frame[nrow(temp)])
    info_stage[nrow(info_stage),4]=tg1
    info_stage[nrow(info_stage),2]=i
    info_stage[nrow(info_stage),1]=temp$ID[1]
    rm(a2)
  } else if(exists("a2")==TRUE && a2=="m"){
    print("enter birth and mitosis")
    out <- sapply(list(temp$Frame,temp$Cell_Intensisty),"[",identify(temp$Frame,temp$Cell_Intensisty))
    info_stage[nrow(info_stage)+1,]=NA
    info_stage[nrow(info_stage),c(3,5)]=c(out[1,1],out[2,1])
    info_stage[nrow(info_stage),4]=tg1
    info_stage[nrow(info_stage),2]=i
    info_stage[nrow(info_stage),1]=temp$ID[1]
    rm(a2)
  } else if(exists("a2")==TRUE && a2=="p"){
    list_redo=rbind(list_redo,temp)
    rm(a2)
  }
  df=rbind(df,temp)
}
}
write.csv(info_stage,paste0(dirout,"df1_info_stage.csv"),row.names = FALSE)


# for uncomplete trajectories ####
d=read.csv("E:/Analysis/00_GROWTH-CURVES/00_SINGLE CURVES/00-DATA/151028-traj_HeLa-rosco_170201/results_df/df2_uncomplete-tracks.csv",stringsAsFactors = FALSE)
mat=matrix(nrow=0,ncol=6)
mat=as.data.frame(mat)
colnames(mat)=c("ID","index_c","tbirth","tg1","tm","unidentified")
info_stage<-mat
rm(mat)
df=NULL
list_redo=NULL
list_check=NULL
x11()
{
  for (i in unique(d$index_c)){#[c(7:length(unique(d$index_c)))]
    temp=d[which(d$index_c==i),]
    if(nrow(temp)>40){
    plot.new()
    plot(y=temp$sum_im.median,x=temp$Frame,main=temp$event[1])
    a=readline("check fucci signla: c=keep and cut, n=exclude, v=NA, y=check the raw movie, p= made a mistake with the previous plot")
    if(a=="p"){
      list_redo=rbind(list_redo,d[which(d$index_c==unique(d$index_c)[(i-1)])])
      a=readline(" now for this plot check fucci signla: c=keep and cut, n=exclude, v=NA, p= made a mistake with the previous plot")
    }
    
    if(a=="c"){
      print("enter G1/S transition")
      out <- sapply(list(temp$Frame,temp$sum_im.median),"[",identify(temp$Frame,temp$sum_im.median))
      if(length(out)==2){
        tg1=out[1]
        plot.new()
        plot(y=temp$Cell_Intensisty,x=temp$Frame)
        abline(v=tg1)
        a2=readline("plot ok? c=keep, n=exclude,p=redo")
      } else {
        list_redo=rbind(list_redo,temp)
      }
    } else if (a=="v"){
      tg1=NA
      plot.new()
      plot(y=temp$Cell_Intensisty,x=temp$Frame)
      a2=readline("plot ok? c=keep, n=exclude,p=redo")
    } else if(a=="y"){
      list_check=rbind(list_check,temp)
    }
    if(exists("a2")==TRUE && a2=="c"){
      # print("enter birth and mitosis")
      # out <- sapply(list(temp$Frame,temp$Cell_Intensisty),"[",identify(temp$Frame,temp$Cell_Intensisty))
      info_stage[(nrow(info_stage)+1),]=NA
      if(temp$event[1]=="s"){
        info_stage[nrow(info_stage),3]=temp$Frame[1]
      } else if(temp$event[1]=="m"){
        info_stage[nrow(info_stage),5]=max(temp$Frame)
      } else if (temp$event[1]=="p"){
        info_stage[nrow(info_stage),6]="p"
      }
      info_stage[nrow(info_stage),4]=tg1
      info_stage[nrow(info_stage),2]=i
      info_stage[nrow(info_stage),1]=temp$ID[1]
      rm(a2)
    } else if(exists("a2")==TRUE && a2=="p"){
      list_redo=rbind(list_redo,temp)
      rm(a2)
    }
    df=rbind(df,temp)
    }
  }
}
write.csv(info_stage,paste0(dirout,"df2_info_stage.csv"),row.names = FALSE)







################### 3) ADD the transition point in early G1 ####
df<-read.csv("E:/Analysis/00_GROWTH-CURVES/00_SINGLE CURVES/170201_HeLa-single-curves-analysis/151028_rosco/151028-rosco-uncomplete-v2-QRfiltering_calcul-gr_dtav3_wo-ol.csv")
table<-read.csv("E:/Analysis/00_GROWTH-CURVES/00_SINGLE CURVES/170201_HeLa-single-curves-analysis/151028_rosco/151028-rosco-uncomplete-v2-QRfiltering_info_stage_keypoints_dtav3.csv")

table$teg1=NA
table$veg1=NA
x11()
{
  for (i in unique(table$index_c)){
    if(is.na(table$tg1[which(table$index_c==i)])==FALSE & (i%in% unique(df$index_c)==TRUE)){
      temp=df[which(df$index_c==i),]
      plot.new()
      plot(y=temp$av,x=temp$Frame)
      abline(v=table$tg1[which(table$index_c==i)])
      a<-readline("transition? c=yes,n=no")
      if(a=="c"){
        out <- sapply(list(temp$Frame,temp$av),"[",identify(temp$Frame,temp$av))
        table$veg1[which(table$index_c==i)]=out[2]
        table$teg1[which(table$index_c==i)]=out[1]
      } else {}
      }
    }
  }
dir.create("E:/Analysis/00_GROWTH-CURVES/00_SINGLE CURVES/170222_new-tests-for-figures/")
dirsave="E:/Analysis/00_GROWTH-CURVES/00_SINGLE CURVES/170222_new-tests-for-figures/"
write.csv(table,"E:/Analysis/00_GROWTH-CURVES/00_SINGLE CURVES/170222_new-tests-for-figures/151028-rosco-uncomplete-v2-QRfiltering_info_stage_keypoints_dtav3_earlyg1.csv")
