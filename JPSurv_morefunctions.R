#### More functions requested for JPSurv web app new features (discussed in meeting 06/05/2019)


#########################################################################################################
# download.data: a function that returns a merged data set for download including 
# - the selected cohort input data
# - the accompanying data for plot.surv.year/plot.dying.year
# arguments:
# input - input dataset read in by function joinpoint.seerdata.
# fit - joinpoint object containing the model output.
# nJP - the number of joinpoints in the model.
# yearvar - the variable name for year of diagnosis used in argument 'year' of the function joinpoint.
# type - two options for this argument: "graph" for graph data and "full" for full data.
# subset - an optional string specifying a subset of observations used in the fitting process.
# interval - the variable name for year since diagnosis. The default is 'Interval'.
# int.col - the interval values selected for the plot if the downloadtype="graph". The default is NULL.
#########################################################################################################

download.data<-function(input,fit,nJP,yearvar,downloadtype,subset=NULL,interval="Interval",int.col=NULL){
  input.sub<-subset(input,eval(parse(text=subset)))
  if(downloadtype=="graph"){
    output.sub<-fit$FitList[[nJP+1]]$predicted
    if(!is.null(int.col)){
      output.sub<-output.sub[which(output.sub[,interval] %in% int.col),]
    }
  }else if(downloadtype=="full"){
    output.sub<-fit$FitList[[nJP+1]]$fullpredicted
  }else{
    print("There is no such type of data for download. Please define the type using either 'graph' or 'full'.")
    break
  }
  yearint.input<-paste(input.sub[,yearvar],input.sub[,interval],sep="_")
  yearint.output<-paste(output.sub[,yearvar],output.sub[,interval],sep="_")
  rows.match<-match(yearint.output,yearint.input)
  cols.input<-colnames(input.sub)
  del.cols<-c("Page_type","Observed_SE_Interval","Observed_SE_Cum")
  cols.input<-cols.input[-which(cols.input %in% del.cols)]
  cols.output<-colnames(output.sub)
  merge.sub<-output.sub
  merge.sub[,cols.input]<-input.sub[rows.match,cols.input]
  colnames(merge.sub)[which(colnames(merge.sub)=="pred_int")]<-"Predicted_Survival_Int"
  colnames(merge.sub)[which(colnames(merge.sub)=="pred_cum")]<-"Predicted_Survival_Cum"
  colnames(merge.sub)[which(colnames(merge.sub)=="pred_int_se")]<-"Predicted_Survival_Int_SE"
  colnames(merge.sub)[which(colnames(merge.sub)=="pred_cum_se")]<-"Predicted_Survival_Cum_SE"
  merge.sub[,"Predicted_ProbDeath_Int"]<-1-merge.sub$Predicted_Survival_Int
  merge.sub[,"Predicted_ProbDeath_Int_SE"]<-merge.sub$Predicted_Survival_Int_SE
  cols.pred<-c("Predicted_Survival_Int","Predicted_ProbDeath_Int","Predicted_Survival_Cum",
              "Predicted_Survival_Int_SE","Predicted_Survival_Cum_SE","Predicted_ProbDeath_Int_SE")
  merge.data<-merge.sub[,c(cols.input,cols.pred)]
  if(downloadtype=="full"){
    merge.data[,yearvar]<-output.sub[,yearvar]
    merge.data[,interval]<-output.sub[,interval]
    subset.list<-strsplit(subset, " & ")
    subset.list[[1]]<-subset.list[[1]][which(grepl("==", subset.list[[1]])==T)]
    cohort.list<-strsplit(subset.list[[1]], "==")
    cohort.vars<-sapply(cohort.list, "[", 1)
    cohort.vars<-gsub(" ","",cohort.vars)
    cohort.values<-sapply(cohort.list, "[", 2)
    cohort.values<-gsub(" ","",cohort.values)
    cohort.vars<-cohort.vars[which(!cohort.vars %in% c(yearvar,interval))]
    cohort.values<-cohort.values[which(!cohort.vars %in% c(yearvar,interval))]
    merge.data[,cohort.vars]<-t(replicate(dim(merge.data)[1],cohort.values))
  }
  merge.data[,yearvar]<-as.numeric(merge.data[,yearvar])
  return(merge.data)
}

#########################################################################################################
# plot.dying.year.annotate: a function that returns a plot for Percent Change in the Annual Probability of
# Dying of Cancer by Diagnosis year using ggplot
# The annotation feature is available for nJP<=3 and the number of multiple intervals selected <=3
# arguments:
# plotdata - the graph data returned by function download.data with downloadtype="graph".
# fit - joinpoint object containing the model output.
# nJP - the number of joinpoints in the model.
# yearvar - the variable name for year of diagnosis used in argument 'year' of the function joinpoint.
# obsintvar - the variable name for observed interval survival. The default is "Relative_Survival_Interval"
#             for relative survival data. For cause-specific data, it needs to be changed accordingly.
# predintvar - the variable name for predicted interval survival. The default is "Predicted_ProbDeath_Int".
# interval - the variable name for year since diagnosis. The default is 'Interval'.
# annotation - the indicator for the annotation feature. The default is 0 (no annotation on the plot).
# topanno - the indicator for showing the top curve annotation. The default is 1 (annotation for the top curve).
#########################################################################################################

### define a plot function for Percent Change in the Annual Probability of Dying of Cancer by Diagnosis year with annotations
plot.dying.year.annotate<-function(plotdata,fit,nJP,yearvar,obsintvar="Relative_Survival_Interval",predintvar="Predicted_ProbDeath_Int",interval="Interval",annotation=0,topanno=1){
  title.rch<-"Percent Change in the Annual Probability of Dying of Cancer \n by Diagnosis Year"
  interval.values<-as.numeric(unique(plotdata[,"Interval"]))
  interval.labels<-paste((interval.values-1)," to ",interval.values," years since diagnosis",sep="")
  if(interval.values[1]==1){
    interval.labels[1]<-"0 to 1 year since diagnosis"
  }
  
  if(annotation==1){  
    if(length(interval.values)<=3 & nJP<=3){
      ### haz results are the same for all interval values
      jp.loc<-fit$FitList[[nJP+1]]$jp
      haz.apc<-aapc(fit$FitList[[nJP+1]], type="RelChgHaz", interval=interval.values[1])
      annot.strs<-paste(sprintf("%.1f",100*haz.apc$estimate),"%",sep="")
      x.values<-list()
      y.values<-list()
      for(i in 1:length(interval.values)){
        int.i<-interval.values[i]
        plotdata.i<-plotdata[which(plotdata[,interval]==int.i),]
        if(max(plotdata.i[,yearvar])==max(haz.apc$end.year) | nJP==0){
          end.year<-haz.apc$end.year
          end.year[length(haz.apc$end.year)]<-max(plotdata.i[,yearvar])
          y.values.i<-1/2*(plotdata.i[which(plotdata.i[,yearvar] %in% haz.apc$start.year),predintvar]+
                             plotdata.i[which(plotdata.i[,yearvar] %in% end.year),predintvar])
          x.values.i<-1/2*(haz.apc$start.year+end.year)
        }
        if(max(plotdata.i[,yearvar])<max(haz.apc$end.year)){
          if(nJP==1){
            if(max(plotdata.i[,yearvar])>jp.loc){
              end.year<-haz.apc$end.year
              end.year[length(haz.apc$end.year)]<-max(plotdata.i[,yearvar])
              y.values.i<-1/2*(plotdata.i[which(plotdata.i[,yearvar] %in% haz.apc$start.year),predintvar]+
                               plotdata.i[which(plotdata.i[,yearvar] %in% end.year),predintvar])
              x.values.i<-1/2*(haz.apc$start.year+end.year)
            }else{
              haz.apc<-haz.apc[1:nJP,]
              y.values.i<-1/2*(plotdata.i[which(plotdata.i[,yearvar] %in% haz.apc$start.year),predintvar]+
                               plotdata.i[which(plotdata.i[,yearvar] %in% haz.apc$end.year),predintvar])
              x.values.i<-1/2*(haz.apc$start.year+haz.apc$end.year)
            }
          } ## nJP=1 end
          if(nJP==2){
            if(max(plotdata.i[,yearvar])>max(jp.loc)){
              end.year<-haz.apc$end.year
              end.year[length(haz.apc$end.year)]<-max(plotdata.i[,yearvar])
              y.values.i<-1/2*(plotdata.i[which(plotdata.i[,yearvar] %in% haz.apc$start.year),predintvar]+
                               plotdata.i[which(plotdata.i[,yearvar] %in% end.year),predintvar])
              x.values.i<-1/2*(haz.apc$start.year+end.year)
            } else if(max(plotdata.i[,yearvar])<=max(jp.loc) & max(plotdata.i[,yearvar])>min(jp.loc)){
              haz.apc<-haz.apc[1:nJP,]
              y.values.i<-1/2*(plotdata.i[which(plotdata.i[,yearvar] %in% haz.apc$start.year),predintvar]+
                               plotdata.i[which(plotdata.i[,yearvar] %in% haz.apc$end.year),predintvar])
              x.values.i<-1/2*(haz.apc$start.year+haz.apc$end.year)
            } else{
              haz.apc<-haz.apc[1:(nJP-1),]
              end.year<-haz.apc$end.year
              end.year[length(haz.apc$end.year)]<-max(plotdata.i[,yearvar])
              y.values.i<-1/2*(plotdata.i[which(plotdata.i[,yearvar] %in% haz.apc$start.year),predintvar]+
                               plotdata.i[which(plotdata.i[,yearvar] %in% end.year),predintvar])
              x.values.i<-1/2*(haz.apc$start.year+end.year)
            }
          } ## nJP=2 end
          if(nJP==3){
            if(max(plotdata.i[,yearvar])>max(jp.loc)){
              end.year<-haz.apc$end.year
              end.year[length(haz.apc$end.year)]<-max(plotdata.i[,yearvar])
              y.values.i<-1/2*(plotdata.i[which(plotdata.i[,yearvar] %in% haz.apc$start.year),predintvar]+
                               plotdata.i[which(plotdata.i[,yearvar] %in% end.year),predintvar])
              x.values.i<-1/2*(haz.apc$start.year+end.year)
            } else if(max(plotdata.i[,yearvar])<=max(jp.loc) & max(plotdata.i[,yearvar])>jp.loc[2]){
              haz.apc<-haz.apc[1:nJP,]
              end.year<-haz.apc$end.year
              end.year[length(haz.apc$end.year)]<-max(plotdata.i[,yearvar])
              y.values.i<-1/2*(plotdata.i[which(plotdata.i[,yearvar] %in% haz.apc$start.year),predintvar]+
                               plotdata.i[which(plotdata.i[,yearvar] %in% end.year),predintvar])
              x.values.i<-1/2*(haz.apc$start.year+end.year)
            } else if(max(plotdata.i[,yearvar])<=jp.loc[2] & max(plotdata.i[,yearvar])>=jp.loc[1]){
              haz.apc<-haz.apc[1:(nJP-1),]
              end.year<-haz.apc$end.year
              end.year[length(haz.apc$end.year)]<-max(plotdata.i[,yearvar])
              y.values.i<-1/2*(plotdata.i[which(plotdata.i[,yearvar] %in% haz.apc$start.year),predintvar]+
                               plotdata.i[which(plotdata.i[,yearvar] %in% end.year),predintvar])
              x.values.i<-1/2*(haz.apc$start.year+end.year)
            } else {
              haz.apc<-haz.apc[1:(nJP-2),]
              end.year<-haz.apc$end.year
              end.year[length(haz.apc$end.year)]<-max(plotdata.i[,yearvar])
              y.values.i<-1/2*(plotdata.i[which(plotdata.i[,yearvar] %in% haz.apc$start.year),predintvar]+
                               plotdata.i[which(plotdata.i[,yearvar] %in% end.year),predintvar])
              x.values.i<-1/2*(haz.apc$start.year+end.year)
            }
          } ## nJP=3 end
        } 
        x.values[[i]]<-x.values.i
        y.values[[i]]<-y.values.i
      }## interval.values iteration end
      if(length(interval.values)==1){
        hues=seq(15,375,length=length(interval.values)+1)
        gg_color_hue<-hcl(h=hues,l=65,c=100)[1:length(interval.values)]
        plot<-ggplot(data=plotdata, aes(x=plotdata[,yearvar], y=plotdata[,predintvar],group=as.factor(plotdata[,interval]), colour=as.factor(plotdata[,interval]))) + 
          geom_line(data=plotdata, aes(x=plotdata[,yearvar], y=plotdata[,predintvar],group=as.factor(plotdata[,interval]), colour=as.factor(plotdata[,interval])),linetype="solid")+
          geom_point(data=plotdata, aes(x=plotdata[,yearvar], y=1-plotdata[,obsintvar], group=as.factor(plotdata[,interval]), colour=as.factor(plotdata[,interval])))+
          xlab('Year at diagnosis') + ylab('Annual Probability of Cancer Death (%)')+
          theme_bw()+
          ggtitle(title.rch) +
          coord_cartesian(ylim=c(0,1))+
          scale_y_continuous(breaks=seq(0,1,0.1),labels = scales::percent)+
          scale_x_continuous(breaks=seq(min(plotdata[,yearvar],na.rm=T),max(plotdata[,yearvar],na.rm=T),5))+
          scale_color_hue(labels = interval.labels)+
          theme(legend.position="bottom")+
          theme(legend.title=element_blank())+
          theme(plot.title = element_text(hjust = 0.5))+
          annotate("text", x = x.values[[1]][1], y = y.values[[1]][1]+0.03, label = annot.strs[1], colour=gg_color_hue[1])+
          annotate("text", x = x.values[[1]][2], y = y.values[[1]][2]+0.03, label = annot.strs[2], colour=gg_color_hue[1])+
          annotate("text", x = x.values[[1]][3], y = y.values[[1]][3]+0.03, label = annot.strs[3], colour=gg_color_hue[1])+
          annotate("text", x = x.values[[1]][4], y = y.values[[1]][4]+0.03, label = annot.strs[4], colour=gg_color_hue[1])
      }
      if(length(interval.values)==2){
        hues=seq(15,375,length=length(interval.values)+1)
        gg_color_hue<-hcl(h=hues,l=65,c=100)[1:length(interval.values)]
        if(topanno==1){
          ymeans.byint<-aggregate(plotdata[,"Predicted_ProbDeath_Int"], list(plotdata[,interval]), mean)
          topint<-which(ymeans.byint[,2]==max(ymeans.byint[,2]))
          plot<-ggplot(data=plotdata, aes(x=plotdata[,yearvar], y=plotdata[,predintvar],group=as.factor(plotdata[,interval]), colour=as.factor(plotdata[,interval]))) + 
            geom_line(data=plotdata, aes(x=plotdata[,yearvar], y=plotdata[,predintvar],group=as.factor(plotdata[,interval]), colour=as.factor(plotdata[,interval])),linetype="solid")+
            geom_point(data=plotdata, aes(x=plotdata[,yearvar], y=1-plotdata[,obsintvar], group=as.factor(plotdata[,interval]), colour=as.factor(plotdata[,interval])))+
            xlab('Year at diagnosis') + ylab('Annual Probability of Cancer Death (%)')+
            theme_bw()+
            ggtitle(title.rch) +
            coord_cartesian(ylim=c(0,1))+
            scale_y_continuous(breaks=seq(0,1,0.1),labels = scales::percent)+
            scale_x_continuous(breaks=seq(min(plotdata[,yearvar],na.rm=T),max(plotdata[,yearvar],na.rm=T),5))+
            scale_color_hue(labels = interval.labels)+
            theme(legend.position="bottom")+
            theme(legend.title=element_blank())+
            theme(plot.title = element_text(hjust = 0.5))+
            annotate("text", x = x.values[[topint]][1], y = y.values[[topint]][1]+0.03, label = annot.strs[1], colour=gg_color_hue[topint])+
            annotate("text", x = x.values[[topint]][2], y = y.values[[topint]][2]+0.03, label = annot.strs[2], colour=gg_color_hue[topint])+
            annotate("text", x = x.values[[topint]][3], y = y.values[[topint]][3]+0.03, label = annot.strs[3], colour=gg_color_hue[topint])+
            annotate("text", x = x.values[[topint]][4], y = y.values[[topint]][4]+0.03, label = annot.strs[4], colour=gg_color_hue[topint])
        }else{
          plot<-ggplot(data=plotdata, aes(x=plotdata[,yearvar], y=plotdata[,predintvar],group=as.factor(plotdata[,interval]), colour=as.factor(plotdata[,interval]))) + 
            geom_line(data=plotdata, aes(x=plotdata[,yearvar], y=plotdata[,predintvar],group=as.factor(plotdata[,interval]), colour=as.factor(plotdata[,interval])),linetype="solid")+
            geom_point(data=plotdata, aes(x=plotdata[,yearvar], y=1-plotdata[,obsintvar], group=as.factor(plotdata[,interval]), colour=as.factor(plotdata[,interval])))+
            xlab('Year at diagnosis') + ylab('Annual Probability of Cancer Death (%)')+
            theme_bw()+
            ggtitle(title.rch) +
            coord_cartesian(ylim=c(0,1))+
            scale_y_continuous(breaks=seq(0,1,0.1),labels = scales::percent)+
            scale_x_continuous(breaks=seq(min(plotdata[,yearvar],na.rm=T),max(plotdata[,yearvar],na.rm=T),5))+
            scale_color_hue(labels = interval.labels)+
            theme(legend.position="bottom")+
            theme(legend.title=element_blank())+
            theme(plot.title = element_text(hjust = 0.5))+
            annotate("text", x = x.values[[1]][1], y = y.values[[1]][1]+0.03, label = annot.strs[1], colour=gg_color_hue[1])+
            annotate("text", x = x.values[[1]][2], y = y.values[[1]][2]+0.03, label = annot.strs[2], colour=gg_color_hue[1])+
            annotate("text", x = x.values[[1]][3], y = y.values[[1]][3]+0.03, label = annot.strs[3], colour=gg_color_hue[1])+
            annotate("text", x = x.values[[1]][4], y = y.values[[1]][4]+0.03, label = annot.strs[4], colour=gg_color_hue[1])+
            annotate("text", x = x.values[[2]][1], y = y.values[[2]][1]+0.03, label = annot.strs[1], colour=gg_color_hue[2])+
            annotate("text", x = x.values[[2]][2], y = y.values[[2]][2]+0.03, label = annot.strs[2], colour=gg_color_hue[2])+
            annotate("text", x = x.values[[2]][3], y = y.values[[2]][3]+0.03, label = annot.strs[3], colour=gg_color_hue[2])+
            annotate("text", x = x.values[[2]][4], y = y.values[[2]][4]+0.03, label = annot.strs[4], colour=gg_color_hue[2])
        }
      }
      if(length(interval.values)==3){
        hues=seq(15,375,length=length(interval.values)+1)
        gg_color_hue<-hcl(h=hues,l=65,c=100)[1:length(interval.values)]
        if(topanno==1){
          ymeans.byint<-aggregate(plotdata[,"Predicted_ProbDeath_Int"], list(plotdata[,interval]), mean)
          topint<-which(ymeans.byint[,2]==max(ymeans.byint[,2]))
          plot<-ggplot(data=plotdata, aes(x=plotdata[,yearvar], y=plotdata[,predintvar],group=as.factor(plotdata[,interval]), colour=as.factor(plotdata[,interval]))) + 
            geom_line(data=plotdata, aes(x=plotdata[,yearvar], y=plotdata[,predintvar],group=as.factor(plotdata[,interval]), colour=as.factor(plotdata[,interval])),linetype="solid")+
            geom_point(data=plotdata, aes(x=plotdata[,yearvar], y=1-plotdata[,obsintvar], group=as.factor(plotdata[,interval]), colour=as.factor(plotdata[,interval])))+
            xlab('Year at diagnosis') + ylab('Annual Probability of Cancer Death (%)')+
            theme_bw()+
            ggtitle(title.rch) +
            coord_cartesian(ylim=c(0,1))+
            scale_y_continuous(breaks=seq(0,1,0.1),labels = scales::percent)+
            scale_x_continuous(breaks=seq(min(plotdata[,yearvar],na.rm=T),max(plotdata[,yearvar],na.rm=T),5))+
            scale_color_hue(labels = interval.labels)+
            theme(legend.position="bottom")+
            theme(legend.title=element_blank())+
            theme(plot.title = element_text(hjust = 0.5))+
            annotate("text", x = x.values[[topint]][1], y = y.values[[topint]][1]+0.03, label = annot.strs[1], colour=gg_color_hue[topint])+
            annotate("text", x = x.values[[topint]][2], y = y.values[[topint]][2]+0.03, label = annot.strs[2], colour=gg_color_hue[topint])+
            annotate("text", x = x.values[[topint]][3], y = y.values[[topint]][3]+0.03, label = annot.strs[3], colour=gg_color_hue[topint])+
            annotate("text", x = x.values[[topint]][4], y = y.values[[topint]][4]+0.03, label = annot.strs[4], colour=gg_color_hue[topint])
        }else{
          plot<-ggplot(data=plotdata, aes(x=plotdata[,yearvar], y=plotdata[,predintvar],group=as.factor(plotdata[,interval]), colour=as.factor(plotdata[,interval]))) + 
            geom_line(data=plotdata, aes(x=plotdata[,yearvar], y=plotdata[,predintvar],group=as.factor(plotdata[,interval]), colour=as.factor(plotdata[,interval])),linetype="solid")+
            geom_point(data=plotdata, aes(x=plotdata[,yearvar], y=1-plotdata[,obsintvar], group=as.factor(plotdata[,interval]), colour=as.factor(plotdata[,interval])))+
            xlab('Year at diagnosis') + ylab('Annual Probability of Cancer Death (%)')+
            theme_bw()+
            ggtitle(title.rch) +
            coord_cartesian(ylim=c(0,1))+
            scale_y_continuous(breaks=seq(0,1,0.1),labels = scales::percent)+
            scale_x_continuous(breaks=seq(min(plotdata[,yearvar],na.rm=T),max(plotdata[,yearvar],na.rm=T),5))+
            scale_color_hue(labels = interval.labels)+
            theme(legend.position="bottom")+
            theme(legend.title=element_blank())+
            theme(plot.title = element_text(hjust = 0.5))+
            annotate("text", x = x.values[[1]][1], y = y.values[[1]][1]+0.03, label = annot.strs[1], colour=gg_color_hue[1])+
            annotate("text", x = x.values[[1]][2], y = y.values[[1]][2]+0.03, label = annot.strs[2], colour=gg_color_hue[1])+
            annotate("text", x = x.values[[1]][3], y = y.values[[1]][3]+0.03, label = annot.strs[3], colour=gg_color_hue[1])+
            annotate("text", x = x.values[[1]][4], y = y.values[[1]][4]+0.03, label = annot.strs[4], colour=gg_color_hue[1])+
            annotate("text", x = x.values[[2]][1], y = y.values[[2]][1]+0.03, label = annot.strs[1], colour=gg_color_hue[2])+
            annotate("text", x = x.values[[2]][2], y = y.values[[2]][2]+0.03, label = annot.strs[2], colour=gg_color_hue[2])+
            annotate("text", x = x.values[[2]][3], y = y.values[[2]][3]+0.03, label = annot.strs[3], colour=gg_color_hue[2])+
            annotate("text", x = x.values[[2]][4], y = y.values[[2]][4]+0.03, label = annot.strs[4], colour=gg_color_hue[2])+
            annotate("text", x = x.values[[3]][1], y = y.values[[3]][1]+0.03, label = annot.strs[1], colour=gg_color_hue[3])+
            annotate("text", x = x.values[[3]][2], y = y.values[[3]][2]+0.03, label = annot.strs[2], colour=gg_color_hue[3])+
            annotate("text", x = x.values[[3]][3], y = y.values[[3]][3]+0.03, label = annot.strs[3], colour=gg_color_hue[3])+
            annotate("text", x = x.values[[3]][4], y = y.values[[3]][4]+0.03, label = annot.strs[4], colour=gg_color_hue[3])
        }
      }
    }
    if(length(interval.values)>3 | nJP>3){
      print("This annotation feature is not available when the number of joinpoints>3 or the number of multiple intervals selected>3.")
      break
    }
  }
  if(annotation==0){
    plot<-ggplot(data=plotdata, aes(x=plotdata[,yearvar], y=plotdata[,predintvar],group=as.factor(plotdata[,interval]), colour=as.factor(plotdata[,interval]))) + 
      geom_line(data=plotdata, aes(x=plotdata[,yearvar], y=plotdata[,predintvar],group=as.factor(plotdata[,interval]), colour=as.factor(plotdata[,interval])),linetype="solid")+
      geom_point(data=plotdata, aes(x=plotdata[,yearvar], y=1-plotdata[,obsintvar], group=as.factor(plotdata[,interval]), colour=as.factor(plotdata[,interval])))+
      xlab('Year at diagnosis') + ylab('Annual Probability of Cancer Death (%)')+
      ggtitle(title.rch) +
      coord_cartesian(ylim=c(0,1))+
      scale_y_continuous(breaks=seq(0,1,0.1),labels = scales::percent)+
      scale_x_continuous(breaks=seq(min(plotdata[,yearvar],na.rm=T),max(plotdata[,yearvar],na.rm=T),5))+
      scale_color_hue(labels = interval.labels)+
      theme(legend.position="bottom")+
      theme(legend.title=element_blank())+
      theme(plot.title = element_text(hjust = 0.5))
  }
  
  if(!annotation %in% c(0,1)){
    print("The annontation indicator should be either 0 or 1.")
    break
  }
  if(!topanno %in% c(0,1)){
    print("The annontation indicator for the top curve only should be either 0 or 1.")
    break
  }
  return(plot)
}

#########################################################################################################
# plot.surv.year.annotate: a function that returns a plot for Average Change in Cumulative Survival by diagnosis year
# using ggplot
# The annotation feature is available for nJP<=3 and the number of multiple intervals selected <=3
# arguments:
# plotdata - the graph data returned by function download.data with downloadtype="graph".
# fit - joinpoint object containing the model output.
# nJP - the number of joinpoints in the model.
# yearvar - the variable name for year of diagnosis used in argument 'year' of the function joinpoint.
# obscumvar - the variable name for observed relative cumulative survival. The default is "Relative_Survival_Cum"
#             for relative survival data. For cause-specific data, it needs to be changed accordingly.
# predcumvar - the variable name for predicted cumulative survival. The default is "Predicted_Survival_Cum".
# interval - the variable name for year since diagnosis. The default is 'Interval'.
# annotation - the indicator for the annotation feature. The default is 0 (no annotation on the plot).
#########################################################################################################

### define a plot function for Average Change in Cumulative Survival by diagnosis year with annotations
plot.surv.year.annotate<-function(plotdata,fit,nJP,yearvar,obscumvar="Relative_Survival_Cum",predcumvar="Predicted_Survival_Cum",interval="Interval",annotation=0){
  title.acs<-"Average Change in Cumulative Survival by Diagnosis Year"
  interval.values<-as.numeric(unique(plotdata[,interval]))
  interval.labels<-paste(interval.values,"-year Cum. Survival",sep="")

  if(annotation==1){  
    if(length(interval.values)<=3 & nJP<=3){
      ### haz results are the same for all interval values
      jp.loc<-fit$FitList[[nJP+1]]$jp
      x.values<-list()
      y.values<-list()
      cs.aaac<-list()
      annot.strs<-list()
      for(i in 1:length(interval.values)){
        int.i<-interval.values[i]
        cs.aaac[[i]] <- aapc(fit$FitList[[nJP+1]], type="AbsChgSur", interval=int.i)
        annot.strs[[i]]<-paste(sprintf("%.1f",100*cs.aaac[[i]]$estimate),sep="")
        plotdata.i<-plotdata[which(plotdata[,interval]==int.i),]
        if(max(plotdata.i[,yearvar])==max(cs.aaac[[i]]$end.year) | nJP==0){
          end.year<-cs.aaac[[i]]$end.year
          end.year[length(cs.aaac[[i]]$end.year)]<-max(plotdata.i[,yearvar])
          y.values.i<-1/2*((plotdata.i[which(plotdata.i[,yearvar] %in% cs.aaac[[i]]$start.year),predcumvar])+
                             (plotdata.i[which(plotdata.i[,yearvar] %in% end.year),predcumvar]))
          x.values.i<-1/2*(cs.aaac[[i]]$start.year+end.year)
        }
        if(max(plotdata.i[,yearvar])<max(cs.aaac[[i]]$end.year)){
          if(nJP==1){
            if(max(plotdata.i[,yearvar])>jp.loc){
              end.year<-cs.aaac[[i]]$end.year
              end.year[length(cs.aaac[[i]]$end.year)]<-max(plotdata.i[,yearvar])
              y.values.i<-1/2*((plotdata.i[which(plotdata.i[,yearvar] %in% cs.aaac[[i]]$start.year),predcumvar])+
                                 (plotdata.i[which(plotdata.i[,yearvar] %in% end.year),predcumvar]))
              x.values.i<-1/2*(cs.aaac[[i]]$start.year+end.year)
            }else{
              cs.aaac[[i]]<-cs.aaac[[i]][1:nJP,]
              y.values.i<-1/2*((obscumvarplotdata.i[which(plotdata.i[,yearvar] %in% cs.aaac[[i]]$start.year),predcumvar])+
                                 (plotdata.i[which(plotdata.i[,yearvar] %in% cs.aaac[[i]]$end.year),predcumvar]))
              x.values.i<-1/2*(cs.aaac[[i]]$start.year+cs.aaac[[i]]$end.year)
            }
          } ## nJP=1 end
          if(nJP==2){
            if(max(plotdata.i[,yearvar])>max(jp.loc)){
              end.year<-cs.aaac[[i]]$end.year
              end.year[length(cs.aaac[[i]]$end.year)]<-max(plotdata.i[,yearvar])
              y.values.i<-1/2*((plotdata.i[which(plotdata.i[,yearvar] %in% cs.aaac[[i]]$start.year),predcumvar])+
                                 (plotdata.i[which(plotdata.i[,yearvar] %in% end.year),predcumvar]))
              x.values.i<-1/2*(cs.aaac[[i]]$start.year+end.year)
            } else if(max(plotdata.i[,yearvar])<=max(jp.loc) & max(plotdata.i[,yearvar])>min(jp.loc)){
              cs.aaac[[i]]<-cs.aaac[[i]][1:nJP,]
              y.values.i<-1/2*((plotdata.i[which(plotdata.i[,yearvar] %in% cs.aaac[[i]]$start.year),predcumvar])+
                                 (plotdata.i[which(plotdata.i[,yearvar] %in% cs.aaac[[i]]$end.year),predcumvar]))
              x.values.i<-1/2*(cs.aaac[[i]]$start.year+cs.aaac[[i]]$end.year)
            } else{
              cs.aaac[[i]]<-cs.aaac[[i]][1:(nJP-1),]
              end.year<-cs.aaac[[i]]$end.year
              end.year[length(cs.aaac[[i]]$end.year)]<-max(plotdata.i[,yearvar])
              y.values.i<-1/2*((plotdata.i[which(plotdata.i[,yearvar] %in% cs.aaac[[i]]$start.year),predcumvar])+
                                 (plotdata.i[which(plotdata.i[,yearvar] %in% end.year),predcumvar]))
              x.values.i<-1/2*(cs.aaac[[i]]$start.year+end.year)
            }
          } ## nJP=2 end
          if(nJP==3){
            if(max(plotdata.i[,yearvar])>max(jp.loc)){
              end.year<-cs.aaac[[i]]$end.year
              end.year[length(cs.aaac[[i]]$end.year)]<-max(plotdata.i[,yearvar])
              y.values.i<-1/2*((plotdata.i[which(plotdata.i[,yearvar] %in% cs.aaac[[i]]$start.year),predcumvar])+
                                 (plotdata.i[which(plotdata.i[,yearvar] %in% end.year),predcumvar]))
              x.values.i<-1/2*(cs.aaac[[i]]$start.year+end.year)
            } else if(max(plotdata.i[,yearvar])<=max(jp.loc) & max(plotdata.i[,yearvar])>jp.loc[2]){
              cs.aaac[[i]]<-cs.aaac[[i]][1:nJP,]
              end.year<-cs.aaac[[i]]$end.year
              end.year[length(cs.aaac[[i]]$end.year)]<-max(plotdata.i[,yearvar])
              y.values.i<-1/2*((plotdata.i[which(plotdata.i[,yearvar] %in% cs.aaac[[i]]$start.year),predcumvar])+
                                 (plotdata.i[which(plotdata.i[,yearvar] %in% end.year),predcumvar]))
              x.values.i<-1/2*(cs.aaac[[i]]$start.year+end.year)
            } else if(max(plotdata.i[,yearvar])<=jp.loc[2] & max(plotdata.i[,yearvar])>=jp.loc[1]){
              cs.aaac[[i]]<-cs.aaac[[i]][1:(nJP-1),]
              end.year<-cs.aaac[[i]]$end.year
              end.year[length(cs.aaac[[i]]$end.year)]<-max(plotdata.i[,yearvar])
              y.values.i<-1/2*((plotdata.i[which(plotdata.i[,yearvar] %in% cs.aaac[[i]]$start.year),predcumvar])+
                                 (plotdata.i[which(plotdata.i[,yearvar] %in% end.year),predcumvar]))
              x.values.i<-1/2*(cs.aaac[[i]]$start.year+end.year)
            } else {
              cs.aaac[[i]]<-cs.aaac[[i]][1:(nJP-2),]
              end.year<-cs.aaac[[i]]$end.year
              end.year[length(cs.aaac[[i]]$end.year)]<-max(plotdata.i[,yearvar])
              y.values.i<-1/2*((plotdata.i[which(plotdata.i[,yearvar] %in% cs.aaac[[i]]$start.year),predcumvar])+
                                 (plotdata.i[which(plotdata.i[,yearvar] %in% end.year),predcumvar]))
              x.values.i<-1/2*(cs.aaac[[i]]$start.year+end.year)
            }
          } ## nJP=3 end
        } 
        x.values[[i]]<-x.values.i
        y.values[[i]]<-y.values.i
      }## interval.values iteration end
      if(length(interval.values)==1){
        hues=seq(15,375,length=length(interval.values)+1)
        gg_color_hue<-hcl(h=hues,l=65,c=100)[1:length(interval.values)]
        plot<-ggplot(data=plotdata, aes(x=plotdata[,yearvar], y=plotdata[,predcumvar],group=as.factor(plotdata[,interval]), colour=as.factor(plotdata[,interval]))) + 
          geom_line(data=plotdata, aes(x=plotdata[,yearvar], y=plotdata[,predcumvar],group=as.factor(plotdata[,interval]), colour=as.factor(plotdata[,interval])),linetype="solid")+
          geom_point(data=plotdata, aes(x=plotdata[,yearvar], y=plotdata[,obscumvar], group=as.factor(plotdata[,interval]), colour=as.factor(plotdata[,interval])))+
          xlab('Year at diagnosis') + ylab('Cumulative Survival (%)')+
          theme_bw()+
          ggtitle(title.acs) +
          coord_cartesian(ylim=c(0,1))+
          scale_y_continuous(breaks=seq(0,1,0.1),labels = scales::percent)+
          scale_x_continuous(breaks=seq(min(plotdata[,yearvar],na.rm=T),max(plotdata[,yearvar],na.rm=T),5))+
          scale_color_hue(labels = interval.labels)+
          theme(legend.position="bottom")+
          theme(legend.title=element_blank())+
          theme(plot.title = element_text(hjust = 0.5))+
          annotate("text", x = x.values[[1]][1], y = y.values[[1]][1]+0.03, label = annot.strs[[1]][1], colour=gg_color_hue[1])+
          annotate("text", x = x.values[[1]][2], y = y.values[[1]][2]+0.03, label = annot.strs[[1]][2], colour=gg_color_hue[1])+
          annotate("text", x = x.values[[1]][3], y = y.values[[1]][3]+0.03, label = annot.strs[[1]][3], colour=gg_color_hue[1])+
          annotate("text", x = x.values[[1]][4], y = y.values[[1]][4]+0.03, label = annot.strs[[1]][4], colour=gg_color_hue[1])
      }
      if(length(interval.values)==2){
        hues=seq(15,375,length=length(interval.values)+1)
        gg_color_hue<-hcl(h=hues,l=65,c=100)[1:length(interval.values)]
        plot<-ggplot(data=plotdata, aes(x=plotdata[,yearvar], y=plotdata[,predcumvar],group=as.factor(plotdata[,interval]), colour=as.factor(plotdata[,interval]))) + 
          geom_line(data=plotdata, aes(x=plotdata[,yearvar], y=plotdata[,predcumvar],group=as.factor(plotdata[,interval]), colour=as.factor(plotdata[,interval])),linetype="solid")+
          geom_point(data=plotdata, aes(x=plotdata[,yearvar], y=plotdata[,obscumvar], group=as.factor(plotdata[,interval]), colour=as.factor(plotdata[,interval])))+
          xlab('Year at diagnosis') + ylab('Cumulative Survival (%)')+
          theme_bw()+
          ggtitle(title.acs) +
          coord_cartesian(ylim=c(0,1))+
          scale_y_continuous(breaks=seq(0,1,0.1),labels = scales::percent)+
          scale_x_continuous(breaks=seq(min(plotdata[,yearvar],na.rm=T),max(plotdata[,yearvar],na.rm=T),5))+
          scale_color_hue(labels = interval.labels)+
          theme(legend.position="bottom")+
          theme(legend.title=element_blank())+
          theme(plot.title = element_text(hjust = 0.5))+
          annotate("text", x = x.values[[1]][1], y = y.values[[1]][1]+0.03, label = annot.strs[[1]][1], colour=gg_color_hue[1])+
          annotate("text", x = x.values[[1]][2], y = y.values[[1]][2]+0.03, label = annot.strs[[1]][2], colour=gg_color_hue[1])+
          annotate("text", x = x.values[[1]][3], y = y.values[[1]][3]+0.03, label = annot.strs[[1]][3], colour=gg_color_hue[1])+
          annotate("text", x = x.values[[1]][4], y = y.values[[1]][4]+0.03, label = annot.strs[[1]][4], colour=gg_color_hue[1])+
          annotate("text", x = x.values[[2]][1], y = y.values[[2]][1]+0.03, label = annot.strs[[2]][1], colour=gg_color_hue[2])+
          annotate("text", x = x.values[[2]][2], y = y.values[[2]][2]+0.03, label = annot.strs[[2]][2], colour=gg_color_hue[2])+
          annotate("text", x = x.values[[2]][3], y = y.values[[2]][3]+0.03, label = annot.strs[[2]][3], colour=gg_color_hue[2])+
          annotate("text", x = x.values[[2]][4], y = y.values[[2]][4]+0.03, label = annot.strs[[2]][4], colour=gg_color_hue[2])    
      }
      if(length(interval.values)==3){
        hues=seq(15,375,length=length(interval.values)+1)
        gg_color_hue<-hcl(h=hues,l=65,c=100)[1:length(interval.values)]
        plot<-ggplot(data=plotdata, aes(x=plotdata[,yearvar], y=plotdata[,predcumvar],group=as.factor(plotdata[,interval]), colour=as.factor(plotdata[,interval]))) + 
          geom_line(data=plotdata, aes(x=plotdata[,yearvar], y=plotdata[,predcumvar],group=as.factor(plotdata[,interval]), colour=as.factor(plotdata[,interval])),linetype="solid")+
          geom_point(data=plotdata, aes(x=plotdata[,yearvar], y=plotdata[,obscumvar], group=as.factor(plotdata[,interval]), colour=as.factor(plotdata[,interval])))+
          xlab('Year at diagnosis') + ylab('Cumulative Survival (%)')+
          theme_bw()+
          ggtitle(title.acs) +
          coord_cartesian(ylim=c(0,1))+
          scale_y_continuous(breaks=seq(0,1,0.1),labels = scales::percent)+
          scale_x_continuous(breaks=seq(min(plotdata[,yearvar],na.rm=T),max(plotdata[,yearvar],na.rm=T),5))+
          scale_color_hue(labels = interval.labels)+
          theme(legend.position="bottom")+
          theme(legend.title=element_blank())+
          theme(plot.title = element_text(hjust = 0.5))+
          annotate("text", x = x.values[[1]][1], y = y.values[[1]][1]+0.03, label = annot.strs[[1]][1], colour=gg_color_hue[1])+
          annotate("text", x = x.values[[1]][2], y = y.values[[1]][2]+0.03, label = annot.strs[[1]][2], colour=gg_color_hue[1])+
          annotate("text", x = x.values[[1]][3], y = y.values[[1]][3]+0.03, label = annot.strs[[1]][3], colour=gg_color_hue[1])+
          annotate("text", x = x.values[[1]][4], y = y.values[[1]][4]+0.03, label = annot.strs[[1]][4], colour=gg_color_hue[1])+
          annotate("text", x = x.values[[2]][1], y = y.values[[2]][1]+0.03, label = annot.strs[[2]][1], colour=gg_color_hue[2])+
          annotate("text", x = x.values[[2]][2], y = y.values[[2]][2]+0.03, label = annot.strs[[2]][2], colour=gg_color_hue[2])+
          annotate("text", x = x.values[[2]][3], y = y.values[[2]][3]+0.03, label = annot.strs[[2]][3], colour=gg_color_hue[2])+
          annotate("text", x = x.values[[2]][4], y = y.values[[2]][4]+0.03, label = annot.strs[[2]][4], colour=gg_color_hue[2])+
          annotate("text", x = x.values[[3]][1], y = y.values[[3]][1]+0.03, label = annot.strs[[3]][1], colour=gg_color_hue[3])+
          annotate("text", x = x.values[[3]][2], y = y.values[[3]][2]+0.03, label = annot.strs[[3]][2], colour=gg_color_hue[3])+
          annotate("text", x = x.values[[3]][3], y = y.values[[3]][3]+0.03, label = annot.strs[[3]][3], colour=gg_color_hue[3])+
          annotate("text", x = x.values[[3]][4], y = y.values[[3]][4]+0.03, label = annot.strs[[3]][4], colour=gg_color_hue[3])
      }
    }
    if(length(interval.values)>3 | nJP>3){
      print("This annotation feature is not available when the number of joinpoints>3 or the number of multiple intervals selected>3.")
      break
    }
  }
  if(annotation==0){
    plot<-ggplot(data=plotdata, aes(x=plotdata[,yearvar], y=plotdata[,predcumvar],group=as.factor(plotdata[,interval]), colour=as.factor(plotdata[,interval]))) + 
      geom_line(data=plotdata, aes(x=plotdata[,yearvar], y=plotdata[,predcumvar],group=as.factor(plotdata[,interval]), colour=as.factor(plotdata[,interval])),linetype="solid")+
      geom_point(data=plotdata, aes(x=plotdata[,yearvar], y=plotdata[,obscumvar], group=as.factor(plotdata[,interval]), colour=as.factor(plotdata[,interval])))+
      xlab('Year at diagnosis') + ylab('Cumulative Survival (%)')+
      theme_bw()+
      ggtitle(title.acs) +
      coord_cartesian(ylim=c(0,1))+
      scale_y_continuous(breaks=seq(0,1,0.1),labels = scales::percent)+
      scale_x_continuous(breaks=seq(min(plotdata[,yearvar],na.rm=T),max(plotdata[,yearvar],na.rm=T),5))+
      scale_color_hue(labels = interval.labels)+
      theme(legend.position="bottom")+
      theme(legend.title=element_blank())+
      theme(plot.title = element_text(hjust = 0.5))
  }
  if(!annotation %in% c(0,1)){
    print("The annontation indicator should be either 0 or 1.")
    break
  }
  return(plot)
}


#########################################################################################################
# plot.surv.int.multiyears: a function that returns a plot for Cumulative Survival by Interval supporting
# multiple selected years
# arguments:
# plotdata - 1. the graph data returned by function download.data with downloadtype="graph" and all intervals needed for int.col.
#            2. the full data returned by function download.data with downloadtype="full".
#            Either of the above would work.
# fit - joinpoint object containing the model output.
# nJP - the number of joinpoints in the model.
# yearvar - the variable name for year of diagnosis used in argument 'year' of the function joinpoint.
# obscumvar - the variable name for observed relative cumulative survival. The default is "Relative_Survival_Cum"
#             for relative survival data. For cause-specific data, it needs to be changed accordingly.
# predcumvar - the variable name for predicted cumulative survival. The default is "Predicted_Survival_Cum".
# interval - the variable name for year since diagnosis. The default is 'Interval'.
# year.col - the year values selected for the plot. The default is NULL.
#########################################################################################################

### define a plot function for Cumulative Survival by Interval
plot.surv.int.multiyears<-function(plotdata,fit,nJP,yearvar,obscumvar="Relative_Survival_Cum",predcumvar="Predicted_Survival_Cum",interval="Interval",year.col=NULL){
  year.col.str<- paste(year.col,collapse=", ",sep="")
  if(length(year.col)<=3){
    title.survint<-paste("Cumulative Survival by Interval for ",year.col.str,sep="")
  }else{
    title.survint<-"Cumulative Survival by Interval"
  }
  x.combo<-function(x){
    paste(x,collapse="_",sep="")
  }
  cohort.combo <- rep(NA,dim(plotdata)[1])
  cohort.combo <- apply(data.frame(plotdata[,1:(which(colnames(plotdata)==interval)-1)]), 1, x.combo)
  cohort.combo.int0<-(paste(unique(cohort.combo),"_00",sep=""))
  int.2d<-sprintf("%02d",plotdata[,interval])
  cohort.combo.int<-paste(cohort.combo,"_",int.2d,sep="")
  plotdata.add<-plotdata[which(plotdata[,interval]==1),]
  plotdata.add[,interval]<-0
  plotdata.add[,c(obscumvar,predcumvar)]<-1
  cols.keep<-c(colnames(plotdata)[1:which(colnames(plotdata)==interval)],obscumvar,predcumvar)
  plotdata.add[,colnames(plotdata)[which(!colnames(plotdata) %in% cols.keep)]]<-NA
  plotdata[,"cohort.combo.int"]<-cohort.combo.int
  plotdata.add[,"cohort.combo.int"]<-cohort.combo.int0
  plotdata.merge<-rbind(plotdata,plotdata.add)
  plotdata.sort<-plotdata.merge[order(plotdata.merge$cohort.combo.int),]
  
  plotdata.sub<-plotdata.sort[which((plotdata.sort[,yearvar] %in% year.col) & !is.na(plotdata.sort[,obscumvar])),]
  int.max<-max(plotdata.sub[,interval],na.rm=T)
  if(int.max<=12){
    xbreaks.by<-1
  }else{
    xbreaks.by<-2
  }
  plot<-ggplot(data=plotdata.sub, aes(x=plotdata.sub[,interval], y=plotdata.sub[,predcumvar],group=as.factor(plotdata.sub[,yearvar]), colour=as.factor(plotdata.sub[,yearvar]))) + 
    geom_line(data=plotdata.sub, aes(x=plotdata.sub[,interval], y=plotdata.sub[,predcumvar],group=as.factor(plotdata.sub[,yearvar]), colour=as.factor(plotdata.sub[,yearvar])),linetype="solid")+
    geom_point(data=plotdata.sub, aes(x=plotdata.sub[,interval], y=plotdata.sub[,obscumvar], group=as.factor(plotdata.sub[,yearvar]), colour=as.factor(plotdata.sub[,yearvar])))+
    xlab('Interval') + ylab('Cumulative Survival (%)')+
    theme_bw()+
    ggtitle(title.survint) +
    coord_cartesian(ylim=c(0,1))+
    scale_y_continuous(breaks=seq(0,1,0.1),labels = scales::percent)+
    scale_x_continuous(breaks=seq(min(plotdata.sub[,interval],na.rm=T),max(plotdata.sub[,interval],na.rm=T),xbreaks.by))+
    scale_color_hue(labels = year.col)+
    theme(legend.position="bottom")+
    theme(legend.title=element_blank())+
    theme(plot.title = element_text(hjust = 0.5))
}

##################################################################################################################################################
### test using example data sets
library(JPSurv)
library(ggplot2)
setwd("P:/srab/surv/JPSurvival/Fanni_work/06052019.Rpackage/")

# example1 relative survival data
DICsample <- dictionary.overview("SEER9_Survival_6CancerSitesByStage_1975_2007.dic")
head(DICsample)

input <- joinpoint.seerdata(seerfilename="SEER9_Survival_6CancerSitesByStage_1975_2007", 
                            newvarnames=c("Site", "Year_of_diagnosis_7507_individual"), NoFit=T,
                            UseVarLabelsInData=c("Sites: CR LB B O P T","Year_of_diagnosis_7507_individual"))
head(input)
colnames(input)[2]<-"Site"
colnames(input)[3]<-"Sex"
colnames(input)[4]<-"Stage"
head(input)

subsetStr<-"Site==0 & Sex==0 & Stage==0"
fit <- joinpoint(input, subsetStr,
                 year="Year_of_diagnosis_7507_individual", observedrelsurv="Relative_Survival_Cum",
                 model.form = NULL, maxnum.jp = 3, proj.year.num = 5)

yearvar<-"Year_of_diagnosis_7507_individual"
obsintvar<-"Relative_Survival_Interval"
predintvar<-"Predicted_ProbDeath_Int"
obscumvar<-"Relative_Survival_Cum"
predcumvar<-"Predicted_Survival_Cum"
interval<-"Interval"
annotation<-1

#nJP<-1
nJP<-3

# test download.data
data.graph5<-download.data(input,fit,nJP,yearvar,"graph",subsetStr,interval="Interval",int.col=c(5))
data.graph15<-download.data(input,fit,nJP,yearvar,"graph",subsetStr,interval="Interval",int.col=c(1,5))
data.graph135<-download.data(input,fit,nJP,yearvar,"graph",subsetStr,interval="Interval",int.col=c(1,3,5))
data.graph1to20<-download.data(input,fit,nJP,yearvar,"graph",subsetStr,interval="Interval",int.col=c(1:20))
data.graph1tomax<-download.data(input,fit,nJP,yearvar,"graph",subsetStr,interval="Interval",int.col=c(1:max(input$Interval,na.rm=T)))
data.full<-download.data(input,fit,nJP,yearvar,"full",subsetStr,interval="Interval")

# test plot.dying.year.annotate
pp5.anno<-plot.dying.year.annotate(data.graph5,fit,nJP,yearvar,obsintvar,predintvar,interval,annotation)
pp5<-plot.dying.year.annotate(data.graph5,fit,nJP,yearvar,obsintvar,predintvar,interval)
pp15.anno<-plot.dying.year.annotate(data.graph15,fit,nJP,yearvar,obsintvar,predintvar,interval,annotation)
pp15<-plot.dying.year.annotate(data.graph15,fit,nJP,yearvar,obsintvar,predintvar,interval)
pp135.anno<-plot.dying.year.annotate(data.graph135,fit,nJP,yearvar,obsintvar,predintvar,interval,annotation)
pp135<-plot.dying.year.annotate(data.graph135,fit,nJP,yearvar,obsintvar,predintvar,interval)
pp.1to20<-plot.dying.year.annotate(data.graph1to20,fit,nJP,yearvar,obsintvar,predintvar,interval,annotation)
pp.1to20<-plot.dying.year.annotate(data.graph1to20,fit,nJP,yearvar,obsintvar,predintvar,interval)

# test plot.surv.year.annotate
pp5.anno<-plot.surv.year.annotate(data.graph5,fit,nJP,yearvar,obscumvar,predcumvar,interval,annotation)
pp5<-plot.surv.year.annotate(data.graph5,fit,nJP,yearvar,obscumvar,predcumvar,interval)
pp15.anno<-plot.surv.year.annotate(data.graph15,fit,nJP,yearvar,obscumvar,predcumvar,interval,annotation)
pp15<-plot.surv.year.annotate(data.graph15,fit,nJP,yearvar,obscumvar,predcumvar,interval)
pp135.anno<-plot.surv.year.annotate(data.graph135,fit,nJP,yearvar,obscumvar,predcumvar,interval,annotation)
pp135<-plot.surv.year.annotate(data.graph135,fit,nJP,yearvar,obscumvar,predcumvar,interval)
pp.1to20<-plot.surv.year.annotate(data.graph1to20,fit,nJP,yearvar,obscumvar,predcumvar,interval,annotation)
pp.1to20<-plot.surv.year.annotate(data.graph1to20,fit,nJP,yearvar,obscumvar,predcumvar,interval)

# test plot.surv.int.multiyears
pp.int.test1<-plot.surv.int.multiyears(data.full,fit,nJP,yearvar,obscumvar,predcumvar,interval,year.col=c(1985,1990))
pp.int.test2<-plot.surv.int.multiyears(data.graph1tomax,fit,nJP,yearvar,obscumvar,predcumvar,interval,year.col=seq(1975,1990,5))
pp.int.test3<-plot.surv.int.multiyears(data.graph1tomax,fit,nJP,yearvar,obscumvar,predcumvar,interval,year.col=seq(1998,2008,4))


#save.f<-"P:/srab/surv/JPSurvival/Fanni_work/06052019.Rpackage/test_download_graph135_example1.csv"
#write.csv(data.graph135,save.f,row.names=F,quote=F)
#save.f<-"P:/srab/surv/JPSurvival/Fanni_work/06052019.Rpackage/test_download_full_example1.csv"
#write.csv(data.full,save.f,row.names=F,quote=F)
#save.f<-"P:/srab/surv/JPSurvival/Fanni_work/06052019.Rpackage/test_example1_pp.jpg"
#ggsave(save.f,pp15.anno)

# example2  cause-specific data
DICsample = dictionary.overview("Breast_causespecific.dic")
head(DICsample)
seerdata_cause = joinpoint.seerdata(seerfilename="Breast_causespecific", 
                                    newvarnames=c("Year_of_diagnosis_1975"),NoFit=T,
                                    UseVarLabelsInData=c("Year_of_diagnosis_1975"))
subsetStr="Year_of_diagnosis_1975 >= 1975 & Breast_stage==1 & Age_groups == 2"
fit.cs <- joinpoint(seerdata_cause, 
                    subsetStr,
                    year="Year_of_diagnosis_1975",observedrelsurv="CauseSpecific_Survival_Cum",
                    model.form = NULL,
                    maxnum.jp = 3,proj.year.num=0)


nJP<-2
yearvar<-"Year_of_diagnosis_1975"
obsintvar<-"CauseSpecific_Survival_Interval"
predintvar<-"Predicted_ProbDeath_Int"
obscumvar<-"CauseSpecific_Survival_Cum"
predcumvar<-"Predicted_Survival_Cum"
interval<-"Interval"
annotation<-1

# test download.data
data.graph5<-download.data(seerdata_cause,fit.cs,nJP,yearvar,"graph",subsetStr,interval="Interval",int.col=c(5))
data.graph15<-download.data(seerdata_cause,fit.cs,nJP,yearvar,"graph",subsetStr,interval="Interval",int.col=c(1,5))
data.graph135<-download.data(seerdata_cause,fit.cs,nJP,yearvar,"graph",subsetStr,interval="Interval",int.col=c(1,3,5))
data.graph1to20<-download.data(seerdata_cause,fit.cs,nJP,yearvar,"graph",subsetStr,interval="Interval",int.col=c(1:20))
data.graph1tomax<-download.data(seerdata_cause,fit.cs,nJP,yearvar,"graph",subsetStr,interval="Interval",int.col=c(1:max(input$Interval,na.rm=T)))
data.full<-download.data(seerdata_cause,fit.cs,nJP,yearvar,"full",subsetStr)

# test plot.dying.year.annotate
pp5.anno<-plot.dying.year.annotate(data.graph5,fit.cs,nJP,yearvar,obsintvar,predintvar,interval,annotation)
pp5<-plot.dying.year.annotate(data.graph5,fit.cs,nJP,yearvar,obsintvar,predintvar,interval)
pp15.anno<-plot.dying.year.annotate(data.graph15,fit.cs,nJP,yearvar,obsintvar,predintvar,interval,annotation)
pp15<-plot.dying.year.annotate(data.graph15,fit.cs,nJP,yearvar,obsintvar,predintvar,interval)
pp135.anno<-plot.dying.year.annotate(data.graph135,fit.cs,nJP,yearvar,obsintvar,predintvar,interval,annotation)
pp135<-plot.dying.year.annotate(data.graph135,fit.cs,nJP,yearvar,obsintvar,predintvar,interval)
pp.1to20.anno<-plot.dying.year.annotate(data.graph1to20,fit.cs,nJP,yearvar,obsintvar,predintvar,interval,annotation)
pp.1to20<-plot.dying.year.annotate(data.graph1to20,fit.cs,nJP,yearvar,obsintvar,predintvar,interval)

# test plot.surv.year.annotate
pp5.anno<-plot.surv.year.annotate(data.graph5,fit.cs,nJP,yearvar,obscumvar,predcumvar,interval,annotation)
pp5<-plot.surv.year.annotate(data.graph5,fit.cs,nJP,yearvar,obscumvar,predcumvar,interval)
pp15.anno<-plot.surv.year.annotate(data.graph15,fit.cs,nJP,yearvar,obscumvar,predcumvar,interval,annotation)
pp15<-plot.surv.year.annotate(data.graph15,fit.cs,nJP,yearvar,obscumvar,predcumvar,interval)
pp135.anno<-plot.surv.year.annotate(data.graph135,fit.cs,nJP,yearvar,obscumvar,predcumvar,interval,annotation)
pp135<-plot.surv.year.annotate(data.graph135,fit.cs,nJP,yearvar,obscumvar,predcumvar,interval)
pp.1to20<-plot.surv.year.annotate(data.graph1to20,fit.cs,nJP,yearvar,obscumvar,predcumvar,interval,annotation)
pp.1to20<-plot.surv.year.annotate(data.graph1to20,fit.cs,nJP,yearvar,obscumvar,predcumvar,interval)

# test plot.surv.int.multiyears
pp.int.test1<-plot.surv.int.multiyears(data.full,fit,nJP,yearvar,obscumvar,predcumvar,interval,year.col=c(1985,1990))
pp.int.test2<-plot.surv.int.multiyears(data.graph1tomax,fit,nJP,yearvar,obscumvar,predcumvar,interval,year.col=seq(1975,1990,3))
pp.int.test3<-plot.surv.int.multiyears(data.graph1tomax,fit,nJP,yearvar,obscumvar,predcumvar,interval,year.col=seq(1998,2008,4))

#save.f<-"P:/srab/surv/JPSurvival/Fanni_work/06052019.Rpackage/test_example2_pp.jpg"
#ggsave(save.f,pp135)
#save.f<-"P:/srab/surv/JPSurvival/Fanni_work/06052019.Rpackage/test_download_graph135_example2.csv"
#write.csv(data.graph15,save.f,row.names=F,quote=F)
#save.f<-"P:/srab/surv/JPSurvival/Fanni_work/06052019.Rpackage/test_download_full_example2.csv"
#write.csv(data.full,save.f,row.names=F,quote=F)


