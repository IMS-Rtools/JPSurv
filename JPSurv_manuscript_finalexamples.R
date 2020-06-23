#################################################################################################
### P:\srab\surv\JPSurvival\Fanni_work\Paper\JPSurv_manuscript_finalexamples.R
### Created by Fanni Zhang 06/17/2020
### - JPSurv runs for final examples in paper P:/srab/surv/JPSurvival/Fanni_work/Paper/JPsurvpaper_2020.06.16.docx
### - Create plot functions plot.dying.year.full and plot.surv.year.full with projection feature. 
### Request: Tables and figures in P:/srab/surv/JPSurvival/Fanni_work/Paper/Tables_Figures_2020.06.16.docx 
###          P:/srab/surv/JPSurvival/Fanni_work/Paper/Request_06152020_paper and finalizing Tables and Figures
#################################################################################################

library(JPSurv)
library(ggplot2)
library(ggrepel)
setwd("P:/srab/surv/JPSurvival/Fanni_work/Paper/")
#################################################################################################
# paper example relative survival data
DICsample <- dictionary.overview("Melanoma_Pancreas.dic")
head(DICsample)

input <- joinpoint.seerdata(seerfilename="Melanoma_Pancreas", 
                            newvarnames=c("Year_of_diagnosis_1975"),
                            UseVarLabelsInData=c("Year_of_diagnosis_1975"))
head(input)

subsetStr<-list()
subsetStr[[1]]<-"Stage_AllLRD==0 & Pancreas_Melanoma ==1"
subsetStr[[2]]<-"Stage_AllLRD==1 & Pancreas_Melanoma ==1"
subsetStr[[3]]<-"Stage_AllLRD==2 & Pancreas_Melanoma ==1"
subsetStr[[4]]<-"Stage_AllLRD==0 & Pancreas_Melanoma ==0"
subsetStr[[5]]<-"Stage_AllLRD==1 & Pancreas_Melanoma ==0"
subsetStr[[6]]<-"Stage_AllLRD==2 & Pancreas_Melanoma ==0"
subsetStr[[7]]<-"Stage_AllLRD==3 & Pancreas_Melanoma ==1"
subsetStr[[8]]<-"Stage_AllLRD==3 & Pancreas_Melanoma ==0"

titleStr<-list()
titleStr[[1]]<-"Melanoma - Local"
titleStr[[2]]<-"Melanoma - Regional"
titleStr[[3]]<-"Melanoma - Distant"
titleStr[[4]]<-"Pancreas - Local"
titleStr[[5]]<-"Pancreas - Regional"
titleStr[[6]]<-"Pancreas - Distant"
titleStr[[7]]<-"Melanoma - All Stages"
titleStr[[8]]<-"Pancreas - All Stages"

out.sum<-list()
out.surv<-list()
out.dying<-list()

yearvar<-"Year_of_diagnosis_1975"
obscumvar<-"Relative_Survival_Cum"
predcumvar<-"Predicted_Survival_Cum"
obsintvar<-"Relative_Survival_Interval"
predintvar<-"Predicted_ProbDeath_Int"
interval<-"Interval"
site<-c(rep("Melanoma",3),rep("Pancreas",3),"Melanoma","Pancreas")
stage<-c(rep(c("Local","Regional","Distant"),2),"All stages","All stages")

### first run using maxjp=3 for Melanoma and Pancreas
### figures created using plot.dying.year.annotate and plot.surv.year.annotate without projection
maxjp<-3

for(i in 1:length(subsetStr)){
  subset<-subsetStr[[i]]
  print(c(i,subset))
  fit <- joinpoint(input, subset,
                   year=yearvar, observedrelsurv="Relative_Survival_Cum",
                   model.form = NULL, maxnum.jp = maxjp, proj.year.num = 5)
  #Cancer site	Stage 	No. JPs. Final Model 	No. Alive	Start Year	End Year 	AACS(1) 	95% C.I.	AACS(5) 	95% C.I.	PCAPD	95% C.I.
  cnames<-c("Cancer site","Stage","No. JPs Final Model","No. Alive","Start Year","End Year",
            "AACS(1)","AACS(1) 95% C.I.","AACS(5)","AACS(5) 95% C.I.","PCAPD","PCAPD 95% C.I.")
  sum<-array("",dim=c((length(fit$jp)+1),length(cnames)))
  colnames(sum)<-cnames
  
  sum[1,1]<-site[i]
  sum[1,2]<-stage[i]
  if(is.null(fit$jp)){
    sum[1,3]<-0
  }else{ 
    sum[1,3]<-length(fit$jp)
  }
  nJP<-length(fit$jp)
  data.graph<-download.data(input,fit,nJP,yearvar,"graph",subset,interval,int.select=c(1,5,10))
  title.rch<-paste(titleStr[[i]]," \n Annual Probability of Dying of Cancer by Diagnosis Year",sep="")
  title.acs<-paste(titleStr[[i]]," \n Relative Survival by Diagnosis Year",sep="")
  dying.out.anno<-plot.dying.year.annotate(data.graph,fit,nJP,yearvar,obsintvar,predintvar,interval,annotation=1,topanno=1,trend=1,title=title.rch)
  surv.out.anno<-plot.surv.year.annotate(data.graph,fit,nJP,yearvar,obscumvar,predcumvar,interval,annotation=1,trend=1,title=title.acs)
  dyingplot.f<-paste("P:/srab/surv/JPSurvival/Fanni_work/Paper/plot/dying_plot_anno_",titleStr[[i]],".jpg",sep="")
  survplot.f<-paste("P:/srab/surv/JPSurvival/Fanni_work/Paper/plot/surv_plot_anno_",titleStr[[i]],".jpg",sep="")
  ggsave(dyingplot.f,dying.out.anno$plot_anno, width=6.76, height=5.32)
  ggsave(survplot.f,surv.out.anno$plot_anno, width=6.76, height=5.32)
  
  input.sub<-subset(input,eval(parse(text=subset)))
  sum[1,4]<-input.sub$Alive_at_Start[1]
  sum[,5]<-dying.out.anno$trends[[1]]$start.year
  sum[,6]<-dying.out.anno$trends[[1]]$end.year
  
  ### surv 1-year
  sum[,7]<-surv.out.anno$trends[[1]]$estimate*100
  survci.l<-round(surv.out.anno$trends[[1]]$lowCI*100,3)
  survci.u<-round(surv.out.anno$trends[[1]]$upCI*100,3)
  survci.l<-sprintf("%.3f",survci.l)  
  survci.u<-sprintf("%.3f",survci.u)
  survci.str<-paste("(",survci.l,", ",survci.u,")",sep="")
  sum[,8]<-survci.str
  
  ### surv 5-year
  sum[,9]<-surv.out.anno$trends[[2]]$estimate*100
  survci.l<-round(surv.out.anno$trends[[2]]$lowCI*100,3)
  survci.u<-round(surv.out.anno$trends[[2]]$upCI*100,3)
  survci.l<-sprintf("%.3f",survci.l)  
  survci.u<-sprintf("%.3f",survci.u)
  survci.str<-paste("(",survci.l,", ",survci.u,")",sep="")
  sum[,10]<-survci.str
  
  ### dying
  sum[,11]<-dying.out.anno$trends[[1]]$estimate*100
  ci.l<-round(dying.out.anno$trends[[1]]$lowCI*100,3)
  ci.u<-round(dying.out.anno$trends[[1]]$upCI*100,3)
  ci.l<-sprintf("%.3f",ci.l)  
  ci.u<-sprintf("%.3f",ci.u)
  ci.str<-paste("(",ci.l,", ",ci.u,")",sep="")
  sum[,12]<-ci.str
  
  out.sum[[i]]<-sum
  out.surv[[i]]<-surv.out.anno
  out.dying[[i]]<-dying.out.anno
}

out.sum.combine<-rbind(out.sum[[1]],out.sum[[2]],out.sum[[3]],out.sum[[4]],out.sum[[5]],out.sum[[6]])
save.f<-"P:/srab/surv/JPSurvival/Fanni_work/Paper/JPsurv_manuscript_Melanoma_Pancreas_trend_summary.txt.xls"
write.table(out.sum.combine,save.f,row.names=F,quote=F,sep="\t")


##################################################################################################
### FemaleBreast_Liver
# paper example relative survival data
DICsample <- dictionary.overview("FemaleBreast_Liver.dic")
head(DICsample)

input <- joinpoint.seerdata(seerfilename="FemaleBreast_Liver", 
                            newvarnames=c("Year_of_diagnosis_1975"),
                            UseVarLabelsInData=c("Year_of_diagnosis_1975"))
head(input)

subsetStr<-list()
subsetStr[[1]]<-"Stage_AllLRD==0 & Liver_FBreast ==1"
subsetStr[[2]]<-"Stage_AllLRD==1 & Liver_FBreast ==1"
subsetStr[[3]]<-"Stage_AllLRD==2 & Liver_FBreast ==1"
subsetStr[[4]]<-"Stage_AllLRD==0 & Liver_FBreast ==0"
subsetStr[[5]]<-"Stage_AllLRD==1 & Liver_FBreast ==0"
subsetStr[[6]]<-"Stage_AllLRD==2 & Liver_FBreast ==0"
subsetStr[[7]]<-"Stage_AllLRD==3 & Liver_FBreast ==1"
subsetStr[[8]]<-"Stage_AllLRD==3 & Liver_FBreast ==0"

titleStr<-list()
titleStr[[1]]<-"Female Breast - Local"
titleStr[[2]]<-"Female Breast - Regional"
titleStr[[3]]<-"Female Breast - Distant"
titleStr[[4]]<-"Liver - Local"
titleStr[[5]]<-"Liver - Regional"
titleStr[[6]]<-"Liver - Distant"
titleStr[[7]]<-"Female Breast - All Stages"
titleStr[[8]]<-"Liver - All Stages"

out.sum<-list()
out.surv<-list()
out.dying<-list()

yearvar<-"Year_of_diagnosis_1975"
obscumvar<-"Relative_Survival_Cum"
predcumvar<-"Predicted_Survival_Cum"
obsintvar<-"Relative_Survival_Interval"
predintvar<-"Predicted_ProbDeath_Int"
interval<-"Interval"
site<-c(rep("Female Breast",3),rep("Liver",3),"Female Breast","Liver")
stage<-c(rep(c("Local","Regional","Distant"),2),"All stages","All stages")  

maxjp<-3
### figures created using plot.dying.year.annotate and plot.surv.year.annotate without projection

for(i in 1:length(subsetStr)){
  subset<-subsetStr[[i]]
  print(c(i,subset))
  fit <- joinpoint(input, subset,
                   year=yearvar, observedrelsurv="Relative_Survival_Cum",
                   model.form = NULL, maxnum.jp = maxjp, proj.year.num = 5)
  #Cancer site	Stage 	No. JPs. Final Model 	No. Alive	Start Year	End Year 	AACS(1) 	95% C.I.	AACS(5) 	95% C.I.	PCAPD	95% C.I.
  cnames<-c("Cancer site","Stage","No. JPs Final Model","No. Alive","Start Year","End Year",
            "AACS(1)","AACS(1) 95% C.I.","AACS(5)","AACS(5) 95% C.I.","PCAPD","PCAPD 95% C.I.")
  sum<-array("",dim=c((length(fit$jp)+1),length(cnames)))
  colnames(sum)<-cnames
  
  sum[1,1]<-site[i]
  sum[1,2]<-stage[i]
  if(is.null(fit$jp)){
    sum[1,3]<-0
  }else{ 
    sum[1,3]<-length(fit$jp)
  }
  nJP<-length(fit$jp)
  data.graph<-download.data(input,fit,nJP,yearvar,"graph",subset,interval,int.select=c(1,5,10))
  title.rch<-paste(titleStr[[i]]," \n Annual Probability of Dying of Cancer by Diagnosis Year",sep="")
  title.acs<-paste(titleStr[[i]]," \n Relative Survival by Diagnosis Year",sep="")
  dying.out.anno<-plot.dying.year.annotate(data.graph,fit,nJP,yearvar,obsintvar,predintvar,interval,annotation=1,topanno=1,trend=1,title=title.rch)
  surv.out.anno<-plot.surv.year.annotate(data.graph,fit,nJP,yearvar,obscumvar,predcumvar,interval,annotation=1,trend=1,title=title.acs)
  dyingplot.f<-paste("P:/srab/surv/JPSurvival/Fanni_work/Paper/plot/dying_plot_anno_",titleStr[[i]],".jpg",sep="")
  survplot.f<-paste("P:/srab/surv/JPSurvival/Fanni_work/Paper/plot/surv_plot_anno_",titleStr[[i]],".jpg",sep="")
  ggsave(dyingplot.f,dying.out.anno$plot_anno, width=6.76, height=5.32)
  ggsave(survplot.f,surv.out.anno$plot_anno, width=6.76, height=5.32)
  
  input.sub<-subset(input,eval(parse(text=subset)))
  sum[1,4]<-input.sub$Alive_at_Start[1]
  sum[,5]<-dying.out.anno$trends[[1]]$start.year
  sum[,6]<-dying.out.anno$trends[[1]]$end.year
  
  ### surv 1-year
  sum[,7]<-surv.out.anno$trends[[1]]$estimate*100
  survci.l<-round(surv.out.anno$trends[[1]]$lowCI*100,3)
  survci.u<-round(surv.out.anno$trends[[1]]$upCI*100,3)
  survci.l<-sprintf("%.3f",survci.l)  
  survci.u<-sprintf("%.3f",survci.u)
  survci.str<-paste("(",survci.l,", ",survci.u,")",sep="")
  sum[,8]<-survci.str
  
  ### surv 5-year
  sum[,9]<-surv.out.anno$trends[[2]]$estimate*100
  survci.l<-round(surv.out.anno$trends[[2]]$lowCI*100,3)
  survci.u<-round(surv.out.anno$trends[[2]]$upCI*100,3)
  survci.l<-sprintf("%.3f",survci.l)  
  survci.u<-sprintf("%.3f",survci.u)
  survci.str<-paste("(",survci.l,", ",survci.u,")",sep="")
  sum[,10]<-survci.str
  
  ### dying
  sum[,11]<-dying.out.anno$trends[[1]]$estimate*100
  ci.l<-round(dying.out.anno$trends[[1]]$lowCI*100,3)
  ci.u<-round(dying.out.anno$trends[[1]]$upCI*100,3)
  ci.l<-sprintf("%.3f",ci.l)  
  ci.u<-sprintf("%.3f",ci.u)
  ci.str<-paste("(",ci.l,", ",ci.u,")",sep="")
  sum[,12]<-ci.str
  
  out.sum[[i]]<-sum
  out.surv[[i]]<-surv.out.anno
  out.dying[[i]]<-dying.out.anno
}

out.sum.combine<-rbind(out.sum[[1]],out.sum[[2]],out.sum[[3]],out.sum[[4]],out.sum[[5]],out.sum[[6]])
save.f<-"P:/srab/surv/JPSurvival/Fanni_work/Paper/JPsurv_manuscript_FemaleBreast_Liver_trend_summary.txt.xls"
write.table(out.sum.combine,save.f,row.names=F,quote=F,sep="\t")

##################################################################################################
### NHL_CML
# paper example relative survival data
DICsample <- dictionary.overview("NHL_CML.dic")
head(DICsample)

input <- joinpoint.seerdata(seerfilename="NHL_CML", 
                            newvarnames=c("Year_of_diagnosis_1975"),
                            UseVarLabelsInData=c("Year_of_diagnosis_1975"))

subsetStr<-list()
subsetStr[[1]]<-"Site_recode_NHL_and_CML==0"
subsetStr[[2]]<-"Site_recode_NHL_and_CML==1"

titleStr<-list()
titleStr[[1]]<-"Non-Hodgkin Lymphoma"
titleStr[[2]]<-"Chronic Myeloid Leukemia"

out.sum<-list()
out.surv<-list()
out.dying<-list()
out.fit<-list()

yearvar<-"Year_of_diagnosis_1975"
obscumvar<-"Relative_Survival_Cum"
predcumvar<-"Predicted_Survival_Cum"
obsintvar<-"Relative_Survival_Interval"
predintvar<-"Predicted_ProbDeath_Int"
interval<-"Interval"
site<-c("Non-Hodgkin Lymphoma","Chronic Myeloid Leukemia")
stage<-c(NA,NA)

maxjp<-4

for(i in 1:length(subsetStr)){
  subset<-subsetStr[[i]]
  print(c(i,subset))
  fit <- joinpoint(input, subset,
                   year=yearvar, observedrelsurv="Relative_Survival_Cum",
                   model.form = NULL, maxnum.jp = maxjp, proj.year.num = 5)
  #trend.acs<-aapc.multiints(fit, type="AbsChgSur", int.select=interval.values)
  #Cancer site	Stage 	No. JPs. Final Model 	No. Alive	Start Year	End Year 	AACS(1) 	95% C.I.	AACS(5) 	95% C.I.	PCAPD	95% C.I.
  cnames<-c("Cancer site","Stage","No. JPs Final Model","No. Alive","Start Year","End Year",
            "AACS(1)","AACS(1) 95% C.I.","AACS(5)","AACS(5) 95% C.I.","PCAPD","PCAPD 95% C.I.")
  sum<-array("",dim=c((length(fit$jp)+1),length(cnames)))
  colnames(sum)<-cnames
  
  sum[1,1]<-site[i]
  sum[1,2]<-stage[i]
  if(is.null(fit$jp)){
    sum[1,3]<-0
  }else{ 
    sum[1,3]<-length(fit$jp)
  }
  nJP<-length(fit$jp)
  data.graph<-download.data(input,fit,nJP,yearvar,"graph",subset,interval,int.select=c(1,5,10))
  title.rch<-paste(titleStr[[i]]," \n Annual Probability of Dying of Cancer by Diagnosis Year",sep="")
  title.acs<-paste(titleStr[[i]]," \n Relative Survival by Diagnosis Year",sep="")
  if(length(fit$jp)<=3){
    dying.out.anno<-plot.dying.year.annotate(data.graph,fit,nJP,yearvar,obsintvar,predintvar,interval,annotation=1,topanno=1,trend=1,title=title.rch)
    surv.out.anno<-plot.surv.year.annotate(data.graph,fit,nJP,yearvar,obscumvar,predcumvar,interval,annotation=1,trend=1,title=title.acs)
    dying.plot<-dying.out.anno$plot_anno
    surv.plot<-surv.out.anno$plot_anno
  }else{
    dying.out.anno<-plot.dying.year.annotate(data.graph,fit,nJP,yearvar,obsintvar,predintvar,interval,annotation=0,topanno=0,trend=1,title=title.rch)
    surv.out.anno<-plot.surv.year.annotate(data.graph,fit,nJP,yearvar,obscumvar,predcumvar,interval,annotation=0,trend=1,title=title.acs)
    dying.plot<-dying.out.anno$plot_no_anno
    surv.plot<-surv.out.anno$plot_no_anno    
  }
  dyingplot.f<-paste("P:/srab/surv/JPSurvival/Fanni_work/Paper/plot/dying_plot_anno_",titleStr[[i]],".jpg",sep="")
  survplot.f<-paste("P:/srab/surv/JPSurvival/Fanni_work/Paper/plot/surv_plot_anno_",titleStr[[i]],".jpg",sep="")
  #ggsave(dyingplot.f,dying.plot, width=6.76, height=5.32)
  #ggsave(survplot.f,surv.plot, width=6.76, height=5.32)
  
  input.sub<-subset(input,eval(parse(text=subset)))
  sum[1,4]<-input.sub$Alive_at_Start[1]
  sum[,5]<-dying.out.anno$trends[[1]]$start.year
  sum[,6]<-dying.out.anno$trends[[1]]$end.year
  
  ### surv 1-year
  sum[,7]<-surv.out.anno$trends[[1]]$estimate*100
  survci.l<-round(surv.out.anno$trends[[1]]$lowCI*100,3)
  survci.u<-round(surv.out.anno$trends[[1]]$upCI*100,3)
  survci.l<-sprintf("%.3f",survci.l)  
  survci.u<-sprintf("%.3f",survci.u)
  survci.str<-paste("(",survci.l,", ",survci.u,")",sep="")
  sum[,8]<-survci.str
  
  ### surv 5-year
  sum[,9]<-surv.out.anno$trends[[2]]$estimate*100
  survci.l<-round(surv.out.anno$trends[[2]]$lowCI*100,3)
  survci.u<-round(surv.out.anno$trends[[2]]$upCI*100,3)
  survci.l<-sprintf("%.3f",survci.l)  
  survci.u<-sprintf("%.3f",survci.u)
  survci.str<-paste("(",survci.l,", ",survci.u,")",sep="")
  sum[,10]<-survci.str
  
  ### dying
  sum[,11]<-dying.out.anno$trends[[1]]$estimate*100
  ci.l<-round(dying.out.anno$trends[[1]]$lowCI*100,3)
  ci.u<-round(dying.out.anno$trends[[1]]$upCI*100,3)
  ci.l<-sprintf("%.3f",ci.l)  
  ci.u<-sprintf("%.3f",ci.u)
  ci.str<-paste("(",ci.l,", ",ci.u,")",sep="")
  sum[,12]<-ci.str
  
  
  out.fit[[i]]<-fit
  out.sum[[i]]<-sum
  #out.surv[[i]]<-surv.out.anno
  #out.dying[[i]]<-dying.out.anno
}

out.sum.combine<-rbind(out.sum[[1]],out.sum[[2]])
save.f<-"P:/srab/surv/JPSurvival/Fanni_work/Paper/JPsurv_manuscript_NHL_CML_overall_trend_summary.txt.xls"
write.table(out.sum.combine,save.f,row.names=F,quote=F,sep="\t")

### Model fit for CML using up to 2-year survival and 5 -year survival
subset.cml.fup2<-"Site_recode_NHL_and_CML==1 & Interval<=2"
fit.cml.fup2 <- joinpoint(input, subset.cml.fup2,
                 year=yearvar, observedrelsurv="Relative_Survival_Cum",
                 model.form = NULL, maxnum.jp = maxjp, proj.year.num = 5)

subset.cml.fup5<-"Site_recode_NHL_and_CML==1 & Interval<=5"
fit.cml.fup5 <- joinpoint(input, subset.cml.fup5,
                          year=yearvar, observedrelsurv="Relative_Survival_Cum",
                          model.form = NULL, maxnum.jp = maxjp, proj.year.num = 5)
fit.cml<-list()
fit.cml[[1]]<-fit.cml.fup2
fit.cml[[2]]<-fit.cml.fup5

### create table 2 for CML
cnames<-c("Cancer Site","Max Follow-up","No. of JPs","Loc. of JPs","BIC","AIC","Log Likelihood","Converged")
table.fit<-array("",dim=c((maxjp+1),length(cnames)))
colnames(table.fit)<-cnames
table.fit[1,1]<-"Chronic Myeloid Leukemia" ##cancer site

fup<-c("2 years","5 years")  ## max intervals from diagnosis (max follow-up)
table2.fup<-list()
for(j in 1:length(fup)){
  fit.cml.fup<-fit.cml[[j]]
  table.fit[1,2]<-fup[j]
  table.fit[,3]<-c(0:4)  ## no. of jps
  for(k in 1:(maxjp+1)){
    fit.k<-fit.cml.fup$FitList[[k]]
    loc.jps <- toString(fit.k$jp)
    table.fit[k,4]<-loc.jps  ## loc. of jps
    table.fit[k,5]<-fit.k$bic
    table.fit[k,6]<-fit.k$aic
    table.fit[k,7]<-fit.k$ll
    table.fit[k,8]<-fit.k$converged
  }
  table2.fup[[j]]<-table.fit
}
table2<-rbind(table2.fup[[1]],table2.fup[[2]])
save.f<-"P:/srab/surv/JPSurvival/Fanni_work/Paper/JPsurv_manuscript_NHL_CML_overall_modelfit_summary.txt.xls"
write.table(table2,save.f,row.names=F,quote=F,sep="\t")

#########################################################################################################
# plot.dying.year.full: a function that returns plots for Percent Change in the Annual Probability of
# Dying of Cancer by Diagnosis year using ggplot
# The annotation feature is available for nJP<=4 and the number of multiple intervals selected <=3
# arguments:
# fulldata - the full data returned by function download.data with downloadtype="full".
# plotdata - the graph data returned by function download.data with downloadtype="graph".
# fit - joinpoint object containing the model output.
# nJP - the number of joinpoints in the model.
# yearvar - the variable name for year of diagnosis used in argument 'year' of the function joinpoint.
# obsintvar - the variable name for observed interval survival. The default is "Relative_Survival_Interval"
#             for relative survival data. For cause-specific data, it needs to be changed accordingly.
# predintvar - the variable name for predicted interval survival. The default is "Predicted_ProbDeath_Int".
# interval - the variable name for year since diagnosis. The default is 'Interval'.
# annotation - the indicator for the annotation feature. The default is 0 (no annotation on the plot).Two plots
#              with and without annotation will be returned in a list when annotation=1.
# topanno - the indicator for showing the top curve annotation. The default is 1 (annotation for the top curve).
# trend - the indicator for returning the trend measure tables. The default is 0 (no trend tables returned).
# proj - the indicator for showing the trend projection from fulldata
#########################################################################################################

plot.dying.year.full<-function(data.full,plotdata,fit,nJP,yearvar,obsintvar="Relative_Survival_Interval",predintvar="Predicted_ProbDeath_Int",interval="Interval",annotation=0,topanno=1,trend=0,proj=0,title=NULL){
  title.rch<-"Percent Change in the Annual Probability of Dying of Cancer \n by Diagnosis Year"
  if(is.null(title)){
    title.rch<-"Annual Probability of Dying of Cancer by Diagnosis Year"
  }else{
    title.rch<-title
  }
  interval.values<-as.numeric(unique(plotdata[,"Interval"]))
  interval.labels<-paste((interval.values-1)," to ",interval.values," years since diagnosis",sep="")
  if(interval.values[1]==1){
    interval.labels[1]<-"0 to 1 year since diagnosis"
  }
  trends.out<-list()
  
  fulldata<-data.full[which(data.full[,interval] %in% interval.values),]
  yearint.plot<-paste(plotdata[,yearvar],"_",plotdata[,interval],sep="")
  yearint.full<-paste(fulldata[,yearvar],"_",fulldata[,interval],sep="")
  yearend<-plotdata[length(yearint.plot),yearvar]
  yearend.plot<-rep(yearend,length(interval.values))+1-interval.values
  yearint.end<-paste(yearend.plot,"_",interval.values,sep="")
  plotdata.add<-fulldata[which((yearint.full %in% yearint.end) | (!yearint.full %in% yearint.plot)),]
  plotdata.add[,"trend.label"]<-""
  #plotdata0<-plotdata
  #plotdata<-rbind(plotdata0,plotdata.add)
  
  if(trend==1 | annotation==1){ 
    trends.out <- aapc.multiints(fit$FitList[[nJP+1]], type="RelChgHaz", int.select=interval.values)
    if(length(interval.values)<=3 & nJP<=4){
      ### haz results are the same for all interval values
      jp.loc<-fit$FitList[[nJP+1]]$jp
      annot.strs<-paste(sprintf("%.1f",100*trends.out[[1]]$estimate),"%",sep="")
      x.values<-list()
      plotdata[,"trend.label"]<-""
      for(i in 1:length(interval.values)){
        int.i<-interval.values[i]
        haz.apc<-trends.out[[i]]
        plotdata.i<-plotdata[which(plotdata[,interval]==int.i),]
        if(max(plotdata.i[,yearvar])==max(haz.apc$end.year) | nJP==0){
          end.year<-haz.apc$end.year
          end.year[length(haz.apc$end.year)]<-max(plotdata.i[,yearvar])
          x.values.i<-1/2*(haz.apc$start.year+end.year)
        }
        if(max(plotdata.i[,yearvar])<max(haz.apc$end.year)){
          if(nJP==1){
            if(max(plotdata.i[,yearvar])>jp.loc){
              end.year<-haz.apc$end.year
              end.year[length(haz.apc$end.year)]<-max(plotdata.i[,yearvar])
              x.values.i<-1/2*(haz.apc$start.year+end.year)
            }else{
              haz.apc<-haz.apc[1:nJP,]
              x.values.i<-1/2*(haz.apc$start.year+haz.apc$end.year)
            }
          } ## nJP=1 end
          if(nJP==2){
            if(max(plotdata.i[,yearvar])>max(jp.loc)){
              end.year<-haz.apc$end.year
              end.year[length(haz.apc$end.year)]<-max(plotdata.i[,yearvar])
              x.values.i<-1/2*(haz.apc$start.year+end.year)
            } else if(max(plotdata.i[,yearvar])<=max(jp.loc) & max(plotdata.i[,yearvar])>min(jp.loc)){
              haz.apc<-haz.apc[1:nJP,]
              x.values.i<-1/2*(haz.apc$start.year+haz.apc$end.year)
            } else{
              haz.apc<-haz.apc[1:(nJP-1),]
              end.year<-haz.apc$end.year
              end.year[length(haz.apc$end.year)]<-max(plotdata.i[,yearvar])
              x.values.i<-1/2*(haz.apc$start.year+end.year)
            }
          } ## nJP=2 end
          if(nJP==3){
            if(max(plotdata.i[,yearvar])>max(jp.loc)){
              end.year<-haz.apc$end.year
              end.year[length(haz.apc$end.year)]<-max(plotdata.i[,yearvar])
              x.values.i<-1/2*(haz.apc$start.year+end.year)
            } else if(max(plotdata.i[,yearvar])<=max(jp.loc) & max(plotdata.i[,yearvar])>jp.loc[2]){
              haz.apc<-haz.apc[1:nJP,]
              end.year<-haz.apc$end.year
              end.year[length(haz.apc$end.year)]<-max(plotdata.i[,yearvar])
              x.values.i<-1/2*(haz.apc$start.year+end.year)
            } else if(max(plotdata.i[,yearvar])<=jp.loc[2] & max(plotdata.i[,yearvar])>=jp.loc[1]){
              haz.apc<-haz.apc[1:(nJP-1),]
              end.year<-haz.apc$end.year
              end.year[length(haz.apc$end.year)]<-max(plotdata.i[,yearvar])
              x.values.i<-1/2*(haz.apc$start.year+end.year)
            } else {
              haz.apc<-haz.apc[1:(nJP-2),]
              end.year<-haz.apc$end.year
              end.year[length(haz.apc$end.year)]<-max(plotdata.i[,yearvar])
              x.values.i<-1/2*(haz.apc$start.year+end.year)
            }
          } ## nJP=3 end
          if(nJP==4){
            if(max(plotdata.i[,yearvar])>max(jp.loc)){
              end.year<-haz.apc$end.year
              end.year[length(haz.apc$end.year)]<-max(plotdata.i[,yearvar])
              x.values.i<-1/2*(haz.apc$start.year+end.year)
            } else if(max(plotdata.i[,yearvar])<=max(jp.loc) & max(plotdata.i[,yearvar])>jp.loc[3]){
              haz.apc<-haz.apc[1:nJP,]
              end.year<-haz.apc$end.year
              end.year[length(haz.apc$end.year)]<-max(plotdata.i[,yearvar])
              x.values.i<-1/2*(haz.apc$start.year+end.year)
            } else if(max(plotdata.i[,yearvar])<=jp.loc[3] & max(plotdata.i[,yearvar])>jp.loc[2]){
              haz.apc<-haz.apc[1:(nJP-1),]
              end.year<-haz.apc$end.year
              end.year[length(haz.apc$end.year)]<-max(plotdata.i[,yearvar])
              x.values.i<-1/2*(haz.apc$start.year+end.year)
            } else if(max(plotdata.i[,yearvar])<=jp.loc[2] & max(plotdata.i[,yearvar])>=jp.loc[1]){
              haz.apc<-haz.apc[1:(nJP-2),]
              end.year<-haz.apc$end.year
              end.year[length(haz.apc$end.year)]<-max(plotdata.i[,yearvar])
              x.values.i<-1/2*(haz.apc$start.year+end.year)
            } else {
              haz.apc<-haz.apc[1:(nJP-3),]
              end.year<-haz.apc$end.year
              end.year[length(haz.apc$end.year)]<-max(plotdata.i[,yearvar])
              x.values.i<-1/2*(haz.apc$start.year+end.year)
            }
          } ## nJP=4 end
        } 
        x.values[[i]]<-round(x.values.i)
        plotdata$trend.label[which((plotdata[,yearvar] %in% x.values[[i]]) & plotdata[,interval]==interval.values[[i]])]<-annot.strs[1:length(x.values[[i]])]  ### 11/25
      }## interval.values iteration end
      
      plotdata0<-plotdata
      plotdata0[,"ind"]<-"obs"
      plotdata.add[,"ind"]<-"proj"
      plotdata<-rbind(plotdata0,plotdata.add)
      plotdata[,"obsintvar_1"]<-1-plotdata[,obsintvar]
      obsintvar_1minus<-"obsintvar_1"
      if(annotation==1){
        hues=seq(15,375,length=length(interval.values)+1)
        gg_color_hue<-hcl(h=hues,l=65,c=100)[1:length(interval.values)]
        if(topanno==1){
          ymeans.byint<-aggregate(plotdata[,"Predicted_ProbDeath_Int"], list(plotdata[,interval]), mean)
          topint<-which(ymeans.byint[,2]==max(ymeans.byint[,2]))
          plotdata$trend.label[which(plotdata[,interval]!=interval.values[topint])]<-""
          plotdata0$trend.label[which(plotdata0[,interval]!=interval.values[topint])]<-""
          
          plot.proj<-ggplot(data=plotdata, aes(x=plotdata[,yearvar], y=plotdata[,predintvar],group=as.factor(plotdata[,interval]), colour=as.factor(plotdata[,interval]))) + 
            geom_line(data=(subset(plotdata,ind=="obs")), aes(x=eval(parse(text=yearvar)), y=eval(parse(text=predintvar)),
                                                              group=as.factor(eval(parse(text=interval))), colour=as.factor(eval(parse(text=interval)))),linetype="solid")+
            geom_point(data=(subset(plotdata,ind=="obs")), aes(x=eval(parse(text=yearvar)), y=eval(parse(text=obsintvar_1minus)), 
                                                               group=as.factor(eval(parse(text=interval))), colour=as.factor(eval(parse(text=interval)))))+
            geom_line(data=(subset(plotdata,ind=="proj")), aes(x=eval(parse(text=yearvar)), y=eval(parse(text=predintvar)),
                                                               group=as.factor(eval(parse(text=interval))), colour=as.factor(eval(parse(text=interval)))),lty=2)+
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
            geom_text_repel(aes(label = trend.label),nudge_x=0.01,nudge_y=0.01,show.legend = FALSE)
          
          plot<-ggplot(data=plotdata0, aes(x=plotdata0[,yearvar], y=plotdata0[,predintvar],group=as.factor(plotdata0[,interval]), colour=as.factor(plotdata0[,interval]))) + 
            geom_line(data=plotdata0, aes(x=plotdata0[,yearvar], y=plotdata0[,predintvar],group=as.factor(plotdata0[,interval]), colour=as.factor(plotdata0[,interval])),linetype="solid")+
            geom_point(data=plotdata0, aes(x=plotdata0[,yearvar], y=1-plotdata0[,obsintvar], group=as.factor(plotdata0[,interval]), colour=as.factor(plotdata0[,interval])))+
            xlab('Year at diagnosis') + ylab('Annual Probability of Cancer Death (%)')+
            theme_bw()+
            ggtitle(title.rch) +
            coord_cartesian(ylim=c(0,1))+
            scale_y_continuous(breaks=seq(0,1,0.1),labels = scales::percent)+
            scale_x_continuous(breaks=seq(min(plotdata0[,yearvar],na.rm=T),max(plotdata0[,yearvar],na.rm=T),5))+
            scale_color_hue(labels = interval.labels)+
            theme(legend.position="bottom")+
            theme(legend.title=element_blank())+
            theme(plot.title = element_text(hjust = 0.5))+
            geom_text_repel(aes(label = trend.label),nudge_x=0.01,show.legend = FALSE)
        }else{
          plot.proj<-ggplot(data=plotdata, aes(x=plotdata[,yearvar], y=plotdata[,predintvar],group=as.factor(plotdata[,interval]), colour=as.factor(plotdata[,interval]))) + 
            geom_line(data=(subset(plotdata,ind=="obs")), aes(x=eval(parse(text=yearvar)), y=eval(parse(text=predintvar)),
                                                              group=as.factor(eval(parse(text=interval))), colour=as.factor(eval(parse(text=interval)))),linetype="solid")+
            geom_point(data=(subset(plotdata,ind=="obs")), aes(x=eval(parse(text=yearvar)), y=eval(parse(text=obsintvar_1minus)), 
                                                               group=as.factor(eval(parse(text=interval))), colour=as.factor(eval(parse(text=interval)))))+
            geom_line(data=(subset(plotdata,ind=="proj")), aes(x=eval(parse(text=yearvar)), y=eval(parse(text=predintvar)),
                                                               group=as.factor(eval(parse(text=interval))), colour=as.factor(eval(parse(text=interval)))),lty=2)+
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
            geom_text_repel(aes(label = trend.label),nudge_x=0.01,nudge_y=0.01,show.legend = FALSE)
          plot<-ggplot(data=plotdata0, aes(x=plotdata0[,yearvar], y=plotdata0[,predintvar],group=as.factor(plotdata0[,interval]), colour=as.factor(plotdata0[,interval]))) + 
            geom_line(data=plotdata0, aes(x=plotdata0[,yearvar], y=plotdata0[,predintvar],group=as.factor(plotdata0[,interval]), colour=as.factor(plotdata0[,interval])),linetype="solid")+
            geom_point(data=plotdata0, aes(x=plotdata0[,yearvar], y=1-plotdata0[,obsintvar], group=as.factor(plotdata0[,interval]), colour=as.factor(plotdata0[,interval])))+
            xlab('Year at diagnosis') + ylab('Annual Probability of Cancer Death (%)')+
            theme_bw()+
            ggtitle(title.rch) +
            coord_cartesian(ylim=c(0,1))+
            scale_y_continuous(breaks=seq(0,1,0.1),labels = scales::percent)+
            scale_x_continuous(breaks=seq(min(plotdata0[,yearvar],na.rm=T),max(plotdata0[,yearvar],na.rm=T),5))+
            scale_color_hue(labels = interval.labels)+
            theme(legend.position="bottom")+
            theme(legend.title=element_blank())+
            theme(plot.title = element_text(hjust = 0.5))+
            geom_text_repel(aes(label = trend.label),nudge_x=0.01,show.legend = FALSE)
        }
      }
    }
    if(length(interval.values)>3 | nJP>4){
      if(annotation==1){
        stop("This annotation feature is not available when the number of joinpoints>4 or the number of multiple intervals selected>3.")
      }
    }
  }
  plot0.proj<-ggplot(data=plotdata, aes(x=plotdata[,yearvar], y=plotdata[,predintvar],group=as.factor(plotdata[,interval]), colour=as.factor(plotdata[,interval]))) + 
    geom_line(data=(subset(plotdata,ind=="obs")), aes(x=eval(parse(text=yearvar)), y=eval(parse(text=predintvar)),
                                                      group=as.factor(eval(parse(text=interval))), colour=as.factor(eval(parse(text=interval)))),linetype="solid")+
    geom_point(data=(subset(plotdata,ind=="obs")), aes(x=eval(parse(text=yearvar)), y=eval(parse(text=obsintvar_1minus)), 
                                                       group=as.factor(eval(parse(text=interval))), colour=as.factor(eval(parse(text=interval)))))+
    geom_line(data=(subset(plotdata,ind=="proj")), aes(x=eval(parse(text=yearvar)), y=eval(parse(text=predintvar)),
                                                       group=as.factor(eval(parse(text=interval))), colour=as.factor(eval(parse(text=interval)))),lty=2)+
    xlab('Year at diagnosis') + ylab('Annual Probability of Cancer Death (%)')+
    theme_bw()+
    ggtitle(title.rch) +
    coord_cartesian(ylim=c(0,1))+
    scale_y_continuous(breaks=seq(0,1,0.1),labels = scales::percent)+
    scale_x_continuous(breaks=seq(min(plotdata[,yearvar],na.rm=T),max(plotdata[,yearvar],na.rm=T),5))+
    scale_color_hue(labels = interval.labels)+
    theme(legend.position="bottom")+
    theme(legend.title=element_blank())+
    theme(plot.title = element_text(hjust = 0.5))
  
  plot0<-ggplot(data=plotdata0, aes(x=plotdata0[,yearvar], y=plotdata0[,predintvar],group=as.factor(plotdata0[,interval]), colour=as.factor(plotdata0[,interval]))) + 
    geom_line(data=plotdata0, aes(x=plotdata0[,yearvar], y=plotdata0[,predintvar],group=as.factor(plotdata0[,interval]), colour=as.factor(plotdata0[,interval])),linetype="solid")+
    geom_point(data=plotdata0, aes(x=plotdata0[,yearvar], y=1-plotdata0[,obsintvar], group=as.factor(plotdata0[,interval]), colour=as.factor(plotdata0[,interval])))+
    xlab('Year at diagnosis') + ylab('Annual Probability of Cancer Death (%)')+
    theme_bw()+
    ggtitle(title.rch) +
    coord_cartesian(ylim=c(0,1))+
    scale_y_continuous(breaks=seq(0,1,0.1),labels = scales::percent)+
    scale_x_continuous(breaks=seq(min(plotdata0[,yearvar],na.rm=T),max(plotdata0[,yearvar],na.rm=T),5))+
    scale_color_hue(labels = interval.labels)+
    theme(legend.position="bottom")+
    theme(legend.title=element_blank())+
    theme(plot.title = element_text(hjust = 0.5))
  
  if(!annotation %in% c(0,1)){
    stop("The annotation indicator should be either 0 or 1.")
  }
  if(!topanno %in% c(0,1)){
    stop("The annotation indicator for the top curve only should be either 0 or 1.")
  }
  if(!proj %in% c(0,1)){
    stop("The projection indicator should be either 0 or 1.")
  }
  if(trend==0){
    if(annotation==0 | length(interval.values)>3 | nJP>4){
      if(proj==0){
        out<-c(plot0)
        names(out)<-c("plot_no_anno")
      }else if(proj==1){
        out<-c(plot0,plot0.proj)
        names(out)<-c("plot_no_anno","plot_no_anno_proj")
      }
    }else{
      if(proj==0){
        out<-list(plot,plot0)
        names(out)<-c("plot_anno","plot_no_anno")
      }else if(proj==1){
        out<-list(plot,plot.proj,plot0,plot0.proj)
        names(out)<-c("plot_anno","plot_anno_proj","plot_no_anno","plot_no_anno_proj")
      }
    }
  }else if(trend==1){
    if(annotation==0 | length(interval.values)>3 | nJP>4){
      if(proj==0){
        out<-list(trends.out,plot0)
        names(out)<-c("trends","plot_no_anno")
      }else if(proj==1){
        out<-list(trends.out,plot0)
        names(out)<-c("trends","plot_no_anno","plot_no_anno_proj")
      }
    }else{
      if(proj==0){
        out<-list(trends.out,plot,plot0)
        names(out)<-c("trends","plot_anno","plot_no_anno")
      }else if(proj==1){
        out<-list(trends.out,plot,plot.proj,plot0,plot0.proj)
        names(out)<-c("trends","plot_anno","plot_anno_proj","plot_no_anno","plot_no_anno_proj")
      }  
    }
  }else{
    stop("The trend indicator should be either 0 or 1.")
  }
  return(out)
}

##################################################################################################################################################################################################################
# plot.surv.year.full: a function that returns a plot with projection for Average Change in Cumulative Survival by diagnosis year
# using ggplot
# The annotation feature is available for nJP<=4 and the number of multiple intervals selected <=3
# arguments:
# plotdata - the graph data returned by function download.data with downloadtype="graph".
# fit - joinpoint object containing the model output.
# nJP - the number of joinpoints in the model.
# yearvar - the variable name for year of diagnosis used in argument 'year' of the function joinpoint.
# obscumvar - the variable name for observed relative cumulative survival. The default is "Relative_Survival_Cum"
#             for relative survival data. For cause-specific data, it needs to be changed accordingly.
# predcumvar - the variable name for predicted cumulative survival. The default is "Predicted_Survival_Cum".
# interval - the variable name for year since diagnosis. The default is 'Interval'.
# annotation - the indicator for the annotation feature. The default is 0 (no annotation on the plot).Two plots
#              with and without annotation will be returned in a list when annotation=1.
# trend - the indicator for returning the trend measure tables. The default is 0 (no trend tables returned).
# proj - the indicator for showing the trend projection from fulldata
#########################################################################################################

plot.surv.year.full<-function(data.full,plotdata,fit,nJP,yearvar,obscumvar="Relative_Survival_Cum",predcumvar="Predicted_Survival_Cum",interval="Interval",annotation=0,trend=0,proj=0,title=NULL){
  interval.values<-as.numeric(unique(plotdata[,interval]))
  trends.out<-list()
  fulldata<-data.full[which(data.full[,interval] %in% interval.values),]
  yearint.plot<-paste(plotdata[,yearvar],"_",plotdata[,interval],sep="")
  yearint.full<-paste(fulldata[,yearvar],"_",fulldata[,interval],sep="")
  yearend<-plotdata[length(yearint.plot),yearvar]
  yearend.plot<-rep(yearend,length(interval.values))+1-interval.values
  yearint.end<-paste(yearend.plot,"_",interval.values,sep="")
  plotdata.add<-fulldata[which((yearint.full %in% yearint.end) | (!yearint.full %in% yearint.plot)),]
  plotdata.add[,"trend.label"]<-""
  #plotdata0<-plotdata
  #plotdata<-rbind(plotdata0,plotdata.add)
  if(is.null(title)){
    if(grepl("Relative", obscumvar)==T){
      title.acs<-"Relative Survival by Diagnosis Year"
      ylabel<-"Relative Survival (%)"
      interval.labels<-paste(interval.values,"-year Relative Survival",sep="")
    }else if(grepl("Cause", obscumvar)==T){
      title.acs<-"Cause-Specific Survival by Diagnosis Year"
      ylabel<-"Cause-Specific Survival (%)"
      interval.labels<-paste(interval.values,"-year Cause-Specific Survival",sep="")
    }else{
      title.acs<-"Survival by Diagnosis Year"
      ylabel<-"Survival (%)"
      interval.labels<-paste(interval.values,"-year Survival",sep="")
    }
  }else{
    title.acs<-title
    if(grepl("Relative", obscumvar)==T){
      ylabel<-"Relative Survival (%)"
      interval.labels<-paste(interval.values,"-year Relative Survival",sep="")
    }else if(grepl("Cause", obscumvar)==T){
      ylabel<-"Cause-Specific Survival (%)"
      interval.labels<-paste(interval.values,"-year Cause-Specific Survival",sep="")
    }else{
      ylabel<-"Survival (%)"
      interval.labels<-paste(interval.values,"-year Survival",sep="")
    }
  }
  
  if(trend==1 | annotation==1){
    trends.out <- aapc.multiints(fit$FitList[[nJP+1]], type="AbsChgSur", int.select=interval.values)
    if(length(interval.values)<=3 & nJP<=4){
      ### haz results are the same for all interval values
      jp.loc<-fit$FitList[[nJP+1]]$jp
      x.values<-list()
      annot.strs<-list()
      plotdata[,"trend.label"]<-""
      for(i in 1:length(interval.values)){
        int.i<-interval.values[i]
        cs.aaac<-trends.out[[i]]
        annot.strs[[i]]<-paste(sprintf("%.1f",100*cs.aaac$estimate),sep="")
        plotdata.i<-plotdata[which(plotdata[,interval]==int.i),]
        if(max(plotdata.i[,yearvar])==max(cs.aaac$end.year) | nJP==0){
          end.year<-cs.aaac$end.year
          end.year[length(cs.aaac$end.year)]<-max(plotdata.i[,yearvar])
          x.values.i<-1/2*(cs.aaac$start.year+end.year)
        }
        if(max(plotdata.i[,yearvar])<max(cs.aaac$end.year)){
          if(nJP==1){
            if(max(plotdata.i[,yearvar])>jp.loc){
              end.year<-cs.aaac$end.year
              end.year[length(cs.aaac$end.year)]<-max(plotdata.i[,yearvar])
              x.values.i<-1/2*(cs.aaac$start.year+end.year)
            }else{
              cs.aaac<-cs.aaac[1:nJP,]
              x.values.i<-1/2*(cs.aaac$start.year+cs.aaac$end.year)
            }
          } ## nJP=1 end
          if(nJP==2){
            if(max(plotdata.i[,yearvar])>max(jp.loc)){
              end.year<-cs.aaac$end.year
              end.year[length(cs.aaac$end.year)]<-max(plotdata.i[,yearvar])
              x.values.i<-1/2*(cs.aaac$start.year+end.year)
            } else if(max(plotdata.i[,yearvar])<=max(jp.loc) & max(plotdata.i[,yearvar])>min(jp.loc)){
              cs.aaac<-cs.aaac[1:nJP,]
              x.values.i<-1/2*(cs.aaac$start.year+cs.aaac$end.year)
            } else{
              cs.aaac<-cs.aaac[1:(nJP-1),]
              end.year<-cs.aaac$end.year
              end.year[length(cs.aaac$end.year)]<-max(plotdata.i[,yearvar])
              x.values.i<-1/2*(cs.aaac$start.year+end.year)
            }
          } ## nJP=2 end
          if(nJP==3){
            if(max(plotdata.i[,yearvar])>max(jp.loc)){
              end.year<-cs.aaac$end.year
              end.year[length(cs.aaac$end.year)]<-max(plotdata.i[,yearvar])
              x.values.i<-1/2*(cs.aaac$start.year+end.year)
            } else if(max(plotdata.i[,yearvar])<=max(jp.loc) & max(plotdata.i[,yearvar])>jp.loc[2]){
              cs.aaac<-cs.aaac[1:nJP,]
              end.year<-cs.aaac$end.year
              end.year[length(cs.aaac$end.year)]<-max(plotdata.i[,yearvar])
              x.values.i<-1/2*(cs.aaac$start.year+end.year)
            } else if(max(plotdata.i[,yearvar])<=jp.loc[2] & max(plotdata.i[,yearvar])>=jp.loc[1]){
              cs.aaac<-cs.aaac[1:(nJP-1),]
              end.year<-cs.aaac$end.year
              end.year[length(cs.aaac$end.year)]<-max(plotdata.i[,yearvar])
              x.values.i<-1/2*(cs.aaac$start.year+end.year)
            } else {
              cs.aaac<-cs.aaac[1:(nJP-2),]
              end.year<-cs.aaac$end.year
              end.year[length(cs.aaac$end.year)]<-max(plotdata.i[,yearvar])
              x.values.i<-1/2*(cs.aaac$start.year+end.year)
            }
          } ## nJP=3 end
          if(nJP==4){
            if(max(plotdata.i[,yearvar])>max(jp.loc)){
              end.year<-cs.aaac$end.year
              end.year[length(cs.aaac$end.year)]<-max(plotdata.i[,yearvar])
              x.values.i<-1/2*(cs.aaac$start.year+end.year)
            } else if(max(plotdata.i[,yearvar])<=max(jp.loc) & max(plotdata.i[,yearvar])>jp.loc[3]){
              cs.aaac<-cs.aaac[1:nJP,]
              end.year<-cs.aaac$end.year
              end.year[length(cs.aaac$end.year)]<-max(plotdata.i[,yearvar])
              x.values.i<-1/2*(cs.aaac$start.year+end.year)
            } else if(max(plotdata.i[,yearvar])<=max(jp.loc) & max(plotdata.i[,yearvar])>jp.loc[2]){
              cs.aaac<-cs.aaac[1:(nJP-1),]
              end.year<-cs.aaac$end.year
              end.year[length(cs.aaac$end.year)]<-max(plotdata.i[,yearvar])
              x.values.i<-1/2*(cs.aaac$start.year+end.year)
            } else if(max(plotdata.i[,yearvar])<=jp.loc[2] & max(plotdata.i[,yearvar])>=jp.loc[1]){
              cs.aaac<-cs.aaac[1:(nJP-2),]
              end.year<-cs.aaac$end.year
              end.year[length(cs.aaac$end.year)]<-max(plotdata.i[,yearvar])
              x.values.i<-1/2*(cs.aaac$start.year+end.year)
            } else {
              cs.aaac<-cs.aaac[1:(nJP-3),]
              end.year<-cs.aaac$end.year
              end.year[length(cs.aaac$end.year)]<-max(plotdata.i[,yearvar])
              x.values.i<-1/2*(cs.aaac$start.year+end.year)
            }
          } ## nJP=4 end
        } 
        x.values[[i]]<-round(x.values.i)
        plotdata$trend.label[which((plotdata[,yearvar] %in% x.values[[i]]) & plotdata[,interval]==interval.values[[i]])]<-annot.strs[[i]][1:length(x.values[[i]])]
      }## interval.values iteration end
      
      plotdata0<-plotdata
      plotdata0[,"ind"]<-"obs"
      plotdata.add[,"ind"]<-"proj"
      plotdata<-rbind(plotdata0,plotdata.add)
      
      if(annotation==1){
        hues=seq(15,375,length=length(interval.values)+1)
        gg_color_hue<-hcl(h=hues,l=65,c=100)[1:length(interval.values)]
        plot.proj<-ggplot(data=plotdata, aes(x=plotdata[,yearvar], y=plotdata[,predcumvar],group=as.factor(plotdata[,interval]), colour=as.factor(plotdata[,interval]))) + 
          geom_line(data=(subset(plotdata,ind=="obs")), aes(x=eval(parse(text=yearvar)), y=eval(parse(text=predcumvar)),
                                                            group=as.factor(eval(parse(text=interval))), colour=as.factor(eval(parse(text=interval)))),linetype="solid")+
          geom_point(data=(subset(plotdata,ind=="obs")), aes(x=eval(parse(text=yearvar)), y=eval(parse(text=obscumvar)), 
                                                             group=as.factor(eval(parse(text=interval))), colour=as.factor(eval(parse(text=interval)))))+
          geom_line(data=(subset(plotdata,ind=="proj")), aes(x=eval(parse(text=yearvar)), y=eval(parse(text=predcumvar)),
                                                             group=as.factor(eval(parse(text=interval))), colour=as.factor(eval(parse(text=interval)))),lty=2)+
          
          xlab('Year at diagnosis') + ylab(ylabel)+
          theme_bw()+
          ggtitle(title.acs) +
          coord_cartesian(ylim=c(0,1.02))+ ### could do manual adjustment of ylim range for better locations of annotations 
          scale_y_continuous(breaks=seq(0,1,0.1),labels = scales::percent)+
          scale_x_continuous(breaks=seq(min(plotdata[,yearvar],na.rm=T),max(plotdata[,yearvar],na.rm=T),5))+
          scale_color_hue(labels = interval.labels)+
          theme(legend.position="bottom")+
          theme(legend.title=element_blank())+
          theme(plot.title = element_text(hjust = 0.5))+
          geom_text_repel(aes(label = trend.label),nudge_x=0.02,nudge_y=0.01,show.legend = FALSE)
        plot<-ggplot(data=plotdata0, aes(x=plotdata0[,yearvar], y=plotdata0[,predcumvar],group=as.factor(plotdata0[,interval]), colour=as.factor(plotdata0[,interval]))) + 
          geom_line(data=plotdata0, aes(x=plotdata0[,yearvar], y=plotdata0[,predcumvar],group=as.factor(plotdata0[,interval]), colour=as.factor(plotdata0[,interval])),linetype="solid")+
          geom_point(data=plotdata0, aes(x=plotdata0[,yearvar], y=plotdata0[,obscumvar], group=as.factor(plotdata0[,interval]), colour=as.factor(plotdata0[,interval])))+
          xlab('Year at diagnosis') + ylab(ylabel)+
          theme_bw()+
          ggtitle(title.acs) +
          coord_cartesian(ylim=c(0,1))+
          scale_y_continuous(breaks=seq(0,1,0.1),labels = scales::percent)+
          scale_x_continuous(breaks=seq(min(plotdata[,yearvar],na.rm=T),max(plotdata[,yearvar],na.rm=T),5))+
          scale_color_hue(labels = interval.labels)+
          theme(legend.position="bottom")+
          theme(legend.title=element_blank())+
          theme(plot.title = element_text(hjust = 0.5))+
          geom_text_repel(aes(label = trend.label),nudge_x=0.02,nudge_y=0.01,show.legend = FALSE)
      }
    }
    if(length(interval.values)>3 | nJP>4){
      stop("This annotation feature is not available when the number of joinpoints>4 or the number of multiple intervals selected>3.")
    }
  }
  plot0.proj<-ggplot(data=plotdata, aes(x=plotdata[,yearvar], y=plotdata[,predcumvar],group=as.factor(plotdata[,interval]), colour=as.factor(plotdata[,interval]))) + 
    geom_line(data=(subset(plotdata,ind=="obs")), aes(x=eval(parse(text=yearvar)), y=eval(parse(text=predcumvar)),
                                                      group=as.factor(eval(parse(text=interval))), colour=as.factor(eval(parse(text=interval)))),linetype="solid")+
    geom_point(data=(subset(plotdata,ind=="obs")), aes(x=eval(parse(text=yearvar)), y=eval(parse(text=obscumvar)), 
                                                       group=as.factor(eval(parse(text=interval))), colour=as.factor(eval(parse(text=interval)))))+
    geom_line(data=(subset(plotdata,ind=="proj")), aes(x=eval(parse(text=yearvar)), y=eval(parse(text=predcumvar)),
                                                       group=as.factor(eval(parse(text=interval))), colour=as.factor(eval(parse(text=interval)))),lty=2)+
    xlab('Year at diagnosis') + ylab(ylabel)+
    theme_bw()+
    ggtitle(title.acs) +
    coord_cartesian(ylim=c(0,1))+
    scale_y_continuous(breaks=seq(0,1,0.1),labels = scales::percent)+
    scale_x_continuous(breaks=seq(min(plotdata[,yearvar],na.rm=T),max(plotdata[,yearvar],na.rm=T),5))+
    scale_color_hue(labels = interval.labels)+
    theme(legend.position="bottom")+
    theme(legend.title=element_blank())+
    theme(plot.title = element_text(hjust = 0.5))
  plot0<-ggplot(data=plotdata0, aes(x=plotdata0[,yearvar], y=plotdata0[,predcumvar],group=as.factor(plotdata0[,interval]), colour=as.factor(plotdata0[,interval]))) + 
    geom_line(data=plotdata0, aes(x=plotdata0[,yearvar], y=plotdata0[,predcumvar],group=as.factor(plotdata0[,interval]), colour=as.factor(plotdata0[,interval])),linetype="solid")+
    geom_point(data=plotdata0, aes(x=plotdata0[,yearvar], y=plotdata0[,obscumvar], group=as.factor(plotdata0[,interval]), colour=as.factor(plotdata0[,interval])))+
    xlab('Year at diagnosis') + ylab(ylabel)+
    theme_bw()+
    ggtitle(title.acs) +
    coord_cartesian(ylim=c(0,1))+
    scale_y_continuous(breaks=seq(0,1,0.1),labels = scales::percent)+
    scale_x_continuous(breaks=seq(min(plotdata[,yearvar],na.rm=T),max(plotdata[,yearvar],na.rm=T),5))+
    scale_color_hue(labels = interval.labels)+
    theme(legend.position="bottom")+
    theme(legend.title=element_blank())+
    theme(plot.title = element_text(hjust = 0.5))
  if(!annotation %in% c(0,1)){
    stop("The annotation indicator should be either 0 or 1.")
  }
  if(trend==0){
    if(annotation==0 | length(interval.values)>3 | nJP>4){
      if(proj==0){
        out<-c(plot0)
        names(out)<-c("plot_no_anno")
      }else if(proj==1){
        out<-c(plot0,plot0.proj)
        names(out)<-c("plot_no_anno","plot_no_anno_proj")
      }
    }else{
      if(proj==0){
        out<-list(plot,plot0)
        names(out)<-c("plot_anno","plot_no_anno")
      }else if(proj==1){
        out<-list(plot,plot.proj,plot0,plot0.proj)
        names(out)<-c("plot_anno","plot_anno_proj","plot_no_anno","plot_no_anno_proj")
      }
    }
  }else if(trend==1){
    if(annotation==0 | length(interval.values)>3 | nJP>4){
      if(proj==0){
        out<-list(trends.out,plot0)
        names(out)<-c("trends","plot_no_anno")
      }else if(proj==1){
        out<-list(trends.out,plot0)
        names(out)<-c("trends","plot_no_anno","plot_no_anno_proj")
      }
    }else{
      if(proj==0){
        out<-list(trends.out,plot,plot0)
        names(out)<-c("trends","plot_anno","plot_no_anno")
      }else if(proj==1){
        out<-list(trends.out,plot,plot.proj,plot0,plot0.proj)
        names(out)<-c("trends","plot_anno","plot_anno_proj","plot_no_anno","plot_no_anno_proj")
      }  
    }
  }else{
    stop("The trend indicator should be either 0 or 1.")
  }
  return(out)
}

#####################################################################################################
### test functions plot.dying.year.full and plot.surv.year.full (2 examples)
### 1. test plot.dying.year.full and plot.surv.year.full using cohort female breast - regional
DICsample <- dictionary.overview("FemaleBreast_Liver.dic")
head(DICsample)

input <- joinpoint.seerdata(seerfilename="FemaleBreast_Liver", 
                            newvarnames=c("Year_of_diagnosis_1975"),
                            UseVarLabelsInData=c("Year_of_diagnosis_1975"))
head(input)
yearvar<-"Year_of_diagnosis_1975"
obscumvar<-"Relative_Survival_Cum"
predcumvar<-"Predicted_Survival_Cum"
obsintvar<-"Relative_Survival_Interval"
predintvar<-"Predicted_ProbDeath_Int"
interval<-"Interval"
# cohort: female breast - regional
subset<-"Stage_AllLRD==1 & Liver_FBreast ==1"
maxjp<-3
fit.bcr<-joinpoint(input, subset,
                   year=yearvar, observedrelsurv="Relative_Survival_Cum",
                   model.form = NULL, maxnum.jp = maxjp, proj.year.num = 5)

# best model njp=3
nJP<-3

data.graph<-download.data(input,fit.bcr,nJP,yearvar,"graph",subset,interval="Interval",int.select=c(1,5,10))
data.full<-download.data(input,fit.bcr,nJP,yearvar,"full",subset,interval="Interval")
save.f<-"P:/srab/surv/JPSurvival/Fanni_work/Paper/plotdata_Female_Breast_Regional_graph.csv"
write.csv(data.graph,save.f,row.names=F,quote=F)
save.f<-"P:/srab/surv/JPSurvival/Fanni_work/Paper/plotdata_Female_Breast_Regional_full.csv"
write.csv(data.full,save.f,row.names=F,quote=F)

# plots with projections
title<-"Female Breast - Regional"
title.rch<-paste(title," \n Annual Probability of Dying of Cancer by Diagnosis Year",sep="")
title.acs<-paste(title," \n Relative Survival by Diagnosis Year",sep="")
out.dash.anno.dying<-plot.dying.year.full(data.full,data.graph,fit.bcr,nJP,yearvar,obsintvar,predintvar,interval,annotation=1,topanno=1,trend=1,proj=1,title=title.rch)
out.dash.anno.surv<-plot.surv.year.full(data.full,data.graph,fit.bcr,nJP,yearvar,obscumvar,predcumvar,interval,annotation=1,trend=1,proj=1,title=title.acs)
plot.proj.dying.bcr<-out.dash.anno.dying$plot_anno_proj
plot.proj.surv.bcr<-out.dash.anno.surv$plot_anno_proj

dyingplot.f<-paste("P:/srab/surv/JPSurvival/Fanni_work/Paper/plot/test_dying_plot_anno_projection_",title,".jpg",sep="")
survplot.f<-paste("P:/srab/surv/JPSurvival/Fanni_work/Paper/plot/test_surv_plot_anno_projection_",title,".jpg",sep="")
ggsave(dyingplot.f,plot.proj.dying.bcr, width=6.76, height=5.32)
ggsave(survplot.f,plot.proj.surv.bcr, width=6.76, height=5.32)


### 2. test plot.dying.year.full and plot.surv.year.full using NHL cohort
DICsample <- dictionary.overview("NHL_CML.dic")
head(DICsample)

input <- joinpoint.seerdata(seerfilename="NHL_CML", 
                            newvarnames=c("Year_of_diagnosis_1975"),
                            UseVarLabelsInData=c("Year_of_diagnosis_1975"))
yearvar<-"Year_of_diagnosis_1975"
obscumvar<-"Relative_Survival_Cum"
predcumvar<-"Predicted_Survival_Cum"
obsintvar<-"Relative_Survival_Interval"
predintvar<-"Predicted_ProbDeath_Int"
interval<-"Interval"
subset<-"Site_recode_NHL_and_CML==0"
maxjp<-4
fit.nhl<-joinpoint(input, subset,
                   year=yearvar, observedrelsurv="Relative_Survival_Cum",
                   model.form = NULL, maxnum.jp = maxjp, proj.year.num = 5)

# best model njp=4
nJP<-4
data.graph<-download.data(input,fit.nhl,nJP,yearvar,"graph",subset,interval="Interval",int.select=c(1,5,10))
data.full<-download.data(input,fit.nhl,nJP,yearvar,"full",subset,interval="Interval")

save.f<-"P:/srab/surv/JPSurvival/Fanni_work/Paper/plotdata_Non_Hodgkin_Lymphoma_graph.csv"
write.csv(data.graph,save.f,row.names=F,quote=F)
save.f<-"P:/srab/surv/JPSurvival/Fanni_work/Paper/plotdata_Non_Hodgkin_Lymphoma_full.csv"
write.csv(data.full,save.f,row.names=F,quote=F)

# plots with projections
title<-"Non-Hodgkin Lymphoma"
title.rch<-paste(title," \n Annual Probability of Dying of Cancer by Diagnosis Year",sep="")
title.acs<-paste(title," \n Relative Survival by Diagnosis Year",sep="")
out.dash.anno.dying<-plot.dying.year.full(data.full,data.graph,fit.nhl,nJP,yearvar,obsintvar,predintvar,interval,annotation=1,topanno=1,trend=1,proj=1,title=title.rch)
out.dash.anno.surv<-plot.surv.year.full(data.full,data.graph,fit.nhl,nJP,yearvar,obscumvar,predcumvar,interval,annotation=1,trend=1,proj=1,title=title.acs)
plot.proj.dying.nhl<-out.dash.anno.dying$plot_anno_proj
plot.proj.surv.nhl<-out.dash.anno.surv$plot_anno_proj
dyingplot.f<-paste("P:/srab/surv/JPSurvival/Fanni_work/Paper/plot/test_dying_plot_anno_projection_",title,".jpg",sep="")
survplot.f<-paste("P:/srab/surv/JPSurvival/Fanni_work/Paper/plot/test_surv_plot_anno_projection_",title,".jpg",sep="")
ggsave(dyingplot.f,plot.proj.dying.nhl, width=6.76, height=5.32)
ggsave(survplot.f,plot.proj.surv.nhl, width=6.76, height=5.32)
