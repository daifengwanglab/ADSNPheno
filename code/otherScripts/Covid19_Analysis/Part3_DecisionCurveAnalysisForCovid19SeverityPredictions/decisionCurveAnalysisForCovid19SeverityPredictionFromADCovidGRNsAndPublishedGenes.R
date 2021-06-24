
# Please note that the code below on Decision Curve Analysis is from the Memorial
# Sloan Kettering Cancer Center:
dca <- function(data, outcome, predictors, xstart=0.01, xstop=0.99, xby=0.01,
                ymin=-0.05, probability=NULL, harm=NULL,graph=TRUE, intervention=FALSE,
                interventionper=100, smooth=FALSE,loess.span=0.10, context = "") {

  print(paste0(":) Please note that this code on DCA (Decision Curve Analysis) is directly from the Memorial Sloan Kettering Cancer Center (MSKCC): https://www.mskcc.org/departments/epidemiology-biostatistics/biostatistics/decision-curve-analysis"))
  # LOADING REQUIRED LIBRARIES
  require(stats)

  # data MUST BE A DATA FRAME
  if (class(data)!="data.frame") {
    stop("Input data must be class data.frame")
  }

  #ONLY KEEPING COMPLETE CASES
  data=data[complete.cases(data[append(outcome,predictors)]),append(outcome,predictors)]

  # outcome MUST BE CODED AS 0 AND 1
  if (max(data[[outcome]])>1 | min(data[[outcome]])<0) {
    stop("outcome cannot be less than 0 or greater than 1")
  }
  # xstart IS BETWEEN 0 AND 1
  if (xstart<0 | xstart>1) {
    stop("xstart must lie between 0 and 1")
  }

  # xstop IS BETWEEN 0 AND 1
  if (xstop<0 | xstop>1) {
    stop("xstop must lie between 0 and 1")
  }

  # xby IS BETWEEN 0 AND 1
  if (xby<=0 | xby>=1) {
    stop("xby must lie between 0 and 1")
  }

  # xstart IS BEFORE xstop
  if (xstart>=xstop) {
    stop("xstop must be larger than xstart")
  }

  #STORING THE NUMBER OF PREDICTORS SPECIFIED
  pred.n=length(predictors)

  #IF probability SPECIFIED ENSURING THAT EACH PREDICTOR IS INDICATED AS A YES OR NO
  if (length(probability)>0 & pred.n!=length(probability)) {
    stop("Number of probabilities specified must be the same as the number of predictors being checked.")
  }

  #IF harm SPECIFIED ENSURING THAT EACH PREDICTOR HAS A SPECIFIED HARM
  if (length(harm)>0 & pred.n!=length(harm)) {
    stop("Number of harms specified must be the same as the number of predictors being checked.")
  }

  #INITIALIZING DEFAULT VALUES FOR PROBABILITES AND HARMS IF NOT SPECIFIED
  if (length(harm)==0) {
    harm=rep(0,pred.n)
  }
  if (length(probability)==0) {
    probability=rep(TRUE,pred.n)
  }


  #CHECKING THAT EACH probability ELEMENT IS EQUAL TO YES OR NO,
  #AND CHECKING THAT PROBABILITIES ARE BETWEEN 0 and 1
  #IF NOT A PROB THEN CONVERTING WITH A LOGISTIC REGRESSION
  for(m in 1:pred.n) {
    if (probability[m]!=TRUE & probability[m]!=FALSE) {
      stop("Each element of probability vector must be TRUE or FALSE")
    }
    if (probability[m]==TRUE & (max(data[predictors[m]])>1 | min(data[predictors[m]])<0)) {
      stop(paste(predictors[m],"must be between 0 and 1 OR sepcified as a non-probability in the probability option",sep=" "))
    }
    if(probability[m]==FALSE) {
      model=NULL
      pred=NULL
      model=glm(data.matrix(data[outcome]) ~ data.matrix(data[predictors[m]]), family=binomial("logit"))
      pred=data.frame(model$fitted.values)
      pred=data.frame(pred)
      names(pred)=predictors[m]
      data=cbind(data[names(data)!=predictors[m]],pred)
      print(paste(predictors[m],"converted to a probability with logistic regression. Due to linearity assumption, miscalibration may occur.",sep=" "))
    }
  }

  # THE PREDICTOR NAMES CANNOT BE EQUAL TO all OR none.
  if (length(predictors[predictors=="all" | predictors=="none"])) {
    stop("Prediction names cannot be equal to all or none.")
  }

  #########  CALCULATING NET BENEFIT   #########
  N=dim(data)[1]
  event.rate=colMeans(data[outcome])

  # CREATING DATAFRAME THAT IS ONE LINE PER THRESHOLD PER all AND none STRATEGY
  nb=data.frame(seq(from=xstart, to=xstop, by=xby))
  names(nb)="threshold"
  interv=nb

  nb["all"]=event.rate - (1-event.rate)*nb$threshold/(1-nb$threshold)
  nb["none"]=0

  # CYCLING THROUGH EACH PREDICTOR AND CALCULATING NET BENEFIT
  for(m in 1:pred.n){
    for(t in 1:length(nb$threshold)){
      # COUNTING TRUE POSITIVES AT EACH THRESHOLD
      tp=mean(data[data[[predictors[m]]]>=nb$threshold[t],outcome])*sum(data[[predictors[m]]]>=nb$threshold[t])
      # COUNTING FALSE POSITIVES AT EACH THRESHOLD
      fp=(1-mean(data[data[[predictors[m]]]>=nb$threshold[t],outcome]))*sum(data[[predictors[m]]]>=nb$threshold[t])
      #setting TP and FP to 0 if no observations meet threshold prob.
      if (sum(data[[predictors[m]]]>=nb$threshold[t])==0) {
        tp=0
        fp=0
      }

      # CALCULATING NET BENEFIT
      nb[t,predictors[m]]=tp/N - fp/N*(nb$threshold[t]/(1-nb$threshold[t])) - harm[m]
    }
    interv[predictors[m]]=(nb[predictors[m]] - nb["all"])*interventionper/(interv$threshold/(1-interv$threshold))
  }

  # CYCLING THROUGH EACH PREDICTOR AND SMOOTH NET BENEFIT AND INTERVENTIONS AVOIDED
  for(m in 1:pred.n) {
    if (smooth==TRUE){
      lws=loess(data.matrix(nb[!is.na(nb[[predictors[m]]]),predictors[m]]) ~ data.matrix(nb[!is.na(nb[[predictors[m]]]),"threshold"]),span=loess.span)
      nb[!is.na(nb[[predictors[m]]]),paste(predictors[m],"_sm",sep="")]=lws$fitted

      lws=loess(data.matrix(interv[!is.na(nb[[predictors[m]]]),predictors[m]]) ~ data.matrix(interv[!is.na(nb[[predictors[m]]]),"threshold"]),span=loess.span)
      interv[!is.na(nb[[predictors[m]]]),paste(predictors[m],"_sm",sep="")]=lws$fitted
    }
  }

  # PLOTTING GRAPH IF REQUESTED
  if (graph==TRUE) {
    require(graphics)

    # PLOTTING INTERVENTIONS AVOIDED IF REQUESTED
    if(intervention==TRUE) {
      # initialize the legend label, color, and width using the standard specs of the none and all lines
      legendlabel <- NULL
      legendcolor <- NULL
      legendwidth <- NULL
      legendpattern <- NULL

      #getting maximum number of avoided interventions
      ymax=max(interv[predictors],na.rm = TRUE)

      #INITIALIZING EMPTY PLOT WITH LABELS
      plot(x=nb$threshold, y=nb$all, type="n" ,xlim=c(xstart, xstop), ylim=c(ymin, ymax), xlab="Threshold probability", ylab=paste("Net reduction in interventions per",interventionper,"patients"))

      #PLOTTING INTERVENTIONS AVOIDED FOR EACH PREDICTOR
      for(m in 1:pred.n) {
        if (smooth==TRUE){
          lines(interv$threshold,data.matrix(interv[paste(predictors[m],"_sm",sep="")]),col=m,lty=2)
        } else {
          lines(interv$threshold,data.matrix(interv[predictors[m]]),col=m,lty=2)
        }

        # adding each model to the legend
        legendlabel <- c(legendlabel, predictors[m])
        legendcolor <- c(legendcolor, m)
        legendwidth <- c(legendwidth, 1)
        legendpattern <- c(legendpattern, 2)
      }
    } else {
      # PLOTTING NET BENEFIT IF REQUESTED

      # initialize the legend label, color, and width using the standard specs of the none and all lines
      legendlabel <- c("None", "All")
      legendcolor <- c(17, 8)
      legendwidth <- c(2, 2)
      legendpattern <- c(1, 1)

      #getting maximum net benefit
      ymax=max(nb[names(nb)!="threshold"],na.rm = TRUE)

      # inializing new benfit plot with treat all option
      plot(x=nb$threshold, y=nb$all, type="l", col=8, lwd=2 ,xlim=c(xstart, xstop), ylim=c(ymin, ymax),
           xlab="Threshold probability", ylab="Net benefit", main = context)
      # adding treat none option
      lines(x=nb$threshold, y=nb$none,lwd=2)
      #PLOTTING net benefit FOR EACH PREDICTOR
      for(m in 1:pred.n) {
        if (smooth==TRUE){
          lines(nb$threshold,data.matrix(nb[paste(predictors[m],"_sm",sep="")]),col=m,lty=2)
        } else {
          lines(nb$threshold,data.matrix(nb[predictors[m]]),col=m,lty=2)
        }
        # adding each model to the legend
        legendlabel <- c(legendlabel, predictors[m])
        legendcolor <- c(legendcolor, m)
        legendwidth <- c(legendwidth, 1)
        legendpattern <- c(legendpattern, 2)
      }
    }
    # then add the legend
    legend("topright", legendlabel, cex=0.8, col=legendcolor, lwd=legendwidth, lty=legendpattern)

  }

  #RETURNING RESULTS
  results=list()
  results$N=N
  results$predictors=data.frame(cbind(predictors,harm,probability))
  names(results$predictors)=c("predictor","harm.applied","probability")
  results$interventions.avoided.per=interventionper
  results$net.benefit=nb
  results$interventions.avoided=interv

  return(results)

}

# https://www.mskcc.org/departments/epidemiology-biostatistics/biostatistics/decision-curve-analysis

library(MASS)
data.set <- birthwt
model = glm(low ~ age + lwt, family=binomial(link="logit"), data=data.set)
data.set$predlow = predict(model, type="response")
dca(data=data.set, outcome="low", predictors="predlow", smooth="TRUE", xstop=0.50)

#
df = read.csv("F://updatedPredDF_Covid_SVC_LTL_4Fold.csv", header = TRUE)
df$Region[1]
dfCol = grep("CovidICU.1.", colnames(df))
colnames(df)[grep("CovidICU.1.", colnames(df))] = paste0(colnames(df)[grep("CovidICU.1.", colnames(df))], "_", df$Region[1])

ltlNB = dca(data=df, outcome="actual", predictors=colnames(df)[dfCol], smooth="FALSE", xstop=1, context = df$Region[1])

write.csv(ltlNB$net.benefit, "F://ltlNB_svc_new.csv")



df = read.csv("F://updatedPredDF_Covid_SVC_DLPFC_4Fold.csv", header = TRUE)
#df = read.csv("D://predDF_Covid_LTL.csv", header = TRUE)
df$Region[1]
dfCol = grep("CovidICU.1.", colnames(df))
colnames(df)[grep("CovidICU.1.", colnames(df))] = paste0(colnames(df)[grep("CovidICU.1.", colnames(df))], "_", df$Region[1])

dlpfcNB = dca(data=df, outcome="actual", predictors=colnames(df)[dfCol], smooth="FALSE", xstop=1, context = df$Region[1])
write.csv(dlpfcNB$net.benefit, "F://dlpfcNB_svc_new.csv")




df = read.csv("F://updatedPredDF_Covid_SVC_DLPFC_4Fold.csv", header = TRUE)
#df = read.csv("D://predDF_Covid_LTL.csv", header = TRUE)
df$Region[1]
dfCol = grep("CovidICU.1.", colnames(df))
colnames(df)[grep("CovidICU.1.", colnames(df))] = paste0(colnames(df)[grep("CovidICU.1.", colnames(df))], "_", df$Region[1])

dlpfcNB = dca(data=df, outcome="actual", predictors=colnames(df)[dfCol], smooth="FALSE", xstop=1, context = df$Region[1])
write.csv(dlpfcNB$net.benefit, "F://all3NB_svc_new.csv")
# df$Region[1]
#
# dlpfcNB = dca(data=df, outcome="actual", predictors="CovidICU.1.", smooth="FALSE", xstop=1, context = df$Region[1])
# write.csv(dlpfcNB$net.benefit, "D://dlpfcNB.csv")



df = read.csv("F://updatedPredDF_Covid_SVC_published_4Fold.csv", header = TRUE)
df$Region[1]
dfCol = grep("CovidICU.1.", colnames(df))
colnames(df)[grep("CovidICU.1.", colnames(df))] = paste0(colnames(df)[grep("CovidICU.1.", colnames(df))], "_", df$Region[1])

dlpfcNB = dca(data=df, outcome="actual", predictors=colnames(df)[dfCol], smooth="FALSE", xstop=1, context = df$Region[1])
write.csv(dlpfcNB$net.benefit, "F://bmcNB_svc_new.csv")

##########################################################################
df = read.csv("F://updatedPredDF_Covid_SVC_Hippocampus_4Fold.csv", header = TRUE)
df$Region[1]
dfCol = grep("CovidICU.1.", colnames(df))
colnames(df)[grep("CovidICU.1.", colnames(df))] = paste0(colnames(df)[grep("CovidICU.1.", colnames(df))], "_", df$Region[1])

dlpfcNB = dca(data=df, outcome="actual", predictors=colnames(df)[dfCol], smooth="FALSE", xstop=1, context = df$Region[1])
write.csv(dlpfcNB$net.benefit, "F://hippNB_svc_new.csv")


####################################################
covidADNetBenefitSVC
#nb = read.csv("D://covidNB.csv", header = TRUE)
#nb = read.csv("F://covidADNetBenefitSVC.csv", header = TRUE)
nb = read.csv("F://covidADNetBenefitSVCUpdated.csv", header = TRUE)

colnames(nb) = c("threshold", "netbenefit", "info")#, "old")
unique(nb$info)
head(nb)
library("ggplot2")
p<-ggplot(nb, aes(x=threshold, y=netbenefit, group=info)) +
  geom_line(aes(color=info))+
  geom_point(aes(color=info))
p = p + ylim(0, 0.5)
p = p + ggtitle("Decision Curve Analysis for Covid-19 Severity Prediction for\nCovid-19 Positive Patients: ICU (1) versus Non-ICU (0)") +
  xlab("Probability Threshold (for Predicting Severe and Sending Covid-19 Patients to the ICU)") + ylab("Net Benefit")
p
p = p + scale_color_manual(values=c("darkorchid2", "blue2",  "lightseagreen", "violet", "orange", "turquoise", "indianred3"))

p
pdf("F://organizedAlzheimers//dcaCovidAdSVC.pdf", width = 10, height = 8)
p
dev.off()

#p = p + scale_color_manual(values=c("coral1", "darkorchid2", "blue2", "lightseagreen", "violet", "darkseagreen3", "orange"))
p

#p = p + scale_color_manual(values=c("darkorchid2", "blue2", "lightseagreen", "violet", "orange"))


dca(data=df, outcome="actual", predictors="CovidICU.1.", smooth="FALSE", xstop=1, context = df$Region[1])

########################
par(mfrow = c(2, 2))

df = read.csv("D://predDF_Covid_LTL.csv", header = TRUE)
df$Region[1]

dca(data=df, outcome="actual", predictors="CovidICU.1.", smooth="FALSE", xstop=1, context = df$Region[1])



df = read.csv("D://predDF_Covid_DLPFC.csv", header = TRUE)
df$Region[1]

dca(data=df, outcome="actual", predictors="CovidICU.1.", smooth="FALSE", xstop=1, context = df$Region[1])




df = read.csv("D://predDF_Covid_bmc.csv", header = TRUE)
df$Region[1]

dca(data=df, outcome="actual", predictors="CovidICU.1.", smooth="FALSE", xstop=1, context = df$Region[1])


df = read.csv("D://predDF_Covid_hipp.csv", header = TRUE)
df$Region[1]

dca(data=df, outcome="actual", predictors="CovidICU.1.", smooth="FALSE", xstop=1, context = df$Region[1])



dev.off()



#nb = read.csv("D://optimalCovidADSVCNetBenefits.csv", header = TRUE)
nb = read.csv("F://optimalCovidADSVCNetBenefitsNEW.csv", header = TRUE)

unique(df$info)
colnames(nb) = c("threshold", "netbenefit", "info")
unique(nb$info)
head(nb)
library("ggplot2")
# p<-ggplot(nb, aes(x=threshold, y=netbenefit, group=info)) +
p<-ggplot(nb, aes(x=threshold, y=netbenefit)) +

  geom_point(aes(color=info))
p = p + ylim(0, 0.5)
p = p + ggtitle("Optimal Models for Various Probability Thresholds Based on Decision Curve Analysis\nfor Severity Prediction for Covid-19 Positive Patients: ICU (1) versus Non-ICU (0)") +
  xlab("Probability Threshold (for Predicting Severe and Sending Covid-19 Patients to the ICU)") + ylab("Net Benefit")
p
p = p + scale_color_manual(values=c("darkorchid2", "blue2",  "lightseagreen", "violet", "orange", "turquoise", "indianred3"))


pdf("F://organizedAlzheimers//dcaCovidAdSVC_Optimal.pdf", width = 10, height = 8)
p
dev.off()

p


dca(data=df, outcome="actual", predictors="CovidICU.1.", smooth="FALSE", xstop=1, context = df$Region[1])

########################


#
nb = read.csv("F://optimalCovidADSVCNetBenefitsNEW.csv", header = TRUE)


colnames(nb) = c("threshold", "netbenefit", "info")
unique(nb$info)
head(nb)
library("ggplot2")
p<-ggplot(nb, aes(x=threshold, y=netbenefit)) + geom_point(aes(color=info))
p = p + ylim(0, 0.5)
p = p + ggtitle("Optimal Models for Various Probability Thresholds Based on Decision Curve Analysis\nfor Severity Prediction for Covid-19 Positive Patients: ICU (1) versus Non-ICU (0)") +
  xlab("Probability Threshold (for Predicting Severe and Sending Covid-19 Patients to the ICU)") + ylab("Net Benefit")
p

p = p + scale_color_manual(values=c("blue2",  "lightseagreen", "violet")) #"turquoise", "indianred3"
p
pdf("F://organizedAlzheimers//dcaCovidAdSVC_Optimal.pdf", width = 10, height = 8)
p
dev.off()


##########################

## Gain in optimality



nb = read.csv("F://netBenefitOptimalAnalysisVsPublishedUpdated1.csv",
              header = TRUE)

colnames(nb) = c("threshold", "netbenefit", "info")
unique(nb$info)
head(nb)
library("ggplot2")
# p<-ggplot(nb, aes(x=threshold, y=netbenefit, group=info)) +
p<-ggplot(nb, aes(x=threshold, y=netbenefit)) +

  #geom_line(aes(color=info))+
  geom_point(aes(color=info))
p = p + ylim(0, 0.5)
p = p + ggtitle("Gain in Net Benefit from AD-Covid GRN Optimal SVC Linear Kernel Models\nRelative to SVC Linear Kernel Model from Published Covid-19 Genes List\nfor Various Probability Thresholds Based on Decision Curve Analysis\nfor Severity Prediction for Covid-19 Positive Patients: ICU (1) versus Non-ICU (0)") +
  xlab("Probability Threshold (for Predicting Severe and Sending Covid-19 Patients to the ICU)") + ylab("Net Benefit")
p
#p = p + scale_color_manual(values=c("darkorchid2", "blue2",  "red",
# lightseagreen","orange", "violet",  "turquoise", "indianred3"))
p
p = p + scale_color_manual(values=c("blue2", "tomato",  "lightseagreen", "violet", "orange")) #"turquoise", "indianred3"
p
pdf("F://organizedAlzheimers//dcaCovidAdSVC_OptimalComparison.pdf", width = 10, height = 9)
p
dev.off()
