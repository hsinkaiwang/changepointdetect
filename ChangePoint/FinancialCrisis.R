library(readr)
library(lubridate)
library(TSA)
library(rugarch)
library(e1071)
library(WeightSVM)
require(dplyr)
source("Change_point_function.R")

# read data
data_org <- read_csv("./../data/BDI.csv")
colnames(data_org)<-c("Date","Value")
data_org$Date=ymd(data_org$Date)


# Extract the data which between 2003 to 2010
data_org=data_org[which(data_org$Date>ymd(20030101)),]
data_org=data_org[which(data_org$Date<ymd(20100101)),]
data_org=data_org %>% arrange(Date) %>% filter(!(Date%in%ymd(20070101,20080101,20090101)))

Sys.setlocale("LC_TIME", "C")
date_range=range(data_org$Date)
data=data.frame(Date=seq(date_range[1],date_range[2],1)) %>% filter(!weekdays(Date)%in%c("Saturday","Sunday")) %>% 
  left_join(data_org)
data$Value=imputeTS::na.interpolation(data$Value)

# Calculate the log return
data$return=c(NA,diff(log(data$Value)))
data$RW=c(NA,(data$return[-nrow(data)]))
data=data[-1,]

# Partition data into train and test sets
train=data[-which(data$Date>ymd(20060101)),]
test=data[which(data$Date>ymd(20060101)),]

# Function used in building SVR ARMA GARCH
svr.arma.grach=function(
    sequence, # time series sequence
    arma.order, # ARMA order
    epsilon_arma, # ARMA svm parameter epsilon
    cost_arma, # ARMA svm parameter cost
    fgamma_arma, # ARMA svm parameter gamma
    garch.order, # GARCH order
    epsilon_garch, # GARCH  svm parameter epsilon
    cost_garch, # GARCH  svm parameter cost
    fgamma_garch, # GARCH  svm parameter gamma
    stop_rate 
){
  # Original sequence
  svm_data=data.frame(
    sequence=sequence
  )
  # Create AR terms
  tmp_ar=as.data.frame(matrix(NA,nrow(svm_data),arma.order[1]))
  for(i in 1:arma.order[1]){
    tmp_ar[,i]=c(rep(NA,i),sequence[1:(length(sequence)-i)])
  }
  colnames(tmp_ar)=paste0("AR_",c(1:arma.order[1]))
  svm_data=cbind(svm_data,tmp_ar)
  svm_data_train=na.omit(svm_data)
  # Build ARMA using SVR
  arma0=wsvm(sequence~.,data = svm_data_train, weight = rep(1,nrow(svm_data_train)),
             kernel ="radial",
             type = "eps-regression",
             cost = cost_arma,
             gamma=fgamma_arma,#???????????????
             epsilon=epsilon_arma)
  # Using the residual of the ARMA to calculate the MA terms of ARMA
  tmp_ma=as.data.frame(matrix(NA,nrow(svm_data_train),arma.order[2]))
  for(i in 1:arma.order[2]){
    tmp_ma[,i]=c(rep(0,i),arma0$residuals[1:(length(arma0$residuals)-i)])
  }
  colnames(tmp_ma)=paste0("MA_",c(1:arma.order[2]))
  svm_data_train_MA=cbind(svm_data_train,tmp_ma)
  # Initial residuals
  R_0=arma0$residuals
  # Initial a2 of the GARCH
  a2=R_0^2
  garch.order=c(max(garch.order),garch.order[2])
  garch_data=data.frame(
    sequence=a2
  )
  # Create the AR terms of GARCH
  garch_ar=as.data.frame(matrix(NA,nrow(garch_data),garch.order[1]))
  for(i in 1:garch.order[1]){
    garch_ar[,i]=c(rep(NA,i),a2[1:(length(a2)-i)])
  }
  colnames(garch_ar)=paste0("AR_",c(1:garch.order[1]))
  garch_data=cbind(garch_data,garch_ar)
  garch_data=na.omit(garch_data)
  # Build GARCH using SVR
  grach0=wsvm(sequence~.,data = garch_data, weight = rep(1,nrow(garch_data)),
              kernel ="radial",
              type = "eps-regression",
              cost = cost_garch,
              gamma=fgamma_garch,#???????????????
              epsilon=epsilon_garch)
  # Using the residual of the GARCH to calculate the MA terms of GARCH
  garch_ma=as.data.frame(matrix(NA,nrow(garch_data),garch.order[2]))
  for(i in 1:garch.order[2]){
    garch_ma[,i]=c(rep(0,i),grach0$residuals[1:(length(grach0$residuals)-i)])
  }
  colnames(garch_ma)=paste0("MA_",c(1:garch.order[2]))
  garch_data_ma=cbind(garch_data,garch_ma)
  
  # Initial prediction of a2 of GARCH
  hat_a2=c(rep(mean(predict(grach0)),max(garch.order)),predict(grach0))
  hat_a2[which(hat_a2<0)]=min(hat_a2[which(hat_a2>0)])
  V0=grach0$residuals
  Keep=T
  R=0
  
  # Iterative the ARMA and GARCH until the stopping criteria match
  while (Keep == T) {
    R=1+R
    arma1=wsvm(sequence~.,data = svm_data_train_MA, weight = 1/sqrt(hat_a2),
               kernel ="radial",
               type = "eps-regression",
               cost = cost_arma*mean(sqrt(hat_a2)),
               gamma=fgamma_arma,#???????????????
               epsilon=epsilon_arma)
    tmp_ma=as.data.frame(matrix(NA,nrow(svm_data_train_MA),arma.order[2]))
    for(i in 1:arma.order[2]){
      tmp_ma[,i]=c(rep(0,i),arma1$residuals[1:(length(arma1$residuals)-i)])
    }
    colnames(tmp_ma)=paste0("MA_",c(1:arma.order[2]))
    svm_data_train_MA=cbind(svm_data_train,tmp_ma)
    R_1=arma1$residuals
    
    a2=R_1^2
    garch_data=data.frame(
      sequence=a2
    )
    garch_ar=as.data.frame(matrix(NA,nrow(garch_data),garch.order[1]))
    for(i in 1:garch.order[1]){
      garch_ar[,i]=c(rep(NA,i),a2[1:(length(a2)-i)])
    }
    colnames(garch_ar)=paste0("AR_",c(1:garch.order[1]))
    garch_data=cbind(garch_data,garch_ar)
    tmp=matrix(NA,max(garch.order),garch.order[2])
    colnames(tmp)=colnames(garch_ma)
    garch_data=cbind(garch_data,rbind(tmp,garch_ma))
    garch_data=na.omit(garch_data)
    grach1=wsvm(sequence~.,data = garch_data, weight = rep(1,nrow(garch_data)),
                kernel ="radial",
                type = "eps-regression",
                cost = cost_garch,
                gamma=fgamma_garch,#???????????????
                epsilon=epsilon_garch)
    garch_ma=as.data.frame(matrix(NA,nrow(garch_data),garch.order[2]))
    for(i in 1:garch.order[2]){
      garch_ma[,i]=c(rep(0,i),grach1$residuals[1:(length(grach1$residuals)-i)])
    }
    colnames(garch_ma)=paste0("MA_",c(1:garch.order[2]))
    hat_a2=c(rep(mean(predict(grach1)),max(garch.order)),predict(grach1))
    hat_a2[which(hat_a2<0)]=min(hat_a2[which(hat_a2>0)])
    V1=grach1$residuals
    print(paste(R,":","The rate of change of ARMA parameter-",mean(abs((R_1-R_0)/R_0))))
    print(paste(R,":","The rate of change of GARCH parameter-",mean(abs((V1-V0)/V0))))
    print("#########################")
    tmp=(mean(abs((R_1-R_0)/R_0))>stop_rate)|(mean(abs((V1-V0)/V0))>stop_rate)
    Keep=(tmp&(R<20))
    R_0=R_1
    V0=V1
  }
  output=list(arma_model=arma1,
              arma_data=svm_data_train_MA,
              arma_order=arma.order,
              garch_model=grach1,
              garch_data=garch_data,
              garch_order=garch.order
  )
  output
}



##############
m=svr.arma.grach(
  sequence=train$return,
  arma.order = c(2,2),
  epsilon_arma = 0.1,
  cost_arma = 1,
  fgamma_arma = 0.01,
  garch.order = c(1,1),
  epsilon_garch = 0.08,
  cost_garch = 1,
  fgamma_garch = 0.075,
  stop_rate = 0.01
)
pre=predict.svr.arma.grach(newdata=test$return,model=m)


######################
data=data.frame(
  date=test$Date,
  true=pre$mean$sequence,
  predict=pre$mean$sequence,
  residuals=pre$mean$sequence-pre$mean$predict
)
###################
data$date=ymd(data$date)
data$CUSUM=0
data$CUSUM_s=0
for(i in 251:nrow(data)){
  tmp=T_mean(true=data$true[1:i],
             predict = data$predict[1:i],
             residuals = data$residuals[1:i])
  # Original CUMSUM
  data$CUSUM[i]=as.numeric(sum(tmp$T_mean>1.358)>0)
  L=length(tmp$T_mean)+1
  true_sd=sqrt(sapply(c(1:L)/L, tmp_f))[-L]
  UCL_d=3.48*true_sd
  # CUMSUM using modified boundary 
  l=as.integer(L*0.22)
  UCL_h=c(rep(1.358,l),UCL_d[(l+1):(L-1)])
  UCL_max_scale=2.413
  data$CUSUM_s[i]=as.numeric(sum(tmp$T_mean>UCL_h)>0)
}



require(ggplot2)
require(gridExtra)
require(grid)
tmp_cumsum_s_date = data$date[data$CUSUM_s>=1] 
test_cumsum2=test[test$Date%in%tmp_cumsum_s_date,]
two_date = ymd(c(20070518,20080110))


p1 <- ggplot()+
  geom_line(test,mapping=aes(x=Date,y=return))+
  geom_point(test_cumsum2,mapping=aes(x=Date,y=return),color="red")+
  geom_vline(aes(xintercept = two_date),color="red",linetype ="dashed")+
  geom_text(mapping=aes(x = two_date+days(60),y=rep(quantile(test$return,0.999),2),label=as.character(two_date)),color="red")+
  xlab("Date")+
  ylab("Return")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

p2 <- ggplot()+
  geom_line(test,mapping=aes(x=Date,y=Value))+
  geom_point(test_cumsum2,mapping=aes(x=Date,y=Value),color="red")+
  geom_vline(aes(xintercept = two_date),color="red",linetype ="dashed")+
  geom_text(mapping=aes(x = two_date+days(60),y=rep(quantile(test$Value,0.999),2),label=as.character(two_date)),color="red")+
  xlab("Date")+
  ylab("BDI Value")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

grid.arrange(p1,p2,
             top = textGrob(expression("Times of detect change point by"~T[s]),
                            gp=gpar(fontsize=20,font=1)))


# Collect the first time when the test is rejected
detech_time=lapply(data$date[data$CUSUM_s>=1], function(tmp_time){
  i=which(data$date==tmp_time)
  tmp=T_mean(true=data$true[1:i],
             predict = data$predict[1:i],
             residuals = data$residuals[1:i])
  
  L=length(tmp$T_mean)+1
  true_sd=sqrt(sapply(c(1:L)/L, tmp_f))[-L]
  UCL_d=3.48*true_sd
  l=as.integer(L*0.22)
  UCL_h=c(rep(1.358,l),UCL_d[(l+1):(L-1)])
  # ??????????????? ??????
  data$date[min(which(tmp$T_mean>UCL_h))]#"2007-08-09"
  
  
})  %>% do.call(what="c") 

# The first column of the following table shows the time when the change point is detected.
# The second column of the following table shows the first time when the test is rejected
data.frame(detech_time=data$date[data$CUSUM_s>=1],reject_time=detech_time)

