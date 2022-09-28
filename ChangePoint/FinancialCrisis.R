library(readr)
library(lubridate)
library(TSA)
library(rugarch)
library(e1071)
library(WeightSVM)

# read data
data_org <- read_csv("./../data/BDI.csv")
colnames(data_org)<-c("Date","Value")
data_org$Date=ymd(data_org$Date)
require(lubridate)
require(dplyr)


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
             gamma=fgamma_arma,#越大越複雜
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
              gamma=fgamma_garch,#越大越複雜
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
               gamma=fgamma_arma,#越大越複雜
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
                gamma=fgamma_garch,#越大越複雜
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


#################################
# Function used to predict by SVR ARMA GARCH model
predict.svr.arma.grach=function(model,newdata){
  all_sequence=c(model$arma_data$sequence,newdata)
  tmp_data=data.frame(
    sequence=all_sequence
  )
  tmp_ar=as.data.frame(matrix(NA,nrow(tmp_data),model$arma_order[1]))
  for(i in 1:model$arma_order[1]){
    tmp_ar[,i]=c(rep(NA,i),all_sequence[1:(length(all_sequence)-i)])
  }
  colnames(tmp_ar)=paste0("AR_",c(1:model$arma_order[1]))
  tmp_data=cbind(tmp_data,tmp_ar)
  tmp=as.data.frame(matrix(NA,length(newdata),model$arma_order[2]))
  tmp2=as.data.frame(as.matrix(model$arma_data[,c((ncol(model$arma_data)-model$arma_order[2]+1):ncol(model$arma_data))]))
  colnames(tmp)=colnames(tmp2)=paste0("MA_",c(1:model$arma_order[2]))
  tmp_MA=rbind(tmp2,tmp)
  
  tmp_data=cbind(tmp_data,tmp_MA)
  tmp_data$predict=c(model$arma_model$fitted,rep(NA,length(newdata)))
  tmp_data$MA_0=tmp_data$sequence-tmp_data$predict
  all_i=which(is.na(tmp_data$predict))
  all_arma=which(is.na(tmp_data$predict))
  for(i in all_i){
    for(tmpi in 1:model$arma_order[2]){
      tmp_data[i,model$arma_order[1]+1+tmpi]=tmp_data$MA_0[i-tmpi]
    }
    tmp_data$predict[i]=predict(model$arma_model,tmp_data[i,])
    tmp_data$MA_0[i]=tmp_data$sequence[i]-tmp_data$predict[i]
  }
  garch_sequence=c(tmp_data$MA_0^2)
  garch_data=data.frame(
    sequence=garch_sequence
  )
  garch_ar=as.data.frame(matrix(NA,nrow(garch_data),model$garch_order[1]))
  for(i in 1:model$garch_order[1]){
    garch_ar[,i]=c(rep(NA,i),garch_sequence[1:(length(garch_sequence)-i)])
  }
  colnames(garch_ar)=paste0("AR_",c(1:model$garch_order[1]))
  garch_data=cbind(garch_data,garch_ar)
  tmp=as.data.frame(matrix(NA,length(newdata),model$garch_order[2]))
  tmp2=as.matrix(model$garch_data[,c((ncol(model$garch_data)-model$garch_order[2]+1):ncol(model$garch_data))])
  colnames(tmp2)=paste0("MA_",c(1:model$garch_order[2]))
  colnames(tmp)=colnames(tmp2)
  tmp_MA=rbind(tmp2,tmp)
  nrow(garch_data)-nrow(tmp_MA)
  tmp_MA=rbind(matrix(NA,nrow(garch_data)-nrow(tmp_MA),ncol(tmp_MA)),as.matrix(tmp_MA))
  colnames(tmp_MA)=colnames(tmp)
  garch_data=cbind(garch_data,tmp_MA)
  garch_data$predict=c(model$garch_model$fitted,rep(NA,nrow(garch_data)-length(model$garch_model$fitted)))
  garch_data$MA_0=garch_data$sequence-garch_data$predict
  all_i=which(is.na(garch_data$predict))
  i=all_i[1]
  for(i in all_i){
    for(tmpi in 1:model$garch_order[2]){
      garch_data[i,model$garch_order[1]+1+tmpi]=garch_data$MA_0[i-tmpi]
    }
    garch_data$predict[i]=predict(model$garch_model,garch_data[i,])
    garch_data$MA_0[i]=garch_data$sequence[i]-garch_data$predict[i]
  }
  output=list(mean=tmp_data[all_arma,],var=garch_data[all_arma,])
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

#######################
tmp_f=function(x){x*(1-x)}

T_mean=function(true=NA,predict=NA,residuals=NA,K=3.48){
  if(length(true)==1){true=predict+residuals}
  if(length(predict)==1){predict=true-residuals}
  if(length(residuals)==1){residuals=true-predict}
  tau1=sqrt(mean((predict^2)*(residuals^2)))
  tmp1=mean(predict*residuals)
  
  n=length(true)
  LS=function(K){
    sqrt(((abs(sum(((predict)*(residuals))[1:K])-K*tmp1))^2)/(n*(tau1^2)))
  }
  tmp=(sapply(1:n,LS))
  if(sum(is.na(tmp))==length(tmp)){
    tmp=rep(0,length(tmp))
  }
  tmp_f=function(x){x*(1-x)}
  true_sd=sqrt(sapply(c(1:length(tmp))/length(tmp), tmp_f))
  output_CC=data.frame(
    T_mean=tmp,
    sd=true_sd,
    UCL=K*true_sd
  )
  output_CC[-nrow(output_CC),]
}
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

data
min(data$date[data$CUSUM_s>=1])
min(data$date[data$CUSUM>=1])



data[c(which(data$CUSUM_s>=1)-1,which(data$CUSUM_s>=1),which(data$CUSUM_s>=1)+1) %>% unique() %>% sort,] 

#par(mfrow=c(2,1))
plot(data$date,data$true,type="l")
points(data$date[data$CUSUM_s>=1],data$true[data$CUSUM_s>=1],col="red")

plot(data$date,data$true,type="l")
points(data$date[data$CUSUM>=1],data$true[data$CUSUM>=1],col="red")
#plot(data$date,runmean(data$CUSUM_s,k=10))
require(caTools)

#plot(data$date,runmean(data$CUSUM_s,k=10))



## 爆發流動性危機
tmp_time = ymd(20070218)# 發現有改變點的時間
#tmp_time = ymd(20070919)# 發現有改變點的時間
i=which(data$date==tmp_time)
tmp=T_mean(true=data$true[1:i],
           predict = data$predict[1:i],
           residuals = data$residuals[1:i])

L=length(tmp$T_mean)+1
true_sd=sqrt(sapply(c(1:L)/L, tmp_f))[-L]
UCL_d=3.48*true_sd
l=as.integer(L*0.22)
UCL_h=c(rep(1.358,l),UCL_d[(l+1):(L-1)])
plot(tmp$T_mean,type="l")
points(UCL_h,type="l",col="red")
# 判斷改變點 位置
data$date[min(which(tmp$T_mean>UCL_h))]#"2007-08-09"


detech_time=lapply(data$date[data$CUSUM_s>=1], function(tmp_time){
  ## 爆發流動性危機
  #tmp_time = ymd(20070218)# 發現有改變點的時間
  #tmp_time = ymd(20070919)# 發現有改變點的時間
  i=which(data$date==tmp_time)
  tmp=T_mean(true=data$true[1:i],
             predict = data$predict[1:i],
             residuals = data$residuals[1:i])
  
  L=length(tmp$T_mean)+1
  true_sd=sqrt(sapply(c(1:L)/L, tmp_f))[-L]
  UCL_d=3.48*true_sd
  l=as.integer(L*0.22)
  UCL_h=c(rep(1.358,l),UCL_d[(l+1):(L-1)])
  #plot(tmp$T_mean,type="l")
  #points(UCL_h,type="l",col="red")
  # 判斷改變點 位置
  data$date[min(which(tmp$T_mean>UCL_h))]#"2007-08-09"
  
})  %>% do.call(what="c") #%>% unique %>% sort

data.frame(data$date[data$CUSUM_s>=1],detech_time) %>% View

plot(data$date,data$true,type="l")
points(data$date[data$date%in%detech_time],data$true[data$date%in%detech_time],col="red")



plot(data_org$Date,data_org$Value,type="l")
points(data_org$Date[data_org$Date%in%detech_time],data_org$Value[data_org$Date%in%detech_time],col="red")
