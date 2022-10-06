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



# Extract the data after year 2015
data_org=data_org[which(data_org$Date>ymd(20150101)),]
data_org=data_org %>% arrange(Date) %>% filter(!(Date%in%ymd(20160101,20170101,20180101,
                                                             20190101,20200101,20210101)))

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
train=data[-which(data$Date>ymd(20190101)),]
test=data[which(data$Date>ymd(20190101)),]

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
two_date = ymd(c(20200103,20200331))


p1 <- ggplot()+
  geom_line(test,mapping=aes(x=Date,y=return))+
  geom_point(test_cumsum2,mapping=aes(x=Date,y=return),color="red")+
  geom_vline(aes(xintercept = two_date),color="red",linetype ="dashed")+
  geom_text(mapping=aes(x = two_date+days(30),y=rep(quantile(test$return,0.999),2),label=as.character(two_date)),color="red")+
  xlab("Date")+
  ylab("Return")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

p2 <- ggplot()+
  geom_line(test,mapping=aes(x=Date,y=Value))+
  geom_point(test_cumsum2,mapping=aes(x=Date,y=Value),color="red")+
  geom_vline(aes(xintercept = two_date),color="red",linetype ="dashed")+
  geom_text(mapping=aes(x = two_date+days(30),y=rep(quantile(test$Value,0.999),2),label=as.character(two_date)),color="red")+
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
  # 判斷改變點 位置
  data$date[min(which(tmp$T_mean>UCL_h))]#"2007-08-09"
  
  
})  %>% do.call(what="c") 

# The first column of the following table shows the time when the change point is detected.
# The second column of the following table shows the first time when the test is rejected
data.frame(detech_time=data$date[data$CUSUM_s>=1],reject_time=detech_time)

