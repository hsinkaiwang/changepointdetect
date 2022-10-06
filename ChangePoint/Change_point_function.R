
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