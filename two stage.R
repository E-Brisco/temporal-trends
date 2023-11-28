TwoStageCMA<-function(measure,yi,vi,mC,mT,sdC,sdT,nC,nT,level,delta,order){
  if(missing(measure)){
    stop("measure is not specified")
  }
  if(measure != "RAW" && measure != "SMD"){
    stop("measure type is not valid")
  }else if(measure=="SMD"){
    k<-length(yi)
    if (missing(yi)||missing(vi)||missing(nC)||missing(nT)){
      stop("one or more arguments are not specified for the SMD measure")
    }
  }else if(measure=="RAW"){
    if (missing(mC)||missing(mT)||missing(nC)||missing(nT)||missing(sdC)||missing(sdT)){
      stop("one or more arguments are not specified for the RAW measure")
    }
    k<-length(mC)
    SMD<-escalc(measure="SMD",m1i=mC,sd1i=sdC,n1i=nC,m2i=mT,sd2i=sdT,n2i=nT)
    yi<-as.numeric(SMD$yi)
    vi<-as.numeric(SMD$vi)
  }
  if(missing(order)){
    order<-c(1:k)
  }
  yi<-yi[order]
  vi<-vi[order]
  nC<-nC[order]
  nT<-nT[order]
  for(i in 1:k){
    if(nC[i]<0 || nT[i]<0 || nC[i]+nT[i]<=0){
      stop("value of nC or nT have to be greater than 1 and nC+nT have to be greater than 8 to run the test")
    }
  }
  if (missing(level)) {
    stop("level is not specified")
  }
  if(k<10){
    stop("Number of studies is too small for stage 1 (less than 10); a fixed \\tau^2 and \\delta estimates cannot be obtain for stage 2 CMA and the tests for shift")
  }
  if (missing(delta)) {
    delta<-NA
  }
  tilde_n<-(nC*nT)/(nC+nT)
  
  ##################################################################################
  ##STAGE 1: Standard CMA (re-estimated \tau^2) using first 10 studies##
  if(is.na(delta)==TRUE){
    stg1delta<-0
  }else{
    stg1delta<-delta
  }
  ##IV-PM
  s_NA<-1
  standard<-rma(yi[1:10],vi[1:10],level=level*100,method="PM")
  cumul_stand<-cumul(standard)
  s_eiters<- cumul_stand$estimate
  s_eabove<- cumul_stand$ci.ub
  s_ebelow<- cumul_stand$ci.lb
  s_change<-((cumul_stand$ci.ub>=stg1delta)*(cumul_stand$ci.lb<=stg1delta))
  s_titers<- cumul_stand$tau2
  s_tau2CI<-c(NA)
  for (i in 2:10){
    a<-rma.uni(yi[1:i],vi[1:i],level=level*100)
    b<-confint.rma.uni(a,level=level*100)
    c<-b$random
    s_tau2CI<-rbind(s_tau2CI,c[1,])
  }
  s_tabove<- s_tau2CI[,3]
  s_tbelow<- s_tau2CI[,2]
  if(sum(is.na(cumul_stand$estimate)) == 0){
    s_NA <- s_NA-1
  }
  s_number<-1+match(TRUE, s_change[5:10] == 0)
  if(is.na(s_number)==TRUE){
    s_number<-10
    s_numberlist<- s_number
    s_delta_no <- cumul_stand$estimate[s_number]
    s_tau2_no <- cumul_stand$tau2[s_number]
  }else if (is.na(s_number)==FALSE){
    s_numberlist<- s_number
    s_delta_no <- cumul_stand$estimate[s_number]
    s_tau2_no <- cumul_stand$tau2[s_number-1]
  }
  #################################################################################
  ##STAGE 2: Fixed effect model (\tau^2 estimated from stage 1 added to variance)##
  ##IV-PM
  s_NA2<-1
  standardIV<-rma(yi,vi+s_tau2_no,method='FE',level=level*100)
  cumul_standIV<-cumul(standardIV)
  s_eiters2<- cumul_standIV$estimate
  s_eabove2<- cumul_standIV$ci.ub
  s_ebelow2<- cumul_standIV$ci.lb
  s_change2b<-((cumul_standIV$ci.ub>=s_delta_no)*(cumul_standIV$ci.lb<=s_delta_no))
  if(sum(is.na(cumul_standIV$estimate)) == 0){
    s_NA2 <- s_NA2-1
  }
  s_number2b<-1+match(TRUE, s_change2b[5:k] == 0)
  if(is.na(s_number2b)==FALSE){
    s_numberlist2b<-s_number2b
  }else if(is.na(s_number2b)==TRUE){
    s_numberlist2b<-NA
  }
  ##################################################################################
  ###**Produce Output**###
  ##Output of stage 1
  output<-cbind(s_eiters,s_eabove,s_ebelow,s_titers,s_tabove,s_tbelow)
  colnames(output)<-c('PM_cumulative_effect','PM_effect_ub','PM_effect_lb','PM_tau2','PM_tau2_ub','PM_tau2_lb')
  rownames(output)<-c()
  outputb<-cbind(s_delta_no, s_tau2_no, s_numberlist,s_NA)
  colnames(outputb)<-c('PM_delta_estimate', 'PM_tau2_estimate','PM_change_point','PM_NA_case')
  
  ##Output of stage 2 CMA
  output2<-cbind(s_eiters2,s_eabove2,s_ebelow2)
  colnames(output2)<-c('PM2_cumulative_effect','PM2_effect_ub','PM2_effect_lb')
  rownames(output2)<-c()
  output2b<-cbind(s_numberlist2b,s_NA2)
  colnames(output2b)<-c('PM2_change_point','PM2_NA_case')
  
  ##Input
  input<-cbind(delta,k,level)
  colnames(input)<-cbind('True \\delta','Number of studies','Confidence level')
  
  ##Final output
  finaloutput<-list('Input'=input,'Stg1output'=outputb,'Stg1output2'=output,'Stg2output'=output2b,'Stg2output2'=output2)
  return(finaloutput)
}
