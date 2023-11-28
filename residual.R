####  RESIDUALS FROM META REGRESSION AND THEIR CUMULATIVE META_ANALYSIS 
####  function requires a model matrix m from meta-regression
####  and columns with corresponding year, publication, study variables
####  output is a list of four data frames
####  mac[[1]] data from cumulative meta-analysis of residuals at the study level
####  mac[[2]] data from cumulative meta-analysis of residuals at the publication level
####  mac[[3]] data from cumulative meta-analysis of residuals at the year level
####  mac[[4]] residuals from  meta-regression m at the year level
####  output is sorted by year and publication (in alphabetical order) within year
#### 
####  Version 10/06/2020, Elena Kulinskaya

cumul_ma2<-function( year, publication, study, m) {
yy<-resid(m)[complete.cases(resid(m))] ### observed values
res<-resid(m)[complete.cases(resid(m))] ### residuals
year<-year[complete.cases(resid(m))]
publication<-publication[complete.cases(resid(m))]
study<-study[complete.cases(resid(m))]

npubl<-as.integer(publication)
nstudy<-as.integer(study)
nyear<-as.integer(as.factor(year))

 M<-m$M ### variances 
 X<-m$X
 
resdat<-data.frame(yy,res,year, publication, study,nyear,npubl,nstudy)

resdat <-resdat[order(nyear,npubl, nstudy),]

M<-M[order(nyear,npubl, nstudy),order(nyear,npubl, nstudy)]
X<-X[order(nyear,npubl, nstudy),]
H<-X%*%solve(t(X)%*%solve(M)%*%X)%*%t(X)%*%solve(M) ### hat matrix


K<-length(yy)  ###studies
Yrs<-length(unique(year))   ### Number of years
  
Pbs<-length(unique(npubl))  ### Number ofpublications

 theta_cc<-numeric(K)
 vi_cc<-numeric(K)

 for(c in 1:K){
  x_c<-rep(1,c)
H_c<-H[1:c,1:c]
 M_c<-M[1:c,1:c]
 yy_c<-resdat$res[1:c]
W_c<-(t(x_c)%*%solve(M_c)%*%x_c) ##### sum of weights
theta_c<-(1/W_c)*(t(x_c)%*%solve(M_c)%*%yy_c)    ### corrected
 var_c<- W_c^(-2)*(				 ### corrected
t(x_c)%*%solve(M_c)%*%x_c-t(x_c)%*%solve(M_c)%*%H_c%*%x_c)
theta_cc[c]<-theta_c
 vi_cc[c]<-var_c }

ci_l<-theta_cc-qnorm(.975)*sqrt(vi_cc)
ci_u<-theta_cc+qnorm(.975)*sqrt(vi_cc)
 p<-2*pnorm(-abs(theta_cc/sqrt(vi_cc)))

study<-resdat$study
publication<-resdat$publication
year<-resdat$year
mac<-data.frame(theta_cc,vi_cc,ci_l,ci_u,p,study,publication,year)
mac_publ<-mac[!duplicated(resdat$npubl, fromLast = TRUE),][-6]
mac_year<-mac[!duplicated(resdat$nyear, fromLast = TRUE),][-c(6,7)]

year<-year[!duplicated(resdat$nyear,fromLast = TRUE)]

theta_y<-numeric(Yrs)
 vi_y<-numeric(Yrs)
jjU<-jjL<-0
for(j in 1:Yrs){
yy<-resdat$res[nyear==j]
 jjL<-jjU
 jjU<-jjL+length(yy)
  x_j<-rep(1,length(yy))
H_c<-H[(jjL+1):jjU,(jjL+1):jjU]
 M_c<-M[(jjL+1):jjU,(jjL+1):jjU]
 
W_c<-(t(x_j)%*%solve(M_c)%*%x_j) ##### sum of weights
theta_c<-(1/W_c)*(t(x_j)%*%solve(M_c)%*%yy)    ### corrected
 var_c<- W_c^(-2)*(				 ### corrected
t(x_j)%*%solve(M_c)%*%x_j-t(x_j)%*%solve(M_c)%*%H_c%*%x_j)
theta_y[j]<-theta_c
 vi_y[j]<-var_c 
res_y<-data.frame(theta_y,vi_y,year)}

 return(list(mac,mac_publ,mac_year,res_y))
 }



