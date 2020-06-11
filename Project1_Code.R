#데이터 불러오기
x1 = read.csv('판매이력.csv')
broken = read.csv('고장이력.csv')

#library
library(dplyr)

#데이터 중복 여부 확인
length(unique(x1$ID))
length(unique(broken$ID))
#데이터 row수와 unique ID(PK) length가 같으므로 데이터내의 중복은 없다.


#판매이력에서 고장난거들 빼서 안고장난 데이터(unbroken) 추출
unbroken = anti_join(x1,broken,by='ID')
length(unique(unbroken[,1]))==length(unique(x1[,1]))-length(unique(broken[,1]))
#판매이력중 ID가 겹치는 고장이력 부분을 제외한 데이터의 수  = 판매이력의 데이터수 - 고장이력 데이터수
#따라서 판매이력이 누락된 고장데이터는 없음


#(데이터수집날짜 - 판매날짜) 이용해 고장나지 않은 데이터의 사용일수 추출

unbroken$판매일자 = as.Date(unbroken$판매일자)

date_get = as.Date('2017-3-30')

unbroken$사용일수 = date_get - unbroken$판매일자


#우측관측중단 데이터 생성 -> 사일자가 1095이상인 데이터 = 우측관측중단
unbroken1 = unbroken %>%
  filter(사용일수>=1095)

unbroken2 = unbroken %>%
  filter(사용일수<1095)

#3년 이상의 우측중단데이터는 사용일수를 1095일로 변환
unbroken1$사용일수 = 1095

unbroken = rbind(unbroken1,unbroken2)

#broken unbroken 의 우측관측중단 column - cens(0,1)추가후 병합
unbroken = unbroken%>%
  mutate(cens = c(rep(0,length(unbroken[,1]))))

unbroken = unbroken[,-3] #column level 맞춰주기 위해 연산에 필요없는 '판매일자'는 삭제

broken = broken %>%
  mutate(cens=c(rep(1,length(broken[,1]))))

data = rbind(broken,unbroken)

head(data,15)
######################################################################
#####경험적 분포######
#median rank의 경험적분포 사용
x = sort(broken$사용일수)
length(x)
cdf = ((1:length(x)-0.3))/(length(data[,1])+0.4)  #전체데이터중에서 완전고장데이터(약5%)만의 경험적 누적분포함수
plot(x,cdf,main='경험적 누적분포함수',xlab='time')




#####수명분포 추정1. weibull##########################################
##우도함수 정의(weibull)
Lik = function(x,param){
  times = x$사용일수
  cens = x$cens
  
  value = -sum(cens*log(dweibull(times,param[1],param[2]))+
                 (1-cens)*log(1-pweibull(times,param[1],param[2])))
  return(value)
}

#모수추정
fit =optim(c(1,quantile(data[,3],0.632)),Lik,x=data,hessian=TRUE)

fit

param = fit$par
param

##모수의 신뢰구간 구하기
mat = fit$hessian
std_err = sqrt(diag(solve(mat)))
std_err

conf_interval = c(param-qnorm(0.975)*std_err,param+qnorm(0.975)*std_err)
conf_interval #정규분포가정의 근사적 95%신뢰구간

conf_interval_lnorm = c(param*exp(-qnorm(0.975)*std_err/param),param*exp(qnorm(0.975)*std_err/param))
conf_interval_lnorm #대수정규분포가정 근사적 95% 신뢰구간




#Bp 수명, 3년시점의 신뢰도, MTTF, 3년시점의 누적고장(%)
weibull_Bp = qweibull(c(0.01,0.03,0.05,0.1),param[1],param[2])
weibull_Bp
weibull_R = 1 - pweibull(1095,param[1],param[2])
weibull_R
weibull_MTTF = param[2]*gamma(1+1/param[1])
weibull_MTTF
pweibull(1095,param[1],param[2])

#####수명분포 추정2. 대수정규##########################################

##우도함수 정의(lnorm)
Lik2 = function(x,param){
  times = x$사용일수
  cens = x$cens
  
  value = -sum(cens*log(dlnorm(times,param[1],param[2]))+
                 (1-cens)*log(1-plnorm(times,param[1],param[2])))
  return(value)
}

##모수추정
fit2 = optim(c(1,1),Lik2,x=data,hessian=TRUE)
fit2
param2=fit2$par
param2
##모수의 신뢰구간 구하기
mat2 = fit2$hessian
std_err2 = sqrt(diag(solve(mat2)))
std_err2

conf_interval2 = c(param2-qnorm(0.975)*std_err2,param2+qnorm(0.975)*std_err2)
conf_interval2 #정규분포가정의 근사적 95%신뢰구간

conf_interval_lnorm2 = c(param2*exp(-qnorm(0.975)*std_err2/param2),param2*exp(qnorm(0.975)*std_err2/param2))
conf_interval_lnorm2 #대수정규분포가정 근사적 95% 신뢰구간


##Bp수명, 3년시점 신뢰도,MTTF,3년시점 누적고
lnorm_Bp = qlnorm(c(0.01,0.03,0.05,0.1),param2[1],param2[2])
lnorm_Bp
R_lnorm =1- plnorm(1095,param2[1],param2[2])
R_lnorm
lnorm_MTTF = exp(param2[1]+(param2[2])^2/2)
lnorm_MTTF  
plnorm(1095,param2[1],param2[2])
#####수명분포 추정3. 지수분포##########################################
Lik3 = function(x){
  times = x$사용일수
  cens = x$cens
  
  lambda = length(broken[,1])/sum(times)
  
  value = sum(cens*log(dexp(times,lambda))+
                 (1-cens)*log(1-pexp(times,lambda)))
  return(c(lambda,value))
}

fit3 = Lik3(data)
fit3


#Bp수명, 3년시점 신뢰도, MTTF, 3년시점의 누적고장(%)
Bp_exp = qexp(c(0.01,0.03,0.05,0.1),fit3[1])
Bp_exp
R_exp = 1 - pexp(1095,fit3[1])
R_exp
MTTF_exp = 1/fit3[1]
MTTF_exp
pexp(1095,fit3[1])

##########################적합도검정#######################

##도시적 검정
##최우추정을 통해 얻은 이론적 분포(weibull)의 누적분포함수 경험적 누적분포함수 확률지 비교
par(mfrow=c(3,1))

#와이블분포
plot(log(x),log(-log(1-cdf)),type='p',main = '와이블분포',xlab='log(t)')
lines(log(x),log(-log(1-pweibull(x,param[1],param[2]))),col=2,lwd=2)

#대수정규분포
plot(log(x),qnorm(cdf),main = '대수정규분포',xlab='log(t)')
lines(log(x),qnorm(plnorm(x,param2[1],param2[2])),col=2,lwd=2)

#지수분포
plot(fit3[1]*x,-log(1-cdf),main = '지수분포',xlab='lambda*t')
lines(fit3[1]*x,-log(1-pexp(x,fit3[1])),col=2,lwd=2)


#Information Criteria를 이용한 적합도검정
AIC_Weibull = 2*fit$value + 2*length(param)
AIC_lnorm = 2*fit2$value + 2*length(param2)
AIC_exp = -2*fit3[2]+2*1

AIC_Weibull
AIC_lnorm
AIC_exp


BIC_weibull = 2*fit$value + length(param)*log(length(data[,1]))
BIC_lnorm = 2*fit2$value + length(param2)*log(length(data[,1]))
BIC_exp = -2*fit3[2]+ 1*log(length(data[,1]))

BIC_weibull
BIC_lnorm
BIC_exp

#Weibull이 가장 좋은 적합도
