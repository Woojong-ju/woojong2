#������ �ҷ�����
x1 = read.csv('�Ǹ��̷�.csv')
broken = read.csv('�����̷�.csv')

#library
library(dplyr)

#������ �ߺ� ���� Ȯ��
length(unique(x1$ID))
length(unique(broken$ID))
#������ row���� unique ID(PK) length�� �����Ƿ� �����ͳ��� �ߺ��� ����.


#�Ǹ��̷¿��� ���峭�ŵ� ���� �Ȱ��峭 ������(unbroken) ����
unbroken = anti_join(x1,broken,by='ID')
length(unique(unbroken[,1]))==length(unique(x1[,1]))-length(unique(broken[,1]))
#�Ǹ��̷��� ID�� ��ġ�� �����̷� �κ��� ������ �������� ��  = �Ǹ��̷��� �����ͼ� - �����̷� �����ͼ�
#���� �Ǹ��̷��� ������ ���嵥���ʹ� ����


#(�����ͼ�����¥ - �Ǹų�¥) �̿��� ���峪�� ���� �������� ����ϼ� ����

unbroken$�Ǹ����� = as.Date(unbroken$�Ǹ�����)

date_get = as.Date('2017-3-30')

unbroken$����ϼ� = date_get - unbroken$�Ǹ�����


#���������ߴ� ������ ���� -> �����ڰ� 1095�̻��� ������ = ���������ߴ�
unbroken1 = unbroken %>%
  filter(����ϼ�>=1095)

unbroken2 = unbroken %>%
  filter(����ϼ�<1095)

#3�� �̻��� �����ߴܵ����ʹ� ����ϼ��� 1095�Ϸ� ��ȯ
unbroken1$����ϼ� = 1095

unbroken = rbind(unbroken1,unbroken2)

#broken unbroken �� ���������ߴ� column - cens(0,1)�߰��� ����
unbroken = unbroken%>%
  mutate(cens = c(rep(0,length(unbroken[,1]))))

unbroken = unbroken[,-3] #column level �����ֱ� ���� ���꿡 �ʿ���� '�Ǹ�����'�� ����

broken = broken %>%
  mutate(cens=c(rep(1,length(broken[,1]))))

data = rbind(broken,unbroken)

head(data,15)
######################################################################
#####������ ����######
#median rank�� ���������� ���
x = sort(broken$����ϼ�)
length(x)
cdf = ((1:length(x)-0.3))/(length(data[,1])+0.4)  #��ü�������߿��� �������嵥����(��5%)���� ������ ���������Լ�
plot(x,cdf,main='������ ���������Լ�',xlab='time')




#####�������� ����1. weibull##########################################
##�쵵�Լ� ����(weibull)
Lik = function(x,param){
  times = x$����ϼ�
  cens = x$cens
  
  value = -sum(cens*log(dweibull(times,param[1],param[2]))+
                 (1-cens)*log(1-pweibull(times,param[1],param[2])))
  return(value)
}

#�������
fit =optim(c(1,quantile(data[,3],0.632)),Lik,x=data,hessian=TRUE)

fit

param = fit$par
param

##����� �ŷڱ��� ���ϱ�
mat = fit$hessian
std_err = sqrt(diag(solve(mat)))
std_err

conf_interval = c(param-qnorm(0.975)*std_err,param+qnorm(0.975)*std_err)
conf_interval #���Ժ��������� �ٻ��� 95%�ŷڱ���

conf_interval_lnorm = c(param*exp(-qnorm(0.975)*std_err/param),param*exp(qnorm(0.975)*std_err/param))
conf_interval_lnorm #������Ժ������� �ٻ��� 95% �ŷڱ���




#Bp ����, 3������� �ŷڵ�, MTTF, 3������� ��������(%)
weibull_Bp = qweibull(c(0.01,0.03,0.05,0.1),param[1],param[2])
weibull_Bp
weibull_R = 1 - pweibull(1095,param[1],param[2])
weibull_R
weibull_MTTF = param[2]*gamma(1+1/param[1])
weibull_MTTF
pweibull(1095,param[1],param[2])

#####�������� ����2. �������##########################################

##�쵵�Լ� ����(lnorm)
Lik2 = function(x,param){
  times = x$����ϼ�
  cens = x$cens
  
  value = -sum(cens*log(dlnorm(times,param[1],param[2]))+
                 (1-cens)*log(1-plnorm(times,param[1],param[2])))
  return(value)
}

##�������
fit2 = optim(c(1,1),Lik2,x=data,hessian=TRUE)
fit2
param2=fit2$par
param2
##����� �ŷڱ��� ���ϱ�
mat2 = fit2$hessian
std_err2 = sqrt(diag(solve(mat2)))
std_err2

conf_interval2 = c(param2-qnorm(0.975)*std_err2,param2+qnorm(0.975)*std_err2)
conf_interval2 #���Ժ��������� �ٻ��� 95%�ŷڱ���

conf_interval_lnorm2 = c(param2*exp(-qnorm(0.975)*std_err2/param2),param2*exp(qnorm(0.975)*std_err2/param2))
conf_interval_lnorm2 #������Ժ������� �ٻ��� 95% �ŷڱ���


##Bp����, 3����� �ŷڵ�,MTTF,3����� ������
lnorm_Bp = qlnorm(c(0.01,0.03,0.05,0.1),param2[1],param2[2])
lnorm_Bp
R_lnorm =1- plnorm(1095,param2[1],param2[2])
R_lnorm
lnorm_MTTF = exp(param2[1]+(param2[2])^2/2)
lnorm_MTTF  
plnorm(1095,param2[1],param2[2])
#####�������� ����3. ��������##########################################
Lik3 = function(x){
  times = x$����ϼ�
  cens = x$cens
  
  lambda = length(broken[,1])/sum(times)
  
  value = sum(cens*log(dexp(times,lambda))+
                 (1-cens)*log(1-pexp(times,lambda)))
  return(c(lambda,value))
}

fit3 = Lik3(data)
fit3


#Bp����, 3����� �ŷڵ�, MTTF, 3������� ��������(%)
Bp_exp = qexp(c(0.01,0.03,0.05,0.1),fit3[1])
Bp_exp
R_exp = 1 - pexp(1095,fit3[1])
R_exp
MTTF_exp = 1/fit3[1]
MTTF_exp
pexp(1095,fit3[1])

##########################���յ�����#######################

##������ ����
##�ֿ������� ���� ���� �̷��� ����(weibull)�� ���������Լ� ������ ���������Լ� Ȯ���� ��
par(mfrow=c(3,1))

#���̺�����
plot(log(x),log(-log(1-cdf)),type='p',main = '���̺�����',xlab='log(t)')
lines(log(x),log(-log(1-pweibull(x,param[1],param[2]))),col=2,lwd=2)

#������Ժ���
plot(log(x),qnorm(cdf),main = '������Ժ���',xlab='log(t)')
lines(log(x),qnorm(plnorm(x,param2[1],param2[2])),col=2,lwd=2)

#��������
plot(fit3[1]*x,-log(1-cdf),main = '��������',xlab='lambda*t')
lines(fit3[1]*x,-log(1-pexp(x,fit3[1])),col=2,lwd=2)


#Information Criteria�� �̿��� ���յ�����
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

#Weibull�� ���� ���� ���յ�