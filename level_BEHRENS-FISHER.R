rm(list=ls())
library(asht)

simulate = function(R,n1,n2,var1, var2)
{
  
  sigma1 = sqrt(var1); sigma2 = sqrt(var2)
  
  # Sample 1
  s1 = matrix(rnorm(R*n1, 20, sigma1), nrow = n1, ncol = R)
  sm_1 = apply(s1, 2, mean)
  svar_1 = apply(s1, 2, var)
  
  # Sample 2
  s2 = matrix(rnorm(R*n2, 20, sigma2), nrow = n2, ncol = R)
  sm_2 = apply(s2, 2, mean)
  svar_2 = apply(s2, 2, var)
  
  # Pooled variance
  var_pooled=(((n1-1)*svar_1) + ((n2-1)*svar_2))/(n1+n2-2)
  
  # Test Statistic
  Z=(sm_1-sm_2)/(sqrt(var_pooled)*sqrt((1/n1)+(1/n2)))
  
  #Fiducial Test Statistic
  Z1=(sm_1-sm_2)/sqrt((svar_1/n1)+(svar_2/n2))
  z2 = (sm_1-sm_2)
  
  crit = qt(0.025,n1+n2-2,lower.tail=F)
  t1 = qt(0.025,n1-1,lower.tail = F)
  t2 = qt(0.025,n2-1,lower.tail = F)
  fiducial_level = 0; welch_level = 0;banerjee_level = 0
  count = 0
  count1 = 0
  count2= 0
  
  t1 = qt(0.025,n1-1,lower.tail = F)
  t2 = qt(0.025,n2-1,lower.tail = F)
  
  
  for(i in 1:R){
    crit1=qbf(0.025,n1,n2,s1=sqrt(svar_1[i]),s2=sqrt(svar_2[i]))
    crit2=qbf(0.975,n1,n2,s1=sqrt(svar_1[i]),s2=sqrt(svar_2[i]))
    welch.df = (((svar_1[i]/n1)+(svar_2[i]/n2))^2)/((((svar_1[i])^2)/((n1^2)*(n1-1)))+(((svar_2[i])^2)/((n2^2)*(n2-1))))
    welch.crit = qt(0.025,welch.df,lower.tail = F)
    w1 = svar_1[i]/n1
    w2 = svar_2[i]/n2
    c1 = w1*(t1^2)
    c2 = w2*(t2^2)
    crit3 = sqrt(c1+c2)
    
    if((Z1[i]<crit1) | (Z1[i] > crit2)){
      count = count+1
    }
    if(abs(Z1[i])>welch.crit){
      count1 = count1+1
    }
    if(abs(z2[i]) > crit3){
      count2 = count2+1
    }
  }
  fiducial_level = count/R
  welch_level = count1/R
  banerjee_level = count2/R
  
  # Proportion
  fisher_t_level = mean(abs(Z) > crit)
  
  
  return(c(fisher_t_level,fiducial_level,welch_level,banerjee_level))
}
#Reporting the levels 
var1=rep(5,times=9);var1
var2=c(5,10,25,40,60,90,120,150,200)
n1 =  rep(100,times=9)
n2 = rep(100,times = 9)

fisher.t=array(0)
fidu=array(0)
welch=array(0)
banerjee = array(0)
for(i in 1:length(var2)){
  set.seed(seed=987654321)
  level.func = simulate(5000,100,100,5,var2[i])
  fisher.t[i]=level.func[1]
  fidu[i]=level.func[2]
  welch[i]=level.func[3]
  banerjee[i]=level.func[4]
}
data=data.frame(n1,n2,var1,var2,fisher.t,fidu,welch,banerjee)
data