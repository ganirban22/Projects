rm(list = ls())
library(asht)
f1 = function(n1,n2)
{
 power_func = function(R,n1,n2,mu1,mu2,var1, var2)
{
  sigma1 = sqrt(var1); sigma2 = sqrt(var2)
  
  # Sample 1
  s1 = matrix(rnorm(R*n1, mu1, sigma1), nrow = n1, ncol = R)
  sm_1 = apply(s1, 2, mean)
  svar_1 = apply(s1, 2, var)
  
  # Sample 2
  s2 = matrix(rnorm(R*n2, mu2, sigma2), nrow = n2, ncol = R)
  sm_2 = apply(s2, 2, mean)
  svar_2 = apply(s2, 2, var)
  
  # Pooled variance
  var_pooled=(((n1-1)*svar_1) + ((n2-1)*svar_2))/(n1+n2-2)
  
  #test.stat.fisher=test.stat.fidu=array(dim=1)
  #power
  c = qt(0.975,n1+n2-2,lower.tail=T)
  test.stat.fisher=(sm_1-sm_2)/(sqrt(var_pooled) * sqrt((1/n1)+(1/n2)))
  test.stat.fidu = test.stat.welch = (sm_1-sm_2)/(sqrt((svar_1/n1)+(svar_2/n2)))
  test.stat.banerjee = sm_1-sm_2
  #k1 = (c - zi_1)/(sqrt(var_pooled) * sqrt((1/n1)+(1/n2)))
  #k2 = (-c - zi_1)/(sqrt(var_pooled) * sqrt((1/n1)+(1/n2)))
  #p1 = pt(k1 , df = n1+n2-2 , lower.tail = F)
  #p2 = pt(k2 , df = n1+n2-2 , lower.tail = T)
  
  #power_fisher =mean( p1 + p2 )
  
  
  power.fidu = power.welch = 0
  count1 = count3 = count2 = 0
  power.banerjee = 0
  
  t1 = qt(0.025,n1-1,lower.tail = F)
  t2 = qt(0.025,n2-1,lower.tail = F)
  for(i in 1:R){
    d1=qbf(0.025,n1,n2,s1=sqrt(svar_1[i]),s2=sqrt(svar_2[i]))
    d2=qbf(0.975,n1,n2,s1=sqrt(svar_1[i]),s2=sqrt(svar_2[i]))
    welch.df = ((svar_1[i]/n1+svar_2[i]/n2)^2)/((((svar_1[i])^2)/((n1^2)*(n1-1)))+(((svar_2[i])^2)/((n2^2)*(n2-1))))
    welch.crit = qt(0.025,welch.df,lower.tail = F)
    if(test.stat.fidu[i] > d2 | test.stat.fidu[i] < d1 ){
      count1 = count1+1
    }
    if(abs(test.stat.welch[i]) > welch.crit){
      count2 = count2 +1
    }
    w1 = svar_1[i]/n1
    w2 = svar_2[i]/n2
    c1 = w1*(t1^2)
    c2 = w2*(t2^2)
    crit = sqrt(c1+c2)
    if(abs(test.stat.banerjee[i]) > crit){
      count3=count3+1
    }
  }
  power.fidu = count1/R
  power.welch = count2/R
  power.banerjee = count3/R
  
  return(c(power.fidu,power.welch,power.banerjee))
}
var1= rep(1,times=6)
var2= c(1,4,7,10,15,20)

zi_1= seq(0,5,length.out=10)
k = length(zi_1)
emp.power.fiducial = matrix(0,nrow = length(var2),ncol = length(zi_1))
emp.power.welch = matrix(0,nrow = length(var2),ncol = length(zi_1))
emp.power.banerjee = matrix(0,nrow = length(var2),ncol = length(zi_1))
set.seed(seed = 987654321)
for(i in 1:length(var2))
{
  set.seed(seed = 987654321) 
  for(j in 1:k)
  {
    pfunc = power_func(1000,n1,n2,0,0+zi_1[j],1,var2[i])
    emp.power.fiducial[i,j] = pfunc[1]
    emp.power.welch[i,j] = pfunc[2]
    emp.power.banerjee[i,j] = pfunc[3]
  }
}
dt = data.frame(id = gl(6, 10), 
                 mu1 = rep(0,times = 60), 
                 mu2 = rep(0+zi_1,times = 6), 
                 sigma1.sq = rep(1, times = 60), 
                 sigma2.sq = rep(c(1,4,7,10,15,20),each = 10), 
                 Fiducial  = c(t(emp.power.fiducial)), 
                 Welch  = c(t(emp.power.welch)),
                 Banerjee = c(t(emp.power.banerjee))); dt

theme_new=theme(plot.title= element_text(size=16, 
                                         hjust=0.5,
                                         face= "bold"),
                plot.subtitle= element_text(size=13, 
                                            hjust=0.5,
                                            face= "bold"),
                legend.title = element_text(size=14,
                                            face="bold"),
                axis.title = element_text(face="bold"),
                axis.text = element_text(face= "bold"))

dt %>% tidyr::pivot_longer(-(1:5), "Index", "Power") %>% 
  ggplot(aes(mu2-mu1, value, col = Index)) + 
  geom_line(lwd = 1) + 
  facet_wrap(sigma1.sq ~ sigma2.sq,labeller = label_bquote(sigma[1]^2==.(sigma1.sq)~","~sigma[2]^2==.(sigma2.sq)))+
  labs(title="Power Curves",
       subtitle = bquote(n[1]==.(n1)~","~n[2]==.(n2)),
       x=bquote(mu[2]-mu[1]),
       y="Power",
       col="Index")+
  theme_new
}
f1(5,5)
f1(5,10)
f1(10,5)
f1(10,10)
f1(10,25)
f1(25,10)
f1(25,25)
f1(25,50)
f1(50,25)
f1(50,50)
f1(50,100)
f1(100,50)
f1(100,100)
