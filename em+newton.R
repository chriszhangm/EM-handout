##################################################################################
#Newton-Rapson
firstderiv=function(theta){(125/4)/(0.5+0.25*theta)-38/(1-theta)+34/theta}
secondderiv = function(theta){-125/(2+theta)^2-38/(1-theta)^2-34/theta^2}

newtonway = function(){
	options(digits=16)
	iteration = numeric(10)
	pi0 = 0.5
	i = 1
	iteration[1] = 0.5
	while (i<=10){
		iteration[i+1] = pi0-firstderiv(pi0)/secondderiv(pi0)
		pi0 = iteration[i+1]
		i = i +1
	}
	return(iteration)	 
}
##################################################################################

#EM 

#y1=125, y2=18,y3=20,y4=34
emmulti = function(tol=1e-8){
	iteration1 = numeric(100)
	iteration10 = numeric(100)
	iteration100 = numeric(100)
	theta1 = 0.5
	theta10 = 0.5
	theta100 = 0.5
	iteration1[1] = theta1;iteration10[1] = theta10;iteration100[1] = theta100
	i=1
	while(i<=100){
		#x1 = 125*((theta0/4)/(0.5+theta0/4))  #e-step
		
		x1_1 = rbinom(1,125,(theta1/4)/(0.5+theta1/4))
		x1_10 = mean(rbinom(10,125,(theta10/4)/(0.5+theta10/4)))
		x1_100  =mean(rbinom(100,125,(theta100/4)/(0.5+theta100/4)))
		
		theta1 = (x1_1+34)/(x1_1+34+38) #m-step
		theta10 = (x1_10+34)/(x1_10+34+38)
		theta100 = (x1_100+34)/(x1_100+34+38)
		
		iteration1[i+1] = theta1
		iteration10[i+1] = theta10
		iteration100[i+1] = theta100

		i=i+1
	}	
	return(list('iteration1'=iteration1,'iteration10'=iteration10,'iteration100'=iteration100))
}
process1  = emmulti()
plot(1:100,process1$iteration1[1:100],type='b',pch=16,col='red',xlab='numbers of iterations',ylab='Estimated value',main='MCEM for Multinomial Case')
points(1:100,process1$iteration10[1:100],type='b',pch=3,col='blue')
points(1:100,process1$iteration100[1:100],type='b',pch=16,col='green')
legend('bottomright',legend=c('n=1','n=10','n=100'),col=c('red','blue','green'),pch=c(16,3,16))






#EM+GMM

### two component EM 
### pN(0,1)+(1-p)N(4,1)
EM_TwoMixtureNormal = function(p, mu1, mu2, sd1, sd2, X, maxiter=1000, tol=1e-10){
diff=1
iter=1
mu11 = numeric(100)
mu22 = numeric(100)
sigma11 = numeric(100)
sigma22 = numeric(100)
ppp = numeric(100)
mu11[1] = mu1
mu22[1] = mu2
sigma11[1] = sd1
sigma22[1] = sd2
ppp[1] = p
while (diff>tol & iter<maxiter) {
    ## E-step: compute omega:
    	d1=dnorm(X, mean=mu1, sd=sd1)
    	d2=dnorm(X, mean=mu2, sd=sd2)
    	omega=d1*p/(d1*p+d2*(1-p))  
    	#responsibility
    	# compute density in two groups
## M-step: update p, mu and sd
	p.new=mean(omega)
	mu1.new=sum(X*omega) / sum(omega)
	mu2.new=sum(X*(1-omega)) / sum(1-omega)
	resid1=X-mu1
	resid2=X-mu2
	sd1.new=sqrt(sum(resid1^2*omega) / sum(omega))
	sd2.new=sqrt(sum(resid2^2*(1-omega)) / sum(1-omega))
	## calculate diff to check convergence
	diff=sqrt(sum((mu1.new-mu1)^2+(mu2.new-mu2)^2+(sd1.new-sd1)^2+(sd2.new-sd2)^2))
	p=p.new;
	mu1=mu1.new;
	mu2=mu2.new;
	sd1=sd1.new;
	sd2=sd2.new;
	iter=iter+1
	mu11[iter] = mu1
	mu22[iter] = mu2
	sigma11[iter] = sd1
	sigma22[iter] = sd2
	ppp[iter] = p
	;
	cat("Iter", iter, ": mu1=", mu1.new, ", mu2=",mu2.new, ", sd1=", sd1.new,", sd2=",sd2.new, ", p=", p.new, ", diff=", diff, "\n")
	}
	return(list(mu11 = mu11,mu22=mu22,sigma11=sigma11,sigma22=sigma22,probility=ppp))
}

## simulation
p0=0.3;
n=5000;
X1=rnorm(n*p0);             # n*p0 indiviudals from N(0,1)
X2=rnorm(n*(1-p0), mean=4)  # n*(1-p0) individuals from N(4,1)
X=c(X1,X2)                  # observed data
hist(X, 50) 
## initial values for EM
p=0.5
mu1=1
mu2=5
sd1=sd2=sd(X)
c(p, mu1, mu2, sd1, sd2)
#0.5000000 -0.3903964  5.0651073  2.0738555  2.0738555
EM_TwoMixtureNormal(p, mu1, mu2, sd1, sd2, X)

process1 = EM_TwoMixtureNormal(p, mu1, mu2, sd1, sd2, X)
plot(1:60,process1$mu11[1:60],type='b',lty=3,pch=16,col='red',ylim=c(0,5),xlab='Numbers of Iterations',ylab='Estimated Value',main='EM for  N(4,1) and N(0,1) with p = 0.3')
points(1:60,process1$mu22[1:60],type='b',lty=3,pch=16,col='blue')
points(1:60,process1$sigma11[1:60],type='b',lty=1,pch=3,col='red')
points(1:60,process1$sigma22[1:60],type='b',lty=1,pch=3,col='blue')
points(1:60,process1$probility[1:60],type='b',lty=5,pch=16,col='green')

legend('topright', legend=c('mu1','mu2','sigma1','sigma2','p'),col=c('red','blue','red','blue','green'),lty=c(3,3,1,1,5),pch=c(16,16,3,3,16),cex=0.5)







