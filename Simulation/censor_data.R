
### Function generates data according to a weibull distrubution or uniform
### That can be used to compare and censor
censoring_data = function(censoring_p, type,
                          betaT, Treatment, LP, eta, lambda){
  n = length(LP)
  alpha.t = eta
  theta1 <- sqrt(1/lambda)
  lambda_i <- theta1*exp(-(betaT * Treatment + LP)/alpha.t)
  lambda_i_final = lambda_i
  max.lambda_i<-max(lambda_i)
  
  density.fun.lambda<-function(x){
    pred.y <- predict(y.loess, newdata=x)
    return(pred.y)
  }
  ### compute the density function ###
  h = 1.06*min(sd(lambda_i), IQR(lambda_i)/1.34)*(n^(-0.2))
  dens<-density(lambda_i,n=n,bw=h,from=0,to=max.lambda_i,na.rm=TRUE)
  x<-dens$x
  y<-dens$y
  y.loess <-loess(y~x,span=0.1, control=loess.control(surface="direct"))
  #plot(dens,lty=1,lwd=2)
  #lines(y.loess,col=2)
  #integrate(function(u){predict(y.loess, newdata=u)}, 0, max.lambda_i)
  
  censor.prop<-function(theta,args){
    censoring_p<-as.numeric(args[1])
    cen.P<-integrate(function(u){
      alpha.t = as.numeric(args[2])
      type = args[3]
      lambda.i<-u
      part1<-density.fun.lambda(lambda.i)
      if (type == 'weibull'){part2 = 1/(1 + ((theta/lambda.i)^alpha.t))}
      if (type == 'uniform'){
        part2 = (lambda.i/(alpha.t*theta))*pgamma((theta/lambda.i)^alpha.t,1/alpha.t)*gamma(1/alpha.t)}
      return(part1*part2)
    },0,max.lambda_i, subdivisions=1000, rel.tol = 2e-4)$value
    return(cen.P-censoring_p)
  }
  args<-c(censoring_p,alpha.t, type)
  theta<-uniroot(censor.prop,args=args,c(0.0001,10000),tol=0.00000001)$root
  if (type == 'weibull'){censor_values = rweibull(n = n, shape = eta, scale = theta)}
  if (type == 'uniform'){censor_values = runif(n = n, min = 0, max = theta)}
  return(list(censor_values = censor_values, theta = theta, lambda_i = lambda_i))
}
