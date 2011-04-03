Metro_Hastings <-
function(li_func,pars,prop_sigma=diag(1,length(pars)),par_names=NULL,iterations=50000,burn_in=1000,adapt_par=c(100,20,0.5,0.75),quiet=FALSE,...)
{
    if (!is.finite(li_func(pars, ...))) 
        stop("Seed parameter values <pars> are not in the defined parameter space.  Try new starting values for <pars>.")

    if (is.null(par_names))
        par_names<-letters[1:length(pars)]

    prop_sigma<-makePositiveDefinite(prop_sigma)
	mu<-pars
	pi_X<-li_func(pars,...)			        # initial likelihood evaluation
	k_X<-pars
	trace<-array(dim=c(iterations,length(pars)))
    deviance<-array(dim=iterations) 
    announce<-floor(seq(iterations/10,iterations,length.out=10))
	for(i in 1:iterations)
	{
	        k_Y<-rmnorm(1,Sigma=prop_sigma,mu=k_X)	# Draw proposal point
        				
        	pi_Y<-li_func(k_Y,...)		        # evaluate likelihood at proposal point
        		
        	a_X_Y = (pi_Y)-(pi_X)		        # Compare relative likelihoods
        	if(is.nan(a_X_Y))
                    a_X_Y<--Inf                          # never jump outside the range of the parameter space
        	if( log(runif(1,0,1)) <= a_X_Y)	        # Make the jump according to MH probability
              	{
        		k_X = k_Y
        		pi_X = pi_Y
        	}	
        
        	trace[i,]<-k_X	
                deviance[i]<-(-2*pi_X)                      # Store the deviance for calculating DIC
		if(i > adapt_par[1] && i %% adapt_par[2] == 0 && i > (adapt_par[4]*iterations) )	        # adapt the proposal covariance structure
		{   
                    len<-floor(i*adapt_par[3]):i
                    x<-trace[len,]
                    N<-length(len)
                    p_sigma <- (N-1) * var(x)/N
                    p_sigma <-makePositiveDefinite(p_sigma)   # To deal with rounding problems that can de-symmetrize
                    if(!(0 %in% p_sigma) ) 
                        prop_sigma<-p_sigma
         }
         if(!quiet && i %in% announce)
            print(paste('updating: ',i/iterations*100,'%',sep=''))
	}
    trace<-trace[burn_in:iterations,]
    DIC<-NULL
    ## Calculate the DIC
    if(length(pars)>1)
    {
        D_bar<-mean(deviance[burn_in:iterations])
        theta_bar<-sapply(1:length(pars),function(x){mean( trace[,x] )})
        D_hat<-li_func(theta_bar,...)
        pD<-D_bar-D_hat
        DIC<-D_hat + 2*pD
    }
    ##
    if(length(pars)>1)
        accept_rate<-length(unique(trace[,1]))/(iterations-burn_in) else
    accept_rate<-length(unique(trace))/(iterations-burn_in)

	return(list("trace"=trace,"prop_sigma"=prop_sigma,"par_names"=par_names,"DIC"=DIC,'acceptance_rate'=accept_rate))
}

