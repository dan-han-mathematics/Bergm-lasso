  #' Bayesian Adaptive Lasso Method on ERGM models
  #'
  #' Function to fit Bayesian exponential random graphs models
  #' using the approximate exchange algorithm.
  #'
  #' @param formula formula; 
  #' an \code{\link[ergm]{ergm}} formula object,
  #' of the form  <network> ~ <model terms>
  #' where <network> is a \code{\link[network]{network}} object
  #' and <model terms> are \code{ergm-terms}.
  #'
  #' @param mean.prior vector; 
  #' mean vector of the multivariate Normal prior.
  #' By default set to a vector of 0's.
  #'
  #' @param sigma.prior square matrix; 
  #' variance/covariance matrix for the multivariate Normal prior.
  #' By default set to a diagonal matrix with every diagonal entry equal to 100.
  #' 
  #' @param burn.in count; 
  #' number of burn-in iterations for every chain of the population.
  #'
  #' @param main.iters count; 
  #' number of iterations for every chain of the population.
  #'
  #' @param aux.iters count; 
  #' number of auxiliary iterations used for network simulation.
  #'
  #' @param nchains count; 
  #' number of chains of the population MCMC.
  #' By default set to twice the model dimension (number of model terms).
  #'
  #' @param gamma scalar; 
  #' parallel adaptive direction sampling move factor.
  #' 
  #' @param V.proposal count; 
  #' diagonal entry for the multivariate Normal proposal.
  #' By default set to 0.0025.
  #'
  #' @param startVals vector;
  #' optional starting values for the parameter estimation. 
  #'
  #' @param offset.coef vector;
  #' A vector of coefficients for the offset terms.
  #'
  #' @param ... additional arguments, to be passed to lower-level functions.
  #'
  #' @references
  #' Caimo, A. and Friel, N. (2011), "Bayesian Inference for Exponential Random Graph Models,"
  #' Social Networks, 33(1), 41-55. \url{https://arxiv.org/abs/1007.5192}
  #'
  #' Caimo, A. and Friel, N. (2014), "Bergm: Bayesian Exponential Random Graphs in R,"
  #' Journal of Statistical Software, 61(2), 1-25. \url{https://www.jstatsoft.org/article/view/v061i02}
  
  
  
  
  
  library(statnet)
  library(ergm)
  library(Bergm)
  library(mvtnorm)
  # library(VGAM) # rinv.gaussian
  library(matrixStats)
  library(mcmc)
  library(coda)
  library(statmod)
  library(MCMCpack)
  
  BergmLasso_full_bayes <- function(formula, 
                         data = NULL,
                         sigma_prior_diag_value=NULL, 
                         lambda=5,
                         mean_prior=NULL , 
                         invTau2=1/10,
                         aux.iters=100,
                         main.iters=1000,
                         burn.in =200, 
                         V.proposal=0.0025,
                         gamma=0.5, 
                         nchains=NULL,
                         startVals=NULL,
                         offset.coef=NULL,
                         thin = 1,
                         r0=1,
                         delta0=0.5)
  {
  
    net <-data
    
    y <- ergm.getnetwork(formula)
    model <- ergm_model(formula, y)
    sy <- summary(formula)
    dim <- length(sy)
    if (is.null( sigma_prior_diag_value)){
      sigma.prior<-100*diag(1/invTau2,dim)
    } else {
      sigma.prior <- sigma_prior_diag_value*diag(1/invTau2,dim)
      }
    
   #regression results by bayesian addpative lasso
    lambda2<-rep(lambda^(2), dim)
  
    tot.iters <- burn.in + main.iters
  
    
    #initial settings
  
    if (!is.null(offset.coef)) {
      if (length(control$init[model$etamap$offsettheta]) !=
          length(offset.coef)) {
        stop("Invalid offset parameter vector offset.coef: ",
             "wrong number of parameters: expected ",
             length(control$init[model$etamap$offsettheta]),
             " got ", length(offset.coef), ".")}
      control$init[model$etamap$offsettheta] <- offset.coef
    }
    sigma2Samples <- array(NA, c(tot.iters,1, nchains))
    invTau2Samples <- array(NA, c(tot.iters, dim, nchains))
    lambda2Samples <- array(NA, c(tot.iters, dim, nchains))
    thetaSamples <- array(NA, c(tot.iters, dim, nchains))
    
    aux.iters<-aux.iters
    k <- 0
    
    
    
    control <- control.ergm(MCMC.burnin = aux.iters, 
                            MCMC.interval = 1,
                            MCMC.samplesize = 1)
    
    if (!is.null(control$init)) {
      if (length(control$init) != length(model$etamap$offsettheta)) {
        stop(paste("Invalid starting parameter vector control$init:",
                   "wrong number of parameters.", 
                   "If you are passing output from another ergm run as control$init,",
                   "in a model with curved terms, see help(enformulate.curved)."))
      }
    } else {control$init <- rep(NA, length(model$etamap$offsettheta))}
    
    proposal <- ergm_proposal(object = ~., constraints = ~.,
                              arguments = control$MCMC.prop.args, 
                              nw = y)
    
    S.prop <- diag(V.proposal, dim, dim)
    Theta <- array(NA, c(main.iters, dim, nchains))
    
    
    if (is.null(startVals)) {
      suppressMessages(mple <- coef(ergm(formula, estimate = "MPLE",
                                    verbose = FALSE,
                                    offset.coef = offset.coef)))
      # suppressMessages(mple <- ergm(formula, estimate = "MPLE",
      #                               verbose = FALSE,
      #                               offset.coef = offset.coef)$coef)
      theta <- matrix(mple + runif(dim * nchains, min = -0.1, max = 0.1), 
                      dim, 
                      nchains)
    } else {
      theta <- matrix(startVals + 
                        runif(dim * nchains, min = -0.1, max = 0.1), 
                      dim, 
                      nchains)
    }  
    
    theta[model$etamap$offsettheta,] <- offset.coef
    
    mean.prior <- rep(mean_prior,dim)
    if (is.null(mean_prior)){mean.prior<-mple
    nas_in_mean.prior=which(!is.finite(mean.prior))
    mean.prior[nas_in_mean.prior]=0}
    
    
    nas_in_theta=which(!is.finite(theta))
    theta[nas_in_theta]=theta[nas_in_theta-1]
    
    
    y0 <- simulate(formula,
                   coef        = theta[,1],
                   nsim        = 1,
                   control     = control.simulate(MCMC.burnin = 1, # !!!
                                                  MCMC.interval = 1),
                   return.args = "ergm_state")$object
    
    
    
    
    theta1 <- rep(NA, dim)
    
    message(" > Bergm Lasso Full Bayes MCMC start")
    clock.start <- Sys.time()
    for (j in 1:tot.iters) {
      cat("iteration",j,"\n")
      for (h in 1:nchains) {
        # cat("chain",h,"\n")
        
        invD <- diag(invTau2,dim)
        
        theta1 <- theta[, h] + 
          gamma * 
          apply(theta[, sample(seq(1, nchains)[-h], 2)], 1, diff) + 
          rmvnorm(1, sigma =   S.prop)[1,]
        
        
        delta <- ergm_MCMC_sample(y0,
                                  theta   = theta1,
                                  control = control)$stats[[1]][1,] - sy
        
        pr <- dmvnorm(rbind(theta1, theta[, h]), 
                      mean = mean.prior,
                      sigma = sigma.prior, 
                      log = TRUE)
        
        beta <- (theta[, h] - theta1) %*% delta + pr[1] - pr[2]
        
          
        if ((beta >= log(runif(1))) & (!is.na(beta))) theta[, h] <- theta1
        
        thetaSamples[j,,h]=theta[, h]
        
        
        #sample sigma.prior
        shape<-dim/2
        scale<-t(theta[, h])%*%invD%*%theta[, h]/2
        sigma2<-1/rgamma(1, shape=shape, rate=scale)
        # sigma2<-1/rgamma(1, shape=shape, rate=scale)
        # sigma2=rinvgamma(1,shape=2,rate=scale)
        
        # sigma.prior<-sigma2*(invD)
        # sigma.prior<-sigma2*solve(invD)
        sigma.prior<-sigma2*diag(1/invTau2,dim)
        sigma2Samples[j,,h]<-sigma2
        
        
        # sample tau2
        
        muPrime <- sqrt(lambda2* sigma2 /(theta[, h])^2)
        lambdaPrime <- lambda2
        invTau2 <- rep(0, dim)
        for (i in seq(dim)) {
          # invTau2[i] <- rinv.gaussian(1, muPrime[i], lambdaPrime[i])
          invTau2[i] <- rinvgauss(1, mean=muPrime[i], shape=lambdaPrime[i])
          
        }
        #invTau2 is 1/tau_j^2
        invTau2Samples[j,,h] <- invTau2
        
        #sample lambda2
        #the lambda2 has prior Gamma(shape=r00,rate=delta00)
        r00=r0
        delta00=delta0
        lambda2_shape<-r00+1
        #lambda2_scale is 2/tau_j^2
        lambda2_rate<-1/invTau2/2+delta00
        lambda2<-rgamma(dim,shape=lambda2_shape,rate=lambda2_rate)
        lambda2Samples[j,,h] <-lambda2
        
        # if (j %% 10 == 0) {
        #   low <- j - 9
        #   high <- j
        #   delta0 <-  r0/ (colMeans(lambda2Samples[low:high, ,h ]))
        # }

      }
      
      
      if (j > burn.in) Theta[j - burn.in, , ] <- theta
    }
    clock.end <- Sys.time()
    runtime <- difftime(clock.end, clock.start)
    FF <- mcmc(na.omit(apply(Theta, 2, cbind)),thin = thin)
    #acceptance rate
    AR <- round(1 - rejectionRate(FF)[1], 2)
    #ess is the effective sample size
    ess <- round(effectiveSize(FF), 0)
    names(ess) <- param_names(model)
    
    specs=param_names(model)
    # Theta1 = Theta[ , ,1]
    
    # 
    # sink("Bergm Lasso output.txt")
    # cat("Posterior Density Estimate for Model: ") 
    # formula
    # 
    # Theta1 <- as.mcmc(Theta1)
    # quantiles <- c(0.025, 0.25, 0.5, 0.75, 0.975)
    # 
    # statnames <- c("Mean", "SD", "Naive SE", "Time-series SE")
    # varstats <- matrix(nrow = nvar(Theta1), ncol = length(statnames), 
    #                    dimnames = list(varnames(Theta1), statnames))
    # 
    # Thetavar <- apply(Theta1, 2, var)
    # Thetatsvar <- apply(Theta1, 2, function(x) coda::spectrum0.ar(x)$spec)
    # varquant <- t(apply(Theta1, 2, quantile, quantiles))
    # 
    # varstats[, 1] <- apply(Theta1, 2, mean)
    # varstats[, 2] <- sqrt(Thetavar)
    # varstats[, 3] <- sqrt(Thetavar / niter(Theta1))
    # varstats[, 4] <- sqrt(Thetatsvar / niter(Theta1))
    # 
    # table1 <- drop(varstats)
    # table2 <- drop(varquant)
    # 
    # rNames <- paste("theta", seq(1, dim), " (", specs[seq(1, dim)], ")", sep = "")
    # 
    # rownames(table1) <- rownames(table2) <- rNames
    # print(table1); cat("\n"); print(table2)
    # cat("\n", "Acceptance rate:", AR, "\n", "\n", "\n")
    # 
    # 
    # sink() 
    # 
  
    
    out = list(formula = formula, specs = specs,
               dim = dim, Theta = FF, AR = AR,ess=ess, runtime=runtime)
    return(out)
    
  }
  
  # #setwd("~/Thesis/Codes/07122022 kapferer")
  # data("kapferer" )
  # net <-kapferer
  # 
  # 
  # Bergm.Lasso = BergmLasso(formula = net ~ edges + triangle , 
  #                        data = net,
  #                        sigma_prior_diag_value=100, 
  #                        lambda=5,
  #                        mean_prior=0 , 
  #                        invTau2=1/10 , 
  #                        main.iters=500,
  #                        burn.in =100, 
  #                        V.proposal=0.0025,
  #                        gamma=0.5, 
  #                        nchains=4 )
