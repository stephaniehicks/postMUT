#####################################################################
# pred.mat
# n.algos = Number of algorithms
# Function that sets up the disjoint categories for the for each of the algorithms
#####################################################################

pred.mat <- function(n.algos){
	if(n.algos == 2){ dat.mat <- expand.grid(c(0,1),c(0,1))  }
	if(n.algos == 3){ dat.mat <- expand.grid(c(0,1),c(0,1),c(0,1)) }
	if(n.algos == 4){ dat.mat <- expand.grid(c(0,1),c(0,1),c(0,1),c(0,1)) }

	return(dat.mat)
}
# pred.mat(2)



#####################################################################
# group.counts
# dat: Takes in a data set to calculate which group each set predictions falls into					
# n.algos: Number of algorithms
#####################################################################

group.counts <- function(dat,n.algos){
	if(n.algos == 2){
		e.grid <- expand.grid(c(0,1),c(0,1))
		n.r = nrow(dat)
		
		joint <- apply(e.grid,1,function(.row){
			count <- length(which(.row[1] == dat[,1] & .row[2] == dat[,2] ))
		})
	}
	
	if(n.algos == 3){
		e.grid <- expand.grid(c(0,1),c(0,1),c(0,1))
		n.r = nrow(dat)
		
		joint <- apply(e.grid,1,function(.row){
			count <- length(which(.row[1] == dat[,1] & .row[2] == dat[,2] & .row[3] == dat[,3] ))
		})
	}
	
	if(n.algos == 4){
		e.grid <- expand.grid(c(0,1),c(0,1),c(0,1),c(0,1))
		n.r = nrow(dat)
		
		joint <- apply(e.grid,1,function(.row){
			count <- length(which(.row[1] == dat[,1] & .row[2] == dat[,2] & .row[3] == dat[,3] & .row[4] == dat[,4] ))
		})
	}
	
	return(joint)
}		


#####################################################################
# sim.categories.xy
# n.size = Number of mutations to simulate
# n.algos = Number algorithms to make predictions for each of the mutations; the number of categories is determined by the number of algorithms
# lkd = IF TRUE, then sim.categories.xy returns the probability of being in each category p.t given input (a0, b0, p0); 
#       IF FALSE, then function returns simulated categories according to the CR model
#####################################################################

sim.categories.xy <- function(n.size, n.algos, input, lkd=FALSE){	
  n.cats = 2^n.algos # create the number of categories regardless of the number of parameters to estimate (e.g. length(a0) )
  
  xs <- pred.mat(n.algos)
  p.t = array(0,dim=n.cats) # probability of being in each category
  a0 <- input[[1]]; b0 <- input[[2]]; p0 <- input[[3]]

  if(length(a0) == 1){
    if(n.cats == 4){
      p.t = a0^xs[,1] * a0^xs[,2] * (1-a0)^(1-xs[,1]) * (1-a0)^(1-xs[,2]) * (1-p0) + 
            b0^xs[,1] * b0^xs[,2] * (1-b0)^(1-xs[,1]) * (1-b0)^(1-xs[,2]) * p0
      out <- sample(1:n.cats, n.size, prob=p.t,replace=TRUE) }
  
    if(n.cats == 8){		
      p.t = a0^xs[,1] * a0^xs[,2] * a0^xs[,3] * (1-a0)^(1-xs[,1]) * (1-a0)^(1-xs[,2]) * (1-a0)^(1-xs[,3]) * (1-p0) + 
            b0^xs[,1] * b0^xs[,2] * b0^xs[,3] * (1-b0)^(1-xs[,1]) * (1-b0)^(1-xs[,2]) * (1-b0)^(1-xs[,3]) * p0
      out <- sample(1:n.cats, n.size, prob=p.t,replace=TRUE) }

    if(n.cats == 16){
      p.t = a0^xs[,1] * a0^xs[,2] * a0^xs[,3] * a0^xs[,4] * (1-a0)^(1-xs[,1]) * (1-a0)^(1-xs[,2]) * (1-a0)^(1-xs[,3]) * (1-a0)^(1-xs[,4]) * (1-p0) + 
            b0^xs[,1] * b0^xs[,2] * b0^xs[,3] * b0^xs[,4] * (1-b0)^(1-xs[,1]) * (1-b0)^(1-xs[,2]) * (1-b0)^(1-xs[,3]) * (1-b0)^(1-xs[,4]) * p0
      out <- sample(1:n.cats,n.size, prob=p.t,replace=TRUE) }
  }


  if(length(a0) == 2){
	  p.t = a0[1]^xs[,1] * a0[2]^xs[,2] * (1-a0[1])^(1-xs[,1]) * (1-a0[2])^(1-xs[,2]) * (1-p0) + 
          b0[1]^xs[,1] * b0[2]^xs[,2] * (1-b0[1])^(1-xs[,1]) * (1-b0[2])^(1-xs[,2]) * p0
    out <- sample(1:n.cats,n.size, prob=p.t,replace=TRUE) 
  }		
		
  if(length(a0) == 3){
	  p.t = a0[1]^xs[,1] * a0[2]^xs[,2] * a0[3]^xs[,3] * (1-a0[1])^(1-xs[,1]) * (1-a0[2])^(1-xs[,2]) * (1-a0[3])^(1-xs[,3]) * (1-p0) + 
          b0[1]^xs[,1] * b0[2]^xs[,2] * b0[3]^xs[,3] * (1-b0[1])^(1-xs[,1]) * (1-b0[2])^(1-xs[,2]) * (1-b0[3])^(1-xs[,3]) * p0
	  out <- sample(1:n.cats,n.size, prob=p.t,replace=TRUE) 
	}
		
  if(length(a0) == 4){
    p.t = a0[1]^xs[,1] * a0[2]^xs[,2] * a0[3]^xs[,3] * a0[4]^xs[,4] * (1-a0[1])^(1-xs[,1]) * (1-a0[2])^(1-xs[,2]) * (1-a0[3])^(1-xs[,3]) * (1-a0[4])^(1-xs[,4]) * (1-p0) + 
          b0[1]^xs[,1] * b0[2]^xs[,2] * b0[3]^xs[,3] * b0[4]^xs[,4] * (1-b0[1])^(1-xs[,1]) * (1-b0[2])^(1-xs[,2]) * (1-b0[3])^(1-xs[,3]) * (1-b0[4])^(1-xs[,4]) * p0
    out <- sample(1:n.cats,n.size, prob=p.t,replace=TRUE) 
  }		
  
  if(lkd){
    return(p.t)
  } else { return(out) }
  
}        
# sim.categories.xy(N.size, n.algos = 4, true.xy)
# sim.categories.xy(N.size,N.algo,true.xy,lkd=TRUE)




#####################################################################
# sim.categories.xyz
# n.size = Number of mutations to simulate
# n.algos = Number algorithms to make predictions for each of the mutations; the number of categories is determined by the number of algorithms
# lkd = IF TRUE, then sim.categories.xy returns the probability of being in each category p.t given input (d0, e0, delta, gamm, p0); 
#       IF FALSE, then function returns simulated categories according to the CR model
#####################################################################

sim.categories.xyz <- function(n.size, n.algos, input, lkd=FALSE){  	
  n.cats = 2^n.algos # create the number of categories regardless of the number of parameters to estimate (e.g. length(d0) )
  
  xs <- pred.mat(n.algos)
  p.t = array(0,dim=n.cats) # probability of being in each category
  d0 <- input[[1]]; e0 <- input[[2]]; delta <- input[[3]]; gamm <- input[[4]]; p0 <- input[[5]]

  if(length(d0) == 1){
    if(n.cats == 4){
      p.t = d0^xs[,1] * d0^xs[,2] * (1-d0)^(1-xs[,1]) * (1-d0)^(1-xs[,2]) * delta * (1-p0) + 
            d0^xs[,1] * d0^xs[,2] * (1-d0)^(1-xs[,1]) * (1-d0)^(1-xs[,2]) * (1-gamm) * p0 + 
            e0^xs[,1] * e0^xs[,2] * (1-e0)^(1-xs[,1]) * (1-e0)^(1-xs[,2]) * (1-delta) * (1-p0) + 
            e0^xs[,1] * e0^xs[,2] * (1-e0)^(1-xs[,1]) * (1-e0)^(1-xs[,2]) * gamm * p0 
      out <- sample(1:n.cats, n.size, prob=p.t,replace=TRUE) }
    
    if(n.cats == 8){  
      p.t = d0^xs[,1] * d0^xs[,2] * d0^xs[,3] * (1-d0)^(1-xs[,1]) * (1-d0)^(1-xs[,2]) * (1-d0)^(1-xs[,3]) * delta * (1-p0) + 
            d0^xs[,1] * d0^xs[,2] * d0^xs[,3] * (1-d0)^(1-xs[,1]) * (1-d0)^(1-xs[,2]) * (1-d0)^(1-xs[,3]) * (1-gamm) * p0 + 
            e0^xs[,1] * e0^xs[,2] * e0^xs[,3] * (1-e0)^(1-xs[,1]) * (1-e0)^(1-xs[,2]) * (1-e0)^(1-xs[,3]) * (1-delta) * (1-p0) + 
            e0^xs[,1] * e0^xs[,2] * e0^xs[,3] * (1-e0)^(1-xs[,1]) * (1-e0)^(1-xs[,2]) * (1-e0)^(1-xs[,3]) * gamm * p0
      out <- sample(1:n.cats,n.size, prob=p.t,replace=TRUE) }
      
    if(n.cats == 16){
      p.t = d0^xs[,1] * d0^xs[,2] * d0^xs[,3] * d0^xs[,4] * (1-d0)^(1-xs[,1]) * (1-d0)^(1-xs[,2]) * (1-d0)^(1-xs[,3]) * (1-d0)^(1-xs[,4]) * delta * (1-p0) + 
            d0^xs[,1] * d0^xs[,2] * d0^xs[,3] * d0^xs[,4] * (1-d0)^(1-xs[,1]) * (1-d0)^(1-xs[,2]) * (1-d0)^(1-xs[,3]) * (1-d0)^(1-xs[,4]) * (1-gamm) * p0 + 
            e0^xs[,1] * e0^xs[,2] * e0^xs[,3] * e0^xs[,4] * (1-e0)^(1-xs[,1]) * (1-e0)^(1-xs[,2]) * (1-e0)^(1-xs[,3]) * (1-e0)^(1-xs[,4]) * (1-delta) * (1-p0) + 
            e0^xs[,1] * e0^xs[,2] * e0^xs[,3] * e0^xs[,4] * (1-e0)^(1-xs[,1]) * (1-e0)^(1-xs[,2]) * (1-e0)^(1-xs[,3]) * (1-e0)^(1-xs[,4]) * gamm * p0
      out <- sample(1:n.cats,n.size, prob=p.t,replace=TRUE) }
  }
   
  if(length(d0) == 2){
      p.t = d0[1]^xs[,1] * d0[2]^xs[,2] * (1-d0[1])^(1-xs[,1]) * (1-d0[2])^(1-xs[,2]) * delta * (1-p0) + 
            d0[1]^xs[,1] * d0[2]^xs[,2] * (1-d0[1])^(1-xs[,1]) * (1-d0[2])^(1-xs[,2]) * (1-gamm) * p0 + 
            e0[1]^xs[,1] * e0[2]^xs[,2] * (1-e0[1])^(1-xs[,1]) * (1-e0[2])^(1-xs[,2]) * (1-delta) * (1-p0) + 
            e0[1]^xs[,1] * e0[2]^xs[,2] * (1-e0[1])^(1-xs[,1]) * (1-e0[2])^(1-xs[,2]) * gamm * p0 
      out <- sample(1:n.cats, n.size, prob=p.t,replace=TRUE)
  }    
    
  if(length(d0) == 3){	
    p.t = d0[1]^xs[,1] * d0[2]^xs[,2] * d0[3]^xs[,3] * (1-d0[1])^(1-xs[,1]) * (1-d0[2])^(1-xs[,2]) * (1-d0[3])^(1-xs[,3]) * delta * (1-p0) + 
          d0[1]^xs[,1] * d0[2]^xs[,2] * d0[3]^xs[,3] * (1-d0[1])^(1-xs[,1]) * (1-d0[2])^(1-xs[,2]) * (1-d0[3])^(1-xs[,3]) * (1-gamm) * p0 + 
          e0[1]^xs[,1] * e0[2]^xs[,2] * e0[3]^xs[,3] * (1-e0[1])^(1-xs[,1]) * (1-e0[2])^(1-xs[,2]) * (1-e0[3])^(1-xs[,3]) * (1-delta) * (1-p0) + 
          e0[1]^xs[,1] * e0[2]^xs[,2] * e0[3]^xs[,3] * (1-e0[1])^(1-xs[,1]) * (1-e0[2])^(1-xs[,2]) * (1-e0[3])^(1-xs[,3]) * gamm * p0
    out <- sample(1:n.cats,n.size, prob=p.t,replace=TRUE) 
  }
    
  if(length(d0) == 4){
    p.t = d0[1]^xs[,1] * d0[2]^xs[,2] * d0[3]^xs[,3] * d0[4]^xs[,4] * (1-d0[1])^(1-xs[,1]) * (1-d0[2])^(1-xs[,2]) * (1-d0[3])^(1-xs[,3]) * (1-d0[4])^(1-xs[,4]) * delta * (1-p0) + 
          d0[1]^xs[,1] * d0[2]^xs[,2] * d0[3]^xs[,3] * d0[4]^xs[,4] * (1-d0[1])^(1-xs[,1]) * (1-d0[2])^(1-xs[,2]) * (1-d0[3])^(1-xs[,3]) * (1-d0[4])^(1-xs[,4]) * (1-gamm) * p0 + 
          e0[1]^xs[,1] * e0[2]^xs[,2] * e0[3]^xs[,3] * e0[4]^xs[,4] * (1-e0[1])^(1-xs[,1]) * (1-e0[2])^(1-xs[,2]) * (1-e0[3])^(1-xs[,3]) * (1-e0[4])^(1-xs[,4]) * (1-delta) * (1-p0) + 
          e0[1]^xs[,1] * e0[2]^xs[,2] * e0[3]^xs[,3] * e0[4]^xs[,4] * (1-e0[1])^(1-xs[,1]) * (1-e0[2])^(1-xs[,2]) * (1-e0[3])^(1-xs[,3]) * (1-e0[4])^(1-xs[,4]) * gamm * p0
    out <- sample(1:n.cats,n.size, prob=p.t,replace=TRUE) 
  }	

  if(lkd){
  return(p.t)
  } else { return(out) }

}        
# sim.categories.xyz(N.size,N.algo,true.xyz)
# sim.categories.xyz(N.size,N.algo,true.xyz,lkd=TRUE)


#####################################################################
# posterior.xy
# obs_pred = Observed functional predictions
# input_xy = The estimates from capture.em.xy
# n.algos = Number algorithms to make predictions for each of the mutations; the number of categories is determined by the number of algorithms
#####################################################################

posterior.xy <- function(obs_pred, input, n.algos, prob_D = TRUE){  
  n.cats = 2^n.algos # create the number of categories regardless of the number of parameters to estimate (e.g. length(a0) )
  
  xs <- obs_pred
  p.t = array(0,dim=n.cats) # posterior probability
  a0 <- input[1:n.algos]; b0 <- input[(n.algos+1):(2*n.algos)]; p0 <- input[(2*n.algos)+1]
  
  if(length(a0) == 1){
    if(n.cats == 4){
      p.t = b0^xs[,1] * b0^xs[,2] * (1-b0)^(1-xs[,1]) * (1-b0)^(1-xs[,2]) * p0 / 
        ( a0^xs[,1] * a0^xs[,2] * (1-a0)^(1-xs[,1]) * (1-a0)^(1-xs[,2]) * (1-p0) + b0^xs[,1] * b0^xs[,2] * (1-b0)^(1-xs[,1]) * (1-b0)^(1-xs[,2]) * p0 ) 
    }
    
    if(n.cats == 8){
      p.t = b0^xs[,1] * b0^xs[,2] * b0^xs[,3] * (1-b0)^(1-xs[,1]) * (1-b0)^(1-xs[,2]) * (1-b0)^(1-xs[,3]) * p0 / 
        ( a0^xs[,1] * a0^xs[,2] * a0^xs[,3] * (1-a0)^(1-xs[,1]) * (1-a0)^(1-xs[,2]) * (1-a0)^(1-xs[,3]) * (1-p0) + b0^xs[,1] * b0^xs[,2] * b0^xs[,3] * (1-b0)^(1-xs[,1]) * (1-b0)^(1-xs[,2]) * (1-b0)^(1-xs[,3]) * p0 ) 
    }
    
    if(n.cats == 16){
      p.t = b0^xs[,1] * b0^xs[,2] * b0^xs[,3] * b0^xs[,4] * (1-b0)^(1-xs[,1]) * (1-b0)^(1-xs[,2]) * (1-b0)^(1-xs[,3]) * (1-b0)^(1-xs[,4]) * p0 / 
        ( a0^xs[,1] * a0^xs[,2] * a0^xs[,3] * a0^xs[,4] * (1-a0)^(1-xs[,1]) * (1-a0)^(1-xs[,2]) * (1-a0)^(1-xs[,3]) * (1-a0)^(1-xs[,4]) * (1-p0) + b0^xs[,1] * b0^xs[,2] * b0^xs[,3] * b0^xs[,4] * (1-b0)^(1-xs[,1]) * (1-b0)^(1-xs[,2]) * (1-b0)^(1-xs[,3]) * (1-b0)^(1-xs[,4]) * p0 ) 
    } 
  }  
  
  if(length(a0) == 2){
    p.t = b0[1]^xs[,1] * b0[2]^xs[,2] * (1-b0[1])^(1-xs[,1]) * (1-b0[2])^(1-xs[,2]) * p0 / 
      ( a0[1]^xs[,1] * a0[2]^xs[,2] * (1-a0[1])^(1-xs[,1]) * (1-a0[2])^(1-xs[,2]) * (1-p0) + b0[1]^xs[,1] * b0[2]^xs[,2] * (1-b0[1])^(1-xs[,1]) * (1-b0[2])^(1-xs[,2]) * p0 ) 
  }  	
  
  if(length(a0) == 3){
    p.t = b0[1]^xs[,1] * b0[2]^xs[,2] * b0[3]^xs[,3] * (1-b0[1])^(1-xs[,1]) * (1-b0[2])^(1-xs[,2]) * (1-b0[3])^(1-xs[,3]) * p0 / 
      ( a0[1]^xs[,1] * a0[2]^xs[,2] * a0[3]^xs[,3] * (1-a0[1])^(1-xs[,1]) * (1-a0[2])^(1-xs[,2]) * (1-a0[3])^(1-xs[,3]) * (1-p0) + b0[1]^xs[,1] * b0[2]^xs[,2] * b0[3]^xs[,3] * (1-b0[1])^(1-xs[,1]) * (1-b0[2])^(1-xs[,2]) * (1-b0[3])^(1-xs[,3]) * p0 )
  }
  
  if(length(a0) == 4){
    p.t = b0[1]^xs[,1] * b0[2]^xs[,2] * b0[3]^xs[,3] * b0[4]^xs[,4] * (1-b0[1])^(1-xs[,1]) * (1-b0[2])^(1-xs[,2]) * (1-b0[3])^(1-xs[,3]) * (1-b0[4])^(1-xs[,4]) * p0 / 
      ( a0[1]^xs[,1] * a0[2]^xs[,2] * a0[3]^xs[,3] * a0[4]^xs[,4] * (1-a0[1])^(1-xs[,1]) * (1-a0[2])^(1-xs[,2]) * (1-a0[3])^(1-xs[,3]) * (1-a0[4])^(1-xs[,4]) * (1-p0) + b0[1]^xs[,1] * b0[2]^xs[,2] * b0[3]^xs[,3] * b0[4]^xs[,4] * (1-b0[1])^(1-xs[,1]) * (1-b0[2])^(1-xs[,2]) * (1-b0[3])^(1-xs[,3]) * (1-b0[4])^(1-xs[,4]) * p0 )
  }		
  
  if(prob_D){
    return(p.t)
  } else { return(1-p.t) }
  
}        
# posterior.xy(pred.mat(4), output.xy, n.algos = 4) # Posterior Pr(D)
# posterior.xy(pred.mat(4), output.xy, n.algos = 4, prob_D = FALSE) # Posterior Pr(N)


#####################################################################
# posterior.xyz
# obs_pred = Observed functional predictions
# input = The estimates from capture.em.xyz
# n.algos = Number algorithms to make predictions for each of the mutations; the number of categories is determined by the number of algorithms
#####################################################################

posterior.xyz <- function(obs_pred, input, n.algos, prob_D = TRUE, de_1 = FALSE){  
  n.cats = 2^n.algos # create the number of categories regardless of the number of parameters to estimate (e.g. length(d0) )
  
  xs <- obs_pred
  p.t.10 = array(0,dim=n.cats) # (Z = 1, Y = 0)
  p.t.01 = array(0,dim=n.cats) # (Z = 0, Y = 1)
  p.t.11 = array(0,dim=n.cats) # (Z = 1, Y = 1) posterior probability 
  d0 <- input[1:n.algos]; e0 <- input[(n.algos+1):(2*n.algos)]; delta <- input[(2*n.algos)+1]; gamm <- input[(2*n.algos)+2]; p0 <- input[(2*n.algos)+3]
  
  # Sanity check
  if(de_1){ stopifnot(length(d0) == 1)}
  if(length(d0) != 1){ stopifnot(length(d0) == n.algos) }
  # stopifnot(length(groups) == n.cats)
  
  if(length(d0) == 1){
    if(n.cats == 4){
      p.t.01 =  d0^xs[,1] * d0^xs[,2] * (1-d0)^(1-xs[,1]) * (1-d0)^(1-xs[,2]) * (1-gamm) * p0 / 
        ( d0^xs[,1] * d0^xs[,2] * (1-d0)^(1-xs[,1]) * (1-d0)^(1-xs[,2]) * ( delta * (1-p0) + (1-gamm) * p0 )  + 
        e0^xs[,1] * e0^xs[,2] * (1-e0)^(1-xs[,1]) * (1-e0)^(1-xs[,2]) * ( (1-delta) * (1-p0) + gamm * p0 ) ) 
      p.t.11 =  e0^xs[,1] * e0^xs[,2] * (1-e0)^(1-xs[,1]) * (1-e0)^(1-xs[,2]) * gamm * p0 / 
        ( d0^xs[,1] * d0^xs[,2] * (1-d0)^(1-xs[,1]) * (1-d0)^(1-xs[,2]) * ( delta * (1-p0) + (1-gamm) * p0 )  +
        e0^xs[,1] * e0^xs[,2] * (1-e0)^(1-xs[,1]) * (1-e0)^(1-xs[,2]) * ( (1-delta) * (1-p0) + gamm * p0 ) )
    }
    
    if(n.cats == 8){  
      p.t.01 =  d0^xs[,1] * d0^xs[,2] * d0^xs[,3] * (1-d0)^(1-xs[,1]) * (1-d0)^(1-xs[,2]) * (1-d0)^(1-xs[,3]) * (1-gamm) * p0 / 
        ( d0^xs[,1] * d0^xs[,2] * d0^xs[,3] * (1-d0)^(1-xs[,1]) * (1-d0)^(1-xs[,2]) * (1-d0)^(1-xs[,3]) * ( delta * (1-p0) + (1-gamm) * p0 )  +
        e0^xs[,1] * e0^xs[,2] * e0^xs[,3] * (1-e0)^(1-xs[,1]) * (1-e0)^(1-xs[,2]) * (1-e0)^(1-xs[,3]) * ( (1-delta) * (1-p0) + gamm * p0 ) )
      p.t.11 =  e0^xs[,1] * e0^xs[,2] * e0^xs[,3] * (1-e0)^(1-xs[,1]) * (1-e0)^(1-xs[,2]) * (1-e0)^(1-xs[,3]) * gamm * p0 / 
        ( d0^xs[,1] * d0^xs[,2] * d0^xs[,3] * (1-d0)^(1-xs[,1]) * (1-d0)^(1-xs[,2]) * (1-d0)^(1-xs[,3]) * ( delta * (1-p0) + (1-gamm) * p0 )  +
        e0^xs[,1] * e0^xs[,2] * e0^xs[,3] * (1-e0)^(1-xs[,1]) * (1-e0)^(1-xs[,2]) * (1-e0)^(1-xs[,3]) * ( (1-delta) * (1-p0) + gamm * p0 ) )
    }
    
    if(n.cats == 16){
      p.t.01 =  d0^xs[,1] * d0^xs[,2] * d0^xs[,3] * d0^xs[,4] * (1-d0)^(1-xs[,1]) * (1-d0)^(1-xs[,2]) * (1-d0)^(1-xs[,3]) * (1-d0)^(1-xs[,4]) * (1-gamm) * p0 / 
        ( d0^xs[,1] * d0^xs[,2] * d0^xs[,3] * d0^xs[,4] * (1-d0)^(1-xs[,1]) * (1-d0)^(1-xs[,2]) * (1-d0)^(1-xs[,3]) * (1-d0)^(1-xs[,4]) * ( delta * (1-p0) + (1-gamm) * p0 )  +
        e0^xs[,1] * e0^xs[,2] * e0^xs[,3] * e0^xs[,4] * (1-e0)^(1-xs[,1]) * (1-e0)^(1-xs[,2]) * (1-e0)^(1-xs[,3]) * (1-e0)^(1-xs[,4]) * ( (1-delta) * (1-p0) + gamm * p0 ) )
      p.t.11 =  e0^xs[,1] * e0^xs[,2] * e0^xs[,3] * e0^xs[,4] * (1-e0)^(1-xs[,1]) * (1-e0)^(1-xs[,2]) * (1-e0)^(1-xs[,3]) * (1-e0)^(1-xs[,4]) * gamm * p0 / 
        ( d0^xs[,1] * d0^xs[,2] * d0^xs[,3] * d0^xs[,4] * (1-d0)^(1-xs[,1]) * (1-d0)^(1-xs[,2]) * (1-d0)^(1-xs[,3]) * (1-d0)^(1-xs[,4]) * ( delta * (1-p0) + (1-gamm) * p0 )  +
        e0^xs[,1] * e0^xs[,2] * e0^xs[,3] * e0^xs[,4] * (1-e0)^(1-xs[,1]) * (1-e0)^(1-xs[,2]) * (1-e0)^(1-xs[,3]) * (1-e0)^(1-xs[,4]) * ( (1-delta) * (1-p0) + gamm * p0 ) )
    }
  }
  
  if(length(d0) == 2){
    p.t.01 =  d0[1]^xs[,1] * d0[2]^xs[,2] * (1-d0[1])^(1-xs[,1]) * (1-d0[2])^(1-xs[,2]) * (1-gamm) * p0 /
      ( d0[1]^xs[,1] * d0[2]^xs[,2] * (1-d0[1])^(1-xs[,1]) * (1-d0[2])^(1-xs[,2]) * ( delta * (1-p0) + (1-gamm) * p0 )  +
      e0[1]^xs[,1] * e0[2]^xs[,2] * (1-e0[1])^(1-xs[,1]) * (1-e0[2])^(1-xs[,2]) * ( (1-delta) * (1-p0) + gamm * p0 ) )
    p.t.11 =  e0[1]^xs[,1] * e0[2]^xs[,2] * (1-e0[1])^(1-xs[,1]) * (1-e0[2])^(1-xs[,2]) * gamm * p0 /
      ( d0[1]^xs[,1] * d0[2]^xs[,2] * (1-d0[1])^(1-xs[,1]) * (1-d0[2])^(1-xs[,2]) * ( delta * (1-p0) + (1-gamm) * p0 )  +
      e0[1]^xs[,1] * e0[2]^xs[,2] * (1-e0[1])^(1-xs[,1]) * (1-e0[2])^(1-xs[,2]) * ( (1-delta) * (1-p0) + gamm * p0 ) )
  }    
  
  if(length(d0) == 3){ 
    p.t.01 =  d0[1]^xs[,1] * d0[2]^xs[,2] * d0[3]^xs[,3] * (1-d0[1])^(1-xs[,1]) * (1-d0[2])^(1-xs[,2]) * (1-d0[3])^(1-xs[,3]) * (1-gamm) * p0 / 
      ( d0[1]^xs[,1] * d0[2]^xs[,2] * d0[3]^xs[,3] * (1-d0[1])^(1-xs[,1]) * (1-d0[2])^(1-xs[,2]) * (1-d0[3])^(1-xs[,3]) * ( delta * (1-p0) + (1-gamm) * p0 )  +
      e0[1]^xs[,1] * e0[2]^xs[,2] * e0[3]^xs[,3] * (1-e0[1])^(1-xs[,1]) * (1-e0[2])^(1-xs[,2]) * (1-e0[3])^(1-xs[,3]) * ( (1-delta) * (1-p0) + gamm * p0 ) )
    p.t.11 =  e0[1]^xs[,1] * e0[2]^xs[,2] * e0[3]^xs[,3] * (1-e0[1])^(1-xs[,1]) * (1-e0[2])^(1-xs[,2]) * (1-e0[3])^(1-xs[,3]) * gamm * p0 / 
      ( d0[1]^xs[,1] * d0[2]^xs[,2] * d0[3]^xs[,3] * (1-d0[1])^(1-xs[,1]) * (1-d0[2])^(1-xs[,2]) * (1-d0[3])^(1-xs[,3]) * ( delta * (1-p0) + (1-gamm) * p0 )  +
      e0[1]^xs[,1] * e0[2]^xs[,2] * e0[3]^xs[,3] * (1-e0[1])^(1-xs[,1]) * (1-e0[2])^(1-xs[,2]) * (1-e0[3])^(1-xs[,3]) * ( (1-delta) * (1-p0) + gamm * p0 ) )
  }
  
  if(length(d0) == 4){
    p.t.01 =  d0[1]^xs[,1] * d0[2]^xs[,2] * d0[3]^xs[,3] * d0[4]^xs[,4] * (1-d0[1])^(1-xs[,1]) * (1-d0[2])^(1-xs[,2]) * (1-d0[3])^(1-xs[,3]) * (1-d0[4])^(1-xs[,4]) * (1-gamm) * p0 / 
      ( d0[1]^xs[,1] * d0[2]^xs[,2] * d0[3]^xs[,3] * d0[4]^xs[,4] * (1-d0[1])^(1-xs[,1]) * (1-d0[2])^(1-xs[,2]) * (1-d0[3])^(1-xs[,3]) * (1-d0[4])^(1-xs[,4]) * ( delta * (1-p0) + (1-gamm) * p0 )  +
      e0[1]^xs[,1] * e0[2]^xs[,2] * e0[3]^xs[,3] * e0[4]^xs[,4] * (1-e0[1])^(1-xs[,1]) * (1-e0[2])^(1-xs[,2]) * (1-e0[3])^(1-xs[,3]) * (1-e0[4])^(1-xs[,4]) * ( (1-delta) * (1-p0) + gamm * p0 ) )
    p.t.11 =  e0[1]^xs[,1] * e0[2]^xs[,2] * e0[3]^xs[,3] * e0[4]^xs[,4] * (1-e0[1])^(1-xs[,1]) * (1-e0[2])^(1-xs[,2]) * (1-e0[3])^(1-xs[,3]) * (1-e0[4])^(1-xs[,4]) * gamm * p0 / 
      ( d0[1]^xs[,1] * d0[2]^xs[,2] * d0[3]^xs[,3] * d0[4]^xs[,4] * (1-d0[1])^(1-xs[,1]) * (1-d0[2])^(1-xs[,2]) * (1-d0[3])^(1-xs[,3]) * (1-d0[4])^(1-xs[,4]) * ( delta * (1-p0) + (1-gamm) * p0 )  +
      e0[1]^xs[,1] * e0[2]^xs[,2] * e0[3]^xs[,3] * e0[4]^xs[,4] * (1-e0[1])^(1-xs[,1]) * (1-e0[2])^(1-xs[,2]) * (1-e0[3])^(1-xs[,3]) * (1-e0[4])^(1-xs[,4]) * ( (1-delta) * (1-p0) + gamm * p0 ) )
  }
  
  
  if(prob_D){
    return(p.t.01 + p.t.11)
  } else { return(1 - p.t.01 - p.t.11) }
  
}        
# posterior.xyz(pred.mat(4), output.xy, n.algos = 4) # Posterior Pr(D)
# posterior.xyz(pred.mat(4), output.xy, n.algos = 4, prob_D = FALSE) # Posterior Pr(N)





#####################################################################
# capture.em.step.xy
# groups = the observed data (i.e. which category each mutation fell into)
# input = input parameters for the expectation step (a0, b0)
# n.algos = Number algorithms to make predictions for each of the mutations; the number of categories is determined by the number of algorithms
# ab_1 = IF TRUE, then the EM algorithm is estimated assuming the a0, b0, p0 is the same for all the algorithm
#        IF FALSE, then the EM algorithm is estimated assuming each algorithm has its own a0, b0
#####################################################################

capture.em.step.xy <- function(groups, input, n.algos, ab_1 = FALSE, expect = FALSE){	
  a0 <- input[[1]]; b0 <- input[[2]]; p0 <- input[[3]]
  n.cats = 2^n.algos # create the number of categories regardless of the number of parameters to estimate (e.g. length(a0) )
  xs <- pred.mat(n.algos)
  
  # Sanity check
  if(ab_1){ stopifnot(length(a0) == 1)}
  if(length(a0) != 1){ stopifnot(length(a0) == n.algos) }
  stopifnot(length(groups) == n.cats)
  
  ####### Expectation Step
  p.t = array(0,dim=n.cats)
  
  if(length(a0) == 1){
      if(n.cats == 4){
        p.t = b0^xs[,1] * b0^xs[,2] * (1-b0)^(1-xs[,1]) * (1-b0)^(1-xs[,2]) * p0 / 
              ( a0^xs[,1] * a0^xs[,2] * (1-a0)^(1-xs[,1]) * (1-a0)^(1-xs[,2]) * (1-p0) + b0^xs[,1] * b0^xs[,2] * (1-b0)^(1-xs[,1]) * (1-b0)^(1-xs[,2]) * p0 ) 
        }
    
      if(n.cats == 8){
        p.t = b0^xs[,1] * b0^xs[,2] * b0^xs[,3] * (1-b0)^(1-xs[,1]) * (1-b0)^(1-xs[,2]) * (1-b0)^(1-xs[,3]) * p0 / 
              ( a0^xs[,1] * a0^xs[,2] * a0^xs[,3] * (1-a0)^(1-xs[,1]) * (1-a0)^(1-xs[,2]) * (1-a0)^(1-xs[,3]) * (1-p0) + b0^xs[,1] * b0^xs[,2] * b0^xs[,3] * (1-b0)^(1-xs[,1]) * (1-b0)^(1-xs[,2]) * (1-b0)^(1-xs[,3]) * p0 ) 
        }
    
      if(n.cats == 16){
        p.t = b0^xs[,1] * b0^xs[,2] * b0^xs[,3] * b0^xs[,4] * (1-b0)^(1-xs[,1]) * (1-b0)^(1-xs[,2]) * (1-b0)^(1-xs[,3]) * (1-b0)^(1-xs[,4]) * p0 / 
              ( a0^xs[,1] * a0^xs[,2] * a0^xs[,3] * a0^xs[,4] * (1-a0)^(1-xs[,1]) * (1-a0)^(1-xs[,2]) * (1-a0)^(1-xs[,3]) * (1-a0)^(1-xs[,4]) * (1-p0) + b0^xs[,1] * b0^xs[,2] * b0^xs[,3] * b0^xs[,4] * (1-b0)^(1-xs[,1]) * (1-b0)^(1-xs[,2]) * (1-b0)^(1-xs[,3]) * (1-b0)^(1-xs[,4]) * p0 ) 
        } 
  }  

  if(length(a0) == 2){
    p.t = b0[1]^xs[,1] * b0[2]^xs[,2] * (1-b0[1])^(1-xs[,1]) * (1-b0[2])^(1-xs[,2]) * p0 / 
          ( a0[1]^xs[,1] * a0[2]^xs[,2] * (1-a0[1])^(1-xs[,1]) * (1-a0[2])^(1-xs[,2]) * (1-p0) + b0[1]^xs[,1] * b0[2]^xs[,2] * (1-b0[1])^(1-xs[,1]) * (1-b0[2])^(1-xs[,2]) * p0 ) 
  }		
  
  if(length(a0) == 3){
    p.t = b0[1]^xs[,1] * b0[2]^xs[,2] * b0[3]^xs[,3] * (1-b0[1])^(1-xs[,1]) * (1-b0[2])^(1-xs[,2]) * (1-b0[3])^(1-xs[,3]) * p0 / 
          ( a0[1]^xs[,1] * a0[2]^xs[,2] * a0[3]^xs[,3] * (1-a0[1])^(1-xs[,1]) * (1-a0[2])^(1-xs[,2]) * (1-a0[3])^(1-xs[,3]) * (1-p0) + b0[1]^xs[,1] * b0[2]^xs[,2] * b0[3]^xs[,3] * (1-b0[1])^(1-xs[,1]) * (1-b0[2])^(1-xs[,2]) * (1-b0[3])^(1-xs[,3]) * p0 )
  }
  
  if(length(a0) == 4){
    p.t = b0[1]^xs[,1] * b0[2]^xs[,2] * b0[3]^xs[,3] * b0[4]^xs[,4] * (1-b0[1])^(1-xs[,1]) * (1-b0[2])^(1-xs[,2]) * (1-b0[3])^(1-xs[,3]) * (1-b0[4])^(1-xs[,4]) * p0 / 
          ( a0[1]^xs[,1] * a0[2]^xs[,2] * a0[3]^xs[,3] * a0[4]^xs[,4] * (1-a0[1])^(1-xs[,1]) * (1-a0[2])^(1-xs[,2]) * (1-a0[3])^(1-xs[,3]) * (1-a0[4])^(1-xs[,4]) * (1-p0) + b0[1]^xs[,1] * b0[2]^xs[,2] * b0[3]^xs[,3] * b0[4]^xs[,4] * (1-b0[1])^(1-xs[,1]) * (1-b0[2])^(1-xs[,2]) * (1-b0[3])^(1-xs[,3]) * (1-b0[4])^(1-xs[,4]) * p0 )
  }		
  
  ####### Maximization Step
  # Estimate p_temp
  p_temp <- sum(p.t*groups) / sum(groups)
  
  # Estimate a_temp, b_temp
  if(ab_1) { 
    a_temp = b_temp = array(0,dim = c(1,1))
    row.algo <- rowSums(pred.mat(n.algos))
    a_temp <- sum((1-p.t)*groups*row.algo) / sum((1-p.t)*groups*ncol(pred.mat(n.algos))) 
    b_temp <- sum(p.t*groups*row.algo) / sum(p.t*groups*ncol(pred.mat(n.algos)))
    
  } else { 
    a_temp = b_temp = array(0,dim = c(n.algos,1)) 
    for(i in 1:n.algos){
      ind <- which(pred.mat(n.algos)[,i] == 1)
      a_temp[i] <- sum((1-p.t[ind])*groups[ind]) / sum((1-p.t)*groups) 
      b_temp[i] <- sum(p.t[ind]*groups[ind]) / sum(p.t*groups)
    }
  }

  if(expect){ return(p.t)
    } else { list(a_temp, b_temp, p_temp) }
  

}		
# capture.em.step.xy(group, list(c(.10,.10,.10,.10),c(.70,.70,.70,.70),.3), n.algos = N.algo, ab_1 = FALSE) # should work
# capture.em.step.xy(group, list(c(.10,.10,.10),c(.70,.70,.70),.3), n.algos = 3, ab_1 = FALSE) # only works if the number of categories in group data is 2^N.algo
  
# capture.em.step.xy(group, list(c(.10,.10,.10,.10),c(.70,.70,.70,.70),.3), n.algos = N.algo, ab_1 = FALSE) # should work


#####################################################################
# capture.em.step.xyz
# groups = the observed data (i.e. which category each mutation fell into)
# input = input parameters for the expectation step (d0, e0, delta, gamm, p0)
# n.algos = Number algorithms to make predictions for each of the mutations; the number of categories is determined by the number of algorithms
# de_1 = IF TRUE, then the EM algorithm is estimated assuming the d0, e0 is the same for all the algorithm
#        IF FALSE, then the EM algorithm is estimated assuming each algorithm has its own d0, e0
#####################################################################

# capture.em.step.xyz is a function that computes the expectation and maximization step
capture.em.step.xyz <- function(groups, input, n.algos, de_1 = FALSE, expect = FALSE){  	
  d0 <- input[[1]]; e0 <- input[[2]]; delta <- input[[3]]; gamm <- input[[4]]; p0 <- input[[5]]
  n.cats = 2^n.algos # create the number of categories regardless of the number of parameters to estimate (e.g. length(a0) )
  xs <- pred.mat(n.algos)
  
  # Sanity check
  if(de_1){ stopifnot(length(d0) == 1)}
  if(length(d0) != 1){ stopifnot(length(d0) == n.algos) }
  stopifnot(length(groups) == n.cats)
  
  ####### Expectation Step
  p.t.10 = p.t.01 = p.t.11 = array(0,dim=n.cats)
  if(length(d0) == 1){
    if(n.cats == 4){
      p.t.01 =  d0^xs[,1] * d0^xs[,2] * (1-d0)^(1-xs[,1]) * (1-d0)^(1-xs[,2]) * (1-gamm) * p0 / 
              ( d0^xs[,1] * d0^xs[,2] * (1-d0)^(1-xs[,1]) * (1-d0)^(1-xs[,2]) * ( delta * (1-p0) + (1-gamm) * p0 )  + 
                e0^xs[,1] * e0^xs[,2] * (1-e0)^(1-xs[,1]) * (1-e0)^(1-xs[,2]) * ( (1-delta) * (1-p0) + gamm * p0 ) ) 
      p.t.10 =  e0^xs[,1] * e0^xs[,2] * (1-e0)^(1-xs[,1]) * (1-e0)^(1-xs[,2]) * (1-delta) * (1-p0) / 
              ( d0^xs[,1] * d0^xs[,2] * (1-d0)^(1-xs[,1]) * (1-d0)^(1-xs[,2]) * ( delta * (1-p0) + (1-gamm) * p0 )  +
                e0^xs[,1] * e0^xs[,2] * (1-e0)^(1-xs[,1]) * (1-e0)^(1-xs[,2]) * ( (1-delta) * (1-p0) + gamm * p0 ) )
      p.t.11 =  e0^xs[,1] * e0^xs[,2] * (1-e0)^(1-xs[,1]) * (1-e0)^(1-xs[,2]) * gamm * p0 / 
              ( d0^xs[,1] * d0^xs[,2] * (1-d0)^(1-xs[,1]) * (1-d0)^(1-xs[,2]) * ( delta * (1-p0) + (1-gamm) * p0 )  +
                e0^xs[,1] * e0^xs[,2] * (1-e0)^(1-xs[,1]) * (1-e0)^(1-xs[,2]) * ( (1-delta) * (1-p0) + gamm * p0 ) )
      }
    
    if(n.cats == 8){  
      p.t.01 =  d0^xs[,1] * d0^xs[,2] * d0^xs[,3] * (1-d0)^(1-xs[,1]) * (1-d0)^(1-xs[,2]) * (1-d0)^(1-xs[,3]) * (1-gamm) * p0 / 
              ( d0^xs[,1] * d0^xs[,2] * d0^xs[,3] * (1-d0)^(1-xs[,1]) * (1-d0)^(1-xs[,2]) * (1-d0)^(1-xs[,3]) * ( delta * (1-p0) + (1-gamm) * p0 )  +
                e0^xs[,1] * e0^xs[,2] * e0^xs[,3] * (1-e0)^(1-xs[,1]) * (1-e0)^(1-xs[,2]) * (1-e0)^(1-xs[,3]) * ( (1-delta) * (1-p0) + gamm * p0 ) )
      p.t.10 =  e0^xs[,1] * e0^xs[,2] * e0^xs[,3] * (1-e0)^(1-xs[,1]) * (1-e0)^(1-xs[,2]) * (1-e0)^(1-xs[,3]) * (1-delta) * (1-p0) / 
              ( d0^xs[,1] * d0^xs[,2] * d0^xs[,3] * (1-d0)^(1-xs[,1]) * (1-d0)^(1-xs[,2]) * (1-d0)^(1-xs[,3]) * ( delta * (1-p0) + (1-gamm) * p0 )  +
                e0^xs[,1] * e0^xs[,2] * e0^xs[,3] * (1-e0)^(1-xs[,1]) * (1-e0)^(1-xs[,2]) * (1-e0)^(1-xs[,3]) * ( (1-delta) * (1-p0) + gamm * p0 ) )
      p.t.11 =  e0^xs[,1] * e0^xs[,2] * e0^xs[,3] * (1-e0)^(1-xs[,1]) * (1-e0)^(1-xs[,2]) * (1-e0)^(1-xs[,3]) * gamm * p0 / 
              ( d0^xs[,1] * d0^xs[,2] * d0^xs[,3] * (1-d0)^(1-xs[,1]) * (1-d0)^(1-xs[,2]) * (1-d0)^(1-xs[,3]) * ( delta * (1-p0) + (1-gamm) * p0 )  +
                e0^xs[,1] * e0^xs[,2] * e0^xs[,3] * (1-e0)^(1-xs[,1]) * (1-e0)^(1-xs[,2]) * (1-e0)^(1-xs[,3]) * ( (1-delta) * (1-p0) + gamm * p0 ) )
      }
    
    if(n.cats == 16){
      p.t.01 =  d0^xs[,1] * d0^xs[,2] * d0^xs[,3] * d0^xs[,4] * (1-d0)^(1-xs[,1]) * (1-d0)^(1-xs[,2]) * (1-d0)^(1-xs[,3]) * (1-d0)^(1-xs[,4]) * (1-gamm) * p0 / 
              ( d0^xs[,1] * d0^xs[,2] * d0^xs[,3] * d0^xs[,4] * (1-d0)^(1-xs[,1]) * (1-d0)^(1-xs[,2]) * (1-d0)^(1-xs[,3]) * (1-d0)^(1-xs[,4]) * ( delta * (1-p0) + (1-gamm) * p0 )  +
                e0^xs[,1] * e0^xs[,2] * e0^xs[,3] * e0^xs[,4] * (1-e0)^(1-xs[,1]) * (1-e0)^(1-xs[,2]) * (1-e0)^(1-xs[,3]) * (1-e0)^(1-xs[,4]) * ( (1-delta) * (1-p0) + gamm * p0 ) )
      p.t.10 =  e0^xs[,1] * e0^xs[,2] * e0^xs[,3] * e0^xs[,4] * (1-e0)^(1-xs[,1]) * (1-e0)^(1-xs[,2]) * (1-e0)^(1-xs[,3]) * (1-e0)^(1-xs[,4]) * (1-delta) * (1-p0) / 
              ( d0^xs[,1] * d0^xs[,2] * d0^xs[,3] * d0^xs[,4] * (1-d0)^(1-xs[,1]) * (1-d0)^(1-xs[,2]) * (1-d0)^(1-xs[,3]) * (1-d0)^(1-xs[,4]) * ( delta * (1-p0) + (1-gamm) * p0 )  +
                e0^xs[,1] * e0^xs[,2] * e0^xs[,3] * e0^xs[,4] * (1-e0)^(1-xs[,1]) * (1-e0)^(1-xs[,2]) * (1-e0)^(1-xs[,3]) * (1-e0)^(1-xs[,4]) * ( (1-delta) * (1-p0) + gamm * p0 ) )
      p.t.11 =  e0^xs[,1] * e0^xs[,2] * e0^xs[,3] * e0^xs[,4] * (1-e0)^(1-xs[,1]) * (1-e0)^(1-xs[,2]) * (1-e0)^(1-xs[,3]) * (1-e0)^(1-xs[,4]) * gamm * p0 / 
              ( d0^xs[,1] * d0^xs[,2] * d0^xs[,3] * d0^xs[,4] * (1-d0)^(1-xs[,1]) * (1-d0)^(1-xs[,2]) * (1-d0)^(1-xs[,3]) * (1-d0)^(1-xs[,4]) * ( delta * (1-p0) + (1-gamm) * p0 )  +
                e0^xs[,1] * e0^xs[,2] * e0^xs[,3] * e0^xs[,4] * (1-e0)^(1-xs[,1]) * (1-e0)^(1-xs[,2]) * (1-e0)^(1-xs[,3]) * (1-e0)^(1-xs[,4]) * ( (1-delta) * (1-p0) + gamm * p0 ) )
      }
  }
  
  if(length(d0) == 2){
    p.t.01 =  d0[1]^xs[,1] * d0[2]^xs[,2] * (1-d0[1])^(1-xs[,1]) * (1-d0[2])^(1-xs[,2]) * (1-gamm) * p0 /
            ( d0[1]^xs[,1] * d0[2]^xs[,2] * (1-d0[1])^(1-xs[,1]) * (1-d0[2])^(1-xs[,2]) * ( delta * (1-p0) + (1-gamm) * p0 )  +
              e0[1]^xs[,1] * e0[2]^xs[,2] * (1-e0[1])^(1-xs[,1]) * (1-e0[2])^(1-xs[,2]) * ( (1-delta) * (1-p0) + gamm * p0 ) )
    p.t.10 =  e0[1]^xs[,1] * e0[2]^xs[,2] * (1-e0[1])^(1-xs[,1]) * (1-e0[2])^(1-xs[,2]) * (1-delta) * (1-p0) /
            ( d0[1]^xs[,1] * d0[2]^xs[,2] * (1-d0[1])^(1-xs[,1]) * (1-d0[2])^(1-xs[,2]) * ( delta * (1-p0) + (1-gamm) * p0 )  +
              e0[1]^xs[,1] * e0[2]^xs[,2] * (1-e0[1])^(1-xs[,1]) * (1-e0[2])^(1-xs[,2]) * ( (1-delta) * (1-p0) + gamm * p0 ) )
    p.t.11 =  e0[1]^xs[,1] * e0[2]^xs[,2] * (1-e0[1])^(1-xs[,1]) * (1-e0[2])^(1-xs[,2]) * gamm * p0 /
            ( d0[1]^xs[,1] * d0[2]^xs[,2] * (1-d0[1])^(1-xs[,1]) * (1-d0[2])^(1-xs[,2]) * ( delta * (1-p0) + (1-gamm) * p0 )  +
              e0[1]^xs[,1] * e0[2]^xs[,2] * (1-e0[1])^(1-xs[,1]) * (1-e0[2])^(1-xs[,2]) * ( (1-delta) * (1-p0) + gamm * p0 ) )
  }    
  
  if(length(d0) == 3){ 
    p.t.01 =  d0[1]^xs[,1] * d0[2]^xs[,2] * d0[3]^xs[,3] * (1-d0[1])^(1-xs[,1]) * (1-d0[2])^(1-xs[,2]) * (1-d0[3])^(1-xs[,3]) * (1-gamm) * p0 / 
            ( d0[1]^xs[,1] * d0[2]^xs[,2] * d0[3]^xs[,3] * (1-d0[1])^(1-xs[,1]) * (1-d0[2])^(1-xs[,2]) * (1-d0[3])^(1-xs[,3]) * ( delta * (1-p0) + (1-gamm) * p0 )  +
              e0[1]^xs[,1] * e0[2]^xs[,2] * e0[3]^xs[,3] * (1-e0[1])^(1-xs[,1]) * (1-e0[2])^(1-xs[,2]) * (1-e0[3])^(1-xs[,3]) * ( (1-delta) * (1-p0) + gamm * p0 ) )
    p.t.10 =  e0[1]^xs[,1] * e0[2]^xs[,2] * e0[3]^xs[,3] * (1-e0[1])^(1-xs[,1]) * (1-e0[2])^(1-xs[,2]) * (1-e0[3])^(1-xs[,3]) * (1-delta) * (1-p0) / 
            ( d0[1]^xs[,1] * d0[2]^xs[,2] * d0[3]^xs[,3] * (1-d0[1])^(1-xs[,1]) * (1-d0[2])^(1-xs[,2]) * (1-d0[3])^(1-xs[,3]) * ( delta * (1-p0) + (1-gamm) * p0 )  +
              e0[1]^xs[,1] * e0[2]^xs[,2] * e0[3]^xs[,3] * (1-e0[1])^(1-xs[,1]) * (1-e0[2])^(1-xs[,2]) * (1-e0[3])^(1-xs[,3]) * ( (1-delta) * (1-p0) + gamm * p0 ) )
    p.t.11 =  e0[1]^xs[,1] * e0[2]^xs[,2] * e0[3]^xs[,3] * (1-e0[1])^(1-xs[,1]) * (1-e0[2])^(1-xs[,2]) * (1-e0[3])^(1-xs[,3]) * gamm * p0 / 
            ( d0[1]^xs[,1] * d0[2]^xs[,2] * d0[3]^xs[,3] * (1-d0[1])^(1-xs[,1]) * (1-d0[2])^(1-xs[,2]) * (1-d0[3])^(1-xs[,3]) * ( delta * (1-p0) + (1-gamm) * p0 )  +
              e0[1]^xs[,1] * e0[2]^xs[,2] * e0[3]^xs[,3] * (1-e0[1])^(1-xs[,1]) * (1-e0[2])^(1-xs[,2]) * (1-e0[3])^(1-xs[,3]) * ( (1-delta) * (1-p0) + gamm * p0 ) )
  }
  
  if(length(d0) == 4){
    p.t.01 =  d0[1]^xs[,1] * d0[2]^xs[,2] * d0[3]^xs[,3] * d0[4]^xs[,4] * (1-d0[1])^(1-xs[,1]) * (1-d0[2])^(1-xs[,2]) * (1-d0[3])^(1-xs[,3]) * (1-d0[4])^(1-xs[,4]) * (1-gamm) * p0 / 
            ( d0[1]^xs[,1] * d0[2]^xs[,2] * d0[3]^xs[,3] * d0[4]^xs[,4] * (1-d0[1])^(1-xs[,1]) * (1-d0[2])^(1-xs[,2]) * (1-d0[3])^(1-xs[,3]) * (1-d0[4])^(1-xs[,4]) * ( delta * (1-p0) + (1-gamm) * p0 )  +
              e0[1]^xs[,1] * e0[2]^xs[,2] * e0[3]^xs[,3] * e0[4]^xs[,4] * (1-e0[1])^(1-xs[,1]) * (1-e0[2])^(1-xs[,2]) * (1-e0[3])^(1-xs[,3]) * (1-e0[4])^(1-xs[,4]) * ( (1-delta) * (1-p0) + gamm * p0 ) )
    
    p.t.10 =  e0[1]^xs[,1] * e0[2]^xs[,2] * e0[3]^xs[,3] * e0[4]^xs[,4] * (1-e0[1])^(1-xs[,1]) * (1-e0[2])^(1-xs[,2]) * (1-e0[3])^(1-xs[,3]) * (1-e0[4])^(1-xs[,4]) * (1-delta) * (1-p0) / 
            ( d0[1]^xs[,1] * d0[2]^xs[,2] * d0[3]^xs[,3] * d0[4]^xs[,4] * (1-d0[1])^(1-xs[,1]) * (1-d0[2])^(1-xs[,2]) * (1-d0[3])^(1-xs[,3]) * (1-d0[4])^(1-xs[,4]) * ( delta * (1-p0) + (1-gamm) * p0 )  +
              e0[1]^xs[,1] * e0[2]^xs[,2] * e0[3]^xs[,3] * e0[4]^xs[,4] * (1-e0[1])^(1-xs[,1]) * (1-e0[2])^(1-xs[,2]) * (1-e0[3])^(1-xs[,3]) * (1-e0[4])^(1-xs[,4]) * ( (1-delta) * (1-p0) + gamm * p0 ) )
    
    p.t.11 =  e0[1]^xs[,1] * e0[2]^xs[,2] * e0[3]^xs[,3] * e0[4]^xs[,4] * (1-e0[1])^(1-xs[,1]) * (1-e0[2])^(1-xs[,2]) * (1-e0[3])^(1-xs[,3]) * (1-e0[4])^(1-xs[,4]) * gamm * p0 / 
            ( d0[1]^xs[,1] * d0[2]^xs[,2] * d0[3]^xs[,3] * d0[4]^xs[,4] * (1-d0[1])^(1-xs[,1]) * (1-d0[2])^(1-xs[,2]) * (1-d0[3])^(1-xs[,3]) * (1-d0[4])^(1-xs[,4]) * ( delta * (1-p0) + (1-gamm) * p0 )  +
              e0[1]^xs[,1] * e0[2]^xs[,2] * e0[3]^xs[,3] * e0[4]^xs[,4] * (1-e0[1])^(1-xs[,1]) * (1-e0[2])^(1-xs[,2]) * (1-e0[3])^(1-xs[,3]) * (1-e0[4])^(1-xs[,4]) * ( (1-delta) * (1-p0) + gamm * p0 ) )
  }

  ####### Maximization Step
  # Estimate p_temp
  p_temp <- sum(groups*(p.t.01 + p.t.11)) / sum(groups)
  
  # Estimate delta_temp, gamm_temp
  delta_temp <- sum(groups*(1 - p.t.01 - p.t.10 - p.t.11)) / sum(groups*(1 - p.t.01 - p.t.11))
  gamm_temp  <- sum(groups*p.t.11) / sum(groups*(p.t.01 + p.t.11))
  
  # Estimate d0_temp, e0_temp
  if(de_1) { 
    d0_temp = e0_temp = array(0,dim = c(1,1))
    row.algo <- rowSums(pred.mat(n.algos))
    d0_temp <- sum(groups*(1-p.t.10 - p.t.11)*row.algo) / sum(groups*(1-p.t.10-p.t.11)*ncol(pred.mat(n.algos))) 
    e0_temp <- sum(groups*(p.t.10 + p.t.11)*row.algo) / sum(groups*(p.t.10+p.t.11)*ncol(pred.mat(n.algos)))
    
  } else { 
    d0_temp = e0_temp = array(0,dim = c(n.algos,1))
    for(i in 1:n.algos){
      ind <- which(pred.mat(n.algos)[,i] == 1)
      d0_temp[i] <- sum(groups[ind]*(1-p.t.10[ind]-p.t.11[ind])) / sum(groups*(1-p.t.10-p.t.11)) 
      e0_temp[i] <- sum(groups[ind]*(p.t.10[ind]+p.t.11[ind])) / sum(groups*(p.t.10+p.t.11))
    }
  }
  
  if(expect){ list(p.t.01, p.t.10, p.t.11) 
  } else { list(d0_temp, e0_temp, delta_temp, gamm_temp, p_temp) }

} 												




#####################################################################
# capture.em.xy
# groups = the observed data (i.e. which category each mutation fell into)
# input = initial parameters
# n.algos = Number algorithms to make predictions for each of the mutations; the number of categories is determined by the number of algorithms
# ab_1 = IF TRUE, then the EM algorithm is estimated assuming the a0, b0, p0 is the same for all the algorithm
#        IF FALSE, then the EM algorithm is estimated assuming each algorithm has its own a0, b0, p0\
# N.random.start = Number of times to randomly start the EM 
# maxit: max iterations  												
# tol: tolerance for determining convergence
#####################################################################

capture.em.xy <- function(groups, random.start = FALSE, input = NA, n.algos, ab_1 = FALSE, N.random.start = 100, maxit=1000,tol=1e-6){
  n.cats = 2^n.algos # create the number of categories regardless of the number of parameters to estimate (e.g. length(a0) )
  stopifnot(length(groups) == n.cats)
  
  if(random.start) { 
    flag <- 0
    lkd_iter = array(0, dim = N.random.start)
    
    if(ab_1){
      par_iter = array(0, dim = c(N.random.start,3))
    } else { par_iter = array(0, dim = c(N.random.start, 2*n.algos + 1)) }
    
    for(lkd in 1:N.random.start){
      if(ab_1){
        a_cur <- runif(1, min = 0, max = 0.5)
        b_cur <- runif(1, min = 0.5, max = 1)
        p_cur <- runif(1, min = 0, max = 1)
        } else {  a_cur <- runif(n.algos, min = 0, max = 0.5)
                  b_cur <- runif(n.algos, min = 0.5, max = 1)
                  p_cur <- runif(1, min = 0, max = 1) 
        }
    
      for(i in 1:maxit){
        cur <- c(a_cur,b_cur,p_cur)
        new <- capture.em.step.xy(groups,list(a_cur,b_cur,p_cur), n.algos, ab_1); a_new <- c(new[[1]]); b_new <- c(new[[2]]); p_new <- c(new[[3]])
        new_step <- c(a_new,b_new,p_new)
      
        if( all(abs(cur - new_step) < tol) ) { 
          flag <- 1; break
          } else { a_cur <- a_new; b_cur <- b_new; p_cur <- p_new }
      }
      
      par_iter[lkd,] = c(a_new, b_new, p_new)
      
      p.t = sim.categories.xy(n.size = 10, n.algos, list(a_new, b_new, p_new), lkd=TRUE)
      lkd_iter[lkd] = lgamma(sum(groups) + 1) - sum(lgamma(groups + 1)) + sum(groups*log(p.t) )
    }
    max.lkd_N = which(lkd_iter == max(lkd_iter), arr.ind = TRUE)
    
    if(ab_1){
      list(par_iter[max.lkd_N,1],par_iter[max.lkd_N,2],par_iter[max.lkd_N,3]) 
    } else {
      list(par_iter[max.lkd_N,1:n.algos],par_iter[max.lkd_N,(n.algos+1):(2*n.algos)],par_iter[max.lkd_N,(2*n.algos + 1)])
    }
      
  } else { 
    flag <- 0
    a_cur <- input[[1]]; b_cur <- input[[2]]; p_cur <- input[[3]]
    
    # Sanity check
    if(ab_1){ stopifnot(length(a_cur) == 1)}
    if(length(a_cur) != 1){ stopifnot(length(a_cur) == n.algos) }
    
    for(i in 1:maxit){
      cur <- c(a_cur,b_cur,p_cur)
      new <- capture.em.step.xy(groups,list(a_cur,b_cur,p_cur),n.algos, ab_1); a_new <- c(new[[1]]); b_new <- c(new[[2]]); p_new <- c(new[[3]])
      new_step <- c(a_new,b_new,p_new)
      
      if( all(abs(cur - new_step) < tol) ) { 
        flag <- 1; break
        } else { a_cur <- a_new; b_cur <- b_new; p_cur <- p_new }
    }
    
    list(a_cur,b_cur,p_cur)
  }
  
}




#####################################################################
# capture.em.xyz
# groups = the observed data (i.e. which category each mutation fell into)
# input = initial parameters
# n.algos = Number algorithms to make predictions for each of the mutations; the number of categories is determined by the number of algorithms
# ab_1 = IF TRUE, then the EM algorithm is estimated assuming the d0, e0 is the same for all the algorithm
#        IF FALSE, then the EM algorithm is estimated assuming each algorithm has its own d0, e0
# N.random.start = Number of times to randomly start the EM 
# maxit: max iterations    											
# tol: tolerance for determining convergence
#####################################################################

capture.em.xyz <- function(groups, random.start = FALSE, input = NA, n.algos, de_1 = FALSE, N.random.start = 100, maxit=1000,tol=1e-6){
  n.cats = 2^n.algos # create the number of categories regardless of the number of parameters to estimate (e.g. length(a0) )
  stopifnot(length(groups) == n.cats)
  
  if(random.start) { 
    flag <- 0
    lkd_iter = array(0, dim = N.random.start)
    
    if(de_1){
      par_iter = array(0, dim = c(N.random.start,5))
    } else { par_iter = array(0, dim = c(N.random.start, 2*n.algos + 3)) }
    
    for(lkd in 1:N.random.start){
      if(de_1){
        d0_cur <- runif(1, min = 0, max = 0.5)
        e0_cur <- runif(1, min = 0.5, max = 1)
        delta_cur <- runif(1, min = 0, max = 0.5)
        gamm_cur <- runif(1, min = 0, max = 0.5)
        p_cur <- runif(1, min = 0, max = 1)
      } else {  d0_cur <- runif(n.algos, min = 0, max = 0.5)
                e0_cur <- runif(n.algos, min = 0.5, max = 1)
                delta_cur <- runif(1, min = 0.8, max = 1)
                gamm_cur <- runif(1, min = 0.8, max = 1)
                p_cur <- runif(1, min = 0, max = 1)
    }    
      
      for(i in 1:maxit){
        cur <- c(d0_cur,e0_cur,delta_cur,gamm_cur,p_cur)
        new <- capture.em.step.xyz(groups,list(d0_cur,e0_cur,delta_cur,gamm_cur,p_cur), n.algos, de_1)
        d0_new <- c(new[[1]]); e0_new <- c(new[[2]]); delta_new <- c(new[[3]]); gamm_new <- c(new[[4]]); p_new <- c(new[[5]])
        new_step <- c(d0_new,e0_new,delta_new,gamm_new,p_new)
        
        if( all(abs(cur - new_step) < tol) ) {
          flag <- 1; break 
        } else { d0_cur <- d0_new; e0_cur <- e0_new; delta_cur <- delta_new; gamm_cur <- gamm_new; p_cur <- p_new }
      }
      
      par_iter[lkd,] = c(d0_new, e0_new, delta_new, gamm_new, p_new)
      
      p.t = sim.categories.xyz(n.size = 10, n.algos, list( d0_new, e0_new, delta_new, gamm_new, p_new ), lkd=TRUE)
      lkd_iter[lkd] = lgamma(sum(groups) + 1) - sum(lgamma(groups + 1)) + sum(groups*log(p.t) )
    }
    max.lkd_N = which(lkd_iter == max(lkd_iter), arr.ind = TRUE)

    if(de_1){
      list(par_iter[max.lkd_N,1],par_iter[max.lkd_N,2],par_iter[max.lkd_N,3], par_iter[max.lkd_N,4],par_iter[max.lkd_N,5] ) 
    } else {
      list(par_iter[max.lkd_N,1:n.algos],par_iter[max.lkd_N,(n.algos+1):(2*n.algos)],par_iter[max.lkd_N,(2*n.algos + 1)],par_iter[max.lkd_N,(2*n.algos + 2)],par_iter[max.lkd_N,(2*n.algos + 3)])
    }
    
  } else { 

    d0_cur <- input[[1]]; e0_cur <- input[[2]]; delta_cur <- input[[3]]; gamm_cur <- input[[4]]; p_cur <- input[[5]]
    flag <- 0
    
    # Sanity check
    if(de_1){ stopifnot(length(d0_cur) == 1)}
    if(length(d0_cur) != 1){ stopifnot(length(d0_cur) == n.algos) }

    for(i in 1:maxit){
      cur <- c(d0_cur,e0_cur,delta_cur,gamm_cur,p_cur)
      new <- capture.em.step.xyz(groups,list(d0_cur,e0_cur,delta_cur,gamm_cur,p_cur), n.algos, de_1)
      d0_new <- c(new[[1]]); e0_new <- c(new[[2]]); delta_new <- c(new[[3]]); gamm_new <- c(new[[4]]); p_new <- c(new[[5]])
      new_step <- c(d0_new,e0_new,delta_new,gamm_new,p_new)
      
      if( all(abs(cur - new_step) < tol) ) {
        flag <- 1; break 
      } else { d0_cur <- d0_new; e0_cur <- e0_new; delta_cur <- delta_new; gamm_cur <- gamm_new; p_cur <- p_new }
    }
  list(d0_cur,e0_cur,delta_cur,gamm_cur,p_cur)  	
  }
}   




# Computes dj and ej given aj, bj, delta and gamma
d_e <- function(a0,b0,delta,gamm){
  d0 = e0 = array(0,dim=c(length(a0)))
  for(i in 1:length(a0)){
    d0[i] = (b0[i] - b0[i]*delta - a0[i]*gamm) / (1- delta - gamm)
    e0[i] = (a0[i] - a0[i]*gamm - b0[i]*delta) / (1- delta - gamm)
  }
  list(d0,e0)
}
# d_e(a.true,b.true,delta.true,gamm.true)

# Computes aj and bj given dj, ej, delta and gamma
a_b <- function(em.output){
  em.output <- unlist(em.output)
  if(length(em.output) == 7){
    d0 <- em.output[1:2]; e0 <- em.output[3:4]; delta <- em.output[5]; gamm <- em.output[6]; p0 <- em.output[7]
    
    a0 = b0 = array(0,dim=c(2,1))
    a0[1] = d0[1]*delta + e0[1]*(1-delta)
    a0[2] = d0[2]*delta + e0[2]*(1-delta)
    
    b0[1] = d0[1]*(1-gamm) + e0[1]*gamm
    b0[2] = d0[2]*(1-gamm) + e0[2]*gamm
  }
  
  if(length(em.output) == 9){
    d0 <- em.output[1:3]; e0 <- em.output[4:6]; delta <- em.output[7]; gamm <- em.output[8]; p0 <- em.output[9]
    
    a0 = b0 = array(0,dim=c(3,1))
    a0[1] = d0[1]*delta + e0[1]*(1-delta)
    a0[2] = d0[2]*delta + e0[2]*(1-delta)
    a0[3] = d0[3]*delta + e0[3]*(1-delta)
    
    b0[1] = d0[1]*(1-gamm) + e0[1]*gamm
    b0[2] = d0[2]*(1-gamm) + e0[2]*gamm
    b0[3] = d0[3]*(1-gamm) + e0[3]*gamm
  }
  
  if(length(em.output) == 11){
      d0 <- em.output[1:4]; e0 <- em.output[5:8]; delta <- em.output[9]; gamm <- em.output[10]; p0 <- em.output[11]
  
      a0 = b0 = array(0,dim=c(4,1))
      a0[1] = d0[1]*delta + e0[1]*(1-delta)
      a0[2] = d0[2]*delta + e0[2]*(1-delta)
      a0[3] = d0[3]*delta + e0[3]*(1-delta)
      a0[4] = d0[4]*delta + e0[4]*(1-delta)
  
      b0[1] = d0[1]*(1-gamm) + e0[1]*gamm
      b0[2] = d0[2]*(1-gamm) + e0[2]*gamm
      b0[3] = d0[3]*(1-gamm) + e0[3]*gamm
      b0[4] = d0[4]*(1-gamm) + e0[4]*gamm
  }
  
  list(a0,b0)  
}
# a_b(true.xyz)



