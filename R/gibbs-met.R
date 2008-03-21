## a generic function for Gibbs sampling with Metropolis steps
## iters_mc         --- iterations of Gibbs sampling
## iters_met        --- iterations of Metropolis for each 1-dimensional sampling
## log_f            --- the log of the density function
## p                --- the dim of variables sampled
## x0               --- the initial value
## stepsizes_met    --- a vector of length p, with
##                  stepsizes_met[i] being the standard deviation of Gaussian proposal
##                  for updating 'i'th parameter
## iters_per.iter   --- for each transition specified by `iters', 
##                  the Markov chain sampling is run 'iters_per.iter' times
##                  this argument is used to avoid saving the whole Markov chain
## ...          --- extra arguments needed to compute log_f

## a matrix with of (iters + 1) * p is returned, with each row for an iteration

gibbs_met <- function(log_f,p,x0,iters_mc,iters_met,stepsizes_met,
                      iters_per.iter=1,...)
{
   if(p != length(x0))
      stop("The number of variables in initial values does NOT match p")
   chain <- matrix(0, iters_mc + 1, p)
   chain[1,] <- x0

   ## update only 'i_par' th parameter in 'iter' iteration
   gibbs_update_per.iter_per.par <- function(iter, i_par)
   {   log_f_condition <- function(x_i)
       {   x <- chain[iter,]
           x[i_par] <- x_i
           log_f(x,...)
       }
       ## updating parameer with index i_par with Metropolis method
       chain[iter,i_par] <<- 
               met_gaussian(log_f_condition,iters_met,1,
                       chain[iter,i_par],stepsizes_met[i_par])[iters_met+1]
   }
   ## gibbs sampling for iteration `iter' for all parameters
   gibbs_update_per.iter <- function(iter) {
        ## copy the states in 'iter-1' iteration to this iteration
        chain[iter,] <<- chain[iter-1,]
        replicate(iters_per.iter,
                  sapply(1:p,gibbs_update_per.iter_per.par,iter=iter) )
   }
   ## perform iters_mc gibbs sampling for all parameters
   sapply(seq(2,iters_mc+1), gibbs_update_per.iter)
   chain
}

## a generic function for Metropolis sampling with Gaussian proposal
## iters            --- length of Markov chain
## log_f            --- the log of the density function
## p                --- the number of variables being sampled
## x0               --- the initial value
## stepsizes        --- the standard deviations of Gaussian proposal
##                  stepsizes[i] for the `i'th  variable 
## iters_per.iter   --- for each transition specified by `iters', 
##                  the Markov chain sampling is run 'iters_per.iter' times
##                  this argument is used to avoid saving the whole Markov chain
## ...              --- other arguments needed to compute log_f

## a matrix with of (iters + 1) * p is returned, with each row for an iteration
 
met_gaussian <- function(log_f,iters, p, x0, stepsizes, iters_per.iter=1, ...)
{   
    if(p != length(x0))
      stop("The number of variables in initial values does NOT match p")

    ## creating a matrix to save the Markov chain, with each row for an iteration
    chain <- matrix(0,iters+1,p)
    chain[1,] <- x0
    old_log_f <- log_f(x0,...)

    ## doing one transition
    one_transition <- function(i)
    {   chain[i,] <<- chain[i-1,]
        i_inside <- 0
        repeat{
                i_inside <- i_inside + 1
                ## propose a point
                x_prop <- rnorm(p) * stepsizes + chain[i,]
                ## decide whether to accept it
                new_log_f <- log_f(x_prop,...)
                if(log(runif(1)) < new_log_f - old_log_f) 
                {  chain[i,] <<- x_prop
                   old_log_f <<- new_log_f
                }
                if(i_inside == iters_per.iter) break
        }
    }
    ## perform iters Metropolis updating
    sapply(seq(2,iters+1),one_transition)
    ## return the whole chain 
    chain
}
