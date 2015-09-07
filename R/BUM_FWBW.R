#' @title Function BUM_FWBW
#' @param pVals a matrix of empirical p-values at each position, 1 row per each experiment.
#' @param trans the transition matrix (2x2, we assume only 2 hidden states, state 1 is always the unbound state). By default, we expect unbound stretches of average length 20 and bound stretches of average length 5.
#' @param beta the beta parameter (the Beta component is assumed to be Beta(1,beta)). The default value, 10, is chosen heuristically so that a p-value of 0.2 has roughly equal probability of coming from each mixture.
#' @param in_prob vector with initial probabilities. By default, we assume we always start in unbound state, so the vector has value c(1,0).
#' @return a matrix containing posteriors in Beta-Uniform mixture HMM (possibly from multiple exps).
#' @export

#' @description Computes posteriors in Beta-Uniform mixture HMM (possibly from multiple exps). pVals is a matrix of empirical p-values at each position, 1 row per each experiment (only treatment).  trans is the transition matrix (2x2, we assume only 2 hidden states,state 1 is always the unbound state).beta is the beta parameter (the Beta component is assumed to be Beta(1,beta))

BUM_FWBW<-function(pVals,trans=matrix(c(0.95, 0.2, 0.05, 0.8),nrow=2,ncol=2),beta=10,in_prob=c(1,0)){
  nexp=length(pVals[,1])
  nBins=length(pVals[1,])
  nStates=length(trans[,1])
  beta=c(1,beta) #create vector of beta params; the first entry =1 implies the control distribution is uniform
  
  #Variables in FWBW
  obsLike=matrix(1,ncol=nBins,nrow=nStates)#log likelihood of observations given state
  fwdMessage=matrix(0,ncol=nBins,nrow=nStates)
  bwdMessage=matrix(0,ncol=nBins,nrow=nStates)
  
  #Calculation of likelihoods (sum over replicates of likelihoods of each experiment)
  for (index in 1:nexp){
    for (index2 in 1:nBins){
			obsLike[,index2]=obsLike[,index2]*dbeta(pVals[index,index2],1,beta)
		}
		# obsLike <- obsLike * exp( (alpha * log(beta/(1+beta)) - lgamma(alpha))%*%matrix(1, ncol=nBins)
				# -  log(1+beta)%*%pVals[index,] + lgamma( matrix(1, nrow=nStates)%*%pVals[index,] +   alpha%*%matrix(1, ncol=nBins))
				# -  matrix(1, nrow=nStates)%*%densNorm[index,] )
	}
	#Calculation of the forward messages
	fwdMessage[,1]=in_prob*obsLike[,1]
	fwdMessage[,1]=fwdMessage[,1]/sum(fwdMessage[,1])
	for(index in 2:nBins){
		fwdMessage[,index]=(trans%*%fwdMessage[,index-1])*obsLike[,index]
		fwdMessage[,index]=fwdMessage[,index]/sum(fwdMessage[,index])
	}
	#Calculation of the backward message
	bwdMessage[,nBins]=1
	for(index in (nBins-1):1){
		bwdMessage[,index]=(trans%*%bwdMessage[,index+1])*obsLike[,index]
				bwdMessage[,index]=bwdMessage[,index]/sum(bwdMessage[,index])

	}
	#Calculation of posteriors
	posterior=fwdMessage*bwdMessage
	posterior=posterior/(matrix(1,nrow=length(posterior[,1]))%*%colSums(posterior))
	return(posterior)
}