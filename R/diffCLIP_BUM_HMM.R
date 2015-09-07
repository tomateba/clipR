#' @title Function BUM_FWBW
#' @param filename a matrix containing the number of mapped reads in each location, for each sample.
#' @param Nc number of replicates in condition 1 (or control)
#' @param Nt number of replicates in condition 2 (or treatment)
#' @return a matrix containing posteriors in Beta-Uniform mixture HMM
#' @description diffCLIP_BUM_HMM compares CLIP data using HMM with Beta-Uniform Mixture emission (BUM model). Takes as input a matrix, two integers (number of replicates in condition 1 and condition 2 respectively), and generates a vector containing all the log ratios in read counts between replicates at every location.
#' @export

diffCLIP_BUM_HMM<-function(filename,Nc,Nt){
  dataNum<-filename
  totalN=dim(dataNum)[2]
  numNucleotides=dim(dataNum)[1]
  if(totalN!=Nc+Nt){stop('inconsistent data table and number of replicates')}
  dataNum[which(dataNum==0,arr.ind=T)]=1 #gets rid of read counts equal to zero, setting them to 1
  
  dataControl=dataNum[,1:Nc]
  dataTreatment=dataNum[,(Nc+1):(Nc+Nt)]
  numComp1= (Nc*(Nc-1))/2 #number of pairwise comparisons between controls
  numComp2= (Nt*(Nt-1))/2 #number of pairwise comparisons between treatments
  
  # For stratification, I don't execute the original code, logRatioRepsCont is the input
  
  logRatioRepsCont1=matrix(0,nrow=numNucleotides,ncol=numComp1)#big matrix containing all the log ratios in read counts between controls and treatments at every location
  
  logRatioRep=matrix(0,nrow=numNucleotides,ncol=1)
  for(j in 1:(Nc-1)){
    for(n in (j+1):Nc){
      logRatioRep[,1]=log(dataControl[,j])-log(dataControl[,n])
      index=(j-1)*Nc+n-j*(j+1)/2 #transforms indices in the more compact representation
      logRatioRepsCont1[,index]=logRatioRep
    }
    
  }	
  
  logRatioRepsCont1<-as.vector(logRatioRepsCont1)
  
  logRatioRepsCont2=matrix(0,nrow=numNucleotides,ncol=numComp2)#big matrix containing all the log ratios in read counts between controls and treatments at every location
  
  logRatioRep2=matrix(0,nrow=numNucleotides,ncol=1)
  for(j in 1:(Nt-1)){
    for(n in (j+1):Nt){
      logRatioRep2[,1]=log(dataTreatment[,j])-log(dataTreatment[,n])
      index=(j-1)*Nc+n-j*(j+1)/2 #transforms indices in the more compact representation
      logRatioRepsCont2[,index]=logRatioRep2
    }
    
  }	
  
  logRatioRepsCont2<-as.vector(logRatioRepsCont2)
  #dim(logRatioRepsCont)=c(1,numNucleotides*Nc*(Nc-1)/2)
  logRatioRepsCont<-c(logRatioRepsCont1,logRatioRepsCont2)
  
  logRatioTC=matrix(0,nrow=numNucleotides,ncol=Nc*Nt)
  
  #another big matrix containing the log ratios in read counts between each treatment sample and each control sample
  
  logRatio=matrix(0,nrow=numNucleotides,ncol=1)
  for(j in 1:Nc){
    for(n in 1:Nt){
      logRatio[,1]=log(dataTreatment[,n])-log(dataControl[,j])
      index=Nt*(j-1)+n #transforms indices in the more compact representation
      logRatioTC[,index]=logRatio
    }
    
  }
  
  
  empPvals=matrix(0,nrow=numNucleotides,ncol=Nc*Nt)
  
  #this contains the empirical p-values that a read count ratio as high as the observed
  for(i in 1:numNucleotides){#computes empirical p-values (one sided, only interested in larger read counts than normal)
    for(j in 1:(Nc*Nt)){
      empPvals[i,j]=sum(abs(logRatioRepsCont)>=abs(logRatioTC[i,j]))/length(logRatioRepsCont)
    }
  }
  empPvals=t(empPvals)
  
  #RUN THE BUM-HMM INFERENCE

  posterior=BUM_FWBW(empPvals)
  return(posterior)	
}
