## Author : Ingrid Vilmus
##

library(reshape2)

##########################################################################################################################################################
# function to simulate genotype frequencies

# nhab corresponds to the number of habitats
# shared is a logical variable, TRUE if the genotypes can belong to several habitat and FALSE if genotypes are specific to habitat
# if shared=FALSE, the user must furnish:
# * ngeno_hab which can be an integer representing the same number of genotype per habitat or it can be a vector of integer with length = nhab for which
# each element represent the number of genotype in each habitat
# * balanced which is a logical value, set to TRUE if genotype frequencies are balanced per habitat and to FALSE for non balanced frequencies
# in this case, the user can furnish alpha wich is a list containing nhab element and for each is a weight vector with length = ngeno_hab for the habitat concerned
# example : ngeno_hab = 5, nhab = 2 means 5 genotypes for each habitat. alpha can be for example list(c(1,2,1,3,5),c(4,8,7,1,2))
# without alpha, the weights are sampled in a uniform distribution from 1 to n (default = 10)
# if shared = TRUE, the user can modify:
# * m parameter which is used to sample genotypes frequencies per habitat using a negative binomial distribution
# default for m is set to 1 meaning for each habitat, most of the genotypes will have null or very small frequencies and a few genotypes will have high frequencies
# large values for m means genotype frequencies for each habitat will follow a Normal distribution
# * ngeno_hab which corresponds to the global number of genotypes shared accross habitats and must be an integer

sim_genofreq<-function(nhab,shared=TRUE,ngeno_hab,m=1,balanced=FALSE,alpha=NULL,n=10) {
  if (shared == TRUE) {
    freq_geno<-matrix(0,nhab,ngeno_hab)
    for (i in 1:nhab) {
      tmp<-rnbinom(ngeno_hab,m,0.5)
      freq_geno[i,]<-tmp/sum(tmp)
    }
    
  } else {
    if (length(ngeno_hab)==1) {
      ngeno_hab<-rep(ngeno_hab,nhab)
    } 
    
    freq_geno<-matrix(0,nhab,0)
    if (is.null(alpha)) {
      alpha<-list()
      for (i in 1:nhab) {
        alpha[[i]]<-sample(x=1:n,size=ngeno_hab[i])
      }
    }
    if (balanced==FALSE) {

      for (i in 1:nhab) {
        tmp<-matrix(0,nhab,ngeno_hab[i])
        tmp[i,]<-alpha[[i]]/sum(alpha[[i]])
        freq_geno<-cbind(freq_geno,tmp)
      }
      
    } else {
      for (i in 1:nhab) {
        tmp<-matrix(0,nhab,ngeno_hab[i])
        tmp[i,]<-1/ngeno_hab[i]
        freq_geno<-cbind(freq_geno,tmp)
      }
    }
  }
  return(freq_geno)
}



#########################################################################################################################################################
# function to simulate data at the genotype scale without accounting for deme entity to simulate survival probabilities with or without local adaptation
# params
# size correspond to the size of each observation, default value is set to 1 and correspond to a Bernouilli trial
# can be an integer if the size is the same or a vector if the size depend on the observation
# ngeno chatracterize the number of genotypes in the metapopulation
# nhab characterize the number of habitat ; one habitat correspond to location from where a subset of the metapopulation come from
# m is the parameter used to sample genotype frequencies per habitat with a negative binomial distribution
# m represent the number of success to obtain
# if m is high, the law converge to a Poisson
# if m is small (m<1), the data are overdispersed, many genotypes having very low frequencies and a few genotypes having high frequencies
# if m = 1, the law becomes geometric, the variance increasing with the mean value
# nrep is the number of replication per habitat

# values of standard deviation must be defined in a list called sdlist with those names:
# * sdhab is the standard deviation at the habitat-level
# * sdrep is the standard deviation at the replicat:habitat interaction level
# * sdint is the standard deviation around the fixed intercept (correspond to a genotype variation, pure geno effect)
# * sdslop is the standard deviation around the fixed slope (correspond to a genotype variation, response to the local adaptation)
# * sdgenohb is the standard deviation at the geno:hab interaction level
# * sdres is the standard deviation at the intra replicat level
# values for average slope and intercept characterizing local adaptation at the metapopulation scale
# a correspond to the fixed average slope (across all the observations)
# b correspond to the fixed average intercept (across all the observations)
# to simulate with or witout local adaptation at the metapopulation scale, adjust a and b values using the shiny app dedicated leaving others values to default
# default values indicates high local adaptation

# linkfamily indicates the link function used to transform probabilities at a continous scale
# freq correspond to the genotype frequencies in each habitat, must be a matrix (nhab,ngeno), one can use the function sim_genofreq to parameterize the sampling
# by default, the function is used with shared parameter



LA_sim<-function(size=1,ngeno,nhab,nrep,sdlist=list(sdhab=2,sdrep=2,sdint=2,sdslop=2,sdres=0),ris=0,a=10,b=-5,linkfamily="logit",freq=NULL,m=1) {

  #### design
  mat_design<-expand.grid(geno=seq(1,ngeno),replicat=seq(1,nrep),hab=seq(1,nhab))
  for (i in 1:ncol(mat_design)) {mat_design[,i]<-as.factor(mat_design[,i])}
  
  #### frequencies
  if (is.null(freq)) {
    freq<-sim_genofreq(nhab=nhab,shared=TRUE,ngeno_hab=ngeno)
  }
  
  X<-melt(freq)
  colnames(X)<-c("hab","geno","frequency")
  data_sim<-merge(mat_design,X,by=c("geno","hab"),sort=FALSE)
  data_sim<-data_sim[order(data_sim$geno,data_sim$hab,data_sim$replicat),]
  # add the size of each observation
  data_sim$size<-size
  
  #### effects
  # sample the habitat effects
  hab_eff<-rnorm(nhab,0,sdlist$sdhab)
  data_sim$hab_eff<-model.matrix(~-1+hab,data_sim)%*%hab_eff
  # sample the replicat effect
  rep_eff<-rnorm(nhab*nrep,0,sdlist$sdrep)
  data_sim$rep_eff<-model.matrix(~-1+replicat:hab,data_sim)%*%rep_eff
  # sample the random intercepts and the random slopes according to the value of the correlation between both vectors effects
  Sigma<-matrix(c(sdlist$sdint^2,ris*sdlist$sdint*sdlist$sdslop,ris*sdlist$sdint*sdlist$sdslop,sdlist$sdslop^2),2,2)
  geno_eff<-mvrnorm(ngeno,c(0,0),Sigma)

  data_sim$int_effect<-model.matrix(~-1+geno,data_sim)%*%geno_eff[,1]
  data_sim$slop_effect<-model.matrix(~-1+geno,data_sim)%*%geno_eff[,2]
  # sample the residual effect at the observation scale (correspond to an intra replicat residual variation)
  resid_eff<-rnorm(nrow(data_sim),0,sdlist$sdres)
  data_sim$resid_eff<-resid_eff
  # add fixed intercept and slope for all the observations
  data_sim$slope<-a
  data_sim$intercept<-b
  
  # calculate the response variable at the latent scale
  data_sim$Y<-data_sim$frequency*(data_sim$slope+data_sim$slop_effect)+data_sim$intercept+data_sim$int_effect+data_sim$hab_eff+data_sim$rep_eff+data_sim$resid_eff
  if (linkfamily=="logit") {
    data_sim$P<-plogis(data_sim$Y)
  } else if (linkfamily=="probit") {
    data_sim$P<-pnorm(data_sim$Y)
  }
  
  # account for overdispersion at the observation level
  data_sim$survivors<-rbinom(n=1:nrow(data_sim),size=data_sim$size,prob=data_sim$P)
  data_sim$obs<-seq(1,nrow(data_sim))
  
}

