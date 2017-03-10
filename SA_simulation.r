# adapted from Blanquart et al 2013
#################################################################################
################################
# FUNCTION TO GENERATE DATA THAT FOLLOWS PERFECTLY THE ASSUMPTIONS OF THE LINEAR MODEL #
#################################################################################
library(gtools)


################################
# to generate fitness response (number of descendant)
make_data <- function(np, d, ind,geno, sd_ind, sd_pop, sd_hab, sd_geno, SA, constant,useDir=FALSE,alpha)
{
  #PARAMETER DEFINITIONS:
  # np is the number of populations tested
  # d indicates design, from d = 2 (twice as many allopatric as sympatric transplant), to np - 1 (full factorial design)
  # ind is the number of individuals tested by transplant (correspond to number of replicat per genotype)
  # geno is the number of genotypes in a deme
  # sd_ind is the individual - level standard deviation of error (GxE)
  # sd_pop is the population - level standard deviation of error (overdispersion sd term)
  # sd_hab is the standard deviation of habitat effects
  # sd_geno is the standard deviation of geno quality effects --> used to calculate demeffects
 
  # SA is the sympatric vs. allopatric contrast (i.e., local adaptation)
  # constant is the global mean
  # alpha is the vector of alpha values used to sample the genotype frequencies per deme using a Dirichlet distribution if useDir=TRUE
  # and it correspond to the relative weight of genotypes for frequency calculation
  # else the genotype frequencies are not allowed to vary over simulations and corresponds to alpha/sum(alpha)
  # if (length(alpha)=geno) the sum of the genotype frequencies per habitat are equal to 1 and if (length(alpha)>geno), virtual genotypes 
  # are created for the calcul of the frequencies (and thus the sum is < 1) but not used for the rest of the data making
  
  Xik<-matrix(0,np,geno*np)
  for (i in 1:np) {
    if (useDir==TRUE) {
      
      Xik[i,(i*geno-geno+1):(i*geno)]<-rdirichlet(1, alpha )[1:geno]
      
    } else {
      Xik[i,(i*geno-geno+1):(i*geno)]<-(alpha/sum(alpha))[1:geno]
    }
  }
  
  #habitat and deme quality effects are drawn in a normal distribution with sd as defined above:
  habeffects<-rnorm(np,m = 0,sd = sd_hab)
  genoeffects<-rnorm(np*geno,m=0, sd = sd_geno)
  demeffects<-as.numeric(Xik%*%genoeffects)
  #now generate the transplant matrix data:
  fitness <- c(); habs <- c(); dems <- c(); symp <- c();
  for(i in 1:np){ #loop on deme of destination
    for(j in 1:np){ #loop on deme of origin
      #draw the population-level error and add the "SA" contrast tosympatric transplants:
      if(i == j){error_pop <- rnorm(1,m = SA,sd = sd_pop)}
      if(i != j){error_pop <- rnorm(1,m = 0,sd = sd_pop)}
      for(k in 1:ind){
        #draw error for this individual:
        error_ind = rnorm(1,m = 0,sd = sd_ind);
        #Generate fitness data for sympatric transplant:
        if(i == j){
          fitness = append(fitness,constant + habeffects[i] +
                             demeffects[j] + error_pop + error_ind)
          symp = append(symp,1)
          habs=append(habs,i)
          dems=append(dems,j)
        }
        #Generate fitness data for allopatric transplant,depending on the experimental design d
        if((i-j >= 1 - np & i-j <= d - np) || (i-j >= 1 & i-j <=
                                               d)){
          fitness=append(fitness,constant + habeffects[i] +
                           demeffects[j] + error_pop + error_ind)
          symp = append(symp,0)
          habs=append(habs,i)
          dems=append(dems,j)
        }
      }
    }
  }
  #return the transplant experiment data
  return(data.frame(habs, dems, symp, fitness))
}


##########################################################################################################################

# to generate data at a latent scale (logit transformation) --> integrate geno and quadrat in the loop k for case where d!=np-1 --> TO DO
latent_data1 <- function(np, d, ind,geno, sd_ind, sd_pop, habeffects, sd_geno, SA, constant,useDir=FALSE,alpha)
{
  #PARAMETER DEFINITIONS:
  # np is the number of populations tested
  # d indicates design, from d = 2 (twice as many allopatric as sympatric transplant), to np - 1 (full factorial design)
  # ind is the number of quadrat (number of pop replicat per habitat * nb geno)
  
  # geno is the number of genotypes in a deme --> number total of geno = np * geno
  # sd_pop is the quadrat - level standard deviation of error 
  # sd_ind is the observation-level standard deviation of error
  # sd_geno is the standard deviation of geno quality effects --> used to calculate demeffects
  
  # constant is the global mean
  # habeffects is a vector of fixed coef for each zone
  
  # NB : the values must be set at the latent scale --> log(p/(1-p))
  # SA is the sympatric vs. allopatric contrast (i.e., local adaptation)
  
  # alpha is the vector of alpha values used to sample the genotype frequencies per deme using a Dirichlet distribution if useDir=TRUE
  # and it correspond to the relative weight of genotypes for frequency calculation
  # else the genotype frequencies are not allowed to vary over simulations and corresponds to alpha/sum(alpha)
  # if (length(alpha)=geno) the sum of the genotype frequencies per habitat are equal to 1 and if (length(alpha)>geno), virtual genotypes 
  # are created for the calcul of the frequencies (and thus the sum is < 1) but not used for the rest of the data making
  
  Xik<-matrix(0,np,geno*np)
  for (i in 1:np) {
    if (useDir==TRUE) {
      
      Xik[i,(i*geno-geno+1):(i*geno)]<-rdirichlet(1, alpha )[1:geno]
      
    } else {
      Xik[i,(i*geno-geno+1):(i*geno)]<-(alpha/sum(alpha))[1:geno]
    }
  }
  # vector of genotypes
  vec_geno<-paste0("g",1:(geno*np))
  colnames(Xik)<-vec_geno
  rownames(Xik)<-paste0("hab",1:np)
  # deme quality effects are drawn in a normal distribution with sd as defined above:
  
  genoeffects<-rnorm(np*geno,m=0, sd = sd_geno)
  demeffects<-as.numeric(Xik%*%genoeffects)
  #now generate the transplant matrix data:
  fitness <- c(); habs <- c(); dems <- c(); symp <- c(); quads<-c(); genot<-c()
  for(i in 1:np){ #loop on deme of destination
    for(j in 1:np){ #loop on deme of origin
      #draw the population-level error and add the "SA" contrast tosympatric transplants:
      if(i == j){error_pop <- rnorm(1,m = SA,sd = sd_pop)}
      if(i != j){error_pop <- rnorm(1,m = 0,sd = sd_pop)}
      
      # name of the genotype considered for this deme
      genot<-c(genot,rep(vec_geno[(j*geno-geno+1):(j*geno)],ind/np))
      # quadrat for this zone
      tmp1<-i*ind/np-ind/np+1
      tmp2<-i*ind/np
      quads<-c(quads,rep(paste0("q",seq(tmp1,tmp2)),np))
      
      # loop on individuals (nb of geno * nb of quadrat)
      for(k in 1:ind){
        #draw error for this individual:
        error_ind = rnorm(1,m = 0,sd = sd_ind);
        #Generate fitness data for sympatric transplant:
        if(i == j){
          fitness = append(fitness,constant + habeffects[i] +
                             demeffects[j] + error_pop + error_ind)
          symp = append(symp,1)
          habs=append(habs,i)
          dems=append(dems,j)
          
        }
        #Generate fitness data for allopatric transplant,depending on the experimental design d
        if((i-j >= 1 - np & i-j <= d - np) || (i-j >= 1 & i-j <=
                                               d)){
          fitness=append(fitness,constant + habeffects[i] +
                           demeffects[j] + error_pop + error_ind)
          symp = append(symp,0)
          habs=append(habs,i)
          dems=append(dems,j)
         
        }
      }
    }
  }
  #return the transplant experiment data
  W<-data.frame(habs, dems, symp, genot, quads, fitness)
  
  return(list(W=W,freq=Xik))
}



#########################################################################################################################################################################

# to generate data at a latent scale (logit transformation), only work for a full reciprocal transplant (d=np-1)

#PARAMETER DEFINITIONS:

# replicat is the number of replications per habitat for the whole population
# geno is a vector containing the number of genotype in each deme and thus the length of geno indicates the number of deme (and consequently habitat)

# sd_pop is the population - level standard deviation of error ()
# sd_error is the observation-level standard deviation of error
# sd_geno is the standard deviation of geno quality effects --> used to calculate demeffects

# mu is the global mean

# habeffects is an otpionnal vector of fixed coef for each habitat
# sd_habs is an optionnal habitat-level standard deviation if one want to set habitat as a random effect (if nhab > 5 for example)
# habeffects or sd_habs must be given

# NB : the values must be set at the latent scale --> log(p/(1-p))
# SA is the sympatric vs. allopatric contrast (i.e., local adaptation)

# alpha is a list containing as many elements as length(geno) to sample the genotype frequencies per deme 
# using a Dirichlet distribution if useDir=TRUE
# and it correspond to the relative weight of genotypes for frequency calculation
# else the genotype frequencies are not allowed to vary over simulations and corresponds to alpha/sum(alpha)
# if each element of alpha correspond to the number of geno for a given, the sum of the genotype frequencies per habitat are equal to 1 and 
# else, virtual genotypes are created for the calcul of the frequencies (and thus the sum is < 1) but not used for the rest of the data making

latent_data<-function(replicat,geno, sd_error, sd_pop, habeffects=NULL, sd_habs=NULL, sd_geno, SA, mu,useDir=FALSE,alpha) {
  # test
  if (is.null(habeffects)&is.null(sd_habs)) {
    stop("habeffects or sd_habs must be given")
  }
  
  
  # number of demes
  nd<-length(geno)
  # number of geno
  ngeno<-sum(geno)
  
  # output
  res<-matrix(NA,0,6)
  colnames(res)<-c("habitat","deme","genotype","replicat","symp","Y")
  res<-as.data.frame(res)
  
  # geno frequencies
  Xik<-matrix(0,nd,0)
  geno_names<-list()
  m<-0
  # loop on deme
  for (i in 1:nd) {
    # loop on geno within deme
    ni<-geno[i]
    Xik_tmp<-matrix(0,nd,ni)
    for (j in 1:ni) {
      if (useDir==TRUE) {
        Xik_tmp[i,1:ni]<-rdirichlet(1, alpha[[i]] )[1:ni]
      } else {
        Xik_tmp[i,1:ni]<-(alpha[[i]]/sum(alpha[[i]]))[1:ni]
      }
    }
    Xik<-cbind(Xik,Xik_tmp)
    geno_names[[i]]<-seq((m+1),(ni+m))
    m<-m+ni
  }

  #colnames(Xik)<-paste0("g",unlist(geno_names))
  #rownames(Xik)<-paste0("hab",1:np)
  
  # deme quality effects are drawn in a normal distribution with sd as defined above:
  
  genoeffects<-rnorm(ngeno, m=0, sd = sd_geno)
  demeffects<-as.numeric(Xik%*%genoeffects)
  
  
  
  if (is.null(habeffects)) {
    habeffects<-rnorm(nd, m=0, sd = sd_habs)
  } 
  
  # loop on habitat of destination
  for (i in 1:nd) {
    # loop on deme of origin
    for (j in 1:nd) {
      
      #draw the population-level error and add the "SA" contrast tosympatric transplants:
      if(i == j){error_pop <- rnorm(replicat,m = SA,sd = sd_pop)}
      if(i != j){error_pop <- rnorm(replicat,m = 0,sd = sd_pop)}
      # loop on the replicat
      for (k in 1:replicat) {
        
        tmp<-matrix(NA,geno[j],6)
        colnames(tmp)<-c("habitat","deme","genotype","replicat","symp","Y")
        tmp<-as.data.frame(tmp)
        tmp$genotype<-geno_names[[j]]
        tmp$deme<-j
        tmp$habitat<-i
        tmp$replicat<-k+i*10-10
        
        if(i == j){tmp$symp<-1}
        if(i != j){tmp$symp<-0}
        
        #draw error for one observation for the deme in the replicat
        error_ind = rnorm(geno[j],m = 0,sd = sd_error)
        
        tmp$Y<-mu + habeffects[i] + demeffects[j] + error_pop[k] + error_ind
        
        res<-rbind(res,tmp)

      }
 
    }
  }
  return(list(W=res,freq=Xik))
}
  
  
  
  
