# adapted from Blanquart et al 2013
#################################################################################
################################
# FUNCTION TO GENERATE DATA THAT FOLLOWS PERFECTLY THE ASSUMPTIONS OF THE LINEAR MODEL #
#################################################################################
library(gtools)

################################
make_data <- function(np, d, ind, sd_ind, sd_pop, sd_hab, sd_geno, SA, constant,useDir=FALSE,alpha)
{
  #PARAMETER DEFINITIONS:
  # np is the number of populations tested
  # d indicates design, from d = 2 (twice as many allopatric as sympatric transplant), to np - 1 (full factorial design)
  # ind is the number of individuals tested by transplant
  # sd_ind is the individual - level standard deviation of error
  # sd_pop is the population - level standard deviation of error
  # sd_hab is the standard deviation of habitat effects
  # sd_geno is the standard deviation of geno quality effects --> used to calculate demeffects
  # SA is the sympatric vs. allopatric contrast (i.e., local adaptation)
  # constant is the baseline fitness
  # alpha is the vector of alpha values used to sample the genotype frequencies per habitat using a Dirichlet distribution if useDir=TRUE
  # it correspond to the relative weight of genotypes for frequency calculation
  # else the genotype frequencies are not allowed to vary over simulations and corresponds to alpha/sum(alpha)
  # if (length(alpha)=ind) the sum of the genotype frequencies per habitat are equal to 1 and if (length(alpha)>ind), virtual genotypes 
  # are created for the calcul of the frequencies (and thus the sum is < 1) but not used for the rest of the data making
  
  Xik<-matrix(0,np,ind*np)
  for (i in 1:np) {
    if (useDir==TRUE) {
      
      Xik[i,(i*ind-ind+1):(i*ind)]<-rdirichlet(1, alpha )[1:ind]
      
    } else {
      Xik[i,(i*ind-ind+1):(i*ind)]<-(alpha/sum(alpha))[1:ind]
    }
  }
  
  #habitat and deme quality effects are drawn in a normal distribution with sd as defined above:
  habeffects<-rnorm(np,m = 0,sd = sd_hab)
  genoeffects<-rnorm(np*ind,m=0, sd = sd_geno)
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


