
#-----------------------------------------------------------------------------------------------------------------------#
# Functions for bootstrap estimation of Confidence Interval on a mixed model


mySumm <- function(.) {
  c(beta=fixef(.), sig01=sqrt(unlist(VarCorr(.))))
}

bCI.tab <- function(b,ind=length(b$t0), type="perc", conf=0.95) {
  btab0 <- t(sapply(as.list(seq(ind)),
                    function(i)
                      boot.ci(b,index=i,conf=conf, type=type)$percent))
  btab <- btab0[,4:5]
  rownames(btab) <- names(b$t0)
  a <- (1 - conf)/2
  a <- c(a, 1 - a)
  pct <- stats:::format.perc(a, 3)
  colnames(btab) <- pct
  return(btab)
}



#-----------------------------------------------------------------------------------------------------------------------#
# Generic Function to simulate the data 
# onlu works for glmer(Y~X1+...+Xn+(1|a1)+...+(1|ar)+(1|X1:a1)+...+(1|Xn:ar))
# ie no nested random effect within fixed effect eg. (0+X1|a1) --> TO DO
# return simulated dataset
# has to be call inside a sapply 
# parameters to set:  
  
# - npods --> number of grains per pods, can be an integer or a vector of integers for each genotype (in the genotype levels order)  
# - varpods --> variance of the number of seeds per pods, if npods is a single integer, model an inter-genotypic variation and if npods is a vector of integers for each genotype, model an intra genotype variation. Default value is 0.
# The user must choose between extract the parameters from the model run on real data or give parameters:  
#   - mod --> model used to extract coefficients and variance  
# - scenario --> list of x elements :  
# * $design --> dataframe corresponding to the design of the data  
# * $fixed --> list of coefficients for fixed effects, (can be set as log of odd or odd ratio if the link function used is logit for example if the values are at the data scale). The first element must be called "(Intercept)", 
# and the other coef must be in the order of levels(), see examples.  
# ex fixed.ef<-list("(Intercept)"=-2.554329103,"zone"=c("3"=-0.489371384 ,"8"=0.006808232,"interpop"=-0.533948813 ))
# * $random --> vector with variances for random effects names for simple effect must correspond to names in the design
# matrix and names for interaction effect must corresponds to eg. X1:a1 or a1:a2, see example on scenario
# * $linkfamily --> the link family --> "logit"

# function to be call inside a apply function
simgerm_apply<-function(npods=6,varpods=0,factpods=NULL,mod=NULL,scenario=NULL) {
  
  # take coefficients and variance from the model run on real dataset
  if (!is.null(mod)) {
    data_fact<-model.frame(mod)[,2:ncol(model.frame(mod))]
    fixed.ef<-fixef(mod)
    data_fact$fixed<-model.matrix(mod)%*%fixed.ef
    rand.var<-as.data.frame(VarCorr(mod))$vcov
    names(rand.var)<-as.data.frame(VarCorr(mod))$grp
    linkfunction<-family(mod)$link # must be logit
    
  } else if (!is.null(scenario)){
    data_fact<-scenario$design
    fixed.ef<-scenario$fixed
    if (length(fixed.ef)>1){
      fix.names<-names(fixed.ef)[2:length(fixed.ef)]
      form<-as.formula(paste("~",paste(fix.names,collapse="+")))
    } else {
      form<-as.formula("~1")
    }
    don<-model.matrix(form,data_fact)
    # check if it works on sevral fixed effects
    data_fact$fixed<-don%*%unlist(fixed.ef)
    rand.var<-scenario$random
    linkfunction<-scenario$linkfamily # must be logit
  } else {
    stop("A model or a scenario must be given")
  }
  # random effect
  if (!is.null(rand.var)) {
    rand.names<-names(rand.var)
    # identify interaction factors
    rand.var_interact<-rand.names[!rand.names%in%colnames(data_fact)]
    # create a list whith for each element the simple effect forming each interaction term
    if (length(rand.var_interact)>1) {
      colsup<-lapply(strsplit(rand.var_interact,":"),function(x) {
        as.character(Reduce(function(x,y)interaction(x,y,sep=":"),as.list(data_fact[,colnames(data_fact)%in%x])))
      })
      names(colsup)<-rand.var_interact
      
      data_fact<-cbind(data_fact,do.call(cbind,colsup))
    }
    
    # sample random effects : produce a list, each element correspond to a sample of random effect (blup) rn with a length = nlevels(data_fact[,rn])
    rand.list<-lapply(rand.names,function(rn) {
      data_fact[,rn]<-factor(data_fact[,rn])
      nr<-nlevels(data_fact[,rn])
      rnorm(nr,0,sqrt(rand.var[rn]))
    })
    names(rand.list)<-rand.names
  # set the random effect for each observation based on the design matrix obtained with model.matrix %*% blup 
  rand.obs<-lapply(rand.names,function(x) {
    frand<-as.formula(paste("~-1+",x))
    don<-model.matrix(frand,data_fact)
    # for interaction, needs suppress non existing interaction levels
    if (length(which(colSums(don)==0)>0)) {
      don<-don[,-which(colSums(don)==0)]
    }
    don%*%rand.list[[x]]
  })
  names(rand.obs)<-paste0("r",rand.names) 
  # 3 steps necessary because the colnames are not ok in one-step
  tmp<-do.call(cbind,rand.obs)
  colnames(tmp)<-names(rand.obs)
  data_fact<-cbind(data_fact,tmp)
  }

  # calculate the log odd of probabilities
  data_fact$logodd<-apply(data_fact[,c("fixed",paste0("r",rand.names) )],MARGIN=1,FUN=sum)
  
  # conditions on the number of pods
  if (length(npods)==1) {
    if (varpods==0) {
      data_fact$nseeds<-rep(npods,nrow(data_fact))
    } else if (!is.null(factpods)){
      nseeds<-rnorm(n=nlevels(data_fact[,factpods]),mean=npods,sqrt(varpods))
      
      data_fact$nseeds<-don2%*%nseeds
    } else {
      stop("the factor used to sample random numbers of seeds per pod around the mean is needed")
    }
  } else if ((!is.null(factpods))&(length(npods)==nlevels(data_fact[,factpods]))){
    if (varpods==0) {
      don_pods<-model.matrix(as.formula(paste("~-1+",factpods)),data_fact)
      data_fact$nseeds<-don_pods%*%npods[levels(data_fact[,factpods])]
    } else {
      ni<-tapply(data_fact[,factpods],list(data_fact[,factpods]),length)
      nseeds<-unlist(mapply(FUN=rnorm,ni,mean=npods,MoreArgs = list(sd=sqrt(varpods))))
      data_fact$nseeds<-nseeds
    }
  } else {
    stop("npods must be a single integer or a vector if integers corresponds to the values for each factpods level.")
  }
  # survivors
  data_fact$Y<-rbinom(n=1:nrow(data_fact),size=round(data_fact$nseeds),prob=plogis(data_fact$logodd))
  # deads
  data_fact$D<-round(data_fact$nseeds)-data_fact$Y
  
  return(data_fact)
}




#-----------------------------------------------------------------------------------------------------------------------#

# TO DO : function to simulate binomial data from a design (drop the part with npods)

