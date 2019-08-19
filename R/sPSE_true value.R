#' General approach of causal mediation analysis for survival outcome under sequential mediators.
#'
#' Derive the true iPSEsv based on real parameters.
#'
#' @name iPSEsv_t
#' @author An-Shun Tai \email{daansh13@gmail.com}, Pei-Hsuan Lin \email{a52012232@gmail.com}, and Sheng-Hsuan Lin \email{shenglin@nctu.edu.tw}
#' @param mediators The name of mediators.
#' @param parameters The real values of parameters
#' @return A list of true iPSEsv and partPSEsv
#' @export
#' @examples
#' TrueV <- iPSEsv_t(mediators=c("x1","x2"),
#' alpha_Y=matrix(1,1,1),beta_Y=matrix(1),gamma_Y=matrix(c(0.5,0.6),1,2),
#' delta_Y=matrix(c(0.5,0.6),1,2),
#' alpha_M=matrix(1,2,2),beta_M=matrix(1,2,1),gamma_M=matrix(1,2,2),
#' delta_M=matrix(1,2,1),sigma2_M=matrix(1,2,1),
#' alpha_C=matrix(1,2,2),beta_C=matrix(1,2,1),gamma_C=matrix(1,2,2),
#' delta_C=matrix(1,2,1),sigma2_C=matrix(1,2,1),
#' C_hat = c(1,0.5),
#' a1=1,a0=0)

iPSEsv_t <- function(mediators,
                   alpha_Y,beta_Y,gamma_Y,delta_Y,
                   alpha_M,beta_M,gamma_M,delta_M,sigma2_M,
                   alpha_C,beta_C,gamma_C,delta_C,sigma2_C,
                   C_hat = c(0.5),
                   a1=1,a0=0){


  ########################
  # create the matrix of the exposure status
  Ematrix <- t(apply(as.matrix(1:2^length(mediators)),1,function(x)c(rep(a1,x),
                                                                     rep(a0,2^length(mediators)-x))))
  Ematrix <- rbind(0,Ematrix)

  ########################
  #
  parH <- c()
  parH$alpha_M <- alpha_M; parH$alpha_C <- alpha_C
  parH$C_hat <- C_hat
  parH$beta_M <- beta_M; parH$beta_C <- beta_C
  parH$gamma_M <- gamma_M; parH$gamma_C <- gamma_C
  parH$delta_M <- delta_M; parH$delta_C <- delta_C


  ############ a section of generating function
  NK <- length(mediators)

  my.new.env <- new.env()
  assign("mu_M_1",function(A,i,parH) {parH$alpha_M[1,]%*%parH$C_hat+parH$beta_M[1,1]*A[1]+
      parH$gamma_M[1,1]*(parH$alpha_C[1,]%*%parH$C_hat+parH$beta_C[1,1]*A[1])},my.new.env)
  assign("mu_C_1",function(A,i,parH) {parH$alpha_C[1,]%*%parH$C_hat+parH$beta_C[1,1]*A[1]},my.new.env)
  for(p in 2:NK){
    assign(paste("mu_M_",p,sep=""),function(A,i,parH){
      s1 <- 0
      for(h in 1:i){
        s1 <- parH$gamma_M[i,h]*getFunction(paste("mu_C_",h,sep=""),where = my.new.env)(A[1:2^(h-1)],i=h,parH) +s1
      }

      s2 <- 0
      for(h in 1:(i-1)){
        s2 <- parH$delta_M[i,h]*getFunction(paste("mu_M_",h,sep=""),where = my.new.env)(A[(2^(h-1)+1):2^h],i=h,parH) +s2
      }
      s <- parH$alpha_M[i,]%*%parH$C_hat+parH$beta_M[i,1]*A[1] + s1+s2
      return(s)
    },my.new.env)
  }# for M

  for(p in 2:NK){
    assign(paste("mu_C_",p,sep=""),function(A,i,parH){
      s1 <- 0
      for(h in 1:(i-1)){
        s1 <- parH$gamma_C[i,h]*getFunction(paste("mu_C_",h,sep=""),where = my.new.env)(A[1:2^(h-1)],i=h,parH) +s1
      }

      s2 <- 0
      for(h in 1:(i-1)){
        s2 <- parH$delta_M[i,h]*getFunction(paste("mu_M_",h,sep=""),where = my.new.env)(A[(2^(h-1)+1):2^h],i=h,parH) +s2
      }
      s <- parH$alpha_C[i,]%*%parH$C_hat+parH$beta_C[i,1]*A[1] + s1+s2
      return(s)
    },my.new.env)
  }#for C

  ############ END: a section of generating function

  lambda_funtion <- function(E.status){
    #tau2_M_1 <- function(A,i) {sigma2_M[1,1]+gamma_M[1,1]^2*sigma2_C[1,1]}
    #Because the term of variance is absent in the calculation of causal effect, we ignore it here.

    R <- matrix(NA,length(mediators),1)
    R[length(mediators),1] <- gamma_Y[1,length(mediators)]
    for(j in (length(mediators)-1) : 1 ){
      dd <- (j+1):length(mediators)
      R[j,1] <- gamma_Y[1,j]+sum(R[dd,1]*gamma_C[dd,j])
    }

    # a temporary version
    K <- length(mediators); Z <- matrix(NA,K,1)
    if(K==1){
      Z[K,1] <- delta_Y[1,K]
    }
    if(K==2){
      Z[K,1] <- delta_Y[1,K]
      Z[K-1,1] <- delta_Y[1,K-1]+gamma_Y[1,K]*delta_C[K,K-1]
    }
    if(K==3){
      Z[K,1] <- delta_Y[1,K]
      Z[K-1,1] <- delta_Y[1,K-1]+gamma_Y[1,K]*delta_C[K,K-1]
      Z[K-2,1] <- delta_Y[1,K-2]+gamma_Y[1,K]*delta_C[K,K-2]+
        gamma_Y[1,K]*gamma_C[K,K-1]*delta_C[K-1,K-2]+gamma_Y[1,K-1]*delta_C[K-1,K-2]
    }
    if(K==4){
      Z[K,1] <- delta_Y[1,K]
      Z[K-1,1] <- delta_Y[1,K-1]+gamma_Y[1,K]*delta_C[K,K-1]
      Z[K-2,1] <- delta_Y[1,K-2]+gamma_Y[1,K]*delta_C[K,K-2]+
        gamma_Y[1,K]*gamma_C[K,K-1]*delta_C[K-1,K-2]+gamma_Y[1,K-1]*delta_C[K-1,K-2]
      Z[K-3,1] <- delta_Y[1,K-3]+gamma_Y[1,K]*delta_C[K,K-3]+
        gamma_Y[1,K]*gamma_C[K,K-1]*delta_C[K-1,K-3]+gamma_Y[1,K]*gamma_C[K,K-2]*delta_C[K-2,K-3]+
        gamma_Y[1,K]*gamma_C[K,K-1]*gamma_C[K-1,K-2]*delta_C[K-2,K-3]+
        gamma_Y[1,K-1]*delta_C[K-1,K-3]+gamma_Y[1,K-1]*gamma_C[K-1,K-2]*delta_C[K-2,K-3]+
        gamma_Y[1,K-2]*delta_C[K-2,K-3]
    }

    Zmu <- apply(matrix(1:K),1,
                 function(x)Z[x,]*getFunction(paste("mu_M_",x,sep=""),where = my.new.env)(E.status[(2^(x-1)+1):(2^x)],x,parH))

    lambda <- beta_Y*E.status[1] + sum(R*beta_C)*E.status[1] +
      sum(( c(0,alpha_Y) + alpha_C%*%R )*C_hat) + sum(Zmu)
    return(lambda)
  }

  iPSEsv <- as.matrix(apply(matrix(1:2^length(mediators)),1,
                            function(i)lambda_funtion(Ematrix[i+1,])-lambda_funtion(Ematrix[i,])))

  convert_to_binary <- function(n) {
    bin <- ''
    if(n > 1) {
      bin <- convert_to_binary(as.integer(n/2))
    }
    bin <- paste0(bin, n %% 2)
    return(as.character(bin))
  }
  NN <- c(); NN[1] <- "A->Y"
  for(u in 2:2^length(mediators)){
    aaa <- as.numeric(strsplit(convert_to_binary(u-1),"")[[1]])
    NN[u] <- paste("A",paste(mediators[which(aaa[length(aaa):1]==1)],collapse = "->"),"Y",sep="->")[1]
  }
  rownames(iPSEsv) <- NN
  colnames(iPSEsv) <- "PSE"
  #######################################################
  ########################################################################################
  part_lambda_funtion <- function(){
    R0 <- matrix(NA,length(mediators),1)
    R0[length(mediators),1] <- delta_Y[1,length(mediators)]
    for(j in (length(mediators)-1) : 1 ){
      dd <- (j+1):length(mediators)
      R0[j,1] <- delta_Y[1,j]+sum(R0[dd,1]*delta_C[dd,j])
    }
    return(R0*beta_M)
  }
  partPSEsv <- rbind( beta_Y[1,],part_lambda_funtion() )*(a1-a0)

  NN0 <- c(); NN0[1] <- "A->Y"
  for(u in 1:length(mediators)){
    NN0[u+1] <- paste("A",paste(c(mediators[u],"..."),collapse = "->"),"Y",sep="->")[1]
  }
  rownames(partPSEsv) <- NN0
  colnames(partPSEsv) <- "true PSE"

  #output
  return(structure(list(iPSEsv=iPSEsv,partPSEsv=partPSEsv)))
}
