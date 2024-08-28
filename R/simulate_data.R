#' Simulate data with gene sets
#'
#' @param  pPos Probability of positive coefficient for indicator variables
#' @param  n  Total number of samples.
#' @param  pfh The number of features in the functional hub
#' @param  pnh The number of features in each non-functional hub.
#' @param  mnh The number of non-functional hubs
#' @param  pfi The number of features that are functional independents.
#' @param  pni The number of features that are nonfunctional independents.
#'         Note: The total number of features will be denoted by p where
#'             p=pfh+mnh*pnh+pfi+pni.
#' @param  SNRh Hub signal to noise ratio. Number between 0 and 1  (1 indicates high signal)
#' @param  SNRo Signal to noise ratio for a functional independent. Number between 0 and 1 (1 indicates high signal)
#' @param  rFHub Correlation parameter between features within the functional hub
#' @param  rNFHub Correlation between features corresponding to the nonfunctional
#'              hub and each of their respective hubs (see xNFHub below)
#'              note: for the time being, betaNFHub is the same for all
#'                    non functional hubsalpha = correlation parameter between xhub values (see below) and Y.
#' @param sdNoise Standard deviation for the non-functional independents
#' @param vrb Verbose flag to enable printing of intermediate steps
#'
#' @author Lorin Towle-Miller
#'
#' @description
#' \code{simulate_data} simulates data with gene sets for different relationships. A detailed description is given in Miecznikowski et.al.
#'
#' @details
#' Simulates the following gene sets:
#' - Functional hub (genes prefixed with "fh"): Genes that are correlated with each other and the outcome
#' - Non-functional hubs (genes prefixed with "nh"): Genes that are correlated with each other but not the outcome
#' - Functional independents (genes prefixed with "fi"): Genes that are correlated with the outcome but not each other
#' - Non-functional independents (genes prefixed with "ni"): Genes that are not correlated with each other or the outcome. Independent noise.
#'
#' @return{ object of class "simHub" containing
#'  x = p x n data matrix (see below).
#'  y = n x 1 vector.
#' featureID=px1 character vector with values ('fh','nh','fi','ni').
#'
#'
#'
#'
#'
#'NOTE: the rows of x are ordered as follows:
#'
#' functional hub features      1:pfh
#'
#' non-functional hub features  (pfh+1):((pfh+1)+mnh*pnh)
#'
#' functional independents      ((pfh+1)+mnh*pnh+1):(((pfh+1)+mnh*pnh+1)+pfi)
#'
#' non-functional independents  (((pfh+1)+mnh*pnh+1)+pfi+1):p
#'
#' standard deviation of background noise in y is 1.0}
#'
#'
#'@references Miecznikowski JC, Gaile DP, Chen X, Tritchler DL.
#'  Identification of consistent functional genetic modules.
#'  \emph{Stat Appl Genet Mol Biol.} 2016 Mar;15(1):1-18. doi:
#'  10.1515/sagmb-2015-0026.
#'
#'
#' @export
simulate_data <- function(pPos=1, #probabily of positive x, hub association
                      n=100,  # total number of samples
                      pfh=10, # the number of features in the one functional hub
                      pnh=10, # the number of features in each non-functional hub
                      mnh=5, # the number of non-functional hubs
                      pfi=10, # the number of features that are functional independents
                      pni=100, # the number of features that are nonfunctional independents
                      #         Note: So, the total number of features will be denoted by p where
                      #             p=pfh+mnh*pnh+pfi+pni
                      SNRh=.5, # hub signal to noise ratio
                      SNRo=.5, # signal to noise ratio for a functional independent
                      rFHub=.3, # correlation parameter  between x rows corresponding to the functional
                      #              hub and xHub
                      rNFHub=.3, # correlation between x rows corresponding to the nonfunctional
                      #              hub and each of their respective hubs (see xNFHub below)
                      #              note: for the time being, betaNFHub is the same for all
                      #                    non functional hubs
                      sdNoise = 1, # Initial standard deviation for x
                      vrb=F # a verbose flag to enable printing of intermediate steps


)
{
  alphaHub = SNRh   # regression parameter between xhub value (see below) and Y

  alphaInd = SNRo/sqrt(pfi)  # regression parameter between a functional independent and Y
  # alphaInd = SNRo  # regression parameter between a functional independent and Y

  betaFHub = sqrt(rFHub/(1 - rFHub))  # regression parameter  between x in same functional module

  betaNFHub = sqrt(rNFHub/(1 - rNFHub))
  #             independent feature and y
  datout <- sim1typeData( pPos, n,  # total number of samples
                          pfh, # the number of features in the one functional hub
                          pnh, # the number of features in each non-functional hub
                          mnh, # the number of non-functional hubs
                          pfi, # the number of features that are functional independents
                          pni, # the number of features that are nonfunctional independents
                          #         Note: So, the total number of features will be denoted by p where
                          #             p=pfh+mnh*pnh+pfi+pni
                          alphaHub, # regression parameter between xhub value (see below) and Y
                          alphaInd, # regression parameter between a functional
                          #             independent feature and y
                          betaFHub, # regression parameter  between x in same functional module
                          betaNFHub, # regression parameter between x in same non-functional module
                          sdNoise,
                          vrb=F # a verbose flag to enable printing of intermediate steps
  )
}

#===============================================================================================#
#
# sim1typeData - Function to simulate data under the Hub models
#                   based on correlations between functional genes and y
#  {add documentation later}
#
# Arguments :
#      pPos = probability of positive coeficient for module indicator variables
#      n  = total number of samples
#      pfh = the number of features in the one functional hub
#      pnh = the number of features in each non-functional hub
#      mnh = the number of non-functional hubs
#      pfi = the number of features that are functional independents
#      pni = the number of features that are nonfunctional independents
#         Note: So, the total number of features will be denoted by p where
#             p=pfh+mnh*pnh+pfi+pni
#
#   alphaHub = correlation parameter between xhub values (see below) and Y
#   alphaInd = correlation parameter between x rows corresponding to the functional
#             independent features and y
#   betaFHub = correlation parameter  between x rows corresponding to the functional
#              hub and xHub
#   betaNFHub = correlation between x rows corresponding to the nonfunctional
#              hub and each of their respective hubs (see xNFHub below)
#              note: for the time being, betaNFHub is the same for all
#                    non functional hubs
#
#
# Output : List object of class "FMDsimHub" containing
#      xHub = 1 x n vector of the values for the latent hub across samples
#      xNFHub=  mnh x n matric of values for the latent nonfunctional hubs
#      x = p x n data matrix (see below)
#      y = n x 1 vector
#      simpars= list of arguments given above
#      featureID=px1 character vector with values ('fh','nh','fi','ni')
#
#  comments regarding the data matrix x:
#     the rows of x are ordered as follows:
#            functional hub features      1:pfh
#            non-functional hub features  (pfh+1):((pfh+1)+mnh*pnh)
#            functional independents      ((pfh+1)+mnh*pnh+1):(((pfh+1)+mnh*pnh+1)+pfi)
#            non-functional independents  (((pfh+1)+mnh*pnh+1)+pfi+1):p
#
#    standard deviation of background noise in y is 1.0
#===============================================================================================#
# Non-export
sim1typeData=function(
    pPos, # probability of positive coeficient for module indicator variables
    n,  # total number of samples
    pfh, # the number of features in the one functional hub
    pnh, # the number of features in each non-functional hub
    mnh, # the number of non-functional hubs
    pfi, # the number of features that are functional independents
    pni, # the number of features that are nonfunctional independents
    #         Note: So, the total number of features will be denoted by p where
    #             p=pfh+mnh*pnh+pfi+pni
    alphaHub, # correlation parameter between xhub values (see below) and Y
    alphaInd, # correlation parameter between x rows corresponding to the functional
    #             independent features and y
    betaFHub, # correlation parameter  between x rows corresponding to the functional
    #              hub and xHub
    betaNFHub, # correlation between x rows corresponding to the nonfunctional
    #              hub and each of their respective hubs (see xNFHub below)
    #              note: for the time being, betaNFHub is the same for all
    #                    non functional hubs
    #
    sdNoise = 1, # Initial standard deviation for x
    vrb=F # a verbose flag to enable printing of intermediate steps

    # standard deviation of background noise in y is 1.0
){


  p=pfh+mnh*pnh+pfi+pni

  x=matrix(rnorm(n*p,0,sdNoise),p,n)
  xNFHub=matrix(NA,mnh,n)
  y=rep(NA,n)
  featureID=rep("",p)
  featureID[1:pfh]="fh"
  featureID[(pfh+1):(pfh+mnh*pnh)]="nh"
  fi.index <- (pfh+mnh*pnh+1):(pfh+mnh*pnh+pfi)
  featureID[fi.index]="fi"
  featureID[(pfh+mnh*pnh+pfi+1):p]="ni"


  # simulate functional hub data

  betaSign = ifelse(rbinom(pfh,1, pPos),1,-1)
  xHub=rnorm(n,0,1)
  for(i in 1:pfh) x[i,]=betaSign[i]*betaFHub*xHub+rnorm(n,0,1)

  # simulate nonfunctional hub data

  if ( mnh > 0) {
    for(h in 1:mnh){
      betaSign = ifelse(rbinom(pnh,1, pPos),1,-1)
      xNFHub[h,]=rnorm(n,0,1)
      for(i in 1:pnh) x[(pfh+i+(h-1)*pnh),]=betaSign[i]*betaNFHub*xNFHub[h,]+rnorm(n,0,1)
    }
  }

  # generate y values based on correlations
  # suppose y = alphaHub*xHub + e where e~N(0,1)
  # let cor.xHub: correlation between xHub and y
  # suppose y = alphaInd*xInd + e
  # let cor.xInd: correlation between functional independents and y
  # then
  cor.xHub <- alphaHub/sqrt(alphaHub^2+1)
  cor.xInd <- alphaInd/sqrt(alphaInd^2+1)
  # and we can simulate y as:
  y <- cor.xHub*xHub + sqrt(1-cor.xHub^2)*rnorm(n)

  # simulate functional independents based on cor.xInd

  for(i in 1:pfi) x[fi.index[i],] <- cor.xInd*y + sqrt(1-cor.xInd^2)*rnorm(n) + rnorm(n)

  x_rownames <- featureID
  x_rownames[featureID=="nh"][1:pnh] <- "nh1"
  for(hub in 1:mnh){
    if(hub==1){
      indices <- 1:pnh
    } else {
      indices <- (indices[length(indices)] + 1):(indices[length(indices)] + pnh)
    }
    x_rownames[featureID=="nh"][indices] <- paste0("nh", hub)
  }
  x_rownames <- paste0(x_rownames, "_", 1:length(x_rownames))

  rownames(x) <- x_rownames
  colnames(x) <- paste0("s", 1:ncol(x))


  simpars <- list(n,pfh,pnh,mnh,pfi,pni,alphaHub,alphaInd,betaFHub,betaNFHub)
  lists <- list(x = x, y = y, featureID=featureID)

  class(lists) <- "simHub"
  return(lists)
}
