#' Build model for multiple systems estimation
#'
#' For multiple systems estimation model corresponding to a specified set of two-list effects,
#'    set up the GLM model formula and data matrix.
#'
#'
#'@param zdat  Data matrix with \eqn{t+1} columns. The first \eqn{t} columns, each corresponding to a particular list,
#' are 0s and 1s defining the capture histories
#' observed. The last column is the count of cases with that particular capture history.
#' List names A, B, ... are constructed if not supplied. Where a capture history is not explicitly listed,
#' it is assumed that it has observed count zero.
#'
#'@param mX A \eqn{2 \times k} matrix giving the \eqn{k} two-list interactions to be included in the model.
#' Each column of \code{mX} contains the numbers of the corresponding pair of lists.
#' If \code{mX = 0}, then all two-list interactions are included. If \code{mX = NULL}, no interactions are included and
#'  the main effects model is fitted.
#'
#' @return A list with components as below.
#'
#'  \code{datamatrix} A matrix with all possible capture histories, other than those corresponding to empty overlaps
#'  within the model.  An empty overlap is a pair of lists \eqn{(i,j)} such that no case is observed in both lists,
#'  regardless of whether it is present on any other lists.   If \eqn{(i,j)} is within the model specified by \code{mX},
#'  all capture histories containing both \eqn{i} and \eqn{j} are
#'  then excluded.
#'
#' \code{modelform} The model formula suitable to be called by the Generalized Linear Model function \code{glm}.
#'
#' \code{emptyoverlaps} A matrix with two rows, whose columns give the indices of pairs of lists for which there are empty overlaps and where
#' the pair is within the specified model. The column names give the names of the lists
#' corresponding to each pair.
#'
#' @examples
#' data(NewOrl)
#' buildmodel(NewOrl, mX=NULL)
#'
#' @importFrom stats as.formula
#'
#'@export


buildmodel <- function(zdat, mX) {
  m = dim(zdat)[2] - 1
  #Recover the data matrix with combination of lists which have count zero
  zdfull = tidylists(zdat, includezerocounts = T)
  vn = dimnames(zdfull)[[2]]
  #count
  cn = vn[m+1]
  #names of list excluding count
  vn = vn[ -(m+1)]
  modf = paste(cn, "~ .")
  if (is.null(mX)) return(list(datamatrix = zdfull, modelform = modf, emptyoverlaps = matrix(nrow=2,ncol=0)))
  if (sum(mX) == 0) {
    mX = t(expand.grid(1:m, 1:m))
    mX = mX[ , mX[1,]<mX[2,]]
  }
  mX=as.matrix(mX)
  # ---  find if any of the overlaps specified within the model are empty
  neff = dim(mX)[2]
  effc = rep(T,neff)
  for (j in (1:neff)) {
    i1 = mX[1,j]
    i2 = mX[2,j]
    joverlap = sum(zdfull[,i1]*zdfull[,i2]*zdfull[,m+1])
    if (joverlap == 0) {
      effc[j] = F
      zdfull = zdfull[zdfull[,i1]*zdfull[,i2]==0, ]
    }
  }
  #----construct matrix of interaction parameters where overlap is empty
  mempty = as.matrix(mX[, !effc])
  dimnames(mempty)[[1]] = c("","")
  menames =  apply(matrix(vn[mempty],nrow=2),2, paste, collapse=":")
  dimnames(mempty)[[2]] = menames
  #----construct matrix of interaction parameters where overlap is not empty
  mX = as.matrix(mX[,effc])
  mXc = vn[mX]
  dim(mXc) = dim(mX)
  #----construct the model formula to include
  #      the interactions in the model for which the overlap is not empty
  mXmod = paste(apply(mXc,2, paste, collapse=":"), collapse=" + ")
  if (nchar(mXmod) > 0) modf = paste(modf, mXmod, sep= " + ")
  modf = as.formula(modf)
  return(list(datamatrix = zdfull, modelform = modf, emptyoverlaps = mempty))
}

#' Check all possible models
#'
#'This routine checks every possible model for existence and identifiability of the maximum likelihood estimates.
#'
#'  The routine calls the routine \code{\link{checkident}} for every model.
#'  If there are \eqn{t} lists then there are \eqn{t(t-1)/2} pairs of lists, and hence \eqn{2^{t(t-1)/2}} possible models, because
#'   the models correspond to subsets of the set of all pairs of lists.   If \eqn{t = 7} there are 2,097,152 models to check, which would take
#'   several hours.   If \eqn{t} is equal to 8 or more, the routine terminates with a statement of the number of models and an explanation that
#'   checking all of these is not possible in a reasonable time.
#'
#'@param zdata Data matrix with \eqn{t+1} columns. The first \eqn{t} columns, each corresponding to a particular list,
#' are 0s and 1s defining the capture histories
#' observed. The last column is the count of cases with that particular capture history.
#' List names A, B, ... are constructed if not supplied. Where a capture history is not explicitly listed,
#' it is assumed that it has observed count zero.
#'
#'@param verbose Specifies whether all possible models are listed in the output, or only those which generate a non-zero error code.
#'
#'@return
#' If \code{verbose=F}, it gives a matrix each of whose rows specifies a model, with last entry equal to the error code.  Only those models yielding a non-zero error code are included, so if no models lead to an error the matrix is empty. Each of the first \eqn{t(t-1)/2} columns corresponds to a pair of lists, and for each row, presence in or absence from the corresponding model is indicated by the value 1 or 0 respectively.
#'
#' If \code{verbose=T}, it gives the full matrix of models together with the error codes they generate.
#'
#'
#' @examples
#'data(Artificial_3)
#'data(Western)
#'checkallmodels(Artificial_3, verbose= TRUE)
#'checkallmodels(Western)
#'@export
checkallmodels <- function (zdata, verbose=F){
    m = dim(zdata)[2] - 1
    if (m == 6)
      print("The program is checking 32728 models.")
    if (m == 7)
      print("There are 2,097,152 models to check. This may take several hours.")
    if (m >= 8) {
      cat("The number of models is", 2^{
        m*(m - 1)/2
        }, ". It is impossible to check them in a reasonable time.")
      return(invisible())
      }
    mX0 = rbind(rep(1:m, each = m), rep(1:m, times = m))
    mX0 = mX0[, (mX0[1, ] < mX0[2, ])]
    varnames = dimnames(zdata)[[2]][1:m]
    mvnames = rbind(varnames[mX0[1,]], varnames[mX0[2,]])
    pairnames = apply(mvnames, 2, paste, collapse=":")
    npairs = m * (m - 1)/2
    nmodels = 2^npairs
    modmat = matrix(0, nrow=nmodels, ncol=npairs)
    dimnames(modmat)[[2]] = pairnames
    ierr = rep(NA, nmodels)
    for (j in (1:nmodels)) {
      jkeep = (intToBits(j)[1:m] == 1)
      mX = mX0[, jkeep]
      modmat[j,jkeep] = 1
      ierr[j] = checkident(zdata, mX)
      }
    resmat = cbind(modmat, ierr)
    if (!verbose) resmat= resmat[ierr>0,]
    return(resmat)
}




#' Produce a data matrix with a unique row for each capture history
#'
#' This routine finds rows with the same capture history and consolidates them into a single row whose count is the sum of counts of
#' the relevant rows.  If \code{includezerocounts = T} then it also includes rows for all the capture histories with zero counts; otherwise
#' these are all removed.
#'
#' @param zdat Data matrix with \eqn{t+1} columns. The first \eqn{t} columns, each corresponding to a particular list,
#' are 0s and 1s defining the capture histories
#' observed. The last column is the count of cases with that particular capture history.
#' List names A, B, ... are constructed if not supplied. Where a capture history is not explicitly listed,
#' it is assumed that it has observed count zero.
#'
#' @param includezerocounts  If \code{F} then remove rows corresponding to capture histories with zero count.
#' If \code{T} then include all possible capture histories including those with zero count,
#' excluding the all-zero row corresponding to the dark figure.
#'
#' @return A data matrix in the form specified above, including all capture histories with zero counts if  \code{includezerocounts=T}.
#'
#' @examples
#' data(NewOrl)
#' zdat<-tidylists(NewOrl,includezerocounts=TRUE)
#'@export
tidylists <- function(zdat, includezerocounts = F) {
  zdat = as.matrix(zdat)
  m = dim(zdat)[2] - 1
  # ----construct full capture history matrix
  zm = NULL
  #-----produce an unordered matrix of all possible capture history including the one corresponds to dark figure
  for (j in (1:m)) zm = rbind(cbind(1, zm), cbind(0, zm))
  # ----calculate the number of 1s in each row of the unordered matrix
  ztot = apply(zm, 1, sum)
  #----order the rows of the matrix according to the number of 1s appearing on each row in ascending order
  zm = zm[order(ztot), ]
  #---remove the row that has 0s for all its entries and add a column of 0s
  zm = cbind(zm[-1, ], 0)
  #---supply column names if necessary
  vn = dimnames(zdat)[[2]]
  if(is.null(vn)) {
    vn = c(LETTERS[1:m], "count")
  }
  dimnames(zm)[[2]] = vn
  # find row of zm corresponding to each row of zmse and update count
  bcode = zdat[, -(m + 1)] %*% (2^(1:m))
  bc = zm[, -(m + 1)] %*% (2^(1:m))
  for (j in (1:dim(zdat)[1])) {
    ij = (1:dim(zm)[1])[bc == bcode[j]]
    zm[ij, m + 1] = zm[ij, m + 1] + zdat[j, m + 1]
  }
  # remove rows with zero counts if includezerocounts is F
  if (!includezerocounts)
    zm = zm[zm[, m + 1] > 0, ]
  zdatr = as.data.frame(zm)
  return(zdatr)
}

#' Build the model matrix based on particular data, as required to check for identifiability and existence of the maximum likelihood estimate
#'
#' This routine builds a model matrix as required by the linear program check \code{\link{checkident}} and checks if the matrix is of full rank.
#'  In addition, for each individual list, and for each pair of lists included in the model,
#'   it returns the total count of individuals appearing on the specific list or lists whether or not in combination with other lists.
#'
#' @param zdat Data matrix with \eqn{t+1} columns. The first \eqn{t} columns, each corresponding to a particular list,
#' are 0s and 1s defining the capture histories
#' observed. The last column is the count of cases with that particular capture history.
#' List names A, B, ... are constructed if not supplied. Where a capture history is not explicitly listed,
#' it is assumed that it has observed count zero.
#'
#' @param mX  A \eqn{2 \times k} matrix giving the \eqn{k} two-list interactions to be included in the model.
#' Each column of \code{mX} contains the numbers of the corresponding pair of lists.
#' If \code{mX = 0}, then all two-list interactions are included. If \code{mX = NULL}, no interactions are included and
#'  the main effects model is fitted.
#'
#' @return A list with components as below
#'
#' \code{modmat} The matrix that maps the parameters in the model (excluding any corresponding to non-overlapping lists) to the log expected value
#'    of the counts of capture histories that do not contain non-overlapping pairs in the data.
#'
#' \code{tvec} A vector indexed by the parameters in the model, excluding those corresponding to non-overlapping pairs of lists.  For each parameter
#'    the vector contains the total count of individuals in all the capture histories that dominate that parameter.
#'
#' \code{rankdef} The column rank deficiency of the matrix \code{modmat}. If \code{rankdef = 0}, the matrix has full column rank.
#'
#'
#'@examples
#' data(NewOrl)
#' buildmodelmatrix(NewOrl, mX=NULL)
#'
#' @export
#'
buildmodelmatrix = function(zdat, mX=NULL) {
  m = dim(zdat)[2] - 1
  #Recover the data matrix
  zdfull = tidylists(zdat, T)
  if (length(mX)==0) { mX= matrix(nrow=2,ncol=0)
  }
  else {
    if (sum(mX) == 0) {
      mX = t(expand.grid(1:m, 1:m))
      mX = mX[ , mX[1,]<mX[2,]]
    }
    mX=as.matrix(mX)
    # ---  find if any of the overlaps specified within the model are empty
    neff = dim(mX)[2]
    effc = rep(T,neff)
    for (j in (1:neff)) {
      i1 = mX[1,j]
      i2 = mX[2,j]
      joverlap = sum(zdfull[,i1]*zdfull[,i2]*zdfull[,m+1])
      if (joverlap == 0) {
        effc[j] = F
        zdfull = zdfull[zdfull[,i1]*zdfull[,i2]==0, ]
      }
    }
    #----construct matrix of interaction parameters where overlap is not empty
    mX = as.matrix(mX[,effc])
  }
  counts = zdfull[,m+1]
  nobs = dim(zdfull)[1]
  vn = c(0,dimnames(zdfull)[[2]][-(m+1)])
  modmat = cbind(rep(1,nobs), zdfull[, -(m+1)])
  #----construct the model matrix to include
  #      the interactions in the model for which the overlap is not empty
  neffs = dim(mX)[2]
  if (neffs > 0) {for ( j in (1:dim(mX)[2])) {
    i1 = mX[1,j]
    i2 = mX[2,j]
    modmat = cbind(modmat, zdfull[,i1]*zdfull[,i2])
    vn = c(vn, paste(vn[i1+1],vn[i2+1],sep=":"))
  }}
  dimnames(modmat)[[2]] =vn
  tvec = t(modmat)%*%counts
  rankdef = dim(modmat)[2] - qr(modmat)$rank
  return(list(modmat=modmat, tvec=tvec, rankdef=rankdef))
}
#' Check a model for the existence and identifiability of the maximum likelihood estimate
#'
#' Apply the linear programming test as derived by Fienberg and Rinaldo (2012), and a calculation of the rank of the design
#' matrix, to check whether a particular model yields an identifiable maximum likelihood estimate
#' based on the given data.  The particular algorithm applied is described on page 3 of the supplementary material, with a typographical error corrected.
#'
#' @param zdat Data matrix with \eqn{t+1} columns. The first \eqn{t} columns, each corresponding to a particular list,
#' are 0s and 1s defining the capture histories
#' observed. The last column is the count of cases with that particular capture history.
#' List names A, B, ... are constructed if not supplied. Where a capture history is not explicitly listed,
#' it is assumed that it has observed count zero.
#'
#' @param mX  A \eqn{2 \times k} matrix giving the \eqn{k} two-list interactions to be included in the model.
#' Each column of \code{mX} contains the numbers of the corresponding pair of lists.
#' If \code{mX = 0}, then all two-list interactions are included. If \code{mX = NULL}, no interactions are included and
#'  the main effects model is fitted.
#'
#' @param verbose Specifies the output.   If \code{F} then the error code is returned.  If \code{T} then
#' in addition the routine prints an error message if the model/data fail either of the two tests, and also
#' returns both the error code and the \code{lp} object.
#'
#' @return
#'
#'If \code{verbose=F}, then return the error code \code{ierr} which is 1 if the linear program test shows that the maximum likelihood
#'  estimate does not exist, 2 if it is not identifiable, and 3 if both tests are failed.
#'
#'If \code{verbose=T}, then return a list with components as below
#'
#'\code{ierr} As described above.
#'
#'\code{zlp} Linear programming object, in particular giving the value of the objective function at optimum.
#'
#' @references Fienberg, S. E. and Rinaldo, A. (2012). Maximum likelihood estimation in log-linear
#' models. Ann. Statist. 40, 996-1023.  Supplementary material: Technical report, Carnegie Mellon University. Available from \url{http://www.stat.cmu.edu/~arinaldo/Fienberg_Rinaldo_Supplementary_Material.pdf}.
#' @export
#'
checkident <- function(zdat, mX=0, verbose=F)
{
  #---use LP to check if ML estimate exists
  #--- construct model matrix
  zb = buildmodelmatrix(zdat, mX)
  amat = t(zb$modmat)
  tt = zb$tvec
  npar = dim(amat)[1]
  nobs = dim(amat)[2]
  #----construct objective and constraint matrix
  f.obj = c(rep(0,nobs),1)
  const.rhs = c( tt, rep(0,nobs))
  amat1 = cbind(amat, rep(0,npar))
  amat2 = cbind(diag(nobs), rep(-1, nobs))
  #amat3 = c(rep(0,nobs), 1)
  const.mat = rbind(amat1, amat2)
  const.dir = c( rep("=",npar), rep(">=",nobs))
  #--------------carry out the LP step
  zlp = lpSolve::lp("max",f.obj,const.mat,const.dir,const.rhs)
  ierr = (zlp$objval == 0) + 2*(zb$rankdef > 0 )
  if (!verbose) return(ierr)
  if (zlp$objval == 0) cat ("The estimate is not finite \n")
  if (zb$rankdef > 0 ) cat ("The estimate is not identifiable \n")
  return(list(ierr=ierr, lp= zlp))
}

#' Estimate the total population including the dark figure
#'
#' This routine calculates the estimate of the total population, including the dark figure, together with confidence intervals as specified.  It also returns the details of the fitted model.
#'
#' @param zdat  Data matrix with \eqn{t+1} columns. The first \eqn{t} columns, each corresponding to a particular list,
#' are 0s and 1s defining the capture histories
#' observed. The last column is the count of cases with that particular capture history.
#' List names A, B, ... are constructed if not supplied. Where a capture history is not explicitly listed,
#' it is assumed that it has observed count zero.
#'
#' @param method If \code{method = "stepwise"} the stepwise method implemented in \code{\link{stepwisefit}} is used.  If \code{method = "fixed"} then a specified fixed model is used; the model is then given by \code{mX}.  If \code{method = "main"} then main effects only are fitted.
#'
#' @param quantiles Quantiles of interest for confidence intervals.
#' @param mX  A \eqn{2 \times k} matrix giving the \eqn{k} two-list interactions to be included in the model if \code{method = "fixed"}.
#' Each column of \code{mX} contains the numbers of the corresponding pair of lists.
#' If \code{mX = 0}, then all two-list interactions are included. If \code{mX = NULL}, no interactions are included and
#'  the main effects model is fitted.
#' If only one interaction is to be fitted, it is ok to specify it as a vector of length 2, e.g \code{mX=c(1,3)}
#' for interactions of list 1 and 3.   If \code{method} is equal to \code{"stepwise"} or \code{"fixed"} then \code{mX} is ignored.
#'
#' @param pthresh Threshold p-value used if \code{method = "stepwise"}.
#'
#' @return A list of components as below
#'
#' \code{estimate} Point estimate and confidence interval estimates corresponding to specified quantiles.
#'
#' \code{MSEfit} The model fitted to the data in the format described in \code{\link{modelfit}}.
#'
#' @examples
#' data(NewOrl)
#' data(NewOrl_5)
#' estimatepopulation(NewOrl, method="stepwise", quantiles=c(0.025,0.975))
#' estimatepopulation(NewOrl_5, method="main", quantiles=c(0.01, 0.05,0.95, 0.99))
#'
#' @importFrom stats qnorm
#'
#'@export
estimatepopulation <-function(zdat,method="stepwise", quantiles=c(0.025,0.975),mX=NULL, pthresh=0.001){
  if (method=="stepwise"){
    zMSE=stepwisefit(zdat, pthresh)
  }
  else if (method=="fixed"){
    zMSE=modelfit(zdat,mX)
  }
  else if (method=="main") {
    zMSE=modelfit(zdat, mX=NULL)
  }
  else {stop ("Unknown method")}
  zfit=zMSE$fit
  zzdat=zfit$data
  totobs = sum(zzdat[, dim(zzdat)[2]])#summing all observed history
  estdarklogmean = zfit$coefficients[1]#intercept of the fit
  estdarklogsd   = sqrt(summary(zfit)$cov.unscaled[1,1])
  quantiles_with_pt.estimate<-sort(unique(c(quantiles,0.5)))
  estdarkci = exp(qnorm(quantiles_with_pt.estimate, estdarklogmean, estdarklogsd))#estimated dark figure confidence interval
  totci = totobs + estdarkci
  names(totci)<-paste( quantiles_with_pt.estimate*100,"%",sep="")
  for (i in 1:length(totci)){
    if (names(totci)[i]=="50%"){
      names(totci)[i]<-"point est."}
  }
  return(list(estimate=totci, MSEfit = zMSE))
}

#' Fit a specified model to Multiple Systems Estimation data
#'
#' This routine fits a specified model to multiple systems estimation data, taking account of the possibility of empty overlaps between observed lists.
#'
#' @param zdat Data matrix with \eqn{t+1} columns. The first \eqn{t} columns, each corresponding to a particular list,
#' are 0s and 1s defining the capture histories
#' observed. The last column is the count of cases with that particular capture history.
#' List names A, B, ... are constructed if not supplied. Where a capture history is not explicitly listed,
#' it is assumed that it has observed count zero.

#' @param mX  A \eqn{2 \times k} matrix giving the \eqn{k} two-list interactions to be included in the model.
#' Each column of \code{mX} contains the numbers of the corresponding pair of lists.
#' If \code{mX = 0}, then all two-list interactions are included. If \code{mX = NULL}, no interactions are included and
#'  the main effects model is fitted.
#' If only one interaction is to be fitted, it may be specified as a vector of length 2, e.g \code{mX=c(1,3)}
#' for interactions of list 1 and 3.
#'
#' @param check If \code{check = T} check first of all if the maximum likelihood
#' estimate exists and is identifiable, using the routine \code{\link{checkident}}.  If either condition fails, print an appropriate error message
#' and return the error code.
#'
#'
#' @return A list of components as below
#'
#' \code{fit} Details of the fit of the specified model as output by \code{glm}.  The Akaike information criterion is adjusted to take account
#'   of the number of parameters corresponding to empty overlaps.
#'
#' \code{emptyoverlaps} Matrix with two rows, giving the list pairs within the model for which no cases are observed in common.
#'   Each column gives the indices of a pair of lists, with the names of the lists in the column name.
#'
#' \code{poisspempty} the Poisson p-values of the empty overlaps.
#'
#' @examples
#' data(NewOrl)
#' modelfit(NewOrl,c(1,3), check=TRUE)
#'
#' @importFrom stats glm poisson
#'
#' @export

modelfit <- function(zdat, mX=NULL, check = T) {
  #--- perform check if asked for
  if (check) {
    zcheck = checkident(zdat, mX, verbose = T)
    if (zcheck$ierr > 0) return(zcheck$ierr)
  }
  zz = buildmodel(zdat, mX)
  eo = zz$emptyoverlaps
  nover= dim(eo)[2]
  fit = glm(zz$modelform, family=poisson, data=zz$datamatrix, x=TRUE)
  fit$aic = fit$aic + 2*nover
  #----------deal with empty overlaps if necessary
  probzero = vector("numeric",nover)
  if (nover > 0) {
    names(probzero) = dimnames(eo)[[2]]
    for (iover in (1:nover)) {
      lambda = sum(fit$coefficients[c(1,1 + eo[,iover])])
      probzero[iover] = exp( -exp(lambda))
    }
  }
  return( list (fit=fit, emptyoverlaps = eo, poisspempty = probzero))
}

#'Stepwise fit using Poisson pvalues.
#'
#'Starting with a model with main effects only, pairwise interactions are added one by one.
#'At each stage the interaction with the lowest p-value is added, provided that p-value is lower than \code{pthresh}, and provided that the resulting
#' model does not fail either of the tests in \code{\link{checkident}}.
#'
#' For each candidate interaction for possible addition to the model, the p-value is calculated as follows.
#' The total of cases occurring on both lists indexed by the interaction (regardless of
#' whether or not they are on any other lists) is calculated.
#' On the null hypothesis that the effect is not included in the model, this statistic has a Poisson
#' distribution whose mean depends on the parameters within the model.    The one-sided Poisson p-value of the observed statistic is calculated.
#'
#'@param zdat Data matrix with \eqn{t+1} columns. The first \eqn{t} columns, each corresponding to a particular list,
#' are 0s and 1s defining the capture histories
#' observed. The last column is the count of cases with that particular capture history.
#' List names A, B, ... are constructed if not supplied. Where a capture history is not explicitly listed,
#' it is assumed that it has observed count zero.
#'
#'@param pthresh this is the threshold below which the p-value of the newly added parameter needs to be in order to be included in the model.
#'  If \code{pthresh = 0} then the model with main effects only is returned.
#'
#'@return A list with components as follow
#'
#'\code{fit} Details of the fit of the specified model as output by \code{glm}.  The Akaike information criterion is adjusted to take account
#'   of the parameters corresponding to empty overlaps.
#'
#'\code{emptyoverlaps} Matrix with two rows, each column of which gives the list pairs within the model for which empty overlaps are observed.
#'
#' \code{poisspempty} the Poisson p-value of the empty overlaps.
#'
#'@examples
#'data(NewOrl)
#'stepwisefit(NewOrl, pthresh=0.001)
#'
#'@importFrom stats ppois predict
#'
#'
#'@export


stepwisefit<- function(zdat, pthresh=0.001) {
  m = dim(zdat)[2] - 1
  mX = NULL
  zfit = modelfit(zdat, mX=NULL, check=F)
  if (pthresh == 0) return(zfit)
  # ------  Set up a list of interactions to be included in the model with interactions
  ints = NULL
  for (i in (1:(m - 1))) for (j in ((i + 1):m)) {
    ints = cbind(ints,c(i, j))
  }
  nints = dim(ints)[2]
  nstar = rep(0,nints)
  for (j in (1:nints)) {
    nstar[j] = sum( zdat[,ints[1,j]]*zdat[,ints[2,j]]*zdat[,m+1])
  }
  #--------cycle through and find best interaction to add.
  intsincluded = rep(F, nints)
  for (jcycle in (1:nints)) {
    pred0 = predict(zfit$fit)
    dat0 = zfit$fit$data
    pvalnew = rep(1, nints)
    for (j in (1:nints)) {
      if (!intsincluded[j]) {
        #----- find the predicted Poisson intensity of this overlap
        pstar = sum(exp(pred0) * dat0[,ints[1,j]]*dat0[,ints[2,j]] )
        pvalnew[j] = min( ppois(nstar[j],pstar), ppois(nstar[j]-1, pstar, F) )
        #----- test the potential model and set the p-value to 1 if the model
        #---      would yield an infinite result or be non-identifiable
        ierr = checkident(zdat, cbind(mX, ints[,j]))
        if (ierr > 0) pvalnew[j] = 1
      }
    }
    #---------if none of them have p value below threshold, finish, otherwise incorporate and proceed
    pvmin = min(pvalnew)
    if (pvmin >= pthresh) return(zfit)
    jintmax = min((1:nints)[pvalnew == pvmin])
    mX1 = cbind(mX, ints[, jintmax])
    zf1 = modelfit(zdat, mX = mX1, check=F)
    zfit = zf1
    mX = mX1
    intsincluded[jintmax] = T
  }
  return(zfit)
}

#' Plot of simulation study
#'
#' This routine produces Figure 1 of Chan, Silverman and Vincent (2019).
#'
#'   Simulations are carried out for two different three-list models.
#'   In one model, the probabilities of capture are 0.01, 0.04 and 0.2 for the three lists respectively, while in the other the probability
#'   is 0.3 on all three lists.  In both cases, there are no interaction
#'   effects, so that captures on the lists occur independently of each other.
#'   The first model is chosen to be somewhat more typical of the sparse capture
#'   case, of the kind which often occurs in the human trafficking context,
#'   while the second is a more classical multiple systems estimate.
#'
#'   The probability of an individual having each possible capture history is first evaluated.
#'   Then these probabilities are multiplied by \code{Nsamp = 1000} and, for each simulation
#'   replicate, Poisson random values with expectations equal to these values are generated
#'   to give a full set of observed capture histories;
#'   together with the null capture history the expected number of counts
#'   (population size) is equal to \code{Nsamp}.
#' Inference was carried out both for the model with main effects only, and for the model with the addition of an interaction effect between the first two lists.
#' The reduction in deviance between the two models was determined.
#'
#' Checking for compliance with the conditions for existence and identifiability of the
#' estimates shows that a very small number of the simulations for the sparse model (two out of
#' ten thousand) fail the checks for existence even within the extended maximum likelihood context.
#' Detailed investigation shows that in neither of these cases is the dark figure itself not estimable;
#' although the parameters themselves cannot all be estimated, there is a maximum likelihood estimate
#' of the expected capture frequencies, and hence the deviance can still be calculated.
#'
#' The routine produces QQ-plots
#'   of the resulting deviance reductions against quantiles of the \eqn{\chi^2_1} distribution,
#'   for \code{nsim} simulation replications.
#'
#' @param nsim The number of simulation replications
#' @param Nsamp The expected value of the total population within each simulation
#' @param seed  The random number seed
#'
#' @return
#'
#'  An \code{nsim} \eqn{\times} 2 matrix giving the changes in deviance for each replication for each of the two models.
#'
#'
#'
#' @references
#'Chan, L., Silverman, B. W., and Vincent, K. (2019).
#'  Multiple systems estimation for Sparse Capture Data: Inferential Challenges when there are Non-Overlapping Lists. Available from \url{https://arxiv.org/abs/1902.05156}.
#'
#'@importFrom stats qchisq rpois ppoints
#'@importFrom graphics matplot abline
#'
#' @export
investigateAIC <- function(nsim=10000, Nsamp= 1000, seed = 1001) {
  #
  # ---- define functions to be called
  #
  simulate_deviances <- function(N, plist=0.1, nsim=100, seed= NULL ) {
    if (length(plist)==1) plist= rep(plist, 3)
    set.seed(seed)
    devdiffs = rep(NA, nsim)
    for (j in (1:nsim)) {
      zdat1 = simulate_capture(N, plist)
      devdiffs[j] = deviance_compare(zdat1)
    }
    return(devdiffs)
  }
  #
  #
  simulate_capture <- function(N=1000, plist=rep(0.1,3)) {
    zdata = rbind(diag(3), 1-diag(3), rep(1,3))
    meancounts = N* exp( zdata %*% log(plist) + (1-zdata)%*% log(1-plist))
    simcounts = rpois(7, meancounts)
    zdata = cbind(zdata,simcounts)
    dimnames(zdata)[[2]] = c("A", "B", "C", "count")
    return(zdata)
  }
  #
  #
  deviance_compare <- function(zdat) {
    if(sum(zdat[4:7,4])==0) return(0)
    devnull =  modelfit(zdat)$fit$deviance
    devalt =   modelfit(zdat, mX=c(1,2),check=F)$fit$deviance
    devdiff = devnull - devalt
    return(devdiff)
  }
  #
  # ----carry out the simulations and plot the results
  #
  zdev0 = simulate_deviances(Nsamp, plist=c(0.01,0.04,0.2),nsim, seed=seed )
  zdev1 = simulate_deviances(Nsamp, plist=0.3,nsim, seed=seed )
  zdev= cbind(sort(zdev0), sort(zdev1))
  zchisq = qchisq(ppoints(nsim), df=1)
  matplot(zchisq, zdev, type="l", xlim=c(0,8), ylim=c(0,8), xaxs="i", yaxs="i",
          xlab="Chi-Squared distribution", ylab="Simulated deviances", main="QQ plot")
  abline( 0,1, lty=3)
  return(zdev)
}



