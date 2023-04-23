#' @title Calculate 95\% CI for NNT based on separable effects
#'
#' @description To calculate the 95\% confidence intervals for Number Needed to Treat (NNT) based on the separable effects.
#'
#' @param dat The data set used for the final analysis.
#' @param Y The outcome with 3 levels. By default, \code{Y=0} for censoring, \code{Y=1} for an event of interest, and \code{Y=2} for competing events.
#' @param dTime The time to events, including censoring, an event of interest, and competing events.
#' @param cutTimes The discrete time intervals for the pooled logistic regression models.
#' @param eoiValue The value of \code{Y=1} for the event of interest by default.
#' @param crsValue The value of \code{Y=2} for the competing events by default.
#' @param A The indicator for exposure or treatment with 2 levels. By default, \code{A=1} for the treatment group and \code{A=0} for the control group.
#' @param X The minimal sufficient of confounders considered in the final analysis,
#' @param id The unique ID for each observed included in \code{dat}.
#' @param nboot The number of bootstrapped samples. By default, \code{nboot=50}.
#'
#' @return A matrix list contains three data.frames
#' \describe{
#'   \item{datRessult}{The estimated NNT and relevant measurements}
#'   \item{datCumulativeIncidece}{The estimated cumulative incidence}
#'   \item{datBoots}{The bootstrapped distributions of NNT and relevant measurements}
#' }
#'
#' @import survival
#' @import stats
#' @import parallel
#'
#' @references
#' \enumerate{
#'   \item Yang et al. (2023) Disentangling the inverse LDL-C-hemorrhagic stroke association in Chinese adults with hypertension: Findings from the Chinese Multi-provincial Cohort Study. Am J Epidemiol. (Revision)
#'   \item Stensrud el al. (2022) Separable effects for causal inference in the presence of competing events. J Am Stat Assoc. 117(537):175-83. <doi: https://doi.org/10.1080/01621459.2020.1765783>
#' }
#'
#' @export
#'
separableNNTCIs <- function(dat,
                            Y,eoiValue,crsValue,dTime,cutTimes,
                            A,X,id,nboot) {
  dat=dat; Y=Y;eoiValue=eoiValue;crsValue=crsValue
  dTime=dTime;cutTimes=cutTimes
  A=A;X=X;id=id;
  # Point estimate
  pe <- separableNNT(dat=dat, Y=Y,
                     eoiValue=eoiValue,
                     crsValue=crsValue,
                     dTime=dTime,
                     cutTimes=cutTimes,
                     A=A, X=X, id=id)
  n <- nrow(dat)

  # Parallel computing
  cl <- parallel::makeCluster(parallel::detectCores())
  parallel::clusterEvalQ(cl, "separableNNT")
  parallel::clusterExport(cl,list("dat","Y","dTime","cutTimes",
                                  "nboot","A","separableNNT","eoiValue","crsValue",
                                  "X"),
                          envir = environment())
  peBoots <- stats::setNames(parallel::parLapply(
    cl,1:nboot,
    fun=function(i) {
      n <- nrow(dat)
      idx <- sample(1:n,size=n,replace=TRUE)
      datBoot <- dat[idx,]
      datBoot$ID <- 1:nrow(datBoot)
      separableNNT(dat=datBoot,
                   Y=Y,
                   eoiValue=eoiValue,
                   crsValue=crsValue,
                   dTime=dTime,
                   cutTimes=cutTimes,
                   A=A, X=X, id="ID")$datResult
    }),1:nboot)
  parallel::stopCluster(cl)
  peBoots <- do.call("rbind",peBoots)


  peLCITE <- apply(
    peBoots[peBoots$Outcome %in% c("Event of interest") &
              peBoots$Effect %in% c("TE=Pr[Y(aY=1,aD=1)=1]-Pr[Y(aY=0,aD=0)=1]")
            ,4:ncol(peBoots)],2,
    FUN = function(i) quantile(i,probs=0.025))
  peUCITE <- apply(
    peBoots[peBoots$Outcome %in% c("Event of interest") &
              peBoots$Effect %in% c("TE=Pr[Y(aY=1,aD=1)=1]-Pr[Y(aY=0,aD=0)=1]")
            ,4:ncol(peBoots)],2,
    FUN = function(i) quantile(i,probs=0.975))
  peLCIDE <- apply(
    peBoots[peBoots$Outcome %in% c("Event of interest") &
              peBoots$Effect %in% c("DE=Pr[Y(aY=1,aD=1)=1]-Pr[Y(aY=0,aD=1)=1]")
            ,4:ncol(peBoots)],2,
    FUN = function(i) quantile(i,probs=0.025))
  peUCIDE <- apply(
    peBoots[peBoots$Outcome %in% c("Event of interest") &
              peBoots$Effect %in% c("DE=Pr[Y(aY=1,aD=1)=1]-Pr[Y(aY=0,aD=1)=1]")
            ,4:ncol(peBoots)],2,
    FUN = function(i) quantile(i,probs=0.975))
  peLCIIDE <- apply(
    peBoots[peBoots$Outcome %in% c("Event of interest") &
              peBoots$Effect %in% c("IDE=Pr[Y(aY=0,aD=1)=1]-Pr[Y(aY=0,aD=0)=1]")
            ,4:ncol(peBoots)],2,
    FUN = function(i) quantile(i,probs=0.025))
  peUCIIDE <- apply(
    peBoots[peBoots$Outcome %in% c("Event of interest") &
              peBoots$Effect %in% c("IDE=Pr[Y(aY=0,aD=1)=1]-Pr[Y(aY=0,aD=0)=1]")
            ,4:ncol(peBoots)],2,
    FUN = function(i) quantile(i,probs=0.975))


  ## total effect
  peResTE <- rbind(peLCITE,
                   pe$datResult[pe$datResult$Effect=="TE=Pr[Y(aY=1,aD=1)=1]-Pr[Y(aY=0,aD=0)=1]",4:ncol(pe$datResult)],
                   peUCITE)
  rownames(peResTE) <- c("LCI", "Est", "UCI")
  datResTE <- as.data.frame(t(peResTE))
  datResTE$Measures <- rownames(datResTE)
  datResTE$Outcome <- "Event of interest"
  datResTE$Effect <- "TE=Pr[Y(aY=1,aD=1)=1]-Pr[Y(aY=0,aD=0)=1]"
  datResTE$Method <- "g-formula"
  rownames(datResTE) <- NULL

  ## separable direct effect
  peResDE <- rbind(peLCIDE,
                   pe$datResult[pe$datResult$Effect=="DE=Pr[Y(aY=1,aD=1)=1]-Pr[Y(aY=0,aD=1)=1]",4:ncol(pe$datResult)],
                   peUCIDE)
  rownames(peResDE) <- c("LCI", "Est", "UCI")
  datResDE <- as.data.frame(t(peResDE))
  datResDE$Measures <- rownames(datResDE)
  datResDE$Outcome <- "Event of interest"
  datResDE$Effect <- "DE=Pr[Y(aY=1,aD=1)=1]-Pr[Y(aY=0,aD=1)=1]"
  datResDE$Method <- "g-formula"
  rownames(datResDE) <- NULL

  ## separable indirect effect
  peResIDE <- rbind(peLCIIDE,
                    pe$datResult[pe$datResult$Effect=="IDE=Pr[Y(aY=0,aD=1)=1]-Pr[Y(aY=0,aD=0)=1]",4:ncol(pe$datResult)],
                    peUCIIDE)
  rownames(peResIDE) <- c("LCI", "Est", "UCI")
  datResIDE <- as.data.frame(t(peResIDE))
  datResIDE$Measures <- rownames(datResIDE)
  datResIDE$Outcome <- "Event of interest"
  datResIDE$Effect <- "IDE=Pr[Y(aY=0,aD=1)=1]-Pr[Y(aY=0,aD=0)=1]"
  datResIDE$Method <- "g-formula"
  rownames(datResIDE) <- NULL

  datRes <- rbind(datResTE,datResDE,datResIDE)
  datRes$Est <- round(datRes$Est, 4)
  datRes$LCI <- round(datRes$LCI, 4)
  datRes$UCI <- round(datRes$UCI, 4)
  peBoots$NNT <- ceiling(peBoots$NNT)
  return(list(datResult=datRes[,c("Outcome","Effect","Measures","Method","LCI","Est","UCI")],
              datCumulativeIncidence=pe$datCumulativeIncidence,
              datBoots=peBoots))
}
