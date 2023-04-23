#' @title Calculate the NNT based on separable effects
#'
#' @description To calculate the Number Needed to Treat (NNT) based on the separable effects.
#'
#' @param dat The data set used for the final analysis
#' @param Y The outcome with 3 levels. By default, \code{Y=0} for censoring, \code{Y=1} for an event of interest, and \code{Y=2} for competing events
#' @param dTime The time to events, including censoring, an event of interest, and competing events
#' @param cutTimes The discrete time intervals for the pooled logistic regression models
#' @param eoiValue The value of \code{Y=1} for the event of interest by default
#' @param crsValue The value of \code{Y=2} for the competing events by default
#' @param A The indicator for exposure or treatment with 2 levels. By default, \code{A=1} for the treatment group and \code{A=0} for the control group
#' @param X The minimal sufficient of confounders considered in the final analysis
#' @param id The unique ID for each observed included in \code{dat}
#'
#' @return A matrix list contains two data.frames
#' \describe{
#'   \item{datRessult}{The estimated NNT and relevant measurements}
#'   \item{datCumulativeIncidece}{The estimated cumulative incidence}
#'   \item{datBoots}{The bootstrapped distributions of NNT and relevant measurements}
#' }
#'
#' @import survival
#' @import stats
#'
#' @references
#' \enumerate{
#'   \item Yang et al. (2023) Disentangling the inverse LDL-C-hemorrhagic stroke association in Chinese adults with hypertension: Findings from the Chinese Multi-provincial Cohort Study. Am J Epidemiol. (Revision)
#'   \item Stensrud el al. (2022) Separable effects for causal inference in the presence of competing events. J Am Stat Assoc. 117(537):175-83. <doi: https://doi.org/10.1080/01621459.2020.1765783>
#' }
#'
#' @export
#'
separableNNT <- function(dat,         ## wide-format data set
                         Y,           ## outcome with 3 levels: 0,1,2
                         dTime,       ## survival time
                         cutTimes,    ## time interval for follow-up
                         eoiValue=1,  ## value for event of interest
                         crsValue=2,  ## value for competing events
                         A,           ## exposure with 2 levels: 0,1
                         X,           ## covariates
                         id           ## id for each subject
) {

  varList <- names(dat)
  if(sum(!c(Y,dTime,A,X,id) %in% varList)>0) {
    stop("Please check the variable names in the original dataset.")
  }

  ## For event of interest
  eValue <- c(eoiValue,crsValue)
  datRes <- datCIF <- data.frame()

  datGFormula <- dat
  datGFormula$Yall <- datGFormula[[Y]] %in% eValue
  datGFormula$Ycensor <- !datGFormula[[Y]] %in% eValue
  datGFormula$dtime_old <- datGFormula[[dTime]]

  # Wide format data -> Long format data
  datGFormula$Tstart <- -0.01
  datGFormulaYall <- survival::survSplit(data=datGFormula,
                                         cut=cutTimes,
                                         start="Tstart",
                                         end=dTime,
                                         event="Yall")
  datGFormulaYcensor <- survival::survSplit(data=datGFormula,
                                            cut=cutTimes,
                                            start="Tstart",
                                            end=dTime,
                                            event="Ycensor")
  datGFormulaYall$Ycensor <- datGFormulaYcensor$Ycensor

  calcCIF <- function(
    dat,
    id=id,
    dTime=dTime,
    timepts=cutTimes,
    competing=FALSE){
    # Create the matrix for saving CIF at each time intervals
    cumulativeIncidence <- matrix(NA,
                                  ncol = length(unique(dat[[id]])),
                                  nrow = length(timepts))
    # Insert the event probabilities at the first time interval:
    cumulativeIncidence[1,] <-
      dat[dat[[dTime]]==0,]$hazardP *
      dat[dat[[dTime]]==0,]$hazardO
    # Create a matrix compatible with 'cumulativeIncidence' with survival probabilities at each time
    survivalProb <- t(aggregate(s~dat[[id]],
                                data = dat,
                                FUN = cumprod)$s) #split the long data into subsets per subject
    for(i in 2:length(timepts)){
      subInputDataP <- dat[dat[[dTime]]==(i-1),]$hazardP     #OBS: dtime starts at 0
      subInputDataO <- (1-dat[dat[[dTime]]==(i-1),]$hazardO) #OBS: dtime starts at 0
      if(!competing) {
        # OBS: survivalProb entry i is time point i-1
        cumulativeIncidence[i,] <-
          subInputDataP * subInputDataO *  survivalProb[(i-1),]
      } else {
        cumulativeIncidence[i,] <-
          (1-subInputDataO) * survivalProb[(i-1),]
      }
    }
    meanCumulativeIncidence <-
      rowMeans(apply(cumulativeIncidence, MARGIN = 2,cumsum))
    return(meanCumulativeIncidence)
  }

  # Event of interest
  datGFormulaYall$Yeoi <- datGFormulaYall$Yall==1 & datGFormulaYall[[Y]]==eoiValue
  datGFormulaYall$Ycrs <- datGFormulaYall$Yall==1 & datGFormulaYall[[Y]]==crsValue
  # Temporal order (Censoring, CompetingRisks, EventOfInterest)
  # datGFormulaYall$Yeoi[datGFormulaYall$Ycensor==1] <- NA
  # datGFormulaYall$Ycrs[datGFormulaYall$Ycensor==1] <- NA
  # datGFormulaYall$Yeoi[datGFormulaYall$Ycrs==1] <- NA

  # sample size
  n <- length(unique(datGFormulaYall[[id]])); n

  # Create 'baseline' - data collected at visit 0
  baseline <- datGFormulaYall[datGFormulaYall[[dTime]]==0,];nrow(baseline)

  # datGFormulaYall <- datGFormulaYall[datGFormulaYall[[dTime]]<length(cutTimes),]

  # Create powers of time to allow time dependent baseline hazard
  datGFormulaYall$dtime2 <- datGFormulaYall[[dTime]] * datGFormulaYall[[dTime]]
  datGFormulaYall$dtime3 <- datGFormulaYall$dtime2 * datGFormulaYall[[dTime]]
  datGFormulaYall$dtime4 <- datGFormulaYall$dtime3 * datGFormulaYall[[dTime]]

  # Create interactions with time and treatment
  datGFormulaYall$A1 <- datGFormulaYall[[dTime]] * datGFormulaYall[[A]]
  datGFormulaYall$A2 <- datGFormulaYall$dtime2 * datGFormulaYall[[A]]
  datGFormulaYall$A3 <- datGFormulaYall$dtime3 * datGFormulaYall[[A]]
  datGFormulaYall$A4 <- datGFormulaYall$dtime4 * datGFormulaYall[[A]]

  # Create time-dependent effect used in the hazard of competing events
  datGFormulaYall$OA  <- datGFormulaYall[[A]]
  datGFormulaYall$OA1 <- datGFormulaYall$A1
  datGFormulaYall$OA2 <- datGFormulaYall$A2
  datGFormulaYall$OA3 <- datGFormulaYall$A3
  datGFormulaYall$OA4 <- datGFormulaYall$A4

  # ----
  # Conditional pooled logistic regression models used for parametric g-formula
  ## Outcome of interest
  fitP <- as.formula(
    paste("Yeoi~",A,"+A1+A2+",dTime,"+dtime2+dtime3+",
          paste(X,collapse="+"),sep=""));fitP
  plrFitP <- glm(fitP,family=binomial(),data=datGFormulaYall)
  ## Competing events
  fitO <- as.formula(
    paste("Ycrs~OA+OA1+OA2+",dTime,"+dtime2+dtime3+",
          paste(X,collapse="+"),sep=""));fitO
  plrFitO <- glm(fitO,family=binomial(),data=datGFormulaYall)

  # Create simulated data where everyone is treated
  # -- For the treated group
  treated <- baseline[rep(1:n,each=length(cutTimes)),]
  treated[[dTime]] <- rep(cutTimes,n)
  treated[[A]] <- TRUE
  treated$OA <- TRUE

  ##
  treated$dtime2 <- treated[[dTime]] * treated[[dTime]]
  treated$dtime3 <- treated$dtime2 * treated[[dTime]]
  treated$dtime4 <- treated$dtime3 * treated[[dTime]]
  ##
  treated$A1 <- treated[[dTime]] * treated[[A]]
  treated$A2 <- treated$dtime2 * treated[[A]]
  treated$A3 <- treated$dtime3 * treated[[A]]
  treated$A4 <- treated$dtime4 * treated[[A]]
  ##
  treated$OA1 <- treated[[dTime]] * treated$OA
  treated$OA2 <- treated$dtime2 * treated$OA
  treated$OA3 <- treated$dtime3 * treated$OA
  treated$OA4 <- treated$dtime4 * treated$OA
  ## Estimate conditional cause-specific discrete hazards
  treated$hazardP <- predict(plrFitP,newdata=treated,type="response")
  treated$hazardO <- predict(plrFitO,newdata=treated,type="response")
  treated$s <- with(treated, (1-hazardP)*(1-hazardO))
  sum(treated$hazardO) < 0

  # -- For the placebo group
  placebo <- baseline[rep(1:n,each=length(cutTimes)),]
  placebo[[dTime]] <- rep(cutTimes,n)
  placebo[[A]] <- FALSE
  placebo$OA <- FALSE
  ##
  placebo$dtime2 <- placebo[[dTime]] * placebo[[dTime]]
  placebo$dtime3 <- placebo$dtime2 * placebo[[dTime]]
  placebo$dtime4 <- placebo$dtime3 * placebo[[dTime]]
  ##
  placebo$A1 <- placebo[[dTime]] * placebo[[A]]
  placebo$A2 <- placebo$dtime2 * placebo[[A]]
  placebo$A3 <- placebo$dtime3 * placebo[[A]]
  placebo$A4 <- placebo$dtime4 * placebo[[A]]
  ##
  placebo$OA1 <- placebo[[dTime]] * placebo$OA
  placebo$OA2 <- placebo$dtime2 * placebo$OA
  placebo$OA3 <- placebo$dtime3 * placebo$OA
  placebo$OA4 <- placebo$dtime4 * placebo$OA
  ## Estimate conditional cause-specific discrete hazards
  placebo$hazardP <- predict(plrFitP,newdata=placebo,type="response")
  placebo$hazardO <- predict(plrFitO,newdata=placebo,type="response")
  placebo$s <- with(placebo, (1-hazardP)*(1-hazardO))

  ## Ay=0, Ad=1
  # -- For the treated group
  treatedAy <- baseline[rep(1:n,each=length(cutTimes)),]
  treatedAy[[dTime]] <- rep(cutTimes,n)
  treatedAy[[A]] <- FALSE
  treatedAy$OA <- TRUE

  ##
  treatedAy$dtime2 <- treatedAy[[dTime]] * treatedAy[[dTime]]
  treatedAy$dtime3 <- treatedAy$dtime2 * treatedAy[[dTime]]
  treatedAy$dtime4 <- treatedAy$dtime3 * treatedAy[[dTime]]
  ##
  treatedAy$A1 <- treatedAy[[dTime]] * treatedAy[[A]]
  treatedAy$A2 <- treatedAy$dtime2 * treatedAy[[A]]
  treatedAy$A3 <- treatedAy$dtime3 * treatedAy[[A]]
  treatedAy$A4 <- treatedAy$dtime4 * treatedAy[[A]]
  ##
  treatedAy$OA1 <- treatedAy[[dTime]] * treatedAy$OA
  treatedAy$OA2 <- treatedAy$dtime2 * treatedAy$OA
  treatedAy$OA3 <- treatedAy$dtime3 * treatedAy$OA
  treatedAy$OA4 <- treatedAy$dtime4 * treatedAy$OA
  ## Estimate conditional cause-specific discrete hazards
  treatedAy$hazardP <- predict(plrFitP,newdata=treatedAy,type="response")
  treatedAy$hazardO <- predict(plrFitO,newdata=treatedAy,type="response")
  treatedAy$s <- with(treatedAy, (1-hazardP)*(1-hazardO))
  sum(treated$hazardO) < 0

  CIFtreatedGFormula <- calcCIF(dat=treated,
                                id=id,dTime=dTime,
                                timepts=cutTimes,
                                competing=FALSE)
  CIFtreatAyGFormula <- calcCIF(dat=treatedAy,
                                id=id,dTime=dTime,
                                timepts=cutTimes,
                                competing=FALSE)
  CIFplaceboGFormula <- calcCIF(dat=placebo,
                                id=id,dTime=dTime,
                                timepts=cutTimes,
                                competing=FALSE)

  # CIFtreatedGFormula[36];CIFtreatedGFormula[40]
  # CIFtreatAyGFormula[36];CIFtreatAyGFormula[40]
  # CIFplaceboGFormula[36];CIFplaceboGFormula[40]

  treatedAy1Ad1RMST <- length(cutTimes)-sum(CIFtreatedGFormula);treatedAy1Ad1RMST
  placeboAy0Ad0RMST <- length(cutTimes)-sum(CIFplaceboGFormula);placeboAy0Ad0RMST
  treatedAy0Ad1RMST <- length(cutTimes)-sum(CIFtreatAyGFormula);treatedAy0Ad1RMST

  ## Total effect
  TERR <- CIFtreatedGFormula[length(cutTimes)]/CIFplaceboGFormula[length(cutTimes)];TERR
  TERD <- CIFtreatedGFormula[length(cutTimes)]-CIFplaceboGFormula[length(cutTimes)];TERD
  TERMSTR <- treatedAy1Ad1RMST/placeboAy0Ad0RMST; TERMSTR
  TERMSTD <- treatedAy1Ad1RMST-placeboAy0Ad0RMST; TERMSTD
  TENNT <- placeboAy0Ad0RMST/TERMSTD; TENNT

  ## Separable direct effect
  DERR <- CIFtreatedGFormula[length(cutTimes)]/CIFtreatAyGFormula[length(cutTimes)];DERR
  DERD <- CIFtreatedGFormula[length(cutTimes)]-CIFtreatAyGFormula[length(cutTimes)];DERD
  DERMSTR <- treatedAy1Ad1RMST/treatedAy0Ad1RMST; DERMSTR
  DERMSTD <- treatedAy1Ad1RMST-treatedAy0Ad1RMST; DERMSTD
  DENNT <- placeboAy0Ad0RMST/DERMSTD; DENNT

  ## Separable indirect effect
  IDERR <- CIFtreatAyGFormula[length(cutTimes)]/CIFplaceboGFormula[length(cutTimes)];IDERR
  IDERD <- CIFtreatAyGFormula[length(cutTimes)]-CIFplaceboGFormula[length(cutTimes)];IDERD
  IDERMSTR <- treatedAy0Ad1RMST/placeboAy0Ad0RMST;IDERMSTR
  IDERMSTD <- treatedAy0Ad1RMST-placeboAy0Ad0RMST;IDERMSTD
  IDENNT <- placeboAy0Ad0RMST/IDERMSTD;IDENNT

  datGFormulaRes <- rbind(
    data.frame(Outcome="Event of interest",
               Method="g-formula",
               Effect="TE=Pr[Y(aY=1,aD=1)=1]-Pr[Y(aY=0,aD=0)=1]",
               RR=TERR, RD=TERD,
               treatedRMST=treatedAy1Ad1RMST,
               placeboRMST=placeboAy0Ad0RMST,
               RMSTR=TERMSTR, RMSTD=TERMSTD, NNT=TENNT),
    data.frame(Outcome="Event of interest",
               Method="g-formula",
               Effect="DE=Pr[Y(aY=1,aD=1)=1]-Pr[Y(aY=0,aD=1)=1]",
               RR=DERR, RD=DERD,
               treatedRMST=treatedAy1Ad1RMST,
               placeboRMST=treatedAy0Ad1RMST,
               RMSTR=DERMSTR, RMSTD=DERMSTD, NNT=DENNT),
    data.frame(Outcome="Event of interest",
               Method="g-formula",
               Effect="IDE=Pr[Y(aY=0,aD=1)=1]-Pr[Y(aY=0,aD=0)=1]",
               RR=IDERR, RD=IDERD,
               treatedRMST=treatedAy0Ad1RMST,
               placeboRMST=placeboAy0Ad0RMST,
               RMSTR=IDERMSTR, RMSTD=IDERMSTD, NNT=IDENNT))
  datGFormulaCIFRes <- rbind(
    data.frame(Outcome="Event of interest",
               Method="g-formula",
               Effect="Pr[Y(aY=1,aD=1)=1]",
               Group="aY=1,aD=1",
               Time=c(0,cutTimes+1),
               CIF=c(0,CIFtreatedGFormula)),
    data.frame(Outcome="Event of interest",
               Method="g-formula",
               Effect="Pr[Y(aY=0,aD=1)=1]",
               Group="aY=0,aD=1",
               Time=c(0,cutTimes+1),
               CIF=c(0,CIFtreatAyGFormula)),
    data.frame(Outcome="Event of interest",
               Method="g-formula",
               Effect="Pr[Y(aY=0,aD=0)=1]",
               Group="aY=0,aD=0",
               Time=c(0,cutTimes+1),
               CIF=c(0,CIFplaceboGFormula)))
  datGFormulaRes$RR <- round(datGFormulaRes$RR,digits=4)
  datGFormulaRes$RD <- round(datGFormulaRes$RD,digits=4)
  datGFormulaRes$treatedRMST <- round(datGFormulaRes$treatedRMST,digits=4)
  datGFormulaRes$placeboRMST <- round(datGFormulaRes$placeboRMST,digits=4)
  datGFormulaRes$RMSTR <- round(datGFormulaRes$RMSTR,digits=4)
  datGFormulaRes$RMSTD <- round(datGFormulaRes$RMSTD,digits=4)
  datGFormulaRes$NNT <- ceiling(datGFormulaRes$NNT)
  datGFormulaCIFRes$CIF <- round(datGFormulaCIFRes$CIF,digits=4)
  return(list(datResult=datGFormulaRes[,c("Outcome","Effect","Method","RR","RD",
                                          "treatedRMST","placeboRMST","RMSTR","RMSTD","NNT")],
              datCumulativeIncidence=datGFormulaCIFRes))
}
