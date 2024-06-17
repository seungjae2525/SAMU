#' @title Sensitivity analysis for marginal exposure effects.
#'
#' @description \code{MarginalSens()} is the function for estimating the sensitivity range(s) of marginal exposure effect.
#'
#' @param data A data frame in which contains outcome, measured confounders, exposures. The order of variables (outcome, measured confounders, exposures) that make up the data is important. See Examples.
#' @param fittedmodel Model output fitted using one of "lm", "randomforest", and "gbm".
#' @param k The number of measured confounders.
#' @param p The number of exposures.
#' @param delta The value of \eqn{\delta}.
#' @param bound The range of \eqn{(\phi_{1}, \ldots, \phi_{k})} and \eqn{(\rho_{1}, \ldots, \rho_{p})}. The order of inputs is (\eqn{\phi_{1,max}, \ldots, } \eqn{ \phi_{k,max},\rho_{1,max}, \ldots, \rho_{p,max}, \phi_{1,min}, \ldots, \phi_{k,min}, \rho_{1,min}, \ldots, \rho_{p,min}}). See Examples.
#' @param n.col Number of columns to display plot.
#' @param xlim x-axis range.
#' @param ylim y-axis range.
#' @param text.which Coordinates to display slopes.
#' @param text.size Size of text displaying slopes.
#' @param line.size Line size of slopes.
#'
#' @details See Section 3 in Jeong et al. (2024) for details.
#'
#' @examples
#' \dontrun{
#' ## Import data and NA omit
#' library (readr)
#' ul="https://raw.githubusercontent.com/lizzyagibson/SHARP.Mixtures.Workshop/master/Data/studypop.csv"
#' data <- as.data.frame(read_csv(url(ul))); data <- na.omit(data)
#'
#' ## Log transformation for outcome
#' data$TELOMEAN <- log(data$TELOMEAN)
#'
#' ## Change the order of variables:: 1 outcome - 11 confounders - 18 exposures
#' data <- data[,c("TELOMEAN",
#'                 # confounders
#'                 "bmi_cat3","edu_cat","race_cat","male","ln_lbxcot","age_cent",
#'                 "LBXWBCSI","LBXLYPCT","LBXMOPCT","LBXEOPCT","LBXBAPCT",
#'                 # exposures
#'                 "LBX074LA","LBX099LA","LBX138LA","LBX153LA","LBX170LA","LBX180LA",
#'                 "LBX187LA","LBX194LA","LBXPCBLA","LBXHXCLA","LBX118LA",
#'                 "LBXD03LA","LBXD05LA","LBXD07LA",
#'                 "LBXF03LA","LBXF04LA","LBXF05LA","LBXF08LA"
#' )]
#' colnames(data) <- c("TELOMEAN",
#'                     "bmi_cat3","edu_cat","race_cat","male","ln_lbxcot","age_cent",
#'                     "LBXWBCSI","LBXLYPCT","LBXMOPCT","LBXEOPCT","LBXBAPCT",
#'                     "PCB74","PCB99","PCB138","PCB153","PCB170","PCB180",
#'                     "PCB187","PCB194","PCB126","PCB169","PCB118",
#'                     "Dioxin1","Dioxin2","Dioxin3",
#'                     "Furan1","Furan2","Furan3","Furan4")
#'
#' ## Change bmi to binary // 0: < 30 1: >= 30
#' data$bmi_cat3 <- ifelse(data$bmi_cat3 >= 3, 1, 0)
#'
#' ## Change race to binary // 4 (Non-Hispanic White) vs others
#' data$race_cat <- ifelse(data$race_cat == 4, 1, 0)
#'
#' ## Change education to binary // after college (3 and 4) vs before college (1 and 2)
#' data$edu_cat <- ifelse(data$edu_cat %in% c(3,4), 1, 0)
#'
#' ## Log transformation of exposures
#' names.logtrans <- c("PCB74","PCB99","PCB138","PCB153","PCB170","PCB180",
#'                     "PCB187","PCB194","PCB126","PCB169","PCB118",
#'                     "Dioxin1","Dioxin2","Dioxin3",
#'                     "Furan1","Furan2","Furan3","Furan4")
#' data[,names.logtrans] <- log(data[,names.logtrans])
#'
#' ## Standardized variables
#' data_r3 <- data
#' data_r3[,-c(1)] <- scale(data_r3[,-c(1)])
#'
#' # Sensitivity analysis for conditional exposure effect
#' ## Working model: Linear regression (confounders adjusted model)
#' ## Number of confounders
#' k <- 11
#' ## Number of exposures
#' p <- 18
#'
#' # Sensitivity analysis for marginal exposure effect
#' ## Working model: Random forest (confounders adjusted model)
#' ## 10-fold cross-validation for tuning parameter (mtry)
#' set.seed(231111)
#' library(caret)
#' library(randomForest)
#' metric <- "RMSE"
#' control <- trainControl(method="cv", number=10, search="grid")
#' tunegrid <- expand.grid(.mtry=seq(1, 29, 2))
#' rf_gridsearch <- train(TELOMEAN~., data=data_r3, method="rf",
#'                        metric=metric, tuneGrid=tunegrid, trControl=control)
#'
#' ## model fitting
#' mtry <- as.numeric(rf_gridsearch$bestTune)
#' rfmodel <- randomForest(TELOMEAN~., data=data_r3, mtry=mtry)
#'
#' ### PDP-based sensitivity analysis results for marginal exposure effect
#' ### Supplementary Figure 22
#' ## L-U is uncorrelated, and X-U correlated with ranges of (0.18, 0.95).
#' mrst1 <- MarginalSens(data_r3, rfmodel, k=k, p=p, delta=-0.02,
#'                       bound=c(rep(0, 11), rep(0.95, 18),  # upper bounds
#'                               rep(0, 11), rep(0.18, 18)), # lower bounds
#'                       n.col=6, xlim=c(-4, 4), ylim=c(-0.2, 0.2),
#'                       text.which=c(-1.0, 4), text.size=2, line.size=0.8)
#' }
#'
#' @import ggplot2
#'
#' @references
#' Jeong B, Lee S, Ye S, Lee D, Lee W (2024):
#' Sensitivity analysis for effects of multiple exposures in the presence of unmeasured confounding
#' \emph{xxx}. DOI: xxx.
#'
#' @keywords Methods
#'
#' @export
MarginalSens <- function(data, fittedmodel, k, p, delta, bound, n.col, xlim, ylim, text.which, text.size, line.size){
  plotsave <- list()
  for (expindx in 1:p) {

    print(paste("Exposure index ---", expindx))
    # grid
    pdvalue <- seq(min(data[-c(1:(k+1))][,expindx]), max(data[-c(1:(k+1))][,expindx]), length.out=50)

    final_bias_band <- c()
    for (i in 1:length(pdvalue)) {
      invcor <- as.matrix(solve(cor(data[,-1])))
      predicttempdata <- data[,-1]
      predicttempdata[,(k+expindx)] <- rep(pdvalue[i], dim(predicttempdata)[1])
      pdvalue_i <- pdvalue[i]
      x_input <- data[,-1]
      if (is.null(fittedmodel$n.tree)) {
        # ex: linear regression, random forest
        if (delta > 0) {
          bias_lower <- mean(predict(fittedmodel, predicttempdata)) -
            delta*(biascal_new_forpdp(data=data, x=x_input, pdvalue_i=pdvalue_i, i=i, k=k, p=p, expindx=expindx, bound=bound)[2])
          bias_upper <- mean(predict(fittedmodel, predicttempdata)) -
            delta*(biascal_new_forpdp(data=data, x=x_input, pdvalue_i=pdvalue_i, i=i, k=k, p=p, expindx=expindx, bound=bound)[1])
        } else {
          bias_lower <- mean(predict(fittedmodel, predicttempdata)) -
            delta*(biascal_new_forpdp(data=data, x=x_input, pdvalue_i=pdvalue_i, i=i, k=k, p=p, expindx=expindx, bound=bound)[1])
          bias_upper <- mean(predict(fittedmodel, predicttempdata)) -
            delta*(biascal_new_forpdp(data=data, x=x_input, pdvalue_i=pdvalue_i, i=i, k=k, p=p, expindx=expindx, bound=bound)[2])
        }
        biasrange <- c(mean(predict(fittedmodel, predicttempdata)), bias_lower, bias_upper)
      } else {
        # ex: gbm
        if (delta > 0) {
          bias_lower <- mean(predict(fittedmodel, data.frame(predicttempdata), n.trees=fittedmodel$n.trees)) -
            delta*(biascal_new_forpdp(data=data, x=x_input, pdvalue_i=pdvalue_i, i=i, k=k, p=p, expindx=expindx, bound=bound)[2])
          bias_upper <- mean(predict(fittedmodel, data.frame(predicttempdata), n.trees=fittedmodel$n.trees)) -
            delta*(biascal_new_forpdp(data=data, x=x_input, pdvalue_i=pdvalue_i, i=i, k=k, p=p, expindx=expindx, bound=bound)[1])
        } else {
          bias_lower <- mean(predict(fittedmodel, data.frame(predicttempdata), n.trees=fittedmodel$n.trees)) -
            delta*(biascal_new_forpdp(data=data, x=x_input, pdvalue_i=pdvalue_i, i=i, k=k, p=p, expindx=expindx, bound=bound)[1])
          bias_upper <- mean(predict(fittedmodel, data.frame(predicttempdata), n.trees=fittedmodel$n.trees)) -
            delta*(biascal_new_forpdp(data=data, x=x_input, pdvalue_i=pdvalue_i, i=i, k=k, p=p, expindx=expindx, bound=bound)[2])
        }
        biasrange <- c(mean(predict(fittedmodel, data.frame(predicttempdata), n.trees=fittedmodel$n.trees)), bias_lower, bias_upper)
      }
      final_bias_band <- rbind(final_bias_band, biasrange)
    }
    biasmldata <- data.frame(final_bias_band)
    colnames(biasmldata) <- c("mean", "low", "upper")
    biasmldata$xlab <- pdvalue
    pd <- data.frame(cbind(biasmldata$xlab, biasmldata$mean))
    colnames(pd) <- c("exposure", "yhat")

    slopelow <- (biasmldata[dim(biasmldata)[1],]$upper - biasmldata[1,]$low) / (biasmldata$xlab[dim(biasmldata)[1]] - biasmldata$xlab[1])
    slopehigh <- (biasmldata[dim(biasmldata)[1],]$low - biasmldata[1,]$upper) / (biasmldata$xlab[dim(biasmldata)[1]] - biasmldata$xlab[1])

    if (slopelow < slopehigh) {
      textdata <- data.frame("sloperange"=paste("(",as.character(round(slopelow,3)),", ",as.character(round(slopehigh,3)),")",sep=""))
    } else {
      textdata <- data.frame("sloperange"=paste("(",as.character(round(slopehigh,3)),", ",as.character(round(slopelow,3)),")",sep=""))
    }

    if (missing(xlim)) {
      xlim <- c(min(data[-c(1:(k+1))][,expindx]), max(data[-c(1:(k+1))][,expindx]))
    }

    if (missing(ylim)) {
      ylim <- c(min(biasmldata$low), max(biasmldata$upper))
    }

    xx <- text.which[1]
    yy <- text.which[2]

    tplot <- ggplot(biasmldata, aes_string(x='xlab', y='mean')) +
      geom_ribbon(data=biasmldata, mapping=aes_string(ymax='upper', ymin='low'), fill="lightblue", alpha=0.6) +
      theme_forest() +
      ylab("yhat") +
      xlab(colnames(data)[expindx+k+1]) +
      geom_line(color="blue3", size=line.size) +
      geom_text(x=xx, y=yy, label=as.character(textdata[1,1]), size=text.size) +
      coord_cartesian(xlim=xlim, ylim=ylim) +
      ggplot2::theme(axis.text.y=element_text(colour="black", size=10),
                     axis.text.x=element_text(colour="black", size=10),
                     axis.title.x=element_text(colour="black", size=10),
                     axis.title.y=element_text(colour="black", size=10),
                     strip.text.x=element_text(colour="black", size=10),
                     strip.text.y=element_text(colour="black", size=10))
    plotsave[[expindx]] <- tplot
  }
  do.call(grid.arrange, c(plotsave, ncol=n.col))
  return(plotsave)
}

