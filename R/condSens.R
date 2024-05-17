condSens <- function(...) UseMethod("condSens")

#' @title Sensitivity analysis for conditional exposure effects.
#'
#' @description \code{condSens()} is the main function of condSens and
#' performs the sensitivity analysis for estimating the sensitivity range(s).
#'
#' @param data A data frame in which contains outcome, measured confounders, exposures. The order of variables (outcome, measured confounders, exposures) that make up the data is important.
#' @param outcome The name of variable for outcome.
#' @param fitmodel Model output fitted using one of "lm", "randomforest", and "gbm".
#' @param model What model was fitted (one of "lm", "randomforest", and "gbm").
#' @param k The number of measured confounders.
#' @param p The number of exposures.
#' @param bound The range of \eqn{(\phi_{1}, \ldots, \phi_{k})} and \eqn{(\rho_{1}, \ldots, \rho_{p})}. The order of inputs is c(\eqn{\phi_{1,max}, \ldots, \phi_{k,max}, } \eqn{\rho_{1,max}, \ldots, \rho_{p,max}, \phi_{1,min}, \ldots, \phi_{k,min}, \rho_{1,min}, \ldots, \rho_{p,min}}). See Examples.
#' @param delta.range The range of \eqn{(\delta_{min}, \delta_{max})}. See Examples.
#' @param delta.diff The increment of the sequence of \eqn{\delta}. Default: 0.1.
#' @param decimal.p Number of digits after the decimal point of numbers appearing in the result table.
#' @param report.result Whether to output the result table. Default: TRUE.
#' @param only.sig Whether to output results only for significant variables. Default: TRUE.
#' @param n.visual.delta The Number of \eqn{\delta} in the plot to be printed. This is only used when only.sig=FALSE. Default: 4.
#'
#' @return An object of class \code{condSens}. The object is a data.frame with the following components: 1. only.sig=FALSE
#' \item{label}{Exposures}
#' \item{model_output}{The coefficients of exposures}
#' \item{cond_min}{The minimum of conditional single- and joint-exposure effect}
#' \item{cond_max}{The maximum of conditional single- and joint-exposure effect}
#' \item{delta}{The values of \eqn{\delta}}
#' \item{joint}{Whether the value is a conditional joint-exposure effect}
#' \item{bs}{The coefficient of significant exposure}
#' \item{bias_min}{Minimum value of bias}
#' \item{bias_max}{Maximum value of bias}
#' \item{delta}{The values of \eqn{\delta}}
#' \item{cond_min}{The minimum of conditional single-exposure effect}
#' \item{cond_max}{The maximum of conditional single-exposure effect}
#'
#' To generate the plot of results for the \code{condSens}, use the \code{\link{SensPlot}} functions.
#'
#' @details See Jeong et al. (2024) for details.
#'
#' @examples
#' ## Import data in out paper and NA omit
#' library(readr)
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
#' ## Change bmi (binary) 0:< 30 1: >=30
#' data$bmi_cat3 <- ifelse(data$bmi_cat3 >= 3, 1, 0)
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
#' ## Fitting the linear model
#' fitmodel <- lm(TELOMEAN ~., data=data_r3)
#'
#' ## Sensitivity analysis results
#' ### Figure 6
#' ## L-U is uncorrelated, and X-U is correlated with ranges of (0.18, 95).
#' ## 1. All exposures
#' rst1 <- condSens(data=data_r3, outcome="TELOMEAN", fitmodel=fitmodel, model="lm",
#'                  k=k, p=p,
#'                  bound=c(rep(0, 11), rep(0.95, 18),
#'                          rep(0, 11), rep(0.18, 18)),
#'                  delta.range=c(-0.03, -0.01), delta.diff=0.01, decimal.p=3,
#'                  report.result=TRUE, only.sig=FALSE, n.visual.delta=3)
#'
#' ## 2. Only significant exposure
#' rst1_sig <- condSens(data=data_r3, outcome="TELOMEAN", fitmodel=fitmodel, model="lm",
#'                      k=k, p=p,
#'                      bound=c(rep(0, 11), rep(0.95, 18),
#'                              rep(0, 11), rep(0.18, 18)),
#'                      delta.range=c(-0.03, -0.01), delta.diff=0.01, decimal.p=3,
#'                      report.result=TRUE, only.sig=TRUE, n.visual.delta=3)
#'
#' @seealso
#'  \code{\link[SAMU]{SensPlot}}
#'
#' @references
#' Jeong B, Lee S, Ye S, Lee D, Lee W (2024):
#' Sensitivity analysis for effects of multiple exposures in the presence of unmeasured confounding
#' \emph{xxx}. DOI: xxx.
#'
#' @keywords methods
#'
#' @export
condSens <- function(data=NULL, outcome=NULL, fitmodel, model="lm",
                     k, p, bound,
                     delta.range, delta.diff=0.1, decimal.p=3,
                     report.result=TRUE, only.sig=TRUE, n.visual.delta=4){

  ##############################################################################
  if(!(model %in% c("lm", "randomforest", "gbm"))){
    stop("For other models except for lm, randomforest, and gbm, condSens() is currently undergoing development.")
  }

  ##############################################################################
  ## Rearrange data
  dataX <- data[,colnames(data) != outcome]
  dataY <- data[,colnames(data) == outcome]
  data <- data.frame("Y"=dataY, dataX)

  if (sum(round(apply(dataX, 2, mean), 5) != 0) > 0 | sum(round(apply(dataX, 2, sd), 5) != 1) > 0) {
    stop("Data for confounders and exposures must be scaled with mean 0 and variance 1!")
  }

  ## calculate correlation matrix
  corrmatrix <- cor(data[,-1])

  ##############################################################################
  if (model == "lm") {
    if (sum(is(fitmodel) %in% "lm") == 0) {
      stop("model is not linear model.")
    }

    ## h*(L,xs+1,X-s) - h*(L,xs,X-s)
    bs <- fitmodel$coef[-1] # k+p vector

  } else if (model == "randomforest") {
    if (sum(is(fitmodel)[2] %in% "randomForest") == 0) {
      stop("model is not random forest.")
    }
    if (is.null(data) | is.null(outcome) | is.null(k) | is.null(p)) {
      stop("'data', 'outcome', 'k', and 'p' must be specified for 'random forest' or 'gbm'.")
    }

    ## Make prediction (Orig)
    prediction <- mean(predict(fitmodel, data))

    ## h*(L,xs+1,X-s) - h*(L,xs,X-s)
    bs <- c()
    for (i in 1:p) {
      predict_data_temp <- data
      predict_data_temp[,1+k+i] <- predict_data_temp[,1+k+i] + 1

      ## Make predictions
      prediction_temp <- mean(predict(fitmodel, predict_data_temp))
      bs[i] <- prediction - prediction_temp # p vector
    }
    bs <- c(rep(0, k), bs) # k+p vector
    names(bs) <- colnames(dataX)

  } else if (model == "gbm") {
    if (sum(is(fitmodel) %in% "gbm") == 0) {
      stop("model is not generalized boosted model.")
    }
    if (is.null(data) | is.null(outcome) | is.null(k) | is.null(p)) {
      stop("'data', 'outcome', 'k', and 'p' must be specified for 'random forest' or 'gbm'.")
    }

    ## Make prediction (Orig)
    prediction <- mean(predict(fitmodel, data, n.trees=fitmodel$n.trees))

    ## h*(L,xs+1,X-s) - h*(L,xs,X-s)
    bs <- c()
    for (i in 1:p) {
      predict_data_temp <- data
      predict_data_temp[,1+k+i] <- predict_data_temp[,1+k+i] + 1

      ## Make predictions
      prediction_temp <- mean(predict(fitmodel, predict_data_temp, n.trees=fitmodel$n.trees))
      bs[i] <- prediction - prediction_temp # p vector
    }
    bs <- c(rep(0, k), bs) # k+p vector
    names(bs) <- colnames(dataX)

  }

  bs_result <- bs

  ##############################################################################
  ##
  invcor <- as.matrix(solve(corrmatrix)) # inverse
  val <- bound
  dir <- c(rep("<=", (k+p)), rep(">=", (k+p)))
  Amat <- as.matrix(rbind(diag(1,(k+p)), diag(1,(k+p))))

  myQ3 <- invcor

  ##
  singlebias <- apply(invcor, 1, function(x) {
    biascal_new(invcor=x, Amat=Amat, dir=dir, val=val, k=k, p=p, Q=myQ3)
  })

  ## Joint bias
  if (k < 0) {
    stop("'k' must be non-negative value.")
  } else if (k == 0) {
    forjointmyq3 <- myQ3[,1:p]
  } else {
    forjointmyq3 <- myQ3[,-c(1:k)]
  }

  # Max
  mycop_max <- cop(f=linfun(a=rowSums(forjointmyq3), name="Obj.function"),
                   lc=lincon(A=Amat, dir=dir, val=val, name=seq(1, (k+p)*2)),
                   qc=quadcon(Q=myQ3, d=0, dir="<=", val=1, use=TRUE),
                   max=TRUE)
  # Min
  mycop_min <- cop(f=linfun(a=rowSums(forjointmyq3), name="Obj.function"),
                   lc=lincon(A=Amat, dir=dir, val=val, name=seq(1, (k+p)*2)),
                   qc=quadcon(Q=myQ3, d=0, dir="<=", val=1, use=TRUE),
                   max=FALSE)

  res_max <- solvecop(op=mycop_max, solver="alabama", quiet=TRUE)
  res_min <- solvecop(op=mycop_min, solver="alabama", quiet=TRUE)
  maxbias <- validate(op=mycop_max, sol=res_max, quiet=TRUE)$obj.fun
  minbias <- validate(op=mycop_min, sol=res_min, quiet=TRUE)$obj.fun

  Joint <- c(minbias, maxbias)

  resultdf <- as.data.frame(t(cbind(singlebias, Joint)))

  ##############################################################################
  if (only.sig == TRUE) {
    ##
    if ((sum(is(fitmodel) %in% "lm") == 0)) {
      stop("If 'only.sig == TRUE', then model must be 'lm'.")
    }

    ##
    resultdf <- resultdf[rownames(resultdf) != "Joint", ] # remove joint effect
    new_dat <- data.frame(bs=bs_result, resultdf); colnames(new_dat) <- c("bs", "bias.min", "bias.max")

    ## Only results for exposure
    if (k==0) {
      new_dat <- new_dat
    } else {
      new_dat <- new_dat[-c(1:k),]
    }

    ## Only significant exposure
    is.sig <- (summary(fitmodel)$coefficients[,4] < 0.05)[-c(1:(1+k))]
    sig.exposure <- names(is.sig)[is.sig]
    sensresultonly <- new_dat[rownames(new_dat) %in% sig.exposure,]

    ## Plot
    plotdat <- list()
    for(iii in 1:length(sig.exposure)){

      plot.dat <- sensresultonly[rownames(sensresultonly) == sig.exposure[iii], ]

      ##
      sign.delta <- unique(sign(delta.range))

      if (sum(sign.delta) %in% c(0,1)) {
        if (delta.range[2] >= 0) {
          delta.range <- sort(delta.range, decreasing=FALSE)
        } else {
          delta.range <- sort(delta.range, decreasing=TRUE)
        }

        ## Find touch zero or Not
        intercept <- plot.dat$bs[1]
        slope1 <- -plot.dat$bias.min[1]
        slope2 <- -plot.dat$bias.max[1]

        y_value1 <- slope1 * delta.range + intercept
        y_value2 <- slope2 * delta.range + intercept

        if (length(unique(sign(y_value1))) != 1) {
          xofy0_1 <- -intercept / slope1
          label.zero_1 <- paste0("(", round(xofy0_1, 3), ", ", 0, ")")
        }
        if (length(unique(sign(y_value2))) != 1) {
          xofy0_2 <- -intercept / slope2
          label.zero_2 <- paste0("(", round(xofy0_2, 3), ", ", 0, ")")
        }

        ## Data for plot
        delta.comb <- unique(sort(c(0, delta.range[1], seq(delta.range[1], delta.range[2], by=delta.diff), delta.range[2])))
        plot.dat.comb <- data.frame(bs=rep(plot.dat$bs, length(delta.comb)),
                                    bias_min=rep(plot.dat$bias.min, length(delta.comb)),
                                    bias_max=rep(plot.dat$bias.max, length(delta.comb)),
                                    delta=delta.comb,
                                    cond_min=pmin(plot.dat$bs - delta.comb * plot.dat$bias.min, plot.dat$bs - delta.comb * plot.dat$bias.max),
                                    cond_max=pmax(plot.dat$bs - delta.comb * plot.dat$bias.min, plot.dat$bs - delta.comb * plot.dat$bias.max))

      } else if (sum(sign.delta) == -1) {
        if (delta.range[2] >= 0) {
          delta.range <- sort(delta.range, decreasing=FALSE)
        } else {
          delta.range <- sort(delta.range, decreasing=TRUE)
        }

        ## Find touch zero or Not
        intercept <- plot.dat$bs[1]
        slope1 <- -plot.dat$bias.min[1]
        slope2 <- -plot.dat$bias.max[1]

        y_value1 <- slope1 * delta.range + intercept
        y_value2 <- slope2 * delta.range + intercept

        if (length(unique(sign(y_value1))) != 1) {
          xofy0_1 <- -intercept / slope1
          label.zero_1 <- paste0("(", round(xofy0_1, 3), ", ", 0, ")")
        }
        if (length(unique(sign(y_value2))) != 1) {
          xofy0_2 <- -intercept / slope2
          label.zero_2 <- paste0("(", round(xofy0_2, 3), ", ", 0, ")")
        }

        ## Data for plot
        delta.comb <- unique(sort(c(0, delta.range[1], seq(delta.range[1], delta.range[2], by=-delta.diff), delta.range[2])))
        plot.dat.comb <- data.frame(bs=rep(plot.dat$bs, length(delta.comb)),
                                    bias_min=rep(plot.dat$bias.min, length(delta.comb)),
                                    bias_max=rep(plot.dat$bias.max, length(delta.comb)),
                                    delta=delta.comb,
                                    cond_min=pmin(plot.dat$bs - delta.comb * plot.dat$bias.min, plot.dat$bs - delta.comb * plot.dat$bias.max),
                                    cond_max=pmax(plot.dat$bs - delta.comb * plot.dat$bias.min, plot.dat$bs - delta.comb * plot.dat$bias.max))
      }

      ##
      plot.dat.comb <- plot.dat.comb[order(plot.dat.comb$delta),]; rownames(plot.dat.comb) <- NULL
      plotdat[[iii]] <- plot.dat.comb
    }

    names(plotdat) <- sig.exposure

    if(report.result == TRUE){
      temp.result <- plotdat
      temp.result <- lapply(temp.result, function(x) round(x, decimal.p))
      print(temp.result)
    }

    invisible(list(result=plotdat))

  } else {

    if (length(n.visual.delta) > 4) {
      stop("The maximum length of 'n.visual.delta' is allowed up to 4. Please enter a vector smaller than or equal to 4.")
    }

    if (k == 0) {
      coef <- bs_result
    } else {
      coef <- bs_result[-c(1:k)]
    }
    label <- c(names(coef), "Joint effect")

    jointeffcoef <- sum(coef)

    if (k==0) {
      resultdf <- resultdf
    } else {
      resultdf <- resultdf[-c(1:k),]
    }

    df.init <- data.frame(cbind(label, resultdf, c(coef, jointeffcoef)))
    colnames(df.init) <- c("label", "bias_low", "bias_upper", "coef")

    if (delta.range[2] >= 0) {
      delta.range <- sort(delta.range, decreasing=FALSE)
    } else {
      delta.range <- sort(delta.range, decreasing=TRUE)
    }

    delta <- round(seq(delta.range[1], delta.range[2], length.out=n.visual.delta), 4)
    for (i in 1:n.visual.delta) {
      tempc0 <- delta[i]

      df_temp <- df.init

      c_1 <- df_temp$coef - tempc0*df_temp$bias_upper
      c_2 <- df_temp$coef - tempc0*df_temp$bias_low

      df_temp$cond_min <- pmin(c_1, c_2)
      df_temp$cond_max <- pmax(c_1, c_2)

      df_temp$c0 <- rep(tempc0, dim(df_temp)[1])

      if (i == 1) {
        df <- df_temp
      } else {
        df <- rbind(df, df_temp)
      }
    }

    colnames(df) <- c("label", "bias_low", "bias_upper", "model_output", "cond_min", "cond_max", "delta")
    sensresult <- data.frame(df[,c("label", "model_output", "cond_min", "cond_max", "delta")])
    sensresult$joint <- ifelse(sensresult$label=="Joint effect", 1, 0)
    rownames(sensresult) <- NULL


    if(report.result == TRUE){
      temp.result <- sensresult
      temp.result[,c(2:4)] <- round(temp.result[,c(2:4)], decimal.p)
      print(temp.result)
    }

    invisible(list(result=sensresult))
  }
}
