#' @title Sensitivity analysis for conditional exposure effects based on an alternative assumption.
#'
#' @description \code{condSens.alt()} is an alternative sensitivity analysis method using assumption for the conditional distribution of U | (L, X) rather than for the joint distribution of (U, L, X).
#'
#' @param data A data frame in which contains outcome, measured confounders, exposures. The order of variables (outcome, measured confounders, exposures) that make up the data is important. See Examples.
#' @param outcome The name of outcome variable.
#' @param fitmodel Model output fitted using one of "lm", "randomforest", and "gbm".
#' @param model What model was fitted (one of "lm", "randomforest", and "gbm").
#' @param k The number of measured confounders.
#' @param p The number of exposures.
#' @param exposure.interest Vector of names for exposures of interest.
#' @param bound The range of \eqn{(\rho_{1}, \ldots, \rho_{p})}. The order of inputs is (\eqn{\rho_{1,max}, \ldots, \rho_{p,max}, } \eqn{\rho_{1,min}, \ldots, } \eqn{ \rho_{p,min}}). See Examples.
#' @param delta.range The range of \eqn{(\delta_{min}, \delta_{max})}. See Examples.
#' @param delta.diff The increment of the sequence of \eqn{\delta}. Default: 0.1.
#' @param decimal.p Number of digits after the decimal point of numbers appearing in the result table.
#' @param report.result Whether to output the result table. Default: TRUE.
#' @param only.sig Whether to output results only for significant exposure(s). Default: TRUE.
#' @param n.visual.delta The number of \eqn{\delta} in the plot to be printed. This is only used when \code{only.sig=FALSE}. Default: 3.
#'
#' @return The object is a data.frame with the following components:
#' \itemize{
#'   \item When \code{only.sig=FALSE}:
#'   \itemize{
#'     \item \code{label}: Exposures
#'     \item \code{model_output}: The coefficients of exposures
#'     \item \code{cond_min}: The minimum of conditional single- and joint-exposure effect
#'     \item \code{cond_max}: The maximum of conditional single- and joint-exposure effect
#'     \item \code{delta}: The values of \eqn{\delta}
#'     \item \code{joint}: Whether the value is a conditional joint-exposure effect
#'   }
#'   \item When \code{only.sig=TRUE}:
#'   \itemize{
#'     \item \code{bs}: The coefficient of significant exposure
#'     \item \code{bias_min}: Minimum value of bias
#'     \item \code{bias_max}: Maximum value of bias
#'     \item \code{delta}: The values of \eqn{\delta}
#'     \item \code{cond_min}: The minimum of conditional single-exposure effect
#'     \item \code{cond_max}: The maximum of conditional single-exposure effect
#'   }
#' }
#'
#' To generate the plot of results for the \code{condSens.alt}, use the \code{\link{SensPlot}} function.
#'
#' @details See Section 2.4.3 in Jeong et al. (2024) for details.
#'
#' @examples
#' # if (!require("MASS", quietly=TRUE)) install.packages("MASS")
#' library(MASS)
#'
#' ## Our simple running example
#' set.seed(231113)
#' n.sample <- 10000
#' delta <- -0.2
#' rho1 <- 0.1
#' rho2 <- 0.2
#' r11 <- 0.6
#' dat <- as.data.frame(mvrnorm(n=n.sample, c(0, 0, 0),
#'                              matrix(c(1,rho1,rho2,rho1,1,r11,rho2,r11,1),3,3),
#'                              empirical=TRUE))
#' colnames(dat) <- c("U", "X1", "X2")
#' dat$Y <- 0.5 * dat$X1 + 0.7 * dat$X2 + delta * dat$U + rnorm(n.sample)
#'
#' dat_noU <- dat[,c("Y", "X1", "X2")]
#' fit.model <- lm(Y ~ X1 + X2, data=dat_noU)
#'
#' rho1_star <- seq(0.05, 0.15, 0.01)
#' rho2_star <- seq(0.15, 0.25, 0.01)
#' rho_grid <- expand.grid(rho1_star, rho2_star)
#'
#' rho_star_f <- function(i, j){
#'   t(c(i, j)) %*% solve(matrix(c(1, 0.6, 0.6, 1), 2, 2))
#' }
#' rho_star <- t(sapply(1:nrow(rho_grid), function(x) {
#'   rho_star_f(rho_grid[x,1], rho_grid[x,2])
#' }))
#'
#' rho1_star_min <- min(rho_star[,1])
#' rho1_star_max <- max(rho_star[,1])
#' rho2_star_min <- min(rho_star[,2])
#' rho2_star_max <- max(rho_star[,2])
#'
#' ## All exposures
#' result.alt <- condSens.alt(data=dat_noU, outcome="Y", fitmodel=fit.model, model="lm",
#'                            k=0, p=2, exposure.interest=c("X1", "X2"),
#'                            bound=c(rho1_star_max, rho2_star_max,  # upper bounds
#'                                    rho1_star_min, rho2_star_min), # lower bounds
#'                            delta.range=c(-0.3, -0.1), delta.diff=0.1, decimal.p=3,
#'                            report.result=TRUE, only.sig=FALSE, n.visual.delta=3)
#' SensPlot(result.alt$result, myxlim=c(0.4, 1.3))
#'
#' ## Only significant exposures
#' result.alt.sig <- condSens.alt(data=dat_noU, outcome="Y", fitmodel=fit.model, model="lm",
#'                                k=0, p=2, exposure.interest=c("X1", "X2"),
#'                                bound=c(rho1_star_max, rho2_star_max,  # upper bounds
#'                                        rho1_star_min, rho2_star_min), # lower bounds
#'                                delta.range=c(-0.3, -0.1), delta.diff=0.1, decimal.p=3,
#'                                report.result=TRUE, only.sig=TRUE, n.visual.delta=3)
#' SensPlot(result.alt.sig$result$X1)
#' SensPlot(result.alt.sig$result$X2)
#'
#' @seealso
#'  \code{\link[SAMU]{SensPlot}}
#'
#' @references
#' Jeong B, Lee S, Ye S, Lee D, Lee W (2024):
#' Sensitivity analysis for effects of multiple exposures in the presence of unmeasured confounding
#' \emph{xxx}. DOI: xxx.
#'
#' @keywords Methods
#'
#' @export
condSens.alt <- function(data=NULL, outcome=NULL, fitmodel, model="lm",
                         k, p, exposure.interest, bound,
                         delta.range, delta.diff=0.1, decimal.p=3,
                         report.result=TRUE, only.sig=TRUE, n.visual.delta=3){
  ##############################################################################
  if(!(model %in% c("lm", "randomforest", "gbm"))){
    stop("For other models except for 'lm', 'randomforest', and 'gbm', condSens() is currently undergoing development.")
  }

  if (is.null(data) | is.null(outcome) | is.null(k) | is.null(p)) {
    stop("'data', 'outcome', 'k', and 'p' must be specified.")
  }

  length_e <- length(exposure.interest)
  length_bound <- length(bound)

  if (p != length_e) {
    stop("The input of 'p' or 'bound' is incorrect. See Examples for how to input 'p' or 'bound'.")
  }

  if (length_e != length_bound/2) {
    warning("The sensitivity intervals of the conditional joint-exposure effect is not the sensitivity intervals of effect when all exposures increase by one unit, but the sensitivity intervals of effect when exposures of interest increase by one unit.")
  }

  upper_bound <- bound[1:length_e]
  lower_bound <- bound[(length_e+1):length_bound]

  if (sum(lower_bound > upper_bound) >= 1) {
    stop("Some 'lower_bound' values are larger than 'upper_bound' values. Please check again.")
  }

  ##############################################################################
  ## Rearrange data
  dataX <- data[,colnames(data) != outcome]
  dataY <- data[,colnames(data) == outcome]
  data <- data.frame("Y"=dataY, dataX)

  if (sum(round(apply(dataX, 2, mean), 5) != 0) > 0 | sum(round(apply(dataX, 2, sd), 5) != 1) > 0) {
    stop("Data for confounders and exposures must be scaled with mean 0 and variance 1!")
  }

  if (identical(intersect(colnames(dataX), exposure.interest), exposure.interest) == FALSE) {
    stop("Some of 'exposure.interest' is specified incorrectly!")
  }

  ##############################################################################
  if (model == "lm") {
    if (sum(is(fitmodel) %in% "lm") == 0) {
      stop("'model' is not linear model.")
    }

    ## h*(L,xs+1,X-s) - h*(L,xs,X-s)
    bs <- fitmodel$coef[-1] # k+p vector

  } else if (model == "randomforest") {
    if (sum(is(fitmodel)[2] %in% "randomForest") == 0) {
      stop("'model' is not random forest.")
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
      stop("'model' is not generalized boosted model.")
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

  ## Only results for exposure
  if (k==0) {
    bs_result <- bs
  } else {
    bs_result <- bs[-c(1:k)]
  }

  ##############################################################################
  resultdf <- matrix(c(lower_bound, sum(lower_bound), upper_bound, sum(upper_bound)),
                     nrow=length_e+1, ncol=2, byrow=FALSE)

  rownames(resultdf) <- c(exposure.interest, "Joint")

  ##############################################################################
  if (only.sig == TRUE) {
    ##
    if ((sum(is(fitmodel) %in% "lm") == 0)) {
      stop("If 'only.sig == TRUE', then 'model' must be 'lm'.")
    }

    ##
    resultdf <- resultdf[rownames(resultdf) != "Joint", ] # remove joint effect
    new_dat <- data.frame(bs=bs_result, resultdf); colnames(new_dat) <- c("bs", "bias.min", "bias.max")

    ## Only significant exposure
    if (k==0) {
      is.sig <- (summary(fitmodel)$coefficients[,4] < 0.05)[-c(1)]
    } else {
      is.sig <- (summary(fitmodel)$coefficients[,4] < 0.05)[-c(1:(1+k))]
    }
    sig.exposure <- names(is.sig)[is.sig]

    if (sum(exposure.interest %in% sig.exposure) == 0) {
      stop("'exposure.interest' does not contain any significant exposures.")
    }

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
      stop("The maximum length of 'n.visual.delta' is allowed up to 4. Please enter a value smaller than or equal to 4.")
    }

    coef <- bs_result
    label <- c(names(coef), "Joint effect")

    jointeffcoef <- sum(coef)

    df.init <- data.frame(label, resultdf, c(coef, jointeffcoef))
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

