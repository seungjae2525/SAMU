#' @title Plot for sensitivity analysis results of conditional exposure effect
#'
#' @description Plot for sensitivity analysis results of conditional exposure effect.
#' When \code{only.sig=FALSE} in \code{condSens()} or \code{condSens.alt()}, the resulting plot for all exposures appears.
#' When \code{only.sig=TRUE} in \code{condSens()} or \code{condSens.alt()}, the resulting plot(s) for only significant exposure(s) appear.
#'
#' @param condSens.result Object from \code{condSens()} or \code{condSens.alt()}.
#' @param title Plot title name. This is only used when \code{only.sig=FALSE} in \code{condSens()} or \code{condSens.alt()}. Default: "Sensitivity analysis".
#' @param myxlim The minimum and maximum values of x-axis. This is only used when \code{only.sig=FALSE} in \code{condSens()} or \code{condSens.alt()}.
#' @param exposure.n The variable name(s) of significant exposure(s). This is only used when \code{only.sig=TRUE} in \code{condSens()} or \code{condSens.alt()}. Default: NULL.
#' @param alpha.range A opacity of sensitivity interval. Values of \code{alpha.range} range from 0 to 1, with lower values corresponding to more transparent colors. This is only used when \code{only.sig=TRUE} in \code{condSens()} or \code{condSens.alt()}. Default: 0.4.
#' @param alpha.ci A opacity of confidence interval. Values of \code{alpha.ci} range from 0 to 1, with lower values corresponding to more transparent colors. This is only used when \code{only.sig=TRUE} in \code{condSens()} or \code{condSens.alt()}. Default: 0.9.
#' @param ytickdiff Distance between y-axis tick. This is only used when \code{only.sig=TRUE} in \code{condSens()} or \code{condSens.alt()}. Default: 0.01.
#' @param point.size Point size of the the lower and upper bound of the sensitivity interval. This is only used when \code{only.sig=TRUE} in \code{condSens()} or \code{condSens.alt()}. Default: 1.4.
#' @param h.width Width of horizon line which represents the effect size. This is only used when \code{only.sig=TRUE} in \code{condSens()} or \code{condSens.alt()}. Default: 1.
#' @param axis.title.size Size of x and y axis title. This is only used when \code{only.sig=TRUE} in \code{condSens()} or \code{condSens.alt()}. Default: 15.
#' @param axis.text.size Size of x and y axis text. This is only used when \code{only.sig=TRUE} in \code{condSens()} or \code{condSens.alt()}. Default: 12.
#' @param label.size label text size. This is only used when \code{only.sig=TRUE} in \code{condSens()} or \code{condSens.alt()}. Default: 4.
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
#' ## Fitting the linear model
#' fitmodel <- lm(TELOMEAN ~., data=data_r3)
#'
#' ## Sensitivity analysis results
#' ### For Figure 6 (a)
#' ## L-U is uncorrelated, and X-U is correlated with ranges of (0.18, 95).
#' ## 1. All exposures
#' rst1 <- condSens(data=data_r3, outcome="TELOMEAN", fitmodel=fitmodel, model="lm",
#'                  k=k, p=p,
#'                  bound=c(rep(0, 11), rep(0.95, 18),  # upper bounds
#'                          rep(0, 11), rep(0.18, 18)), # lower bounds
#'                  delta.range=c(-0.03, -0.01), delta.diff=0.01, decimal.p=3,
#'                  report.result=TRUE, only.sig=FALSE, n.visual.delta=3)
#' SensPlot(condSens.result=rst1$result, myxlim=c(-0.2, 0.2),
#'          title="Sensitivity analysis\nCase 1")
#'
#' ## 2. Only significant exposure
#' ### For Figure 6 (b)
#' rst1_sig <- condSens(data=data_r3, outcome="TELOMEAN", fitmodel=fitmodel, model="lm",
#'                      k=k, p=p,
#'                      bound=c(rep(0, 11), rep(0.95, 18),  # upper bounds
#'                              rep(0, 11), rep(0.18, 18)), # lower bounds
#'                      delta.range=c(-0.03, -0.01), delta.diff=0.01, decimal.p=3,
#'                      report.result=TRUE, only.sig=TRUE, n.visual.delta=3)
#' SensPlot(condSens.result=rst1_sig$result$Furan1, exposure.n="Furan1",
#'          alpha.range=0.4, alpha.ci=0.9, ytickdiff=0.01, point.size=1.4, h.width=1,
#'          axis.title.size=13, axis.text.size=10, label.size=3.5)
#'
#' @seealso
#'  \code{\link[SAMU]{condSens}} and \code{\link[SAMU]{condSens.alt}}
#'
#' @import ggplot2
#'
#' @references
#' Jeong B, Lee S, Ye S, Lee D, Lee W (2024):
#' Sensitivity analysis for effects of multiple exposures in the presence of unmeasured confounding
#' \emph{xxx}. DOI: xxx.
#'
#' @keywords Plot
#'
#' @export
SensPlot <- function(condSens.result,
                     title="Sensitivity analysis", myxlim=NULL,
                     exposure.n=NULL, alpha.range=0.4, alpha.ci=0.9,
                     ytickdiff=0.01, point.size=1.4, h.width=1,
                     axis.title.size=15, axis.text.size=12, label.size=4){
  if (colnames(condSens.result)[6] == "cond_max") {
    #### For results when only.sig=TRUE ####
    plot.dat <- condSens.result
    delta.range <- c(min(plot.dat$delta), max(plot.dat$delta))
    ## Find touch zero or Not
    intercept <- plot.dat$bs[1]
    slope1 <- -plot.dat$bias_min[1]
    slope2 <- -plot.dat$bias_max[1]

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

    ##
    yscale <- 10^(floor(log10(ytickdiff)))
    ylabel <- seq(floor(min(plot.dat$cond_min)/yscale)*yscale,
                  ceiling(max(plot.dat$cond_max)/yscale)*yscale, ytickdiff)
    ylabel <- round(ylabel, 5)


    if (plot.dat$delta[nrow(plot.dat)] == 0) { # all negative delta
      xlabel <- -plot.dat$delta
      plot.dat$delta <- -plot.dat$delta
      ##
      g <- ggplot(plot.dat) +
        geom_ribbon(data=plot.dat, mapping=aes(x=.data$delta, ymin=.data$cond_min, ymax=.data$cond_max),
                    alpha=alpha.range, inherit.aes=F, fill="red") +
        geom_line(aes(x=delta, y=intercept), colour="red", linewidth=h.width) +
        { if(min(plot.dat$cond_min) < 0 | plot.dat$bs[1] < 0) {
          geom_line(aes(x=delta, y=0), colour="black", linewidth=h.width, linetype="dashed", alpha=0.5)
        } } +
        xlab(expression(bold(delta))) + ylab("Effect")  +
        theme_bw() +
        scale_x_continuous(breaks=xlabel, labels=-xlabel, expand=c(0.0015,0.0015)) +
        scale_y_continuous(breaks=ylabel, labels=ylabel) +
        theme(title=element_text(face="bold", size=axis.title.size-1),
              axis.title.x=element_text(face="bold", size=axis.title.size),
              axis.title.y=element_text(face="bold", size=axis.title.size, angle=90),
              axis.text=element_text(face="bold", size=axis.text.size),
              panel.grid.minor.x=element_blank(),
              panel.grid.minor.y=element_blank())

      if (length(unique(sign(y_value1))) != 1 & length(unique(sign(y_value2))) != 1) {
        g <- g +
          annotate("point", x=-xofy0_1, y=0, size=point.size, colour="blue") +
          annotate("label", label=label.zero_1, x=-xofy0_1 - 0.006, y=0 + 0.003, size=label.size, colour="black",
                   hjust=0, vjust="inward") +
          annotate("point", x=-xofy0_2, y=0, size=point.size, colour="blue") +
          annotate("label", label=label.zero_2, x=-xofy0_2 - 0.006, y=0 + 0.003, size=label.size, colour="black",
                   hjust=0, vjust="inward")
      } else if (length(unique(sign(y_value1))) != 1) {
        g <- g +
          annotate("point", x=-xofy0_1, y=0, size=point.size, colour="blue") +
          annotate("label", label=label.zero_1, x=-xofy0_1 - 0.006, y=0 + 0.003, size=label.size, colour="black",
                   hjust=0, vjust="inward")
      } else if (length(unique(sign(y_value2))) != 1) {
        g <- g +
          annotate("point", x=-xofy0_2, y=0, size=point.size, colour="blue") +
          annotate("label", label=label.zero_2, x=-xofy0_2 - 0.006, y=0 + 0.003, size=label.size, colour="black",
                   hjust=0, vjust="inward")
      } else {
        g <- g
      }
    } else {
      xlabel <- plot.dat$delta
      ##
      g <- ggplot(plot.dat) +
        geom_ribbon(data=plot.dat, mapping=aes(x=.data$delta, ymin=.data$cond_min, ymax=.data$cond_max),
                    alpha=alpha.range, inherit.aes=F, fill="red") +
        geom_line(aes(x=delta, y=intercept), colour="red", linewidth=h.width) +
        { if(min(plot.dat$cond_min) < 0 | plot.dat$bs[1] < 0) {
          geom_line(aes(x=delta, y=0), colour="black", linewidth=h.width, linetype="dashed", alpha=0.5)
        } } +
        xlab(expression(bold(delta))) + ylab("Effect")  +
        theme_bw() +
        scale_x_continuous(breaks=xlabel, labels=xlabel, expand=c(0.0015,0.0015)) +
        scale_y_continuous(breaks=ylabel, labels=ylabel) +
        theme(title=element_text(face="bold", size=axis.title.size-1),
              axis.title.x=element_text(face="bold", size=axis.title.size),
              axis.title.y=element_text(face="bold", size=axis.title.size, angle=90),
              axis.text=element_text(face="bold", size=axis.text.size),
              panel.grid.minor.x=element_blank(),
              panel.grid.minor.y=element_blank())

      if (length(unique(sign(y_value1))) != 1 & length(unique(sign(y_value2))) != 1) {
        g <- g +
          annotate("point", x=xofy0_1, y=0, size=point.size, colour="blue") +
          annotate("label", label=label.zero_1, x=xofy0_1 - 0.006, y=0 + 0.003, size=label.size, colour="black",
                   hjust=0, vjust="inward") +
          annotate("point", x=xofy0_2, y=0, size=point.size, colour="blue") +
          annotate("label", label=label.zero_2, x=xofy0_2 - 0.006, y=0 + 0.003, size=label.size, colour="black",
                   hjust=0, vjust="inward")
      } else if (length(unique(sign(y_value1))) != 1) {
        g <- g +
          annotate("point", x=xofy0_1, y=0, size=point.size, colour="blue") +
          annotate("label", label=label.zero_1, x=xofy0_1 - 0.006, y=0 + 0.003, size=label.size, colour="black",
                   hjust=0, vjust="inward")
      } else if (length(unique(sign(y_value2))) != 1) {
        g <- g +
          annotate("point", x=xofy0_2, y=0, size=point.size, colour="blue") +
          annotate("label", label=label.zero_2, x=xofy0_2 - 0.006, y=0 + 0.003, size=label.size, colour="black",
                   hjust=0, vjust="inward")
      } else {
        g <- g
      }
    }

    ## Title
    if (is.null(exposure.n)) {
      g <- g
    } else {
      g <- g +
        ggtitle(paste0("Sensitivity analysis [", exposure.n,"]"))
    }

    suppressWarnings(print(g))

    return(list(plot=g))

  } else if (colnames(condSens.result)[6] == "joint") {
    #### For results when only.sig=FALSE #####
    c0upper <- condSens.result$delta[1]
    delta <- unique(condSens.result$delta)

    if (length(delta) == 1) {
      myColors <- c("#01579B")
    } else if (length(delta) == 2) {
      if (c0upper > 0) {
        myColors <- c("#01579B","#81D4FA")
      } else {
        myColors <- c("#81D4FA","#01579B")
      }
    } else if (length(delta) == 3) {
      if (c0upper > 0) {
        myColors <- c("#01579B","#29B6F6","#81D4FA")
      } else {
        myColors <- c("#81D4FA","#29B6F6","#01579B")
      }
    } else if (length(delta) == 4) {
      if (c0upper > 0) {
        myColors <- c("#01579B","#0288D1","#29B6F6","#81D4FA")
      } else {
        myColors <- c("#81D4FA","#29B6F6","#0288D1","#01579B")
      }
    } else {
      stop("Warnings!")
    }

    names(myColors) <- rev(levels(factor(condSens.result$delta)))
    trev <- levels(factor(condSens.result$delta))
    if (c0upper > 0) {
      condSens.result$delta <- factor(condSens.result$delta, levels=rev(trev))
    } else {
      condSens.result$delta <- factor(condSens.result$delta, levels=trev)
    }

    break.levels <- rev(levels(factor(condSens.result$delta)))

    condSens.result$label <- factor(condSens.result$label, levels=rev(unique(condSens.result$label)))

    if (is.null(myxlim)) {
      myxlim <- c(min(condSens.result$cond_min), max(condSens.result$cond_max))
    }

    suppressWarnings({
      fp <- ggplot(condSens.result) +
        geom_pointrange(data=condSens.result,
                        mapping=aes(x=.data$model_output, y=.data$label,
                                    xmin=.data$cond_min, xmax=.data$cond_max, color=.data$delta),
                        size=0.4, position=position_dodge(width=0.7)) +
        coord_cartesian(xlim=myxlim) +
        geom_vline(xintercept=0, linetype="solid", size=0.4, colour="black") +
        theme_forest() +
        ggtitle(title) +
        xlab("Effect") +
        ylab(" ") +
        guides(shape=FALSE) +
        scale_colour_manual(name=expression(bold(italic(delta)~"(delta)")), values=myColors,
                            breaks=break.levels)
    })

    suppressWarnings(print(fp))

    return(list(plot=fp))
  } else {
    stop("'condSens.result' is incorrect!")
  }
}
