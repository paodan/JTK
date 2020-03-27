
#' Identify data with rhythms using JTK algorithm
#'
#' This is a non-parametric algorithm for detecting rhythmic components
#' in genome-scale datasets.
#' @param data expression matrix, each column represents a sample at
#' certain time, each row represents a gene/protein
#' @param timepoints total time points
#' @param reps replicates per time point
#' @param normal normal approximation
#' @param alt alternative distributions
#' @param periods rhythms of how many time points (i.e. between 5 and 7 time
#' points per cycle if periods = 5:7)
#' @param interval the number of time units (hours) between time points
#' @param ampci calculating amplitude confidence if TRUE, default is FALSE
#' @param conf confidence level
#' @return a data frame showing information of genes/proteins rhythms.
#' "BH.Q" is adjusted p value using BH method, "ADJ.P" is the rhythms
#' significance p value, "PER" is periods (in time unit), "LAG" is the lag of
#' rhythms, "AMP" is the amplitudes.
#' @export
#' @examples
#' time = seq(1, 72, by = 3)
#' set.seed(123)
#' expr1 = 2 + 3 * cos(time/(6) - 4) + rnorm(length(time), 0, 0.2)
#' expr2 = 3 + 2 * cos(time/(3) - 3) + rnorm(length(time), 0, 0.2)
#' expr3 = 3 + 2 * cos(time/(6) - 3) + rnorm(length(time), 0, 2)
#' plot(expr3, col = "black", ylim = c(-1, 8))
#' lines(expr3, col = "black")
#' points(expr2, col = "blue")
#' lines(expr2, col = "blue")
#' points(expr1, col= "red")
#' lines(expr1, col = "red")
#'
#' data = as.data.frame(rbind(expr1, expr2, expr3))
#'
#' JTK(data, timepoints = 24, periods = 2:24, interval = 3)
JTK = function(data, timepoints, reps=1, normal=FALSE, alt=FALSE,
               periods, interval, ampci=FALSE, conf=0.8){

  res = jtkdist(timepoints = timepoints, reps = reps,
                normal = normal, alt = alt, res = list())

  res = jtk.init(periods = periods, interval = interval, res)

  cat("JTK analysis started on",date(),"\n")
  flush.console()

  res1 <- apply(data, 1, function(z) {
    res_tmp = jtkx(z, res = res)
    c(res_tmp$JTK.ADJP, res_tmp$JTK.PERIOD, res_tmp$JTK.LAG, res_tmp$JTK.AMP)
  })
  res1 <- as.data.frame(t(res1))
  bhq <- p.adjust(unlist(res1[,1]),"BH")
  res1 <- cbind(bhq,res1)
  colnames(res1) <- c("BH.Q","ADJ.P","PER","LAG","AMP")
  results <- cbind(res1, data)
  results <- results[order(res1$ADJ.P,-res1$AMP),]
  return(results)
}

#' Fit a cosine curve for rhythmic data
#'
#' @param expr numeric vector, expression of a gene/protein
#' @param Time numeric vector, all time points (time unit)
#' @param n number of best fitted plots using different periods.
#' @param periodHours numeric vector, fit cosine curve (time unit)
#' @return a ggplot object of cosine curves with different periods.
#' @import ggplot2
#' @export
#' @examples
#' time = seq(1, 72, by = 3)
#' set.seed(123)
#' expr1 = 2 + 3 * cos(time/(6) - 4) + rnorm(length(time), 0, 0.2)
#' expr2 = 3 + 2 * cos(time/(3) - 3) + rnorm(length(time), 0, 0.2)
#' expr3 = 3 + 2 * cos(time/(6) - 3) + rnorm(length(time), 0, 2)
#' plot(expr3, col = "black", ylim = c(-1, 8))
#' lines(expr3, col = "black")
#' points(expr2, col = "blue")
#' lines(expr2, col = "blue")
#' points(expr1, col= "red")
#' lines(expr1, col = "red")
#'
#' data = as.data.frame(rbind(expr1, expr2, expr3))
#'
#' result = JTK(data, timepoints = 24, periods = 2:24, interval = 3)
#' plotJTK(as.numeric(result[1, 6:ncol(result)]), time)
#' plotJTK(as.numeric(result[1, 6:ncol(result)]), time, 1, 42)
#' plotJTK(as.numeric(result[2, 6:ncol(result)]), time)
#' plotJTK(as.numeric(result[3, 6:ncol(result)]), time)
plotJTK = function(expr, Time, n = 6,
                   periodHours = seq(1, 60, by = 0.5)){
  fitData = data.frame(expr = expr, Time = Time)
  re = vapply(setNames(periodHours, periodHours), function(mi){
    formula = expr ~ cos(2*pi*Time/mi)
    fit.lm <- lm(formula, data = fitData)
    fit <- fitted(fit.lm)
    sum(fit.lm$residuals^2)
  }, FUN.VALUE = 1)

  timeMany = seq(min(Time), max(Time), by = 0.5)
  n = min(length(periodHours), n)
  predData = lapply(setNames(as.numeric(names(sort(re)))[seq_len(n)],
                             names(sort(re))[seq_len(n)]),
                    function(mi){
                      formula = expr ~ cos(2*pi*Time/mi)
                      fit.lm <- lm(formula, data = fitData)
                      # find predictions for original time series
                      pred <- predict(fit.lm, newdata=data.frame(Time = timeMany),
                                      interval = "confidence")
                      colnames(pred)[1] = "expr"
                      return(data.frame(Time = timeMany, pred, mi = mi))
                    })
  predData = do.call("rbind", predData)

  ggplot(predData, aes(Time, expr))+
    geom_point(data = fitData)+
    geom_line()+
    facet_wrap(~mi, nrow = 3)+
    geom_ribbon(mapping = aes(ymin = lwr, ymax = upr),
                fill = "blue", alpha = 0.3)
}