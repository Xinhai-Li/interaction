#' Fractional-power interaction regression (FPIR): y ~  x1 + x2 + x3 + x1*x2*x3
#'
#' @description This function estimates the optimum values of exponents for x1, x2, and x3 in FPIR,
#'  return model parameters such as regression coefficients, sum of squares of the interaction term, R square, AIC.
#'
#' @author Xinhai Li (Xinhai_li_edu@126.com)
#'
#' @param y the dependent variable in regression
#' @param x1 the first independent variables in regression y ~  x1 + x2 + x3 + x1*x2*x3
#' @param x2 the second independent variable in regression y ~  x1 + x2 + x3 + x1*x2*x3
#' @param x3 the third independent variable in regression y ~  x1 + x2 + x3 + x1*x2*x3
#'
#' @return return a data.frame with columns the optimum values of exponents for x1, x2, and x3 in FPIR,
#' R square value of the whole model, proportion of the variance explained by the interaction term,
#' coefficients of all terms (the interception, coefficients of x1, x2, x3, and the interaction term).
#'
#' @examples
#'
#'  x1 = runif(100); x2 = runif(100); x3 = runif(100); y = x1 * x2^2.5 * x3^-12 + 0.1*x3; y = jitter(y) # generate a dataset
#'  out = FPIR1threeway(y, x1, x2, x3)
#'  out[[2]] # first round exponents
#'  out[[4]] # second round exponents (tuned results)
#'  out[[5]] # final model coefficients
#'  out[[6]] # adjusted R square
#'
#' @export


# Fractional-power interaction regression (FPIR), y ~  x1 + x2 + x3 + x1*x2*x3
# one dependent variable y, three explanatory variables x1, x2 and x3

FPIR1threeway = function(y, x1, x2, x3) {
  out = list()
  # The first round comparison
  # FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
  Exp1 = Exp2 = Exp3 = c(-50, -40, -30, -20, -14, -9, -5, -2, -1, -0.5, 0, 0.5, 1, 2, 5, 9, 14, 20, 30, 40, 50)
  Exps.2 = expand.grid(Exp1, Exp2); Exps.2 = cbind(ID=1:nrow(Exps.2), Exps.2)
  Exps.3 = expand.grid(seq(1, nrow(Exps.2)), Exp3)
  Exps = merge(Exps.2, Exps.3, by.x = 'ID', by.y = 'Var1')
  names(Exps) = c('ID','e1', 'e2','e3') #9261

  Len = nrow(Exps) #9261
  out[[1]] = cbind(Exps, Rsq_T=NA, Rsq_I=NA,  AIC=NA)
  for(i in 1:Len){
    inter = x1^Exps[i,2] * x2^Exps[i,3] * x3^Exps[i,4]
    fit = lm(y ~ x1 + x2 + x3 + inter)
    reg = summary(fit)
    out[[1]]$Rsq_T[i]  = reg$adj.r.squared
    out[[1]]$Rsq_I[i]  = anova(fit)[3,2]/sum(anova(fit)[1:4,2])
    out[[1]]$AIC[i]    =AIC(fit)
    print(paste("The first round:", round(i/Len*100,2), "% finished", sep=' '))
  }
  # FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
  (out[[2]] = out[[1]][out[[1]]$Rsq_T == max(out[[1]]$Rsq_T),]) # the best model
  if (nrow(out[[2]] > 1))  out[[2]] <- out[[2]][1,]

  # The second round comparison
  # SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
  e1 = out[[2]]$e1; e2 = out[[2]]$e2; e3 = out[[2]]$e3
  LIST <- c(-50, -40, -30, -20, -14, -9, -5, -2, -1, -0.5, 0, 0.5, 1, 2, 5, 9, 14, 20, 30, 40, 50)

  if (abs(e1 < 50)) {
    order_1 = LIST==e1;
    order_1_L = order_1[-1];              order_1_L = c(order_1_L, FALSE); Exp1_L = LIST[order_1_L]
    order_1_U = order_1[-length(order_1)];order_1_U = c(FALSE, order_1_U); Exp1_U = LIST[order_1_U]
    Exp1_s = seq(Exp1_L+(e1-Exp1_L)/2, e1+(Exp1_U-e1)/2, length=10)
  }
  if (e1 == -50) {Exp1_s = seq(-55, -45, length=10)}
  if (e1 ==  50) {Exp1_s = seq(55, 45,   length=10)}

  if (abs(e2 < 50)) {
    order_2 = LIST==e2;
    order_2_L = order_2[-1];              order_2_L = c(order_2_L, FALSE); Exp2_L = LIST[order_2_L]
    order_2_U = order_2[-length(order_2)];order_2_U = c(FALSE, order_2_U); Exp2_U = LIST[order_2_U]
    Exp2_s = seq(Exp2_L+(e2-Exp2_L)/2, e2+(Exp2_U-e2)/2, length=10)
  }
  if (e2 == -50) {Exp2_s = seq(-55, -45, length=10)}
  if (e2 ==  50) {Exp2_s = seq(55, 45,   length=10)}

  if (abs(e3 < 50)) {
    order_3 = LIST==e3;
    order_3_L = order_3[-1];              order_3_L = c(order_3_L, FALSE); Exp3_L = LIST[order_3_L]
    order_3_U = order_3[-length(order_3)];order_3_U = c(FALSE, order_3_U); Exp3_U = LIST[order_3_U]
    Exp3_s = seq(Exp3_L+(e3-Exp3_L)/2, e3+(Exp3_U-e3)/2, length=10)
  }
  if (e3 == -50) {Exp3_s = seq(-55, -45, length=10)}
  if (e3 ==  50) {Exp3_s = seq(55, 45,   length=10)}

  exps.2 = expand.grid(Exp1_s, Exp2_s); exps.2 = cbind(ID=1:nrow(exps.2), exps.2)
  exps.3 = expand.grid(seq(1, nrow(exps.2)), Exp3_s)
  exps = merge(exps.2, exps.3, by.x = 'ID', by.y = 'Var1')
  names(exps) = c('ID','e1', 'e2','e3')

  Len = nrow(exps) #
  out[[3]] = cbind(exps, Rsq_T=NA, Rsq_I=NA,  AIC=NA)
  for(i in 1:Len){
    inter = x1^exps[i,2] * x2^exps[i,3] * x3^exps[i,4]
    fit = lm(y ~ x1 + x2 + x3 + inter)
    reg = summary(fit)
    out[[3]]$Rsq_T[i]  = reg$adj.r.squared
    out[[3]]$Rsq_I[i]  = anova(fit)[3,2]/sum(anova(fit)[1:4,2])
    out[[3]]$AIC[i]    =AIC(fit)
    print(paste("The second round:", i/Len*100, "% finished", sep=' '))
  }
  # SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
  (out[[4]] = out[[3]][out[[3]]$Rsq_T == max(out[[3]]$Rsq_T),]) # the best model

  inter = x1^out[[4]]$e1 * x2^out[[4]]$e2 * x3^out[[4]]$e3
  fit = lm(y ~ x1 + x2 + x3 + inter)
  sum.fit = summary(fit)
  out[[5]] = fit$coefficients
  out[[6]] = sum.fit$adj.r.squared
  print(out[[4]][2:4])
  return(out)
}

