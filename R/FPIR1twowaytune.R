#' Fractional-power interaction regression (FPIR) for regression y ~ x1 + x2 + x1^M *x2^N
#'
#' @description This function estimates the tuned values of exponents (M and N) in FPIR,
#'  return model parameters such as regression coefficients, sum of squares of the interaction term, R square, AIC.
#'
#'
#' @author Xinhai Li (Xinhai_li_edu@126.com)
#'
#' @param y the dependent variable in regression
#' @param x1 one of the independent variables in regression
#' @param x2 another independent variable in regression
#' @param Exp1 one of the independent variables in regression
#' @param Exp2 another independent variable in regression
#'
#' @return return a data.frame with columns of M, N, R square value of the whole model,
#' proportion of the variance explained by the interaction term, coefficients of all terms
#' (the interception, x1, x2, and the interaction).
#'
#' @examples
#'
#'  x1 = runif(100); x2 = runif(100);
#'  y = x1 * x2^10; y = y + rnorm(length(y), 0, 0.01) # generate a dataset
#'  results = FPIR.1twoway (y, x1, x2)
#'  Exp1 = results[[2]]$Ex1; Exp2 = results[[2]]$Ex2
#'  results.2 = FPIR1twowaytune (y, x1, x2, Exp1, Exp2)
#'  Ex1.tune   = as.numeric(results.2[[2]][1])
#'  Ex2.tune   = as.numeric(results.2[[2]][2])
#'  results.2[[2]] # parameters for the best model
#'  results.2[[3]] # regression coefficients
#'  results.2[[4]] # adjusted R square
#' @export

FPIR1twowaytune <- function(y, x1, x2, Exp1, Exp2) { # two explanatory variables x1 and x2 and one dependent variable y
  out <- list()
  LIST <- c(-52.5,-48.6,-44.85,-41.25,-37.8,-34.5,-31.35,-28.35,-25.5,-22.8,-20.25,-17.85,-15.6,-13.5,
            -11.55,-9.75,-8.1,-6.6,-5.25,-4.05,-3,-2.1,-1.35,-1,-0.75,-0.3,-0.1,0,0.1,0.3,0.75,1,1.35,2.1,3,4.05,5.25,
            6.6,8.1,9.75,11.55,13.5,15.6,17.85,20.25,22.8,25.5,28.35,31.35,34.5,37.8,41.25,44.85,48.6,52.5)

  if (abs(Exp1) < 52.5) {
    order_1 = LIST==Exp1;
    order_1_L = order_1[-1];              order_1_L = c(order_1_L, FALSE); Exp1_L = LIST[order_1_L]
    order_1_U = order_1[-length(order_1)];order_1_U = c(FALSE, order_1_U); Exp1_U = LIST[order_1_U]
    Exp1_s = seq(Exp1_L+(Exp1-Exp1_L)/2, Exp1+(Exp1_U-Exp1)/2, length=10)
  }
  if (Exp1 == -52.5) {Exp1_s = seq(-56, -50.5, length=10)}
  if (Exp1 ==  52.5) {Exp1_s = seq(50.5, 56,   length=10)}

  if (abs(Exp2) < 52.5) {
    order_2 = LIST==Exp2;
    order_2_L = order_2[-1];              order_2_L = c(order_2_L, FALSE); Exp2_L = LIST[order_2_L]
    order_2_U = order_2[-length(order_2)];order_2_U = c(FALSE, order_2_U); Exp2_U = LIST[order_2_U]
    Exp2_s = seq(Exp2_L+(Exp2-Exp2_L)/2, Exp2+(Exp2_U-Exp2)/2, length=10)
  }
  if (Exp2 == -52.5) {Exp2_s = seq(-56, -50.5, length=10)}
  if (Exp2 ==  52.5) {Exp2_s = seq(50.5, 56,   length=10)}


  Exps <- expand.grid(Exp1_s, Exp2_s)
  names(Exps) <- c('Ex1','Ex2')
  Len <- nrow(Exps) # 100
  out[[1]] <- cbind(Exps, Rsq_T=NA, Rsq_I=NA,  coef_In=NA, coef_x1=NA, coef_x2=NA, coef_I=NA)
  for(K in 1:Len){
    inter <- x1^Exps[K, 1] * x2^Exps[K, 2]
    fit <- lm(y ~ x1 + x2 + inter)
    reg <- summary(fit)
    out[[1]]$Rsq_T[K] <- reg$adj.r.squared
    out[[1]]$Rsq_I[K] <- anova(fit)[3,2]/sum(anova(fit)[1:4,2])
    out[[1]]$coef_In[K]<-coef(fit)[1]
    out[[1]]$coef_x1[K]<-coef(fit)[2]
    out[[1]]$coef_x2[K]<-coef(fit)[3]
    out[[1]]$coef_I[K ]<-coef(fit)[4]
  }
  out[[2]] <- out[[1]][out[[1]]$Rsq_T == max(out[[1]]$Rsq_T),] # the best model
  if (nrow(out[[2]])>1) out[[2]] <- out[[2]][1,]
  if(is.na(unique(out[[2]][1]))) {
    out[[2]] <- out[[2]][1, ]
    out[[3]] = out[[4]] = out[[5]] <- "NA"
  }
  else
  {inter <- x1^out[[2]]$Ex1 * x2^out[[2]]$Ex2
  fit <- lm(y ~ x1 + x2 + inter)
  sum.fit <- summary(fit)
  out[[3]] <- fit$coefficients
  out[[4]] <- sum.fit$adj.r.squared
  if(is.na(fit$coefficients[4]))  out[[5]] <- 1 else out[[5]] <- sum.fit[[4]][4,4]}
  print(out[[2]][1:2])
  return(out)
}

