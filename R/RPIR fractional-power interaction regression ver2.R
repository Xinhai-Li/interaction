#' Fit fractional-power interaction regression (FPIR) for regression y ~ x1 + x2 + x1^M * x2^N
#'
#' @description This function estimates the values of exponents (M and N) in Fractional-power interaction regression
#'  (FPIR)，the formula is: y ~ x1 + x2 + x1^M * x2^N.  Return model parameters such as regression coefficients,
#'  sum of squares of the interaction term, and total R squares of the model.
#'
#' @author Xinhai Li (Xinhai_li_edu@126.com)
#'
#' @param y the dependent variable in regression y ~ x1 + x2 + x1^M *x2^N
#' @param x1 one of the independent variables in regression y ~ x1 + x2 + x1^M *x2^N
#' @param x2 another independent variable in regression y ~ x1 + x2 + x1^M *x2^N
#'
#' @return return a list including a data.frame with columns of M, N, R square value of the whole model,
#' proportion of the variance explained by the interaction term, coefficients of all terms
#' (the interception, x1, x2, and the interaction).
#'
#' @examples
#'
#' # Simulated data
#' x1 = runif(100); x2 = runif(100);
#' y = x1 * x2^10; y = y + rnorm(length(y), 0, 0.01) # generate a dataset
#' results = FPIR1twoway (y, x1, x2)
#' Exp1 = results[[2]]$Ex1; Exp2 = results[[2]]$Ex2
#' results.2 = FPIR1twowaytune (y, x1, x2, Exp1, Exp2)
#' results.2[[1]] # parameters for 10000 models
#' results.2[[2]] # parameters for the best model
#' results.2[[3]] # regression coefficient
#' results.2[[4]] # adjusted R square#'
#'
#' attach(trees)
#' # First round fitting
#' results = FPIR1twoway(trees$Volume, trees$Girth, trees$Height) # 1.35, 0.3
#' # Second round fitting based on the first round fitting
#' results2 = FPIR1twowaytune(trees$Volume, trees$Girth, trees$Height, 1.35, 0.3)
#' @export
#=============================================================================================================


#=============================================================================================================
FPIR1twoway <- function(y, x1, x2) { # two explanatory variables x1 and x2 and one dependent variable y
  out <- list()
  Exp1 = Exp2 <- c(-52.5,-48.6,-44.85,-41.25,-37.8,-34.5,-31.35,-28.35,-25.5,-22.8,-20.25,-17.85,-15.6,-13.5,
                   -11.55,-9.75,-8.1,-6.6,-5.25,-4.05,-3,-2.1,-1.35,-1,-0.75,-0.3,-0.1,0,0.1,0.3,0.75,1,1.35,2.1,3,4.05,5.25,
                   6.6,8.1,9.75,11.55,13.5,15.6,17.85,20.25,22.8,25.5,28.35,31.35,34.5,37.8,41.25,44.85,48.6,52.5)
  Exps <- expand.grid(Exp1, Exp2)
  names(Exps) <- c('Ex1','Ex2')
  Len <- nrow(Exps) # 10000
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
    print(paste(K/Len*100, "% finished", sep=' '))
    # cat(K, " ")
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
  return(out)
}
#=============================================================================================================




#' Tune fractional-power interaction regression (FPIR) for regression y ~ x1 + x2 + x1^M * x2^N
#'
#' @description This function tune the estimated values of exponents (M and N) from FPIR1twoway(y, x1, x2).
#'
#' @author Xinhai Li (Xinhai_li_edu@126.com)
#'
#' @param y the dependent variable in regression y ~ x1 + x2 + x1^M *x2^N
#' @param x1 one of the independent variables in regression y ~ x1 + x2 + x1^M *x2^N
#' @param x2 another independent variable in regression y ~ x1 + x2 + x1^M *x2^N
#' @param Exp1 the initial value of M in regression y ~ x1 + x2 + x1^M *x2^N
#' @param Exp2 the initial value of N in regression y ~ x1 + x2 + x1^M *x2^N
#'
#' @return return a list including tuned values of M and N.
#'
#' @examples
#'
#' # Simulated data
#' x1 = runif(100); x2 = runif(100);
#' y = x1 * x2^10; y = y + rnorm(length(y), 0, 0.01) # generate a dataset
#' results = FPIR1twoway (y, x1, x2)
#' Exp1 = results[[2]]$Ex1; Exp2 = results[[2]]$Ex2
#' results.2 = FPIR1twowaytune (y, x1, x2, Exp1, Exp2)
#' results.2[[1]] # parameters for 10000 models
#' results.2[[2]] # parameters for the best model
#' results.2[[3]] # regression coefficient
#' results.2[[4]] # adjusted R square#'
#'
#' attach(trees)
#' # First round fitting
#' results = FPIR1twoway(trees$Volume, trees$Girth, trees$Height) # 1.35, 0.3
#' # Second round fitting based on the first round fitting
#' results2 = FPIR1twowaytune(trees$Volume, trees$Girth, trees$Height, 1.35, 0.3)
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
  return(out)
}
#=============================================================================================================



# IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
# how to define initial exponents (more values around 1)

# two-way interaction
# P1 = numeric(); P1[1]=0; step=0.03
# for (i in 2:26)   {P1[i] =  P1[i-1]+step*i*5}
# P2 = P1*-1
# POWER = unique(sort(c(P1, P2, -0.1, 0.1, -1, 1)));   plot(POWER)
# POWER=t(POWER)
# write.csv(POWER, "d:/out.csv", row.names=F)

# three way interaction (apart more)
# P1 = numeric(); P1[1]=0; step=0.1
# for (i in 2:10)   {P1[i] =  P1[i-1]+step*i*10}
# P2 = P1*-1
# POWER = unique(sort(c(P1, P2, -0.5, 0.5, -1, 1)));   plot(POWER)
# IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII





#' Fractional-power interaction regression (FPIR) for regression y ~ x1 + x2 + x1^M * x2^N using Bayesian method
#'
#' @description This function estimates the values of exponents (M and N) in
#'  Fractional-power interaction regression (FPIR) for regression y ~ x1 + x2 + x1^M * x2^N,
#'  return regression coefficients.
#'
#' @author Xinhai Li (Xinhai_li_edu@126.com)
#'
#' @param y the dependent variable in regression y ~ x1 + x2 + x1^M *x2^N
#' @param x1 one of the independent variables in regression y ~ x1 + x2 + x1^M *x2^N
#' @param x2 another independent variable in regression y ~ x1 + x2 + x1^M *x2^N
#' @param WinBUGS the install direction of software WinBUGS
#' @param digits controls the number of significant digits to print when printing numeric values (intercept and regression coefficients in formula).
#'
#' @return return a list including a data.frame with estimated model coefficients, and the formula of y ~ x1 + x2 + x1^M * x2^N
#'
#' @examples
#' data(trees)
#' results = interaction_bayes(trees$Volume, trees$Girth, trees$Height, WinBUGS = "d:/softwares/WinBUGS14/", digits=2)
#' results
#'
#' @importFrom R2WinBUGS bugs
#'
#' @export
##=======================================================================================
interaction_bayes <- function(y, x1, x2, WinBUGS, digits=3){
  n <- length(y)

  # WinBUGS model ##===##
  interaction <- function(){
    for(i in 1:n){
      y[i]  ~ dnorm(mu[i], tau)
      mu[i]  <- b0 + b1*x1[i] + b2*x2[i]  + b3 * pow(x1[i], b4) * pow(x2[i], b5)
    }
    b0 ~ dnorm(0,.01)
    b1 ~ dnorm(0,.01)
    b2 ~ dnorm(0,.01)
    b3 ~ dnorm(0,.01)
    b4 ~ dnorm(0,.01)
    b5 ~ dnorm(0,.01)
    tau ~ dgamma(0.01, 0.01)
  }
  interaction_file <- file.path(tempdir(), "interaction.txt")
  ## write model file:
  write.model(interaction, interaction_file)
  ##===##

  data <- list ( "y","x1","x2","n")
  inits <- function()  list(b0=1, b1=1, b2=1, b3=1, b4=1, b5=0.5, tau=1)
  parameters <- c("b0", "b1", "b2", "b3","b4", "b5", "tau")

  out <- bugs(data, inits, parameters, model.file = interaction_file,
              n.chain=3, n.burnin=1000, n.iter=10000, debug=T,
              bugs.directory = WinBUGS)
  MEANS <- round(as.numeric(out$mean), digits)
  Formula <- paste("y = ", MEANS[1], " + ", MEANS[2], "x1", " + ", MEANS[3], "x2", " + ", MEANS[4], "x1^", MEANS[5], "x2^", MEANS[6], sep="")
  Formu <- sub("+ -", "- ", Formula, fixed = T) # fixed = T make sure it works
  OUT <- list()
  OUT[[1]] <- out$summary
  OUT[[2]] <- Formu
  names(OUT) <- c("Summary table for model coefficients", "Formula")
  return(OUT)
}
##=======================================================================================
