#' A class for linear regression with many methods for corresponding values calculations
#'
#' @field y vector The orginial labels of observations in data
#' @field X matrix The features of all observations in data
#' @field beta vector The coefficients in the sample
#' @field y_hat matrix The estimated labels
#' @field e_hat matrix The standard errors between orginial and estimated labels
#' @field df numeric The degree of freedom of data
#' @field var_resid numeric The residual variance
#' @field var_hat matrix The variance of the regression coefficients
#' @field t_val vector The t-values for each coefficient
#' @field p_val vector The p-values for each coefficient
#' @field formula character The formula used in the model
#' @field data data.frame The sample data
#'
#' @import methods
#' @importFrom ggplot2 ggplot
#' @importFrom plyr is.formula
#' @importFrom gridExtra grid.arrange
#' @export linreg
#' @exportClass linreg

linreg <- setRefClass("linreg",
  fields = list(
    y = "numeric",
    X = "matrix",
    y_hat = "matrix",
    beta = "vector",
    e_hat = "matrix",
    df = "numeric",
    var_resid = "numeric",
    var_hat = "matrix",
    t_val = "vector",
    p_val = "vector",
    formula = "formula",
    data = "data.frame",
    name = "character",
    f = "character"
  ),
  methods = list(
    # The initialize methods
    initialize = function(formula, data){
      stopifnot(plyr::is.formula(formula),is.data.frame(data) )
      cat("User::initialize")
      f <<- Reduce(paste, deparse(formula))
      formula <<- formula
      name <<- deparse(substitute(data))
      data <<- data
      X <<- model.matrix(formula,data)
      y <<- data[,all.vars(formula)[1]]
      beta <<- round(as.vector(solve(t(X)%*%X)%*%t(X)%*%y),3)
      y_hat <<- round(as.matrix(X%*%beta),3)
      e_hat <<- y-y_hat
      df <<- length(y) - ncol(X)
      var_resid <<- as.numeric((t(e_hat)%*%e_hat)/df)
      var_hat <<- var_resid * solve(t(X) %*% X)
      t_val <<- beta/(sqrt(diag(var_hat)))
      p_val <<- 2*pt(abs(t_val), df,lower.tail = FALSE)
      names(p_val) <<- colnames(X)
      names(beta) <<- colnames(X)
    },
    #Estimation of beta and their variances by QR decomposition
    qr_beta = function(){
      Q <- qr.Q(qr(X))
      R <- qr.R(qr(X))
      beta <<- solve(R)%*%t(Q)%*%y
      var_hat <<- var_resid*solve(R)%*%t(solve(R))
    },
    #Regressions coefficients
    coef = function(){
      return(beta)
    },
    #The fitted values
    pred = function(){
      return(y_hat)
    },
    #The residuals
    resid = function(){
      return(e_hat)
    },
    #The degrees of freedom
    freedomdegree = function(){
      return(df)
    },
    #The residual variance
    residualvariance = function(){
      return(var_resid)
    },
    #The variance of the regression coefficients
    coefvariance = function(){
      return(var_hat)
    },
    #The t-values for each coefficient
    t_values = function(){
      return(t_val)
    },
    #The p-values for each coefficient
    p_values = function(){
      return(p_Val)
    },
    #The print function
    print = function(){
      cat(paste("linreg(formula = ",f,", data = ", name,")\n",sep=""))
      base::print(beta)
    },
    #The summary function
    summary = function(){
      std_e <- round((sqrt(diag(var_hat))),3)
      m <- as.data.frame(matrix(c(round(beta,3), std_e, round(t_val,3), format(p_val, scientific = 2)), ncol = 4))
      for(i in 1:ncol(X)){
        if(p_val[i]>=0&&p_val[i]<0.001){
          m[i,5]<-"***"
        }else if(
          p_val[i]>=0.001&&p_val[i]<0.01){
          m[i,5]<-"**"
        }else if(
          p_val[i]>=0.01&&p_val[i]<0.05){
          m[i,5]<-"*"
        }else if(
          p_val[i]>=0.05&&p_val[i]<0.1){
          m[i,5]<-"."
        }else if(
          p_val[i]>=0.1&&p_val[i]<1){
          m[i,5]<-""
        }
      }
      colnames(m) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)","")
      rownames(m) <- colnames(X)
      base::print(m)
      base::print(paste("Residual standard error:", round(sqrt(var_resid),3), "on", df, "degrees of freedom"))
    },

    #The plot function
    plot=function(){
      d <- as.data.frame(cbind(y_hat,e_hat, sqrt(abs((e_hat-mean(e_hat))/sd(e_hat)))))

      g1<-ggplot2::ggplot(data=d,ggplot2::aes(x=d[,1],y=d[,2])) +
        ggplot2::geom_point(shape=1, size=4)+
        ggplot2::stat_summary(ggplot2::aes(x=y_hat,group=1),fun.y=median, color="red", geom="line")+
        ggplot2::geom_hline(yintercept = 0, col="grey", linetype="dotted")+
        ggplot2::labs(x="Fitted valies or Predictions",y="Residuals")+
        ggplot2::ggtitle("Residuals vs Fitted Plot")+
        ggplot2::scale_y_continuous()+
        ggplot2::scale_color_discrete("Regression")+
        ggplot2::theme(
          panel.background  = ggplot2::element_blank(),
          plot.background = ggplot2::element_rect(fill="lightskyblue1", color=NA),
          legend.background = ggplot2::element_rect(fill="transparent", color=NA),
          legend.key = ggplot2::element_rect(fill="transparent", color=NA)
        )


      g2<-ggplot2::ggplot(d,ggplot2::aes(y=d[,3],x=d[,1]))+
        ggplot2::geom_point(shape=1, size=4)+
        ggplot2::stat_summary(ggplot2::aes(x=y_hat,group=1),fun.y=median, color="red", geom="line")+
        ggplot2::labs(x="Fitted valies or Predictions",y=" sqrt(|Standardized Residuals|)")+
        ggplot2::ggtitle("sqrt(|Standardized Residuals|) vs Fitted Plot")+
        ggplot2::scale_y_continuous()+
        ggplot2::scale_color_discrete("Regression")+
        ggplot2::theme(
          panel.background  = ggplot2::element_blank(),
          plot.background = ggplot2::element_rect(fill="lightskyblue1", color=NA),
          legend.background = ggplot2::element_rect(fill="transparent", color=NA),
          legend.key = ggplot2::element_rect(fill="transparent", color=NA)
        )

      gg<-gridExtra::grid.arrange(g1, g2)

    }
  )

)

