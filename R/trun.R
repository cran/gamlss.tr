#----------------------------------------------------------------------------------------
trun <-function (par = c(0), 
              family = "NO", 
                name = "tr", 
                type = c("left", "right", "both"),
               local = TRUE,
               delta = NULL, 
                ...)
{
#------------------------------------------
     type <- match.arg(type)
     fam  <- as.gamlss.family(family)
    fname <- fam$family[[1]] 
   family <- c("None", "None")
     dfun <- paste(paste("d",fname,sep=""), name, sep="")
   dorfun <- paste("d",fname,sep="")
   porfun <- paste("p",fname,sep="")
     pfun <- paste(paste("p",fname,sep=""), name, sep="")
   #  qfun <- paste(paste("q",fname,sep=""), name, sep="")
   #  rfun <- paste(paste("r",fname,sep=""), name, sep="")           
 # fname <- family # as.name(family)
if (local)
 {
#--trying to get gamlss sys.frame--  
     rexpr<-regexpr("gamlss",sys.calls())
for (i in 1:length(rexpr)){ 
    position <- i 
    if (rexpr[i]==1) break}
gamlss.environment <- sys.frame(position)      
#--end here------------------------
 }
 else gamlss.environment <- sys.frame(0)
#   generate d within gamlss
    eval(dummy <- trun.d(par, family = fname, type = type, ...))
    eval(call("<-",as.name(dfun),dummy), envir=gamlss.environment)# parent.frame(n = 1)
# eval(call("<-",as.name(dfun),dummy), envir=sys.frame(1)) # to put it in gamlss
# generate p within gamlss
    eval(dummy <- trun.p(par, family = fname, type = type, ...))
    eval(call("<-",as.name(pfun),dummy), envir=gamlss.environment)# parent.frame(n = 1)
# rename the family 
   family[[1]] <- paste(paste(fname, name, sep=""))
   family[[2]] <- paste(type, "truncated",fam$family[[2]])
    fam$family <- family
# Global deviance increment  
           sGD <- gsub(dorfun, dfun, deparse(body(fam$G.dev.incr)))
  body(fam$G.dev.incr) <- parse(text=sGD)
# get the no of parameters  
        nopar <- fam$nopar
# check for the delta
 if (length(delta)==0) delta <- rep(NA,nopar) 
 if (length(delta)==1) delta <- rep(delta,nopar)
 if (length(delta)!=nopar)  stop("delta should be the same length the parameters in the family ") 
# now change the first derivatives
  switch(nopar,  
          { 
      fam$dldm <- function() as.vector(attr(numeric.deriv(TEST(y, mu, log=TRUE), "mu", delta=NULL), "gradient")) 
           sMU <- sub("TEST", dfun, body(fam$dldm))
if (!is.na(delta[1])) sMU <- sub("NULL",  as.character(delta[1]), sMU) 
body(fam$dldm) <- parse(text=sMU[length(sMU)])
          sres <- gsub(porfun, pfun,  deparse(fam$rqres))
          if  (fam$type=="Discrete") sres <- gsub("ymin = 0",  paste("ymin =",par[1]+1),  sres) 
          sres <- gsub("expression", "",  sres)
     fam$rqres <- parse(text=sres)  
          },
          {   
      fam$dldm <- function() as.vector(attr(numeric.deriv(TEST(y, mu, sigma, log=TRUE), "mu", delta=NULL), "gradient"))
      fam$dldd <- function() as.vector(attr(numeric.deriv(TEST(y, mu, sigma, log=TRUE), "sigma", delta=NULL), "gradient"))
           sMU <- sub("TEST", dfun, body(fam$dldm))
if (!is.na(delta[1])) sMU <- sub("NULL",  as.character(delta[1]), sMU) 
body(fam$dldm) <- parse(text=sMU[length(sMU)])   
        sSIGMA <- sub("TEST", dfun, body(fam$dldd))
if (!is.na(delta[2])) sSIGMA <- sub("NULL",  as.character(delta[2]), sSIGMA)
body(fam$dldd) <- parse(text=sSIGMA[length(sSIGMA)]) 
          sres <- gsub(porfun, pfun,  deparse(fam$rqres))
          if  (fam$type=="Discrete") sres <- gsub("ymin = 0",  paste("ymin =",par[1]+1),  sres) 
          sres <- gsub("expression", "",  sres)
     fam$rqres <- parse(text=sres)    
           },
           {   
      fam$dldm <- function() as.vector(attr(numeric.deriv(TEST(y, mu, sigma, nu, log=TRUE), "mu", delta=NULL), "gradient"))
      fam$dldd <- function() as.vector(attr(numeric.deriv(TEST(y, mu, sigma, nu, log=TRUE), "sigma", delta=NULL), "gradient"))
      fam$dldv <- function() as.vector(attr(numeric.deriv(TEST(y, mu, sigma, nu, log=TRUE), "nu", delta=NULL), "gradient"))
           sMU <- sub("TEST", dfun, body(fam$dldm))
if (!is.na(delta[1]))sMU <- sub("NULL",  as.character(delta[1]), sMU)          
body(fam$dldm) <- parse(text=sMU[length(sMU)])  
        sSIGMA <- sub("TEST", dfun, body(fam$dldd))
if (!is.na(delta[2])) sSIGMA <- sub("NULL",  as.character(delta[2]), sSIGMA) 
body(fam$dldd) <- parse(text=sSIGMA[length(sSIGMA)])  
           sNU <- sub("TEST", dfun, body(fam$dldv))
if (!is.na(delta[3])) sNU <- sub("NULL",  as.character(delta[3]), sNU)
body(fam$dldv) <- parse(text=sNU[length(sNU)])
          sres <- gsub(porfun, pfun,  deparse(fam$rqres))
          if  (fam$type=="Discrete") sres <- gsub("ymin = 0",  paste("ymin =",par[1]+1),  sres) 
          sres <- gsub("expression", "",  sres)
     fam$rqres <- parse(text=sres)    
           },
           {
      fam$dldm <- function() as.vector(attr(numeric.deriv(TEST(y, mu, sigma, nu, tau, log=TRUE), "mu", delta=NULL), "gradient"))
      fam$dldd <- function() as.vector(attr(numeric.deriv(TEST(y, mu, sigma, nu, tau, log=TRUE), "sigma", delta=NULL), "gradient"))
      fam$dldv <- function() as.vector(attr(numeric.deriv(TEST(y, mu, sigma, nu, tau, log=TRUE), "nu", delta=NULL), "gradient"))
      fam$dldt <- function() as.vector(attr(numeric.deriv(TEST(y, mu, sigma, nu, tau, log=TRUE), "tau", delta=NULL), "gradient"))
           sMU <- sub("TEST", dfun, body(fam$dldm))
if (!is.na(delta[1])) sMU <- sub("NULL",  as.character(delta[1]), sMU)      
body(fam$dldm) <- parse(text=sMU[length(sMU)])   
        sSIGMA <- sub("TEST", dfun, body(fam$dldd))
if (!is.na(delta[2])) sSIGMA <- sub("NULL",  as.character(delta[2]), sSIGMA) 
body(fam$dldd) <- parse(text=sSIGMA[length(sSIGMA)])  
           sNU <- sub("TEST", dfun, body(fam$dldv))
if (!is.na(delta[3])) sNU <- sub("NULL",  as.character(delta[3]), sNU)           
body(fam$dldv) <- parse(text=sNU[length(sNU)]) 
          sTAU <- sub("TEST", dfun, body(fam$dldt))
if (!is.na(delta[4])) sTAU <- sub("NULL",  as.character(delta[4]), sTAU)
body(fam$dldt) <- parse(text=sTAU[length(sTAU)]) 
          sres <- gsub(porfun, pfun,  deparse(fam$rqres))
          if  (fam$type=="Discrete") sres <- gsub("ymin = 0",  paste("ymin =",par[1]+1),  sres) 
          sres <- gsub("expression", "",  sres)
     fam$rqres <- parse(text=sres)    
           })
      fam 
}
