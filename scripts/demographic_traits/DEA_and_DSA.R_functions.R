# Limma function####
limma_lm <- function(fit, covariate, covariate_data){
    v.contrast <- rep(0,ncol(fit$design))
    if(is.factor(covariate_data[,covariate])){
        if(covariate == "Sex"){
            if(table(covariate_data$Sex)["1"] >= 10  & table(covariate_data$Sex)["2"] >= 10 ){   
                v.contrast[ which( colnames(fit$design) == "Sex2") ] <- 1 # contrast againts basal level which is male
                contrast.matrix <- cbind( "C1" = v.contrast)
            }else{
                return(NA)
            }
        }else if(covariate=="Ancestry"){
            if( table(covariate_data$Ancestry)["AFR"] >= 10  & table(covariate_data$Ancestry)["EUR"] >= 10 ){
                v.contrast[ which( colnames(fit$design) == "AncestryAFR") ] <- 1 # contrast againts basal level which is EUR
                contrast.matrix <- cbind( "C1" = v.contrast)
            }else{
                return(NA)
            }
        }else{
            return(NULL)
        }
    }else{ # Continuous variable
        v.contrast[ which(colnames(fit$design) == covariate)] <- 1 #
        contrast.matrix <- cbind( "C1" = v.contrast)
    }
    fitConstrasts <- contrasts.fit(fit,contrast.matrix)
    eb = eBayes(fitConstrasts)
    tt.smart.sv <- topTable(eb,adjust.method = "BH",number=Inf)
    return(tt.smart.sv)
}  

limma_lm.n3 <- function(fit, covariate, covariate_data){
    v.contrast <- rep(0,ncol(fit$design))
    if(is.factor(covariate_data[,covariate])){
        if(covariate == "Sex"){
            if(table(covariate_data$Sex)["1"] >= 3  & table(covariate_data$Sex)["2"] >= 3 ){   
                v.contrast[ which( colnames(fit$design) == "Sex2") ] <- 1 # contrast againts basal level which is male
                contrast.matrix <- cbind( "C1" = v.contrast)
            }else{
                return(NA)
            }
        }else if(covariate=="Ancestry"){
            if( table(covariate_data$Ancestry)["AFR"] >= 3  & table(covariate_data$Ancestry)["EUR"] >= 3 ){
                v.contrast[ which( colnames(fit$design) == "AncestryAFR") ] <- 1 # contrast againts basal level which is EUR
                contrast.matrix <- cbind( "C1" = v.contrast)
            }else{
                return(NA)
            }
        }else{
            return(NULL)
        }
    }else{ # Continuous variable
        v.contrast[ which(colnames(fit$design) == covariate)] <- 1 #
        contrast.matrix <- cbind( "C1" = v.contrast)
    }
    fitConstrasts <- contrasts.fit(fit,contrast.matrix)
    eb = eBayes(fitConstrasts)
    tt.smart.sv <- topTable(eb,adjust.method = "BH",number=Inf)
    return(tt.smart.sv)
}  

limma_lm.sum2zero <- function(fit, covariate, covariate_data){
    v.contrast <- rep(0,ncol(fit$design))
    # When sum to 0, we need to multiply times 2 the factors
    if(is.factor(covariate_data[,covariate])){
        if(covariate == "Sex"){
            if(table(covariate_data$Sex)["1"] >= 10  & table(covariate_data$Sex)["2"] >= 10 ){   
                v.contrast[ which( colnames(fit$design) == "Sex1") ] <- 2 # contrast againts basal level which is male
                contrast.matrix <- cbind( "C1" = v.contrast)
            }else{
                return(NA)
            }
        }else if(covariate=="Ancestry"){
            if( table(covariate_data$Ancestry)["AFR"] >= 10  & table(covariate_data$Ancestry)["EUR"] >= 10 ){
                v.contrast[ which( colnames(fit$design) == "Ancestry1") ] <- 2 # contrast against basal level which is EUR 
                contrast.matrix <- cbind( "C1" = v.contrast)
            }else{
                return(NA)
            }
        }else{
            return(NULL)
        }
    }else{ # Continuous variable
        v.contrast[ which(colnames(fit$design) == covariate)] <- 1 #
        contrast.matrix <- cbind( "C1" = v.contrast)
    }
    fitConstrasts <- contrasts.fit(fit,contrast.matrix)
    eb = eBayes(fitConstrasts)
    tt.smart.sv <- topTable(eb,adjust.method = "BH",number=Inf)
    return(tt.smart.sv)
}  

limma_lm.interactions <- function(fit, interaction.term){
    print(interaction.term)
    v.contrast <- rep(0,ncol(fit$design))
    if(interaction.term == "Age:Ancestry"){ # Continuous variable -> Age is multiple by 0 if AFR or Male. So basically considering effect of Age in EUR or Female
        v.contrast[ which(colnames(fit$design) == "Age:Ancestry1")] <- 1 #
        contrast.matrix <- cbind( "C1" = v.contrast)
    }else if(interaction.term == "Age:Sex"){ # Continuous variable -> Age is multiple by 0 if AFR or Male. So basically considering effect of Age in EUR or Female
        v.contrast[ which(colnames(fit$design) == "Age:Sex1")] <- 1 #
        contrast.matrix <- cbind( "C1" = v.contrast)
    }else if(interaction.term == "Age:BMI"){ # Continuous variable -> Age is multiple by 0 if AFR or Male. So basically considering effect of Age in EUR or Female
        v.contrast[ which(colnames(fit$design) == "Age:BMI")] <- 1 #
        contrast.matrix <- cbind( "C1" = v.contrast)
    }else if(interaction.term == "Ancestry:BMI"){ # Continuous variable -> Age is multiple by 0 if AFR or Male. So basically considering effect of Age in EUR or Female
        v.contrast[ which(colnames(fit$design) == "Ancestry1:BMI")] <- 1 #
        contrast.matrix <- cbind( "C1" = v.contrast)
    }else if(interaction.term == "Ancestry:Sex"){ # Continuous variable -> Age is multiple by 0 if AFR or Male. So basically considering effect of Age in EUR or Female
        if("Ancestry1:Sex1" %in% colnames(fit$design)){
            v.contrast[ which(colnames(fit$design) == "Ancestry1:Sex1")] <- 1 #
            contrast.matrix <- cbind( "C1" = v.contrast)
        }else{
           v.contrast[ which(colnames(fit$design) == "Sex1:Ancestry1")] <- 1 #
           contrast.matrix <- cbind( "C1" = v.contrast)
        }
    }else{
        if("Sex1:BMI" %in% colnames(fit$design)){
            v.contrast[ which( colnames(fit$design) == "Sex1:BMI") ] <- 1 # contrast against basal level which not EUR not female
            contrast.matrix <- cbind( "C1" = v.contrast)
        }else{
            v.contrast[ which( colnames(fit$design) == "BMI:Sex1") ] <- 1 # contrast against basal level which not EUR not female
            contrast.matrix <- cbind( "C1" = v.contrast)
        }
    }
    fitConstrasts <- contrasts.fit(fit,contrast.matrix)
    eb = eBayes(fitConstrasts)
    tt.smart.sv <- topTable(eb,adjust.method = "BH",number=Inf)
    return(tt.smart.sv)
}  

# Filtering functions ####
#nearZeroVar_fun <- function(i,samples, psi){nearZeroVar(as.numeric(psi[i,samples]), freqCut = 80/20, uniqueCut = 0)}

# Residuals functions ####
mountBeta <- function(glmBatch,batchReg,mydata) {
    # mount beta coefficients
    
    nbeta <- 1
    ncoef <- 2
    beta  <- numeric(1)
    
    for (i in batchReg) {
        if (!is.factor(mydata[,i])) {
            beta[nbeta] <- glmBatch$coefficients[ncoef]
            nbeta <- nbeta + 1
            ncoef <- ncoef + 1
        }
        if (is.factor(mydata[,i])) {
            nlev        <- nlevels(mydata[,i])
            b1          <- -sum(glmBatch$coefficients[ncoef:(ncoef+nlev-2)])/nlev
            beta[nbeta] <- b1
            beta[(nbeta+1):(nbeta+nlev-1)] <- b1+glmBatch$coefficients[ncoef:(ncoef+nlev-2)] 
            nbeta <- nbeta + nlev
            ncoef <- ncoef + (nlev-1)
        }
    }
    return(beta)
}

computeResiduals <- function(mydata, batchReg, traitReg) {
    
    #  traitRegressors <- regressors[!regressors %in% batchReg]
    regressors <- c(batchReg,traitReg)
    myformeq <- paste("y", paste(colnames(mydata)[regressors], collapse=" + "), sep=" ~ ")
    glmFull  <- glm(myformeq, data = mydata,family = quasibinomial('logit'),control = list(maxit = 100) )
    
    ylogit    <- log(mydata$y/(1-mydata$y))
    myformeqB <- paste("~ ", paste(colnames(mydata)[batchReg], collapse=" + ")," - 1")
    #  AUZ       <- model.matrix(as.formula(myformeqB), data=mydata[,regressors])
    if (length(batchReg) > 1) index   <- sapply(mydata[,batchReg],is.factor)
    if (length(batchReg) == 1) index  <- is.factor(mydata[,batchReg])
    fact      <- batchReg[index]
    if (length(fact) > 1) AUZ  <- model.matrix(as.formula(myformeqB), data=mydata[,regressors],
                                               contrasts.arg = lapply(mydata[,fact], contrasts, contrasts=FALSE))
    if (length(fact) <2 ) AUZ  <- model.matrix(as.formula(myformeqB), data=mydata[,regressors])
    
    #-- from model with only batch effects, taking out these batch effects
    
    myformeq <- as.formula(paste("y", paste(colnames(mydata)[batchReg], collapse=" + "), sep=" ~ "))
    glmBatch <- glm(myformeq , data = mydata,family = quasibinomial('logit'),control = list(maxit = 100) )
    
    # only glm with batch effects
    beta <- mountBeta(glmBatch,batchReg,mydata) 
    predBatch   <- AUZ %*% matrix(beta)
    
    ylogitRes1  <- ylogit - predBatch
    yRes1       <- 1/(1+exp(-ylogitRes1))
    myformeq    <- as.formula(paste("yRes1", paste(colnames(mydata)[traitReg], collapse=" + "), sep=" ~ "))
    glmRes1     <- glm(myformeq, data = mydata,family = quasibinomial('logit'),control = list(maxit = 100) )
    
    # -- from full model, taking out the batch effects
    
    # beta2       <- mountBeta(glmFull,batchReg,mydata) 
    # predBatch2  <- AUZ %*% matrix(beta2)
    # ylogitRes2  <-  ylogit - predBatch2
    # yRes2       <-  1/(1+exp(-ylogitRes2))
    # myformeq    <- as.formula(paste("yRes2", paste(colnames(mydata)[traitReg], collapse=" + "), sep=" ~ "))
    # glmRes2     <- glm(myformeq, data = mydata,family = quasibinomial('logit'),control = list(maxit = 100) )
    
    SStot   <- sum((mydata$y-mean(mydata$y))^2)
    SStotR1 <- sum((yRes1-mean(yRes1))^2)
    #SStotR2 <- sum((yRes2-mean(yRes2))^2)
    R2Full  <- round((SStot - sum((glmFull$fitted.values- mydata$y)^2))/SStot,digits=5)
    R2Batch <- round((SStot - sum((glmBatch$fitted.values- mydata$y)^2))/SStot,digits=5)
    R2Res1  <- round((SStotR1 - sum((glmRes1$fitted.values- yRes1)^2))/SStotR1,digits=5)
    #R2Res2  <- round((SStotR2 - sum((glmRes2$fitted.values- yRes2)^2)) /SStotR2,digits=5)
    
    #tableCoefs <- data.frame(Full=glmFull$coefficients,Batch=NA,Res1=NA, Res2=NA)
    tableCoefs <- data.frame(Full=glmFull$coefficients,Batch=NA,Res1=NA)
    index   <- match(names(glmBatch$coefficients),names(glmFull$coefficients))
    tableCoefs[index,2] <- glmBatch$coefficients
    
    index   <- match(names(glmRes1$coefficients),names(glmFull$coefficients))
    tableCoefs[index,3] <- glmRes1$coefficients
    
    #index   <- match(names(glmRes2$coefficients),names(glmFull$coefficients))
    #tableCoefs[index,4] <- glmRes2$coefficients
    
    difGlobalMod1 <- sqrt(sum((glmFull$fitted.values-glmRes1$fitted.values)^2))/length(glmFull$fitted.values)
    #difGlobalMod2 <- sqrt(sum((glmFull$fitted.values-glmRes2$fitted.values)^2))/length(glmFull$fitted.values)
    
    # retObj <- list(c(R2Full=R2Full,R2Batch=R2Batch,R2Res1=R2Res1,R2Res2=R2Res2),
    #                c(difGlobalMod1=difGlobalMod1,difGlobalMod2=difGlobalMod2),
    #                tableCoefs,
    #                cleanedData=data.frame(cleanMod1=yRes1,cleanMod2=yRes2))
    retObj <- list(c(R2Full=R2Full,R2Batch=R2Batch,R2Res1=R2Res1),
                   c(difGlobalMod1=difGlobalMod1),
                   tableCoefs,
                   cleanedData=data.frame(cleanMod1=yRes1))
    
    return(retObj)  
}

get_residuals <- function(event_id, mdata){
    # Model one event at a time
    y <- pmin(pmax(as.numeric(psi[event_id,]),0),1) 
    glm_data <- cbind.data.frame(mdata, y)
    
    obj <- computeResiduals(glm_data, which(colnames(mdata) %in% covariates), which(colnames(mdata) %in% individual_traits))
    return(obj$cleanedData$cleanMod1)
    #return(obj$cleanedData$cleanMod2)
}

# hier.part functions ####
# ---- methods derived from hier.part
# ---- original methods renamed with ".mod" suffix
current.model.mod <- function (y, current.comb, xcan, SStot=0,family = c("gaussian","quasibinomial"), 
                               link = c("logit"), gof = c("Rsqu","RMSPE"), ...)  {
    comb.data <- data.frame(xcan[, current.comb])
    colnames(comb.data) <- colnames(xcan)[current.comb]
    data <- data.frame(y, comb.data)
    depv <- names(data)[1]
    n.comb <- dim(comb.data)[2]
    xs <- vector("character", n.comb)
    for (i in 1:(n.comb - 1)) xs[i] <- paste(names(comb.data)[i], 
                                             "+", sep = "")
    xs[n.comb] <- names(comb.data)[n.comb]
    xss <- paste(xs, collapse = " ", sep = "")
    formu <- stats::formula(paste(depv, "~", xss, sep = ""))
    if (gof == "RMSPE") gf <- sqrt(sum((stats::glm(formu, data = data, family = family,...)$fitted.values - y)^2))
    if (gof == "Rsqu") {
        if (family == "quasibinomial") 
            gf <- (SStot-sum((stats::glm(formu, data = data, family = family,...)$fitted.values - y)^2))/SStot
        if (family == "gaussian") 
            gf <- summary(stats::lm(formu, data = data))$r.squared
    }
    gf
}

all.regs.mod <- function (y, xcan, family = c("gaussian", "quasibinomial"), link = c("logit"), gof = c("Rsqu","RMSPE"),...) { 
    if (sum(is.na(xcan)) > 0) {
        missing <- is.na(apply(xcan, 1, FUN = sum))
        xcan <- xcan[!missing, ]
        y <- y[!missing]
        warning(paste(sum(missing), "observations deleted due to missingness in xcan\n"), 
                call. = FALSE)
    }
    if (sum(is.na(y)) > 0) {
        missing <- is.na(y)
        xcan <- xcan[!missing, ]
        y <- y[!missing]
        warning(paste(sum(missing), "observations deleted due to missingness in y\n"), 
                call. = FALSE)
    }
    pcan <- dim(xcan)[2]
    n <- (2^pcan) - 1
    combs <- combos1(pcan)$ragged
    SStot <- sum((y-mean(y))^2)
    
    if (gof == "RMSPE")  gfs <- sqrt(sum((stats::glm(y ~ 1, family = family,...)$fitted.values - y)^2))
    if (gof == "Rsqu")   gfs <- 0
    
    for (i in 1:n) {
        if (i%%500 == 0) 
            cat(i, "regressions calculated:", n - i, "to go...\n")
        current.comb <- as.vector(combs[i, ][combs[i, ] > 0]) 
        combn <- paste(names(data.frame(xcan)[current.comb]), "", collapse = "")
        if (gof == "RMSPE") new.line <- current.model.mod(y, current.comb, xcan,family = family, gof = "RMSPE",...)
        if (gof == "Rsqu")  new.line <- current.model.mod(y, current.comb, xcan,family = family, SStot=SStot,gof = "Rsqu",...)
        gfs <- c(gfs, new.line)
    }
    gfs
}

hier.part.mod <- function(y,xcan,family='gaussian',gof = "Rsqu", link = "",...) {
    pcan <- dim(xcan)[2]
    gfs <- all.regs.mod(y, xcan, family = family, gof = gof, link = link, ...)
    hp <- partition(gfs, pcan, var.names = names(data.frame(xcan)))
    
    params <- list(full.model = paste("y ~", paste(names(xcan),collapse = " + ")), 
                   family = family, link = link, gof = gof)
    if(sum(hp$IJ$I<0)>0){
        NA
    }else{
        list(gfs = gfs, IJ = hp$IJ, I.perc = hp$I.perc, params = params)
    }
}

# Fractional regression functions ####
# Function to obtain robust standard error for the coefficients
get_tables <- function(event_id, glm_model){
    df <-  as.data.frame(coef(summary(glm_model)))
    df <- df[contrast_names,]
    df$Trait <- individual_traits
    df$Ensembl_id <- unlist(strsplit(event_id, split = ";"))[[1]]
    df$Event_id <- event_id
    df$Event <- unlist(strsplit(event_id, split = ";"))[[2]]
    df$Type <- unlist(strsplit(df$Event, split = ":"))[[1]]
    df$`t value`<- NULL
    df$`Std. Error` = NULL
    df$glm.p.value = df$`Pr(>|t|)` 
    # ----- mcalvo use robust sandwich error #
    cft <- coeftest(glm_model, vcov.=vcovHC(glm_model, type="HC0"))
    df$Std.Error <- cft[contrast_names, 2]
    df$Z.Value <- cft[contrast_names, 3]
    df$P.Value <- cft[contrast_names, 4]
    df$Iter <- glm_model$iter
    
    # -------------------------------------- #
    # If warning, event might be modelled properly and function might return an NA (NaN)
    # Do this to keep track of warnigns. 
    # Note how events with warnings are associated with extremely large betas
    cft <- tryCatch({
        coeftest(glm_model, vcov.=vcovHC(glm_model, type="HC0"))
    }, warning = function(w) {
        return(NA)
        # "Warning message:
        #   In sqrt(diag(se)) : NaNs produced"
    }  
    )
    # Subset data used in downstream analysis
    if(!is.na(cft[1])){
        df$coeftest_Warning <- rep("0",length(contrast_names)) # 0 are events with no warnings
    }else{
        df$coeftest_Warning <- rep("1", length(contrast_names)) # 1 are events with  warnings
    }
    df <- df[,c("Event_id",
                "Trait",
                "Ensembl_id",
                "Type",  
                "Event",
                "Estimate",
                "glm.p.value",
                "Iter",
                "Std.Error",
                "Z.Value", 
                "P.Value",
                "coeftest_Warning")]
    return(df)
}

# Function to fit a glm model per event
model_psi <- function(event_id, mdata){
    print(event_id)
    # Model one event at a time
    psi_values <- pmin(pmax(as.numeric(psi[event_id,]),0),1) # forzar los límites de PSI a los valores teóricos de mínimo y máximo
    glm_data <- cbind.data.frame(mdata, psi_values)
    
    # Model formula: psi ~ covariates + traits
    mod_formula <- as.formula(paste("psi_values ~ ", paste(c(covariates, individual_traits), collapse = " + "), collapse = " ") )
    myglm <- glm(mod_formula, 
                 data = glm_data, 
                 family = quasibinomial('logit'),
                 control = list(maxit = 100) )
    
    # If there is a glm warning cause algorithm did not converge 
    my_warning <- tryCatch({
        glm(mod_formula, 
            data = glm_data, 
            family = quasibinomial('logit'),
            control = list( maxit = 100) )
    }, warning = function(w) {
        return(NA)
        # Warning message:
        #  glm.fit: algorithm did not converge 
    }  
    )
    
    # Create table
    event_results <- get_tables(event_id, myglm)
    
    if(!is.na(my_warning[1])){
        event_results$glm_Warning <- rep("0",length(contrast_names)) # 0 are events with no warnings
    }else{
        event_results$glm_Warning <- rep("1", length(contrast_names)) # 1 are events with  warnings
    }
    
    return(list('Event_id' = event_id,
                'glm' = myglm,
                'res' = event_results))  # table with estimates and P.Values per event
}
