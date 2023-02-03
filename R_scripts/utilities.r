## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 
##                     ----  Utilities ----
##
## --- A comparison of machine learning and statistical species distribution models: ---
##                -- when overfitting hurts interpretation -- 
## 
##                          --- December 2022 ---
##
## --- Emma Chollet, Andreas Scheidegger, Jonas Wydler and Nele Schuwirth ---
##
##                      --- emma.chollet@eawag.ch ---  
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ---- Small functions ----

# Capitalize first letter of a string
Cap <- function(x) {
  paste(toupper(substring(x, 1, 1)), substring(x, 2), sep = "", collapse = " ")
}

## ---- Data preprocessing ----

# Split data for CV and ODG
split.data <- function(data, list.taxa, dir.workspace, info.file.name.data, select.temp, CV, ODG, ODG.info = c(), bottom = T){
  
    info.file.name <- paste0(info.file.name.data, "_",
                             ifelse(CV, "CV_",
                                    ifelse(ODG, paste0(paste(c("ODG_", ODG.info["training.ratio"], ODG.info["variable"]), collapse = ""), "_"),
                                           "FIT_")),
                             "")
    
    if(CV){
      file.name <- paste0(dir.workspace, info.file.name, "Splits.rds")
    } else if(ODG){
      file.name <- paste0(dir.workspace, info.file.name, "SplitBy", ODG.info["model"] ,".rds")
    }
    
    if (file.exists(file.name)){
     splits <- readRDS(file = file.name)
     cat("\nFile with data splits already exists, we read it from", file.name, "and save it in object 'splits'.\n")
    } else {
      cat("\nNo data splits exist yet, we produce it and save it in", file.name, ".\n")
      if(CV){
        set.seed(2021)  
        
        folds <- groupKFold(data$SiteId, 3) # Keep same sites in same split to avoid data leakage
        
        train1 <- data[folds$Fold1,]
        test1 <- data[-folds$Fold1,]
        
        train2 <- data[folds$Fold2,]
        test2 <- data[-folds$Fold2,]
        
        train3 <- data[folds$Fold3,]
        test3 <- data[-folds$Fold3,]
        
        splits <- list("Split1" = list("Training data" = train1, "Testing data" = test1), 
                       "Split2" = list("Training data" = train2, "Testing data" = test2), 
                       "Split3" = list("Training data" = train3, "Testing data" = test3))
      } else if(ODG){
        training.ratio <- as.numeric(ODG.info["training.ratio"])
        variable <- ODG.info["variable"]
        model <- ODG.info["model"]

        if(variable == "random"){
          sample.size <- floor(training.ratio*nrow(data)) # size of training set is 80% of the whole dataset
          set.seed(2021) # used for reproducing the same output of simulations
          rind.train  <- sort(sample(nrow(data), size = sample.size))
          data.train  <- data[rind.train,]
          data.test   <- data[-rind.train,]
        } else {
          if(variable == "temperature"){
            variable <- paste(variable, model, sep = ".")
          }
          v <- data[,variable]
          q <- rank(v, na.last = F) / length(v) # put the NA in the testing set
          rind.train <- if(bottom == T){ 
            which(q < training.ratio)
          } else { 
            which(q > 1-training.ratio) 
          }
          data.train  <- data[rind.train,]
          data.test   <- data[-rind.train,]
        }
        
        splits <- list("Split1" = list("Training data" = data.train, "Testing data" = data.test))
        
      }
      saveRDS(splits, file = file.name)
    }
    return(splits)
}

# Standardize input environmental factors    
standardize.data <- function(data, splits, env.fact.full, dl, CV, ODG){
  
  # Add other temperatures to env fact list
  all.temp <- colnames(data)[which(grepl("temperature.lm", colnames(data)))]
  env.fact.full2 <- c(env.fact.full, all.temp)
  
  # Calculate mean and sd for env data for normalization with data leakage
  mean.dl <- apply(select(data, all_of(env.fact.full2)), 2, function(k){
      mean(k, na.rm = TRUE)
  })
  sd.dl <- apply(select(data, all_of(env.fact.full2)), 2, function(k){
      sd(k, na.rm = TRUE)
  })
  
  if(CV|ODG){
    
    temp.standardized.data <- lapply(splits, function(split){
      
      # Standardize training data
      # split <- splits[[1]]
      training.data <- split[[1]]
      
      mean.env.cond <- apply(select(training.data, all_of(env.fact.full2)), 2, function(k){
        mean(k, na.rm = TRUE)
      })
      sd.env.cond <- apply(select(training.data, all_of(env.fact.full2)), 2, function(k){
        sd(k, na.rm = TRUE)
      })
      
      # Standardize environmental factors
      if(dl){
        training.data[,env.fact.full2] <- as.data.frame(sweep(sweep(training.data[,env.fact.full2], 2, mean.dl, FUN="-"), 2, sd.dl, FUN = "/"))
      } else {
        training.data[,env.fact.full2] <- as.data.frame(sweep(sweep(training.data[,env.fact.full2], 2, mean.env.cond, FUN="-"), 2, sd.env.cond, FUN = "/"))
      }
      
      # Re-calculate temp2 and velocity2 with scaled variables
      training.data$temperature2 <- training.data$temperature^2
      training.data$velocity2 <- training.data$velocity^2
      
      # Standardize testing data
      testing.data <- split[[2]] # load testing data 

      if(dl){
        testing.data[,env.fact.full2] <- as.data.frame(sweep(sweep(testing.data[,env.fact.full2], 2, mean.dl, FUN="-"), 2, sd.dl, FUN = "/"))
      } else {
        testing.data[,env.fact.full2] <- as.data.frame(sweep(sweep(testing.data[,env.fact.full2], 2, mean.env.cond, FUN="-"), 2, sd.env.cond, FUN = "/"))
      }
      
      # Re-calculate temp2 and velocity2 with scaled variables
      testing.data$temperature2 <- testing.data$temperature^2
      testing.data$velocity2 <- testing.data$velocity^2
      
      return(list("Training data" = training.data, "Testing data" = testing.data, "Mean" = mean.env.cond, "SD" = sd.env.cond))
    })
    
    # Extract necessary information
    standardized.data <- lapply(temp.standardized.data,"[", 1:2) # only the splits without the mean, sd info
    normalization.data <- lapply(temp.standardized.data,"[", 3:4) # the mean and sd of the splits
    
  } else {
    
    training.data <- splits[[1]]
    training.data[,env.fact.full2] <- as.data.frame(sweep(sweep(training.data[,env.fact.full2], 2, mean.dl, FUN="-"), 2, sd.dl, FUN = "/"))
    
    # Re-calculate temp2 and velocity2 with scaled variables
    training.data$temperature2 <- training.data$temperature^2
    training.data$velocity2 <- training.data$velocity^2
    
    standardized.data <- list( "Entire dataset" = training.data)
    normalization.data <- list("Mean" = mean.dl, "SD" = sd.dl)
  }
  
  # Replace '0' and '1' by factors in standardized folds
  standardized.data.factors <- lapply(standardized.data, function(split){
    
    if(!CV & !ODG){
      cind.taxa <- which(grepl("Occurrence.",colnames(split)))
      for (i in cind.taxa ) {
        split[which(split[,i] == 0),i] <- "absent"
        split[which(split[,i] == 1),i] <- "present"
        split[,i] = as.factor(split[,i])
      }
      return(split)
      
    } else {
      
    return(lapply(split, function(fold){
      cind.taxa <- which(grepl("Occurrence.",colnames(fold)))
      for (i in cind.taxa ) {
        fold[which(fold[,i] == 0),i] <- "absent"
        fold[which(fold[,i] == 1),i] <- "present"
        fold[,i] = as.factor(fold[,i])
      }
      return(fold)
    }))
    }
    
  })
  
  stand.norm.data <- list("standardized.data" = standardized.data, "standardized.data.factors" = standardized.data.factors, 
                          "normalization.data" = normalization.data)
  return(stand.norm.data)
}

## ---- Models utilities ----

# Apply NULL model
apply.null.model <- function(data, list.taxa){
  
  # Make a list with the "outputs" of null model
  list.outputs <- vector(mode = 'list', length = length(list.taxa))
  names(list.outputs) <- list.taxa
  
  for (j in list.taxa){
    # j <- list.taxa[1]
    temp.output <- c("Likelihood", "Performance", "AUC")
    temp.list <- vector(mode = 'list', length = length(temp.output))
    names(temp.list) <- temp.output
    
    no.pres <- sum(data[, j] == 1, na.rm = TRUE)
    no.abs  <- sum(data[, j] == 0, na.rm = TRUE)
    no.obs  <- no.pres + no.abs
    prev <- no.pres / (no.obs)
    likeli <- rep(c(prev, 1-prev),c(no.pres, no.abs))
    temp.list[["Likelihood"]] <- likeli
    
    st.dev <- -2 * sum(log(likeli)) / no.obs
    temp.list[["Performance"]] <- st.dev
    
    auc <- suppressMessages(auc(data[,j], rep(prev, nrow(data))))
    temp.list[["AUC"]] <- auc
    
    list.outputs[[j]] <- temp.list
  }
  return(list.outputs)
}

## ---- Process output from stat models to fit structure of ml models (makes plotting easier)
transfrom.stat.outputs <- function(stat.outputs, list.taxa, CV, ODG){
    # CV = F
    if ( !CV & !ODG ){ # Model performance (FIT)
        stat.fit.res <- lapply(stat.outputs, function(models){
          #models <- stat.outputs[[1]]
          temp.list.st.dev <- vector(mode = "list", length = length(list.taxa))
          
          for(j in 1:length(list.taxa)){
            #j = 1
            dev.temp <- models[[2]]$deviance
            dev.temp <- dev.temp[c("Taxon", "std.deviance")]
            dev.temp <- subset(dev.temp, Taxon == list.taxa[[j]])
         
            dev.temp$Performance.train <- as.numeric(dev.temp$std.deviance)

            prop.temp <- models[[2]]$probability
            prop.temp <- subset(prop.temp, Taxon == list.taxa[[j]])
            
            prop.temp$Likelihood.train <- ifelse(prop.temp$Obs == 1, prop.temp$Pred, 1 - prop.temp$Pred)
            
            temp.list.st.dev[[j]] <- list("Trained model" = models[[1]],
              
                                          "Prediction factors training set" = ifelse(prop.temp$Pred >= 0.5,"present","absent"),
                                          "Prediction probabilities training set" = data.frame("present" = prop.temp$Pred, "absent" = 1 - prop.temp$Pred),
                                          "Likelihood training set" = prop.temp$Likelihood.train,
                                          "Performance training set" = dev.temp$Performance.train
            )
          }
          names(temp.list.st.dev) <-  list.taxa
          #stat.output.list[[n]] <- temp.list.st.dev
          return(temp.list.st.dev)
        })
        
        return(stat.fit.res) # this return statments are not really needed
        
    } else { # Prediction (CV or ODG)
        stat.cv.res <- lapply(stat.outputs, function(models){
            # models <- stat.outputs[[2]]
            
            stat.cv.res.splits <- lapply(models, function(split){
              # split <- models[[1]]
              temp.list.st.dev <- vector(mode = "list", length = length(list.taxa))
              
              for(j in 1:length(list.taxa)){
                # j = 1
                dev.temp <- split[[2]]$deviance
                dev.temp <- dev.temp[c("Taxon", "Type", "std.deviance")]
                dev.temp <- subset(dev.temp, Taxon == list.taxa[[j]])
                dev.temp.train <- subset(dev.temp, Type == "Training")
                dev.temp.test <- subset(dev.temp, Type == "Testing")
                
                
                dev.temp.train$Performance.train <- as.numeric(dev.temp.train$std.deviance)
                dev.temp.test$Performance.test <- as.numeric(dev.temp.test$std.deviance)

                prop.temp <- split[[2]]$probability
                prop.temp <- subset(prop.temp, Taxon == list.taxa[[j]])
                prop.temp.train <- subset(prop.temp, Type == "Training")
                prop.temp.test <- subset(prop.temp, Type == "Testing")
                
                prop.temp.train$Likelihood.train <- ifelse(prop.temp.train$Obs == 1, prop.temp.train$Pred, 1 - prop.temp.train$Pred)
                prop.temp.test$Likelihood.test <- ifelse(prop.temp.test$Obs == 1, prop.temp.test$Pred, 1 - prop.temp.test$Pred)
                
                
                temp.list.st.dev[[j]] <- list("Trained model" = split[[1]],

                                              "Prediction factors training set" = ifelse(prop.temp.train$Pred >= 0.5,"present","absent"),
                                              "Prediction probabilities training set" = data.frame("present" = prop.temp.train$Pred, "absent" = 1 - prop.temp.train$Pred),
                                              "Likelihood training set" = prop.temp.train$Likelihood.train,
                                              "Performance training set" = dev.temp.train$Performance.train,
                                              
                                              "Prediction factors testing set" = ifelse(prop.temp.test$Pred >= 0.5,"present","absent"),
                                              "Prediction probabilities testing set" = data.frame("present" = prop.temp.test$Pred,"absent" = 1 - prop.temp.test$Pred),
                                              "Likelihood testing set" = prop.temp.test$Likelihood.test,
                                              "Performance testing set" = dev.temp.test$Performance.test
                                              )
                
                }
                names(temp.list.st.dev) <-  list.taxa
                return(temp.list.st.dev)})
            return(stat.cv.res.splits)
            })
    
    return(stat.cv.res)
    }
}

# Function to use stat models results to apply to new matrix of predictors
pred.stat.models <- function(res.extracted, matrix.predictors){
  # res.extracted   <- rstan::extract(res,permuted=TRUE,inc_warmup=FALSE)
  
  # Extract inputs (x), observations (y), and parameters at maximum posterior
  ind.maxpost <- which.max(res.extracted[["lp__"]])
  mu.alpha.comm.maxpost <- res.extracted[["mu_alpha_comm"]][ind.maxpost]
  sigma.alpha.comm.maxpost <- res.extracted[["sigma_alpha_comm"]][ind.maxpost]
  mu.beta.comm.maxpost  <- res.extracted[["mu_beta_comm"]][ind.maxpost,]
  sigma.beta.comm.maxpost  <- res.extracted[["sigma_beta_comm"]][ind.maxpost,]
  alpha.taxa.maxpost <- res.extracted[["alpha_taxa"]][ind.maxpost,]
  beta.taxa.maxpost  <- res.extracted[["beta_taxa"]][ind.maxpost,,]
  
  # Calibration results
  # Check if site effects AND latent variables are disabled
  z <- matrix(rep(alpha.taxa.maxpost,nrow(matrix.predictors)),nrow=nrow(matrix.predictors),byrow=TRUE) + 
    matrix.predictors%*%beta.taxa.maxpost
  
  p.maxpost <- 1/(1+exp(-z))
  
  return(p.maxpost)
  
}

# make.final.outputs.cv <- function(outputs.cv, list.models, list.taxa){
#     
#     outputs <- vector(mode = "list", length = length(list.models))
#     names(outputs) <- list.models
#     
#     out <- c("Observation", #1
#              "Prediction factors", #2 
#              "Prediction probabilities", #3 
#              "Likelihood", #4
#              "Performance", #5
#              "Performance splits") #6
#     output.names <- paste(out, "testing set")
#     
#     for (l in list.models){
#         
#         temp.list.taxa <- vector(mode = "list", length = length(list.taxa))
#         names(temp.list.taxa) <- list.taxa
# 
#         for(j in list.taxa){
#             
#             temp.output <- vector(mode = "list", length = length(output.names))
#             names(temp.output) <- output.names
#             
#             for (m in output.names[3]) {
#                 temp.output[[m]] <- bind_rows(outputs.cv[[1]][[l]][[j]][[m]],
#                                               outputs.cv[[2]][[l]][[j]][[m]],
#                                               outputs.cv[[3]][[l]][[j]][[m]], .id = "Split")
#             }
#             for (m in output.names[c(2,4)]) {
#                 temp.output[[m]] <- c(outputs.cv[[1]][[l]][[j]][[m]],
#                                       outputs.cv[[2]][[l]][[j]][[m]],
#                                       outputs.cv[[3]][[l]][[j]][[m]])
#             }
#             
#             temp.vect <- vector(mode ="numeric", length = length(outputs.cv)) 
#             names(temp.vect) <- names(outputs.cv)
#             for (n in 1:length(outputs.cv)) {
#                 #n = 1
#                 perf <- outputs.cv[[n]][[l]][[j]][["Performance testing set"]]
#                 temp.vect[n] <- ifelse(is.numeric(perf), perf, NA)
#             }
#             temp.output[[5]] <- mean(temp.vect, na.rm = T)
#             temp.output[[6]] <- temp.vect
#             
#             temp.list.taxa[[j]] <- temp.output
#         }
#         outputs[[l]] <- temp.list.taxa
#     }
#     return(outputs)
# }

# Make dataframe with models outputs
make.df.outputs <- function(outputs, list.models, list.taxa, 
                            list.splits = c("Split1", "Split2", "Split3"), null.model, prev.inv, CV, ODG){
    
    if(CV | ODG){
        outputs.cv <- outputs
        no.splits <- length(list.splits)
        no.taxa <- length(list.taxa)
        no.models <- length(list.models)
        
        # make table with perf for each split
        df.pred.perf.cv <- data.frame(matrix(nrow = no.taxa, ncol = no.splits*no.models))
        colnames(df.pred.perf.cv) <- apply(expand.grid(list.splits,list.models), 1, paste, collapse="_")
        df.pred.perf.cv$Taxa <- list.taxa
        df.fit.perf.cv <- df.pred.perf.cv
        df.pred.auc.cv <- df.pred.perf.cv
        df.fit.auc.cv <- df.pred.perf.cv
        
        for (s in list.splits) {
          # s <- list.splits[1]
          # print(s)
            for (l in list.models) {
              # l <- list.models[2]
              # print(l)
                list.taxa.temp <- names(outputs.cv[[s]][[l]])
                for (j in list.taxa.temp) {
                  # j <- list.taxa[1]
                    rind.taxa <- which(df.pred.perf.cv$Taxa == j)
                    if(grepl("GLM", l)){
                      l.obs <- "iGLM"
                    } else {
                      l.obs <- l
                    }
                    fit.perf <- outputs.cv[[s]][[l]][[j]][["Performance training set"]]
                    fit.auc <- suppressMessages(as.numeric(auc(outputs.cv[[s]][[l.obs]][[j]][["Observation training set"]][,j], outputs.cv[[s]][[l]][[j]][["Prediction probabilities training set"]][,"present"])))
                    
                    pred.perf <- outputs.cv[[s]][[l]][[j]][["Performance testing set"]]
                    pred.auc <- suppressMessages(as.numeric(auc(outputs.cv[[s]][[l.obs]][[j]][["Observation testing set"]][,j], outputs.cv[[s]][[l]][[j]][["Prediction probabilities testing set"]][,"present"])))
                    
                    df.fit.perf.cv[rind.taxa,paste(s,l, sep="_")] <- ifelse(is.numeric(fit.perf), fit.perf, NA)
                    df.pred.perf.cv[rind.taxa,paste(s,l, sep="_")] <- ifelse(is.numeric(pred.perf), pred.perf, NA)
                    df.fit.auc.cv[rind.taxa,paste(s,l, sep="_")] <- ifelse(is.numeric(fit.auc), fit.auc, NA)
                    df.pred.auc.cv[rind.taxa,paste(s,l, sep="_")] <- ifelse(is.numeric(pred.auc), pred.auc, NA)
                    
                }
            }
        }
    }
    
    # make table with mean perf across splits
    df.pred.perf <- data.frame(matrix(nrow = no.taxa, ncol = no.models))
    colnames(df.pred.perf) <- list.models
    df.pred.perf$Taxa <- list.taxa
    df.pred.perf$Prevalence <- prev.inv[which(prev.inv$Occurrence.taxa %in% list.taxa), "Prevalence"]
    df.pred.perf[, "Taxonomic level"] <- prev.inv[which(prev.inv$Occurrence.taxa %in% list.taxa), "Taxonomic.level"]
    
    df.pred.auc <- df.pred.perf
    df.fit.auc <- df.pred.perf
    
    # Add columns for explanatory power
    df.pred.perf[,"Null_model"] <- NA
    expl.pow <- paste0("expl.pow_", list.models)
    df.pred.perf[, expl.pow] <- NA
    df.fit.perf <- df.pred.perf
    
    for (j in list.taxa) {
      # print(j)
      # j <- list.taxa[1]
        rind.taxa <- which(df.pred.perf$Taxa == j)
        df.pred.perf[rind.taxa, "Null_model"] <- null.model[[j]][["Performance"]]
        df.fit.perf[rind.taxa, "Null_model"] <- null.model[[j]][["Performance"]]
        df.pred.auc[rind.taxa, "Null_model.auc"] <- null.model[[j]][["AUC"]]
        df.fit.auc[rind.taxa, "Null_model.auc"] <- null.model[[j]][["AUC"]]
        for(l in list.models){
          # print(l)
          # l <- list.models[1]
            if(CV | ODG){
                splits.model <- apply(expand.grid(list.splits,l), 1, paste, collapse="_")
                # For testing/prediction
                # standardized deviance
                mean.dev <- mean(as.matrix(df.pred.perf.cv[rind.taxa, splits.model]), na.rm = T)
                df.pred.perf[rind.taxa,l] <- mean.dev
                val.expl.pow <- (null.model[[j]][["Performance"]] - mean.dev) / null.model[[j]][["Performance"]]
                df.pred.perf[rind.taxa, paste0("expl.pow_",l)] <- val.expl.pow
                # auc
                mean.auc <- mean(as.matrix(df.pred.auc.cv[rind.taxa, splits.model]), na.rm = T)
                df.pred.auc[rind.taxa,l] <- mean.auc

                # For training/fitting
                # standardized deviance
                mean.dev <- mean(as.matrix(df.fit.perf.cv[rind.taxa, splits.model]), na.rm = T)
                df.fit.perf[rind.taxa,l] <- mean.dev
                val.expl.pow <- (null.model[[j]][["Performance"]] - mean.dev) / null.model[[j]][["Performance"]]
                df.fit.perf[rind.taxa, paste0("expl.pow_",l)] <- val.expl.pow
                # auc
                mean.auc <- mean(as.matrix(df.fit.auc.cv[rind.taxa, splits.model]), na.rm = T)
                df.fit.auc[rind.taxa,l] <- mean.auc
            } else {
              # standardized deviance
              dev <- outputs[[l]][[j]][["Performance training set"]]
              df.fit.perf[rind.taxa,l] <- dev
              val.expl.pow <- (null.model[[j]][["Performance"]] - dev) / null.model[[j]][["Performance"]]
              df.fit.perf[rind.taxa, paste0("expl.pow_",l)] <- val.expl.pow
              
              # auc
              if(grepl("GLM", l)){
                l.obs <- "iGLM"
              } else {
                l.obs <- l
              }
              suppressMessages(
                if(j == "Occurrence.Perlodidae" & l == "SVM"){ # convergence problem for SVM applied to Perlodidae
                  auc <- NA
                } else {
                  auc <- as.numeric(auc(outputs[[l.obs]][[j]][["Observation training set"]][,j], outputs[[l]][[j]][["Prediction probabilities training set"]][,"present"]))
                }
              )
              df.fit.auc[rind.taxa,l] <- auc
            }

        }
    }
    
    # Make one df with performance during training AND testing AND expl.pow
    common.vect <- c("Taxa", "Prevalence", "Taxonomic level", "Null_model")
    common.vect2 <- c("Taxa", "Prevalence", "Taxonomic level", "Null_model.auc")
    
    if(CV | ODG){
      # Merge dataframes for comparison
      merged.dev.df <- left_join(df.fit.perf[,unique(c(common.vect, list.models, all_of(expl.pow)))], df.pred.perf[,unique(c(common.vect, list.models, all_of(expl.pow)))], 
                             by = c(common.vect), suffix = c(".dev.fit", ".dev.pred"))
      merged.auc.df <- left_join(df.fit.auc, df.pred.auc, 
                                 by = c(common.vect2), suffix = c(".auc.fit", ".auc.pred"))
      # colnames(merged.auc.df)[which(colnames(merged.auc.df) == "Null_model")] <- "Null_model.auc"
      merged.df <- left_join(merged.dev.df, merged.auc.df, by = common.vect[-which(common.vect == "Null_model")])
      merged.df[, paste0(list.models, ".likelihood.ratio")] <- NA
      merged.df[, paste0(list.models, ".auc.ratio")] <- NA
      merged.df$Big.model.diff <- NA
      merged.df$Big.pred.expl.pow.diff <- NA
      if(any(grepl("chGLM", list.models) == TRUE)){
        merged.df$chGLM.model.diff <- NA
        merged.df$chGLM.pred.expl.pow.diff <- NA
      }
      
      # Compute biggest difference in expl. pow. of prediction
      temp.df <- select(merged.df, "Taxa", contains("expl.pow") & contains(".pred") & !contains("diff"))
      for(j in list.taxa){
        # j <- list.taxa[1]
        # Find min and max expl. pow
        min <- min(temp.df[which(temp.df$Taxa == j), 2:(no.models+1)])
        max <- max(temp.df[which(temp.df$Taxa == j), 2:(no.models+1)])
        model.min <- sub("expl.pow_", "", sub(".pred", "", colnames(temp.df)[which(temp.df[which(temp.df$Taxa == j), ] == min)]))
        model.max <- sub("expl.pow_", "", sub(".pred", "", colnames(temp.df)[which(temp.df[which(temp.df$Taxa == j), ] == max)]))
        merged.df[which(merged.df$Taxa == j), "Big.model.diff"] <- paste(model.max, model.min, sep = "-")
        merged.df[which(merged.df$Taxa == j), "Big.pred.expl.pow.diff"] <- max - min
        
        # Compare with chGLM
        if(any(grepl("chGLM", colnames(temp.df)) == TRUE)){
          expl.pow.chGLM <- temp.df[which(temp.df$Taxa == j), which(grepl("chGLM", colnames(temp.df)))]
          merged.df[which(merged.df$Taxa == j), "chGLM.model.diff"] <- paste(model.max, "chGLM", sep = "-")
          merged.df[which(merged.df$Taxa == j), "chGLM.pred.expl.pow.diff"] <- max - expl.pow.chGLM
        }
        
        # Compute ratio
        for (l in list.models) {
          # likelihood ratio
          pred <- merged.df[which(merged.df$Taxa == j), paste0(l, ".dev.pred")]
          fit <- merged.df[which(merged.df$Taxa == j), paste0(l, ".dev.fit")]
          merged.df[which(merged.df$Taxa == j), paste0(l, ".likelihood.ratio")] <- exp(-(pred - fit) / 2)
          # auc ratio
          auc.pred <- merged.df[which(merged.df$Taxa == j), paste0(l, ".auc.pred")]
          auc.fit <- merged.df[which(merged.df$Taxa == j), paste0(l, ".auc.fit")]
          merged.df[which(merged.df$Taxa == j), paste0(l, ".auc.ratio")] <- auc.pred / auc.fit
        }
      }
    } else {
      merged.df <- left_join(df.fit.perf, df.fit.auc, 
                             by = c(common.vect[-which(common.vect == "Null_model")]), suffix = c(".dev.pred", ".auc.pred"))
    }
    
    if(CV | ODG){
        result <- list("Table predictive performance CV" = df.pred.perf.cv, "Table predictive performance" = df.pred.perf,
                       "Table fit performance CV" = df.fit.perf.cv, "Table fit performance" = df.fit.perf, "Table merged" = merged.df)
    } else {
        result <- merged.df
    }
    
    return(result)
}

# Make dataframe of models performance results for plotting
perf.plot.data <- function(list.df.perf){
  
  names.appcase <- names(list.df.perf)
  plot.data <- data.frame()
  
  for(n in 1:length(list.df.perf)){  
    # n = 1
    
    plot.data0 <- list.df.perf[[n]]
    col.dev.fit <- c(colnames(plot.data0)[grepl(".dev.fit", colnames(plot.data0))], "Null_model")
    col.dev.pred <- c(colnames(plot.data0)[grepl(".dev.pred", colnames(plot.data0))], "Null_model")
    col.auc.fit <- c(colnames(plot.data0)[grepl(".auc.fit", colnames(plot.data0))], "Null_model")
    col.auc.pred <- c(colnames(plot.data0)[grepl(".auc.pred", colnames(plot.data0))], "Null_model")
    
    
    plot.data1 <- gather(plot.data0, key = model, value = performance.fit, all_of(col.dev.fit))
    plot.data2 <- gather(plot.data0, key = model, value = performance.pred, all_of(col.dev.pred))
    plot.data3 <- gather(plot.data0, key = model, value = auc.fit, all_of(col.auc.fit))
    plot.data4 <- gather(plot.data0, key = model, value = auc.pred, all_of(col.auc.pred))
    
    plot.data1$model <- sub(".dev.fit", "", plot.data1$model)
    plot.data1 <- plot.data1[,c("Taxa", "Prevalence", # "Taxonomic level", 
                                "model", "performance.fit")]
    plot.data2$model <- sub(".dev.pred", "", plot.data2$model)
    plot.data2 <- plot.data2[,c("Taxa", "Prevalence", # "Taxonomic level", 
                                "model", "performance.pred")]
    plot.data3$model <- sub(".auc.fit", "", plot.data3$model)
    plot.data3 <- plot.data3[,c("Taxa", "Prevalence", # "Taxonomic level", 
                                "model", "auc.fit")]
    plot.data4$model <- sub(".auc.pred", "", plot.data4$model)
    plot.data4 <- plot.data4[,c("Taxa", "Prevalence", # "Taxonomic level", 
                                "model", "auc.pred")]
    
    if("FIT" %in% names.appcase ){
      plot.data0 <- left_join(plot.data1, plot.data3)
    } else {
      plot.data0 <- left_join(plot.data1, plot.data2, by = c("Taxa", "Prevalence", # "Taxonomic level", 
                                                             "model"))
      plot.data0 <- left_join(plot.data0, plot.data3, by = c("Taxa", "Prevalence", # "Taxonomic level", 
                                                             "model"))
      plot.data0 <- left_join(plot.data0, plot.data4, by = c("Taxa", "Prevalence", # "Taxonomic level", 
                                                             "model"))
    }
      
    plot.data5 <- gather(plot.data0, key = dataset, value = performance, -c("Taxa", "Prevalence", # "Taxonomic level", 
                                                                            "model") )
    plot.data5$dataset <- sub("performance.fit", "Calibration", plot.data5$dataset)
    plot.data5$dataset <- sub("performance.pred", "Prediction", plot.data5$dataset)
    plot.data5$dataset <- sub("auc.fit", "AUC Calibration", plot.data5$dataset)
    plot.data5$dataset <- sub("auc.pred", "AUC Prediction", plot.data5$dataset)
    plot.data5$appcase <- names.appcase[n]
    
    # list.plot.data[[n]] <- plot.data3
    plot.data <- bind_rows(plot.data, plot.data5)
  }
  
  return(plot.data)
}

# Make list of models prediction during CV for each taxon
make.list.models.pred <- function(outputs.cv, list.taxa, list.models){
  
  list.df.models.pred <- vector(mode = "list", length = length(list.taxa))
  names(list.df.models.pred) <- list.taxa
  
  for(j in list.taxa){
    # j <- list.taxa[1]
    df.models.pred <- data.frame()
    for (l in list.models) {
      # l <- list.models[1]
      if(l == "hGLM" | l == "chGLM"){
        temp.l <- "iGLM"
      } else {
        temp.l <- l
      }
      list.split.df.pred <- list(bind_cols(outputs.cv[["Split1"]][[l]][[j]][["Prediction probabilities testing set"]],
                                           outputs.cv[["Split1"]][[temp.l]][[j]][["Observation testing set"]][,c("X","Y", j)]),
                                 bind_cols(outputs.cv[["Split2"]][[l]][[j]][["Prediction probabilities testing set"]],
                                           outputs.cv[["Split2"]][[temp.l]][[j]][["Observation testing set"]][,c("X","Y", j)]),
                                 bind_cols(outputs.cv[["Split3"]][[l]][[j]][["Prediction probabilities testing set"]],
                                           outputs.cv[["Split3"]][[temp.l]][[j]][["Observation testing set"]][,c("X","Y", j)]))
      temp.df.pred <- bind_rows(list.split.df.pred, .id = "Split")
      temp.df.pred$Model <- l
      df.models.pred <- bind_rows(df.models.pred, temp.df.pred)
    }
    df.models.pred$Observation <- df.models.pred[,j]
    list.df.models.pred[[j]] <- df.models.pred
  }
  return(list.df.models.pred)
}

# Make data frame for plotting ICE
ice.plot.data <- function(taxon, outputs, ref.data, list.models, list.taxa, env.fact, select.env.fact, normalization.data, ODG, no.samples, no.steps, vect.seeds = 2021){
  
  # ref.data = standardized.data[[1]]
  # taxon <- select.taxa[1]
  name.taxon <- gsub("Occurrence.", "", taxon)
  cat("\nProducing dataframes for ICE for taxa", taxon)
  
  no.models <- length(list.models)
  # no.subselect <- length(subselect)
  
  if(ODG){
    normalization.data <- normalization.data[[1]]
    
    standardized.training.data <- ref.data$`Training data`
    standardized.training.data$dataset <- "Calibration"
    
    standardized.testing.data <- ref.data$`Testing data`
    standardized.testing.data$dataset <- "Prediction"
    
    ref.data <- bind_rows(standardized.training.data, standardized.testing.data)
  }
  
  # Initialize final list with all arranged plots
  list.ice.plot.data <- vector(mode = "list", length = length(select.env.fact))
  names(list.ice.plot.data) <- names(select.env.fact)
  
  for(k in select.env.fact){
    # taxon <- subselect.taxa[1]
    # k <- select.env.fact[1]
    # seed <- 2021
    name.fact <- names(select.env.fact)[which(select.env.fact == k)]
    cat("\nFor env. fact.", name.fact)
    
    plot.data <- data.frame()
    plot.data.means <- data.frame()
    plot.data.mean.pdp <- data.frame()
    plot.data.min.max.pdp <- data.frame()
    plot.data.rug <- data.frame()
    
    for(l in list.models){
      # l <- list.models[8]
      cat("\nfor model", l)
      
      # Extract trained model
      # trained.mod <- ann.outputs.cv$Split1$ANN$Occurrence.Gammaridae$`Trained model`
      trained.mod <- outputs[[l]][[taxon]][["Trained model"]]
      # trained.mod <- stat.outputs.transformed$hGLM$Occurrence.Gammaridae$`Trained model`
      
      # Prepare environmental data
      if(grepl("GLM", l) & !("temperature2" %in% env.fact)){
        vect.env.fact <- c(env.fact, "temperature2", "velocity2")
      } else if(!grepl("GLM", l)) {
        vect.env.fact <- env.fact
      }
      
      m <- min(ref.data[,k])
      M <- max(ref.data[,k])
      range.test <- seq(m, M, length.out = no.steps)
      
      if(ODG){
        ref.env.df <- ref.data[which(ref.data$dataset == "Calibration"), vect.env.fact]
      } else {
        ref.env.df <- ref.data[ , vect.env.fact]
      }
      
      # Make range of backward normalized values for labelling x axis
      # !! Might be mathematically false 
      m2 <- (m * normalization.data$SD[name.fact]) + normalization.data$Mean[name.fact]
      M2 <- (M * normalization.data$SD[name.fact]) + normalization.data$Mean[name.fact]
      range.orig.fact <- round(seq(m2, M2, length.out = no.steps), digits = 4)
      
      for(s in vect.seeds){
        # Subselect a number of samples/observations to plot
        # no.samples <- 3
        # s <- vect.seeds[2]
        set.seed(s)
        cat(" with seed", s)
        
        env.df <- ref.env.df[sample(nrow(ref.env.df), size = no.samples),]
        # env.df <- env.df[1,]
        
        # Make a dataframe for predicted values for each sample
        pred.df <- data.frame(matrix(nrow = no.samples, ncol = no.steps))
        colnames(pred.df) <- range.orig.fact
        
        for(n in 1:no.samples){
          # n <- 1
          env.fact.test <- env.df[n,]
          env.fact.test <- env.fact.test[rep(seq_len(1), each = no.steps), ]
          env.fact.test[,k] <- range.test
          
          # Recalculate temp2 or velocity2 for whole range
          if( (k == "temperature" | k == "velocity") & grepl("GLM", l)){
            env.fact.test[,paste0(k,2)] <- env.fact.test[,k]^2
            # env.fact.test[13,paste0(k,2)] <- env.fact.test[13,k]^2
          }
          
          # Make predictions
          if(l == "ANN"){
            env.fact.test <- as.matrix(env.fact.test)
            pred.df[n,] <- predict(trained.mod, env.fact.test)[ , which(names(outputs[[l]]) == taxon)]
          } else if (l == "hGLM" | l == "chGLM"){
            res.extracted   <- rstan::extract(trained.mod,permuted=TRUE,inc_warmup=FALSE)
            pred.stat <- pred.stat.models(res.extracted = res.extracted, matrix.predictors = as.matrix(env.fact.test))
            colnames(pred.stat) <- list.taxa
            pred.df[n,] <- pred.stat[,taxon]
          } else {
            pred.df[n,] <- predict(trained.mod, env.fact.test, type = 'prob')[,"present"]
          }
        } 
        
        # Subselect predictions for plot resolution
        # n <- subselect
        # temp.range.fact <- range.orig.fact[seq(1,200, by = subselect[n])]
        temp.pred.df <- pred.df
        temp.pred.df$observation <- rownames(pred.df)
        
        # Prepare plot data
        temp.plot.data <- gather(temp.pred.df, key = variable, value = value, -observation)
        temp.plot.data$observation <- as.factor(temp.plot.data$observation)
        temp.plot.data$Model <- l
        temp.plot.data[, "variable"] <- as.numeric(temp.plot.data[, "variable"])
        temp.plot.data$Seed <- s
        
        plot.data <- bind_rows(plot.data, temp.plot.data)
        
        # Make dataframe for PDP (mean of ICE)
        mean.pdp <- data.frame(colMeans(temp.pred.df[,1:(length(temp.pred.df)-1)]))
        colnames(mean.pdp) <- "mean.pdp"
        mean.pdp$variable <- range.orig.fact
        mean.pdp$Model <- as.factor(l)
        mean.pdp$Seed <- s
        plot.data.mean.pdp <- bind_rows(plot.data.mean.pdp, mean.pdp)
        
        # Make dataframe for min and max value of pdp
        temp.min.max <- data.frame("Model" = l, "min" = min(mean.pdp$mean), "max" = max(mean.pdp$mean))
        temp.min.max$Seed <- s
        plot.data.min.max.pdp <- bind_rows(plot.data.min.max.pdp, temp.min.max)
        
        
        # Make dataframe for predicted values when others are at their mean
        means <- data.frame(variable = range.orig.fact)
        means$Model <- l
        means$Seed <- s
        env.fact.mean.test <- data.frame(matrix(rep(0, times = length(vect.env.fact)*no.steps), nrow = no.steps))
        colnames(env.fact.mean.test) <- vect.env.fact
        env.fact.mean.test[,k] <- range.test
        if( (k == "temperature" | k == "velocity") & grepl("GLM", l)){
          env.fact.mean.test[,paste0(k,2)] <- env.fact.mean.test[,k]^2
          # env.fact.test[13,paste0(k,2)] <- env.fact.test[13,k]^2
        }
        if(l == "ANN"){
          env.fact.mean.test <- as.matrix(env.fact.mean.test)
          means$mean.average <- predict(trained.mod, env.fact.mean.test)[ , which(names(outputs[[l]]) == taxon)]
        } else if (l == "hGLM" | l == "chGLM"){
          res.extracted   <- rstan::extract(trained.mod,permuted=TRUE,inc_warmup=FALSE)
          pred.stat <- pred.stat.models(res.extracted = res.extracted, matrix.predictors = as.matrix(env.fact.mean.test))
          colnames(pred.stat) <- list.taxa
          means$mean.average <- pred.stat[,taxon]
        } else {
          means$mean.average <- predict(trained.mod, env.fact.mean.test, type = 'prob')[,"present"]
        }
        plot.data.means <- bind_rows(plot.data.means, means)
        
        # Make dataframe for rug (datapoints of observations)
        observations <- as.data.frame((env.df[,k] * normalization.data$SD[name.fact]) + normalization.data$Mean[name.fact])
        colnames(observations) <- "variable"
        observations$Model <- as.factor(l)
        observations$Seed <- s
        plot.data.rug <- bind_rows(plot.data.rug, observations)
      }
    }
    
    list.ice.plot.data[[name.fact]] <- list("name.taxon" = name.taxon, "env.fact" = name.fact, "plot.data.ice" = plot.data, 
                                            "plot.data.mean.pdp" = plot.data.mean.pdp, "plot.data.min.max.pdp" = plot.data.min.max.pdp,
                                            "plot.data.means" = plot.data.means, "plot.data.rug" = plot.data.rug)
  }
  return(list.ice.plot.data)
}

# Table 3: Summary statitistics likelihood ratio ####
table.likelihood.ratio <- function(list.df.perf, list.models, models.analysis){
  
  # [min, mean, median, max]x[list.models]
  df.likeli.ratio <- data.frame()
  for (n in 1:length(list.df.perf)) {
    # n = 2
    temp.df <- data.frame("Case" = rep(names(list.df.perf)[n], 4))
    temp.df$Stat <- c("min", "mean", "median", "max")
    for (l in list.models) {
      # l = list.models[1]
      vect.likeli.ratio <- list.df.perf[[n]][,paste0(l, ".likelihood.ratio")]
      temp.df[,l] <- round(c(min(vect.likeli.ratio),
                           mean(vect.likeli.ratio),
                           median(vect.likeli.ratio),
                           max(vect.likeli.ratio)), digits = 2)
    }
    if(any(models.analysis == T) & n == 1){
      temp.df <- temp.df[, order(temp.df[which(temp.df$Stat == "median"), ], decreasing = T)] # order columns by decreasing order
    }
  df.likeli.ratio <- bind_rows(df.likeli.ratio, temp.df)
  }
  return(df.likeli.ratio)
}