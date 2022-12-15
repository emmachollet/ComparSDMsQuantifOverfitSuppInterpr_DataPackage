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

preprocess.data <- function(data.env, data.inv, prev.inv, taxa.ibch = NA, remove.na.env.fact = T, remove.na.taxa = T, prev.restrict, env.fact.full, vect.info, dir.workspace, BDM){
    
  # Print information
  cind.taxa <- which(grepl("Occurrence.", colnames(data.inv)))
  list.taxa <- colnames(data.inv)[cind.taxa]
  cat("\nSummary information of original datasets before preprocessing:\n",
      length(unique(data.inv$SampId)), "samples,\n",
      length(unique(data.inv$SiteId)), "sites,\n",
      length(list.taxa), "taxa (including missing values) with taxonomic level:\n")
  
  print(summary(as.factor(prev.inv[which(prev.inv$Occurrence.taxa %in% list.taxa), "Taxonomic.level"])))
  
  if(remove.na.env.fact){
    # Drop rows with incomplete influence factors
    rind <- !apply(is.na(data.env[,env.fact.full]), 1, FUN=any)
    rind <- ifelse(is.na(rind), FALSE, rind)
    cat(paste("\nExcluding samples because of incomplete environmental factor:", sum(!rind),"\n"))
    data.env <- data.env[rind,]
    dim(data.env)
  }
  data.inv <- data.inv[data.inv$SampId %in% data.env$SampId, ]
  
  if(remove.na.taxa){
    # Drop taxa with too many NAs
    too.many.na <- c()
    threshold <- ifelse(BDM, 100, 200)
    for(i in cind.taxa){
      if(sum(is.na(data.inv[,i])) > threshold){ too.many.na <- c(too.many.na, i)}
    }
    cat("\nExcluding following", length(too.many.na), "taxa because they have more than", 
        threshold, "NAs:\n", gsub("Occurrence.", "", colnames(data.inv)[too.many.na]), "\n")
    if(!is.null(too.many.na)){
      data.inv <- data.inv[, -too.many.na]
    }
    
    # Drop rows with incomplete taxa
    cind.taxa <- which(grepl("Occurrence.", colnames(data.inv)))
    rind <- !apply(is.na(data.inv[,cind.taxa]), 1, FUN=any)
    rind <- ifelse(is.na(rind), FALSE, rind)
    cat(paste("\nExcluding samples because of incomplete taxa:", sum(!rind),"\n"))
    data.inv <- data.inv[rind,]
    dim(data.inv)
    
    # Delete taxa with only presence or absence (happens because we deleted rows)
    # Delete taxa with less than 5% of only present or absent
    cind.taxa <- which(grepl("Occurrence.",colnames(data.inv)))
    no.samples <- nrow(data.inv)
    prev.perc <- prev.restrict*no.samples
    cind.rem <- c()
    for (j in cind.taxa) {
        if(sum(data.inv[,j]) < prev.perc | sum(data.inv[,j]) > no.samples - prev.perc){ cind.rem <- c(cind.rem, j) }
    }
    cat("\nExcluding",  length(cind.rem),"taxa because they have less than", 100*prev.restrict, "% or more than", 
        100 - 100*prev.restrict, "% of prevalence:\n", gsub("Occurrence.", "", colnames(data.inv)[cind.rem]), "\n")
    if(length(cind.rem) != 0){ 
      data.inv <- data.inv[,-cind.rem] 
    }
    dim(data.inv)
  }
  
  # Update prevalence dataframe with new list of taxa and number of samples
  cind.taxa <- which(grepl("Occurrence.", colnames(data.inv))) # update list.taxa
  list.taxa <- colnames(data.inv)[cind.taxa]
  prev.inv <- prev.inv[which(prev.inv$Occurrence.taxa %in% list.taxa),]
  for (j in list.taxa) { 
    prev <- sum(data.inv[,j]) / dim(data.inv)[1]
    prev.inv[which(prev.inv$Occurrence.taxa == j), "Prevalence"] <- prev
  }
  prev.inv[which(grepl("Occurrence.group", prev.inv$Occurrence.taxa)), "Taxonomic.level"] <- "group"
  
  # Reorder taxa by prevalence
  list.taxa <- list.taxa[order(match(list.taxa, prev.inv$Occurrence.taxa))]
  
  # Add other temperatures to env fact list
  all.temp <- colnames(data.env)[which(grepl("temperature.lm", colnames(data.env)))]
  env.fact.full2 <- c(env.fact.full, all.temp)
  
  # Merge data sets
  data <- data.env[, c("SiteId", "SampId", "X", "Y", vect.info, env.fact.full2)] %>%
    right_join(data.inv[, c(which(colnames(data.inv) %in% c("SiteId", "SampId", list.taxa)))], by = c("SiteId", "SampId"))
  # Reorder taxa columns by prevalence
  data <- data[,c(which(!(colnames(data) %in% list.taxa)), match(list.taxa, colnames(data)))]
  dim(data)
  
  suppressWarnings(
  if(!is.na(taxa.ibch)){
  
    # Analyze which taxa needed to calculate the IBCH index
    # Clean taxa names in our list
    temp.list.taxa <- gsub("Occurrence.", "", list.taxa)
    temp.list.taxa <- gsub("group.", "", temp.list.taxa)
    temp.list.taxa <- sort(temp.list.taxa)
    
    # Clean taxa names in IBCH list
    list.taxa.ibch <- taxa.ibch[,which(grepl("IBCH_2019", colnames(taxa.ibch)))]
    list.taxa.ibch <- gsub("Tachet", "", list.taxa.ibch)
    list.taxa.ibch <- gsub(" ", "", list.taxa.ibch)
    list.taxa.ibch <- unlist(sapply(list.taxa.ibch, FUN = strsplit, split = "[^[:alnum:]]", USE.NAMES = F)) # split taxa name with "\"
    list.taxa.ibch <- tolower(list.taxa.ibch) # all character to lower cap
    list.taxa.ibch <- sapply(list.taxa.ibch, FUN =  Cap, USE.NAMES = F) # capitalize first letter
    list.taxa.ibch <- sort(list.taxa.ibch)
    
    ind <- which(list.taxa.ibch %in% temp.list.taxa)
    ind2 <- which(temp.list.taxa %in% list.taxa.ibch)
    cat("\nFollowing", length(ind), "taxa needed for IBCH are present in our dataset:\n",
        list.taxa.ibch[ind],
        "\n\nFollowing", length(list.taxa.ibch) - length(ind),"taxa needed for IBCH are missing in our dataset:\n",
        list.taxa.ibch[-ind],
        "\n\nFollowing", length(temp.list.taxa) - length(ind2), "taxa not needed for IBCH are present in our dataset:\n",
        temp.list.taxa[-ind2], "\n"
    )
  }
  )
  
  # Print information
  cat("\nSummary information of final datasets after preprocessing:\n",
      length(unique(data$SampId)), "samples,\n",
      length(unique(data$SiteId)), "sites,\n",
      length(list.taxa), "taxa", ifelse(remove.na.taxa, "(without missing values)",""), "with prevalence between", 
      100*prev.restrict, "% and", 
      100 - 100*prev.restrict, ", with taxonomic level:\n")
  print(summary(as.factor(prev.inv[which(prev.inv$Occurrence.taxa %in% list.taxa), "Taxonomic.level"])))
  
  preprocessed.data <- list("data" = data, "list.taxa" = list.taxa, "prev.inv" = prev.inv)
  
  return(preprocessed.data)
}

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
    # stat.outputs = stat.outputs
    
    # stat.output.list <- vector(mode = "list", length = length(stat.outputs))
    
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

make.final.outputs.cv <- function(outputs.cv, list.models, list.taxa){
    
    outputs <- vector(mode = "list", length = length(list.models))
    names(outputs) <- list.models
    
    out <- c("Observation", #1
             "Prediction factors", #2 
             "Prediction probabilities", #3 
             "Likelihood", #4
             "Performance", #5
             "Performance splits") #6
    output.names <- paste(out, "testing set")
    
    for (l in list.models){
        
        temp.list.taxa <- vector(mode = "list", length = length(list.taxa))
        names(temp.list.taxa) <- list.taxa

        for(j in list.taxa){
            
            temp.output <- vector(mode = "list", length = length(output.names))
            names(temp.output) <- output.names
            
            for (m in output.names[3]) {
                temp.output[[m]] <- bind_rows(outputs.cv[[1]][[l]][[j]][[m]],
                                              outputs.cv[[2]][[l]][[j]][[m]],
                                              outputs.cv[[3]][[l]][[j]][[m]], .id = "Split")
            }
            for (m in output.names[c(2,4)]) {
                temp.output[[m]] <- c(outputs.cv[[1]][[l]][[j]][[m]],
                                      outputs.cv[[2]][[l]][[j]][[m]],
                                      outputs.cv[[3]][[l]][[j]][[m]])
            }
            
            temp.vect <- vector(mode ="numeric", length = length(outputs.cv)) 
            names(temp.vect) <- names(outputs.cv)
            for (n in 1:length(outputs.cv)) {
                #n = 1
                perf <- outputs.cv[[n]][[l]][[j]][["Performance testing set"]]
                temp.vect[n] <- ifelse(is.numeric(perf), perf, NA)
            }
            temp.output[[5]] <- mean(temp.vect, na.rm = T)
            temp.output[[6]] <- temp.vect
            
            temp.list.taxa[[j]] <- temp.output
        }
        outputs[[l]] <- temp.list.taxa
    }
    return(outputs)
}

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

make.table <- function(df.pred.perf, df.fit.perf, list.models){
  list.models <- append(list.models, "Null_model")
  names(list.models) <- c()

  # calculate mean standardized deviance for the training set (i.e. quality of fit)
  table.fit.mean <- apply(df.fit.perf[list.models],2, FUN = mean)
  table.fit.mean <-t(table.fit.mean)
  table.fit.mean <- as.data.frame(table.fit.mean)
  rownames(table.fit.mean) <- "Mean std. dev. during training"
  
  # calculate corresponding standart deviation
  table.fit.sd <- apply(df.fit.perf[list.models],2, FUN = sd)
  table.fit.sd <-t(table.fit.sd)
  table.fit.sd <- as.data.frame(table.fit.sd)
  rownames(table.fit.sd) <- "SD std. dev. during training"
  
  # calculate mean standardized deviance for the testing set (i.e. quality of fit)
  table.pred.mean <- apply(df.pred.perf[list.models],2, FUN = mean)
  table.pred.mean <-t(table.pred.mean)
  table.pred.mean <- as.data.frame(table.pred.mean)
  rownames(table.pred.mean) <- "Mean std. dev. during testing"
  
  # calculate corresponding standart deviation
  table.pred.sd <- apply(df.pred.perf[list.models],2, FUN = sd)
  table.pred.sd <-t(table.pred.sd)
  table.pred.sd <- as.data.frame(table.pred.sd)
  rownames(table.pred.sd) <- "SD std. dev. during testing"
  
  # calculate mean standardized deviance for the testing set (i.e. predictive perfomance)
  names.expl.pow <- paste("expl.pow_", list.models, sep="")
  
  table.mean.exp <- apply(df.pred.perf[names.expl.pow],2, FUN = mean)
  table.mean.exp <-t(table.mean.exp)
  table.mean.exp <- as.data.frame(table.mean.exp)
  colnames(table.mean.exp) <- list.models
  rownames(table.mean.exp) <- "Mean expl. power"
  
  # calculate corresponding standart deviation
  table.sd.exp <- apply(df.pred.perf[names.expl.pow],2, FUN = sd)
  table.sd.exp <-t(table.sd.exp)
  table.sd.exp <- as.data.frame(table.sd.exp)
  colnames(table.sd.exp) <- list.models
  rownames(table.sd.exp) <- "SD expl. power"
  
  # add performance ratio
  perf.ratio <- (table.pred.mean/table.pred.sd) * table.pred.mean
  rownames(perf.ratio) <- "performance ratio"
  
  # row bind results
  table <- rbind(table.fit.mean, table.fit.sd, table.pred.mean, table.pred.sd, table.mean.exp, table.sd.exp, perf.ratio)
  
  tab1 <- table %>% gt(rownames_to_stub = T) %>% tab_header(
    title = md("**Mean predictive performance across models**") # make bold title
  ) %>%
    fmt_number(
      columns = list.models, # round numbers
      decimals = 2
    ) %>% # remove uneccessary black lines
    tab_options(
      table.border.top.color = "white",
      heading.border.bottom.color = "black",
      row_group.border.top.color = "black",
      row_group.border.bottom.color = "white",
      #stub.border.color = "transparent",
      table.border.bottom.color = "white",
      column_labels.border.top.color = "black",
      column_labels.border.bottom.color = "black",
      table_body.border.bottom.color = "black",
      table_body.hlines.color = "white")
  return(tab1)
}
    
make.table.species <- function(df.merged.perf, list.models){
  
  names(list.models) <- c()
  
  common.vect <- c("Taxa", "Prevalence", "Null_model")
  expl.pow.vect.pred <- paste("expl.pow_", list.models, ".pred", sep="")
  expl.pow.vect.fit <- paste("expl.pow_", list.models, ".fit", sep="")
  likeli.ratio.vect <- colnames(df.merged.perf)[which(grepl("likelihood.ratio", colnames(df.merged.perf)))]
  
  # Make tables
  tmp.table <- df.merged.perf
  tmp.table$Taxa <- sub("Occurrence.", "", tmp.table$Taxa)
  # colnames(tmp.table) <- c("Taxa", "Prevalence", list.models)
  #tmp.table <- tmp.table %>% mutate((across(is.numeric, round, digits=3)))
  tab2 <- tmp.table %>% gt() %>%
    tab_header(
      title = md("**Different performance accross models**") # make bold title
    ) %>%
    tab_spanner(
      label = "Training",
      columns = paste0(list.models, ".fit")
    ) # %>%
    
    fit.vect <- list.models
    names(fit.vect) <- paste0(list.models, ".fit")
      
    tab2 <- tab2 %>% cols_label(
      .list = fit.vect
    ) %>%
    tab_spanner(
      label = "Testing",
      columns = paste0(list.models, ".pred")
    ) %>%
    tab_spanner(
      label = "Explanatory power for prediction",
      columns = all_of(expl.pow.vect.pred)
    ) %>%
    tab_spanner(
      label = "Explanatory power during training",
      columns = all_of(expl.pow.vect.fit)
    ) %>%
    tab_spanner(
      label = "Likelihood ratio",
      columns = all_of(likeli.ratio.vect)
    ) %>%
    tab_spanner(
      label = "Biggest expl. pow. difference",
      columns = c(Big.model.diff, Big.pred.expl.pow.diff, chGLM.model.diff, chGLM.pred.expl.pow.diff)
    ) %>%
    fmt_number(
      columns = c("Prevalence", colnames(tmp.table)[-which(colnames(tmp.table) %in% c("Taxa", "Taxonomic level", "Big.model.diff", "chGLM.model.diff"))]), # round numbers
      decimals = 2
    ) %>% # remove uneccessary black lines
    tab_options(
      table.border.top.color = "white",
      heading.border.bottom.color = "black",
      row_group.border.top.color = "black",
      row_group.border.bottom.color = "white",
      #stub.border.color = "transparent",
      table.border.bottom.color = "white",
      column_labels.border.top.color = "black",
      column_labels.border.bottom.color = "black",
      table_body.border.bottom.color = "black",
      table_body.hlines.color = "white")
  
  return(tab2)
}

make.table.species.rearranged <- function(df.merged.perf, list.models){
  
  names(list.models) <- c()
  no.models <- length(list.models)
  
  # Make tables
  # tmp.table <- df.merged.perf[,-which(grepl("expl.pow",colnames(df.merged.perf)) & grepl(".fit",colnames(df.merged.perf)))]
  tmp.table <- df.merged.perf[,c(1, 2, which((grepl("expl.pow_",colnames(df.merged.perf)) & grepl(".pred",colnames(df.merged.perf))) | grepl("likelihood",colnames(df.merged.perf)))) ]
  
  tmp.table$Taxa <- sub("Occurrence.", "", tmp.table$Taxa)
  # colnames(tmp.table) <- c("Taxa", "Prevalence", list.models)
  #tmp.table <- tmp.table %>% mutate((across(is.numeric, round, digits=3)))
  
  tab3 <- tmp.table %>% gt() %>%
    tab_header(
      title = md("**Different performance accross models**")) # make bold title
    # ) %>%
    # cols_label(
    #   Null_model = "Null model",
    #   Big.model.diff = "Models",
    #   Big.pred.expl.pow.diff = "Value"
    # ) %>%
    # tab_spanner(
    #   label = "Biggest expl. pow. difference in pred. ",
    #   columns = c(Big.model.diff, Big.pred.expl.pow.diff, chGLM.model.diff, chGLM.pred.expl.pow.diff)
    # )
  for (l in 0:(no.models-1)) {
    col.group <- colnames(tmp.table)[which(grepl(list.models[no.models-l], colnames(tmp.table)) & !grepl("diff", colnames(tmp.table)) )]
    col.names <- c("Expl. pow.", "Likelihood ratio")
    names(col.names) <- col.group
    
    tab3 <- tab3 %>%
      cols_move(
        columns = all_of(col.group),
        after = "Prevalence"
        ) %>%
      tab_spanner(
        label = list.models[no.models-l],
        columns = all_of(col.group)
      ) %>%
      cols_label( .list = col.names
      )
  }
  
  tab3 <- tab3 %>%
    fmt_number(
      columns = c("Prevalence", colnames(tmp.table)[-which(colnames(tmp.table) %in% c("Taxa", "Taxonomic level", "Big.model.diff", "chGLM.model.diff"))]), # round numbers
      decimals = 2
    ) %>% # remove uneccessary black lines
    tab_options(
      table.border.top.color = "white",
      heading.border.bottom.color = "black",
      row_group.border.top.color = "black",
      row_group.border.bottom.color = "white",
      #stub.border.color = "transparent",
      table.border.bottom.color = "white",
      column_labels.border.top.color = "black",
      column_labels.border.bottom.color = "black",
      table_body.border.bottom.color = "black",
      table_body.hlines.color = "white")
  
  return(tab3)
}

make.table.species.rearranged.order <- function(df.merged.perf, list.models){
  
  names(list.models) <- c()
  no.models <- length(list.models)
  
  # Make tables
  tmp.table <- df.merged.perf[,-which(grepl("expl.pow",colnames(df.merged.perf)) & grepl(".fit",colnames(df.merged.perf)))]
  
  tmp.table$Taxa <- sub("Occurrence.", "", tmp.table$Taxa)
  tmp.table  <- arrange(tmp.table, desc(chGLM.pred.expl.pow.diff))
  
  tab3 <- tmp.table %>% gt() %>%
    tab_header(
      title = md("**Different performance accross models**") # make bold title
    ) %>%
    cols_label(
      Null_model = "Null model",
      Big.model.diff = "Models",
      Big.pred.expl.pow.diff = "Value"
    ) %>%
    tab_spanner(
      label = "Biggest expl. pow. difference in pred.",
      columns = c(Big.model.diff, Big.pred.expl.pow.diff, chGLM.model.diff, chGLM.pred.expl.pow.diff)
    )
  for (l in 0:(no.models-1)) {
    col.group <- colnames(tmp.table)[which(grepl(list.models[no.models-l], colnames(tmp.table)) & !grepl("diff", colnames(tmp.table)))]
    col.names <- c("Fit", "Prediction", "Expl. pow.", "Likelihood ratio")
    names(col.names) <- col.group
    
    tab3 <- tab3 %>%
      cols_move(
        columns = all_of(col.group),
        after = "Null_model"
      ) %>%
      tab_spanner(
        label = list.models[no.models-l],
        columns = all_of(col.group)
      ) %>%
      cols_label( .list = col.names
      )
  }
  
  tab3 <- tab3 %>%
    fmt_number(
      columns = c("Prevalence", colnames(tmp.table)[-which(colnames(tmp.table) %in% c("Taxa", "Taxonomic level", "Big.model.diff", "chGLM.model.diff"))]), # round numbers
      decimals = 2
    ) %>% # remove uneccessary black lines
    tab_options(
      table.border.top.color = "white",
      heading.border.bottom.color = "black",
      row_group.border.top.color = "black",
      row_group.border.bottom.color = "white",
      #stub.border.color = "transparent",
      table.border.bottom.color = "white",
      column_labels.border.top.color = "black",
      column_labels.border.bottom.color = "black",
      table_body.border.bottom.color = "black",
      table_body.hlines.color = "white")
  
  return(tab3)
}




