## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 
##                      ----  Functions ML models ----
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


# Define metric function for trainControl within caret package
stand.dev <- function(data, lev = c("present","absent"), model = NULL){
    
    no.obs <- dim(data)[1]
    
    likeli <- 1:no.obs
    
    likeli[which(data$obs == lev[1])] <- data[which(data$obs == lev[1]), lev[1]]
    likeli[which(data$obs == lev[2])] <- data[which(data$obs == lev[2]), lev[2]]
    
    likeli[which(likeli < 0.0001)] <- 0.0001
    
    st.dev <- -2 * sum(log(likeli)) / no.obs
    names(st.dev) <- "StandardizedDeviance"
    return(st.dev)
}


# Function to select the best performance of stand dev (i.e the lowest one)
lowest <- function (x, metric, maximize = F){
    
    best <- which.min(x[, metric])
    return(best)
}

# Function to apply ML algorithms
apply.ml.model <- function(splitted.data, tune.grid = NULL, list.algo, list.taxa, env.fact, env.fact.full, selec.metric = "StandardizedDeviance", CV = T, ODG, prev.inv, seed = 2021,...){
    
    data.train <- splitted.data[["Training data"]]
    data.test <- splitted.data[["Testing data"]]
    
    tune.grid <- tune.grid
    no.algo <- length(list.algo)
    
    # Adapt list taxa to taxa actually present in data.train
    list.taxa <- list.taxa[list.taxa %in% colnames(data.train)]
    
    # Make a list to store the outputs of each model
    outputs <- vector(mode = 'list', length = no.algo)
    names(outputs) <- list.algo
    
    for(k in 1:no.algo){
        # k = 3
        algorithm = list.algo[k]
        
        if(algorithm == "glm"){
          expl.var <- env.fact.full # add squared terms for glms
        }else{
          expl.var <- env.fact
        }
        # 
        # Make a list with the outputs of the algorithm for each taxon in list.taxa
        list.outputs <- vector(mode = 'list', length = length(list.taxa))
        names(list.outputs) <- list.taxa
        
        # if testing data set exists, create outputs for results on training and testing data sets
        # else only list outputs for training data set
        if(CV|ODG){which.set <- c("training set", "testing set")
        } else {which.set <- c("training set")}
        out <- c("Observation", #1
                 "Prediction factors", #2 
                 "Prediction probabilities", #3 
                 "Likelihood", #4
                 "Performance")#, #5
                 # "Confusion matrix") #6
        output.names <- c("Trained model", c(outer(out, which.set, FUN = paste)))
        
        for (j in 1:length(list.taxa)){
            # j = 1
            temp.list <- vector(mode = 'list', length = length(output.names))
            names(temp.list) <- output.names
            
            temp.train <- na.omit(data.train[, c("SiteId", "SampId",
                                                 "X", "Y", 
                                                 list.taxa[j], expl.var)]) # create a temporary training dataset with the taxon and env fact, to 
            if(CV|ODG){temp.test <- na.omit(data.test[, c("SiteId", "SampId",
                                                   "X", "Y", 
                                                   list.taxa[j], expl.var)])
                        temp.sets <- list(temp.train, temp.test)
            } else {    temp.sets <- list(temp.train)     } 
            
            f <- reformulate(expl.var, list.taxa[j]) # write formula (target variable ~ explanatory variables) to apply the model
            
            # Why is train() better than using the algorithm's fct directly ?
            # Because train() does other things like:
            # 1. Cross validating the model
            # 2. Tune the hyper parameters for optimal model performance
            # 3. Choose the optimal model based on a given evaluation metric
            # 4. Preprocess the predictors (what we did so far using preProcess())
            
            set.seed(seed) # for reproducibility of the folds
            folds <- groupKFold(temp.train$SiteId, 3) # create 3 folds, grouped by SiteId
                
            cat(paste("\nApplying", algorithm, "to", j, list.taxa[j], "\n"))

            # Define the training control
            train.control <- trainControl(
                method = 'cv',                   # k-fold cross validation
                number = 3,                      # number of folds
                index = folds,                   # provide indices computed with groupKFold for the k-fold CV
                classProbs = T,                  # should class probabilities be returned
                summaryFunction = stand.dev ,
                selectionFunction = lowest       # we want to minimize the metric
            )
            
            model <- train(f, data = temp.train, metric = selec.metric, method = algorithm, trControl = train.control, tuneGrid = tune.grid)
            # model <- train(x = temp.train[,expl.var], y = temp.train[,list.taxa[j]], metric = selec.metric, method = algorithm, trControl = train.control)
            
            
            temp.list[["Trained model"]] <- model
            
            for(n in 1:length(which.set)){
              # n <- 1
                # Observation
                temp.list[[paste(out[1],which.set[n])]] <- temp.sets[[n]]
                n.obs <- dim(temp.sets[[n]])[1]
                
                # Prediction factors
                temp.list[[paste(out[2],which.set[n])]] <- predict(model, temp.sets[[n]])

                # Prediction probabilities
                temp.list[[paste(out[3],which.set[n])]] <- predict(model, temp.sets[[n]], type = 'prob')
 
                # Likelihood
                likeli <- 1:nrow(temp.sets[[n]])
                for(i in 1:nrow(temp.sets[[n]])){
                    if(temp.sets[[n]][i,list.taxa[j]] == "present"){
                        likeli[i] <- temp.list[[paste(out[3],which.set[n])]][i, "present"]
                    } else if (temp.sets[[n]][i,list.taxa[j]] == "absent" ){
                        likeli[i] <- temp.list[[paste(out[3],which.set[n])]][i, "absent"]
                    }
                }
                likeli[which(likeli < 0.0001)] <- 0.0001 # avoid problems when likelihood too small
                temp.list[[paste(out[4],which.set[n])]] <- likeli
                # Performance
                temp.list[[paste(out[5],which.set[n])]] <- -2 * sum(log(likeli)) / nrow(temp.sets[[n]])
                
            }
            
            list.outputs[[j]] <- temp.list
        }
    
        outputs[[k]] <- list.outputs
    }
    
    return(outputs)
}