## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 
##                     ----  Functions ANN models ----
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


# Function to build and train ANN Multilayer Perceptron

build_and_train_model <- function (hyper.param = hyper.param,
                                   split = split,
                                   env.fact = env.fact,
                                   list.taxa = list.taxa,
                                   learning.rate = learning.rate,
                                   batch.size = batch.size,
                                   CV = CV,
                                   ODG = ODG){
  # Set hyperparameters
  num.layers <- hyper.param[,1]
  num.units.per.layer <- hyper.param[,2]
  act.fct <- hyper.param[,3]
  num.epochs <- hyper.param[,4]
  
  # Adapt list taxa to taxa actually present in data.train
  list.taxa <- list.taxa[list.taxa %in% colnames(split[[1]])]
  
  # Set training and testing sets
  Xtrain <- as.matrix(split[[1]][ ,env.fact])
  Ytrain <- as.matrix(split[[1]][ ,list.taxa])
  
  # Build ANN
  model <- keras_model_sequential()
  
  if(act.fct == "leakyrelu"){
    model <- model %>% layer_dense(units = num.units.per.layer,
                                   input_shape = ncol(Xtrain)) %>% 
      layer_activation_leaky_relu()
    
    if (num.layers>1){
      for (i in 1:(num.layers-1)){
        model <- model %>% layer_dense(units = num.units.per.layer) %>% 
          layer_activation_leaky_relu()
      }
    }
  } else {
    model <- model %>% layer_dense(units = num.units.per.layer,
                                   input_shape = ncol(Xtrain),
                                   activation = act.fct)
    
    if (num.layers>1){
      for (i in 1:(num.layers-1)){
        model <- model %>% layer_dense(units = num.units.per.layer,
                                       activation = act.fct)
      }
    }
  }
  
  # Add the output layer, need a sigmoid activation function because outputs probability
  model <- model %>% layer_dense(units = ncol(Ytrain), activation = "sigmoid")
  
  summary(model)
  
  # Specify the learning rate for stochastic gradient descent
  opt <- optimizer_adam(learning_rate = learning.rate)
  
  model %>% compile(optimizer = opt,
                    loss = "binary_crossentropy", # loss_binary_crossentropy(),
                    metrics = list('accuracy')) # categorical_accuracy ?
  
  # Fit the model
  history <- model %>% fit(x = Xtrain,
                           y = Ytrain,
                           epochs = num.epochs,
                           batch_size = batch.size)
  
  # Make a list with the outputs of the algorithm for each taxon in list.taxa
  list.outputs <- vector(mode = 'list', length = length(list.taxa))
  names(list.outputs) <- list.taxa
  
  if(CV | ODG){
    which.set <- c("training set", "testing set")
    Xtest <- as.matrix(split[[2]][ ,env.fact])
    Ytest <- as.matrix(split[[2]][ ,list.taxa])
  } else { 
    which.set <- c("training set") 
    Xtest <- as.matrix(split[[1]][ ,env.fact])
    Ytest <- as.matrix(split[[1]][ ,list.taxa])
  }
  out <- c("Observation", #1
           "Prediction factors", #2 
           "Prediction probabilities", #3 
           "Likelihood", #4
           "Performance") #5
  output.names <- c("Trained model", "Training history", c(outer(out, which.set, FUN = paste)))
  
  pred <- list()
  likeli <- list()
  
  for (n in 1:length(which.set)) {
    # Predict 
    pred[[n]] <- model %>% predict(if(n == 1){ Xtrain } else { Xtest })
    
    # Compute standardized deviance
    likeli[[n]] <- if(n == 1){ Ytrain } else { Ytest }
    for (j in 1:ncol(likeli[[n]])) {
      for (i in 1:nrow(likeli[[n]])) {
        obs <- likeli[[n]][i,j]
        likeli[[n]][i,j] <- ifelse(obs == 1, pred[[n]][i,j], 1 - pred[[n]][i,j])
      }
    }
  }
  
  for (j in 1:length(list.taxa)) {
    temp.list <- vector(mode = 'list', length = length(output.names))
    names(temp.list) <- output.names
    
    temp.sets <- if(CV | ODG){ list(cbind(Xtrain, Ytrain[,j]), cbind(Xtest, Ytest[,j])) } else { list(cbind(Xtrain, Ytrain[,j])) }
    
    temp.list[[1]] <- model
    temp.list[[2]] <- history
    
    for(n in 1:length(which.set)){
      pred.prob <- pred[[n]][,j]
      pred.fact <- ifelse(pred.prob > 0.5, "present", "absent")
      temp.list[[paste(out[1],which.set[n])]] <- split[[n]][,c("SiteId", "SampId", "X", "Y", env.fact)]
      temp.list[[paste(out[1],which.set[n])]][,list.taxa[j]] <- ifelse(split[[n]][,list.taxa[j]] == 1, "present", "absent")
      temp.list[[paste(out[1],which.set[n])]][,list.taxa[j]] <- as.factor(temp.list[[paste(out[1],which.set[n])]][,list.taxa[j]])
      temp.list[[paste(out[2],which.set[n])]] <-  pred.fact
      temp.list[[paste(out[3],which.set[n])]] <-  data.frame("absent" = 1 - pred.prob, "present" = pred.prob)
      temp.list[[paste(out[4],which.set[n])]] <-  likeli[[n]][,j]
      perf <- -2 * sum(log(likeli[[n]][,j])) / nrow(temp.sets[[n]])
      temp.list[[paste(out[5],which.set[n])]] <- perf
    }
    list.outputs[[j]] <- temp.list
  }
  
  # Return the model and the training history
  return(list.outputs)
}