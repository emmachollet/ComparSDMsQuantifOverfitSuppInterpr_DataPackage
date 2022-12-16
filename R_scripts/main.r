## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                    
# PRELIMINARIES ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

getwd() # show working directory, ! set working directory to source file ('main.r') location
rm(list=ls()) # free workspace
graphics.off() # clean graphics display

## Main options ####

BDM <- F # select dataset, "All" or only "BDM" monitoring programs
file.prefix <- ifelse(BDM, "BDM_", "All_")

CV <- F              # train for cross-validation (CV)
ODG <- ifelse(CV, F, # if CV = T, no out-of-domain generalization (ODG)
                  F  # train for out-of-domain generalization (ODG)
                  )  # if CV = F and ODG = F, train on whole dataset (FIT)

ODG.info <- c(training.ratio = 0.8,     # ratio of data used for calibration in ODG
              variable = "temperature", # factor by which data is split for ODG, can be "random", "temperature", "IAR", etc
              model = "lme.area.elev")  # if variable = 'temperature', choose by which one to split, can be "lme.area.elev", "lm.area.elev", etc

train.info <- ifelse(CV, "CV_",
                     ifelse(ODG, paste0(paste(c("ODG_", ODG.info["training.ratio"], ODG.info["variable"]), collapse = ""), "_"),
                            "FIT_"))

dl <- F                  # allow data leakage in standardization between training and testing set
if(!CV){ dl <- F }       # if it's only fitting, no dataleakage

models.analysis <- c(    # comparison analysis of different models structures:
  "rf.hyperparam"   = F, # RF regularization
  "ann.hyperparam"  = F, # ANN hyperparameter tuning
  "rf.random"       = F, # sensitivity analysis RF randomness
  "ann.random"      = F  # sensitivity analysis ANN randomness
)

case.compar <- c(     # comparison analysis of models trained on different datasets:
  "cv.odg"       = ifelse(!CV & !ODG, F, T), # CV vs ODG
  "dl"           = F, # data leakage (dl) vs no data leakage (in data standardization)
  "temp.model"   = F  # different temperature models
)

server <- F # run the script on the server (changes number of cores)
plot.all.ICE <- F # produce ICE/PDP plots for paper

## Libraries ####

# Set a checkpoint to use same library versions, which makes the code repeatable over time
if ( !require("checkpoint") ) { install.packages("checkpoint"); library("checkpoint") }
checkpoint("2022-01-01", r_version = "4.1.1") # replace with desired date and R version

if ( !require("parallel") ) { install.packages("parallel"); library("parallel") }           # to run things in parallel

# Data management
if ( !require("dplyr") ) { install.packages("dplyr"); library("dplyr") }                    # to sort, join, merge data
if ( !require("tidyr") ) { install.packages("tidyr"); library("tidyr") }                    # to sort, join, merge data
if ( !require("splitTools") ) { install.packages("splitTools"); library("splitTools") }     # to split the data
if ( !require("vtable") ) { install.packages("vtable"); library("vtable") }                 # to make table with summary statistics
if ( !require("pROC") ) { install.packages("pROC"); library("pROC") }                       # to compute AUC

# Plots
if ( !require("ggplot2") ) { install.packages("ggplot2"); library("ggplot2") }              # to do nice plots
if ( !require("gridExtra") ) { install.packages("gridExtra"); library("gridExtra") }        # to arrange multiple plots on a page
if ( !require("plot.matrix") ) { install.packages("plot.matrix"); library("plot.matrix") }  # to plot nice tables
if ( !require("scales") ) { install.packages("scales"); library("scales") }                 # to look at colors
if ( !require("reshape2") ) { install.packages("reshape2"); library("reshape2") }           # to reshape dataframes
if ( !require("DiagrammeR")) { install.packages("DiagrammeR"); library("DiagrammeR") }      # to plot trees of BCT
if ( !require("skimr") ) { install.packages("skimr"); library("skimr") }                    # to show key descriptive stats
if ( !require("corrplot") ) { install.packages("corrplot"); library("corrplot") }           # to plot correlation matrix
if ( !require("gt") ) { install.packages("gt"); library("gt") }                             # to make tables
if ( !require("sf") ) { install.packages("sf"); library("sf") }                             # to read layers for plotting maps
if ( !require("ggpubr") ) { install.packages("ggpubr"); library("ggpubr") }                 # to arrange multiple plots on a page

# Statistical models
if ( !require("rstan") ) { install.packages("rstan"); library("rstan") }                    # to run hierarchical models written in Stan

# Artificial Neural Networks (ANN)
if ( !require("reticulate") ) { install.packages("reticulate"); library("reticulate") }
# install_miniconda()              # run this the very first time reticulate is installed
# install.packages("tensorflow")
library("tensorflow")
# install_tensorflow()             # run this line only when opening new R session
# install.packages("keras")
library("keras")
# install_keras()                  # run this line only when opening new R session
# use_condaenv()

# Machine Learning (ML) models
if ( !require("mgcv") ) { install.packages("mgcv"); library("mgcv") }                       # to run Generalized Additive Model (GAM) algorithm
if ( !require("gam") ) { install.packages("gam"); library("gam") }                          # to run Generalized Additive Model (GAM) algorithm
if ( !require("kernlab") ) { install.packages("kernlab"); library("kernlab") }              # to run Support Vector Machine (SVM) algorithm
if( !require("xgboost") ) { install.packages("xgboost"); library("xgboost") }               # to run Boosted Classification Trees (BCT)
if ( !require("randomForest")) {install.packages("randomForest"); library("randomForest") } # to run Random Forest (RF)

# have to be loaded at the end to not cache function 'train'
if ( !require("caret") ) { install.packages("caret"); library("caret") }                    # comprehensive framework to build machine learning models


## Functions ####

source("stat_model_functions.r")
source("ml_model_functions.r")
source("ann_model_functions.r")
source("plot_functions.r")
source("utilities.r")


## Data ####

# Define directory and files
dir.input.data        <- "../Input_data/"
dir.workspace         <- "../Output_data/Tables/"
dir.models.output     <- "../Output_data/Trained_models/"
dir.plots.output      <- "../Plots/Models_analysis_plots/"
dir.expl.plots.output <- "../Plots/Explorative_plots/"
file.input.data       <- "All_2729samples_9envfact_lme.area.elev_ModelInputs.csv"
file.prev.taxa        <- "All_2729samples_9envfact_lme.area.elev_PrevalenceTaxa.csv"

# Load input datasets
data      <- read.csv(paste0(dir.input.data, file.input.data), header = T, sep = ",", stringsAsFactors = F)
prev.inv  <- read.csv(paste0(dir.input.data, file.prev.taxa), header = T, sep = ",", stringsAsFactors = F)

# Prepare inputs for geographic plots
inputs    <- map.inputs(directory = paste0(dir.input.data,"Swiss.map.gdb"), data = data)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DATA WRANGLING ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Select env. factors ####

# Select temperature of interest
print(colnames(data)[which(grepl("temperature", colnames(data)))])
select.temp <- "lme.area.elev"
# select.temp <- "lm.area.elev"
# data$temperature <- data[,which(colnames(data) == paste0("temperature.", select.temp))]

# Add quadratic terms for temperature and flow velocity
# data$temperature2 <- data$temperature^2
# data$velocity2 <- data$velocity^2

# Select environmental factors
env.fact <- c("Temperature"                        = "temperature",       # Temp
              "Flow velocity"                      = "velocity",          # FV
              "Riparian agriculture"               = "A10m",              # A10m
              "Livestock unit density"             = "cow.density",       # LUD
              "Insecticide application rate"       = "IAR",               # IAR
              "Urban area"                         = "urban.area",        # Urban
              "Forest-river intersection"          = "FRI",               # FRI
              "Forest-river intersection buffer"   = "bFRI",              # bFRI
              "Width variability"                  = "width.variability") # WV

env.fact.full <- c(env.fact,
                   "Temperature2" = "temperature2",
                   "Velocity2"    = "velocity2")

no.env.fact <- length(env.fact)

# Select sites (categorical) information to keep for analysis
vect.info <- c("Region", "RiverBasin", "BIOGEO")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Preprocess data ####

# Select list of taxa
list.taxa <- colnames(data)[which(grepl("Occurrence.", colnames(data)))]
names.taxa <- gsub("Occurrence.", "", list.taxa)
no.taxa <- length(list.taxa)
select.taxa <- c("Occurrence.Gammaridae", "Occurrence.Psychodidae", "Occurrence.Nemouridae") # select taxa for further analysis
list.taxa.int <- prev.inv[which(prev.inv[, "Prevalence"] < 0.75 & prev.inv[,"Prevalence"] > 0.25), # list taxa with intermediate prevalence
                          "Occurrence.taxa"]
no.taxa.int <- length(list.taxa.int)

cat("\nSummary information of input dataset:\n",
    length(unique(data$SampId)), "samples,\n",
    length(unique(data$SiteId)), "sites,\n",
    length(list.taxa), "taxa.")
print(summary(as.factor(prev.inv[which(prev.inv$Occurrence.taxa %in% list.taxa), "Taxonomic.level"])))

# Write info for files name
info.file.name.data <- paste0(file.prefix,
                              dim(data)[1], "samples_",
                              no.env.fact, "envfact_",
                              select.temp)

info.file.name.models <- paste0(file.prefix,
                                no.taxa, "taxa_",
                                train.info,
                                ifelse(dl, "DL_", "no_DL_"),
                                select.temp)
cat("\nMain info: ", info.file.name.models, "\n\n")

# Split data
if(CV|ODG){
  splits <- split.data(data = data, list.taxa = list.taxa, dir.workspace = dir.workspace, 
                       info.file.name.data = info.file.name.data, select.temp = select.temp,
                       CV = CV, ODG = ODG, ODG.info = ODG.info, bottom = T)
} else {
  splits <- list("Entire dataset" = data)
}
list.splits <- names(splits)
no.splits <- length(list.splits)
n.cores.splits <- ifelse(server, no.splits, 1) # to use one core per split on the server

# Standardize data
stand.norm.data <- standardize.data(data = data, splits = splits, env.fact.full = env.fact.full, dl = dl, CV = CV, ODG = ODG)
standardized.data <- stand.norm.data$standardized.data
standardized.data.factors <- stand.norm.data$standardized.data.factors
normalization.data <- stand.norm.data$normalization.data
remove(stand.norm.data)

# Produce correlation matrix
file.name <- paste0(dir.expl.plots.output, info.file.name.data, "_CorrelationMatrix")
if(!file.exists(paste0(file.name,".pdf"))){
  pdf.corr.mat.env.fact(data = data, env.fact = env.fact, file.name = file.name)
}

# Analyze splits
file.name <- paste0("_AnalysisSplits_", train.info, ifelse(ODG, paste0("_by", ODG.info["model"],""), ""))
if(!file.exists(paste0(dir.expl.plots.output, info.file.name.data, file.name, ".pdf"))){
  list.plots <- analysis.splits(inputs = inputs, splits = splits, env.fact = env.fact, vect.info = vect.info)
  cat("Printing PDF and PNG\n")
  print.pdf.plots(list.plots = list.plots, width = 15, height = 8, dir.output = dir.expl.plots.output, info.file.name = info.file.name.data, file.name = file.name, png = T, png.ratio = 0.8)
}

# Analyze prevalence
file.name <- paste0("_AnalysisPrevalence_", train.info, ifelse(ODG, paste0("_by", ODG.info["model"],""), ""))
if(!file.exists(paste0(dir.expl.plots.output, info.file.name.data, file.name, ".pdf")) & !CV){
  list.plots <- analysis.prevalence(prev.inv = prev.inv, splits = splits, ODG = ODG)
  cat("Printing PDF\n")
  print.pdf.plots(list.plots = list.plots, width = 15, height = 15, dir.output = dir.expl.plots.output, info.file.name = info.file.name.data, file.name = file.name, 
                  png = T,png.vertical = T,  png.ratio = 1)
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Select models ####
# Already select the colors assigned to each model for the plots

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Statistical models
list.stat.mod <-  c(
  "#256801" = "hGLM", 
  "#59AB2D" = "chGLM"
)
stat.iterations <- 10000 # set iterations for stat models training (needs to be an even number)
n.chain  <- 2 # number of sample chains ran in parallel
n.cores.stat.models <- ifelse(server, length(list.stat.mod), 1) # to use one core per model on the server
comm.corr.options <- c(F,T) # flags for correlation matrix
names(comm.corr.options) <- list.stat.mod

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Machine learning algorithms (! their packages have to be installed first)
list.algo <- c( 
  "deepskyblue4"  = 'glm',       # Generalized Linear Model
  "deepskyblue"   = 'gamLoess',  # Generalized Additive Model
  "#7B1359"       = 'svmRadial', # Support Vector Machine
  "hotpink1"      = "ada" ,      # Boosted Classification Tree
  "hotpink3"      = 'rf'         # Random Forest
)
no.algo <- length(list.algo)

list.ml <- c(
  "iGLM", # Generalized Linear Model
  "GAM",  # Generalized Additive Model
  "SVM",  # Support Vector Machine
  "BCT",  # Boosted Classification Tree
  "RF"    # Random Forest
)
names(list.ml) <- names(list.algo)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Artificial Neural Network

# Hyperparameters with fixed value
learning.rate <- 0.01 # only one value
batch.size <- 64 # only one value

# Hyperparameters to tune
if(models.analysis["ann.hyperparam"]){
  grid.hyperparam <- expand.grid(
    layers  = c(3, 
                5, 10),
    units   = c(4, 8, 
                32),
    act.fct = c( 
      "tanh",
      # "swish",
      "relu",
      "leakyrelu"),
    no.epo  = c(50, 100)
  )
} else {
  grid.hyperparam <- expand.grid(
    layers  = c(5),
    units   = c(4),
    act.fct = c("leakyrelu"),
    no.epo  = c(50)
  )
}
grid.hyperparam$act.fct <- as.character(grid.hyperparam$act.fct)
no.hyperparam <- nrow(grid.hyperparam)
list.hyper.param <- vector("list", no.hyperparam)

if(no.hyperparam == 1){ # fixed hyperparameters
  if(models.analysis["ann.random"]){
    no.ann.runs = 10 # number of runs to do randomness sensitivity analysis
    for (n in 1:no.ann.runs) {
      list.hyper.param[[n]] <- grid.hyperparam[1,]
      names(list.hyper.param)[n] <- paste0("ANN_rand", n)
    }
  } else {
    list.hyper.param[[1]] <- grid.hyperparam[1,]
    names(list.hyper.param) <- c("ANN")
  }
} else { # hyperparameters tuning
  for (n in 1:no.hyperparam) {
    list.hyper.param[[n]] <- grid.hyperparam[n,]
    names(list.hyper.param)[n] <- paste(paste0(grid.hyperparam[n,], c("L", "U_", "", "epo")), collapse = "")
  }
  names(list.hyper.param) <- paste("ANN_", names(list.hyper.param), sep = "")
}
list.ann <- names(list.hyper.param)
no.ann <- length(list.ann)
if(no.ann == 1){
  names(list.ann) <- "#FFB791"
} else {
  names(list.ann) <- hcl.colors(no.ann, palette = "Oranges")
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Final model list
if(models.analysis["ann.hyperparam"]){
  list.models <- list.ann
} else {
  list.models <- c(list.ml[1], list.stat.mod, list.ml[2:no.algo], list.ann)
  list.models.orig <- list.models
}
no.models <- length(list.models)
show_col(names(list.models))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                    
# APPLY MODELS ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Define null model, outputs taxon prevalence everywhere
null.model <- apply.null.model(data = data, list.taxa = list.taxa)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if(!models.analysis["ann.hyperparam"]){

## Statistical models ####

ptm <- proc.time() # to calculate time of simulation

stat.outputs <- mclapply(comm.corr.options, mc.cores = n.cores.stat.models, function(comm.corr){
  
  # comm.corr <- comm.corr.options[[1]]
  model.name <- ifelse(comm.corr, list.stat.mod[2], list.stat.mod[1])
  info.file.stat.name <- paste0("Stat_model_",
                                model.name, "_",
                                stat.iterations,"iterations_",
                                info.file.name.models)

  file.name <- paste0(dir.models.output, info.file.stat.name, ".rds")
  # cat(file.name)

  # If the file with the output already exist, just read it
  if (file.exists(file.name)){
      if(!exists("stat.output")){
          cat("\nFile with statistical model output already exists, we read it from:\n", file.name, "\nand save it in object 'stat.output'.\n")
          stat.output <- readRDS(file = file.name)
          }
      else{
          cat("\nList with statistical model output already exists as object 'stat.output' in this environment.\n")
      }
  } else {
      cat("\nNo statistical model output exist yet, we produce it and save it in", file.name)
      if(CV|ODG){
        stat.output <- mclapply(standardized.data, mc.cores = n.cores.splits, FUN = stat_mod_cv, CV, ODG, model.name, comm.corr, stat.iterations, n.chain, list.taxa, env.fact.full)
        cat("\nSaving output of statistical models in:\n", file.name)
        saveRDS(stat.output, file = file.name, version = 2) #version two here is to ensure compatibility across R versions

      } else {
        stat.output <- stat_mod_cv(standardized.data, CV, ODG, model.name, comm.corr, stat.iterations, n.chain, list.taxa, env.fact.full)
        cat("\nSaving output of statistical models in\n", file.name)
        saveRDS(stat.output, file = file.name, version = 2)
        }

  }
  return(stat.output)
})

# Process output from stat models to fit structure of ml models (makes plotting easier)
stat.outputs.transformed <- transfrom.stat.outputs(stat.outputs = stat.outputs, list.taxa = list.taxa, CV = CV, ODG = ODG)

print(paste("\nSimulation time of statistical model "))
print(proc.time()-ptm)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ML ####

ptm <- proc.time() # to calculate time of simulation

# Split computation of ML algorithm, because it' too heavy to read in
if(CV|ODG){ temp.list.algo <- list(list.algo[1:3], list.algo[4:length(list.algo)])
} else { temp.list.algo <- list(list.algo) }
temp.outputs <- list()

for (n in 1:length(temp.list.algo)) {
  # n = 1
  info.file.ml.name <-  paste0("ML_models_",
                               paste0(paste(temp.list.algo[[n]], collapse = "_"), "_"),
                               info.file.name.models
                               )
  
  file.name <- paste0(dir.models.output, info.file.ml.name, ".rds")

  if(file.exists(file.name)){
    cat("\nThe file already exists:\n", file.name, "\nReading it and uploading it in the environment.\n")
    if(CV | ODG){ temp.outputs[[n]] <- readRDS(file = file.name)
    } else { temp.outputs[[n]] <- readRDS(file = file.name) }
  } else {
    cat("\nNo ML outputs exist yet, we run the models and save the results in:\n", file.name)
    if(CV|ODG){
        if(server){
          # Compute three splits in parallel (should be run on the server)
          temp.outputs[[n]] <- mclapply(standardized.data.factors, mc.cores = n.cores.splits,
                                   FUN = apply.ml.model, tune.grid = NULL, temp.list.algo[[n]], list.taxa, env.fact, env.fact.full, CV, ODG, prev.inv = prev.inv)
        } else {
          # Compute one split after the other
          temp.outputs[[n]] <- lapply(standardized.data.factors, FUN = apply.ml.model, tune.grid = NULL, temp.list.algo[[n]], list.taxa, env.fact, env.fact.full, CV, ODG, prev.inv = prev.inv)
        }
        cat("\nSaving outputs of algorithms in:\n", file.name)
        saveRDS(temp.outputs[[n]], file = file.name, version = 2)
      } else {
        splitted.data <- list("Training data" =  standardized.data.factors[[1]], "Testing data" = data.frame())
        temp.outputs[[n]] <- apply.ml.model(splitted.data = splitted.data, tune.grid = NULL, list.algo = temp.list.algo[[n]], list.taxa = list.taxa,
                                      env.fact = env.fact, env.fact.full = env.fact.full, CV = CV, ODG = ODG, prev.inv = prev.inv)
        cat("\nSaving outputs of algorithms in:\n", file.name)
        saveRDS(temp.outputs[[n]], file = file.name, version = 2)
      }
  }
}

print(paste("Simulation time of different models ", info.file.ml.name))
print(proc.time()-ptm)

} # only ann analysis bracket

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ANN ####

info.file.ann.name <-  paste0("ANN_models_",
                              no.ann, "ann_",
                              ifelse(models.analysis["ann.hyperparam"], "HyperparamTuning_", ""),
                              ifelse(models.analysis["ann.random"], "RandomAnalysis_", ""),
                              info.file.name.models)
file.name <- paste0(dir.models.output, info.file.ann.name, ".rds")
cat(file.name)

if(file.exists(file.name) & plot.all.ICE == F){
  cat("The file already exists. Reading it", file.name)
  if(CV | ODG){
    ann.outputs.cv <- readRDS(file = file.name)
  } else {
    ann.outputs <- readRDS(file = file.name)
  }
} else {
  cat("This ANN output doesn't exist yet, we produce it and save it in", file.name)
  if(CV|ODG){
    ann.outputs.cv <- lapply(standardized.data, function(split){
      lapply(list.hyper.param, FUN = build_and_train_model, split = split,
             env.fact = env.fact, list.taxa = list.taxa,
             learning.rate = learning.rate, batch.size = batch.size,
             CV = CV, ODG = ODG)
    })
    saveRDS(ann.outputs.cv, file = file.name)

  } else {
    ann.outputs <- lapply(list.hyper.param, FUN = build_and_train_model, split = standardized.data,
                          env.fact = env.fact, list.taxa = list.taxa,
                          learning.rate = learning.rate, batch.size = batch.size,
                          CV = CV, ODG = ODG)
    saveRDS(ann.outputs, file = file.name)
  }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Analysis RF ####

if(CV|ODG){
  if(models.analysis["rf.hyperparam"]){ # hyperparameter grid search
    
    list.tuned.grid <- list() # list of dataframes with columns for each hyperparameter
    
    # For 'rf'
    tuned.algo <- "rf"
    tuned.grid <- expand.grid(mtry = c(1,2,4,8))
    # tuned.grid <- list(1,2,4,8)
    
    # # For 'RRF'
    # tuned.algo <- "RRF"
    # tuned.grid <- expand.grid(mtry = c(1,2), coefReg = c(1), coefImp = c(0, 0.5))
    
    for (n in 1:nrow(tuned.grid)) {
      list.tuned.grid[[n]] <- tuned.grid[.(n),]
      names(list.tuned.grid)[[n]] <- paste(c(tuned.algo, "_", paste(tuned.grid[n,], colnames(tuned.grid), sep="")), collapse = "")
    }
    
    list.tuned.algo <- names(list.tuned.grid)
    no.analyzed.algo <- length(list.tuned.algo)
    names(list.tuned.algo) <- hcl.colors(no.analyzed.algo, palette = "Magenta")
    print(list.tuned.algo)
    
    info.file.ml.name <-  paste0("ML_models_",
                                 no.analyzed.algo, "tuned_",
                                 tuned.algo, "_",
                                 info.file.name.models)
    file.name <- paste0(dir.models.output, info.file.ml.name, ".rds")
    
    if(file.exists(file.name)){
      analysis.ml.outputs.cv <- readRDS(file.name)
    } else {
      analysis.ml.outputs.cv <- lapply(standardized.data.factors, function(split){
        # split <- standardized.data.factors[[1]]
        lapply(list.tuned.grid, FUN = apply.ml.model, splitted.data = split, list.algo = tuned.algo, list.taxa = list.taxa,
               env.fact = env.fact, CV = CV, ODG = ODG, prev.inv = prev.inv)
        })
      cat("\nSaving models outputs in:", file.name, "\n")
      saveRDS(analysis.ml.outputs.cv, file = file.name, version = 2)
    }
    
    list.models <- c(list.models, list.tuned.algo)
    
  } else if(models.analysis["rf.random"]){ # randomness analysis
    
    temp.list.algo <- list.algo[4:5] # select BCT and RF
    print(temp.list.algo)
    seeds <- c(
      2020,
      2021,
      2022)
    
    grid.algo.seeds <- expand.grid("algorithm" = temp.list.algo, "seed" = seeds)
    grid.algo.seeds$algorithm <- as.character(grid.algo.seeds$algorithm)
    grid.algo.seeds$seed <- as.integer(grid.algo.seeds$seed)
    no.algo.seeds <- nrow(grid.algo.seeds)
    list.algo.seeds <- vector("list", no.algo.seeds)
    for (n in 1:no.algo.seeds) {
      list.algo.seeds[[n]] <- grid.algo.seeds[n,]
      names(list.algo.seeds)[n] <- paste(grid.algo.seeds[n,], collapse = "")
    }
    list.seed.algo <- names(list.algo.seeds) # make list of new algo names
    for (l in temp.list.algo) { # assign colors
      names(list.seed.algo)[which(grepl(l, list.seed.algo))] <- colorRampPalette(c("white", names(temp.list.algo[which(temp.list.algo == l)])))(length(seeds))
    }
    list.models <- c(list.models, list.seed.algo) # append to list.models
    
    info.file.ml.name <-  paste0("ML_models_",
                                 no.algo.seeds, "RandomAnalysis_",
                                 info.file.name.models)
    file.name <- paste0(dir.models.output, info.file.ml.name, ".rds")
    
    if(file.exists(file.name)){
      analysis.ml.outputs.cv <- readRDS(file.name)
    } else {
      analysis.ml.outputs.cv <- lapply(standardized.data.factors, function(split){
        # split <- standardized.data.factors[[1]]
        lapply(list.algo.seeds, function(algo.seed){apply.ml.model(splitted.data = split, list.algo = algo.seed[,"algorithm"], list.taxa = list.taxa,
                                                                   env.fact = env.fact, CV = CV, ODG = ODG, prev.inv = prev.inv, seed = algo.seed[,"seed"])
        })
      })
      cat("\nSaving models outputs in:", file.name, "\n")
      saveRDS(analysis.ml.outputs.cv, file = file.name, version = 2)
    }
  }
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                    
# PROCESS OUTPUTS ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## Merge outputs ####

if(CV | ODG){
  # Merge all CV outputs in one
  outputs.cv <- vector(mode = "list", length = no.splits)
  names(outputs.cv) <- list.splits

  for (s in list.splits) {
    #s = "Split3"
    if(exists("temp.outputs")){
      outputs.cv[[s]][[list.ml[1]]] <- temp.outputs[[1]][[s]][[list.algo[1]]] # iGLM
    }
    if(exists("stat.outputs.transformed")){
        outputs.cv[[s]][[list.stat.mod[1]]] <- stat.outputs.transformed[[1]][[s]] # hGLM
        outputs.cv[[s]][[list.stat.mod[2]]] <- stat.outputs.transformed[[2]][[s]] # chGLM
    }
    if(exists("temp.outputs")){ # rest of ml algo
      for (n in 1:length(temp.outputs)) {
        if(n == 1){
          for (l in 2:length(temp.outputs[[n]][[s]]) ) {
            outputs.cv[[s]][[list.ml[l]]] <- temp.outputs[[n]][[s]][[list.algo[l]]]
          }
        } else {
          for (l in (length(temp.outputs[[n-1]][[s]])+1):(length(temp.outputs[[n-1]][[s]])+length(temp.outputs[[n]][[s]])) ) {
            outputs.cv[[s]][[list.ml[l]]] <- temp.outputs[[n]][[s]][[list.algo[l]]]
          }
        }
      }
    }
    if(exists("ann.outputs.cv")){
      outputs.cv[[s]] <- append(outputs.cv[[s]], ann.outputs.cv[[s]])
    }
    if(exists("analysis.ml.outputs.cv")){
      if(models.analysis["rf.random"]){ # correct the weird listing for this output
        for (l in list.seed.algo) {
          analysis.ml.outputs.cv[[s]][[l]] <- analysis.ml.outputs.cv[[s]][[l]][[1]]
        }
      }
      if(models.analysis["rf.hyperparam"]){ # correct the weird listing for this output
        for (l in list.tuned.algo) {
          analysis.ml.outputs.cv[[s]][[l]] <- analysis.ml.outputs.cv[[s]][[l]][[1]]
        }
      }
      outputs.cv[[s]] <- append(outputs.cv[[s]], analysis.ml.outputs.cv[[s]])
    }
  }
  
} else {
  # Make final outputs as list
  outputs <- list()
  
  if(exists("temp.outputs")){
    outputs[[list.ml[1]]] <- temp.outputs[[1]][[1]]
  }
  if(exists("stat.outputs.transformed")){
    outputs[[list.stat.mod[1]]] <- stat.outputs.transformed[[1]]
    outputs[[list.stat.mod[2]]] <- stat.outputs.transformed[[2]]
  } 
  if(exists("temp.outputs")){
    for (l in 2:no.algo) {
      outputs[[list.ml[l]]] <- temp.outputs[[1]][[list.algo[l]]]
    } 
  }
  if (exists("ann.outputs")){
    outputs <- append(outputs, ann.outputs)
  }
}

# remove(temp.outputs)
if(exists("outputs")){  
  # Make sure list.models match with models outpur
  if(CV|ODG){
    vec1 <- names(outputs.cv$Split1)
    vec2 <- list.models
    list.models <- vec2[match(vec1, vec2)]
  } else {
    vec1 <- names(outputs)
    vec2 <- list.models
    list.models <- vec2[match(vec1, vec2)]
  }
} 
print(list.models)
no.models <- length(list.models)
show_col(names(list.models))

info.file.name <- paste0(no.models, "models_",
                         ifelse(models.analysis["rf.hyperparam"], paste0(tuned.algo, "_HyperparamTuning_"), ""),
                         ifelse(models.analysis["rf.random"], "RF_RandomAnalysis_", ""),
                         ifelse(models.analysis["ann.hyperparam"], "ANN_HyperparamTuning_", ""),
                         ifelse(models.analysis["ann.random"], "ANN_RandomAnalysis_", ""),
                         info.file.name.models, "_")
cat(info.file.name)

## Results tables ####

# Produce final outputs with mean performance across splits
if(CV | ODG){
  # Make final outputs as tables
  df.cv <- make.df.outputs(outputs = outputs.cv, 
                           list.models = list.models,
                           list.taxa = list.taxa, list.splits = list.splits,
                           null.model = null.model, prev.inv = prev.inv, CV = CV, ODG = ODG)
  df.pred.perf.cv <- df.cv$`Table predictive performance CV`
  df.pred.perf <- df.cv$`Table predictive performance`
  df.fit.perf.cv <- df.cv$`Table fit performance CV`
  df.fit.perf <- df.cv$`Table fit performance`
  df.perf <- df.cv$`Table merged`
  remove(df.cv)
    
} else {
  # Make final outputs as tables
  df.perf <- make.df.outputs(outputs = outputs, list.models = list.models,
                           list.taxa = list.taxa, list.splits = list.splits,
                           null.model = null.model, prev.inv = prev.inv, CV = CV, ODG = ODG)
}

# Make csv file with performance tables
df.results <- df.perf
df.results <- df.results %>% 
  mutate_if(is.numeric, round, digits=2)
df.results$Taxa <- gsub("Occurrence.", "", df.results$Taxa)

# All results
file.name <- paste(dir.workspace, info.file.name, "TableAllResults", ".csv", sep="")
cat(file.name)
write.table(df.results, file.name, sep=",", row.names=F, col.names=TRUE)

# Summary statistics
file.name <- paste(dir.workspace, info.file.name, "TableSummaryStat", ".csv", sep="")
cat(file.name)
df.summary <- sumtable(df.results, add.median = TRUE, out='return')
write.table(df.summary, file.name, sep=",", row.names=F, col.names=TRUE)

# Standardized deviance
df.stand.dev <-  df.results[,c(which(colnames(df.results) %in% c("Taxa", "Prevalence")), which(grepl("dev.pred", colnames(df.results)) & !grepl("expl.pow", colnames(df.results))))]
colnames(df.stand.dev)[which(grepl("dev.pred", colnames(df.stand.dev)))] <- list.models
file.name <- paste(dir.workspace, info.file.name, "TableStandDev", ".csv", sep="")
cat(file.name)
write.table(df.stand.dev, file.name, sep=",", row.names=F, col.names=TRUE)

# Likelihood ratio
df.likelihood.ratio <- df.results[,c(which(colnames(df.results) %in% c("Taxa", "Prevalence")), which(grepl("likeli", colnames(df.perf))))]
colnames(df.likelihood.ratio)[which(grepl("likeli", colnames(df.likelihood.ratio)))] <- list.models
file.name <- paste(dir.workspace, info.file.name, "TableLikelihoodRatio", ".csv", sep="")
cat(file.name)
write.table(df.likelihood.ratio, file.name, sep=",", row.names=F, col.names=TRUE)

# AUC
df.auc.pred <- df.results[,c(which(colnames(df.results) %in% c("Taxa", "Prevalence")), which(grepl("auc.pred", colnames(df.perf))))]
colnames(df.auc.pred)[which(grepl("auc", colnames(df.auc.pred)))] <- list.models
file.name <- paste(dir.workspace, info.file.name, "TableAUC", ".csv", sep="")
cat(file.name)
write.table(df.auc.pred, file.name, sep=",", row.names=F, col.names=TRUE)


# Explore likelihood ratio
if((CV|ODG) & "iGLM" %in% list.models){
  temp.df <- arrange(df.perf, desc(iGLM.likelihood.ratio))
  temp.df <- filter(temp.df, iGLM.likelihood.ratio > 1)
  cat("Following", length(temp.df$Taxa), "taxa have a likelihood ratio above 1 for iGLM:\n",
      gsub("Occurrence.", "", temp.df$Taxa), "summarized like:\n"
      )
  print(summary(temp.df$iGLM.likelihood.ratio))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PLOTS ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## CV vs ODG ####

if(case.compar["cv.odg"]){
  names.appcase <- c("CV", "ODG")
  if(CV){ df.perf1 <- df.perf } else if(ODG){ df.perf2 <- df.perf }
  CV <- !CV
  ODG <- !ODG
} else if(case.compar["temp.model"]){
  names.appcase <- c("lm", "lme")
  df.perf1 <- df.perf
  lme.temp <- !lme.temp
}

temp.info.file.name <- paste0(no.models, "models_",
                              ifelse(models.analysis["rf.hyperparam"], "RF_HyperparamTuning_", ""),
                              ifelse(models.analysis["rf.random"], "RF_RandomAnalysis_", ""),
                              ifelse(models.analysis["ann.hyperparam"], "ANN_HyperparamTuning_", ""),
                              ifelse(models.analysis["ann.random"], "ANN_RandomAnalysis_", ""),
                              file.prefix,
                              no.taxa, "taxa_",
                              ifelse(CV, "CV_",
                                     ifelse(ODG, paste0(paste(c("ODG_", ODG.info["training.ratio"], ODG.info["variable"]), collapse = ""), "_"),
                                            "FIT_")),
                              ifelse(dl, "DL_", "no_DL_"),
                              select.temp,"_")
cat(temp.info.file.name)

if(case.compar["cv.odg"]){
  if(CV){
    df.perf1 <- read.csv(paste(dir.workspace, temp.info.file.name, "TableAllResults", ".csv", sep=""), header=TRUE, sep=",", stringsAsFactors=FALSE)
  } else if(ODG){
    df.perf2 <- read.csv(paste(dir.workspace, temp.info.file.name, "TableAllResults", ".csv", sep=""), header=TRUE, sep=",", stringsAsFactors=FALSE)
  }
  CV <- !CV
  ODG <- !ODG
} else if(case.compar["temp.model"]){
  df.perf2 <- read.csv(paste(dir.workspace, temp.info.file.name, "TableAllResults", ".csv", sep=""), header=TRUE, sep=",", stringsAsFactors=FALSE)
  lme.temp <- !lme.temp
}

if(exists("df.perf1")){
  list.df.perf <- list(df.perf1, df.perf2)
  names(list.df.perf) <- names.appcase
} else {
  list.df.perf <- list(df.perf)
  names(list.df.perf) <- ifelse(CV, "CV",
                                ifelse(ODG, paste(c("ODG_", ODG.info["training.ratio"], ODG.info["variable"]), collapse = ""),
                                       "FIT"))
}

cat("Plots will be produce for the following application cases:\n", names(list.df.perf))
cat("List of selected taxa:", sub("Occurrence.", "", select.taxa))

# Prepare data for plots 
plot.data <- perf.plot.data(list.df.perf = list.df.perf)

# Boxplots standardized deviance
list.plots <- plot.boxplots.compar.appcase(plot.data = plot.data, list.models = list.models, models.analysis = models.analysis)
file.name <- "ModelCompar_Boxplots"
if(length(list.df.perf) != 1){
  temp.info.file.name <- gsub(
    ifelse(CV, "CV_",
           ifelse(ODG, paste0(paste(c("ODG_", ODG.info["training.ratio"], ODG.info["variable"]), collapse = ""), "_"),
                  "FIT_")), "", info.file.name)
} else { 
  temp.info.file.name <- info.file.name 
}
print.pdf.plots(list.plots = list.plots, width = 15, height = ifelse(any(models.analysis == TRUE), 21, 8), 
                dir.output = dir.plots.output, info.file.name = temp.info.file.name, file.name = file.name, 
                png = TRUE, png.vertical = ifelse(any(models.analysis == TRUE), T, F), png.ratio = 1.2)

if(!any(models.analysis == T)){
  # Standardized deviance vs prevalence
  list.plots <- plot.perfvsprev.compar.appcase(plot.data = plot.data, list.models = list.models,
                                               list.taxa = list.taxa)
  file.name <- "ModelCompar_PerfVSPrev"
  print.pdf.plots(list.plots = list.plots, width = 15, height = ifelse(length(list.df.perf) == 1, 10, 15), dir.output = dir.plots.output, info.file.name = temp.info.file.name, file.name = file.name, 
                  png = TRUE, png.square = TRUE, png.ratio = 0.8)
}

if(length(list.df.perf) != 1){
  # Table likelihood ratio
  table.likeli.ratio <- table.likelihood.ratio(list.df.perf = list.df.perf, list.models = list.models, models.analysis = models.analysis)
  file.name <- paste(dir.workspace, temp.info.file.name, "TableStatLikeliRatio", ".csv", sep="")
  cat(file.name)
  write.table(table.likeli.ratio, file.name, sep=",", row.names=F, col.names=TRUE)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Models comparison

## GLM parameters comparison ####

list.glm <- list.models[1:3]

if(CV){
  list.plots <- plot.glm.param(outputs.cv, df.perf, list.glm, list.taxa, env.fact.full, CV)
} else if(!ODG) {
  list.plots <- plot.glm.param(outputs, df.perf, list.glm, list.taxa, env.fact.full, CV)
}

if(CV | !ODG){
  file.name <- "GLMParametersComparison"
  cat(file.name)
  print.pdf.plots(list.plots = list.plots, width = 8.4, height = 11.88, # A4 proportions  
                  dir.output = dir.plots.output, info.file.name = info.file.name, file.name = file.name, png = TRUE, png.vertical = TRUE)
}

# Plots specifically related to trained models (and not to CV)

# If CV or ODG, take just the first split for trained models analysis
# Doesn't work for CV for now
if(!CV | ODG){
  outputs <- outputs.cv[[1]]
  normalization.data.cv <- normalization.data
}

if(!CV){
  if(plot.all.ICE == T){ # to produce all ICE/PDP provided in manuscript and Supp. Inf.
    temp.select.taxa <- list(select.taxa, list.taxa[which(list.taxa != "Occurrence.Perlodidae")])
    temp.select.env.fact <- list(env.fact, env.fact[1])
  } else { # 
    temp.select.taxa <- list(select.taxa[1], select.taxa[2])
    temp.select.env.fact <- list(env.fact[1], env.fact[2])
  }
  
  for (n in 1:length(temp.select.taxa)) {
    
    subselect.taxa <- temp.select.taxa[[n]]
    select.env.fact <- temp.select.env.fact[[n]]
    
    ## ICE/PDP ####
    
    no.samples <- 100
    no.steps <- 200
    subselect <- 1
    
    list.list.plots <- lapply(subselect.taxa, FUN= plot.ice.per.taxa, outputs = outputs, ref.data = standardized.data[[1]], 
                              list.models = list.models, list.taxa = list.taxa,
                              env.fact = env.fact, select.env.fact = select.env.fact,
                              normalization.data = normalization.data, ODG = ODG, 
                              no.samples = no.samples, no.steps = no.steps, subselect = subselect)
    
    for (j in 1:length(subselect.taxa)) {
      taxon <- sub("Occurrence.", "", subselect.taxa[j])
      file.name <- paste0("ICE_", no.samples, "samp_", taxon, "_", length(select.env.fact), "envfact")
      print.pdf.plots(list.plots = list.list.plots[[j]], width = 8.3, height = 11.7, # A4 format in inches 
                      dir.output = paste0(dir.plots.output, "ICE/"), 
                      info.file.name = info.file.name, file.name = file.name) #,
                      # png = T, png.vertical = T, png.ratio = 0.7)
    }
  
    ## Response shape ####
    
    list.list.plots <- lapply(subselect.taxa, FUN = plot.rs.taxa, outputs = outputs, list.models = list.models, 
                              normalization.data = normalization.data, env.fact = env.fact, CV = CV, ODG = ODG)
    
    
    for (j in 1:length(subselect.taxa)) {
      taxon <- sub("Occurrence.", "", subselect.taxa[j])
      file.name <- paste0("RS_", taxon)
      print.pdf.plots(list.plots = list.list.plots[[j]], width = 8.3, height = 11.7, # A4 format in inches
                      dir.output = paste0(dir.plots.output, "ICE/"),
                      info.file.name = info.file.name, file.name = file.name)
    }
    
  }
}
