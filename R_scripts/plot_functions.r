## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 
##                     ----  Functions plots ----
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


## ---- Print PDF and PNG ----

print.pdf.plots <- function(list.plots, width = 12, height = width*3/4, dir.output, info.file.name = "", file.name, png = FALSE, png.square = F, png.vertical = F, png.ratio = 1){
  
  pdf(paste0(dir.output, info.file.name, file.name, ".pdf"), paper = 'special', width = width, height = height, onefile = TRUE)
  
  for (n in 1:length(list.plots)) {
    plot(list.plots[[n]])
  }
  dev.off()
  
  if(png){
  # Print a png file
    if(png.square){
      w = 1280
      h = 960
    } else if(png.vertical){
      w = 960
      h = 1280
      } else {
      w = 1280
      h = 720
    }
    
    for (n in 1:length(list.plots)) {
      # n = 1
      png(paste0(dir.output, info.file.name, file.name, n, ".png"), width = png.ratio*w, height = png.ratio*h) # size in pixels (common 16:9 --> 1920�1080 or 1280�720)
      plot(list.plots[[n]])
      dev.off()
      
    }    
  }
}

## ---- Inputs for Swiss map plots ----

map.inputs <- function(directory, data){
    
    inputs = list()
    
    # Obtain simple feature for borders of Switzerland, and major lakes and rivers
    inputs$ch <- st_read(directory, layer="switzerland", stringsAsFactors = F)
    inputs$ch <- filter(inputs$ch, NAME=="Schweiz")
    
    inputs$rivers.major <- st_read(directory, layer = "major_rivers", stringsAsFactors = F)
    inputs$lakes.major <- st_read(directory, layer = "major_lakes", stringsAsFactors = F)
    
    # Get coordinates of all sites 
    inputs$xy <- select(data, SiteId, X, Y)
    inputs$xy <- distinct(inputs$xy)
    
    return(inputs)
}

## ---- Explorative plots -----

# Figure SI A 2: analysis splits ####
analysis.splits <- function(inputs, splits, env.fact, vect.info){
  
  no.splits <- length(splits)
  plot.data <- data.frame()
  list.plots <- list()
  
  if(no.splits == 1){
    if(grepl("Split", names(splits))){
      plot.data <- bind_rows(splits[[1]][["Training data"]], splits[[1]][["Testing data"]], .id = "Split")
      plot.data$Split <- ifelse(plot.data$Split == 1, "Calibration", "Prediction")
    } else {
      plot.data <- splits[[1]]
      plot.data$Split <- names(splits)
    }
  } else {
    for (n in 1:no.splits) {
      # n=1
      temp.plot.data <- splits[[n]][["Testing data"]]
      temp.plot.data$Split <- paste0("Split",n)
      plot.data <- bind_rows(plot.data, temp.plot.data)
      
    }
  }
  plot.data$Split <- as.factor(plot.data$Split)
  plot.data <- na.omit(plot.data)
  for (v in vect.info ) {
    plot.data[which(plot.data[,v] == "RHEIN"),v] <- "RHINE"
  }  
  
  # map distribution splits
  cat("Plotting distribution splits on Swiss map\n")
  temp.list.plots <- list()
  for(v in vect.info[which(vect.info != "Region")]){
    # v <- vect.info[3]
    g <- ggplot()
    g <- g + geom_sf(data = inputs$ch, fill="#E8E8E8", color="black")
    g <- g + geom_sf(data = inputs$rivers.major, fill=NA, color="lightblue", show.legend = FALSE)
    g <- g + geom_sf(data = inputs$lakes.major, fill="lightblue", color="lightblue", show.legend = FALSE)
    g <- g + geom_point(data = plot.data, aes_string(x="X", y="Y", color= "Split", shape = v), size= 3, alpha= 0.6)
    g <- g + theme_void(base_size = 18)
    if( v == "RiverBasin"){
      g <- g + labs(shape = "Main catchment")
    }
    g <- g + theme(# plot.title = element_text(hjust = 0.5),
      panel.grid.major = element_line(colour="transparent"),
      plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines"))
    # g
    temp.list.plots <- append(temp.list.plots, list(g))
  }
  list.plots <- append(list.plots, temp.list.plots)
  
  # barplots distribution splits
  cat("Plotting distribution splits in barplots\n")
  temp.list.plots <- list()
  for(v in vect.info){
    # v <- vect.info[1]
    p <- ggplot(data = plot.data, aes_string(x=v)) 
    p <- p + geom_histogram(aes(fill = Split), stat="count", position = "dodge")
    p <- p + theme_bw(base_size = 18)
    p <- p + labs(y = "Number of samples"#,
                  # title = "Distribution of the samples in the splits across regions"
    )
    # p
    temp.list.plots <- append(temp.list.plots, list(p))
  }
  list.plots <- append(list.plots, temp.list.plots)
  
  plot.data1 <- gather(plot.data[,c("Split", vect.info, env.fact)], key = factor, value = value, -c("Split", all_of(vect.info)))
  names.env.fact <- names(env.fact)
  names(names.env.fact) <- env.fact
  env.fact_labeller <- function(variable,value){
    return(names.env.fact[value])
  }
  cat("Plotting distribution env. fact per splits/region\n")
  temp.list.plots <- list()
  temp.vect.info <- c("Split", vect.info)
  suppressWarnings(
    for(v in temp.vect.info){
      # v <- vect.info[1]
      q <- ggplot(data = plot.data1, aes(x=value)) 
      q <- q + geom_density(aes_string(colour = v, fill = v), alpha = 0.3)
      q <- q + theme_bw(base_size = 18)
      q <- q + facet_wrap(~factor, scales = "free",
                          labeller = env.fact_labeller, 
                          strip.position="top")
      q <- q + theme(legend.text = element_text(size=22),
                     strip.background = element_rect(fill = "white"))
      q <- q + labs(x = "Values of environmental factor",
                    y = "Density"#,
                    # title = "Distribution of the environmental factors in splits"
      )
      # q
      temp.list.plots <- append(temp.list.plots, list(q))
    }
  )
  list.plots <- append(list.plots, temp.list.plots)
  
  return(list.plots)
}

# Figure SI X X: spatial distribution env. fact. ####
maps.env.fact <- function(inputs, list.env.fact, data){
  
  plot_maps_variable = function (pair.fact.name, inputs, data) {
    variable <- pair.fact.name[1]
    print(variable)
    name.variable <- pair.fact.name[2]
    print(name.variable)
    ggplot(data = data[,c("X","Y", variable)]) + 
      geom_sf(data = inputs$ch, fill="#E8E8E8", color="black") + 
      geom_sf(data = inputs$rivers.major, fill=NA, color="lightblue", show.legend = FALSE) + 
      geom_sf(data = inputs$lakes.major, fill="lightblue", color="lightblue", show.legend = FALSE) + 
      geom_point(aes(x=X, y=Y, color = data[, variable]), size= 3, alpha= 0.6) + 
      scale_colour_gradient2(name = name.variable, 
                                    high = "firebrick3") + 
      theme_void(base_size = 18) + 
      theme(panel.grid.major = element_line(colour="transparent"),
                   plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines"))
  }
  
  # !! almost there, just missing labels 
  list.plots <- lapply(list.env.fact, plot_maps_variable, inputs  = inputs, data = data)
  
  return(list.plots)
}


# Figure SI A 3: correlation matrix ####
pdf.corr.mat.env.fact <- function(data, env.fact, file.name){
  
  plot.data <- na.omit(data[,env.fact])
  plot.data <- as.matrix(plot.data)
  plot.data <- cor(plot.data)
  colnames(plot.data) <- names(env.fact)
  rownames(plot.data) <- names(env.fact)
  
  pdf(paste0(file.name, ".pdf"), paper = 'special', width = 15, height = 15, onefile = TRUE)
  corrplot(plot.data, 
           method = "number", 
           type="upper", 
           tl.col = "black",
           tl.srt=45)
  corrplot.mixed(plot.data, 
                 tl.col = "black",
                 tl.srt=45,
                 tl.pos = "lt",
                 order="hclust", 
                 lower = 'number', upper = "circle", 
                 lower.col = "black", 
                 number.cex = .9)
  dev.off()
  
  png(paste0(file.name, ".png"),  width = 1280*2/3, height = 960*2/3)
  corrplot.mixed(plot.data, 
                 tl.col = "black",
                 tl.srt=45,
                 tl.pos = "lt",
                 order="hclust", 
                 lower = 'number', upper = "circle", 
                 lower.col = "black", 
                 number.cex = .9)
  dev.off()
  cat("PDF and PNG printed")
}

# Figure SI A 9: prevalence analysis ####
analysis.prevalence <- function(prev.inv, splits, ODG){
  
  list.plots <- list()
  
  plot.data1 <- arrange(prev.inv, -desc(Prevalence))
  colnames(plot.data1)[which(colnames(plot.data1) == "Prevalence")] <- "Entire dataset"
  plot.data1$Taxa.ind <- 1:length(plot.data1$Occurrence.taxa)
  
  p <- ggplot(data=plot.data1, aes( x = plot.data1[,"Entire dataset"], y = 1:dim(prev.inv)[1], colour = Taxonomic.level))+ 
    geom_point() +
    labs(title = "Taxa prevalence in entire dataset",
         x = "Prevalence", y = "Taxa",
         color = "Taxonomic level") +
    scale_y_discrete(limits = gsub("Occurrence.", "", plot.data1$Occurrence.taxa)) +
    theme_bw(base_size = 18)
  # p
  
  list.plots[[1]] <- p
  
  if(ODG){
    
    training.set <- splits$Split1$`Training data`
    testing.set <- splits$Split1$`Testing data`
    cind.taxa <- which(grepl("Occurrence.",colnames(testing.set)))
    Occurrence.taxa <- colnames(testing.set)[cind.taxa]
    data.prev <- data.frame(Occurrence.taxa)
    data.prev[,"Training set"] <- NA
    data.prev[,"Testing set"] <- NA

    for(i in cind.taxa){
      # i <- cind.taxa[1]
      n.obs.train <- sum(!is.na(training.set[,i]))
      prev.train <- sum(na.omit(training.set[,i]))/n.obs.train # calculate the of prevalence of this taxon
      data.prev[which(data.prev[,1] == colnames(training.set)[i]), "Training set"] <- prev.train
      
      n.obs.test <- sum(!is.na(testing.set[,i]))
      prev.test <- sum(na.omit(testing.set[,i]))/n.obs.test # calculate the of prevalence of this taxon
      data.prev[which(data.prev[,1] == colnames(testing.set)[i]), "Testing set"] <- prev.test
    }
    
    plot.data2 <- data.prev
    # colnames(plot.data2)[which(colnames(plot.data2) == "Prevalence")] <- "Testing dataset"
    
    plot.data <- left_join(plot.data1, plot.data2, by = c("Occurrence.taxa"))
    plot.data <- gather(plot.data, key = Dataset, value = Prevalence, -c("Occurrence.taxa", "Taxonomic.level", "Missing.values", "Taxa.ind"))
    plot.data$Prevalence <- plot.data$Prevalence * 100
    plot.data3 <- filter(plot.data, Dataset != "Entire dataset")
    
    q <- ggplot(data=plot.data3, aes(x = Prevalence, y = Taxa.ind, colour = Dataset))+ 
      geom_point(size=4) +
      labs( # title = "Change in taxa prevalence",
           x = "Prevalence", y = "Taxa",
           color = "Dataset") +
      scale_y_discrete(limits = gsub("Occurrence.", "", plot.data1$Occurrence.taxa)) +
      theme_bw(base_size = 20)
    q
    list.plots[[2]] <- q
  }
  
  return(list.plots)
}

## ---- Models analysis plots ----

# Paper plots ####

# Extract GLMs parameters
df.glm.param <- function(outputs, df.perf, list.glm, list.taxa, env.fact.full){
  
  temp.temp.df.coef <- data.frame()
  temp.temp.comm.coef <- data.frame()
  
  for (l in list.glm) {
    # l = list.models[2]
    temp.df.coef <- data.frame(Taxa = list.taxa)
    temp.df.coef$Prevalence <- df.perf$Prevalence
    temp.df.coef[,c("Intercept", env.fact.full)] <- NA
    
    temp.comm.coef <- data.frame(variable = env.fact.full)
    
    if(l == "iGLM"){
      
      for(j in list.taxa){
        # j = list.taxa[1]
        trained.mod <- outputs[[l]][[j]][["Trained model"]]
        coef <- trained.mod[["finalModel"]][["coefficients"]]
        temp.df.coef[which(temp.df.coef$Taxa == j), c("Intercept", env.fact.full)] <- coef
      }
      
    } else {
      trained.mod <- outputs[[l]][[1]][["Trained model"]]
      res.extracted   <- rstan::extract(trained.mod,permuted=TRUE,inc_warmup=FALSE)
      ind.maxpost <- which.max(res.extracted[["lp__"]])
      alpha.taxa.maxpost <- res.extracted[["alpha_taxa"]][ind.maxpost,]
      beta.taxa.maxpost  <- res.extracted[["beta_taxa"]][ind.maxpost,,]
      names(alpha.taxa.maxpost) <- list.taxa
      colnames(beta.taxa.maxpost) <- list.taxa
      for(j in list.taxa){
        temp.df.coef[which(temp.df.coef$Taxa == j), "Intercept"] <- alpha.taxa.maxpost[j]
        temp.df.coef[which(temp.df.coef$Taxa == j), env.fact.full] <- beta.taxa.maxpost[,j]
      }
      temp.comm.coef$mu.beta.comm.maxpost  <- res.extracted[["mu_beta_comm"]][ind.maxpost,]
      temp.comm.coef$sigma.beta.comm.maxpost  <- res.extracted[["sigma_beta_comm"]][ind.maxpost,]
      temp.comm.coef$Model <- l
      temp.temp.comm.coef <- bind_rows(temp.temp.comm.coef, temp.comm.coef)
    }
    
    temp.df.coef <- gather(temp.df.coef, key = variable, value = value, -c("Taxa", "Prevalence"))
    temp.df.coef$Model <- l
    temp.temp.df.coef <- bind_rows(temp.temp.df.coef, temp.df.coef)
    
  }
  
  return(list(temp.temp.df.coef, temp.temp.comm.coef))

}

# Analysis GLMs parameters distribution
plot.glm.param <- function(outputs, df.perf, list.glm, list.taxa, env.fact.full, CV){
  
  col.vect <- names(list.glm)
  df.coef <- data.frame()
  
  # Make df coefficients
  if(CV){
    list.splits <- names(outputs)
    
    for (s in list.splits) {
      # s <- list.splits[2]
      temp.temp.df.coef <- df.glm.param(outputs[[s]], df.perf, list.glm, list.taxa, env.fact.full)[[1]]
      temp.temp.df.coef$Split <- s
      df.coef <- bind_rows(df.coef, temp.temp.df.coef)
    }
    
  } else {
    
    df.coef <- df.glm.param(outputs, df.perf, list.glm, list.taxa, env.fact.full)[[1]]
  
  }
  
  if(CV){ 
    p <- ggplot(data = df.coef, aes(value, color = Split, fill = Split))
    p <- p + geom_density(alpha = 0.1)
    p <- p + facet_grid(variable ~ Model, scales = "free")
  } else {
    p <- ggplot(data = df.coef, aes(value, color = Model, fill = Model))
    p <- p + geom_density(alpha = 0.1)
    p <- p + facet_wrap( ~ variable, scales = "free",
                         strip.position="top",
                         ncol = 2)
  }
  p <- p + scale_colour_manual(values=col.vect)
  p <- p + theme(strip.background = element_rect(fill = "white"))
  p <- p + labs(title = "Distribution of parameters", 
                x = "Parameter value",
                y = "Density")
  p <- p + theme_bw(base_size = 10)
  
  
  q <- ggplot(data = df.coef, aes(x = Prevalence, y = value, color = Model, 
                                  shape = Model))
  q <- q + geom_point(alpha = 0.7)
  if(CV){
    q <- q + facet_grid(variable ~ Model, scales = "free")
  } else {
    q <- q + facet_wrap( ~ variable, scales = "free",
                         strip.position="top",
                         ncol = 2)
  }
  q <- q + scale_colour_manual(values=col.vect, name="Model", labels=list.glm) +
    scale_shape_manual(values=c("triangle", "square", "circle"), name="Model", labels=list.glm)
  q <- q + guides(colour = FALSE)
  q <- q + theme(strip.background = element_rect(fill = "white"))
  q <- q + labs(title = "Distribution of parameters across all taxa", 
                x = "Prevalence",
                y = "Parameter value")
  q <- q + theme_bw(base_size = 10)
  
  return(list(p,q))
  
}

# Figure 1: Boxplots models performance ####
plot.boxplots.compar.appcase <- function(plot.data, list.models, models.analysis){
  
  list.models.temp <- c("#000000" = "Null_model", list.models)
  if(models.analysis["ann.hyperparam"] == T){
    list.models.temp <- gsub("ANN_", "", list.models.temp)
    plot.data$model <- gsub("ANN_", "", plot.data$model)
  }
  # Make a vector of colors
  col.vect <- names(list.models.temp)
  names(col.vect) <- list.models.temp
  
  app.case <- unique(plot.data$appcase)
  
  median.null.model <- median(plot.data[which(plot.data$model == "Null_model"), "performance"])
  median.pred.appcase1 <- vector()
  median.pred.appcase2 <- vector()
  
  for (l in list.models.temp[which(list.models.temp != "Null_model")]) {
    median.pred.appcase1[l] <- median(plot.data[which(plot.data$model == l & plot.data$dataset == "Prediction" & plot.data$appcase == app.case[1]), "performance"])
    median.pred.appcase2[l] <- median(plot.data[which(plot.data$model == l & plot.data$dataset == "Prediction" & plot.data$appcase == app.case[2]), "performance"])
  }
  median.pred.appcase1.best <- min(median.pred.appcase1)
  median.pred.appcase2.best <- min(median.pred.appcase2)
  cat("Best median stand. dev. is:\n", 
      names(median.pred.appcase1)[which(median.pred.appcase1 == median.pred.appcase1.best)], median.pred.appcase1.best, " during ", app.case[1], "\n",
      names(median.pred.appcase2)[which(median.pred.appcase2 == median.pred.appcase2.best)], median.pred.appcase2.best, " during ", app.case[2], "\n")
  plot.data.median <- data.frame("median" = c(median.pred.appcase1.best, median.pred.appcase2.best), "appcase" = app.case)
  
  plot.data.labels <- data.frame("label" = c("a)", "b)"),
                                 "dataset" = c("Calibration", "Calibration"),
                                 "model" = c("Null_model", "Null_model"),
                                 "appcase" = app.case)
  
  # Boxplots standardized deviance
  plot.data1 <- filter(plot.data, dataset %in% c("Calibration", "Prediction"))
  plot.data1 <- filter(plot.data1, model %in% list.models.temp)
  if(any(models.analysis == TRUE)){
    lev <- names(median.pred.appcase1[order(median.pred.appcase1, decreasing = T)])
    lev <- c("Null_model", lev)
  }
  p <- ggplot(plot.data1, aes(x = model, y = performance, fill = dataset))
  p <- p + geom_boxplot()
  p <- p + ylim(0, 1.5) # ECR: because perf problems
  p <- p + geom_hline(yintercept = median.null.model, linetype='longdash', col = 'black')
  p <- p + geom_hline(data = plot.data.median, aes(yintercept = median), 
                      linetype='dotted', col = 'black')
  p <- p + scale_fill_manual(values=c(Calibration = "#998ec3", Prediction = "#f1a340"))
  if(any(models.analysis == TRUE)){
    p <- p + scale_x_discrete(limits = lev)
    p <- p + coord_flip()
  } else {
    p <- p + scale_x_discrete(limits = list.models.temp)
  }
  if(length(app.case) != 1){
    p <- p + geom_text(data = plot.data.labels, y = ifelse(any(models.analysis == T), 0.2, 1.4), label = plot.data.labels$label, size = 7, fontface="bold")
  }

  p <- p + theme_bw(base_size = 25)
  p <- p + facet_wrap(~ appcase, nrow = 2, # scales = "free_y", 
                      # labeller=label_parsed, 
                      strip.position= ifelse(any(models.analysis == T), "top", "right"))
  p <- p + theme(legend.text = element_text(size=22),
                 strip.background = element_rect(fill = "white"))
  p <- p + labs(x="Model",
                  y="Standardized deviance",
                  fill = "",
                  # title = "Model performance comparison")
                  title = "")
  # p
  
  # Boxplots auc
  plot.data2 <- filter(plot.data, dataset %in% c("AUC Calibration", "AUC Prediction"))
  plot.data2 <- filter(plot.data2, model %in% list.models.temp)
  q <- ggplot(plot.data2, aes(x = model, y = performance, fill = dataset))
  q <- q + geom_boxplot()
  q <- q + scale_x_discrete(limits = list.models.temp[which(list.models.temp != "Null_model")])
  q <- q + ylim(0.3, 1.1) # ECR: because perf problems
  q <- q + geom_hline(yintercept = 0.5, linetype='longdash', col = 'black')
  q <- q + scale_fill_manual(values=c("#998ec3", "#f1a340"))
  q <- q + theme_bw(base_size = 25)
  q <- q + facet_wrap(~ appcase, nrow = 2, # scales = "free_y", 
                      # labeller=label_parsed, 
                      strip.position = ifelse(any(models.analysis == T), "top", "right"))
  q <- q + theme(legend.text = element_text(size=22),
                 strip.background = element_rect(fill = "white"))
  if(any(models.analysis == TRUE)){
    q <- q + coord_flip()
  }
  q <- q + labs(x="Model",
                y="AUC",
                fill = "",
                # title = "Model performance comparison")
                title = "")
  # q
  
  list.plots <- list(p,q)
  
  return(list.plots)
}

# Figure 2: Models performance vs taxa prevalence ####
plot.perfvsprev.compar.appcase <- function(plot.data, list.models, list.taxa){
  
  list.models.temp <- list.models
  col.vect.temp <- names(list.models)
  names(col.vect.temp) <- list.models
  list.models <- c("#000000" = "Null_model", list.models)
  ex.taxa <- c("Psychodidae", "Gammaridae")
  
  # Make a vector of colors
  col.vect <- names(list.models)
  names(col.vect) <- list.models
  app.case <- unique(plot.data$appcase)
  no.taxa <- length(list.taxa)
  
  plot.data$Prevalence <- plot.data$Prevalence*100
  
  # Prevalence vs stand dev
  plot.data1 <- filter(plot.data, dataset %in% c("Calibration", "Prediction"))
  plot.data1 <- filter(plot.data1, model %in% list.models)
  plot.data.labels <- data.frame("label" = c("a)", "b)", "c)", "d)"),
                                 "dataset" = c("Calibration", "Prediction", "Calibration", "Prediction"),
                                 "appcase" = c(app.case[1], app.case[1], app.case[2], app.case[2]))
  plot.data.select.taxa <- data.frame("taxon" = ex.taxa,
                                      "dataset" = c("Prediction","Prediction"),
                                      "appcase" = c(app.case[1], app.case[1]))
  p <- ggplot()
  p <- p  + geom_point(data = plot.data1, aes_string(x = "Prevalence", y = "performance", 
                                                      colour = "model"), alpha = 0.8,
                       size = 3)
  p <- p + geom_text(data = plot.data.select.taxa,
                     x = c(55, 73), y = c(1.55, 0.45), label = ex.taxa, size = 7, fontface="italic")
  p <- p + xlim(4.5, 95.5)
  p <- p + ylim(0, 1.6) # ECR: only because perf problems
  p <- p + geom_hline(yintercept = plot.data[select.taxa[1], "Prevalence"], linetype='dashed', col = 'grey30')
  p <- p + stat_function(fun=function(x) -2*(x/100*log(x/100) + (1-x/100)*log(1-x/100))) # to plot null model as function line
  p <- p + theme_bw(base_size = 20)
  if(length(app.case) != 1){
    p <- p + facet_grid(appcase ~ dataset)
    p <- p + geom_text(data = plot.data.labels,
                       x = c(10, 90, 10, 90), y = 1.5, label = plot.data.labels$label, size = 7, fontface="bold")
  } else {
    p <- p + facet_wrap(~ dataset)
  }
  p <- p + theme(legend.title = element_text(size=22),
               legend.text = element_text(size=20),
               strip.background = element_rect(fill="white"))
  p <- p + labs(y = "Standardized deviance",
                 x = "Prevalence (%)",
                 # shape = "Taxonomic level",
                 color = "Model",
                 title = "")
  p <- p + scale_colour_manual(values=col.vect)
  p <- p + guides(colour = guide_legend(override.aes = list(size=6)))
  # p
  
  # Prevalence vs AUC
  plot.data2 <- filter(plot.data, dataset %in% c("AUC Calibration", "AUC Prediction"))
  plot.data2 <- filter(plot.data2, model %in% list.models.temp)
  q <- ggplot()
  q <- q  + geom_point(data = plot.data2, aes_string(x = "Prevalence", y = "performance", 
                                                     colour = "model"), alpha = 0.8,
                       size = 3)
  q <- q + xlim(4.5, 95.5)
  q <- q + ylim(0.3, 1.1) # ECR: only because perf problems
  q <- q + geom_hline(yintercept = 0.5, linetype='dashed', col = 'grey30')
  q <- q + theme_bw(base_size = 20)
  q <- q + facet_grid(appcase ~ dataset # , scales = "free"
  )
  q <- q + theme(legend.title = element_text(size=22),
                 legend.text = element_text(size=20),
                 strip.background = element_rect(fill="white"))
  q <- q + labs(y = "AUC",
                x = "Prevalence (%)",
                # shape = "Taxonomic level",
                color = "Model",
                title = "")
  q <- q + scale_colour_manual(values=col.vect.temp)
  q <- q + guides(colour = guide_legend(override.aes = list(size=6)))
  # p
  
  list.plots <- list(p,q)
  
  return(list.plots)
}


# Figure 3, SI B and SI C: ICE and PDP ####
plot.ice.per.taxa <- function(ice.plot.data.taxon, list.models, subselect, ice.random){
  
  # ice.plot.data.taxon <- list.ice.plot.data[[1]]
  plot.ice.per.env.fact = function (plot.data.env.fact, list.models, subselect, ice.random) {
    # plot.data.env.fact <- ice.plot.data.taxon[[1]]
    # subselect <- 1
    # ice.random <- T
    name.taxon <- plot.data.env.fact$name.taxon
    name.fact <- plot.data.env.fact$env.fact
    plot.data <- plot.data.env.fact$plot.data.ice
    plot.data.mean.pdp <- plot.data.env.fact$plot.data.mean.pdp
    plot.data.min.max.pdp <- plot.data.env.fact$plot.data.min.max.pdp
    plot.data.means <- plot.data.env.fact$plot.data.means
    plot.data.rug <- plot.data.env.fact$plot.data.rug
    x.coord.arrow <- max(plot.data$variable)
    
    # Subselect predictions for plot resolution
    temp.range.fact <- unique(plot.data$variable)[seq(1,length(unique(plot.data$variable)), by = subselect)]
    plot.data <- filter(plot.data, plot.data$variable %in% temp.range.fact)
    plot.data.mean.pdp <- filter(plot.data.mean.pdp, plot.data.mean.pdp$variable %in% temp.range.fact)
    plot.data.means <- filter(plot.data.means, plot.data.means$variable %in% temp.range.fact)
    
    # Select seed if no analysis of randomness
    if(!ice.random){
      plot.data <- filter(plot.data, plot.data$Seed == 2021)
      plot.data.mean.pdp <- filter(plot.data.mean.pdp, plot.data.mean.pdp$Seed == 2021)
      plot.data.means <- filter(plot.data.means, plot.data.means$Seed == 2021)
    }
    
    # Reorder model factor levels to have right order on panel
    plot.data <- plot.data %>%
      mutate(across(Model, factor, levels=list.models))
    plot.data.means <- plot.data.means %>%
      mutate(across(Model, factor, levels=list.models))
    plot.data.mean.pdp <- plot.data.mean.pdp %>%
      mutate(across(Model, factor, levels=list.models))
    plot.data.rug <- plot.data.rug %>%
      mutate(across(Model, factor, levels=list.models))
    
    p <- ggplot(plot.data, aes(x = variable, y = value, group=factor(observation))) 
    p <- p + geom_line(aes(color=factor(observation)), alpha = 0.3, show.legend = FALSE) # remove legend that would be the number of samples
    p <- p + geom_line(data = plot.data.mean.pdp,
                       aes(x = variable, y = mean.pdp), 
                       color = "grey20", size = 1.5, # alpha = 0.7, 
                       inherit.aes = F)
    p <- p + geom_line(data = plot.data.means,
                       aes(x = variable, y = mean.average), 
                       color = "grey20", size = 0.5, linetype='dashed',# alpha = 0.7, 
                       inherit.aes = F)
    p <- p + geom_rug(data = plot.data.rug,
                      aes(x = variable), 
                      color = "grey20", alpha = 0.7, inherit.aes = F)
    p <- p + geom_segment(data =  plot.data.min.max.pdp,inherit.aes = FALSE,
                          lineend = "round", linejoin = "round", #linetype = "dashed",
                          aes(x = x.coord.arrow, y = min, xend = x.coord.arrow, yend = max),
                          arrow = arrow(length = unit(0.3, "cm"), ends = "both"))
    if(ODG & name.fact == "Temperature"){
      split.value <- max(plot.data.rug[,"variable"])
      p <- p + geom_vline(xintercept = split.value, linetype='dashed', col = 'grey30')
    }
    if(ice.random){
      p <- p + facet_grid(Model ~ Seed)
    } else {
      p <- p + facet_wrap( ~ Model,# scales = "free_x", 
                           strip.position="top",
                           ncol = 2)
    }
    p <- p + theme_bw(base_size = 10)
    p <- p + theme(strip.background = element_rect(fill = "white"))
    p <- p + labs(title = paste("ICE and PDP for", name.taxon), # paste(l),
                  # subtitle = paste("Resolution:", no.steps/no.subselect, "steps"),
                  x = name.fact,
                  y = "Predicted probability of occurrence")
    # p
    return(p)
  }
  
  list.plots <- lapply(ice.plot.data.taxon, plot.ice.per.env.fact, list.models, subselect, ice.random)
  
  return(list.plots)
}

# Models prediction vs observed presence/absence
plot.rs.taxa <- function(taxa, outputs, list.models, normalization.data, env.fact, CV, ODG){
  
  # taxa <- subselect.taxa[1]
  # taxon <- sub("Occurrence.", "", taxa)
  cat("\nProducing response shape plot for:", sub("Occurrence.", "", taxa), "\n")
  
  # temp.list.models <- c("Observations", list.models)
  temp.list.models <- list.models
  vect.env.fact <- unname(env.fact)
  names.env.fact <- names(env.fact)
  
  if(ODG){
    normalization.data <- normalization.data[[1]]
  }
  # Make list of plots per fact to be returned at the end
  list.plots <- vector(mode = 'list', length = length(env.fact))
  names(list.plots) <- names.env.fact
  
  plot.data0 <- data.frame()
  
  for (l in temp.list.models) {
    # l <- temp.list.models[1]
    # l <- list.models[1]
    
    cat("For model", l, "\n")
    m <- "training set"
    
    if(l == "Observations"){
      plot.data1 <- outputs[["iGLM"]][[taxa]][[paste("Observation", m)]]
      plot.data1$pred <- ifelse(plot.data1[,taxa] == "present", 1, 0)
    } else if(l == "hGLM" | l == "chGLM") { 
      plot.data1 <- outputs[["iGLM"]][[taxa]][[paste("Observation",m)]]
      plot.data1$pred <- outputs[[l]][[taxa]][[paste("Prediction probabilities", m)]][,"present"]
    } else { 
      plot.data1 <- outputs[[l]][[taxa]][[paste("Observation",m)]]
      plot.data1$pred <- outputs[[l]][[taxa]][[paste("Prediction probabilities", m)]][,"present"]
    }
    plot.data1[,vect.env.fact] <- as.data.frame(sweep(sweep(plot.data1[,vect.env.fact], 2, normalization.data$SD[names.env.fact], FUN="*"), 2, normalization.data$Mean[names.env.fact], FUN = "+"))
    plot.data1 <- gather(plot.data1, key = factor, value = value, -SiteId, -SampId, -X, -Y, -taxa, -pred)
    
    if(CV | ODG){
      plot.data1$set <- "Training"
      
      m <- "testing set"
      if(l == "Observations"){
        plot.data2 <- outputs[["iGLM"]][[taxa]][[paste("Observation", m)]]
        plot.data2$pred <- ifelse(plot.data2[,taxa] == "present", 1, 0)
      } else if(l == "hGLM" | l == "chGLM") { 
        plot.data2 <- outputs[["iGLM"]][[taxa]][[paste("Observation",m)]]
        plot.data2$pred <- outputs[[l]][[taxa]][[paste("Prediction probabilities", m)]][,"present"]
      } else { 
        plot.data2 <- outputs[[l]][[taxa]][[paste("Observation", m)]]
        plot.data2$pred <- outputs[[l]][[taxa]][[paste("Prediction probabilities", m)]][,"present"]
      }
      plot.data2[,vect.env.fact] <- as.data.frame(sweep(sweep(plot.data2[,vect.env.fact], 2, normalization.data$SD[names.env.fact], FUN="*"), 2, normalization.data$Mean[names.env.fact], FUN = "+"))
      plot.data2 <- gather(plot.data2, key = factor, value = value, -SiteId, -SampId, -X, -Y, -taxa, -pred)
      plot.data2$set <- "Testing"
      
      plot.data3 <- bind_rows(plot.data1, plot.data2)
    } else {
      plot.data3 <- plot.data1
    }
    plot.data3$model <- l
    
    plot.data0 <- bind_rows(plot.data0, plot.data3)
  }
  
  for (k in env.fact) {
    # k <- env.fact[1]
    env.ind <- which(env.fact == k)
    plot.data <- filter(plot.data0, factor == k)
    if(ODG){
      plot.data$set <- as.factor(plot.data$set)
    } else {
      plot.data$set <- as.factor("Calibration")
    }
    
    plot.data <- plot.data %>%
      mutate(across(.cols = c(model, factor), levels=temp.list.models))
    # plot.data$model <- as.factor(plot.data$model)
    
    g <- ggplot(data = plot.data, aes_string(x = "value", y = "pred", color = taxa)) 
    # g <- ggplot(data = plot.data, aes(x = "value", y = "pred")) + geom_point()
    
    if(CV | ODG){
      g <- g + geom_point(aes(shape = set), alpha = 0.35)
      g <- g + labs(shape = "Dataset")
      g <- g + guides(shape = guide_legend(override.aes = list(size=3)))
    } else {
      g <- g + geom_point(alpha = 0.35)
    }
    g <- g + theme_bw(base_size=10)
    g <- g + facet_wrap( ~ model, scales = "free_x", 
                         #labeller=label_parsed, 
                         strip.position="top",
                         ncol = 2)
    if(ODG & k == "temperature"){
      split.value <- max(plot.data[which(plot.data[,"set"] == "Training"),"value"])
      g <- g + geom_vline(xintercept = split.value, linetype='dashed', col = 'grey30')
    }
    g <- g + theme(strip.background = element_rect(fill = "white"))
    g <- g + scale_color_manual(name = "Observation", values=c(absent = "#c2141b", present = "#007139"), labels = c("Absence", "Presence"))
    g <- g + labs( # title = paste(names(k), "-",paste(taxon)),
      # title = paste("Probability of occurrence vs explanatory variables"),
      # subtitle = paste(k, "-",paste(taxon)),
      x = k,
      y = "Predicted probability of occurrence",
      color = "Observation")
    g <- g + ylim(0,1)
    g <- g + guides(colour = guide_legend(override.aes = list(size=6)))
    # g
    list.plots[[env.ind]] <- g
  }
  
  return(list.plots)
  
}