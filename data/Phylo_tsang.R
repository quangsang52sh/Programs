library(shiny)
library(ape)
library(phangorn)

run_ML <- function(dat, boostrap) {
#Input modeltest Function 
aic.weights <- function(aic){
  diff.aic <- aic-min(aic)
  exp(-0.5 * diff.aic) / sum(exp(-0.5 * diff.aic))
}

modelTest <- function (object, tree = NULL, model = "all", G = TRUE, I = TRUE, FREQ=FALSE, k = 4, 
                       control = pml.control(epsilon = 1e-08, maxit = 100, trace = 1), 
                       multicore = FALSE, mc.cores = NULL) ##c("JC", "F81", "K80","HKY", "SYM", "GTR")
{    
  #    multicore <- mc.cores > 1L
  if(multicore && is.null(mc.cores)){
    mc.cores <- detectCores()
  }
  if (inherits(object,"phyDat")) 
    data <- object
  if (inherits(object,"pml")) {
    data <- object$data
    if (is.null(tree)) 
      tree <- object$tree
  }
  
  if(attr(data, "type")=="DNA") type <- c("JC", "F81", "K80", "HKY", "TrNe", 
                                          "TrN", "TPM1", "K81", "TPM1u", "TPM2", "TPM2u", "TPM3", "TPM3u", 
                                          "TIM1e", "TIM1", "TIM2e", "TIM2", "TIM3e", "TIM3", "TVMe", "TVM", 
                                          "SYM", "GTR")
  if(attr(data, "type")=="AA") type <- .aamodels     
  if( (length(model)==1) && model == "all") model <- type
  model <- match.arg(model, type, TRUE)  
  env <- new.env()
  assign("data", data, envir=env) 
  if (is.null(tree)) 
    tree <- NJ(dist.hamming(data))
  else{
    if(length(tree$tip.label) > 3) tree <- nnls.phylo(tree, dist.ml(data)) 
    # may need something faster for trees > 500 taxa  
  }
  trace <- control$trace
  control$trace <- trace - 1
  fit <- pml(tree, data)
  fit <- optim.pml(fit, control = control)
  l <- length(model)
  if(attr(fit$data, "type")=="DNA")FREQ <- FALSE    
  n <- 1L + sum(I + G + (G & I) + FREQ + (FREQ & I) + (FREQ & G) + 
                  (FREQ & G & I))
  nseq <- sum(attr(data, "weight"))  
  fitPar <- function(model, fit, G, I, k, FREQ) {
    m <- 1
    res <- matrix(NA, n, 6)
    res <- as.data.frame(res)
    colnames(res) <- c("Model", "df", "logLik", "AIC", "AICc", "BIC")
    data.frame(c("Model", "df", "logLik", "AIC", "AICc", "BIC"))
    calls <- vector("list", n)
    trees <- vector("list", n)
    fittmp <- optim.pml(fit, model = model, control = control)
    res[m, 1] <- model
    res[m, 2] <- fittmp$df
    res[m, 3] <- fittmp$logLik
    res[m, 4] <- AIC(fittmp)
    res[m, 5] <- AICc(fittmp)
    res[m, 6] <- AIC(fittmp, k = log(nseq))
    calls[[m]] <- fittmp$call  
    trees[[m]] <- fittmp$tree
    m <- m + 1
    if (I) {
      if(trace>0)print(paste0(model, "+I"))
      fitI <- optim.pml(fittmp, model = model, optInv = TRUE, 
                        control = control)
      res[m, 1] <- paste0(model, "+I")
      res[m, 2] <- fitI$df
      res[m, 3] <- fitI$logLik
      res[m, 4] <- AIC(fitI)
      res[m, 5] <- AICc(fitI)
      res[m, 6] <- AIC(fitI, k = log(nseq))
      calls[[m]] <- fitI$call
      trees[[m]] <- fitI$tree
      m <- m + 1
    }
    if (G) {
      if(trace>0)print(paste0(model, "+G"))
      fitG <- update(fittmp, k = k)
      fitG <- optim.pml(fitG, model = model, optGamma = TRUE, 
                        control = control)
      res[m, 1] <- paste0(model, "+G")
      res[m, 2] <- fitG$df
      res[m, 3] <- fitG$logLik
      res[m, 4] <- AIC(fitG)
      res[m, 5] <- AICc(fitG)
      res[m, 6] <- AIC(fitG, k = log(nseq))
      calls[[m]] <- fitG$call
      trees[[m]] <- fitG$tree
      m <- m + 1
    }
    if (G & I) {
      if(trace>0)print(paste0(model, "+G+I"))
      fitGI <- update(fitI, k = k)
      fitGI <- optim.pml(fitGI, model = model, optGamma = TRUE, 
                         optInv = TRUE, control = control)
      res[m, 1] <- paste0(model, "+G+I")
      res[m, 2] <- fitGI$df
      res[m, 3] <- fitGI$logLik
      res[m, 4] <- AIC(fitGI)
      res[m, 5] <- AICc(fitGI)
      res[m, 6] <- AIC(fitGI, k = log(nseq))
      calls[[m]] <- fitGI$call
      trees[[m]] <- fitGI$tree
      m <- m + 1
    }
    if (FREQ) {
      if(trace>0)print(paste0(model, "+F"))
      fitF <- optim.pml(fittmp, model = model, optBf = TRUE, 
                        control = control)
      res[m, 1] <- paste0(model, "+F")
      res[m, 2] <- fitF$df
      res[m, 3] <- fitF$logLik
      res[m, 4] <- AIC(fitF)
      res[m, 5] <- AICc(fitF)
      res[m, 6] <- AIC(fitF, k = log(nseq))
      calls[[m]] <- fitF$call
      trees[[m]] <- fitF$tree
      m <- m + 1
    }
    if (FREQ & I) {
      if(trace>0)print(paste0(model, "+I+F"))
      fitIF <- update(fitF, inv = fitI$inv)
      fitIF <- optim.pml(fitIF, model=model, optBf = TRUE, optInv = TRUE,
                         control = control)
      res[m, 1] <- paste0(model, "+I+F")
      res[m, 2] <- fitIF$df
      res[m, 3] <- fitIF$logLik
      res[m, 4] <- AIC(fitIF)
      res[m, 5] <- AICc(fitIF)
      res[m, 6] <- AIC(fitIF, k = log(nseq))
      calls[[m]] <- fitIF$call
      trees[[m]] <- fitIF$tree
      m <- m + 1
    }
    if (FREQ & G) {
      if(trace>0)print(paste0(model, "+G+F"))
      fitGF <- update(fitF, k=k, shape=fitG$shape)
      fitGF <- optim.pml(fitGF, model = model, optBf = TRUE, 
                         optGamma = TRUE, control = control)
      res[m, 1] <- paste0(model, "+G+F")
      res[m, 2] <- fitGF$df
      res[m, 3] <- fitGF$logLik
      res[m, 4] <- AIC(fitGF)
      res[m, 5] <- AICc(fitGF)
      res[m, 6] <- AIC(fitGF, k = log(nseq))
      calls[[m]] <- fitGF$call
      trees[[m]] <- fitGF$tree
      m <- m + 1
    }
    if (FREQ & G & I) {
      if(trace>0)print(paste0(model, "+G+I+F"))
      fitGIF <- update(fitIF, k=k)
      fitGIF <- optim.pml(fitGIF, model = model, optBf = TRUE, 
                          optInv = TRUE, optGamma = TRUE, control = control)
      res[m, 1] <- paste0(model, "+G+I+F")
      res[m, 2] <- fitGIF$df
      res[m, 3] <- fitGIF$logLik
      res[m, 4] <- AIC(fitGIF)
      res[m, 5] <- AICc(fitGIF)
      res[m, 6] <- AIC(fitGIF, k = log(nseq))
      calls[[m]] <- fitGIF$call
      trees[[m]] <- fitGIF$tree
      m <- m + 1
    }
    list(res, trees, calls)
  }
  eval.success <- FALSE
  if (!eval.success & multicore) {
    # !require(parallel) ||         
    #        if (.Platform$GUI != "X11") {
    #            warning("package 'parallel' not found or GUI is used, \n      analysis is performed in serial")
    #       }
    #        else {
    RES <- mclapply(model, fitPar, fit, G, I, k, FREQ, mc.cores=mc.cores)
    eval.success <- TRUE
    #        }
  }
  if (!eval.success) 
    RES <- lapply(model, fitPar, fit, G, I, k, FREQ)
  RESULT <- matrix(NA, n * l, 8)
  RESULT <- as.data.frame(RESULT)
  colnames(RESULT) <- c("Model", "df", "logLik", "AIC", "AICw", "AICc", 
                        "AICcw", "BIC")
  
  for (i in 1:l) RESULT[((i - 1) * n + 1):(n * i), c(1,2,3,4,6,8)] <- RES[[i]][[1]]
  RESULT[,5] <- aic.weights(RESULT[,4])
  RESULT[,7] <- aic.weights(RESULT[,6])
  for(i in 1:l){
    for(j in 1:n){
      mo <- RES[[i]][[1]][j,1]
      tname <- paste0("tree_", mo)
      tmpmod <- RES[[i]][[3]][[j]]
      tmpmod["tree"] <- call(tname)
      if(!is.null(tmpmod[["k"]]))tmpmod["k"] <- k
      if(attr(data, "type")=="AA") tmpmod["model"] <- RES[[i]][[1]][1,1]          
      assign(tname, RES[[i]][[2]][[j]], envir=env)
      assign(mo, tmpmod, envir=env) 
    }
  }
  attr(RESULT, "env") <- env 
  RESULT
}

tidy.modelTest <- function(x){
  env <- attr(x, "env")
  l <- nrow(x)
  k <- rep(1L, l)
  shape <- rep(NA_real_, l)
  inv <- rep(0, l)
  for(i in seq_len(l)){
    tmp <- get(x$Model[i], env)
    if(!is.null(tmp[["k"]]))k[i] <- tmp[["k"]]
    if(!is.null(tmp[["shape"]]))shape[i] <- tmp[["shape"]]
    if(!is.null(tmp[["inv"]]))inv[i] <- tmp[["inv"]]
  }
  data.frame(Model = x$Model, k=k, shape=shape, inv=inv)
}


# Model selection using ML
model <- modelTest(dat)
write.csv(model, file = "MLmodel.csv")
env <- attr(model, "env")
ls(env=env)
  #best_model <- model[which.min(model$AIC), "Model"]
  best_model <- model[which.min(model$AIC), "Model"]
  getModel <- get(best_model, env)
  your_model <- eval( getModel, env=env)
  print(your_model)
  save(your_model, file = paste0(best_model,".Rdata"))
  #getModel <- get(best_model, env)
  dm <- dist.ml(dat)
  treeNJ <- NJ(dm)
  fit <- pml(treeNJ, data = dat)
  # Extract the prefix e.g. "GTR" from getModel
  modelPrefix <- substr(best_model, 1, regexpr("\\+", best_model) - 1)
  m1=paste0(modelPrefix, "+G+I")
  m2=paste0(modelPrefix, "+G")
  m3=paste0(modelPrefix, "+I")
  if (best_model == m1) {
    fit1 <- optim.pml(fit, model = modelPrefix, optInv = TRUE, optGamma = TRUE, rearrangement = "none", control = pml.control(trace = 0))
  } else if (best_model == m2) {
    fit1 <- optim.pml(fit, model = modelPrefix, optInv = FALSE, optGamma = TRUE, rearrangement = "none", control = pml.control(trace = 0))
  } else if (best_model == m3) {
    fit1 <- optim.pml(fit, model = modelPrefix, optInv = TRUE, optGamma = FALSE, rearrangement = "none", control = pml.control(trace = 0))
  } else {
    fit1 <- optim.pml(fit, model = best_model, optInv = FALSE, optGamma = FALSE, rearrangement = "none", control = pml.control(trace = 0))
  }
  bs <- bootstrap.pml(fit1, bs = boostrap, optNni = TRUE, control = pml.control(trace = 2))
  tree <- plotBS(midpoint(fit1$tree), bs, p = 50, type = "p")
  write.tree(tree, "MLtree.tre")
  save.image("runML.Rdata")
}

run_ML_withoutModel <- function(dat, modelname, bootstrap) {
  best_model <- modelname  # Get model from input
  dm <- dist.ml(dat)
  treeNJ <- NJ(dm)
  fit <- pml(treeNJ, data = dat)

  # Extract base model name (e.g., "GTR" from "GTR+G+I")
  modelPrefix <- unlist(strsplit(best_model, "\\+"))[1]

  # Define model variations
  m1 <- paste0(modelPrefix, "+G+I")
  m2 <- paste0(modelPrefix, "+G")
  m3 <- paste0(modelPrefix, "+I")

  # Run optimization based on model
  if (best_model == m1) {
    fit1 <- optim.pml(fit, model = modelPrefix, optInv = TRUE, optGamma = TRUE, 
                      rearrangement = "none", control = pml.control(trace = 0))
  } else if (best_model == m2) {
    fit1 <- optim.pml(fit, model = modelPrefix, optInv = FALSE, optGamma = TRUE, 
                      rearrangement = "none", control = pml.control(trace = 0))
  } else if (best_model == m3) {
    fit1 <- optim.pml(fit, model = modelPrefix, optInv = TRUE, optGamma = FALSE, 
                      rearrangement = "none", control = pml.control(trace = 0))
  } else {
    fit1 <- optim.pml(fit, model = modelPrefix, optInv = FALSE, optGamma = FALSE, 
                      rearrangement = "none", control = pml.control(trace = 0))
  }

  # Run bootstrap analysis
  bs <- bootstrap.pml(fit1, bs = bootstrap, optNni = TRUE, control = pml.control(trace = 2))

  # Generate tree
  tree <- plotBS(midpoint(fit1$tree), bs, p = 50, type = "p")
  write.tree(tree, "MLtree_withoutModel.tre")
  save.image("runML_withoutModel.Rdata")
}


run_MrBayes <- function(dat) {
  sys_name <- Sys.info()["sysname"]
  mrbayes_path <- ""
  output_dir <- paste0(getwd(), "/Mrbayes_results")

  # Ensure output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }

  if (sys_name == "Windows") {
    message("Downloading MrBayes for Windows...")
    win_url <- "https://github.com/NBISweden/MrBayes/releases/download/v3.2.7/MrBayes-3.2.7-WIN.zip"
    zip_file <- "MrBayes-3.2.7-WIN.zip"
    extract_folder <- "MrBayes-3.2.7-WIN"

    # Download if not exists
    if (!file.exists(zip_file)) {
      download.file(win_url, zip_file, mode = "wb")
    }

    # Extract only if not extracted
    if (!dir.exists(extract_folder)) {
      unzip(zip_file, exdir = extract_folder)
    }

    # Ensure correct path to executable
    mrbayes_path <- file.path(getwd(), extract_folder, "bin", "mb.3.2.7-win64.exe")

    # Check if executable exists
    if (!file.exists(mrbayes_path)) {
      stop("❌ Error: MrBayes executable not found! Please check extraction.")
    }
    
  } else if (sys_name == "Linux") {
    message("Installing MrBayes for Linux...")
    if (!file.exists("./MrBayes/src/mb")) {
      system("git clone --depth=1 https://github.com/NBISweden/MrBayes.git")
      system("cd MrBayes && ./configure && make")
    }
    mrbayes_path <- "./MrBayes/src/mb"

  } else if (sys_name == "Darwin") { # macOS
    message("Downloading MrBayes for macOS...")
    mac_url <- "https://github.com/NBISweden/MrBayes/releases/download/v3.2.6/MrBayes-3.2.6_MACx64.zip"
    zip_file <- "MrBayes-3.2.6_MACx64.zip"

    if (!file.exists(zip_file)) {
      download.file(mac_url, zip_file, mode = "wb")
    }

    if (!dir.exists("MrBayes")) {
      unzip(zip_file, exdir = "MrBayes")
    }

    mrbayes_path <- file.path(getwd(), "MrBayes", "mb")

  } else {
    stop("❌ Unsupported operating system.")
  }

  # Ensure MrBayes can run
  if (!file.exists(mrbayes_path)) {
    stop("❌ Error: MrBayes executable not found! Please check installation.")
  }

  # Convert input data to Nexus format
  bin <- as.DNAbin(dat)
  nexus_file <- "mrbayes_input.nexus"
  write.nexus.data(bin, nexus_file)

  # Download autoRunMrbayes.txt if missing
  auto_file <- "autoRunMrbayes.txt"
  auto_url <- "https://raw.githubusercontent.com/quangsang52sh/Programs/main/data/autoRunMrbayes.txt"

  if (!file.exists(auto_file)) {
    message("Downloading MrBayes script...")
    download.file(auto_url, auto_file, mode = "wb")
  }

  # Run MrBayes
  message("Running MrBayes...")
  system(paste(shQuote(mrbayes_path), "<", auto_file, "> mrbayes_input_log.txt"))

  # Verify if MrBayes ran successfully
  log_file <- "mrbayes_input_log.txt"
  if (file.exists(log_file)) {
    log_contents <- readLines(log_file, n = 20)
    if (any(grepl("error", tolower(log_contents)))) {
      stop("❌ Error: MrBayes encountered an issue. Check 'mrbayes_input_log.txt'.")
    }
  } else {
    stop("❌ Error: MrBayes did not create a log file. Please check execution.")
  }

  # Move output files to result directory
  tree_file <- paste0("mrbayes_input.nexus.con.tre")
  run1_file <- paste0("mrbayes_input.nexus.run1.t")
  run2_file <- paste0("mrbayes_input.nexus.run2.t")

  if (!file.exists(tree_file) || !file.exists(run1_file) || !file.exists(run2_file)) {
    stop("❌ Error: MrBayes output files were not generated.")
  }

  file.copy(c(tree_file, run1_file, run2_file), output_dir, overwrite = TRUE)

  message("✅ MrBayes analysis completed successfully!")
}


run_MP <- function(dat, bootstrap) {
  # Maximum Parsimony (MP) tree
  MegaMPTree <- pratchet(dat)
  MegaMPTree <- acctran(MegaMPTree, dat)
  bs_megaMP <- bootstrap.phyDat(dat, pratchet, bs = bootstrap)
  MPTree <- plotBS(midpoint(MegaMPTree), bs_megaMP, p = 50, type = "p")
  add.scale.bar()
  write.tree(MPTree, "MPtree.tre")
  save.image("runMP.Rdata")
}

run_NJ <- function(dat, bootstrap) {
  tree <- nj(dist.ml(dat))
  MegaNJTree <- tree
  bs_megaNJ <- bootstrap.phyDat(dat, FUN = function(x) NJ(dist.hamming(x)), bs = bootstrap)
  treeNJ <- plotBS(midpoint(MegaNJTree), bs_megaNJ, p = 50, type = "p")
  write.tree(treeNJ, "NJtree.tre")
  save.image("runNJ.Rdata")
}

runPhylogeneticAnalysis <- function(file, NJ = FALSE, MP = FALSE, ML = FALSE, Mrbayes = FALSE, run_ML_withoutModel = FALSE, modelname, bootstrap) {
  withProgress(message = "Running analysis...", value = 0, {
    setProgress(0)
    
    if (NJ) {
      myfile <- read.dna(file$datapath, format = "fasta")
      dat <- as.phyDat(myfile)
      run_NJ(dat, bootstrap)
    }

    if (MP) {
      myfile <- read.dna(file$datapath, format = "fasta")
      dat <- as.phyDat(myfile)
      run_MP(dat, bootstrap)
    }

    if (ML) {
      myfile <- read.dna(file$datapath, format = "fasta")
      dat <- as.phyDat(myfile)
      run_ML(dat, bootstrap)
    }

    if (Mrbayes) {
      myfile <- read.dna(file$datapath, format = "fasta")
      dat <- as.phyDat(myfile)
      run_MrBayes(dat)
    }


    if (run_ML_withoutModel) {
      myfile <- read.dna(file$datapath, format = "fasta")
      dat <- as.phyDat(myfile)
      run_ML_withoutModel(dat, modelname, bootstrap)
    }
    
    setProgress(1)
    Sys.sleep(0.5)  # Add a small delay for visual feedback
  })
}

ui <- fluidPage(
  tags$head(
    # Adjust window size on load
    tags$script(HTML("
      window.onload = function() {
        window.resizeTo(screen.width * 0.8, screen.height * 0.8);
      }
    ")),
    
    # CSS Styling
    tags$style(HTML("
      /* Title Box */
      .title-box {
        font-size: 25px;
        background-color: #f0f0f0;
        padding: 15px;
        border-radius: 10px;
        border: 2px solid #d3d3d3;
        box-shadow: 2px 2px 10px rgba(0, 0, 0, 0.2);
        margin-bottom: 10px;
      }
      /* Title Box Row */
      .title-row {
        display: flex; 
        align-items: center; 
        justify-content: flex-start;
      }
      /* DNA Image in Title */
      .dna-header {
        width: 80px;   /* Adjust size as needed */
        height: auto; 
        margin-right: 15px;
      }
      /* General Body Styling */
      body {
        background-color: white;
        color: blue;
      }
      /* Profile Picture */
      .profile-pic {
        display: flex;
        justify-content: center;
        margin: 10px 0;
      }
      .profile-pic img {
        border-radius: 50%;
        border: 3px solid #dddddd;
        width: 150px;
        height: 150px;
        object-fit: cover;
        box-shadow: 3px 3px 10px rgba(0, 0, 0, 0.2);
      }
      /* Welcome Box */
      .welcome-box {
        background-color: #f5f5f5; 
        padding: 20px; 
        border: 2px solid #d3d3d3; 
        border-radius: 10px; 
        text-align: center; 
        font-family: monospace; 
        font-size: 16px; 
        line-height: 1.6; 
        width: 90%;
        max-width: 800px;
        margin: auto;
        box-shadow: 2px 2px 10px rgba(0, 0, 0, 0.1);
      }
    "))
  ),

  # Title Section with Fixed DNA Animation
  div(class = "title-box",
    style = "position: relative; 
             font-size: 28px; 
             font-weight: bold; 
             color: #2E86C1; 
             background-color: #f0f0f0; 
             padding: 15px; 
             border-radius: 10px; 
             border: 2px solid #d3d3d3;
             box-shadow: 2px 2px 10px rgba(0, 0, 0, 0.2);
             text-align: center;",  # Keep text centered

    # DNA Animation positioned in the top-left (does not move)
    div(style = "position: absolute; 
                 top: 10px; 
                 left: 10px; 
                 width: 60px; height: auto;",
        img(src = "https://raw.githubusercontent.com/quangsang52sh/Programs/main/data/original-770ce4c43c5396a4e2c979eb302f059a.gif", 
            style = "width: 100%; height: auto;")
    ),
    
    # Keep Title Text Centered
    HTML("<span style='background-color: green; color: white; font-weight: bold;'>Phylogenetic Analysis (Version 5.0)</span>")
),

  sidebarLayout(
    # Sidebar
    sidebarPanel(
      width = 4,
      fileInput("file", "Upload File (Fasta alignment sequence)", accept = c(".fas", ".fasta", ".fa", ".txt", ".tsv")),
      checkboxGroupInput("tree_methods", "Select Phylogenetic Methods:", 
                         choices = list("Neighbor Joining (NJ)" = "NJ",
                                        "Maximum Parsimony (MP)" = "MP",
                                        "Maximum Likelihood (ML)" = "ML",
                                        "Bayesian Inference (BI - MrBayes)" = "Mrbayes")),
      checkboxInput("run_ML_withoutModel", "Run ML Tree (without model test)"),
      textInput("modelname", "Model Name (for ML without modeltest)", value = "GTR+G(4)+I"),
      sliderInput("bootstrapSlider", "Bootstrap Value", min = 0, max = 1000, value = 100),
      actionButton("checkModelButton", "Check Best Model"),
      verbatimTextOutput("bestModelOutput"),
      br(), br(),
      actionButton("runButton", "Run Analysis"),
      actionButton("stopButton", "Stop Analysis",
             style = "background-color: darkred; color: white; font-weight: bold;"),
      br(), br(),
      actionButton("quitButton", "Quit App", 
                   style = "background-color: red; color: white; font-weight: bold;")
    ),

    # Main Panel
    mainPanel(
      width = 8,
      fluidRow(
        column(width = 12,
          div(class = "profile-pic",
              img(src = "https://raw.githubusercontent.com/quangsang52sh/Programs/main/data/51614765.png", 
                  alt = "Profile Picture")
          )
        ),
        column(width = 8, 
          verbatimTextOutput("output"), 
          uiOutput("greeting"), 
          verbatimTextOutput("progressText")
        ),
        column(width = 12,
          div(class = "welcome-box",
    style = "background-color: #f5f5f5; 
             padding: 15px; 
             border: 2px solid #d3d3d3; 
             border-radius: 10px; 
             text-align: left; 
             font-family: monospace; 
             font-size: 17px;
             font-weight: bold;
             white-space: pre-wrap; 
             max-width: 600px;
             word-wrap: break-word;
             overflow: hidden;
             margin: auto;
             box-shadow: 2px 2px 10px rgba(0, 0, 0, 0.1);",
    
    HTML("<pre style='font-family: monospace; font-size: 17px; text-align: left; line-height: 1.4;'>
##############################################
Welcome to Tsang Script!                     
Sang Tran Quang  <span style='background-color: blue; color: white; font-weight: bold;'>(LIMITED VERSION)</span>           
Email: <a href='mailto:sangtq@ntu.edu.vn'>sangtq@ntu.edu.vn</a>                   
Rscript 5.0, Automation Analysis by TSANG              
##############################################
    </pre>")
          )
        )
      )
    )
  )
)



server <- function(input, output, session) {
  analysisDone <- reactiveVal(FALSE)
  bestModel <- reactiveVal(NULL)

  # Model Checking
  observeEvent(input$checkModelButton, {
    if (!is.null(input$file)) {
      myfile <- read.dna(input$file$datapath, format = "fasta")
      dat <- as.phyDat(myfile)

      # Run model test
      model_results <- modelTest(dat)
      write.csv(model_results, file = "MLmodel_checking.csv")
      # Get best model based on minimum AIC
      best_model <- model_results[which.min(model_results$AIC), "Model"]
      bestModel(best_model)  # Store best model
    }
  })

  output$bestModelOutput <- renderText({
    if (!is.null(bestModel())) {
      paste("Best Model Selected:", bestModel())
    } else {
      "Click 'Check Best Model' to determine the best model."
    }
  })

 # Running Analysis with Progress Bar
observeEvent(input$runButton, {
  if (!is.null(input$file)) {
    withProgress(message = "Running Phylogenetic Analysis", value = 0, {
      
      # Read input file
      myfile <- read.dna(input$file$datapath, format = "fasta")
      dat <- as.phyDat(myfile)
      
      incProgress(0.1, detail = "File loaded successfully...")

      # Run NJ safely
      if (isTRUE("NJ" %in% input$tree_methods)) {
        run_NJ(dat, input$bootstrapSlider)
        incProgress(0.2, detail = "Neighbor Joining (NJ) analysis completed...")
      }

      # Run MP safely
      if (isTRUE("MP" %in% input$tree_methods)) {
        run_MP(dat, input$bootstrapSlider)
        incProgress(0.2, detail = "Maximum Parsimony (MP) analysis completed...")
      }

      # Run ML safely
      if (isTRUE("ML" %in% input$tree_methods)) {
        run_ML(dat, input$bootstrapSlider)
        incProgress(0.2, detail = "Maximum Likelihood (ML) analysis completed...")
      }

      # Run ML without model test safely
      if (isTRUE(input$run_ML_withoutModel)) {
        run_ML_withoutModel(dat, input$modelname, input$bootstrapSlider)
        incProgress(0.2, detail = "ML analysis without model test completed...")
      }

      # Run MrBayes safely
      if (isTRUE("Mrbayes" %in% input$tree_methods)) {
        incProgress(0.5, detail = "Running Bayesian Inference (MrBayes)... This may take a while.")
        run_MrBayes(dat)
        incProgress(1, detail = "Bayesian analysis completed!")
      }

      analysisDone(TRUE)
    })
  }
})

  output$progressText <- renderText({
    if (is.null(input$file) || !analysisDone()) {
      "Waiting for input file..."
    } else {
      "Analysis in progress..."
    }
  })

  output$output <- renderText({
    if (analysisDone()) {
      "Analysis completed!"
    }
  })

  output$greeting <- renderUI({
    if (!is.null(input$file) && analysisDone()) {
      HTML("<pre>Have a nice day! ....</pre><pre>Good luck to you ....</pre>")
    }
  })

  observeEvent(input$quitButton, {
    stopApp()  # Stops the Shiny application
  })
}


phyloApp3 <- function() {
  shinyApp(ui = ui, server = server) }

phyloApp3()
