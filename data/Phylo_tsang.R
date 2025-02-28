# Load required packages
if (!require("ape")) {
  install.packages("ape")
}
library(ape)
if (!require("phangorn")) {
  install.packages("phangorn")
}
library(phangorn)
#if (!require("tidyverse")) {
#  install.packages("tidyverse")
#}
#library(tidyverse)
if (!require("zip")) {
  install.packages("zip")
}
library(zip)
  
ID_name="tsang"
#print
options(max.print=1000000)
# Information input
cat("

##################################################################################################
##################################################################################################
##                               Phylogenetic tree (version 3.2)                                ##
##                               Welcome to tsang script !.....                                 ##
##                               This is modified by tsang                                      ##
##                               Email: sangtq@ntu.edu.vn                                       ##
##################################################################################################
##################################################################################################
")
cat("\n________________________Greeting to you", ID_name," ....  ________________________________________\n\n")
cat("
Please! Prior to running this script, following this tutorior
1) open Rscript and type 'source(tsang_script.R)', if already source, try to next step
2) runPhylogeneticAnalysis(file_name, ML = 'FALSE', NJ = 'FALSE', MP = 'FALSE', Mrbayes = 'FALSE', RaxML = 'FALSE', modeltestNG = 'FALSE', without_modeltest='FALSE', Model='GTR+G+I')
if run ML, turn into TRUE 
if run NL, turn into TRUE
if run MP, turn into TRUE
if you are already run model test and has a name like 'GTR+G+I', please use (without_modeltest='TRUE', Model='GTR+I'), if not type modelname, this is automatic run un default GTR+G+I
if run modeltest only, try to (modeltestNG = TRUE)
if run RaxML with modeltest, try to run modeltestNG first ( modeltestNG = 'TRUE') and then run RaxML (RaxML = 'TRUE', Model = 'modeltestNG[name]')
Thank you!
Tsang

")

run_ML <- function(dat) {
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
  bs <- bootstrap.pml(fit1, bs = 1000, optNni = TRUE, control = pml.control(trace = 2))
  tree <- plotBS(midpoint(fit1$tree), bs, p = 50, type = "p")
  write.tree(tree, "MLtree.tre")
  save.image("runML.Rdata")
}

run_ML_withoutModel <- function(dat, modelname) {

  best_model <- modelname
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
  bs <- bootstrap.pml(fit1, bs = 1000, optNni = TRUE, control = pml.control(trace = 2))
  tree <- plotBS(midpoint(fit1$tree), bs, p = 50, type = "p")
  write.tree(tree, "MLtree_withoutModel.tre")
  save.image("runML_withoutModel.Rdata")
}

run_modeltestNG <- function(myfile) {

# Check the operating system
if (Sys.info()["sysname"] == "Linux") {
  #system("rm -r mrbayes-3.2.7/")
  download.file(url = "https://raw.githubusercontent.com/quangsang52sh/Programs/main/data/modeltest-ng-0.1.7-linux.zip", destfile = "modeltest-ng-0.1.7-linux.zip")
  unzip("modeltest-ng-0.1.7-linux.zip")
  mt_path <- paste0(getwd(),"/modeltest-ng-0.1.7-static/modeltest-ng-static")
  if (dir.exists("mtng_results") == "TRUE") {
  mt_output = paste0(getwd(),"/mtng_results")
  } else {
  dir.create("mtng_results")
  mt_output = paste0(getwd(),"/mtng_results")
  }
} else if (Sys.info()["sysname"] == "Darwin") {
  download.file(url="https://raw.githubusercontent.com/quangsang52sh/Programs/main/data/modeltest-ng-MacOS.zip",destfile="modeltest-ng-MacOS.zip")
  unzip("modeltest-ng-MacOS.zip")
  mt_path <- paste0(getwd(),"/modeltest-ng/modeltest-ng-osx")
    if (dir.exists("mtng_results") == "TRUE") {
  mt_output = paste0(getwd(),"/mtng_results")
  } else {
  dir.create("mtng_results")
  mt_output = paste0(getwd(),"/mtng_results")
  }
} else {
  stop("Unsupported operating system.")
}

system(paste0(mt_path, " -i ", myfile, " -o ", mt_output,"/modeltest"," -v"))
if (file.create(paste0(mt_output,"/modeltest.out")) == "TRUE") {
 print(paste0("Your detail model have been saved in ",mt_output,"/modeltest.out directory"))
 print(paste0("Have nice a day!....."))
 print(paste0("Tsang."))
}
else {
 stop("You are running fail modeltestng. Please try again !")
}
}


run_RaxML <- function(myfile,modelname) {

if (Sys.info()["sysname"] == "Linux") {
  dir.create("raxml-ng_v1.2.0_linux")
  download.file(url = "https://raw.githubusercontent.com/quangsang52sh/Programs/main/data/raxml-ng_v1.2.0_linux_x86_64.zip", destfile = "./raxml-ng_v1.2.0_linux/raxml-ng_v1.2.0_linux_x86_64.zip")
  zip::unzip("./raxml-ng_v1.2.0_linux/raxml-ng_v1.2.0_linux_x86_64.zip", exdir = "./raxml-ng_v1.2.0_linux")
  raxml_path <- paste0(getwd(),"/raxml-ng_v1.2.0_linux/raxml-ng")
  if (dir.exists("raxml_results") == "TRUE") {
  raxml_output = paste0(getwd(),"/raxml_results")
  } else {
  dir.create("raxml_results")
  raxml_output = paste0(getwd(),"/raxml_results")
  }
} else if (Sys.info()["sysname"] == "Darwin") {
  dir.create("raxml-ng_v1.2.0_macos")
  download.file(url="https://raw.githubusercontent.com/quangsang52sh/Programs/main/data/raxml-ng_v1.2.0_macos_x86_64.zip",destfile="./raxml-ng_v1.2.0_macos/raxml-ng_v1.2.0_macos_x86_64.zip")
  zip::unzip("./raxml-ng_v1.2.0_macos/raxml-ng_v1.2.0_macos_x86_64.zip", exdir = "./raxml-ng_v1.2.0_macos")
  raxml_path <- paste0(getwd(),"/raxml-ng_v1.2.0_macos/raxml-ng")
    if (dir.exists("raxml_results") == "TRUE") {
  raxml_output = paste0(getwd(),"/raxml_results")
  } else {
  dir.create("raxml_results")
  raxml_output = paste0(getwd(),"/raxml_results")
  }
} else {
  stop("Unsupported operating system.")
}

system(paste0(raxml_path," --all --msa ",myfile," --model ", modelname," --bs-trees 1000"))

source_folder <- paste0(getwd())
pattern <- myfile
destination_folder <- raxml_output
files_to_move <- list.files(path = source_folder, pattern = pattern, full.names = TRUE)
file.rename(files_to_move, file.path(destination_folder, basename(files_to_move)))

RaxMLTree <- read.tree(paste0(raxml_output,"/",pattern,".raxml.bestTree"))
bs_raxml <- read.tree(paste0(raxml_output,"/",pattern,".raxml.bootstraps"))
treeRaxML <- plotBS(midpoint(RaxMLTree), bs_raxml, p = 50, type = "p")
write.tree(treeRaxML, "treeRaxML.tre")
save.image("treeRaxML.Rdata")

}


run_MrBayes <- function(dat) {

# Check the operating system
if (Sys.info()["sysname"] == "Windows") {
  download.file(url = "https://github.com/NBISweden/MrBayes/releases/download/v3.2.7/MrBayes-3.2.7-WIN.zip", destfile = "MrBayes-3.2.7-WIN.zip")
  zip::unzip("MrBayes-3.2.7-WIN.zip")
  mrbayes_path <- paste0(getwd(),"/MrBayes-3.2.7-WIN/bin/mb.3.2.7-win64.exe")
  if (dir.exists("Mrbayes_results") == "TRUE") {
  Mrbayes_output = paste0(getwd(),"/Mrbayes_results")
  }
  else {
  dir.create("Mrbayes_results")
  Mrbayes_output = paste0(getwd(),"/Mrbayes_results")
  }
} else if (Sys.info()["sysname"] == "Linux") {
  #system("rm -r mrbayes-3.2.7/")
  if (!file.exists("./MrBayes/src/mb")) {
  system("git clone --depth=1 https://github.com/NBISweden/MrBayes.git")
  system("cd ./MrBayes && ./configure")
  system("cd ./MrBayes && make")
  mrbayes_path <- "./MrBayes/src/mb" 
  } else { mrbayes_path <- "./MrBayes/src/mb" }
  if (dir.exists("Mrbayes_results") == "TRUE") {
  Mrbayes_output = paste0(getwd(),"/Mrbayes_results")
  }
  else {
  dir.create("Mrbayes_results")
  Mrbayes_output = paste0(getwd(),"/Mrbayes_results")
  }
} else if (Sys.info()["sysname"] == "Darwin") {
  download.file(url="https://github.com/NBISweden/MrBayes/releases/download/v3.2.6/MrBayes-3.2.6_MACx64.zip",destfile="MrBayes-3.2.6_MACx64.zip")
  zip::unzip("MrBayes-3.2.6_MACx64.zip")
  mrbayes_path <- paste0(getwd(),"/MrBayes/mb")
    if (dir.exists("Mrbayes_results") == "TRUE") {
  Mrbayes_output = paste0(getwd(),"/Mrbayes_results")
  }
  else {
  dir.create("Mrbayes_results")
  Mrbayes_output = paste0(getwd(),"/Mrbayes_results")
  }
} else {
  stop("Unsupported operating system.")
}

bin <- as.DNAbin(dat)
write.nexus.data(bin,"mrbayes_input.nexus")
input = paste0("mrbayes_input.nexus")
if (file.exists("autoRunMrbayes.txt") == "TRUE") {
  auto = paste0("./autoRunMrbayes.txt") }
else {
  download.file("https://raw.githubusercontent.com/quangsang52sh/Programs/main/data/autoRunMrbayes.txt","autoRunMrbayes.txt")
  auto = paste0("./autoRunMrbayes.txt")
  }

system(paste0(mrbayes_path, " < ", auto, " > mrbayes_input_log.txt"))
source_folder <- paste0(getwd())
pattern <- "mrbayes_input"
destination_folder <- Mrbayes_output
files_to_move <- list.files(path = source_folder, pattern = pattern, full.names = TRUE)
file.rename(files_to_move, file.path(destination_folder, basename(files_to_move)))
MrbayesTree <- read.nexus(paste0(Mrbayes_output,"/",input,".con.tre"))
run1 <- read.nexus(paste0(Mrbayes_output,"/",input,".run1.t"))
run2 <- read.nexus(paste0(Mrbayes_output,"/",input,".run2.t"))
# burnning 10%
bs_Mrbayes <- c(run1[101:10001],run2[101:10001])
MrbayesTree <- plotBS(midpoint(MrbayesTree), bs_Mrbayes, p = 50, type = "p")
write.tree(MrbayesTree, "MrbayesTree.tre")
save.image("runMrbayes.Rdata")

}


run_MP <- function(dat) {
# Maximum Parsimony (MP) tree
MegaMPTree <- pratchet(dat)
MegaMPTree <- acctran(MegaMPTree, dat)
bs_megaMP <- bootstrap.phyDat(dat, pratchet, bs = 1000)
MPTree <- plotBS(midpoint(MegaMPTree), bs_megaMP, p = 50, type = "p")
add.scale.bar()
write.tree(MPTree, "MPtree.tre")
save.image("runMP.Rdata")
}


run_NJ <- function(dat) {

# Neighbor Joining (NJ) tree
#myfile <- read.dna(file_name, format = "fasta")
#dat <- as.phyDat(myfile)
tree <- nj(dist.ml(dat))
MegaNJTree <- tree
bs_megaNJ <- bootstrap.phyDat(dat, FUN = function(x) NJ(dist.hamming(x)), bs = 1000)
treeNJ <- plotBS(midpoint(MegaNJTree), bs_megaNJ, p = 50, type = "p")
write.tree(treeNJ, "NJtree.tre")
save.image("runNJ.Rdata")
}



runPhylogeneticAnalysis <- function(file_name, ML = FALSE, NJ = FALSE, MP = FALSE, Mrbayes = FALSE, RaxML = FALSE, without_modeltest = FALSE, modeltestNG = FALSE, Model = "GTR+G+I") {

if (ML == "TRUE") {
    myfile <- read.dna(file_name, format = "fasta")
    dat <- as.phyDat(myfile)
    run_ML(dat)

   }

  if (NJ == "TRUE") {
    myfile <- read.dna(file_name, format = "fasta")
    dat <- as.phyDat(myfile)
    run_NJ(dat)
  } 

  if (MP == "TRUE") {
    myfile <- read.dna(file_name, format = "fasta")
    dat <- as.phyDat(myfile)
    run_MP(dat)
  }
  
  if (Mrbayes == "TRUE") {
    myfile <- read.dna(file_name, format = "fasta")
    dat <- as.phyDat(myfile)
    run_MrBayes(dat)
  }
  
    if ((ML == "TRUE") & (MP == "TRUE")) {
    myfile <- read.dna(file_name, format = "fasta")
    dat <- as.phyDat(myfile)
    run_ML(dat)
    run_MP(dat)
   
  } 
  
    if ((ML == "TRUE") & (NJ == "TRUE")) {
    myfile <- read.dna(file_name, format = "fasta")
    dat <- as.phyDat(myfile)
    run_ML(dat)
    run_NJ(dat)
   
  } 
  
    if ((NJ == "TRUE") & (MP == "TRUE")) {
    myfile <- read.dna(file_name, format = "fasta")
    dat <- as.phyDat(myfile)
    run_NJ(dat)
    run_MP(dat)
   
  } 

    if ((ML == "TRUE") & (NJ == "TRUE") & (MP == "TRUE")) {
    myfile <- read.dna(file_name, format = "fasta")
    dat <- as.phyDat(myfile)
    run_ML(dat)
    run_NJ(dat)
    run_MP(dat)
   
  }

    if ((ML == "FASLE") & (NJ == "FALSE") & (MP == "FALSE") & (without_modeltest == "FALSE")) {
    myfile <- read.dna(file_name, format = "fasta")
    dat <- as.phyDat(myfile)
    run_NJ(dat)
   
  }

    if (without_modeltest == "TRUE") {
    myfile <- read.dna(file_name, format = "fasta")
    dat <- as.phyDat(myfile)
    run_ML_withoutModel(dat, Model)
   
  }

    if (modeltestNG == "TRUE") {
    run_modeltestNG(file_name)
   
  }

    if (RaxML == "TRUE") {
    run_RaxML(file_name, Model)
   
  }

}
