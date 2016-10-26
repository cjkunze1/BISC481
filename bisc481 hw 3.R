#############################
# Colin Kunze
# ckunze@usc.edu
# 10.25.2016
# Homework 3
# BISC 481
#############################


#### If the following packages are not installed, uncomment the single '#' lines and install ####

## Bioconductor
# source("https://bioconductor.org/biocLite.R")
# biocLite()

## DNAshapeR
# biocLite("DNAshapeR")

## Caret
#install.packages("caret")
#install.packages("e1071")
#install.packages("ROCR")

## Install and initialize packages
#install.packages("ggplot2")
#install.packages("grid")

#### Initialization ####
library(DNAshapeR)
library(caret)
library(ggplot2)
library(grid)
library(ROCR)
library(Biostrings)
#------- Adjust the working path to YOUR working path-------
workingPath <- "/Users/Colin/Documents/SCHOOLWORK/Senior yearrr/BISC481 structural bio/BISC481-master/gcPBM/"
#-----------------------------------------------------------

# All of the different DNA sequences we want to analyze
sequences <- c("Mad.txt", "Max.txt", "Myc.txt")

#vectors for the results
results <- vector(,3)
names(results) <- sequences

results1mer <- vector(,3)
names(results1mer) <- sequences

for (i in sequences){
  ## Acquire the DNA shape prediction
  fn_fasta <- paste0(workingPath, i, ".fa") #full file name of the DNA sequence being analyzed (including directories)
  pred <- getShape(fn_fasta)         #DNA shape prediction from sequence
  
  ## Encode feature vectors for 1-mer only
  featureVector1mer <- encodeSeqShape(fn_fasta, pred, "1-mer")
  
  ## Encode feature vectors for 1-mer and 1-shape
  featureType <- c("1-mer", "1-shape")
  featureVector <- encodeSeqShape(fn_fasta, pred, featureType)
  
  ## Build MLR model by using Caret
  # Data preparation
  fn_exp <- paste0(workingPath, i) # name of the dna sequence and information in.txt format
  exp_data <- read.table(fn_exp)   # reads the data into an object
  df <- data.frame(affinity=exp_data$V2, featureVector)  #turns the data into a data.frame object// "1-mer" "1-shape" features
  df1mer <- data.frame(affinity=exp_data$V2, featureVector1mer)  #does the same for the "1-mer" features
  
  # Arguments setting for Caret
  trainControl <- trainControl(method = "cv", number = 10, savePredictions = TRUE) #uses 10-fold cross validation in resampling
  
  
  # Prediction with L2-regularized (helps prevent overfitting)
  model2 <- train(affinity~., data = df, trControl=trainControl, 
                  method = "glmnet", tuneGrid = data.frame(alpha = 0, lambda = c(2^c(-15:15)))) #using both 1mer and shape
  
  model21mer <- train(affinity~., data = df1mer, trControl=trainControl,
                      method = "glmnet", tuneGrid = data.frame(alpha = 0, lambda = c(2^c(-15:15)))) #using only 1mer
  
  results[i] <- model2$results$Rsquared[1]
  results1mer[i] <- model21mer$results$Rsquared[1]
  
  print(paste0("R-squared value of using '1-mer' model for ", i))
  print(unname(results1mer[i]))
  
  print(paste0("R-squared value using '1mer' +'1-shape' model for ", i))
  print(unname(results[i]))
}

plottheme <- theme(
  plot.margin = unit(c(0.1, 0.5, 0.1, 0.1), "cm"),
  axis.text = element_text(colour="black", size=12),
  axis.title.x = element_text(colour="black", size=12),
  axis.title.y = element_text(colour="black", size=12),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  axis.line = element_line(colour = "black"),
  axis.text = element_text(colour ="black"),
  axis.ticks = element_line(colour = "black")
)

#plotting the r-squared values of the model with shape against the old model
ggplot() +
  geom_point(aes(x = results1mer, y = results), color = "red", size=1) +
  geom_abline(slope=1) + geom_vline(xintercept=0) + geom_hline(yintercept=0) +
  coord_fixed(ratio = 1, xlim = c(0,1), ylim = c(0,1)) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  plottheme

#################
# Determining the differences in physical features between bound and unbound DNA
##################

plotseq <- c("bound_500.fa", "unbound_500.fa") #the sequences to be analyzed
newworkingpath <- "/Users/Colin/Documents/SCHOOLWORK/Senior yearrr/BISC481 structural bio/BISC481-master/CTCF/"

for (i in plotseq){
  
  #getting the predicted shape of the sequences
  filename <- paste0(newworkingpath, i)
  prediction <- getShape(filename)
  
  #plotting the minior groove width and its heat map
  dev.off()
  plotShape(prediction$MGW)
  dev.print(pdf, paste0(newworkingpath, i, "_MGW_plot.pdf"))
  heatShape(prediction$MGW, 20)
  dev.print(pdf, paste0(newworkingpath, i, "_MGW_HeatMap.pdf"))
  
  #plotting the propeller twist and its heat map
  dev.off()
  plotShape(prediction$ProT)
  dev.print(pdf, paste0(newworkingpath, i, "_ProT_plot.pdf"))
  heatShape(prediction$ProT, 20)
  dev.print(pdf, paste0(newworkingpath, i, "_ProT_HeatMap.pdf"))
  
  #plotting the roll and its heat map
  dev.off()
  plotShape(prediction$Roll)
  dev.print(pdf, paste0(newworkingpath, i, "_Roll_plot.pdf"))
  dummycolumn <- matrix(0, ncol = 1, nrow = 5000)
  rollmatrix <- cbind(prediction$Roll, dummycolumn)
  heatShape(rollmatrix, 20)
  dev.print(pdf, paste0(newworkingpath, i, "_Roll_HeatMap.pdf"))
  
  #plotting the helical twist and its heat map
  dev.off()
  plotShape(prediction$HelT)
  dev.print(pdf, paste0(newworkingpath, i, "_HelT_plot.pdf"))
  heltmatrix <- cbind(prediction$HelT, dummycolumn)
  heatShape(heltmatrix, 20)
  dev.print(pdf, paste0(newworkingpath, i, "_HelT_HeatMap.pdf"))
}

#############################################
# Applying logistic regression
#############################################

## Generate data for the classifcation (assign Y to bound and N to non-bound)
# bound
boundFasta <- readDNAStringSet(paste0(newworkingpath, "bound_500.fa"))
sequences <- paste(boundFasta)
boundTxt <- data.frame(seq=sequences, isBound="Y")

# non-bound
nonboundFasta <- readDNAStringSet(paste0(newworkingpath, "unbound_500.fa"))
sequences <- paste(nonboundFasta)
nonboundTxt <- data.frame(seq=sequences, isBound="N")

# merge two datasets
writeXStringSet( c(boundFasta, nonboundFasta), paste0(newworkingpath, "ctcf.fa"))
exp_data <- rbind(boundTxt, nonboundTxt)


## DNAshapeR prediction
pred <- getShape(paste0(newworkingpath, "ctcf.fa"))

##Encode feature vector for just using 1-mer
featureVector1mer <- encodeSeqShape(paste0(newworkingpath, "ctcf.fa"), pred, "1-mer")
df1mer <- data.frame(isBound = exp_data$isBound, featureVector1mer)


## Encode feature vectors for using shape and 1-mer
featureType <- c("1-mer", "1-shape")
featureVector <- encodeSeqShape(paste0(newworkingpath, "ctcf.fa"), pred, featureType)
df <- data.frame(isBound = exp_data$isBound, featureVector)


## Logistic regression
# Set parameters for Caret
trainControl <- trainControl(method = "cv", number = 10, 
                             savePredictions = TRUE, classProbs = TRUE)
# Perform prediction
model <- train(isBound~ ., data = df, trControl = trainControl,
               method = "glm", family = binomial, metric ="ROC")
model1mer <- train(isBound~ ., data = df1mer, trControl = trainControl,
                   method = "glm", family = binomial, metric = "ROC")

## Plot AUROC
prediction <- prediction( model$pred$Y, model$pred$obs )
prediction1mer <- prediction(model1mer$pred$Y, model1mer$pred$obs )
performance <- performance( prediction, "tpr", "fpr" )
preformance1mer <- performance( prediction1mer, "tpr", "fpr")
plot(performance)
par(new=TRUE)
plot(preformance1mer, col="green")

## Caluculate AUROC
auc <- performance(prediction, "auc")
auc1mer <- performance(prediction1mer, "auc")
auc <- unlist(slot(auc, "y.values"))
auc1mer <- unlist(slot(auc1mer, "y.values"))

print("auc for 1mer model")
print(auc1mer)

print("auc for 1-shape and 1-mer model")
print(auc)


