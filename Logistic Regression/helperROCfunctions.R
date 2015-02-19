#####
# Helper file with functions that construct ROC curves for glm and penfit objects
#####

# Build a function that makes a ROC curve for a given model (GLM or penfit)

constructROC <- function(model, ...) {
  UseMethod("constructROC", model)
}

constructROC.logistf <- function(model, test = TRUE, thr = seq(from = 0, to = 1, by = 0.1)) {
  m <- length(thr)
  roc <- data.frame(x = numeric(m), y = numeric(m))
  
  for(i in seq_along(thr)) {
    # Compute predictions
    # if test = TRUE, then use the test set asthma_test, otherwise use the training set.
    if (test) { 
      preds <- as.matrix(data.frame(b = 1, asthma_test[, c(86, 256, 43)])) %*% model$coefficients
      preds <- ifelse(exp(preds)/(exp(preds) + 1) >= thr[i], "case", "control")
      
      obs <- asthma_test$CASE
    } else {
      preds <- ifelse(model$predict >= thr[i], "case", "control")
      obs <- relevel(factor(model$data$CASE, levels = c(1,0), labels = c("case", "control")), ref = "control")
    }
    
    preds <- factor(preds, levels = levels(asthma$CASE))
    
    # Compute confusion matrix
    confMatrix <- table(obs, preds, useNA = "always")    
    # compute false positive rate.
    roc[i, "x"] = confMatrix[obs = "control", preds = "case"] / sum(confMatrix[obs = "control", ])
    # compute true positive rate.
    roc[i, "y"] = confMatrix[obs = "case", preds = "case"] / sum(confMatrix[obs = "case", ])
  }  
  
  return(dplyr::arrange(roc, x, y)) 
}

constructROC.penfit <- function(model, test = TRUE, thr = seq(from = 0, to = 1, by = 0.1)) {
  m <- length(thr)
  roc <- data.frame(x = numeric(m), y = numeric(m))
  
  for(i in seq_along(thr)) {
    # Compute predictions
    # if test = TRUE, then use the test set asthma_test, otherwise use the training set.
    if (test) { 
      preds <- ifelse(predict(object = model,
                              penalized = asthma_test[, -1],
                              data = asthma_test) >= thr[i], 
                      "case", 
                      "control")
      obs <- asthma_test$CASE
    } else {
      preds <- ifelse(predict(object = model,
                              penalized = asthma_training[, -1],
                              data = asthma_training) >= thr[i], 
                      "case", 
                      "control")
      obs <- asthma_training$CASE
    }
    
    preds <- factor(preds, levels = levels(asthma$CASE))
    
    # Compute confusion matrix
    confMatrix <- table(obs, preds, useNA = "always")    
    # compute false positive rate.
    roc[i, "x"] = confMatrix[obs = "control", preds = "case"] / sum(confMatrix[obs = "control", ])
    # compute true positive rate.
    roc[i, "y"] = confMatrix[obs = "case", preds = "case"] / sum(confMatrix[obs = "case", ])
  }
  
  return(dplyr::arrange(roc, x, y)) 
}

# Build a function that returns a ROC plot
plotROC <- function(roc) {
  p <- ggplot(roc, aes(x = x, y = y)) +
    geom_line(col = "blue") + 
    geom_point(col = "blue") + 
    geom_line(data = data.frame(x = c(0, 1), y = c(0, 1)), aes(x = x, y = y), linetype = 2)
  
  return(p)
}

# Build a function to assess how well a model fits the given data
fitData <- function(model, ...) {
  UseMethod("fitData", model)
}

fitData.penfit <- function(model) {
  # Get the predictions and the observed values
  preds <- ifelse(predict(model, penalized = asthma_training[, -1], data = asthma_training) >= 0.5, "case", "control")
  preds <- factor(preds, levels = levels(asthma_training$CASE))
  obs <- asthma_training$CASE
  
  # Compute the confusion matrix
  confMatrix <- caret::confusionMatrix(data = preds, 
                                       reference = obs,
                                       positive = "case")
  
  # Return interesting values
  list(table = confMatrix$table,
       stats = c(confMatrix$overall[1],
                 confMatrix$byClass[1:4]),
       roc = constructROC(model, test = FALSE))
}

fitData.logistf <- function(model) {
  # Get the predictions and the observed values
  preds <- ifelse(model$predict >= 0.5, "case", "control")
  preds <- factor(preds, levels = levels(asthma$CASE))
  obs <- relevel(factor(model$data$CASE, levels = c(1,0), labels = c("case", "control")), ref = "control")
  
  # Compute the confusion matrix
  confMatrix <- caret::confusionMatrix(data = preds, 
                                       reference = obs,
                                       positive = "case")
  
  # Return interesting values
  list(table = confMatrix$table,
       stats = c(confMatrix$overall[1],
                 confMatrix$byClass[1:4]),
       roc = constructROC(model, test = FALSE))
}