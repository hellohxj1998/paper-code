library("glmnet")
library("survival")
library("pbapply")
library("survivalROC")
library("e1071")
library("VennDiagram")
library("grid")

if(file.exists("nonzeroCoef.R")) source("nonzeroCoef.R")
if(file.exists("msvmRFE.R")) source("msvmRFE.R")

Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)

myexpr <- read.csv("easy_input_exp.csv", header = T, row.names = 1)
mysurv <- read.csv("easy_input_suv.csv", header = T, row.names = 1)

set.seed(1314)
cvfit <- cv.glmnet(t(myexpr), 
                   Surv(mysurv$months, mysurv$status), 
                   nfold = 10,
                   family = "cox")

plot(cvfit)

fit <- glmnet(t(myexpr), 
              Surv(mysurv$months, mysurv$status), 
              family = "cox")

plot(fit, label = TRUE)

coef.min <- coef(cvfit, s = "lambda.min")
active.min <- which(coef.min != 0)
lasso_gene_ids <- rownames(myexpr)[active.min]
lasso_coefs <- coef.min[active.min]

lasso_result <- cbind(Gene = lasso_gene_ids, Coefficient = lasso_coefs)
write.csv(lasso_result, "LASSO_selected_genes.csv", row.names = FALSE)

input <- input_data 

nfold <- 5
nrows <- nrow(input)
folds <- rep(1:nfold, len = nrows)[sample(nrows)]
folds <- lapply(1:nfold, function(x) which(folds == x))

results <- lapply(folds, svmRFE.wrap, input, k = 5, halve.above = 100)

top.features <- WriteFeatures(results, input, save = FALSE)
write.csv(top.features, "SVM_feature_ranking.csv")

featsweep <- lapply(1:5, FeatSweep.wrap, results, input) 

no.info <- min(prop.table(table(input[, 1])))
errors <- sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))