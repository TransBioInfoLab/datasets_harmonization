##################################################################################################
###                                                                                              #      
###        This Script is used to do age prediction based on DNA methylation data(450K)          #
###        Coefficients of the predictor are based on 13,566 training samples                    #
###        Two age predictors based on different methods(Elastic Net and BLUP) can be used       #
###        Qian Zhang 27-03-2018, Email: q.zhang@uq.edu.au                                       #
##################################################################################################




############# for each probe, change to missing value to the mean value across all individuals #############
addna<-function(methy){
  methy[is.na(methy)]<-mean(methy,na.rm=T)
  return(methy)
}

replace_missing <- function(data) {
  print("1.1 Replacing missing values with mean value")
  dataNona<-apply(data,2,function(x) addna(x))   ###############  replace the NA with mean value for each probe 
  
  dataNona<- dataNona[,apply(dataNona,2,function(x) sum(is.na(x)))!=nrow(dataNona)] ############ remove the probe when it has NA across all individuals
  print(paste0(ncol(data) - ncol(dataNona)," probe(s) is(are) removed since it has (they have) NA across all individuals"))
  
  
  print("1.2 Standardizing")
  dataNona.norm<- apply(dataNona,1,scale)        ############### standardize the DNA methylation within each individual, remove the mean and divided by the SD of each individual     Probe * IND
  rownames(dataNona.norm)<-colnames(dataNona)
  
  dataNona.norm
}

load_predictor <- function(refdir, predname = c("en", "blup")) {
  predname <- match.arg(predname)
  
  coef <- read.table(
    file.path(ref_dir, paste0(predname, ".coef")),
    stringsAsFactor = FALSE,
    header = TRUE
  )
  
  pred_int <- coef[1,2]
  pred_coef <- coef[-1,]
  
  row.names(pred_coef) <- pred_coef$probe
  
  list(
    pred_coef = pred_coef, pred_int = pred_int
  )
}

check_missing <- function(dataNona.norm, pred_coef) {
  pred_comm <- intersect(row.names(pred_coef), row.names(dataNona.norm))
  pred_diff <- nrow(pred_coef) - length(pred_comm)
  
  print(
    paste0(pred_diff, " probe(s) in predictor is(are) not in the data")
  )
  
  pred_comm
}

run_age_prediction <- function(
    dataNona.norm, refdir, predname = c("en", "blup")
) {
  predname <- match.arg(predname)
  
  if (predname == "en") {
    message("Performing Elastic Net age prediction")
  } else {
    message("Performing BLUP age prediction")
  }
  
  predictor_ls <- load_predictor(refdir, predname = predname)
  pred_coef <- predictor_ls$pred_coef
  pred_int <- predictor_ls$pred_int
  pred_comm <- check_missing(dataNona.norm, pred_coef)
  
  pred_coef <- pred_coef[pred_comm,]
  pred_res <- pred_coef$coef %*% dataNona.norm[pred_comm,] + pred_int
  
  pred_res
}

predict_age <- function(data, refdir) {
  # adjust data
  if(nrow(data) > ncol(data)){
    print("I guess you are using Probe in the row, data will be transformed!!!")
    data <- t(data)
  }
  
  dataNona.norm <- replace_missing(data)
  
  pred_en <- run_age_prediction(dataNona.norm, refdir, predname = "en")
  pred_blup <- run_age_prediction(dataNona.norm, refdir, predname = "blup")
  
  data.frame(
    sample = row.names(data),
    ent_age = as.numeric(pred_en),
    blup_age = as.numeric(pred_blup)
  )
}
