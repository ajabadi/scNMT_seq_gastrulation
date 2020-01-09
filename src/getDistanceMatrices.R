source("src/utils/utils.R")

metFile <- grep(list.files(recursive = TRUE), pattern = "acc_dt_list.rds", value = TRUE)
accFile <- grep(list.files(recursive = TRUE), pattern = "met_dt_list.rds", value = TRUE)

stopifnot(length(metFile) == 1 & length(accFile) == 1)

met_dt_list <- readRDS(metFile)
acc_dt_list <- readRDS(accFile)

log_file <- "output/get_dist_mats.log"
# subset <- ifelse(Sys.info()["user"] == "alabadi", 10000, NULL)
met_results <- tryCatch({get_dist_mats(met_dt_list)},error = function(e) {
  cat(e$message, file=log_file, append = TRUE)
  cat(e$call, file=log_file, append = TRUE)
})

if(!is.null(met_results)) {
  saveRDS(met_results, file = "output/met_wdist_mats.rds")
}

acc_results <- tryCatch({get_dist_mats(acc_dt_list)},error = function(e) {
  cat(e$message, file=log_file, append = TRUE)
  cat(e$call, file=log_file, append = TRUE)
})

if(!is.null(acc_results)) {
  saveRDS(acc_results, file = "output/acc_wdist_mats.rds")
}
