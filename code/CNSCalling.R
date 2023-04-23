
# CNS calling

library(sigminer)


tally_W_pcawg <- readRDS("./data/tallydata/tally_W_pcawg.rds")
tally_W_560 <- readRDS("./data/tallydata/tally_W_560.rds")

nmfpcawg <- tally_W_pcawg$nmf_matrix
nmf560 <- tally_W_560$nmf_matrix

nmfall <- rbind(nmfpcawg, nmf560)

# Analysis of mutational signatures
sigprofiler_extract(nmfall,
                    output = "sig_all",
                    range = 2:30,
                    nrun = 100,
                    init_method = "random",
                    is_exome = FALSE,
                    use_conda = TRUE
)


# Signature number determination

sig_all <- sigprofiler_import("./data/modelselection/sig_CNS/sig_all/", order_by_expo = TRUE, type = "all")

# saveRDS(sig_all, file = "./data/modelselection/data/all_cn_solutions_sp.rds")


sig_all <- readRDS("./data/modelselection/data/all_cn_solutions_sp.rds")

show_sig_number_survey(
  sig_all$all_stats %>%
    dplyr::rename(
      s = `Stability (Avg Silhouette)`,
      e = `Mean Cosine Distance`
    ) %>%
    dplyr::mutate(
      SignatureNumber = as.integer(gsub("[^0-9]", "", SignatureNumber))
    ),
  x = "SignatureNumber",
  left_y = "s", right_y = "e",
  left_name = "Stability (Avg Silhouette)",
  right_name = "Mean Cosine Distance",
  highlight = 8
)


all_sigs <- sig_all$solution_list$S8
apply(all_sigs$Exposure, 1, mean)
#      Sig1      Sig2      Sig3      Sig4      Sig5      Sig6      Sig7      Sig8
# 237.96494 148.05574 135.60024 125.57896 116.38957 103.03506  96.67036  75.00599

sig_names(all_sigs)
# [1] "Sig1" "Sig2" "Sig3" "Sig4" "Sig5" "Sig6" "Sig7" "Sig8"

colnames(all_sigs$Signature) <- colnames(all_sigs$Signature.norm) <- rownames(all_sigs$Exposure)
rownames(all_sigs$Exposure.norm) <- all_sigs$Stats$signatures$Signatures <- paste0("CNS", 1:8)

sig_names(all_sigs)


# saveRDS(all_sigs, file = "./data/modelselection/data/all_cn_sigs_8_signature.rds")

















































