kinship.table <- read.table("intermediate_sh/merged_qc_unrel_kin.kin0", header = F)[,c(1:2,6)] # Reading kinship table
  names(kinship.table) <- c("ID1", "ID2", "kinship")

### Checking Kinship between individuals with threshold (0.0442)
if(!any(kinship.table[,3] >= 0.0442)) {
  cat("*** Checking Done! All Remained Individuals are Unrelated (Kinship Coefficient < 0.0442) ***",sep = "\n")
} else {
    cat("Checking Done! Related Individuals found (Kinship Coefficient >= 0.0442):", sep = "\n")
    kinship.table[kinship.table$kinship >= 0.0442,]

}
