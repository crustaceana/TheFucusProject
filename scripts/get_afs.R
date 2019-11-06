########################
######### TODO #########
# outgroup: KRI and TJA
# correlation between line1 and line2 derived afs 

rm(list = ls())
.packages = c("ggplot2", "dplyr", "rstan", "tibble", "bayesplot", "purrr", "reshape2", "pracma", "viridis", "data.table",
              "Cairo", "extrafont", "ggthemes", "bbmle", "svglite", "stringi")
.packagesdev = "thomasp85/patchwork"
# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
.instdev <- basename(.packagesdev) %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])
if(length(.packagesdev[!.instdev]) > 0) devtools::install_github(.packagesdev[!.instdev])
# Load packages into session
lapply(.packages, require, character.only=TRUE)
lapply(basename(.packagesdev), require, character.only=TRUE)

csv_call = read.csv("data/Clone_GATK_final_loci_mean_inform.csv")
csv_call[1:5,1:10]
# colnames(csv_call)[grepl(pattern = "scaffold", colnames(csv_call))] =
#   paste0(colnames(csv_call)[grepl(pattern = "scaffold", colnames(csv_call))], "_A1")
colnames(csv_call)[grepl(pattern = "X", colnames(csv_call))] =
  paste0(colnames(csv_call)[grepl(pattern = "scaffold", colnames(csv_call))], "_A2")

write.csv(x = csv_call, file = "data/Clone_GATK_final_loci_edit.csv", row.names = FALSE)
csv_call = read.csv("data/Clone_GATK_final_loci_edit.csv")
csv_call[1:5,1:10]
data.frame(table(factor(csv_call$Pop)))
sum(grepl(pattern = "scaffold", x = colnames(csv_call)))/2
sum(grepl(pattern = "scaffold", x = colnames(csv_call))) * 2
colSums(csv_call[,-1:-2]) < 16383

########################
######### TEST #########
file.create("fixed_var_outgroup.txt")
outg = c("KRI_sex", "TJA_sex")
# csv_call = csv_call[,1:10]
scafID = colnames(csv_call[!grepl(pattern = "Ind|Pop|A2", x = colnames(csv_call))])
get_daf = function(csv, scaf) {
  scaf_df = csv[, colnames(csv)[grepl(pattern = scaf, x = colnames(csv))]]
  scaf_pop_call = cbind(Pop = csv[,2], scaf_df)
  mean_by_pop = data.frame(table(factor(csv$Pop)),
                           aggregate(scaf_pop_call, by = list(scaf_pop_call$Pop), FUN = mean)[,3:4])
  outg_df = mean_by_pop[which(mean_by_pop$Var1 %in% outg), ]
  if (sum(rowSums(outg_df[,3:4]) == 2) == 2) {
    line = paste0(scaf, " variant IS FIXED in the outgroups ", outg[1], " and ", outg[2], " with sum = 2")
    write(line, file = "fixed_var_outgroup.txt", append=TRUE)
  #   anc_df = data.frame(lapply(mean_by_pop[,3:4], function(x) gsub(pattern = 1, replacement = "ANC", x = x)))
  #   anc_df = data.frame(lapply(anc_df, function(x) gsub(pattern = 2, replacement = "DER", x = x)))
  } else if (sum(rowSums(outg_df[,3:4]) == 4) == 2) {
    line = paste0(scaf, " variant IS FIXED in the outgroups ", outg[1], " and ", outg[2], " with sum = 4")
    write(line, file = "fixed_var_outgroup.txt", append=TRUE)
  # } else {
  #   cat(scaf, "variant is NOT FIXED in the outgroups", outg, "\n")
  }
  # return(outg_df)
}
lapply(scafID, function(s) {
  get_daf(csv = csv_call, scaf = s)
})

scaf_call = csv_call[, colnames(csv_call)[grepl(pattern = "scaffold6531_5857", x = colnames(csv_call))]]
scaf_pop_call = cbind(Pop = csv_call[,2], scaf_call)
head(scaf_pop_call)
scaf_pop_01 = data.frame(lapply(scaf_pop_call, function(x) gsub(pattern = 2, replacement = 0, x = x)))
head(scaf_pop_01)
data.frame(table(scaf_pop_01$Pop, scaf_pop_01$scaffold6531_5857_A2))
table(scaf_pop_01$Pop)
aggregate(scaf_pop_01, by = list(scaf_pop_01$Pop), FUN = mean)
freq_by_pop = data.frame(table(factor(scaf_pop_01$Pop)),
                         aggregate(scaf_pop_01, by = list(scaf_pop_01$Pop), FUN = mean)[,3:4])
if (sum(rowSums(outg_df[,3:4]) == 2) == 2) {
  anc_df = data.frame(lapply(mean_by_pop[,3:4], function(x) gsub(pattern = 1, replacement = "ANC", x = x)))
  anc_df = data.frame(lapply(anc_df, function(x) gsub(pattern = 2, replacement = "DER", x = x)))
}



str(scaf_pop_call)
aggregate(scaf_pop_call, by = list(scaf_pop_call$Pop), FUN = mean)
mean_by_pop = data.frame(table(factor(csv_call$Pop)), aggregate(scaf_pop_call, by = list(scaf_pop_call$Pop), FUN = mean)[,3:4])

outg_df = mean_by_pop[which(mean_by_pop$Var1 %in% outg), ]
if (sum(rowSums(outg_df[,3:4]) == 2) == 2) {
  anc_df = data.frame(lapply(mean_by_pop[,3:4], function(x) gsub(pattern = 1, replacement = "ANC", x = x)))
  anc_df = data.frame(lapply(anc_df, function(x) gsub(pattern = 2, replacement = "DER", x = x)))
}

vcf_call = read.table("test/data/Clonal_GATK_subvcf.tsv", header = TRUE)
vcf_call[1:12, 1:12]
rm_col = c("INFO", "ID", "QUAL", "FILTER", "FORMAT")
vcf_call = vcf_call[,-which(names(vcf_call) %in% rm_col)]

sapply(vcf_call$FX14KRI.1, function(x) strsplit(as.character(x), split = ":"))
get_GT = apply(vcf_call[, c(-1:-4)], MARGIN = 2, FUN = function(x) {
  substr(x, start = 1, stop = 3)
})
get_GT[1:5,1:5]
sex_GT = get_GT[, grepl(pattern = "KRI|TJA|STO|LET", x = colnames(get_GT))]
colnames(sex_GT)
sex_GT[1:5,1:5]
sex_GT = gsub(pattern = "1/1", replacement = 2, x = sex_GT)
sex_GT = gsub(pattern = "0/1", replacement = 1, x = sex_GT)
sex_GT = gsub(pattern = "0/0", replacement = 0, x = sex_GT)
sex_GT = apply(sex_GT, MARGIN = 2, FUN = as.integer)
rowSums(x = sex_GT, na.rm = TRUE)

clo_GT = get_GT[, !grepl(pattern = "KRI|TJA|STO|LET", x = colnames(get_GT))]
colnames(clo_GT)
clo_GT[1:5,1:5]
clo_GT = gsub(pattern = "1/1", replacement = 2, x = clo_GT)
clo_GT = gsub(pattern = "0/1", replacement = 1, x = clo_GT)
clo_GT = gsub(pattern = "0/0", replacement = 0, x = clo_GT)
clo_GT = apply(clo_GT, MARGIN = 2, FUN = as.integer)
rowSums(x = clo_GT, na.rm = TRUE)
