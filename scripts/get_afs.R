#########################
######### NOTES #########
# run script where data directory is
# outgroup: KRI and TJA
# correlation between lineage1 and lineage2 and on of derived afs

rm(list = ls())
.packages = c("ggplot2", "dplyr", "rstan", "tibble", "bayesplot", "purrr", "reshape2", "pracma", "viridis", "data.table",
              "Cairo", "extrafont", "ggthemes", "bbmle", "svglite", "stringi", "optparse")
.packagesdev = "thomasp85/patchwork"
# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
.instdev <- basename(.packagesdev) %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])
if(length(.packagesdev[!.instdev]) > 0) devtools::install_github(.packagesdev[!.instdev])
# Load packages into session
lapply(.packages, require, character.only=TRUE)
lapply(basename(.packagesdev), require, character.only=TRUE)

option_list = list(
  make_option(c("-C", "--csv"), type="character", default=NULL,
              help="input CSV genotype table", metavar="character"),
  make_option(c("-L", "--lineages"), type="character", default=NULL,
              help="pattern to identify the filename which contains info about lineages", metavar="character"),
  make_option(c("-E", "--extension"), type="character", default=NULL,
              help="which format to save the figures, e.g., .svg or .pdf", metavar="character"),
  make_option(c("-T", "--outxt"), type="character", default=NULL,
              help="output text file with list of fixed variants in the outgroups", metavar="character"),
  make_option(c("-O", "--outgroups"), type="character", default=NULL,
              help="define the name of the outgroup(s)", metavar="character"))

opt_parser = OptionParser(option_list=option_list,
                          description = "Identify ancestral and derived allele and compare afs between two populations",
                          epilogue = "Example: Rscript scripts/get_afs.R -C data/Clone_GATK_final_loci_edit.csv -L 6_lineages -E .svg -O KRI_sex TJA_sex")
opt = parse_args(opt_parser)

if (is.null(opt$csv) | is.null(opt$extension) | is.null(opt$outgroups)) {
  print_help(opt_parser)
  stop("All arguments must be supplied.\n", call.=FALSE)
}

# csv_gt = read.csv("data/Clone_GATK_final_loci_edit.csv")
cat("Reading input", opt$csv, "...\n")
csv_gt = read.csv(opt$csv)
# fl_line = list.files(path = "data", pattern = "5_lineages", full.names = TRUE)
if (is.null(opt$lineages)) {
  csv_all = csv_gt
} else {
  fl_line = list.files(path = "data", pattern = opt$lineages, full.names = TRUE)
  csv_line = read.csv(fl_line)
}
if (identical(csv_gt$Pop, csv_line$Pop)) {
  csv_call = cbind(csv_gt[,1:2], Line=csv_line$Line, csv_gt[, 3:ncol(csv_gt)])
}
# csv_call[1:5,1:10]
########################
######### TEST #########
# file.create("fixed_var_outgroup.txt")
dir.create("docs")
# file.create(paste0("docs/", opt$outxt))
# file.create("docs/max_val_allele1.txt")
# outg = c("KRI_sex", "TJA_sex")
outg = strsplit(opt$outgroups, split = " ")[[1]]
scaf_id = colnames(csv_call[!grepl(pattern = "Ind|Pop|Line|A2", x = colnames(csv_call))])

# oldtxt = read.table("fixed_var_outgroup.txt")
# newtxt = read.table("docs/fixed_var_outgroup.txt")
# scaf_id = setdiff(newtxt$V1, oldtxt$V1)

# get_afs function finds homozygous variants in the outgroups and compute AF
get_afs = function(csv) {
  # fl = read.table(file = txt)
  # fl = read.csv(file = csv)
  # scaf_id = colnames(csv[!grepl(pattern = "Ind|Pop|A2", x = colnames(csv))])[id_num]
  # scafID_allele = fl[, c(1, dim(fl)[2])]
  # colnames(scafID_allele) = c("scaf", "allele_num")
  # scaf_id = as.character(scafID_allele$scaf)[id_num]
  # print(scaf_id)
  id_fix_csv = lapply(scaf_id, function(x) {
    cat("Calculating frequency of", x, "...\n")
    one_scaf = csv[, colnames(csv)[grepl(pattern = x, x = colnames(csv))]]
    # group by Pop or Line
    one_scaf_pop = cbind(Pop = csv[,"Line"], one_scaf)
    one_freq_tot = data.frame(table(one_scaf_pop$Pop, one_scaf_pop[, 2]),
                              Tot = data.frame(table(one_scaf_pop$Pop))[,2]*2)
    two_freq_tot = data.frame(table(one_scaf_pop$Pop, one_scaf_pop[, 3]),
                              Tot = data.frame(table(one_scaf_pop$Pop))[,2]*2)
    if (!isempty(two_freq_tot[two_freq_tot$Var2==1,]$Var2)) {
      t3one = data.frame(Pop = one_freq_tot[one_freq_tot$Var2==1,]$Var1,
                         Allele = one_freq_tot[one_freq_tot$Var2==1,]$Var2,
                         Count = one_freq_tot[one_freq_tot$Var2==1, "Freq"] + two_freq_tot[two_freq_tot$Var2==1, "Freq"],
                         Tot = one_freq_tot[one_freq_tot$Var2==1, "Tot"])
    } else {
      t3one = data.frame(Pop = one_freq_tot[one_freq_tot$Var2==1,]$Var1,
                         Allele = one_freq_tot[one_freq_tot$Var2==1,]$Var2,
                         Count = one_freq_tot[one_freq_tot$Var2==1, "Freq"],
                         Tot = one_freq_tot[one_freq_tot$Var2==1, "Tot"])
    }
    if (!isempty(one_freq_tot[one_freq_tot$Var2==2,]$Var2)) {
      t3two = data.frame(Pop = one_freq_tot[one_freq_tot$Var2==2,]$Var1,
                         Allele = one_freq_tot[one_freq_tot$Var2==2,]$Var2,
                         Count = one_freq_tot[one_freq_tot$Var2==2, "Freq"] + two_freq_tot[two_freq_tot$Var2==2, "Freq"],
                         Tot = one_freq_tot[one_freq_tot$Var2==2, "Tot"])
      t4 = rbind(t3one, t3two)
    } else {
      colnames(two_freq_tot) = colnames(t3one)
      t4 = rbind(t3one, two_freq_tot[two_freq_tot$Allele==2,])
    }
    # This part is to double check the values of allele1 #
    # idx_max_a1 = which.max(one_scaf_pop[, 2])
    # mx_allele = paste0(x, " variant has allele1 = ", max(one_scaf_pop[, 2]), " and allele2 = ", one_scaf_pop[idx_max_a1, 3])
    # write(mx_allele, file = "docs/max_val_allele1.txt", append=TRUE)
    ######################################################
    one_af = mutate(t4, AF = Count/Tot)
    colnames(one_af) = c("Pop", "Allele", "Count", "Tot", paste0("AF_", x))
    # If AF of outgroups is below 1 remove the AF value
    if (sum(one_af[,'Pop'] %in% "outg") > 0) {
      outg = "outg"
    }
    outg_fr = one_af[which(one_af$Pop %in% outg), dim(one_af)[2]]
    outg_vc = rep(1,length(outg_fr))
    outg_filter = outg_vc[outg_fr == 1 | outg_fr == 0]
    # change sum to 2 if using lineages
    if (sum(outg_filter) == 2) {
      # Add column to be able to filter out dataframes without the extra column
      one_af$Homo_outg = 1
      # Write variant to txt output
      # line = paste0(x, " variant IS FIXED in the outgroups ", outg[1], " and ", outg[2])
      # write(line, file = paste0("docs/", opt$outxt), append=TRUE)
    }
    return(one_af)
    })
  return(id_fix_csv)
}
afs_list = get_afs(csv = csv_call)
# afs_list[[10]]
# colSums(afs_list[[101]][, c("Count", "Tot")])
afs_list_homo_outg = Filter(function(x) ncol(x)==6, afs_list)
cat("The number of polarised variants is", length(afs_list_homo_outg), "\n")

# get_daf function selects which allele is the derived since the ancestral allele has frequency = 1 in the outgroups
get_daf = function(df) {
  outg_df = df[df$Pop==outg[1], c("Count", "Tot")]
  if (outg_df[1,1]==outg_df[1,2]) {
    der_df = df[df$Allele==2, ]
  } else if (outg_df[2,1]==outg_df[2,2]) {
    der_df = df[df$Allele==1, ]
  }
  return(der_df)
}
# get_daf_line function works the same as get_daf but with lineages
get_daf_line = function(df) {
  outg_df = df[df$Pop=="outg", c("Count", "Tot")]
  if (outg_df[1,1]==outg_df[1,2]) {
    der_df = df[df$Allele==2, ]
  } else if (outg_df[2,1]==outg_df[2,2]) {
    der_df = df[df$Allele==1, ]
  }
  return(der_df)
}
if (sum(csv_call[,'Line'] %in% "outg") > 0) {
  outg = "outg"
  der_afs = lapply(afs_list_homo_outg, function(x) get_daf_line(df = x))
} else {
  der_afs = lapply(afs_list_homo_outg, function(x) get_daf(df = x))
}

# der_afs[[140]]

# Bind dataframes by columns and select only columns for the derived frequency and the column for the population names
scafs_der_afs = bind_cols(der_afs)
der_scaf = colnames(scafs_der_afs)[grepl(pattern = "scaffold", x = colnames(scafs_der_afs))]
daf_scaf = scafs_der_afs[, which(colnames(scafs_der_afs) %in% der_scaf)]
# daf_scaf[, 1:10]
fin_daf_scaf = cbind(Pop=scafs_der_afs$Pop, daf_scaf)
# fin_daf_scaf[, 1:5]
rownames(fin_daf_scaf) = fin_daf_scaf$Pop
fin_daf_scaf$Pop = NULL
fin_var_by_pop = t(fin_daf_scaf)
# head(as.data.frame(fin_var_by_pop))
# table(as.data.frame(fin_var_by_pop)[,4])

# Remove frequency of the derived allele of the outgroups because = 0
if (sum(csv_call[,'Line'] %in% "outg") > 0) {
  outg = "outg"
  fin_var_by_pop = fin_var_by_pop[, -which(colnames(fin_var_by_pop) == outg)]
} else {
  fin_var_by_pop = fin_var_by_pop[, -which(colnames(fin_var_by_pop) == outg[1])]
  fin_var_by_pop = fin_var_by_pop[, -which(colnames(fin_var_by_pop) == outg[2])]
}

# Find varinats with non-zero frequency in all populations
all_shared = fin_var_by_pop[rowSums(fin_var_by_pop > 0) == ncol(fin_var_by_pop), ]
cat("There are", length(rownames(all_shared)), "variants that have non-zero frequency in all populations.\n")

# Find varinats with non-zero frequency in sexual populations
sex_lin = fin_var_by_pop[, colnames(fin_var_by_pop)[grepl(pattern = "sex", x = colnames(fin_var_by_pop))]]
sex_shared = sex_lin[rowSums(sex_lin > 0) == ncol(sex_lin), ]
cat("There are", length(rownames(sex_shared)), "variants that have non-zero frequency in the sexual populations.\n")
cat("The correlation between", colnames(sex_shared)[1], "and", colnames(sex_shared)[2], "frequencies is",
    cor(sex_shared[,1], sex_shared[,2]), "\n")
# plot(sex_shared[,1], sex_shared[,2])

# Find varinats with non-zero frequency in all clonal populations
clo_lin = fin_var_by_pop[, colnames(fin_var_by_pop)[!grepl(pattern = "sex", x = colnames(fin_var_by_pop))]]
# head(clo_lin)
clo_shared = clo_lin[rowSums(clo_lin > 0) == ncol(clo_lin), ]
cat("There are", length(rownames(clo_shared)), "variants that have non-zero frequency in all the clonal populations.\n")

# intersect(rownames(all_shared), rownames(clo_shared))
# intersect(rownames(all_shared), rownames(sex_shared))
cat("The variants with non-frequency that are shared between sexual and clonal populations are:\n")
intersect(rownames(sex_shared), rownames(clo_shared))

# sis_lin = fin_var_by_pop[, colnames(fin_var_by_pop)[grepl(pattern = "sis|Sis", x = colnames(fin_var_by_pop))]]
# head(sis_lin)
# sis_shared = sis_lin[rowSums(sis_lin > 0) == ncol(sis_lin), ]

df_grid = expand.grid(popx = seq(from=0,to = 1,by = 0.1), popy = seq(from=0,to = 1,by = 0.1))
pl_comb = combn(x = 1:ncol(fin_var_by_pop), m = 2)
afs_2d_heat = function(dataf, n_comb) {
  df = as.data.frame(dataf)[, pl_comb[, n_comb]]
  df = round(df, 1)
  df$Count = apply(df, MARGIN = 1, FUN = function(x) {
    nrow(df[df[, 1]==x[1] & df[, 2]==x[2],])
  })
  dfu = unique(df)
  colnames(dfu) = c("popx", "popy", "Count")
  df_merge = merge(dfu, df_grid, by = c("popx", "popy"), all = TRUE)
  df_merge[is.na(df_merge)] = 0
  cat("Heatmapping", colnames(dataf)[pl_comb[1, n_comb]], "and", colnames(dataf)[pl_comb[2, n_comb]], "...\n")
  one_fig = ggplot(df_merge, aes(x = popx, y = popy)) +
    geom_tile(data = subset(df_merge, Count!=0), aes(fill = Count)) +
    geom_tile(data = subset(df_merge,  Count==0), aes(colour = "0"),
              linetype = 0, fill = "#330033") +
    scale_fill_continuous(type = "viridis") +
    labs(x = colnames(dataf)[pl_comb[1, n_comb]],
         y = colnames(dataf)[pl_comb[2, n_comb]]) +
    theme(legend.title = element_blank())
  return(one_fig)
}
# ext = ".svg"
ext = opt$extension
lapply(1:ncol(pl_comb), function(x) {
  cat("Saving heatmap",
      paste0("TheFucusProject/figures/heatmap_der_", colnames(fin_var_by_pop)[pl_comb[1, x]], "_VS_", colnames(fin_var_by_pop)[pl_comb[2, x]], ext),
      "...\n")
  ggsave(filename = paste0("TheFucusProject/figures/heatmap_der_", colnames(fin_var_by_pop)[pl_comb[1, x]],
                           "_VS_", colnames(fin_var_by_pop)[pl_comb[2, x]], ext),
         plot = afs_2d_heat(dataf = fin_var_by_pop, n_comb = x), width = 10, height = 8)
})
