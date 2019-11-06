########################
######### TODO #########
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
  make_option(c("-T", "--outxt"), type="character", default=NULL,
              help="output text file with list of fixed variants in the outgroups", metavar="character"),
  make_option(c("-O", "--outgroups"), type="character", default=NULL,
              help="define the name of the outgroup(s)", metavar="character"))

opt_parser = OptionParser(option_list=option_list,
                          description = "Identify ancestral and derived allele and compare afs between two populations",
                          epilogue = "Example: Rscript scripts/get_afs.R -C data/Clone_GATK_final_loci_edit.csv -T fixed_var_outgroup.txt -O KRI_sex TJA_sex")
opt = parse_args(opt_parser)

if (is.null(opt$csv) | is.null(opt$outxt) | is.null(opt$outgroups)) {
  print_help(opt_parser)
  stop("All arguments must be supplied.\n", call.=FALSE)
}

# csv_call = read.csv("data/Clone_GATK_final_loci_edit.csv")
cat("Reading input", opt$csv, "...\n")
csv_call = read.csv(opt$csv)
# csv_call[1:5,1:10]
########################
######### TEST #########
# file.create("fixed_var_outgroup.txt")
dir.create("docs")
file.create(paste0("docs/", opt$outxt))
# file.create("docs/max_val_allele1.txt")
# outg = c("KRI_sex", "TJA_sex")
outg = strsplit(opt$outgroups, split = " ")[[1]]
# csv_call = csv_call[,1:10]

scaf_id = colnames(csv_call[!grepl(pattern = "Ind|Pop|A2", x = colnames(csv_call))])

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
    one_scaf_pop = cbind(Pop = csv[,2], one_scaf)
    one_freq_tot = data.frame(table(one_scaf_pop$Pop, one_scaf_pop[, 3]),
                              Tot = data.frame(table(one_scaf_pop$Pop))[,2])
    # This part is to double check the values of allele1 #
    # idx_max_a1 = which.max(one_scaf_pop[, 2])
    # mx_allele = paste0(x, " variant has allele1 = ", max(one_scaf_pop[, 2]), " and allele2 = ", one_scaf_pop[idx_max_a1, 3])
    # write(mx_allele, file = "docs/max_val_allele1.txt", append=TRUE)
    ######################################################
    one_af = mutate(one_freq_tot, AF = Freq/Tot)
    colnames(one_af) = c("Pop", "Allele", "Count", "Tot", paste0("AF_", x))
    # If AF of outgroups is below 1 remove the AF value
    outg_fr = one_af[which(one_af$Pop %in% outg), dim(one_af)[2]]
    outg_vc = rep(1,length(outg_fr))
    outg_filter = outg_vc[outg_fr == 1 | outg_fr == 0]
    if (sum(outg_filter) == 4) {
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
afs_list_homo_outg = Filter(function(x) ncol(x)==6, afs_list)
cat("The number of polarised variants is", length(afs_list_homo_outg), "\n")

# scaf_call = csv_call[, colnames(csv_call)[grepl(pattern = "scaffold12659_12459", x = colnames(csv_call))]]
# scaf_pop_call = cbind(Pop = csv_call[,2], scaf_call)

get_daf = function(df) {
  outg_df = df[df$Pop==outg[1], c("Count", "Tot")]
  if (outg_df[1,1]==outg_df[1,2]) {
    der_df = df[df$Allele==2, ]
  } else if (outg_df[2,1]==outg_df[2,2]) {
    der_df = df[df$Allele==1, ]
  }
  return(der_df)
}
der_afs = lapply(afs_list_homo_outg, function(x) get_daf(df = x))
der_afs[[100]]

scafs_der_afs = bind_cols(der_afs)
der_scaf = colnames(scafs_der_afs)[grepl(pattern = "scaffold", x = colnames(scafs_der_afs))]
daf_scaf = scafs_der_afs[, which(colnames(scafs_der_afs) %in% der_scaf)]
daf_scaf[, 1:10]
fin_daf_scaf = cbind(Pop=scafs_der_afs$Pop, daf_scaf)
fin_daf_scaf[, 1:5]
