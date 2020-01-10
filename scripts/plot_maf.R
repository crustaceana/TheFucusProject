rm(list = ls())
summ_tsv = read.table("data/summary.stats.markers.pop.tsv", header = TRUE, sep = "\t")
head(summ_tsv)
csv_gt = read.csv("data/Clone_GATK_final_loci_edit.csv")
csv_gt[1:9,1:9]
unique(csv_gt$Ind)
table(csv_gt$Pop)
unique(summ_tsv$POP_ID)

fl_line = list.files(path = "data", pattern = "6_lineages", full.names = TRUE)
csv_line = read.csv(fl_line)
head(csv_line)

table(summ_tsv$POP_ID)
table(csv_line$Pop)
tsv_popid = unique(substr(summ_tsv$POP_ID, start = 1, 3))
csv_popid = unique(substr(csv_line$Pop, start = 1, 3))

summ_tsv$pop_px = substr(summ_tsv$POP_ID, start = 1, 3)
head(summ_tsv)
csv_line$pop_px = substr(csv_line$Pop, start = 1, 3)
head(csv_line)
csv_line$pop_px[csv_line$pop_px=="Sup"] = "VAL"
csv_line$pop_px[csv_line$pop_px=="Sis"] = "KUG"
csv_line$pop_px[csv_line$pop_px=="Mot"] = "SAL"
table(csv_line$pop_px)
table(summ_tsv$pop_px)

pop_n = 11
tsv_popid[pop_n]
(line_nm = as.character(unique(csv_line$Line)))
summ_tsv$LINE[summ_tsv$pop_px==tsv_popid[pop_n]] = unique(as.character(csv_line$Line[csv_line$pop_px==tsv_popid[pop_n]]))
head(summ_tsv[summ_tsv$pop_px==tsv_popid[pop_n], ])
table(summ_tsv$LINE)[1] + table(summ_tsv$LINE)[2] + table(summ_tsv$LINE)[3] + table(summ_tsv$LINE)[4] + table(summ_tsv$LINE)[5] +
  table(summ_tsv$LINE)[6] + table(summ_tsv$LINE)[7]

head(summ_tsv)
table(summ_tsv$LINE)
line_par = c("line_clone_a", "line_clone_b")
summ_tsv_par = summ_tsv[which(summ_tsv$LINE %in% line_par), ]
# test_par = sample_n(summ_tsv_par, size = 50)
test_par = summ_tsv_par
# df_grid = expand.grid(popx = seq(from=0,to = 1,by = 0.1), popy = seq(from=0,to = 1,by = 0.1))
bin_fq = seq(0,1,0.05)
test_par$bin_maf = cut(test_par$MAF_LOCAL, bin_fq, include.lowest = TRUE)
head(test_par)
unique(test_par$bin_maf)
test_x = data.frame(bin_maf=unique(test_par$bin_maf))
test_bin_m = aggregate(x = test_par$MAF_LOCAL, by=list(LINE=test_par$LINE, bin_maf=test_par$bin_maf), FUN=mean)
test_pl = merge(merge(test_x, test_bin_m[test_bin_m[,1]==line_par[1], ], by = "bin_maf", all = TRUE),
                test_bin_m[test_bin_m[,1]==line_par[2], ], by = "bin_maf", all = TRUE)
ggplot(data = test_pl) +
  geom_point(aes(x = x.x, y = x.y))
