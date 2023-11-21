setwd("~/data/mm10_dysf/")
dir.create("scrna_db", showWarnings = F)
dir.create("figs", showWarnings = F)
source("mmdysf_utils_byst.r")
source("transfer_color_chimera_tmp.r")
options("stringsAsFactors" = TRUE)
mmdysf_vars()
rl()

library("reticulate")
use_python("/net/mraid14/export/data/users/eladch/tools/CO7/python3/3.7.5/bin/python", required=TRUE)
np = import("tgutils.numpy")
skl = import("sklearn")
nmf_py = import("sklearn.decomposition._nmf")
#pacmap = import("pacmap")
importlib = import("importlib")
importlib$reload(nmf_py)
source("~/src/mc_local_nmf_cv.R")
source("~/src/unite_metacell.R")
library(RColorBrewer)
library(wordcloud)

source("~/src/unite_metacell.R")
source("~/src/mc_local_nmf_gene_programs.R")
#source("workspace/mc_local_nmf_cv.R")
umap_py = import("umap")

rl_nmf = function(){
  importlib$reload(nmf_py)
}
rl_nmf()

library(ggplot2)
library(ggpubr)
# scdb_flow_init()
# library("devtools")
# load_all("/home/akhiad/src/metacell.flow")
# tgconfig::override_params("config/netflow.yaml","metacell")
# First, generate the metacell models with mmdysf_data_gen function
# Then probably write a custom plotting function, similar to pd1_time_course_plots

set_param(param = "mc_plot_device", value = "svg", package = "metacell")
# build mat
mmdysf_build_mat = function(force_new=T, max_f_mit=0.4, min_umis=450, max_umis=2^15, sample_sheet_fn="config/DysfunctionalSampleIndex_09Sep20.txt", sample_sheet_fn2="config/DysfunctionalSampleIndex_18Jun20bk.txt") 
{
  # load initial mat
  mcell_import_multi_mars(mat_nm=mm_all_id, dataset_table_fn=sample_sheet_fn, base_dir="/net/mraid14/export/data/users/akhiad/mm10_dysf/umi.tab", force=force_new) # md_filter = "!Amp.Batch.ID %in% c('AB10225.txt', 'AB10226.txt', 'AB10227.txt', 'AB10228.txt', 'AB10229.txt')"
  
  mcell_plot_batch_stats(mm_all_id)
  min_umis_cutoff = mcell_plot_umis_per_cell(mm_all_id, min_umis_cutoff=min_umis)
  mat_a = scdb_mat(mm_all_id)
  gene_sums = rowSums(mat_a@mat)
  if (sum(gene_sums == 0) > 0){
    print("zero sum genes")
    mcell_mat_ignore_genes(mm_all_id, mm_all_id, names(gene_sums[gene_sums == 0]))
  }
  
  mcell_add_gene_stat(mm_all_id, mm_all_id, force=force_new)
  
  # add FACS index sort data
  source("mmdysf_csv_to_tabs.r")
  facs_csv_to_tab(idir = "idx_sort/input_3", odir = "idx_sort/processed_3", sample_sheet = sample_sheet_fn, colnames_res=c("Well", "FSC-A", "FSC-W", "FSC-H", "SSC-A", "SSC-W", "SSC-H", "mCherry-A", "Pacific-Blue-A", "APC-A", "PerCP-Cy5-5-A", "PE-Cy7-A", "GFP-A", "bystander-GFP-A", "PD1-APC-A", "CD45-1-APC-Cy7-A", "CD45-2-PE-A", "TIM3-BV421-A", "CD8-PE-Cy7-A", "time", "Time", "TIM3", "PD1", "Tigit", "Sell", "CD25", "CTLA-4", "Fasl", "CCR5", "TNFrsfa", "CCR2", "SlamF6", "IL7R", "CD45.1", "CD45.2", "CCR7", "CD8a"))
  mcell_add_mars_facs_data(new_mat_id = mm_all_id, mat_id = mm_all_id, base_dir = "idx_sort/processed_3")
  
  mat_a = scdb_mat(mm_all_id)
  
  mat_byst = scdb_mat(paste(mm_all_id, "_byst", sep=""))
  byst_md = mat_byst@cell_metadata
  new_md = mat_a@cell_metadata
  new_md[setdiff(colnames(byst_md), colnames(mat_a@cell_metadata))] = NA
  new_md[rownames(byst_md),colnames(byst_md)] = byst_md
  mat_a@cell_metadata = new_md
  levels(mat_a@cell_metadata$mouse_id) = c(levels(mat_a@cell_metadata$mouse_id), "m285")
  mat_a@cell_metadata[mat_a@cell_metadata$amp_batch_id == "AB11219","mouse_id"] = "m285"
  scdb_add_mat(mm_all_id, mat_a)
  
  mat_a@cell_metadata$bystander.GFP.A_Ab[mat_a@cell_metadata$location %in% c('spleen', 'LN')] = mat_a@cell_metadata$GFP.A_Ab[mat_a@cell_metadata$location %in% c('spleen', 'LN')]
  
  table(mat_a@cell_metadata[colnames(mat_a@mat), "cell_type"])
  mat_a@cell_metadata$bystander.GFP.A_Ab[is.na(mat_a@cell_metadata$bystander.GFP.A_Ab)] = 0
  mat_a@cell_metadata$cell_type[mat_a@cell_metadata$bystander.GFP.A_Ab > 4500] = "byst"
  mat_a@cell_metadata$cell_type = plyr::revalue(mat_a@cell_metadata$cell_type, c("spec_or_byst"="spec"))
  table(mat_a@cell_metadata[colnames(mat_a@mat), "cell_type"])
  #scdb_add_mat(mm_all_id, mat_a)
  
  mat_a@cell_metadata$CD45.1_Ab[is.na(mat_a@cell_metadata$CD45.1_Ab)] = 0
  mat_a@cell_metadata$CD45.1_Ab[is.na(mat_a@cell_metadata$CD45.2_Ab)] = 0
  plot(sort(mat_a@cell_metadata$CD45.1_Ab))
  abline(h=1500)
  plot(sort(log10(mat_a@cell_metadata$CD45.2_Ab[log10(mat_a@cell_metadata$CD45.1_Ab)<3])))
  mat_a@cell_metadata$cell_type[mat_a@cell_metadata$CD45.1_Ab > 1500] = "spec"
  #mat_a@cell_metadata$cell_type[mat_a@cell_metadata$CD45.2_Ab > 4500] = "endo"
  mat_a@cell_metadata$cell_type = plyr::revalue(mat_a@cell_metadata$cell_type, c("spec or endo"="endo"))
  
  
  mat_a@cell_metadata$batch_set_id_saved = mat_a@cell_metadata$batch_set_id
  mat_a@cell_metadata$batch_set_id = plyr::revalue(mat_a@cell_metadata$batch_set_id, c("PD1_2"="PD1_5", "PD1_6"="PD1_5", "Ctrl_2"="Ctrl_5", "Ctrl_6"="Ctrl_5", "InVitro_activated2"="Ctrl_5","InVitro_activated"="Ctrl_5"))
  #scdb_add_mat(mm_all_id, mat_a)
  
  mat_a@cell_metadata$days_post_transfer[mat_a@cell_metadata$days_post_transfer==""] = "1"
  mat_a@cell_metadata$days_post_transfer = factor(as.character(mat_a@cell_metadata$days_post_transfer), levels=as.character(levels(mat_a@cell_metadata$days_post_transfer)[levels(mat_a@cell_metadata$days_post_transfer) != ""]))
  mat_a@cell_metadata$days_post_transfer = plyr::revalue(mat_a@cell_metadata$days_post_transfer, c("16"="15", "13"="12"))
  
  plot(sort(log10(as.integer(mat_a@cell_metadata$CD45.1.APC.Cy7.A_Ab[mat_a@cell_metadata$batch_set_id %in% c("Delay_OT1_GFPOT1_ctrl", "Delay_OT1_GFPOT1_PD1")]))))
  mat_a@cell_metadata$CD45.1.APC.Cy7.A_Ab[is.na(mat_a@cell_metadata$CD45.1.APC.Cy7.A_Ab)] = 0
  levels(mat_a@cell_metadata$cell_type) = c(levels(mat_a@cell_metadata$cell_type), "delay_spec")
  mat_a@cell_metadata$cell_type[mat_a@cell_metadata$days_post_transfer == "6 or 10" & mat_a@cell_metadata[,"CD45.1.APC.Cy7.A_Ab"] > 400] = "delay_spec"
  mat_a@cell_metadata$days_post_transfer[mat_a@cell_metadata[,"days_post_transfer"]=="6 or 10" & mat_a@cell_metadata[,"CD45.1.APC.Cy7.A_Ab"] > 400] = "10"
  mat_a@cell_metadata$days_post_transfer = plyr::revalue(mat_a@cell_metadata$days_post_transfer, c("6 or 10"="6"))
  
  mat_a@cell_metadata$mCherry.A_Ab = as.numeric(mat_a@cell_metadata$mCherry.A_Ab)
  mat_a@cell_metadata$mCherry.A_Ab[mat_a@cell_metadata$mCherry.A_Ab<0] = 0
  mat_a@cell_metadata$mCherry.A_Ab[is.na(mat_a@cell_metadata$mCherry.A_Ab)] = 0
  plot(sort(log10(mat_a@cell_metadata$mCherry.A_Ab)))
  abline(h=log10(2000))
  mat_a@cell_metadata$PD1_KO = "NotKO"
  mat_a@cell_metadata$PD1_KO[mat_a@cell_metadata$mCherry.A_Ab > 2000] = "KO"
  levels(mat_a@cell_metadata$cell_type) = c(levels(mat_a@cell_metadata$cell_type), "pd1_ko_spec", "ctrl_ko_spec")
  mat_a@cell_metadata$cell_type[mat_a@cell_metadata$mCherry.A_Ab > 2000 & mat_a@cell_metadata$batch_set_id %in% c("PD1_ko_pd1", "PD1_ko_ctrl")] = "pd1_ko_spec"
  mat_a@cell_metadata$cell_type[mat_a@cell_metadata$mCherry.A_Ab <= 2000 & mat_a@cell_metadata$batch_set_id %in% c("PD1_ko_pd1", "PD1_ko_ctrl")] = "ctrl_ko_spec"
  
  # levels(mat_at@cell_metadata$batch_set_id) = c(levels(mat_at@cell_metadata$batch_set_id), c("PD1_7_cd25+", "PD1_7_cd25-"))
  # mat_at@cell_metadata[mat_at@cell_metadata[,"amp_batch_id"]=="AB10120","batch_set_id"] = "PD1_7_cd25+"
  # mat_at@cell_metadata[mat_at@cell_metadata[,"amp_batch_id"]=="AB10121","batch_set_id"] = "PD1_7_cd25-"
  # 
  # levels(mat_at@cell_metadata$cell_type) = c(levels(mat_at@cell_metadata$cell_type), c("spec_cd25+", "spec_cd25-", "spec_tim3+pd1+"))
  # mat_at@cell_metadata[mat_at@cell_metadata[,"amp_batch_id"]=="AB10120","cell_type"] = "spec_cd25+"
  # mat_at@cell_metadata[mat_at@cell_metadata[,"amp_batch_id"]=="AB10121","cell_type"] = "spec_cd25-"
  # mat_at@cell_metadata[mat_at@cell_metadata[,"amp_batch_id"]=="AB10122","cell_type"] = "spec_tim3+pd1+"
  # 
  scdb_add_mat(mm_all_id, mat_a)
  # 
  # mat_f = scdb_mat(mm_filt_id)
  # levels(mat_f@cell_metadata$batch_set_id) = c(levels(mat_f@cell_metadata$batch_set_id), c("PD1_7_cd25+", "PD1_7_cd25-"))
  # levels(mat_f@cell_metadata$batch_set_id_saved) = c(levels(mat_f@cell_metadata$batch_set_id_saved), c("PD1_7_cd25+", "PD1_7_cd25-"))
  # levels(mat_f@cell_metadata$cell_type) = c(levels(mat_f@cell_metadata$cell_type), c("spec_cd25+", "spec_cd25-", "spec_tim3+pd1+"))
  # 
  # mat_f@cell_metadata[rownames(mat_at@cell_metadata),] = mat_at@cell_metadata
  # scdb_add_mat(mm_filt_id, mat_f)
  
  
  # remove cells with high %mito and mitochondrial genes
  mt_genes = grep('^mt-', mat_a@genes, perl=T, v=T)
  f_mt = colSums(mat_a@mat[mt_genes, ]) / colSums(mat_a@mat)
  mcell_mat_ignore_cells(mm_filt_id, mm_all_id, names(which(f_mt >= max_f_mit)))
  mcell_mat_ignore_genes(mm_filt_id, mm_filt_id, unique(c(mt_genes, "Malat1", "Xist", grep('.*Rik', mat_f@genes, v=T, perl=T), grep('^Gm[0-9].*', mat_f@genes, v=T, perl=T))))
  mcell_mat_ignore_small_cells(mm_filt_id, mm_filt_id, min_umis)
  mat = scdb_mat(mm_filt_id)
  csize = colSums(mat@mat)
  large_c = names(which(csize > max_umis))
  mcell_mat_ignore_cells(mm_filt_id, mm_filt_id, union(mat@ignore_cells, large_c))
  
  
  mcell_add_gene_stat(mm_filt_id, mm_filt_id, force=force_new)
  mcell_gset_filter_varmean(gset_id=mm_filt_id, gstat_id=mm_filt_id, T_vm=0.1, force_new=T)
  mcell_gset_filter_cov(gset_id = mm_filt_id, gstat_id=mm_filt_id, T_tot=100, T_top3=2)
  mcell_plot_gstats(gstat_id=mm_filt_id, gset_id=mm_filt_id)
  mcell_gset_split_by_dsmat(mm_filt_id, mm_filt_id, 150)
  mcell_plot_gset_cor_mats(mm_filt_id, mm_filt_id)
  
  mat = scdb_mat(mm_filt_id)
  
  rl()
  
  md_info = mat@cell_metadata[mat@cells, ] %>% 
    group_by(batch_set_id_saved, Gender, days_post_transfer, location, cell_type) %>% 
    summarize(n_mice=length(unique(mouse_id)), n_plates=length(unique(amp_batch_id)), n_cells=length(amp_batch_id)) %>%
    data.frame 
  write.table(md_info, scfigs_fn(mm_filt_id, "cells_md_counts", ext="txt"), quote=F, sep="\t")
  
}

mmdysf_build_lat_gsets = function() 
{
  # build gset to mask by the data
  select_gene_modules_by_anchor_genes(mm_filt_id, c('Mki67', 'Hist1h1d', 'Pcna', 'Smc4', 'Mcm3'), gset_nm="mmdysf_cc", cor_thresh=0.1, sz_cor_thresh=0.1)
  
  select_gene_modules_by_anchor_genes(mm_filt_id, c('Isg15', 'Wars', 'Ifit1'), gset_nm="mmdysf_ifn", cor_thresh=0.1, sz_cor_thresh=0.1)
  select_gene_modules_by_anchor_genes(mm_filt_id, c("Hsp90ab1", "Hspa1a", "Fos", "Hif1a"), gset_nm="mmdysf_stress", cor_thresh=0.1, sz_cor_thresh=0.1)
  
  # Stopped at this point to manually select the following gene clusters to filter: based on invitro activated and naive transfer plates (amp_batches 5837 to 6243)
  # mcell_gset_remove_clusts("mmdysf_cc", filt_clusts=c(10, 12, 20), new_id = "mmdysf_cc_filt", reverse=T)
  # mcell_gset_remove_clusts("mmdysf_ifn", filt_clusts=17, new_id = "mmdysf_ifn_filt", reverse=T)
  # mcell_gset_remove_clusts("mmdysf_stress", filt_clusts=c(1, 20), new_id = "mmdysf_stress_filt", reverse=T)
  # mcell_gset_remove_clusts("mmdysf_cc", filt_clusts=c(1,8,11,12,18,19,20), new_id = "mmdysf_cc_filt", reverse=T)
  # mcell_gset_remove_clusts("mmdysf_ifn", filt_clusts=c(14,17), new_id = "mmdysf_ifn_filt", reverse=T)
  # mcell_gset_remove_clusts("mmdysf_stress", filt_clusts=c(16,18,19,20), new_id = "mmdysf_stress_filt", reverse=T)
  mcell_gset_remove_clusts("mmdysf_cc", filt_clusts=c(8,9,14), new_id = "mmdysf_cc_filt", reverse=T)
  mcell_gset_remove_clusts("mmdysf_ifn", filt_clusts=c(1:15), new_id = "mmdysf_ifn_filt", reverse=T)
  mcell_gset_remove_clusts("mmdysf_stress", filt_clusts=c(1:15), new_id = "mmdysf_stress_filt", reverse=T)
  
  lat_gsets = lapply(c('mmdysf_cc_filt','mmdysf_ifn_filt', 'mmdysf_stress_filt'), scdb_gset) #'mmdysf_cc_filt', 
  new_lat_genes = unlist(lapply(1:length(lat_gsets), function(i) { gs = lat_gsets[[i]]@gene_set; v = rep(i, length(gs)); names(v) = names(gs); return(v) }))
  
  mat_f = scdb_mat(mm_filt_id)
  non_filt_gset = scdb_gset(mm_filt_id)
  
  # cc_genes = names(non_filt_gset@gene_set)[non_filt_gset@gene_set %in% c(10,18,25,51,74,92)]#c(17,43,46,48,21,12,5,20)] #c(5,12,17,20,21)
  # stress_genes = names(non_filt_gset@gene_set)[non_filt_gset@gene_set %in% c(75,76,78,93)]#c(34,41,42,47,44,49,3)] #c(3,)
  # interf_genes = names(non_filt_gset@gene_set)[non_filt_gset@gene_set %in% c(14,20,23,28,31,61,83,84,85,86,87,88,89,90,94,95,99)]
  # 
  # cc_genes = names(non_filt_gset@gene_set)[non_filt_gset@gene_set %in% c(43,53,82,87,94)]#c(17,43,46,48,21,12,5,20)] #c(5,12,17,20,21)
  # stress_genes = names(non_filt_gset@gene_set)[non_filt_gset@gene_set %in% c(2,9,35,81,95)]#c(34,41,42,47,44,49,3)] #c(3,)
  # interf_genes = names(non_filt_gset@gene_set)[non_filt_gset@gene_set %in% c(7,12,15,20,24,27,40,46,51,58,71,73,75,77,85,86,88,89,90,91,98:100)]
  
  cc_genes = names(non_filt_gset@gene_set)[non_filt_gset@gene_set %in% c(11,35,58,60,110,130,131)]#c(17,43,46,48,21,12,5,20)] #c(5,12,17,20,21)
  stress_genes = names(non_filt_gset@gene_set)[non_filt_gset@gene_set %in% c(72,77,80,119,124)]#c(34,41,42,47,44,49,3)] #c(3,)
  interf_genes = names(non_filt_gset@gene_set)[non_filt_gset@gene_set %in% c(7,68,75,114,127,135,136)]
  
  high_mean_genes = c(grep('^Rp[ls]', mat_f@genes, v=T, perl=T),"Hsp90ab1", "Hspa1a", "Hspa5", "Ybx1, Ldha", "Tuba1c", "Eif5a", "Eif5b", "Npm1", "Rac2", "AC153498.1", "Tubb4b", "H2afx", "Ube2c", "Birc5", "AC153546", "H2afv", "Tmpo", "Hdgf", "Uhrf1", "Cdk1", "Racgap1", "AW112010", "Tagln2", "Hnrnpa2b", "Anp32b", "Ybx1", "Txn1", "AC149090.1","Isg10", "Ly6c1","Ly6c2", names(non_filt_gset@gene_set)[non_filt_gset@gene_set %in% c(2,8,16,21,64,74,78,82,120,144)]) #, "Cd52", "Npm1", "Vim", "Rac2", "Lgals1", "Gm11478", "S100a10", "AC153498.1", "Gm20390", "Lsp1", "Crip1", "Actb", "Tmsb4x"
  rm_genes = unique(c(cc_genes,stress_genes,interf_genes,high_mean_genes))
  add_cc_genes = grep('^Hist|^Cenp|^Smc[0-9]', mat_f@genes, v=T, perl=T)
  add_cc_genes = add_cc_genes[!add_cc_genes %in% rm_genes]
  add_ifn_genes = grep('^Ifi', mat_f@genes, v=T, perl=T)
  add_ifn_genes = add_ifn_genes[!add_ifn_genes %in% rm_genes]
  non_chosen_lat_genes = rm_genes[!rm_genes %in% names(non_filt_gset@gene_set)]
  add_genes = rep(51:53, times=c(length(add_cc_genes), length(add_ifn_genes), length(non_chosen_lat_genes)))
  names(add_genes) = c(add_cc_genes, add_ifn_genes, non_chosen_lat_genes)
  scdb_add_gset(mm_lateral_gset_id, gset_new_gset(c(non_filt_gset@gene_set[names(non_filt_gset@gene_set) %in% rm_genes], add_genes), 'lateral: CC, IFN, stress'))
  lat_gset = scdb_gset(mm_lateral_gset_id)
  # non_filt_gset@gene_set = non_filt_gset@gene_set[! names(non_filt_gset@gene_set) %in% rm_genes]
  # scdb_add_gset(id = paste(mm_filt_id, "_filt", sep=""), new_gset)
  
  # manually add some genes to make sure they're in
  
  add_cc_genes = grep('^Hist|^Cenp|^Smc[0-9]', mat_f@genes, v=T, perl=T)
  add_ifn_genes = grep('^Ifi', mat_f@genes, v=T, perl=T)
  add_genes = rep(1:2, times=c(length(add_cc_genes), length(add_ifn_genes)))
  names(add_genes) = c(add_cc_genes, add_ifn_genes)
  add_genes = add_genes[setdiff(names(add_genes), names(new_lat_genes))]
  scdb_add_gset(mm_lateral_gset_id, gset_new_gset(c(new_lat_genes, add_genes), 'lateral: CC, IFN, stress'))
  
  
}

mmdysf_data_gen = function(rebuild_lat_gsets=F, rebuild_mat=F) 
{
  if (rebuild_mat) {
    mmdysf_build_mat()
  }
  
  if (rebuild_lat_gsets) {
    mmdysf_build_lat_gsets()
  }
  
  mat = scdb_mat(mm_filt_id)
  md = mat@cell_metadata[colnames(mat@mat),]
  
  ### All B16-OVA and MC38##
  ctype = c('spec', 'byst', 'endo', 'cd8_in_vitro', 'pd1_ko_spec', 'ctrl_ko_spec', 'spec_cd25+', 'spec_cd25-', 'spec_tim3+pd1+')
  ttype = c('tumor', 'spleen', 'LN','cd8_in_vitro' )
  pd1_bsets = c("PD1_2", "PD1_6", "Ctrl_2", "Ctrl_6", 'PD1_5', 'Ctrl_5', 'PD1_treatment', 'PD1_ctrl', 'naive_transfer', 'cd8_in_vitro', 'InVitro_activated', 'InVitro_activated2','41bb_8', 'pd1_8', 'Ctrl_8', '41bb+pd1_8', '41bb_9', 'pd1_9', 'Ctrl_9', '41bb+pd1_9','S6_Panel_PD1_1', 'S6_Panel_Ctrl_1', 'PD1_7', 'MC38_41BB+PD1_2', 'MC38_41BB_2', 'MC38_41BB+PD1_3', 'MC38_Ctrl_2', 'MC38_PD1_3', 'MC38_PD1_2', 'MC38_41BB_3', 'MC38_Ctrl_3', 'BFP_pd1', '41BB_ko_pd1', 'Delay_OT1_GFPOT1_PD1', 'Delay_only_OT1_ctrl', 'Delay_only_OT1_PD1', 'Delay_OT1_GFPOT1_ctrl', 'PD1_ko_pd1', 'PD1_ko_ctrl', 'Zbdb32_ko_pd1', 'Zbdb32_ko_Ctrl', 'BFP_Ctrl', '41BB_ko_Ctrl') #, "Delay_only_OT1_PD1", "Delay_only_OT1_ctrl","Delay_OT1_GFPOT1_PD1", "Delay_OT1_GFPOT1_ctrl" 'ID3_control', 'ID3_KO',
  
  cells_all = mat@cell_metadata[mat@cells, ] %>%
    tibble::rownames_to_column("cell_id") %>%
    filter(batch_set_id_saved %in% pd1_bsets & cell_type %in% ctype & location %in% ttype)
  
  
  #cells_all = rbind(cells_all, cells)
  
  # cells_all = cells_all %>%
  #   filter(!batch_set_id_saved %in% c("InVitro_activated2"))
  
  # ctype = c('spec')
  # ttype = 'tumor'
  # pd1_bsets = c('InVitro_activated')
  # 
  # cells = mat@cell_metadata[mat@cells, ] %>%
  #   tibble::rownames_to_column("cell_id") %>%
  #   filter(batch_set_id_saved %in% pd1_bsets & cell_type %in% ctype & location %in% ttype & days_post_transfer %in% c("29"))
  
  # cells_all = rbind(cells_all, cells)
  
  cells_all = cells_all %>%
    filter(!mouse_id %in% c("m76", "m83", "m82", "m81", "m80", "m89", "m88", "m87", "m1", "m2", "m3", "m4", "m36", "m37", "m38", "m39", "m51", "m198", "m195", "m177", "m216"))
  
  # # First time this should be false, after first run of metacells choose cells to filter
  filter_tum_dend = F
  if (filter_tum_dend){
    filt_cells_csv =  read.csv(file=paste("figs", "/filt_cells.csv", sep=""), stringsAsFactors = F)
    filt_cells = filt_cells_csv$x
  
    
    filt_mcs = (1:max(mc@mc))[mc@colors %in% c("brown", "red")]
    #filt_cells = c()
    # filt_cells_prev = filt_cells
    filt_cells = unique(c(filt_cells, names(mc@mc[mc@mc %in% filt_mcs])))
    cells_all = cells_all[! cells_all$cell_id %in% filt_cells,]
    cells_all = cells_all[cells_all$cell_id %in% names(mc@mc),]
  }
  #write.csv(x = filt_cells, file=paste(scfigs_dir, "/filt_cells.csv", sep=""))
  #scdb_add_mat("mmdysf_filt_save", mat)
  # mat = scdb_mat("mmdysf_filt_save")
  # scdb_add_mat(mm_filt_id, mat)
  # 
  cells_filt = colnames(mat@mat)[!colnames(mat@mat) %in% cells_all$cell_id]
  mcell_mat_ignore_cells(mm_filt_id, mm_filt_id, union(mat@ignore_cells, cells_filt))
  mat = scdb_mat(mm_filt_id)
  zero_genes = rownames(mat@mat)[rowSums(mat@mat) == 0]
  mcell_mat_ignore_genes(mm_filt_id, mm_filt_id, union(mat@ignore_genes, zero_genes))
  mat = scdb_mat(mm_filt_id)
  
  mcell_add_gene_stat(mm_filt_id, mm_filt_id, force=T)
  mcell_gset_filter_varmean(gset_id=mm_filt_id, gstat_id=mm_filt_id, T_vm=0.08, force_new=T)
  mcell_gset_filter_cov(gset_id = mm_filt_id, gstat_id=mm_filt_id, T_tot=50, T_top3=2)
  mcell_plot_gstats(gstat_id=mm_filt_id, gset_id=mm_filt_id)
  mcell_gset_split_by_dsmat(mm_filt_id, mm_filt_id, 150)
  mcell_plot_gset_cor_mats(mm_filt_id, mm_filt_id)
  
  dend_genes = names(non_filt_gset@gene_set)[non_filt_gset@gene_set %in% c(31,97, 89)]
  tum_genes = names(non_filt_gset@gene_set)[non_filt_gset@gene_set %in% c(50)]
  mon_genes = names(non_filt_gset@gene_set)[non_filt_gset@gene_set %in% c(46)]
  ery_genes = names(non_filt_gset@gene_set)[non_filt_gset@gene_set %in% c(37)]
  prf_genes = names(non_filt_gset@gene_set)[non_filt_gset@gene_set %in% c(15)]
  cd4_genes = names(non_filt_gset@gene_set)[non_filt_gset@gene_set %in% c(99)]
  
  cell_depths = colSums(mat@mat)
  dend_score = colSums(mat@mat[dend_genes,]) / cell_depths
  tum_score = colSums(mat@mat[tum_genes,]) / cell_depths
  mon_score = colSums(mat@mat[mon_genes,]) / cell_depths
  ery_score = colSums(mat@mat[ery_genes,]) / cell_depths
  prf_score = colSums(mat@mat[prf_genes,]) / cell_depths
  cd4_score = colSums(mat@mat[cd4_genes,]) / cell_depths
  
  plot(sort(prf_score))
  plot(sort(cd4_score))
  cd4_cells = colnames(mat@mat)[cd4_score > 0.005]
  plot(sort(dend_score))
  dend_cells = colnames(mat@mat)[dend_score > 0.01]
  plot(sort(tum_score))
  tum_cells = colnames(mat@mat)[tum_score > 0.01]
  plot(sort(mon_score))
  mon_cells = colnames(mat@mat)[mon_score > 0.01]
  plot(sort(ery_score))
  ery_cells = colnames(mat@mat)[ery_score > 0.01]
  non_t_cells= unique(c(dend_cells, tum_cells, mon_cells, ery_cells, cd4_cells))
  
  cells_all = cells_all[! cells_all$cell_id %in% non_t_cells,]
  
  mel_basic_mat2mc(mm_filt_id, pd1_bsets, mm_lateral_gset_id,
                   T_top3 = 2, T_lfc=3, cgraph_knn=120, min_mc_size = 30, mc_K = 30, T_vm=0.08, T_tot=50, name=paste("combined_atlas", sep="_", collapse = "_"), cells=cells_all$cell_id) #"spleen_LN"
   rl()
  
  
  #### Picseq cd11c cd8 ####
  ctype ='endo'
  ttype = 'tumor'
  pd1_bsets = c('B16_OVA_Picseq_aPD1_1', 'B16_OVA_Picseq_ctrl_1')
  pic_descrip = c('Tumor5_aPD1_tPIC_1', 'Tumor5_aPD1_tPIC_2' ,'Tumor1_aPD1_tPIC_1', 'Tumor1_aPD1_tPIC_2', 'Tumor3_aPD1_tPIC_1', 'Tumor3_aPD1_tPIC_2', 'Tumor8_ctrl_tPIC_1', 'Tumor8_ctrl_tPIC_2', 'Tumor9_ctrl_tPIC_1')
  cd8_descrip = c('Tumor5_aPD1_CD8+_1', 'Tumor5_aPD1_CD8+_2', 'Tumor1_aPD1_CD8+_1', 'Tumor1_aPD1_CD8+_2', 'Tumor3_aPD1_CD8+_1', 'Tumor3_aPD1_CD8+_2', 'Tumor8_ctrl_CD8+_1', 'Tumor8_ctrl_CD8+_2', 'Tumor9_ctrl_CD8+_1')
  cd11c_descrip = c('Tumor5_aPD1_CD11c_1', 'Tumor5_aPD1_CD11c_2' ,'Tumor1_aPD1_CD11c_1', 'Tumor1_aPD1_CD11c_2', 'Tumor3_aPD1_CD11c_1', 'Tumor3_aPD1_CD11c_2', 'Tumor8_ctrl_CD11c_1', 'Tumor8_ctrl_CD11c_2', 'Tumor9_ctrl_CD11c_1', 'Tumor9_ctrl_CD11c_2')
  
  cells = mat@cell_metadata[mat@cells, ] %>%
    tibble::rownames_to_column("cell_id") %>%
    filter(batch_set_id %in% pd1_bsets & cell_type == ctype & location %in% ttype & Description %in% c(cd8_descrip, cd11c_descrip))
  
  mel_basic_mat2mc(mm_filt_id, pd1_bsets, mm_lateral_gset_id,
                   T_top3 = 2, T_lfc=3, cgraph_knn=80, min_mc_size = 20, mc_K = 40, T_vm=0.08, T_tot=50, name=paste("pic_cd8_dend", sep="_", collapse = "_"), cells=cells$cell_id) #"spleen_LN"
  rl()
  
  cells = mat@cell_metadata[mat@cells, ] %>%
    tibble::rownames_to_column("cell_id") %>%
    filter(batch_set_id %in% pd1_bsets & cell_type == ctype & location %in% ttype & Description %in% pic_descrip)
  
  mel_basic_mat2mc(mm_filt_id, pd1_bsets, mm_lateral_gset_id,
                   T_top3 = 2, T_lfc=3, cgraph_knn=80, min_mc_size = 20, mc_K = 40, T_vm=0.08, T_tot=50, name=paste("pic_cd8_dend", "pics", sep="_", collapse = "_"), cells=cells$cell_id) #"spleen_LN"
  rl()
  
  #### Picseq bi-specific ####
  ctype = 'endo'
  ttype = c('tumor', 'LN')
  pd1_bsets = c('B16_Bispecific_XCR1_1', 'B16_Bispecific_PBS_1', 'B16_Bispecific_Syn_1')

  
  cells = mat@cell_metadata[mat@cells, ] %>%
    tibble::rownames_to_column("cell_id") %>%
    filter(batch_set_id %in% pd1_bsets & cell_type == ctype & location %in% ttype)
  
  mel_basic_mat2mc(mm_filt_id, pd1_bsets, mm_lateral_gset_id,
                   T_top3 = 2, T_lfc=3, cgraph_knn=80, min_mc_size = 20, mc_K = 30, T_vm=0.08, T_tot=50, name=paste("sing_bispecific_xcr1", sep="_", collapse = "_"), cells=cells$cell_id) #"spleen_LN"
  rl()
  

  ctype = 'endo_pic'
  cells = mat@cell_metadata[mat@cells, ] %>%
    tibble::rownames_to_column("cell_id") %>%
    filter(batch_set_id %in% pd1_bsets & cell_type == ctype & location %in% ttype)
  
  mel_basic_mat2mc(mm_filt_id, pd1_bsets, mm_lateral_gset_id,
                   T_top3 = 2, T_lfc=3, cgraph_knn=80, min_mc_size = 20, mc_K = 40, T_vm=0.08, T_tot=50, name=paste("pic_bispecific_xcr1", sep="_", collapse = "_"), cells=cells$cell_id) #"spleen_LN"
  rl()
  
  ######## Non combined models #########
  
  ctype ='spec'
  ttype = 'tumor'
  pd1_bsets = c('PD1_5', 'Ctrl_5')
  
  cells = mat@cell_metadata[mat@cells, ] %>%
    tibble::rownames_to_column("cell_id") %>%
    filter(batch_set_id %in% pd1_bsets & cell_type == ctype & location %in% ttype & (!mouse_id %in% c("m76", "m83", "m82", "m81", "m80", "m89", "m88", "m87")))
  
  mel_basic_mat2mc(mm_filt_id, pd1_bsets, mm_lateral_gset_id,
                   T_top3 = 2, T_lfc=3, cgraph_knn=120, min_mc_size = 40, mc_K = 40, T_vm=0.10, T_tot=50, name=paste(ctype, ttype, sep="_", collapse = "_"), cells=cells$cell_id) #"spleen_LN"
  rl()
  
  ctype =c('spec', 'cd8_in_vitro')
  ttype = c('spleen', 'LN', 'cd8_in_vitro')
  pd1_bsets = c('PD1_5', 'Ctrl_5', 'cd8_in_vitro')
  
  cells = mat@cell_metadata[mat@cells, ] %>%
    tibble::rownames_to_column("cell_id") %>%
    filter(batch_set_id %in% pd1_bsets & cell_type %in% ctype & location %in% ttype & (!mouse_id %in% c("m132","m134", "m140", "m144", "m137", "m145")))
  
  mel_basic_mat2mc(mm_filt_id, pd1_bsets, mm_lateral_gset_id,
                   T_top3 = 2, T_lfc=3, cgraph_knn=100, min_mc_size = 40, mc_K = 40, T_vm=0.09, T_tot=50, name=paste('spec', 'spleen_LN', sep="_", collapse = "_"), cells=cells$cell_id) #"spleen_LN"
  rl()
  
  ctype ='byst'
  ttype = c('tumor','spleen', 'LN')
  pd1_bsets = c('PD1_5', 'Ctrl_5')
  
  cells = mat@cell_metadata[mat@cells, ] %>%
    tibble::rownames_to_column("cell_id") %>%
    filter(batch_set_id %in% pd1_bsets & cell_type == ctype & location %in% ttype)
  
  mel_basic_mat2mc(mm_filt_id, pd1_bsets, mm_lateral_gset_id,
                   T_top3 = 2, T_lfc=3, cgraph_knn=80, min_mc_size = 40, mc_K = 40, T_vm=0.09, T_tot=50, name=paste(ctype, 'spleen', 'LN', sep="_", collapse = "_"), cells=cells$cell_id) #"spleen_LN"
  rl()
  
  ctype ='endo'
  ttype = c('tumor', 'spleen')
  pd1_bsets = c('PD1_treatment', 'PD1_ctrl', 'naive_transfer')
  
  cells = mat@cell_metadata[mat@cells, ] %>%
    tibble::rownames_to_column("cell_id") %>%
    filter(batch_set_id %in% pd1_bsets & cell_type == ctype & location %in% ttype)
  
  mel_basic_mat2mc(mm_filt_id, pd1_bsets, mm_lateral_gset_id,
                   T_top3 = 2, T_lfc=3, cgraph_knn=100, min_mc_size = 40, mc_K = 40, T_vm=0.1, T_tot=50, name=paste(ctype, 'tumor', 'spleen', sep="_", collapse = "_"), cells=cells$cell_id) #"spleen_LN"
  rl()
  
  ctype ='endo'
  ttype = c('tumor')
  pd1_bsets = c('MC38_Ctrl_2', 'MC38_PD1_2', 'MC38_41BB_2', 'MC38_41BB+PD1_2')
  
  cells = mat@cell_metadata[mat@cells, ] %>%
    tibble::rownames_to_column("cell_id") %>%
    filter(batch_set_id %in% pd1_bsets & cell_type == ctype & location %in% ttype)
  filter_tum_dend = TRUE
  if (filter_tum_dend){
    filt_mcs = c(1:13)#(1:max(mc@mc))[mc@colors=="brown"]
    filt_cells = c()
    # filt_cells_prev = filt_cells
    filt_cells = unique(c(filt_cells, names(mc@mc[mc@mc %in% filt_mcs])))
    
    cells = cells[! cells$cell_id %in% filt_cells,]
  }
  mel_basic_mat2mc(mm_filt_id, pd1_bsets, mm_lateral_gset_id,
                   T_top3 = 2, T_lfc=3, cgraph_knn=80, min_mc_size = 40, mc_K = 40, T_vm=0.08, T_tot=50, name=paste("MC38", ctype, ttype, sep="_", collapse = "_"), cells=cells$cell_id) #"spleen_LN"
  rl()
  
  
  # for (ttype in c('spleen', "LN")){
  #   cells = mat@cell_metadata[mat@cells, ] %>% 
  #     tibble::rownames_to_column("cell_id") %>% 
  #     filter(batch_set_id %in% pd1_bsets & location == ttype)  
  #   
  #   mel_basic_mat2mc(mm_filt_id, pd1_bsets, mm_lateral_gset_id, 
  #                    T_top3 = 2, T_lfc=3, cgraph_knn=70, min_mc_size = 30, mc_K = 30, T_vm=0.09, T_tot=50, name=paste('spec_or_byst', ttype, sep="_"), cells=cells$cell_id)
  # }
  
  ######## Projected Models #########
  # S6 panel 1
  ctype =c('spec', 'endo')
  ttype = 'tumor'
  pd1_bsets = c('S6_Panel_PD1_1', 'S6_Panel_Ctrl_1')
  
  cells = mat@cell_metadata[mat@cells, ] %>%
    tibble::rownames_to_column("cell_id") %>%
    filter(batch_set_id %in% pd1_bsets & cell_type %in% ctype & location %in% ttype)
  
  mel_basic_mat2mc(mm_filt_id, pd1_bsets, mm_lateral_gset_id,
                   T_top3 = 2, T_lfc=3, cgraph_knn=80, min_mc_size = 40, mc_K = 40, T_vm=0.08, T_tot=50, name=paste("s6_1",sep="_", collapse = "_"), cells=cells$cell_id) ##[!cells$cell_id %in% dend_cells]
  rl()
  # B16-OVA KO Zbdb32 41bb
  ctype ='spec'
  ttype = 'tumor'
  pd1_bsets = c('Zbdb32_ko_Ctrl', 'Zbdb32_ko_pd1', 'BFP_Ctrl', 'BFP_pd1')
  
  cells = mat@cell_metadata[mat@cells, ] %>%
    tibble::rownames_to_column("cell_id") %>%
    filter(batch_set_id %in% pd1_bsets & cell_type == ctype & location %in% ttype)
  
  mel_basic_mat2mc(mm_filt_id, pd1_bsets, mm_lateral_gset_id,
                   T_top3 = 2, T_lfc=3, cgraph_knn=80, min_mc_size = 40, mc_K = 40, T_vm=0.08, T_tot=50, name=paste(ctype, ttype, "zbdb32_41bb_ko_filt",sep="_", collapse = "_"), cells=cells$cell_id) ##[!cells$cell_id %in% dend_cells]
  rl()
  # B16-OVA 41BB + PD1
  ctype ='spec'
  ttype = 'tumor'
  pd1_bsets = c('41bb+pd1_8', '41bb_8', 'pd1_8', 'Ctrl_8')
  
  cells = mat@cell_metadata[mat@cells, ] %>%
    tibble::rownames_to_column("cell_id") %>%
    filter(batch_set_id %in% pd1_bsets & cell_type == ctype & location %in% ttype)
  
  mel_basic_mat2mc(mm_filt_id, pd1_bsets, mm_lateral_gset_id,
                   T_top3 = 2, T_lfc=3, cgraph_knn=80, min_mc_size = 40, mc_K = 40, T_vm=0.08, T_tot=50, name=paste(ctype, ttype, "41bb_pd1",sep="_", collapse = "_"), cells=cells$cell_id) #"spleen_LN"
  rl()
  
  
  i = 13
  id = mm_endo_spec_ids[i]
  gset = scdb_gset(id)
  
  # MC38 Models total cd3+
  ctype ='endo'
  ttype = 'tumor'
  pd1_bsets = c('MC38_Ctrl_1')
  
  cells = mat@cell_metadata[mat@cells, ] %>%
    tibble::rownames_to_column("cell_id") %>%
    filter(batch_set_id %in% pd1_bsets & cell_type == ctype & location %in% ttype)
  
  mel_basic_mat2mc(mm_filt_id, pd1_bsets, mm_lateral_gset_id,
                   T_top3 = 2, T_lfc=3, cgraph_knn=80, min_mc_size = 40, mc_K = 40, T_vm=0.08, T_tot=50, name=paste(ctype, pd1_bsets, sep="_", collapse = "_"), cells=cells$cell_id) #"spleen_LN"
  rl()
  
  
  mat = scdb_mat("mmdysf_filt_save")
  ctype = c('spec')
  ttype = c('tumor')
  pd1_bsets = c('ID3_control','ID3_KO')
  
  cells = mat@cell_metadata[mat@cells, ] %>%
    tibble::rownames_to_column("cell_id") %>%
    filter(batch_set_id_saved %in% pd1_bsets)
  
  mel_basic_mat2mc(mm_filt_id, pd1_bsets, mm_lateral_gset_id,
                   T_top3 = 2, T_lfc=3, cgraph_knn=80, min_mc_size = 40, mc_K = 40, T_vm=0.08, T_tot=50, name=paste("ID3_KO", sep="_", collapse = "_"), cells=cells$cell_id, ref_gset_id = id) #"spleen_LN"
  rl()
  
  ttype = c('tumor')
  pd1_bsets = c("Delay_only_OT1_PD1", "Delay_only_OT1_ctrl","Delay_OT1_GFPOT1_PD1", "Delay_OT1_GFPOT1_ctrl")
  
  cells = mat@cell_metadata[mat@cells, ] %>%
    tibble::rownames_to_column("cell_id") %>%
    filter(batch_set_id_saved %in% pd1_bsets & location %in% ttype)
  
  mel_basic_mat2mc(mm_filt_id, pd1_bsets, mm_lateral_gset_id,
                   T_top3 = 2, T_lfc=3, cgraph_knn=80, min_mc_size = 40, mc_K = 40, T_vm=0.08, T_tot=50, name=paste("delay_spec", sep="_", collapse = "_"), cells=cells$cell_id, ref_gset_id = id) #"spleen_LN"
  rl()
  
  ctype = c('pd1_ko_spec')
  ttype = c('tumor')
  pd1_bsets = c('PD1_ko_pd1','PD1_ko_ctrl')
  
  cells = mat@cell_metadata[mat@cells, ] %>%
    tibble::rownames_to_column("cell_id") %>%
    filter(batch_set_id_saved %in% pd1_bsets & cell_type %in% ctype & location %in% ttype)
  
  mel_basic_mat2mc(mm_filt_id, pd1_bsets, mm_lateral_gset_id,
                   T_top3 = 2, T_lfc=3, cgraph_knn=80, min_mc_size = 40, mc_K = 40, T_vm=0.08, T_tot=50, name=paste("pd1_crispr_ko", sep="_", collapse = "_"), cells=cells$cell_id, ref_gset_id = id) #"spleen_LN"
  rl()
  
  ctype = c('ctrl_ko_spec')
  ttype = c('tumor')
  pd1_bsets = c('PD1_ko_pd1','PD1_ko_ctrl')
  
  cells = mat@cell_metadata[mat@cells, ] %>%
    tibble::rownames_to_column("cell_id") %>%
    filter(batch_set_id_saved %in% pd1_bsets & cell_type %in% ctype & location %in% ttype)
  
  mel_basic_mat2mc(mm_filt_id, pd1_bsets, mm_lateral_gset_id,
                   T_top3 = 2, T_lfc=3, cgraph_knn=80, min_mc_size = 40, mc_K = 40, T_vm=0.08, T_tot=50, name=paste("pd1_crispr_ctrl", sep="_", collapse = "_"), cells=cells$cell_id, ref_gset_id = id) #"spleen_LN"
  rl()
  
  ######## Annotation #########
  confs = list()
  for (i in c(5, 7:12)) { #seq_along(mm_mc_endo_spec_ids)
    id = mm_endo_spec_ids[i]
    mc_id = mm_mc_endo_spec_ids[i]
    message(id)
    
    #12 "#cee6b9" effector-memory
    # first time - run without coloring, then color by supmc and marks files
    confs[[id]] = colorize_by_confusion_mat(mc_id=mc_id, graph_id=id, show_mc_ids=T)
    # 188   "darkblue"  naive-dysf
    confs[[id]] = colorize_by_confusion_mat(mc_id=mc_id, graph_id=id, res=confs[[id]], show_mc_ids=T, supmc_file=sprintf("config/%s_supmc.txt", id), marks_file=sprintf("config/%s_marks.txt", mm_endo_spec_ids[1]))
    
  }
  
  for (i in c(7:10, 12)){
    mf1 = c('mouse_id', 'location', 'cell_type', 'batch_set_id', 'days_post_transfer')
    id = mm_endo_spec_ids[i]
    mc_id = mm_mc_endo_spec_ids[i]
    message(id)
    mc = scdb_mc(mc_id)
    mc@colors = rep("white", max(mc@mc))
    scdb_add_mc(mc_id, mc)
    mel_basic_mc_mc2d_plots(mc_id, mc_id, id, id, mm_lateral_gset_id, metadata_fields_to_export=mf1)
  }
  
  loc_mc_ids = c(7:10, 12, 13, 14, 15, 17)
  atlas_mc_ids = c(1,1,1,1, 1, 1, 1, 1, 1)
  atlas_ids = c("mmdysf_filt_combined_atlas")
  for (ind in 1:length(loc_mc_ids)){
    print(ind)
    print(atlas_ids[atlas_mc_ids[ind]])
    atlas_id = atlas_ids[atlas_mc_ids[ind]]
    atlas_mc_id = paste(atlas_id, "_outClean", sep="")
    query_id = mm_mc_endo_spec_ids[loc_mc_ids[ind]]
    atlas = mcell_gen_atlas(mat_id=id, mc_id=atlas_mc_id, gset_id=atlas_mc_id, mc2d_id=atlas_mc_id, atlas_cols = NULL)
    proj_on_atlas = mcell_proj_on_atlas(mat_id="mmdysf_filt", mc_id=query_id, atlas=atlas, 
                                        fig_cmp_dir= paste("figs/", query_id, '.atlas_proj',sep=""), 
                                        ten2mars=F,
                                        gene_name_map=NULL,
                                        recolor_mc_id = paste(query_id, sep=""), #paste(query_id, "_proj", sep="")
                                        plot_all_mcs = T,
                                        md_field='batch_set_id',
                                        max_entropy=4,
                                        burn_cor=0.8)
  }
  
  
  i = 5
  id = mm_endo_spec_ids[i]
  mc_id = mm_mc_endo_spec_ids[i]
  message(id)
  
  # first time - run without coloring, then color by supmc and marks files
  #conf = colorize_by_confusion_mat(mc_id=mc_id, graph_id=id, show_mc_ids=T)
  sup = confs[[id]]$mc_sup
  length(sup)
  mc = scdb_mc(mc_id)
  col2group = get_mc_col2group(mc)
  #mc@colors[c(12,13)] = "mediumorchid4"
  #mc@colors[c(9)] = "gold"
  #mc@colors[c(6)] = "darkgreen"
  # mc@colors[122] = "#5fb8f4"
  # scdb_add_mc(mc_id, mc)
  lfp = log2(mc@mc_fp)
  lfp_cor = tgs_cor(t(lfp))
  
  # 10 74 204 170 (313) 285 333 442 555 (526) 710 953 1316
  ## MC38 Model annotation
  plt('Cd8a', 'Cd4', lfp, mc@colors)
  plt('Nrgn', 'Tcf7', lfp, mc@colors)
  plt('Mki67', 'Pcna', lfp, mc@colors)
  plt('Nrgn', 'Cd74', lfp, mc@colors)
  plt('Mlana', 'Cd74', lfp, mc@colors)
  plt('Apoe', 'Cd74', lfp, mc@colors)
  plt('Cd8a', 'Cd4', lfp, ifelse(1:ncol(lfp) %in% sup[[56]]$mcs, 'blue', 'white'))
  plt('Tcf7', 'Klf2', lfp, ifelse(1:ncol(lfp) %in% sup[[17]]$mcs, 'blue', 'white'))
  plt('Eomes', 'Klf2', lfp, ifelse(1:ncol(lfp) %in% sup[[20]]$mcs, 'blue', 'white'))
  plt('Prf1', 'Gzmb', lfp, ifelse(1:ncol(lfp) %in% sup[[32]]$mcs, 'blue', 'white'))
  plt('Prf1', 'Havcr2', lfp, ifelse(1:ncol(lfp) %in% sup[[2]]$mcs, 'blue', 'white'))
  plt('Prf1', 'Pdcd1', lfp, ifelse(1:ncol(lfp) %in% sup[[1009]]$mcs, 'blue', 'white'))
  plt('Pdcd1', 'Havcr2', lfp, ifelse(1:ncol(lfp) %in% sup[[464]]$mcs, 'blue', 'white'))
  plt('Il7r', 'Ccl5', lfp, ifelse(1:ncol(lfp) %in% sup[[6]]$mcs, 'blue', 'white'))
  plt('Sell', 'Ccl5', lfp, ifelse(1:ncol(lfp) %in% sup[[14]]$mcs, 'blue', 'white'))
  
  #9V-26V-29V-38-52-64V-75-93V
  # supid color name
  # 4 "brown" tumor-na
  # 47  "#5fb8f4" naive1
  # 91  "#1f78b4" cycling_trans
  # 14  "#cee6b9" em1
  # 59  "#6a8255" em2
  # 31  "#975e9f" trans-dysf2
  # 38  "deeppink"  trans-dysf3
  # 2 "#e2abd6" trans-dysf1
  
  #Spec only
  # supid color name
  # 6 "brown" tumor-na
  # 36  "#1f78b4" cycling_trans
  # 54  "deeppink"  trans-dysf3
  # 8 "darkseagreen3" em3
  # 12  "#5fb8f4" naive1
  # 15  "#6a8255" em2
  # 67  "#975e9f" trans-dysf2
  # 2 "#e2abd6" trans-dysf1
  
  # 1dend 5mon_dend 17tumor 204central_mem 415_trans_eff 614trans_dysf 796resm_mem 890cytox_dys 944dys 1029dys_low 1159nai 1256em
  # 225gold 264nfkbid 340tum-dend 345em 424naive 584trans 854yellowcyc 870resmem 1043centmem
  plt('Tcf7', 'Klf2', lfp, ifelse(1:ncol(lfp) %in% sup[[32]]$mcs, 'blue', 'white'))
  plt('Pcna', 'Mki67', lfp, ifelse(1:ncol(lfp) %in% sup[[101]]$mcs, 'blue', 'white'))
  cor_lfp_genes = tgs_cor(t(lfp))
  tail(sort(cor_lfp_genes["Havcr2",]), 30)
  tail(sort(cor_lfp_genes["Klf2",]), 30)
  plt('Ctla4','Lag3', lfp, mc@colors)
  plt('Pdcd1','Havcr2', lfp, mc@colors)
  plt('Dgka', 'Klf2', lfp, mc@colors)
  plt('Eomes', 'Lag3', lfp, mc@colors)
  plt('Klrk1', 'Cd44', lfp, mc@colors)
  plt('Cx3cr1', 'Gzmb', lfp, mc@colors)
  plt('Gzmb', 'Cd69', lfp, mc@colors)
  plt('Il7r', 'Klf2', lfp, mc@colors)
  plt('Dapl1', 'Vim', lfp, mc@colors)
  plt('Lef1', 'Klf2', lfp, mc@colors)
  plt('Il7r', 'Tcf7', lfp, mc@colors)
  plt('Tigit', 'Eomes', lfp, mc@colors)
  plt('Nfkbid', 'Ccl4', lfp, mc@colors)
  plt('Il7r', 'Ccl5', lfp, mc@colors)
  plt('Tox', 'Slamf6', lfp, mc@colors)
  plt('Lag3', 'Tcf7', lfp, mc@colors)
  plt('Ccl4', 'Ccl3', lfp, mc@colors)
  plt('Gzmb', 'Prf1', lfp, mc@colors)
  plt('Sell', 'S1pr1', lfp, mc@colors)
  plt('Klf2', 'S1pr1', lfp, mc@colors)
  plt('Prf1', 'Havcr2', lfp, mc@colors)
  plt('Prf1', 'Pdcd1', lfp, mc@colors)
  plt('Pdcd1', 'Havcr2', lfp, mc@colors)
  plt('Nrgn', 'Tcf7', lfp, mc@colors)
  plt('Mki67', 'Pcna', lfp, mc@colors)
  plt('Nrgn', 'Cd74', lfp, mc@colors)
  plt('Mlana', 'Cd74', lfp, mc@colors)
  plt('Sell', 'Ccr7', lfp, ifelse(1:ncol(lfp) %in% sup[[353]]$mcs, 'blue', 'white'))
  plt('Tcf7', 'Klf2', lfp, ifelse(1:ncol(lfp) %in% sup[[64]]$mcs, 'blue', 'white'))
  plt('Eomes', 'Klf2', lfp, ifelse(1:ncol(lfp) %in% sup[[20]]$mcs, 'blue', 'white'))
  plt('Prf1', 'Gzmb', lfp, ifelse(1:ncol(lfp) %in% sup[[32]]$mcs, 'blue', 'white'))
  plt('Prf1', 'Havcr2', lfp, ifelse(1:ncol(lfp) %in% sup[[2]]$mcs, 'blue', 'white'))
  plt('Prf1', 'Pdcd1', lfp, ifelse(1:ncol(lfp) %in% sup[[1]]$mcs, 'blue', 'white'))
  plt('Pdcd1', 'Havcr2', lfp, ifelse(1:ncol(lfp) %in% sup[[464]]$mcs, 'blue', 'white'))
  plt('Il7r', 'Ccl5', lfp, ifelse(1:ncol(lfp) %in% sup[[6]]$mcs, 'blue', 'white'))
  plt('Sell', 'Ccl5', lfp, ifelse(1:ncol(lfp) %in% sup[[14]]$mcs, 'blue', 'white'))
  
  plt('Sell', 'Mlana', lfp, ifelse(1:ncol(lfp) %in% sup[[4]]$mcs, 'blue', 'white'))
  plt('Cd74', 'Vim', lfp, ifelse(1:ncol(lfp) %in% sup[[3]]$mcs, 'blue', 'white'))
  plt('Cd74', 'Hbb-bs', lfp, ifelse(1:ncol(lfp) %in% sup[[22]]$mcs, 'blue', 'white'))
  tail(sort(lfp[, 54]), 30)
  plt('Gzma', 'Ccl1', lfp, ifelse(1:ncol(lfp) %in% sup[[227]]$mcs, 'blue', 'white'))
  plt('Gzma', 'Il7r', lfp, ifelse(1:ncol(lfp) %in% sup[[43]]$mcs, 'blue', 'white'))
  plt('Il7r', 'Ccl5', lfp, ifelse(1:ncol(lfp) %in% sup[[73]]$mcs, 'blue', 'white'))
  plt('Havcr2', 'Pdcd1', lfp, ifelse(1:ncol(lfp) %in% sup[[2]]$mcs, 'blue', 'white'))
  plt('Ccl4', 'Ccl3', lfp, ifelse(1:ncol(lfp) %in% sup[[2]]$mcs, 'blue', 'white'))
  plt('Nfkbid', 'Ccl4', lfp, ifelse(1:ncol(lfp) %in% sup[[178]]$mcs, 'blue', 'white'))
  plt('Sell', 'S1pr1', lfp, ifelse(1:ncol(lfp) %in% sup[[29]]$mcs, 'blue', 'white'))
  plt('Nrgn', 'Tcf7', lfp, ifelse(1:ncol(lfp) %in% sup[[2]]$mcs, 'blue', 'white'))
  query_sup_by_mcs(sup, which(lfp['Prf1', ] > 0.8 & lfp['Pdcd1', ] > 0.8))
  query_sup_by_mcs(sup, which(lfp['Cd74', ] > 2 | lfp['Mlana', ] > 1))
  query_sup_by_mcs(sup, which(lfp['Dapl1', ] > 1.3))
  query_sup_by_mcs(sup, which(lfp['Gzmc', ] > 2 ))
  query_sup_by_mcs(sup, c(1729:1732, 616))
  query_sup_by_mcs(sup, as.integer(names(table(mc@mc[rownames(mat@cell_metadata[mat@cell_metadata$batch_set_id=="cd8_in_vitro",])][mc@colors[mc@mc]=="#1f78b4"]))))
  confs[[id]] = colorize_by_confusion_mat(mc_id=mc_id, graph_id=id, res=confs[[id]], show_mc_ids=T, supmc_file=sprintf("config/%s_supmc.txt", id), marks_file=sprintf("config/%s_marks.txt", mm_endo_spec_ids[1]))
  
  # supid color name
  # 2 "brown" tumor-na
  # 271   "#1f78b4" cycling-dysf
  # 202   "blue"  cycling-trans
  # 174   "yellow"  cc
  # 19  "#e2abd6" cytokine-dysf
  # 4 "deeppink"  nfkb-dysf
  # 150 "#975e9f" effector-trans
  # 129 "gold"  cytotoxic-effector
  # 81  "#6a8255" memory-like
  # 112 "#5fb8f4" naive-memory
  # 57  "#cee6b9" em1
  # confs = list()
  # for (i in seq_along(mm_mc_endo_spec_ids)) {
  #   id = mm_endo_spec_ids[i]
  #   mc_id = mm_mc_endo_spec_ids[i]
  #   message(id)
  #   
  #   #12 "#cee6b9" effector-memory
  #   # first time - run without coloring, then color by supmc and marks files
  #   confs[[id]] = colorize_by_confusion_mat(mc_id=mc_id, graph_id=id, show_mc_ids=T)
  #   confs[[id]] = colorize_by_confusion_mat(mc_id=mc_id, graph_id=id, res=confs[[id]], show_mc_ids=T, supmc_file=sprintf("config/%s_supmc.txt", id), marks_file=sprintf("config/%s_marks.txt", id))
  #   
  # }
  mc = scdb_mc(mc_id)
  mc2d = scdb_mc2d(mc_id)
  md = mat@cell_metadata[names(mc@mc),]
  
  ignore_edges = NULL
  mgraph = mc2d_comp_mgraph(mc_id, id, ignore_mismatch=F, symmetrize=F)
  scdb_add_mgraph(id = mc_id, mgraph = tgMCManifGraph(mgraph = mgraph,mc_id = mc_id))
  ignore_edges = mgraph[mgraph$mc1==1777& mgraph$mc2==1785, ]
  ignore_edges = rbind(ignore_edges, mgraph[mgraph$mc1==1451& mgraph$mc2==1152, ])
  ignore_edges = rbind(ignore_edges, mgraph[mgraph$mc1==99& mgraph$mc2==427, ])
  ignore_edges = rbind(ignore_edges, mgraph[mgraph$mc1==1771& mgraph$mc2==38, ])
  ignore_edges = rbind(ignore_edges, mgraph[mgraph$mc1==1468& mgraph$mc2==99, ])
  ignore_edges = rbind(ignore_edges, mgraph[mgraph$mc1==1765& mgraph$mc2==1480, ])
  ignore_edges = rbind(ignore_edges, mgraph[mgraph$mc1==708& mgraph$mc2==1028, ])
  ignore_edges = rbind(ignore_edges, mgraph[mgraph$mc1==690& mgraph$mc2==1774, ])
  ignore_edges = rbind(ignore_edges, mgraph[mgraph$mc1==1771& mgraph$mc2==1730, ])
  ignore_edges = rbind(ignore_edges, mgraph[mgraph$mc1==884& mgraph$mc2==1730, ])
  ignore_edges = rbind(ignore_edges, mgraph[mgraph$mc1==884& mgraph$mc2==1743, ])
  ignore_edges = rbind(ignore_edges, mgraph[mgraph$mc1==665& mgraph$mc2==634, ])
  ignore_edges = rbind(ignore_edges, mgraph[mgraph$mc1==655& mgraph$mc2==636, ])
  
  ignore_edges = rbind(ignore_edges, mgraph[mgraph$mc1==636& mgraph$mc2==606, ])
  ignore_edges = rbind(ignore_edges, mgraph[mgraph$mc1==634& mgraph$mc2==606, ])
  ignore_edges = rbind(ignore_edges, mgraph[mgraph$mc1==605& mgraph$mc2==623, ])
  ignore_edges = rbind(ignore_edges, mgraph[mgraph$mc1==176& mgraph$mc2==572, ])
  ignore_edges = rbind(ignore_edges, mgraph[mgraph$mc1==802& mgraph$mc2==1775, ])
  ignore_edges = rbind(ignore_edges, mgraph[mgraph$mc1==660& mgraph$mc2==183, ])
  ignore_edges = rbind(ignore_edges, mgraph[mgraph$mc1==660& mgraph$mc2==609, ])
  ignore_edges = rbind(ignore_edges, mgraph[mgraph$mc1==1741& mgraph$mc2==1779, ])
  ignore_edges = rbind(ignore_edges, mgraph[mgraph$mc1==1480& mgraph$mc2==1481, ])
  ignore_edges = rbind(ignore_edges, mgraph[mgraph$mc1==800& mgraph$mc2==888, ])
  ignore_edges = rbind(ignore_edges, mgraph[mgraph$mc1==802& mgraph$mc2==1775, ])
  ignore_edges = rbind(ignore_edges, mgraph[mgraph$mc1==817& mgraph$mc2==1780, ])
  ignore_edges = rbind(ignore_edges, mgraph[mgraph$mc1==1044& mgraph$mc2==1469, ])
  ignore_edges = rbind(ignore_edges, mgraph[mgraph$mc1==101& mgraph$mc2==1613, ])
  # ignore_edges = rbind(ignore_edges, mgraph[mgraph$mc1==990& mgraph$mc2==992, ])
  # ignore_edges = rbind(ignore_edges, mgraph[mgraph$mc1==787& mgraph$mc2==795, ])
  # ignore_edges = rbind(ignore_edges, mgraph[mgraph$mc1==787& mgraph$mc2==815, ])
  # ignore_edges = rbind(ignore_edges, mgraph[mgraph$mc1==936& mgraph$mc2==919, ])
  # ignore_edges = rbind(ignore_edges, mgraph[mgraph$mc1==784& mgraph$mc2==108, ])
  # ignore_edges = rbind(ignore_edges, mgraph[mgraph$mc1==577& mgraph$mc2==888, ])
  # ignore_edges = rbind(ignore_edges, mgraph[mgraph$mc1==795& mgraph$mc2==790, ])
  # ignore_edges = NULL
  # "#975e9f"
  rl()
  mcell_mc2d_force_knn(mc_id, mc_id, id, ignore_mismatch=T, ignore_edges = ignore_edges)
  my_mcell_mc2d_plot(mc2d_id=mc_id, plot_single_cells = F, plot_edges = F)
  # Spec comb changes
  # mc@colors[c(433,700,544,695,520,542,672,233,678,565,553,537,552,707,375,357,692,377,659)] = "#1f78b4"
  # mc@colors[c()] = "gold"
  # mc@colors[c(517,659)] = "cyan"
  # mc@colors[c(736)] = "blue"
  # mc@colors[c()] = "magenta"
  # mc@colors[c(311,356,351,281,109,744,734)] = "#6a8255"
  # mc@colors[c(263)] = "#5fb8f4"
  # #mc@colors[c(517)] = "#e2abd6"
  # # mc@colors[c(267,44,75,36,65)] = "#cee6b9"
  # mc@colors[c(560,349,289,657)] = "#975e9f"
  # mc@colors[c(282,257,694)] = "brown"  #Cd4, "Hba-a1 Mlana Cd74
  
  plt('Prf1', 'Pdcd1', lfp, mc@colors)
  plt('Prf1', 'Gzmb', lfp, mc@colors)
  
  ### Combined b16 mc38 changes
  mc@colors[c()] = "#1f78b4"
  mc@colors[c()] = "gold"
  mc@colors[c(1738,885,1480,855,882,733,1741,742,797,836,837,789,787,788,998)] = "cyan"
  mc@colors[c()] = "blue"
  mc@colors[c(949, 1025)] = "magenta"
  mc@colors[c(1449,101,1783,1767,1467)] = "#6a8255"
  mc@colors[c(125,124)] = "#5fb8f4"
  mc@colors[c(838,1600,1661,1693)] = "#e2abd6"
  mc@colors[c()] = "#cee6b9"
  mc@colors[c(1422,1485)] = "#975e9f"
  mc@colors[c()] = "brown"  #Cd4, "Hba-a1 Mlana Cd74
  
  
  ### Combined changes
  mc@colors[c(808,856,927,858)] = "#1f78b4"
  mc@colors[c(561,850,552,560,829,675,845,930,912,913)] = "gold"
  mc@colors[c(846)] = "cyan"
  mc@colors[c()] = "blue"
  mc@colors[c(597,543, 578,592,558,835,568,553,562,565,832,851)] = "magenta"
  mc@colors[c()] = "#6a8255"
  mc@colors[c()] = "#5fb8f4"
  mc@colors[c(905,686)] = "#e2abd6"
  mc@colors[c(191,730,723,402)] = "#cee6b9"
  mc@colors[c()] = "#975e9f"
  mc@colors[c()] = "brown"  #Cd4, "Hba-a1 Mlana Cd74
  # 
  # mcell_mc_plot_hierarchy(mc_id, id, mc_order=c((1:max(mc@mc))[mc@colors=="deeppink"], (1:max(mc@mc))[mc@colors=="#e2abd6"], (1:max(mc@mc))[mc@colors=="#6a8255"], (1:max(mc@mc))[mc@colors=="#5fb8f4"], (1:max(mc@mc))[mc@colors=="#cee6b9"], (1:max(mc@mc))[mc@colors=="magenta"], (1:max(mc@mc))[mc@colors=="cyan"], (1:max(mc@mc))[mc@colors=="gold"], (1:max(mc@mc))[mc@colors=="#975e9f"], (1:max(mc@mc))[mc@colors=="#1f78b4"], (1:max(mc@mc))[mc@colors=="blue"], (1:max(mc@mc))[mc@colors=="brown"]), sup_mcs = sup, width=4000, height = 8000, show_mc_ids = T,min_nmc = 2)
  # 
  # PD1 ko changes
  mc@colors[8] = "#6a8255"
  mc@colors[c(19)] = "#e2abd6"
  mc@colors[c(5,7)] = "#1f78b4"
  mc@colors[c(1,12,13)] = "blue"
  
  # PD1 ko ctrl changes
  mc@colors[12] = "gold"
  mc@colors[c(7)] = "#6a8255"
  mc@colors[c(10,9)] = "#e2abd6"
  
  # delay spec changes
  mc@colors[2] = "#975e9f"
  mc@colors[c(8)] = "gold"
  mc@colors[c(3)] = "cyan"
  mc@colors[c(20)] = "#6a8255"
  mc@colors[c(9)] = "#e2abd6"
  
  # MC38 cd8 changes
  mc@colors[c(5)] = "#5fb8f4"
  mc@colors[c(2:4)] = "brown"
  mc@colors[c(10,26,18,17)] = "#1f78b4"
  mc@colors[c(24)] = "blue"
  mc@colors[c(19,23)] = "cyan"
  
  
  
  
  # Spec tumor changes
  mc@colors[203] = "#1f78b4"
  mc@colors[58] = "#cee6b9"
  mc@colors[320] = "blue"
  mc@colors[172] = "#975e9f"
  mc@colors[305] = "#6a8255"
  mc@colors[167] = "#1f78b4"
  mc@colors[195] = "#1f78b4"
  mc@colors[124] = "#e2abd6"
  mc@colors[123] = "#e2abd6"
  mc@colors[c(205,202)] = "gold"
  
  mc2d@graph = mc2d@graph[c(-677,-901),]
  scdb_add_mc2d(mc_id, mc2d)
  mcell_mc2d_plot(mc2d_id=mc_id)
  mc@color_key = mc@color_key[-7,]
  
  # Spec spleen changes
  mc@colors[27] = "brown"
  
  # Spec LN changes
  mc@colors[30] = "deeppink"
  mc@colors[29] = "#1f78b4"
  
  # Spec spleen LN changes
  mc@colors[c(47)] = "brown"
  
  
  #Byst changes
  mc@colors[c(38)] = "brown"
  mc@colors[c(47)] = "#1f78b4"
  
  #cd8_in_vitro changes
  mc@colors[c(5,6)] = "#5fb8f4"
  
  #endo changes
  mc@colors[c(3)] = "brown"
  mc@colors[c(44)] = "#e2abd6"
  mc@colors[c(61)] = "magenta"
  
  # mc@colors[26] = "brown"
  #Spec changes
  mc@colors[c(9,74,76)] = "brown"
  mc@colors[c(38,43,94,30,56)] = "#1f78b4"
  # mc@colors[c(89, 118, 115, 93, 94, 61)] = "#5fb8f4"
  mc@colors[c(26)] = "#6a8255"
  mc@colors[c(76,92)] = "darkseagreen3"
  #mc@colors[c(80,84,85)] = "#e2abd6"#"#975e9f"
  # mc@colors[c(70,63,66,71,65,67,91,69,76,73,75,64,68,72)] = "#1f78b4"
  scdb_add_mc(mc_id, mc)
  confs[[id]] = colorize_by_confusion_mat(mc_id=mc_id, graph_id=id, res=confs[[id]], show_mc_ids=T)
  
  mel_knn_mat = mc2d_comp_weighted_mgraph(mc_id = mc_id, graph_id = id, ignore_mismatch=F)
  write.csv(mel_knn_mat, "figs/mel_knn_mat.csv")
  write.csv(mc@colors, "figs/mel_mc_colors.csv")
  
  
}


genes2plt = unique(c('Cd3d', 'Cd8a', 'Nkg7', 'Lag3', 'Tigit', 'Mki67', 'Casp3', 'Il7r', 'Gzma', 'Id2', 'Ccl5', 'Cxcr6', 'Ccr2', 'Gzmk', 'Sell', "Id3", 'Tcf7', 'Tox', 'Ccr7', 'Lef1', 'Dapl1', 'Ctla4',  'Ccr7',  'Ctla4', 'Nrgn','Top2a', 'Hspa5', 'Il10ra','Tnfaip3', 'Klrd1','Havcr2', 'Gzmb', 'Pdcd1', 'Prf1', 'S100a6','S100a4', "Nfkbid", 'Tnf', "Tnfrsf9", "Ifng", "Nr4a1", 'Xcl1', 'Il2', 'Il2ra','Il2rb','Klrk1', 'Cd44','Ccl3', 'Ccl4','Il10','Myc', 'Nr4a1', 'Zfpm1', "Myb"))
length(genes2plt)
add_genes_plot = genes2plt
idx = c(1:5,7:10,13, 14, 15)
gsets = sapply(mm_mc_endo_spec_ids[idx], scdb_gset)
all_genes = unique(c(unlist(sapply(gsets, function(x){return(names(x@gene_set))}, simplify = T)), add_genes_plot))
good_genes = all_genes

tfs_genes = c('Tbx21', 'Prdm1', 'Id2', 'Tox', 'Zbtb38', 'Spry2', 'Zfp683', 'Ikzf3', 'Cers4', 'Cebpb', 'Eld1', 'Stat3', 'Runx3', 'Batf', 'Stat4', 'Stat1', 'Tcf7', 'Id3', 'Eomes', 'Klf2', 'Klf3', 'Foxo1', 'Rela', 'Nfkbid', 'Irf8', 'Irf4', 'Zbtb32', 'Spry1', 'Nr4a1', 'Myc', 'Eea1', 'Pou2f2','Zeb2' ,'Rel', 'Myb', 'Zfpm1', 'Nfatc1', 'Ikzf2', 'Ikzf1', 'Pcna', 'Mki67')
marks_genes =c('Pdcd1', 'Havcr2', 'Tigit', 'Lag3', 'Cd244', 'Entpd1', 'Ctla4', 'Cxcr6', 'Cx3cr1', 'Klrc1', 'Klrd1', 'Klrk1', 'Il10ra','Il2rb', 'Prf1', 'Gzmb', 'Tnfrsf9', 'Xcl1', 'Ccl4', 'Ccl3', 'Tnf', 'Ifng', 'Myc', 'Il2', 'Il10', 'Cd44', 'Il2ra','Cd69', 'Sell', 'Slamf6', 'Dapl1', 'S1pr1', 'Lef1', 'Itgae','Cxcr3', 'Il7r', 'Ccr7', 'Ccl5', 'Itga1', 'Gzma', 'Gzmk', 'Nrgn', 'Mki67', 'Pcna', 'Top2a', 'Hmgb2')
tfs_marks_genes = c(tfs_genes, marks_genes)
i = 5
id = mm_endo_spec_ids[i]
mc_id = mm_mc_endo_spec_ids[i]
message(id)

# first time - run without coloring, then color by supmc and marks files
#conf = colorize_by_confusion_mat(mc_id=mc_id, graph_id=id, show_mc_ids=T)
sup = confs[[id]]$mc_sup
length(sup)
mc = scdb_mc(mc_id)
col2group = get_mc_col2group(mc)
group2col = get_mc_group2col(mc)

mmdysf_plots(idx=c(5))
mmdysf_plots = function(idx=c(5, 8:10), main_mc_grps_ls=c('spec', 'byst'), col2group=col2group, group2col=group2col) 
{
  mf1 = c('mouse_id', 'location', 'cell_type', 'batch_set_id', 'days_post_transfer')
  mf2 = c('mouse_id', 'location', 'batch_set_id')
  
  
  for (i in idx) {
    
    id = mm_endo_spec_ids[i]
    mc_id = mm_mc_endo_spec_ids[i]
    message(id)
    mc = scdb_mc(mc_id)
    
    if (i ==5){
      rl()
      mel_basic_mc_mc2d_plots(mc_id, mc_id, id, id, mm_lateral_gset_id, metadata_fields_to_export=mf1, ignore_edges=ignore_edges)
    }
    if (i ==11){
      rl(conf_f = "config/mm_dysf_mc38.yaml")
      mel_basic_mc_mc2d_plots(mc_id, mc_id, id, id, mm_lateral_gset_id, metadata_fields_to_export=mf1, ignore_edges=NULL)
    }
    else{
      rl(conf_f = "config/mm_dysf_6.yaml")
      mel_basic_mc_mc2d_plots(mc_id, mc_id, id, id, mm_lateral_gset_id, metadata_fields_to_export=mf1, ignore_edges=NULL)
    }
    
    
  }
  
  for (i in idx) {
    print(i)
    id = mm_endo_spec_ids[i]
    mc_id = mm_mc_endo_spec_ids[i]
    mc = scdb_mc(mc_id)
    
    main_mc_grps = col2group[names(table(mc@colors))]
    
    if (i == 1){
      nbars=87
      mmdysf_mc_group_composition_barplots(mc_id, id, set='pd1_byst', groups=main_mc_grps, nbars=nbars, exp_batch=FALSE, ann_colors = list(location=c(spleen="burlywood1", tumor="darkgoldenrod4", LN="brown", cd8_in_vitro="purple"), batch_set_id=c(PD1_5="darkblue",Ctrl_5="lightblue", naive_transfer="orange"),cell_type=c(spec="red",byst="green", cd8_in_vitro="purple")))
    }
    else if (i == 5){
      nbars=436
      mmdysf_mc_group_composition_barplots(mc_id, id, set='pd1_all', groups=main_mc_grps[!names(main_mc_grps) %in% c("red", "brown")], nbars=nbars, exp_batch=FALSE, col2group = col2group, group2col=group2col) #[!names(main_mc_grps) %in% c("red", "brown")]
    }
    else if (i == 2){
      nbars=87
      mmdysf_mc_group_composition_barplots(mc_id, id, set='pd1_byst', groups=main_mc_grps, nbars=nbars, exp_batch=FALSE, ann_colors = list(location=c(spleen="burlywood1", tumor="darkgoldenrod4", LN="brown", cd8_in_vitro="purple"), batch_set_id=c(PD1_5="darkblue",Ctrl_5="lightblue", cd8_in_vitro="orange"),cell_type=c(spec="red",byst="green", cd8_in_vitro="purple")))
    }
    else if (i == 4){
      nbars=50
      mmdysf_mc_group_composition_barplots(mc_id, id, set='pd1_byst', groups=main_mc_grps, nbars=nbars, exp_batch=FALSE, ann_colors = list(location=c(spleen="burlywood1", tumor="darkgoldenrod4", LN="brown", cd8_in_vitro="purple"), batch_set_id=c(PD1_treatment="darkblue",PD1_ctrl="lightblue", naive_transfer="orange"),cell_type=c(spec="red",endo="green", cd8_in_vitro="purple")))
    }
    else if (i == 13){
      nbars=50
      mmdysf_mc_group_composition_barplots(mc_id, id, set='pd1_byst', groups=main_mc_grps, nbars=nbars, exp_batch=FALSE, ann_colors = list(location=c(spleen="burlywood1", tumor="darkgoldenrod4", LN="brown", cd8_in_vitro="purple"), batch_set_id=c(MC38_PD1_2="darkblue",MC38_Ctrl_2="lightblue", MC38_41BB_2="purple",'MC38_41BB+PD1_2'="springgreen3"),cell_type=c(spec="red",endo="green", cd8_in_vitro="purple")))
    }
    else if (i == 14){
      nbars=12
      mmdysf_mc_group_composition_barplots(mc_id, id, set='pd1_byst', groups=main_mc_grps, nbars=nbars, exp_batch=FALSE, ann_colors = list(location=c(spleen="burlywood1", tumor="darkgoldenrod4", LN="brown", cd8_in_vitro="purple"), batch_set_id=c(pd1_8="darkblue",Ctrl_8="lightblue", '41bb_8'="purple",'41bb+pd1_8'="springgreen3"),cell_type=c(spec="red",endo="green", cd8_in_vitro="purple")), col2group = col2group, group2col=group2col)
    }
    else if (i == 15){
      nbars=12
      mmdysf_mc_group_composition_barplots(mc_id, id, set='pd1_byst', groups=main_mc_grps, nbars=nbars, exp_batch=FALSE, ann_colors = list(location=c(spleen="burlywood1", tumor="darkgoldenrod4", LN="brown", cd8_in_vitro="purple"), batch_set_id=c(Zbdb32_ko_Ctrl="darkblue",Zbdb32_ko_pd1="lightblue", BFP_Ctrl="purple",BFP_pd1="springgreen3"),cell_type=c(spec="red",endo="green", cd8_in_vitro="purple")), col2group = col2group, group2col=group2col)
    }
    else {
      nbars=20
      mmdysf_mc_group_composition_barplots(mc_id, id, set='pd1_byst', groups=main_mc_grps, nbars=nbars, exp_batch=FALSE)
    }
    
    my_mel_plot_e_gc_barplots(mc_id, id, genes=add_genes_plot[add_genes_plot %in% rownames(mc@e_gc)], ord_first_by_color = T)
    my_mel_plot_e_gc_barplots(mc_id, 'combined_tfs', genes=tfs_genes[tfs_genes %in% rownames(mc@e_gc)], ord_first_by_color = T, ncolumns = 2)
    
    for (ff in mf1) { 
      for (single_plot in c(T, F)) {
        tgconfig::set_param("mcell_mc2d_plot_key", F, "metacell")
        mcell_mc2d_plot_by_factor(mc_id, id, ff, single_plot=single_plot) 
        tgconfig::set_param("mcell_mc2d_plot_key", T, "metacell")
      }
    }
  }
  for (i in idx) {
    id = mm_endo_spec_ids[i]
    mc_id = mm_mc_endo_spec_ids[i]
    
    mc = scdb_mc(mc_id)
    lfp = log2(mc@mc_fp)
    #col2group = get_mc_col2group(mc)
    #group2col = get_mc_group2col(mc)
    
    
    genes_2d = c('Tcf7', 'Il7r', 'Lef1', 'Sell', 'Gzmk', 'Ccl5', 'Ccr2', 'Ccr7', 'Cxcr6', 'Nr4a2', 'Nr4a3', 'Id2', 'Prf1', 'Gzmb', 'Pdcd1', 'Tigit', 'Tnfrsf9', 'Havcr2', 'Il2ra', 'Il2', 'Il10', 'Il10ra', 'Tnf', 'Irf8', 'Ccl3', 'Ccl4', 'Tox', 'Lag3', 'Tcf7', 'Ctla4', 'Nr4a1', 'Myc', 'Gzma')
    genes_2d = all_genes
    for (gene in genes_2d){
      if (gene %in% rownames(lfp)){
        mcell_mc2d_plot_gene(mc_id, gene, show_legend = F, show_mc_ids = F)
      }
    }
    
    mcell_mc_export_tab(mc_id = mc_id, gstat_id = id, mat_id = id, T_fold=2, metadata_fields=mf1)
    
    n_top_genes = 200
    
    
    f = col2group[mc@colors] %in% main_mc_grps
    lfp_f = lfp[, f]
    model_nm = mc_id
    gg = names(tail(sort(apply(apply(lfp_f, 1, range), 2, diff)), n_top_genes))
    cc = cor(t(lfp_f[gg, ]))
    hc = hclust(dist(cc), method='ward.D2')
    pheatmap(cc[hc$order, hc$order], cluster_cols=F, cluster_rows=F, breaks=seq(-1, 1, length=100), main=model_nm, filename=sprintf("%s/%s_top%d_genes_over_%s_lfp_cor.png", .scfigs_base, model_nm, n_top_genes, paste0(mm_main_grps, collapse="-")), width=12, height=12)
    geoms = sapply(c(mc), function(v) v@e_gc)
    
    # Facs plots
    fn = paste("figs/", mc_id,"_mc_table_barplot1.png",sep="")
    png(fn, w=1200, h=500)
    barplot(table(mc@colors), names.arg = col2group[names(table(mc@colors))], col=names(table(mc@colors)))
    dev.off()
    
    # for (grp in names(table(mc@colors))){
    #   mcs = which(mc@colors==grp)
    #   fn = paste("figs/facs/", mc_id, "_",col2group[grp], "_GFP.png",sep="")
    #   png(fn, w=1200, h=500)
    #   plot(sort(log10(mat@cell_metadata[names(mc@mc[mc@mc %in% mcs]), "bystander.GFP.A_Ab"])), col=mc@colors[mc@mc[mc@mc %in% mcs]], ylim=c(0,5.5))
    #   abline(h=log10(1500), col="red")
    #   dev.off()
    #   fn = paste("figs/facs/", mc_id, "_", col2group[grp], "_TIM3.png",sep="")
    #   png(fn, w=1200, h=500)
    #   plot(sort(log10(mat@cell_metadata[names(mc@mc[mc@mc %in% mcs]), "TIM3.BV421.A_Ab"])), col=mc@colors[mc@mc[mc@mc %in% mcs]], ylim=c(0,5.5))
    #   #abline(h=log10(1500), col="red")
    #   dev.off()
    #   fn = paste("figs/facs/", mc_id, "_", col2group[grp], "_PD1.png",sep="")
    #   png(fn, w=1200, h=500)
    #   plot(sort(log10(mat@cell_metadata[names(mc@mc[mc@mc %in% mcs]), "PD1.APC.A_Ab"])), col=mc@colors[mc@mc[mc@mc %in% mcs]], ylim=c(0,5.5))
    #   #abline(h=log10(1500), col="red")
    #   dev.off()
    # }
    # 
    # library(RColorBrewer)
    # rf = colorRampPalette(rev(brewer.pal(11,'Spectral')))
    # r = rf(32)
    # mat@cell_metadata$color = mc@colors[mc@mc[rownames(mat@cell_metadata)]]
    # mat@cell_metadata$annot = col2group[mc@colors[mc@mc[rownames(mat@cell_metadata)]]]
    # facs_plot = mat@cell_metadata[names(mc@mc), ] %>% 
    #   group_by(batch_set_id, days_post_transfer, cell_type, annot) %>% do(plots=ggplot(data=.) + aes(x = log10(TIM3.BV421.A_Ab), y = log10(PD1.APC.A_Ab)) + stat_bin2d() + scale_fill_gradientn(colours=r, trans="log") + facet_grid(. ~ batch_set_id + days_post_transfer + cell_type + annot) + xlim(0,5) + ylim(0,5))
    # for (ii in 1:length(facs_plot$plots)){
    #   ggsave(filename = paste("figs/facs/density/", mc_id, "_", facs_plot$batch_set_id[[ii]], facs_plot$days_post_transfer[[ii]], facs_plot$cell_type[[ii]], facs_plot$annot[[ii]],"with_col.png",sep="_"), plot = facs_plot$plots[[ii]])
    # }
    # facs_plot = mat@cell_metadata[names(mc@mc), ] %>% 
    #   group_by(batch_set_id, days_post_transfer, cell_type, annot) %>% do(plots=ggplot(data=.) + aes(x = log10(TIM3.BV421.A_Ab), y = log10(PD1.APC.A_Ab)) + geom_point() + facet_grid(. ~ batch_set_id + days_post_transfer + cell_type + annot) + xlim(0,5) + ylim(0,5))
    # for (ii in 1:length(facs_plot$plots)){
    #   ggsave(filename = paste("figs/facs/scatter/", mc_id, "_", facs_plot$batch_set_id[[ii]], facs_plot$days_post_transfer[[ii]], facs_plot$cell_type[[ii]], facs_plot$annot[[ii]],"scatter.png",sep="_"), plot = facs_plot$plots[[ii]])
    # }
    
    #load_all("/net/mraid14/export/data/users/atanay/proj/metac/metacell")
  }
  saved_genes = all_genes
  for (i in idx) {
    id = mm_endo_spec_ids[i]
    mc_id = mm_mc_endo_spec_ids[i]
    
    mc = scdb_mc(mc_id)
    lfp = log2(mc@mc_fp)
    saved_genes = intersect(saved_genes, rownames(mc@mc_fp))
  }
  for (i in idx) {
    id = mm_endo_spec_ids[i]
    mc_id = mm_mc_endo_spec_ids[i]
    
    mc = scdb_mc(mc_id)
    lfp = log2(mc@mc_fp)
    #col2group = get_mc_col2group(mc)
    #group2col = get_mc_group2col(mc)
    
    if (i == 5){
      good_genes = mcell_mc_plot_marks(mc_id, mc_id, mat_id = mm_filt_id,
                                       fig_fn = NULL, lateral_gset_id = NULL,
                                       mc_ord = order(mc@colors),#confs[[id]]$mc_hc$order,#confs[[id]]$mc_hc$order
                                       plot_cells = T,
                                       zero_median = F,
                                       add_genes = NULL,
                                       add_metadata=c("location"),
                                       ext_metadata=NULL,
                                       md_level_colors = NULL,
                                       focus_mcs = NULL,
                                       gene_list = saved_genes,
                                       reorder_marks=T)
      col2ord = c(10,8,7,2,5,1,3,4,6,9)
      names(col2ord) = names(table(mc@colors))
      good_genes = mcell_mc_plot_marks(mc_id, mc_id, mat_id = mm_filt_id,
                                       fig_fn = "figs/mcell_marks_with_cell_type.png", lateral_gset_id = NULL,
                                       mc_ord = order(col2ord[mc@colors]),#confs[[id]]$mc_hc$order,#confs[[id]]$mc_hc$order
                                       plot_cells = F,
                                       zero_median = F,
                                       add_genes = NULL,
                                       add_metadata=NULL,#c("cell_type"),
                                       ext_metadata=NULL,
                                       md_level_colors = NULL,
                                       focus_mcs = NULL,
                                       gene_list = saved_genes,
                                       reorder_marks=T, fold_burn=2)
    }
    else{
      mcell_mc_plot_marks(mc_id, mc_id, mat_id = mm_filt_id,
                          fig_fn = NULL, lateral_gset_id = NULL,
                          mc_ord = confs[[id]]$mc_hc$order,
                          plot_cells = T,
                          zero_median = F,
                          add_genes = NULL,
                          add_metadata=c("cell_type"),
                          ext_metadata=NULL,
                          md_level_colors = NULL,
                          focus_mcs = NULL,
                          gene_list = NULL,
                          reorder_marks=F, fold_burn=3)
    }
    
    #rl()
  }
  idx = c(5, 22)
  idx= c(rev(idx))
  mcs = sapply(mm_mc_endo_spec_ids[idx], scdb_mc)
  lfps = sapply(mcs, function(v) log2(v@mc_fp))
  geoms = sapply(mcs, function(v) log2(v@e_gc + 1e-5))
  col2groups = sapply(mcs, get_mc_col2group)
  group2cols = sapply(mcs, get_mc_group2col)
  
  n_top_genes = 100
  
  # mc = scdb_mc(mc_id)
  # lfp = log2(mc@mc_fp)
  # col2group = get_mc_col2group(mc)
  # group2col = get_mc_group2col(mc)
  # f = col2group[mc@colors] %in% mm_main_grps
  # lfp_f = lfp[, f]
  # model_nm = mc_id
  # gg = names(tail(sort(apply(apply(lfp_f, 1, range), 2, diff)), n_top_genes))
  # cc = cor(t(lfp_f[gg, ]))
  # hc = hclust(dist(cc), method='ward.D2')
  # pheatmap(cc[hc$order, hc$order], cluster_cols=F, cluster_rows=F, breaks=seq(-1, 1, length=100), main=model_nm, filename=sprintf("%s/%s_top%d_genes_over_%s_lfp_cor.png", .scfigs_base, model_nm, n_top_genes, paste0(mm_main_grps, collapse="-")), width=8, height=8)
  # geoms = sapply(c(mc), function(v) v@e_gc)
  
  
  for (i in seq_along(lfps)) {
    #col2group = col2groups[[i]]
    mc = mcs[[i]]
    
    main_mc_grps = col2group[names(table(mc@colors))]
    
    
    f = col2group[mc@colors] %in% main_mc_grps
    lfp_f = lfps[[i]][, f]
    model_nm = mm_model_nms[names(lfps)[i]]
    gg = names(tail(sort(apply(apply(lfp_f, 1, range), 2, diff)), n_top_genes))
    cc = cor(t(lfp_f[gg, ]))
    hc = hclust(dist(cc), method='ward.D2')
    pheatmap(cc[hc$order, hc$order], cluster_cols=F, cluster_rows=F, breaks=seq(-1, 1, length=100), main=model_nm, filename=sprintf("%s/%s_top%d_genes_over_%s_lfp_cor.png", .scfigs_base, model_nm, n_top_genes, paste0(mm_main_grps, collapse="-")), width=8, height=8)
  }
  
  .cross_model_plt = function(id, g1, g2) 
  {
    .plot_start(sprintf("%s/geom_2d_plots/plt_geom_%d_%s_%s.png", .scfigs_base, id, g1, g2), 300 * length(geoms), 330)
    par(mfrow=c(1, length(geoms)))
    par(mar=c(4,4,4,1))
    xlim = range(as.vector(sapply(geoms, function(geom) { if (g1 %in% rownames(geom)) { range(geom[g1, ]) } else { NA }})), na.rm=T)
    ylim = range(as.vector(sapply(geoms, function(geom) { if (g2 %in% rownames(geom)) { range(geom[g2, ]) } else { NA }})), na.rm=T)
    for (i in seq_along(lfps)) { #c(1,3,2,4,5)
      geom = geoms[[i]]
      mc = mcs[[i]]
      model_nm = mm_model_nms[names(geoms)[i]]
      if (length(intersect(c(g1, g2), rownames(geom))) == 2) {
        plt(g1, g2, geom, mc@colors, xlim=xlim, ylim=ylim, main=model_nm)
      } else {
        plot(1:2, col=NA, xlab=g1, ylab=g2, main=model_nm, xlim=xlim, ylim=ylim)
      }
    }
    dev.off()
  }
  
  gs1 = c('Cd3d', 'Il7r', 'Sell', 'Gzma', 'Id2',  'Ccl3', 'Ccl5',  'Nkg7', 'Klrd1', 'S100a4', 'Pdcd1',  'Lag3', 'Mki67', 'Gzmk', 'Il7r','Mlana', 'Havcr2', 'Prf1', 'Prf1', 'Prf1', 'Nfkbid', 'Serpinb9', 'Nr4a1', 'Il2', 'Klf2','Tnf','Ccl1','Pdcd1', 'Sell', "Xcr1", "Xcr1", 'Clec9a', 'Cd8a', 'Tigit', 'Xcl1','Ccr7',"Il7r",'Tcf7','Ccl5')
  gs2 = c('Cd8a', 'Tcf7', 'Klf2', 'Gzmb', 'Gzmk', 'Ccl4', 'Cxcr6', 'Prf1', 'Klrk1', 'S100a6', 'Havcr2', 'Tigit', 'Top2a', 'Gzma', 'Ccl5','Cd74', 'Pdcd1', 'Pdcd1', 'Havcr2', 'Gzmb', 'Myc', 'Serpinb6b', 'Nr4a3', 'Nr4a2', 'Itgb1','Csf2','Ccl9','Csf1', 'S1pr1', "Cd40", 'Clec9a', 'Cd8a', 'Cd4', 'Cd4' , 'Xcr1','Sell','Tcf7','Klrk1','Tigit')
  for (i in seq_along(gs1)) {
    message(sprintf("%s vs %s", gs1[i], gs2[i]))
    .cross_model_plt(i, gs1[i], gs2[i])
  }
  
  
}

mmdysf_cluster_genes_over_mcs_1_model = function(groups=mm_main_grps, min_log_geomean=-15, T_enr=0.5, reg=2**-16, zlim=4, filter_genes_re=c('^Rp|^Gm[0-9]'),ids=1) 
{
  mc = scdb_mc(mm_mc_endo_spec_ids[ids])
  geom = mc@e_gc
  col2group = get_mc_col2group(mc)
  group2col = get_mc_group2col(mc)
  
  gg = rownames(geom)
  if (!is.null(filter_genes_re)) {
    gg = grep(filter_genes_re, gg, inver=T, perl=T, v=T)
  }
  
  models = c()
  annts = c()
  orig_ids = c()
  m = NULL
  
  f = col2group[mc@colors] %in% groups
  m = cbind(m, geom[gg, f])
  models = c(models, rep(mm_model_nms[mm_mc_endo_spec_ids[ids]], sum(f)))
  annts = c(annts, col2group[mc@colors[f]])
  orig_ids = c(orig_ids, seq_along(mc@colors)[f])
  
  
  colnames(m) = 1:ncol(m)
  
  # filter genes 
  ug = apply(m, 1, mean)
  ord = order(ug)
  
  m_n = log2( (m + reg) / apply(m + reg, 1, median))
  
  top_enr = apply(m_n, 1, quantile, prob=0.95)
  top_enr_ord = top_enr[ord]
  cmin = median(top_enr_ord[1:101])
  cmax = median(top_enr_ord[(length(top_enr_ord)-101):length(top_enr_ord)])
  enr_trend = zoo::rollmedian(top_enr_ord, 101, fill=c(cmin, NA, cmax))
  
  top_enr_n = top_enr
  top_enr_n[ord] = top_enr[ord] - enr_trend
  
  f = log2(ug)>= min_log_geomean & top_enr_n >= T_enr
  message(sprintf("%d genes passed", sum(f)))
  m_n_f = m_n[f, ]
  
  mc_ann = data.frame(row.names=colnames(m_n_f), model=models, grp=annts, orig_ids=orig_ids)
  png(sprintf("%s/genes_geom_clust_%s.png", .scfigs_base, paste0(groups, collapse="_")), width=200+4*ncol(m_n_f), height=200+8 * nrow(m_n_f))
  rhc = hclust(dist(cor(t(m_n_f))), method="ward.D2")
  chc = hclust(dist(cor(m_n_f)), method="ward.D2")
  pheatmap(pmin(pmax(m_n_f[rhc$order, chc$order], -zlim), zlim), cluster_rows=F, cluster_cols=F, annotation_col=mc_ann[, c('model', 'grp')], cellwidth=4, cellhight=8, annotation_colors=list(grp=group2col))
  dev.off()
  
  
  #list(m=m_n_f, ann=mc_ann)
}

mmdysf_cluster_genes_over_mcs = function(groups=mm_main_grps, min_log_geomean=-15, T_enr=0.5, reg=2**-16, zlim=4, filter_genes_re=c('^Rp|^Gm[0-9]')) 
{
  mcs = sapply(mm_mc_endo_spec_ids, scdb_mc)
  geoms = sapply(mcs, function(v) v@e_gc)
  col2groups = sapply(mcs, get_mc_col2group)
  group2cols = sapply(mcs, get_mc_group2col)
  
  gg = names(which(table(unlist(sapply(geoms, rownames))) == length(mm_mc_endo_spec_ids)))
  if (!is.null(filter_genes_re)) {
    gg = grep(filter_genes_re, gg, inver=T, perl=T, v=T)
  }
  
  models = c()
  annts = c()
  orig_ids = c()
  m = NULL
  for (i in seq_along(mcs)) {
    mc = mcs[[i]]
    f = col2groups[[i]][mc@colors] %in% groups
    m = cbind(m, geoms[[i]][gg, f])
    models = c(models, rep(mm_model_nms[names(mcs)[i]], sum(f)))
    annts = c(annts, col2groups[[i]][mc@colors[f]])
    orig_ids = c(orig_ids, seq_along(mc@colors)[f])
  }
  colnames(m) = 1:ncol(m)
  
  # filter genes 
  ug = apply(m, 1, mean)
  ord = order(ug)
  
  m_n = log2( (m + reg) / apply(m + reg, 1, median))
  
  top_enr = apply(m_n, 1, quantile, prob=0.95)
  top_enr_ord = top_enr[ord]
  cmin = median(top_enr_ord[1:101])
  cmax = median(top_enr_ord[(length(top_enr_ord)-101):length(top_enr_ord)])
  enr_trend = zoo::rollmedian(top_enr_ord, 101, fill=c(cmin, NA, cmax))
  
  top_enr_n = top_enr
  top_enr_n[ord] = top_enr[ord] - enr_trend
  
  f = log2(ug)>= min_log_geomean & top_enr_n >= T_enr
  message(sprintf("%d genes passed", sum(f)))
  m_n_f = m_n[f, ]
  
  mc_ann = data.frame(row.names=colnames(m_n_f), model=models, grp=annts, orig_ids=orig_ids)
  png(sprintf("%s/genes_geom_clust_%s.png", .scfigs_base, paste0(groups, collapse="_")), width=200+4*ncol(m_n_f), height=200+8 * nrow(m_n_f))
  rhc = hclust(dist(cor(t(m_n_f))), method="ward.D2")
  chc = hclust(dist(cor(m_n_f)), method="ward.D2")
  pheatmap(pmin(pmax(m_n_f[rhc$order, chc$order], -zlim), zlim), cluster_rows=F, cluster_cols=F, annotation_col=mc_ann[, c('model', 'grp')], cellwidth=4, cellhight=8, annotation_colors=list(grp=group2cols[[1]]))
  dev.off()
  
  
  list(m=m_n_f, ann=mc_ann)
}

pd1_time_course_plots = function(n_top_genes = 100, min_max_lfp=1, main_mc_grps=mm_main_grps) 
{
  lfp_f_cc = list()
  for (i in seq_along(lfps)) {
    col2group = col2groups[[i]]
    mc = mcs[[i]]
    f = col2group[mc@colors] %in% main_mc_grps
    lfp = lfps[[i]][, f]
    lfp = lfp[apply(lfp, 1, max) >= min_max_lfp, ]
    model_nm = mm_model_nms[names(lfps)[i]]
    gg = names(tail(sort(apply(apply(lfp, 1, range), 2, diff)), n_top_genes))
    cc = cor(t(lfp[gg, ]), method = 'spearman')
    lfp_f_cc[[names(lfps)[i]]] = cc
    
    max_g_lfp = apply(lfp[gg, ], 1, max)
    rownames(cc) = colnames(cc) = sprintf("%s (%.2f)", gg, max_g_lfp)
    
    hc = hclust(dist(cc), method='ward.D2')
    pheatmap(cc[hc$order, hc$order], cluster_cols=F, cluster_rows=F, breaks=seq(-1, 1, length=100), main=model_nm, filename=sprintf("%s/%s_top%d_genes_over_%s_lfp_cor.png", .scfigs_base, model_nm, n_top_genes, paste0(main_mc_grps, collapse="-")), width=n_top_genes/5, height=n_top_genes/5)
  }
  
  # joint endo-spec PD1 genes corr
  le = lfp_f_cc[[mm_mc_iva_5_spec_id]]
  ls = lfp_f_cc[[mm_mc_iva_5_byst_id]]
  endo_spec_g = intersect(rownames(le), rownames(ls))
  njg = length(endo_spec_g)
  
  le_j = le[endo_spec_g, endo_spec_g]
  ls_j = ls[endo_spec_g, endo_spec_g]
  
  le_hc = hclust(dist(le_j), method='ward.D2')
  ls_hc = hclust(dist(ls_j), method='ward.D2')
  pheatmap(le_j[le_hc$order, le_hc$order], cluster_cols=F, cluster_rows=F, breaks=seq(-1, 1, length=100), main="Spec by spec", filename=sprintf("%s/PD1_Spec_joint_byst_%d_genes_by_spec_over_%s_lfp_cor.png", .scfigs_base, njg, paste0(main_mc_grps, collapse="-")), width=max(6, njg/5), height=max(6, njg/5))
  pheatmap(ls_j[le_hc$order, le_hc$order], cluster_cols=F, cluster_rows=F, breaks=seq(-1, 1, length=100), main="Spec by endo", filename=sprintf("%s/PD1_Byst_joint_spec_%d_genes_by_spec_over_%s_lfp_cor.png", .scfigs_base, njg, paste0(main_mc_grps, collapse="-")), width=max(6, njg/5), height=max(6, njg/5))
  pheatmap(le_j[ls_hc$order, ls_hc$order], cluster_cols=F, cluster_rows=F, breaks=seq(-1, 1, length=100), main="Endo by spec", filename=sprintf("%s/PD1_Spec_joint_byst_%d_genes_by_byst_over_%s_lfp_cor.png", .scfigs_base, njg, paste0(main_mc_grps, collapse="-")), width=max(6, njg/5), height=max(6, njg/5))
  pheatmap(ls_j[ls_hc$order, ls_hc$order], cluster_cols=F, cluster_rows=F, breaks=seq(-1, 1, length=100), main="Spec by spec", filename=sprintf("%s/PD1_Byst_joint_spec_%d_genes_by_byst_over_%s_lfp_cor.png", .scfigs_base, njg, paste0(main_mc_grps, collapse="-")), width=max(6, njg/5), height=max(6, njg/5))
  
  
  #list(lfp_f_cc=lfp_f_cc)
}

pd1_time_course_plots_1_g = function(n_top_genes = 100, min_max_lfp=1, main_mc_grps=mm_main_grps, ids=1) 
{
  # mmdysf_mc_group_composition_barplots(mm_mc_endo_spec_ids[ids], mm_endo_spec_ids[ids],  groups=main_mc_grps, set='pd1_days')
  # 
  # for (gene in c('Tcf7', 'Il7r', 'Lef1', 'Sell', 'Gzmk', 'Ccl5', 'Ccr2', 'Ccr7', 'Cxcr6', 'Nr4a2', 'Nr4a3', 'Id2', 'Prf1', 'Gzmb', 'Pdcd1', 'Tigit', 'Tnfrsf9', 'Havcr2')) {
  #   for (i in ids) {
  #     mcell_mc2d_plot_gene(mm_mc_endo_spec_ids[i], gene=gene, show_legend=F)
  #   }
  # }
  
  mc = scdb_mc(mm_mc_endo_spec_ids[ids])
  geom = mc@e_gc
  col2grp = get_mc_col2group(mc)
  group2col = get_mc_group2col(mc)
  lfp = log2(mc@mc_fp)
  
  f = col2grp[mc@colors] %in% main_mc_grps
  lfp = lfp[, f]
  lfp = lfp[apply(lfp, 1, max) >= min_max_lfp, ]
  model_nm = mm_model_nms[mm_mc_endo_spec_ids[ids]]
  gg = names(tail(sort(apply(apply(lfp, 1, range), 2, diff)), n_top_genes))
  cc = cor(t(lfp[gg, ]), method = 'spearman')
  #lfp_f_cc[mm_mc_endo_spec_ids[ids]] = cc
  
  max_g_lfp = apply(lfp[gg, ], 1, max)
  rownames(cc) = colnames(cc) = sprintf("%s (%.2f)", gg, max_g_lfp)
  
  hc = hclust(dist(cc), method='ward.D2')
  pheatmap(cc[hc$order, hc$order], cluster_cols=F, cluster_rows=F, breaks=seq(-1, 1, length=100), main=model_nm, filename=sprintf("%s/%s_top%d_genes_over_%s_lfp_cor.png", .scfigs_base, model_nm, n_top_genes, paste0(main_mc_grps, collapse="-")), width=n_top_genes/5, height=n_top_genes/5)
  
  
}

mmdysf_pd1_gene_diff_expr_by_cond = function(min_max_umi=0, n_genes_to_show=20, min_enr=2, min_cells=20,ids=c(5,8:10)) 
{
  
  mcs = sapply(c(mm_mc_endo_spec_ids[ids]), scdb_mc)
  col2grps = sapply(mcs, get_mc_col2group)
  mats = sapply(c(mm_endo_spec_ids[ids]), scdb_mat)
  mats_ds = sapply(mats, function(mat) { scm_downsamp(mat@mat, scm_which_downsamp_n(mat))})
  
  
  # plot func 
  pd1_inner_diff_expr_plot_func = function(de, min_enr, n_genes_to_show, mod_nm, ofn, xlab, ylab, min_max_umi = 8) 
  {
    png(scfigs_fn(mod_nm, ofn, scfigs_dir(mod_nm, 'diff_expr')), 800, 800)
    par(mar=c(4,4,8,8))
    x1 = log2(de$tot1 + 1)
    x2 = log2(de$tot2 + 1)
    names(x1) = names(x2) = de$gene
    ind = abs(de$enr) > min_enr
    lim = range(c(x1, x2))
    plot(x1, x2, ylab=sprintf('%s, %d cells (log2)', ylab, nrow(x2)), xlab=sprintf('%s, %d cells (log2)', xlab, nrow(x1)), xlim=lim, ylim=lim, pch=19, cex=0.5, col=ifelse(ind, 'red', 'darkgrey'), main=ofn)
    abline(a=0, b=1, col='black')
    
    enr_g = de %>% filter(enr > min_enr & (tot1 >= min_max_umi | tot2 >= min_max_umi)) %>% head(n_genes_to_show) %>% arrange(tot2)
    dpl_g = de %>% filter(enr < -min_enr & (tot1 >= min_max_umi | tot2 >= min_max_umi)) %>% tail(n_genes_to_show) %>% arrange(tot1)
    lab_off = (lim[2] - lim[1]) / 10 
    if (nrow(enr_g) > 0) {
      mtext(enr_g$gene, side=4, line=0.2, at=seq(lim[1] + lab_off, lim[2] - lab_off, length=nrow(enr_g)), las=2, cex=0.9)
      segments(x0=lim[2], y0=seq(lim[1] + lab_off, lim[2] - lab_off, length=nrow(enr_g)), x1=x1[enr_g$gene], y1=x2[enr_g$gene], col='grey')
      
    }
    if (nrow(dpl_g) > 0) {
      text(seq(lim[1] + lab_off, lim[2] - lab_off, length=nrow(dpl_g)), par("usr")[4] + par("csi"), offset=0, dpl_g$gene, srt = 30, xpd = TRUE, pos = 4, cex=0.9)
      
      segments(x0=seq(lim[1] + lab_off, lim[2] - lab_off, length=nrow(dpl_g)), y0=lim[2], x1=x1[dpl_g$gene], y1=x2[dpl_g$gene], col='grey')
    }
    dev.off()
  }
  
  mc = scdb_mc(mm_mc_endo_spec_ids[ids])
  geom = mc@e_gc
  col2grp = get_mc_col2group(mc)
  group2col = get_mc_group2col(mc)
  mat_m = scdb_mat(mm_endo_spec_ids[ids])
  mat_ds = scm_downsamp(mat_m@mat, scm_which_downsamp_n(mat_m))
  
  model_nms = c('Spec')
  grps_filt = list(Spec=mm_main_grps)
  
  # CTRL vs PD1 by model, location, day and mc_grp
  des = list()
  for (i in seq_along(model_nms)) {
    md = mat_m@cell_metadata[names(mc@mc), ]
    md$mc_grp = col2grp[mc@colors[mc@mc]]
    
    mod_nm = model_nms[i]
    for (loc in c('tumor')) {
      for (grp in grps_filt[[mod_nm]]) {
        for (day in c(3,6,9,12,15)) { #as.numeric(unique(as.character(md$days_post_transfer)))
          message(sprintf("%s, %s, %s, %d", mod_nm, loc, grp, day))
          c1 = md %>% tibble::rownames_to_column('cell_name') %>%
            filter(location == loc & mc_grp == grp & days_post_transfer == day & batch_set_id == 'PD1_5') 
          c2 = md %>% tibble::rownames_to_column('cell_name') %>% 
            filter(location == loc & mc_grp == grp & days_post_transfer == day & batch_set_id == 'Ctrl_5')
          if (nrow(c1) >= min_cells && nrow(c2) >= min_cells) {
            de = diff_expr(mc, mat_ds, mcs1=NULL, mcs2=NULL, min_max_umi=min_max_umi, geo_mean=T, nms1=c1$cell_name, nms2=c2$cell_name)
            nm = sprintf("%s_%s_d%d_PD1_vs_ctrl_%s", mod_nm, loc, day, grp)
            des[[nm]] = de
            pd1_inner_diff_expr_plot_func(de, min_enr, n_genes_to_show, mod_nm, nm, xlab="PD1", ylab="Control") 
          }
        }
      }
    }
  }
  
  
  # Day 6 vs day 3 by model, location, PD1/ctrl and mc_grp
  for (i in seq_along(model_nms)) {
    md = mat_m@cell_metadata[names(mc@mc), ]
    md$mc_grp = col2grp[mc@colors[mc@mc]]
    
    mod_nm = model_nms[i]
    for (loc in c('tumor')) {
      for (grp in grps_filt[[mod_nm]]) {
        for (cond in c('PD1_5', 'Ctrl_5')) {
          message(sprintf("%s, %s, %s, %s", mod_nm, loc, grp, cond))
          c1 = md %>% tibble::rownames_to_column('cell_name') %>%
            filter(location == loc & mc_grp == grp & batch_set_id == cond & days_post_transfer == 6) 
          c2 = md %>% tibble::rownames_to_column('cell_name') %>% 
            filter(location == loc & mc_grp == grp & batch_set_id == cond & days_post_transfer == 3)
          if (nrow(c1) >= min_cells && nrow(c2) >= min_cells) {
            de = diff_expr(mc, mat_ds, mcs1=NULL, mcs2=NULL, min_max_umi=min_max_umi, geo_mean=T, nms1=c1$cell_name, nms2=c2$cell_name)
            nm = sprintf("%s_%s_%s_d6_vs_d3_%s", mod_nm, loc, cond, grp)
            des[[nm]] = de
            pd1_inner_diff_expr_plot_func(de, min_enr, n_genes_to_show, mod_nm, nm, xlab='d6', ylab='d3') 
          }
        }
      }
    }
  }
  
  
  
  return(des)
}

mmdysf_cluster_genes_over_mcs_1_model()
#c('naive1', 'cycling_trans','cycling_trans_','trans-dysf1', 'trans-dysf2','na','em')
pd1_time_course_plots(main_mc_grps = c('trans-dysf3'))
pd1_time_course_plots_1_g(main_mc_grps = c('cycling-trans'), ids = 3)
des = mmdysf_pd1_gene_diff_expr_by_cond()



for (de_ind in 1:length(des)){
  de_nm = names(des)[de_ind]
  de_path = paste("figs/Spec.diff_expr/", de_nm, ".csv", sep="")
  de = des[[de_ind]]
  de = de[abs(de$enr) > 1,]
  print(dim(de)[1])
  write.csv(de, de_path)
}
de_path = paste("figs/Spec.diff_expr/", "trans-dysf1_max_lfp_diff_genes.csv", sep="")

i = 5
id = mm_endo_spec_ids[i]
mc_id = mm_mc_endo_spec_ids[i]
mc = scdb_mc(mc_id)

library("anndata")
mat = scdb_mat(id)
sc = import("scanpy")
mtx_fn = paste("scrna_db", "/cd8_cells.mtx", sep="")
h5ad_fn = paste("scrna_db", "/cd8_cells.h5ad", sep="")
Matrix::writeMM(obj = t(mat@mat), file = mtx_fn)
mtx_mat = anndata::read_mtx(mtx_fn)
mtx_mat$var_names = rownames(mat@mat)
mtx_mat$obs_names = colnames(mat@mat)
sc$write(adata = mtx_mat, filename = h5ad_fn)

mc_2 = anndata::read_h5ad('metacells_vignette/metacells_2.h5ad')
mc_2_cells = anndata::read_h5ad('metacells_vignette/cells_2.h5ad')
mc2_cols = read.csv('metacells_vignette/cluster-colors_2.csv')
mc_clusts = mc_2$obs$cluster
mc_colors = mc2_cols$color[mc_clusts]
cells_mc = mc_2_cells$obs$metacell
cells_mc = cells_mc + 1
cell_names = mc_2_cells$obs_names
names(cells_mc) = cell_names
gene_names = mc_2$var_names

outlier_cells = names(cells_mc)[cells_mc == 0]
cells_mc_filt = cells_mc[cells_mc > 0]

umis2 = as.matrix(mc_2$X)
fractions = umis2 / rowSums(umis2)
fractions = t(fractions)

fake_tgMC = mc

filt_fake_egc = fractions
colnames(filt_fake_egc) = 1:max(cells_mc_filt)

fake_tgMC@mc = cells_mc_filt
fake_tgMC@outliers = outlier_cells   #cell_names[mc_per_cell == 0]
fake_tgMC@cell_names = colnames(mat@mat)
fake_tgMC@mc_fp = as.matrix(filt_fake_egc)
fake_tgMC@e_gc = as.matrix(filt_fake_egc)
fake_tgMC@cov_gc = as.matrix(filt_fake_egc)
fake_tgMC@n_bc = as.matrix(t(as.data.frame(colSums(mat@mat))))
fake_tgMC@annots = colnames(filt_fake_egc)
fake_tgMC@colors = as.character(mc_colors)
mc_id_2 = paste(mc_id, "_2", sep="")
scdb_add_mc(mc_id_2, fake_tgMC)

mc_new = scdb_mc(mc_id_2)
mc = scdb_mc(mc_id)
mc_id_2 = mc_id

sp = import("scipy.sparse")
mc2 = import("metacells")
umap_edges_2 = sp$coo_matrix(mc2$ut$get_oo_proper(mc_2, 'obs_outgoing_weights'))

umap_edges_2_sparse = 1*as.matrix((umap_edges_2 > 0))


mc2d = scdb_mc2d(mc_id)
mc2d_2 = mc2d

mc2d_2@mc_id = mc_id_2
umap_x = mc2$ut$get_o_numpy(mc_2, 'umap_x')
umap_y = mc2$ut$get_o_numpy(mc_2, 'umap_y')
mc2d_2@mc_x = umap_x
mc2d_2@mc_y = umap_y
scdb_add_mc2d(mc_id_2,mc2d_2)

mc = scdb_mc(mc_id)
lfp = log2(1e-5 + mc@mc_fp)
egc = mc@e_gc
col2group = get_mc_col2group(mc)
col2grp = col2group
group2col = get_mc_group2col(mc)
mc2d = scdb_mc2d(mc_id)
md = mat@cell_metadata[names(mc@mc),]

tf_tb = read.csv("figs/mouse_tf_list.csv")
tf_list = as.character(tf_tb$TFs) #toupper(tf_tb$Tf)_

tf_tb_human = read.csv("figs/human_TFs.csv")
tf_list_human = as.character(tf_tb_human$Name) #toupper(tf_tb$Tf)_

mouse2human_df = read.csv(paste(scfigs_dir, "/mouse2human.csv", sep=""), stringsAsFactors = F)
mouse2human = mouse2human_df$Human.gene.name
names(mouse2human) = mouse2human_df$Gene.name
human2mouse = mouse2human_df$Gene.name
names(human2mouse) = mouse2human_df$Human.gene.name

tf_list_human_m = human2mouse[tf_list_human]
tf_list_human_m = tf_list_human_m[!is.na(tf_list_human_m)]
tf_list_human_m = tf_list_human_m[tf_list_human_m %in% rownames(egc)]
tf_list_human_m_add = tf_list_human_m[!tf_list_human_m %in% tf_list]

tf_list_all_non_lat = unique(c(tf_list, tf_list_human_m_add))
tf_list_all_non_lat = tf_list_all_non_lat[!tf_list_all_non_lat %in% all_lat_genes]

surfs_tb = read.csv("figs/surfs.csv")
surfs_list = as.character(surfs_tb$Surf)

shades = colorRampPalette(c("blue", "white", "red"))(1001)
day_tab = table(mc@mc, md$days_post_transfer)
#d_ratio = (day_tab[,1])/(day_tab[,1] + day_tab[,2])
d_ratio = (rowSums(day_tab[,c("0","1","2","3")]))/(rowSums(day_tab[,c("0","1","2","3")]) + rowSums(day_tab[,c("4","5","6","9","12","15")]))
d_colors = floor(d_ratio * 1000)
png("figs/spec_day3_day6_mc2d.png", w=400,h=400)
plot(mc2d@mc_x,mc2d@mc_y, pch=21, bg=shades[d_colors], cex=2, main="Day3 (red) vs Day6 (blue)")
dev.off()
#text(mc2d@mc_x,mc2d@mc_y, labels = 1:ncol(lfp))

b_tab = table(mc@mc, md$batch_set_id)
b_ratio = (b_tab[,"PD1_5"])/(b_tab[,"PD1_5"] + b_tab[,"Ctrl_5"])
b_colors = floor(b_ratio * 1000)
png("figs/spec_Pd1_Ctrl_mc2d.png", w=400,h=400)
plot(mc2d@mc_x,mc2d@mc_y, pch=21, bg=shades[b_colors], cex=2, main="Pd1 (red) vs Ctrl (blue)")
dev.off()
#write.csv(tail(sort(td1_diff),50), de_path)

# mmdysf_build_cc_gsets = function() 
# {
#   # build gset to mask by the data
#   select_gene_modules_by_anchor_genes(mm_filt_id, c('Mki67', 'Smc4'), gset_nm="mmdysf_cc_g2_m", cor_thresh=0.1, sz_cor_thresh=0.1, nclusts = 5)
#   
#   select_gene_modules_by_anchor_genes(mm_filt_id, c('Pcna', 'Mcm3'), gset_nm="mmdysf_cc_g1_s", cor_thresh=0.1, sz_cor_thresh=0.1, nclusts = 5)
# 
#   mcell_gset_remove_clusts("mmdysf_cc_g2_m", filt_clusts=c(2,4), new_id = "mmdysf_cc_g2_m_filt", reverse=T)
#   mcell_gset_remove_clusts("mmdysf_cc_g1_s", filt_clusts=c(4), new_id = "mmdysf_cc_g1_s_filt", reverse=T)
#   
#   
# }
# mmdysf_build_cc_gsets()
# 
# cc_gset = scdb_gset("mmdysf_cc_filt")
# cc_genes = names(cc_gset@gene_set)
# cc_genes = cc_genes[cc_genes %in% rownames(lfp)]
# 
# cc_g2_gset = scdb_gset("mmdysf_cc_g2_m_filt")
# cc_g2_genes = names(cc_g2_gset@gene_set)
# cc_g2_genes = cc_g2_genes[cc_g2_genes %in% rownames(lfp)]
# 
# cc_g1_gset = scdb_gset("mmdysf_cc_g1_s_filt")
# cc_g1_genes = names(cc_g1_gset@gene_set)
# cc_g1_genes = cc_g1_genes[cc_g1_genes %in% rownames(lfp)]
# 
# gene_cor_number = 40
# gene_cor = cor(t(lfp))
# # cc_g1_genes = unique(names(tail(sort(gene_cor["Mcm3",]),gene_cor_number)), names(tail(sort(gene_cor["Pcna",]),gene_cor_number)))
# # cc_g2_genes = unique(names(tail(sort(gene_cor["Mki67",]),gene_cor_number)), names(tail(sort(gene_cor["Smc4",]),gene_cor_number)))
# cc_genes = names(non_filt_gset@gene_set)[non_filt_gset@gene_set %in% c(10,23,24,36)]#c(17,43,46,48,21,12,5,20)] #c(5,12,17,20,21)
#cc_g1_genes = names(non_filt_gset@gene_set)[non_filt_gset@gene_set %in% c(24)]#c(17,43,46,48,21,12,5,20)] #c(5,12,17,20,21)
cc_g1_genes = c("Pcna", "Rrm2", "Mcm5", "Mcm6", "Mcm4", "Ung", "Mcm7", "Mcm2","Uhrf1", "Orc6", "Tipin")
cc_g2_genes = c("Mki67","Cenpf","Top2a","Smc4","Ube2c","Ccnb1","Cdk1","Arl6ip1","Ankrd11","Hmmr","Cenpa","Tpx2","Aurka","Kif4", "Kif2c","Bub1b","Ccna2", "Kif23","Kif20a","Sgol2","Smc2", "Kif11", "Cdca2","Incenp","Cenpe") 
cc_g1_genes = intersect(cc_g1_genes, rownames(lfp))
cc_g2_genes = intersect(cc_g2_genes, rownames(lfp))
cc_genes = c(cc_g1_genes, cc_g2_genes)
#cc_g2_genes = names(non_filt_gset@gene_set)[non_filt_gset@gene_set %in% c(36)]#c(17,43,46,48,21,12,5,20)] #c(5,12,17,20,21)

cc_sig = apply(lfp[cc_genes,],2, mean)
cc_g2_sig = apply(lfp[cc_g2_genes,],2, mean)
cc_g1_sig = apply(lfp[cc_g1_genes,],2, mean)
byst_cells = names(mc@mc)[mat@cell_metadata[names(mc@mc), "cell_type"]=="byst"]
spec_cells = names(mc@mc)[mat@cell_metadata[names(mc@mc), "cell_type"]=="spec"]
ccg2_byst = mean(cc_g2_sig[mc@mc[byst_cells]])
ccg1_byst = mean(cc_g1_sig[mc@mc[byst_cells]])
ccg2_spec = mean(cc_g2_sig[mc@mc[spec_cells]])
ccg1_spec = mean(cc_g1_sig[mc@mc[spec_cells]])

plot(lfp["Pcna",], lfp["Mcm3",], pch=21, cex=3, bg=mc@colors, xlab="Pcna", ylab="Mcm3", cex.lab=1, main="Plt 1", xlim=range(lfp["Pcna",]), ylim=range(lfp["Mcm3",]))
text(lfp["Pcna",], lfp["Mcm3",], colnames(lfp), cex=3/4)
plot(lfp["Mki67",], lfp["Smc4",], pch=21, cex=3, bg=mc@colors, xlab="Mki67", ylab="Smc4", cex.lab=1, main="Plt 2", xlim=range(lfp["Mki67",]), ylim=range(lfp["Smc4",]))
text(lfp["Mki67",], lfp["Smc4",], colnames(lfp), cex=3/4)

png("figs/ccg1_vs_ccg2.png", w=600,h=600)
plot(cc_g2_sig, cc_g1_sig, pch=21, cex=3, bg=mc@colors, xlab="cc_g2_sig", ylab="cc_g1_sig", cex.lab=1, main="cc_g2_sig vs cc_g1_sig", xlim=range(cc_g2_sig), ylim=range(cc_g1_sig))
text(cc_g2_sig, cc_g1_sig, colnames(lfp), cex=3/4)
dev.off()

plot(lfp["Il2ra",], cc_sig, pch=21, cex=3, bg=mc@colors, xlab="Il2ra", ylab="cc_sig", cex.lab=1, main="Il2ra vs cc_sig", xlim=range(lfp["Il2ra",]), ylim=range(cc_sig))
text(lfp["Il2ra",], cc_sig, colnames(lfp), cex=3/4)
plot(lfp["Il2ra",], cc_g2_sig, pch=21, cex=3, bg=mc@colors, xlab="Il2ra", ylab="cc_g2_sig", cex.lab=1, main="Il2ra vs cc_g2_sig", xlim=range(lfp["Il2ra",]), ylim=range(cc_g2_sig))
text(lfp["Il2ra",], cc_g2_sig, colnames(lfp), cex=3/4)
plot(lfp["Il2ra",], cc_g1_sig, pch=21, cex=3, bg=mc@colors, xlab="Il2ra", ylab="cc_g1_sig", cex.lab=1, main="Il2ra vs cc_g1_sig", xlim=range(lfp["Il2ra",]), ylim=range(cc_g1_sig))
text(lfp["Il2ra",], cc_g1_sig, colnames(lfp), cex=3/4)

plt('Pcna', 'Mki67', lfp, mc@colors)

Il2_cor = apply(lfp, 1, cor, lfp["Il2",])
cor_path = paste("figs/Spec.diff_expr/", "Il2_top_100_cor_genes.csv", sep="")
write.csv(tail(sort(Il2_cor),100), cor_path)

Il2ra_cor = apply(lfp, 1, cor, lfp["Il2ra",])
cor_path = paste("figs/Spec.diff_expr/", "Il2ra_top_100_cor_genes.csv", sep="")
write.csv(tail(sort(Il2ra_cor),100), cor_path)

mc_td1 = (1:max(mc@mc))[mc@colors=="#e2abd6"]
mc_td2 = (1:max(mc@mc))[mc@colors=="gold"]
mc_cycling = (1:max(mc@mc))[mc@colors=="#1f78b4"]
td1_mean  = apply(lfp[,mc_td1], 1, mean)
td2_mean  = apply(lfp[,mc_td2], 1, mean)
cycling_mean  = apply(lfp[,mc_cycling], 1, mean)
rest_mean_td2  = apply(lfp[,-c(mc_td2)], 1, mean)
rest_mean_central_mem  = apply(lfp[,-c(mc_td1)], 1, mean)
td1_diff = td1_mean - rest_mean_td2
td2_diff = td2_mean - rest_mean_central_mem
de_path = paste("figs/Spec.diff_expr/", "trans-dysf2_max_lfp_diff_genes.csv", sep="")
write.csv(tail(sort(td2_diff),50), de_path)
de_path = paste("figs/Spec.diff_expr/", "em2_max_lfp_diff_genes.csv", sep="")
write.csv(tail(sort(td1_diff),50), de_path)

mc_td1 = (1:max(mc@mc))[mc@colors=="#e2abd6"]
mc_td2 = (1:max(mc@mc))[mc@colors=="gold"]
mc_cycling = (1:max(mc@mc))[mc@colors=="#1f78b4"]
td1_mean  = apply(egc[,mc_td1], 1, sum)/ table(mc@colors)["#e2abd6"]
td2_mean  = apply(egc[,mc_td2], 1, sum) / table(mc@colors)["gold"]
cycling_mean  = apply(egc[,mc_cycling], 1, sum)/ table(mc@colors)["#1f78b4"]

td2_diff = td2_mean - cycling_mean
td2_diff_tf = td2_diff[names(td2_diff) %in% tf_list]
de_path = paste("figs/Spec.diff_expr/", "cycling_cytotoxic_egc_diff_sum_genes.csv", sep="")
write.csv(c(tail(sort(td2_diff),50), head(sort(td2_diff),50)), de_path)

de_path = paste("figs/Spec.diff_expr/", "cycling_cytotoxic_egc_diff_mean_tf_genes.csv", sep="")
write.csv(c(tail(sort(td2_diff_tf),10), head(sort(td2_diff_tf),10)), de_path)


ss = embdyn_cc_plots(mm_filt_id, mc_id, mc2d_id=mc_id, tag=mm_filt_id, thr_score = 0.005,m_score_scaling = 0.7, m_genes=cc_g2_genes, s_genes=cc_g1_genes)
embdyn_cc_plots = function(mat_id, mc_id,  mc2d_id=mc_id, tag=mat_id,thr_score = 0.005,m_score_scaling = 0.7, m_genes=NULL, s_genes=NULL)
{
  # output: function returns names of cells which are growth-arrested, i.e. lye below a threshold line
  
  if (is.null(m_genes)){
    m_genes = c("Mki67","Cenpf","Top2a","Smc4;SMC4","Ube2c","Ccnb1","Cdk1","Arl6ip1","Ankrd11","Hmmr;IHABP","Cenpa;Cenp-a","Tpx2","Aurka","Kif4", "Kif2c","Bub1b","Ccna2", "Kif23","Kif20a","Sgol2","Smc2", "Kif11", "Cdca2","Incenp","Cenpe")  
  }
  else{
    m_genes = unique(c("Mki67","Cenpf","Top2a","Smc4;SMC4","Ube2c","Ccnb1","Cdk1","Arl6ip1","Ankrd11","Hmmr;IHABP","Cenpa;Cenp-a","Tpx2","Aurka","Kif4", "Kif2c","Bub1b","Ccna2", "Kif23","Kif20a","Sgol2","Smc2", "Kif11", "Cdca2","Incenp","Cenpe", m_genes))
  }
  
  if (is.null(s_genes)){
    #s_genes = c("Pcna", "Rrm2", "Mcm5", "Npm1", "Mcm6", "Mcm4", "Ung", "Mcm7", "Mcm2","Uhrf1", "Orc6", "Tipin")
    s_genes = c("Pcna", "Rrm2", "Mcm5", "Mcm6", "Mcm4", "Ung", "Mcm7", "Mcm2","Uhrf1", "Orc6", "Tipin") # Npm1
  }
  else{
    s_genes = unique(c("Pcna", "Rrm2", "Mcm5", "Mcm6", "Mcm4", "Ung", "Mcm7", "Mcm2","Uhrf1", "Orc6", "Tipin", s_genes))
  }
  
  
  m = scdb_mat(mat_id)
  mc = scdb_mc(mc_id)
  mc2d = scdb_mc2d(mc2d_id)
  
  s_genes = intersect(rownames(mc@mc_fp), s_genes)
  m_genes = intersect(rownames(mc@mc_fp), m_genes)
  
  tot  = colSums(m@mat)
  s_tot = colSums(m@mat[s_genes,])
  m_tot = colSums(m@mat[m_genes,])
  
  s_score = s_tot/tot
  m_score = m_tot/tot
  
  if(!dir.exists(sprintf("figs/cc_%s",tag))) {
    dir.create(sprintf("figs/cc_%s", tag))
  }
  
  f = s_score +  m_score_scaling*m_score < thr_score
  
  
  # Next try to fit a smooth spline
  png(sprintf("figs/cc_%s/cc_scores.png", tag), w=600, h=600)
  plot(s_score, m_score, pch=19, cex=0.1)
  points(s_score[f], m_score[f], pch=19, cex=0.1, col="darkred")
  dev.off()
  shades = colorRampPalette(c("white","lightgray","lightblue","blue", "red", "yellow"))
  png(sprintf("figs/cc_%s/cc_sm_scores.png", tag), w=600, h=600)
  smoothScatter(s_score, m_score, colramp=shades, pch=19, cex=0.1)
  abline(a = thr_score/m_score_scaling,b = - 1/m_score_scaling)
  dev.off()
  
  mc_cc_tab = table(mc@mc, f[names(mc@mc)])
  mc_cc = 1+floor(99*mc_cc_tab[,2]/rowSums(mc_cc_tab))
  
  shades = colorRampPalette(c("white","lightblue", "blue", "purple"))(100)
  png(sprintf("figs/cc_%s/2d_cc.png", tag), w=800, h=800)
  plot(mc2d@sc_x, mc2d@sc_y, pch=19, cex=0.4, col=ifelse(f[names(mc2d@sc_x)], "black", "lightgray"))
  points(mc2d@mc_x, mc2d@mc_y, pch=21, cex=2.5, bg=shades[mc_cc])
  dev.off()
  
  png(sprintf("figs/cc_%s/bars_cc.png", tag), w=1200, h=500)
  barplot(mc_cc, col=mc@colors, las=2, cex.names=0.7)
  dev.off()
  
  
  return(colnames(mat@mat)[f])
}



#### --- Supervised analysis of gene programs ---- #####
#load_all("/net/mraid14/export/data/users/atanay/proj/metac/metacell")

mc = scdb_mc(mc_id)
lfp = log2(mc@mc_fp)
egc = mc@e_gc

f_cytox_h = mc@colors=="gold"
c_cytox_h = apply(lfp[,f_cytox_h], 1, cor, lfp["Id2", f_cytox_h])
f_dys_trans = mc@colors=="#975e9f"
c_dys_trans = apply(lfp[,f_dys_trans], 1, cor, lfp["Pdcd1", f_dys_trans])
f_cytokine = mc@colors=="#e2abd6"
c_cytokine = apply(lfp[,f_cytokine], 1, cor, lfp["Ccl3", f_cytokine])
f_cytokine_nfkb = mc@colors=="deeppink"
c_cytokine_nfkb = apply(lfp[,f_cytokine_nfkb], 1, cor, lfp["Nfkbid", f_cytokine_nfkb])
f_cycling_trans = mc@colors=="#1f78b4"
c_cycling_trans = apply(lfp[,f_cycling_trans], 1, cor, lfp["Nrgn", f_cycling_trans])
f_central_mem = mc@colors=="#6a8255"
c_central_mem = apply(lfp[,f_central_mem], 1, cor, lfp["Dapl1", f_central_mem])
f_resident_mem = mc@colors=="blue"
c_resident_mem = apply(lfp[,f_resident_mem], 1, cor, lfp["Tcf7", f_resident_mem])
f_dysf = mc@colors=="magenta"
c_dysf = apply(lfp[,f_dysf], 1, cor, lfp["Tigit", f_dysf])
f_dysf_low = mc@colors=="cyan"
c_dysf_low = apply(lfp[,f_dysf_low], 1, cor, lfp["Tigit", f_dysf_low])
f_naive1 = mc@colors=="#5fb8f4" #
c_naive1 = apply(lfp[,f_naive1], 1, cor, lfp["Il7r", f_naive1])
f_em1 = mc@colors=="#cee6b9"
c_em1 = apply(lfp[,f_em1], 1, cor, lfp["Ccr2", f_em1])
f_tum = mc@colors=="brown"
c_tum = apply(lfp[,f_tum], 1, cor, lfp["Mlana", f_tum])

tail(sort(c_cytox_h),40)
tail(sort(c_dys_trans),40)
tail(sort(c_cycling_trans),40)
tail(sort(c_central_mem),40)
tail(sort(c_resident_mem),40)

e_dys_trans = rowMeans(lfp[,f_dys_trans])
plot(c_dys_trans, e_dys_trans, cex=0.4, pch=19)
abline(h=0.4, col="red")
tail(sort(c_dys_trans[e_dys_trans>0.3])[!names(sort(c_dys_trans[e_dys_trans>0.3])) %in% core_genes ],60) #[!names(sort(c_dys_trans[e_dys_trans>0.4])) %in% core_genes]
max_dys_trans = apply(lfp[,f_dys_trans], 1, max)
which(max_dys_trans>0.8)
g_dys_trans = names(which(max_dys_trans>0.8))
max_dys_trans = apply(lfp[,f_dys_trans], 1, max)
pheatmap::pheatmap(pmax(pmin(lfp[g_dys_trans, f_dys_trans],1),-1))
lcore = lfp[g_dys_trans, f_dys_trans]
lcore = lcore[,order(lcore["Xcl1",])]
pheatmap::pheatmap(pmax(pmin(lcore,1),-1),cluster_cols=F)


mcell_mc_plot_gg(mc_id, "Prf1", "Pdcd1", use_egc=T)

e_cytox_h = rowMeans(lfp[,f_cytox_h])
tail(sort(e_cytox_h)[!names(sort(e_cytox_h)) %in% core_genes],60)
plot(c_cytox_h, e_cytox_h, cex=0.4, pch=19)
abline(h=0.4, col="red")
tail(sort(c_cytox_h[e_cytox_h>0.5])[!names(sort(c_cytox_h[e_cytox_h>0.5])) %in% core_genes],40) #[core_genes_ref_gene=="Prf1"], [!names(sort(c_cytox_h[e_cytox_h>0.1])) %in% core_genes & names(sort(c_cytox_h[e_cytox_h>0.1])) %in% tf_list]
max_cytox_h = apply(lfp[,f_cytox_h], 1, max)
which(max_cytox_h>0.8)
g_cytox_h = names(which(max_cytox_h>0.8))
max_cytox_h = apply(lfp[,f_cytox_h], 1, max)
pheatmap::pheatmap(pmax(pmin(lfp[g_cytox_h, f_cytox_h],1),-1))
lcore = lfp[g_cytox_h, f_cytox_h]
lcore = lcore[,order(lcore["Prf1",])]
pheatmap::pheatmap(pmax(pmin(lcore,1),-1),cluster_cols=F)


mcell_mc_plot_gg(mc_id, "Prf1", "Pdcd1", use_egc=T)
mcell_mc_plot_gg(mc_id, "Cd44", "Pdcd1", use_egc=T)
mcell_mc_plot_gg(mc_id, "Tcf7", "Slamf6", use_egc=T)
mcell_mc_plot_gg(mc_id, "Mki67", "Cd69", use_egc=T)
mcell_mc_plot_gg(mc_id, "Cd69", "Slamf6", use_egc=T)
mcell_mc_plot_gg(mc_id, "Cd69", "Tcf7", use_egc=T)
mcell_mc_plot_gg(mc_id, "Eomes", "Tbx21", use_egc=T)
mcell_mc_plot_gg(mc_id, "Havcr2", "Cxcr5", use_egc=T)
mcell_mc_plot_gg(mc_id, "Havcr2", "Cxcr6", use_egc=T)
mcell_mc_plot_gg(mc_id, "Cd44", "Pdcd1", use_egc=T)
mcell_mc_plot_gg(mc_id, "Prf1", "Pdcd1", use_egc=T)
mcell_mc_plot_gg(mc_id, "Prf1", "Tnfrsf9", use_egc=T)
mcell_mc_plot_gg(mc_id, "Prf1", "Havcr2", use_egc=T)
mcell_mc_plot_gg(mc_id, "Prf1", "Ccl4", use_egc=T)
mcell_mc_plot_gg(mc_id, "Prf1", "Ccl3", use_egc=T)
mcell_mc_plot_gg(mc_id, "Prf1", "Il2ra", use_egc=T)
mcell_mc_plot_gg(mc_id, "Prf1", "Il2rb", use_egc=T)
mcell_mc_plot_gg(mc_id, "Ccl3", "Capg", use_egc=T)
mcell_mc_plot_gg(mc_id, "Prf1", "Irf8", use_egc=T)
mcell_mc_plot_gg(mc_id, "Prf1", "Eomes", use_egc=T)
mcell_mc_plot_gg(mc_id, "Stat3", "Stat1", use_egc=T)
mcell_mc_plot_gg(mc_id, "Sell", "Il7r", use_egc=T)
mcell_mc_plot_gg(mc_id, "Tox", "Pdcd1", use_egc=T)

e_resident_mem = rowMeans(lfp[,f_resident_mem])
plot(c_resident_mem, e_resident_mem, cex=0.4, pch=19)
abline(h=0.3)
tail(sort(e_resident_mem),60)#[!names(sort(e_resident_mem)) %in% core_genes],60)
tail(sort(c_resident_mem[e_resident_mem>0.3])[!names(sort(c_resident_mem[e_resident_mem>0.3])) %in% core_genes],60)
max_resident_mem = apply(lfp[,f_resident_mem], 1, max)
which(max_resident_mem>0.6)
g_resident_mem = names(which(max_resident_mem>0.6))
max_resident_mem = apply(lfp[,f_resident_mem], 1, max)
pheatmap::pheatmap(pmax(pmin(lfp[g_resident_mem, f_resident_mem],1),-1))
lcore = lfp[g_resident_mem, f_resident_mem]
lcore = lcore[,order(lcore["Sell",])]
pheatmap::pheatmap(pmax(pmin(lcore,1),-1),cluster_cols=F)

e_central_mem = rowMeans(lfp[,f_central_mem])
plot(c_central_mem, e_central_mem, cex=0.4, pch=19)
abline(h=0.3)
tail(sort(e_central_mem)[!names(sort(e_central_mem)) %in% core_genes],60)
tail(sort(c_central_mem[e_central_mem>0.3])[!names(sort(c_central_mem[e_central_mem>0.3])) %in% core_genes & names(sort(c_central_mem[e_central_mem>0.3])) %in% tf_list],60)
max_central_mem = apply(lfp[,f_central_mem], 1, max)
which(max_central_mem>0.6)
g_central_mem = names(which(max_central_mem>0.6))
max_central_mem = apply(lfp[,f_central_mem], 1, max)
pheatmap::pheatmap(pmax(pmin(lfp[g_central_mem, f_dem],1),-1))
lcore = lfp[g_central_mem, f_central_mem2]
lcore = lcore[,order(lcore["Sell",])]
pheatmap::pheatmap(pmax(pmin(lcore,1),-1),cluster_cols=F)

mcell_mc_plot_gg(mc_id, "Dapl1", "Lsp1", use_egc=T)

e_cycling_trans = rowMeans(lfp[,f_cycling_trans])
plot(c_cycling_trans, e_cycling_trans, cex=0.4, pch=19)
abline(h=0.3)
tail(sort(e_cycling_trans)[!names(sort(e_cycling_trans)) %in% core_genes],60)
tail(sort(c_cycling_trans[e_cycling_trans>0.2])[!names(sort(c_cycling_trans[e_cycling_trans>0.2])) %in% core_genes ],60) #& names(sort(c_cycling_trans[e_cycling_trans>0.2])) %in% tf_list
max_cycling_trans = apply(lfp[,f_cycling_trans], 1, max)
which(max_cycling_trans>0.8)
g_cycling_trans = names(which(max_cycling_trans>0.8))
max_cycling_trans = apply(lfp[,f_cycling_trans], 1, max)
pheatmap::pheatmap(pmax(pmin(lfp[g_cycling_trans, f_cycling_trans],1),-1))
lcore = lfp[g_cycling_trans, f_cycling_trans]
lcore = lcore[,order(lcore["Tcf7",])]
pheatmap::pheatmap(pmax(pmin(lcore,1),-1),cluster_cols=F)

mcell_mc_plot_gg(mc_id, "Prf1", "Pdcd1", use_egc=T)

e_cytokine = rowMeans(lfp[,f_cytokine])
plot(c_cytokine, e_cytokine, cex=0.4, pch=19)
abline(h=1)
tail(sort(c_cytokine[e_cytokine>0.5])[!names(sort(c_cytokine[e_cytokine>0.5])) %in% core_genes & names(sort(c_cytokine[e_cytokine>0.5])) %in% tf_list],60)
max_cytokine = apply(lfp[,f_cytokine], 1, max)
length(which(max_cytokine>1.2))
g_cytokine = names(which(max_cytokine>1.2))
max_cytokine = apply(lfp[,f_cytokine], 1, max)
pheatmap::pheatmap(pmax(pmin(lfp[g_cytokine, f_cytokine],1),-1))
lcore = lfp[g_cytokine, f_cytokine]
lcore = lcore[,order(lcore["Ccl4",])]
pheatmap::pheatmap(pmax(pmin(lcore,1),-1),cluster_cols=F)

mcell_mc_plot_gg(mc_id, "Prf1", "Pdcd1", use_egc=T)

e_naive1 = rowMeans(lfp[,f_naive1])
plot(c_naive1, e_naive1, cex=0.4, pch=19)
abline(h=0.2)
tail(sort(e_naive1)[!names(sort(e_naive1)) %in% core_genes],60)
tail(sort(c_naive1[e_naive1>0.2])[!names(sort(c_naive1[e_naive1>0.2])) %in% core_genes ],60)
max_naive1 = apply(lfp[,f_naive1], 1, max)
length(which(max_naive1>0.8))
g_naive1 = names(which(max_naive1>0.8))
max_naive1 = apply(lfp[,f_naive1], 1, max)
pheatmap::pheatmap(pmax(pmin(lfp[g_naive1, f_naive1],1),-1))
lcore = lfp[g_naive1, f_naive1]
lcore = lcore[,order(lcore["Ccl4",])]
pheatmap::pheatmap(pmax(pmin(lcore,1),-1),cluster_cols=F)

mcell_mc_plot_gg(mc_id, "Prf1", "Pdcd1", use_egc=T)

e_em1 = rowMeans(lfp[,f_em1])
plot(c_em1, e_em1, cex=0.4, pch=19)
abline(h=0.8)
tail(sort(e_em1)[!names(sort(e_em1)) %in% core_genes],60)
tail(sort(c_em1[e_em1>0.5])[!names(sort(c_em1[e_em1>0.5])) %in% core_genes ],60)
max_em1 = apply(lfp[,f_em1], 1, max)
length(which(max_em1>1.2))
g_em1 = names(which(max_em1>1.2))
max_em1 = apply(lfp[,f_em1], 1, max)
pheatmap::pheatmap(pmax(pmin(lfp[g_em1, f_em1],1),-1))
lcore = lfp[g_em1, f_em1]
lcore = lcore[,order(lcore["Ccr2",])]
pheatmap::pheatmap(pmax(pmin(lcore,1),-1),cluster_cols=F)

mcell_mc_plot_gg(mc_id, "Prf1", "Pdcd1", use_egc=T)

e_tum = rowMeans(lfp[,f_tum])
plot(c_tum, e_tum, cex=0.4, pch=19)
abline(h=0.2)
head(sort(c_tum[e_tum>0.2]),60)
max_tum = apply(lfp[,f_tum], 1, max)
length(which(max_tum>0.8))
g_tum = names(which(max_tum>0.8))
max_tum = apply(lfp[,f_tum], 1, max)
pheatmap::pheatmap(pmax(pmin(lfp[g_tum, f_tum],1),-1))
lcore = lfp[g_tum, f_tum]
lcore = lcore[,order(lcore["Mlana",])]
pheatmap::pheatmap(pmax(pmin(lcore,1),-1),cluster_cols=F)

mcell_mc_plot_gg(mc_id, "Prf1", "Pdcd1", use_egc=T)
mcell_mc_plot_gg(mc_id, "Id2", "Id3", use_egc=T)

scfigs_dir = "figs"
# old_gp_list = read.csv(file=paste(scfigs_dir, "/gp_share_new_f.csv", sep=""))
# old_core_genes = as.character(old_gp_list$Gene)

gp_list = read.csv(file=paste(scfigs_dir, "/gp_full_model_final.csv", sep=""))
core_genes = as.character(gp_list$Gene)
core_genes_ref = as.character(gp_list$Related.Annot)
core_genes_ref = core_genes_ref[core_genes %in% rownames(lfp)]

core_genes_ref_gene = as.character(gp_list$Ref.Gene)
core_genes_ref_gene = core_genes_ref_gene[core_genes %in% rownames(lfp)]
gp_list = gp_list[core_genes %in% rownames(lfp),]
core_genes = core_genes[core_genes %in% rownames(lfp)]
gp_list$Tf = core_genes %in% tf_list
gp_list$Surf = core_genes %in% surfs_list
#write.csv(gp_list, file=paste(scfigs_dir, "/gp_share_tf.csv", sep=""))

gene_inds = 1:length(table(core_genes_ref_gene))
names(gene_inds) = names(table(core_genes_ref_gene))
gene_s = gene_inds[core_genes_ref_gene]
names(gene_s) = core_genes
g_id = "main_prog_unique_n"

scdb_add_gset(id=g_id, tgGeneSets(gene_s))

main_gset = scdb_gset(g_id)
names(main_gset@set_names) = names(table(core_genes_ref_gene))[main_gset@set_names]
main_gset@set_names = sort(main_gset@set_names)
scdb_add_gset(id=g_id, main_gset)
main_gset = scdb_gset(g_id)

filt_main_gset = main_gset@gene_set
names(core_genes_ref) = core_genes
names(core_genes_ref_gene) = core_genes

main_cols = c("#e2abd6", "#cee6b9", "brown", "red", "grey", "#5fb8f4", "orange", "indianred4", "deeppink", "#1f78b4", "yellow", "magenta", "gold", "blue")
names(main_cols) = names(table(core_genes_ref_gene))
#non_lat_gp = c("cycling-trans", "resident-memory", "cytotoxic-effector-high", "cytotoxic-effector-low", "effector-cytokines-trans", "naive-effector-memory-percursor", "central-memory-like", "naive", "dysfunctional", "short-lived-eff-nfkb", "short-lived-eff-cytokines")
non_lat_gp = c("memory", "cytotoxicity", "effector-cytokines", "effector-memory", "effector-nfkb", "immunomodulatory", "naive-like", "transitional")

# Look at lfp scores per gene modules:
filt_main_gset = filt_main_gset[core_genes]
filt_main_gset_ntf = filt_main_gset[!gp_list$Tf]
filt_main_gset_nsurf = filt_main_gset[!gp_list$Surf]
# Prog score based on mean lfp
prog_score = tgs_matrix_tapply(t(lfp[names(filt_main_gset_ntf),]), filt_main_gset_ntf, mean)
rownames(prog_score) = paste(names(sort(main_gset@set_names)), "_prog", sep="")
prog_score_with_tf = tgs_matrix_tapply(t(lfp[names(filt_main_gset),]), filt_main_gset, mean)
rownames(prog_score_with_tf) = paste(names(sort(main_gset@set_names)), "_prog", sep="")
# Prog score based on total egc
prog_score = tgs_matrix_tapply(t(egc[names(filt_main_gset_ntf),]), filt_main_gset_ntf, sum)
rownames(prog_score) = paste(names(sort(main_gset@set_names)), "_prog", sep="")
prog_score_with_tf = tgs_matrix_tapply(t(egc[names(filt_main_gset),]), filt_main_gset, sum)
rownames(prog_score_with_tf) = paste(names(sort(main_gset@set_names)), "_prog", sep="")
prog_score_with_nsurf = tgs_matrix_tapply(t(egc[names(filt_main_gset_nsurf),]), filt_main_gset_nsurf, sum)
rownames(prog_score_with_nsurf) = paste(names(sort(main_gset@set_names)), "_prog", sep="")
gp_cor_dir = paste(scfigs_dir, "/gp_cors_egc", sep="")
dir.create(gp_cor_dir,showWarnings = F)

rownames(prog_score) = c("eff-cytokines-prog", "em-prog", "dend-prog", "stress-prog", "interf1-prog", "naive-prog", "m-cc-prog", "tumor-prog", "eff-nfkb-prog", "trans-mem-prog", "s-cc-prog", "immunomodulatory-prog", "eff-cytotoxic-prog", "memory-prog")
rownames(prog_score_with_tf) = rownames(prog_score)

progs_umis = tgs_matrix_tapply(t(mat@mat[names(filt_main_gset),names(mc@mc)]), filt_main_gset, sum)
progs_umis = t(tgs_matrix_tapply(progs_umis, mc@mc, sum))
mc_tot_umis = colSums(progs_umis)
progs_depth = t(t(progs_umis) / mc_tot_umis)
rownames(progs_umis) = rownames(prog_score)
rownames(progs_depth) = rownames(prog_score)


colnames(prog_score) = 1:max(mc@mc)
colnames(prog_score_with_tf) = 1:max(mc@mc)
colnames(progs_umis) = 1:max(mc@mc)
colnames(progs_depth) = 1:max(mc@mc)
my_mel_plot_e_gc_barplots = function(mc_id, name, genes=NULL, ncolumns=2, panel_height=50, panel_width=300, ord_first_by_color=T, n_ideal_umi=1000, egc=NULL, col2group=NULL, ord_by_id=NULL) 
{
  mc = scdb_mc(mc_id)
  if (is.null(col2group)){
    col2group = get_mc_col2group(mc)
  }

  
  if (is.null(genes)) {
    marks_gset = scdb_gset(mc_id)
    genes = names(marks_gset@gene_set)
  }
  if (is.null(egc)){
    e_gc = mc@e_gc[genes, ] * n_ideal_umi
  }
  else{
    e_gc = egc[genes, ] * n_ideal_umi
    colnames(e_gc) = 1:max(mc@mc)
  }
  
  if (ord_first_by_color & !is.null(ord_by_id) ) {
    e_gc = e_gc[, order(as.numeric(ordered(col2group[mc@colors], levels=ord_by_id)) + as.numeric(colnames(e_gc)) * 1e-6)]
  }
  
  .plot_start(scfigs_fn(mc_id, sprintf("mc_geom_mean_%s", name)), w=ncolumns * panel_width, h=panel_height * ceiling(length(genes) / ncolumns))
  #png(scfigs_fn(mc_id, sprintf("mc_geom_mean_%s", name)), ncolumns * panel_width, panel_height * ceiling(length(genes) / ncolumns))
  layout(matrix(1:(length(genes) + length(genes) %% 2), ncol=ncolumns))
  par(mar=c(0.5, 12, 0.5, 1))
  
  for (g in genes) {
    barplot(e_gc[g, ], border=NA, col=mc@colors[as.numeric(colnames(e_gc))], xaxt='n', yaxt='n', space=0)
    yaxp = par("yaxp")
    axis(2, yaxp=c(yaxp[1], yaxp[2], 1), las=2, cex=1)
    mtext(g, 2, line=1.5, cex=1.5, las=2)
  }
  dev.off()
  
}
my_mel_plot_e_gc_barplots(mc_id, "main_prog_pd1", genes=rownames(prog_score), egc=prog_score,ord_first_by_color=F)
my_mel_plot_e_gc_barplots(mc_id, "main_prog_pd1_with_tfs", genes=rownames(prog_score_with_tf), egc=prog_score_with_tf,ord_first_by_color=F)
my_mel_plot_e_gc_barplots(mc_id, "main_prog_pd1_umis", genes=rownames(progs_depth), egc=progs_depth,ord_first_by_color=F)

my_mel_plot_e_gc_barplots(mc_id, "mc_v2", genes=rownames(progs_depth), egc=progs_depth,ord_first_by_color=F)

for (x_gp in 1:(dim(prog_score)[1]-1)){
  for (y_gp in (x_gp+1):dim(prog_score)[1]){
    
    png(paste(gp_cor_dir, "/gp_cor_", rownames(prog_score)[x_gp], "_", rownames(prog_score)[y_gp], ".png", sep=""), 500, 500)
    plot(prog_score[x_gp,], prog_score[y_gp,], col=mc@colors, pch=19, cex=2, xlab=rownames(prog_score)[x_gp], ylab=rownames(prog_score)[y_gp])
    text(prog_score[x_gp,], prog_score[y_gp,], labels=1:length(prog_score[y_gp,]), cex=0.5)
    dev.off()
  }
}

prog_score_log = log2(prog_score + 1e-6)
for (x_gp in 1:(dim(prog_score)[1]-1)){
  for (y_gp in (x_gp+1):dim(prog_score)[1]){
    
    png(paste(gp_cor_dir, "/gp_cor_", rownames(prog_score)[x_gp], "_", rownames(prog_score)[y_gp], "_log.png", sep=""), 500, 500)
    plot(prog_score_log[x_gp,], prog_score_log[y_gp,], col=mc@colors, pch=19, cex=2, xlab=rownames(prog_score)[x_gp], ylab=rownames(prog_score)[y_gp])
    text(prog_score_log[x_gp,], prog_score_log[y_gp,], labels=1:length(prog_score[y_gp,]), cex=0.5)
    dev.off()
  }
}

gp_hm_cor = cor(t(prog_score))
gp_hm_cor[is.na(gp_hm_cor)] = 0
png(paste(gp_cor_dir, "/gp_cor_hm.png", sep=""), 800, 800)
pheatmap(gp_hm_cor)
dev.off()

for (gp_ind in 1:dim(prog_score)[1]){
  gp_genes = names(filt_main_gset_ntf)[filt_main_gset_ntf==gp_ind]
  gp_genes = gp_genes[gp_genes %in% rownames(lfp)]
  anc_gene = names(main_gset@set_names)[gp_ind]
  gp_col = main_cols[gp_ind]
  if (length(gp_genes)>1){
    print(gp_ind)
    c_gp = sort(apply(lfp[gp_genes,], 1, cor, lfp[anc_gene,]))
    # c_gp = c_gp[c_gp>0.3]
    png(paste(gp_cor_dir, "/gp_gene_cor_", anc_gene, "_", gp_col, ".png", sep=""), 400, 1200)
    barplot(c_gp, col=gp_col, horiz = T, las = 2)
    dev.off()
  }
  
}

# Looks for TFs correlated to each programs
tf_genes = tf_list[tf_list %in% rownames(lfp)]
tf_lfp = lfp[tf_genes,]
tf_pg_cor = tgs_cor(t(tf_lfp), t(prog_score))
for (gp_ind in 1:dim(prog_score)[1]){
  gp_name = names(main_gset@set_names)[gp_ind]
  gp_col = main_cols[gp_ind]
  c_gp = sort(tf_pg_cor[,gp_ind])
  c_gp = c_gp[c_gp>0.6]
  if(length(c_gp > 0)){
    png(paste(gp_cor_dir, "/tfs_gp_cor_", gp_name, "_", gp_col, ".png", sep=""), 200, 600)
    barplot(c_gp, col=gp_col, horiz = T, las = 2)
    dev.off()
  }
  
}

# Looks for Surface proteins correlated to each programs
surf_genes = surfs_list[surfs_list %in% rownames(lfp)]
surf_lfp = lfp[surf_genes,]
surf_pg_cor = tgs_cor(t(surf_lfp), t(prog_score_with_nsurf))
for (gp_ind in 1:dim(prog_score)[1]){
  gp_name = names(main_gset@set_names)[gp_ind]
  gp_col = main_cols[gp_ind]
  c_gp = sort(surf_pg_cor[,gp_ind])
  c_gp = c_gp[c_gp>0.6]
  if(length(c_gp > 0)){
    png(paste(gp_cor_dir, "/surfs_gp_cor_", gp_name, "_", gp_col, ".png", sep=""), 200, 600)
    barplot(c_gp, col=gp_col, horiz = T, las = 2)
    dev.off()
  }
  
}
tf_genes = tf_genes[!tf_genes %in% surf_genes]



gene_ann = data.frame(group=core_genes_ref_gene)
rownames(gene_ann) = core_genes

# paletteLength <- 100
# myColor <- colorRampPalette(c("white","steelblue","red"))(paletteLength)#colorRampPalette(c("blue","green", "white", "orange","red"))(paletteLength)
# myBreaks <- c(0,seq(0.5, max(cor_topics), length.out=floor(paletteLength)-1))

heatmap_core_genes_path = paste(gp_cor_dir,"/core_genes_heatmap_no_lat_nc.png",sep="")
png(heatmap_core_genes_path, 6000, 6000)
cor_core = tgs_cor(t(egc[core_genes[core_genes_ref %in% non_lat_gp],])) #
#cor_core = cor_core[!is.na(cor_core)]
hm = pheatmap(cor_core, cluster_col=F, cluster_row=F, fontsize=25,annotation_col=gene_ann, annotation_row=gene_ann,annotation_colors=list(group=main_cols)) #
dev.off()
# save_pheatmap_pdf(hm, heatmap_cor_topics_path)

mc_sizes = table(mc@mc)
barplot(mc_sizes, col = mc@colors, las=2, cex.names = 0.5)

for (gp_ind in 1:dim(prog_score_with_tf)[1]){
  gp_genes = names(filt_main_gset)[filt_main_gset==gp_ind]
  anc_gene = names(main_gset@set_names)[gp_ind]
  gp_col = main_cols[gp_ind]
  if (length(gp_genes) > 1){
    c_gp = sort(apply(lfp[gp_genes,], 1, cor, lfp[anc_gene,]))
    # c_gp = c_gp[c_gp>0.3]
    genes_names = names(c_gp)
    genes_names[genes_names %in% tf_list] = paste(genes_names[genes_names %in% tf_list], "**", sep="")
    genes_names[genes_names %in% surfs_list] = paste(genes_names[genes_names %in% surfs_list], "++", sep="")
    png(paste(gp_cor_dir, "/with_tf_gp_gene_cor_", anc_gene, "_", gp_col, ".png", sep=""), 200, 600)
    barplot(c_gp, col=gp_col, horiz = T, las = 2, names.arg = genes_names)
    dev.off()
  }
  
}

gp_ind = 8
gp_genes = names(filt_main_gset_ntf)[filt_main_gset_ntf==gp_ind]
anc_gene = names(main_gset@set_names)[gp_ind]
gp_col = main_cols[gp_ind]
c_gp = sort(apply(lfp, 1, cor, lfp[anc_gene,]))

gp_ind = 2
gp_genes = names(filt_main_gset_ntf)[filt_main_gset_ntf==gp_ind]
anc_gene = names(main_gset@set_names)[gp_ind]
gp_col = main_cols[gp_ind]
c_gp_2 = sort(apply(lfp, 1, cor, lfp[anc_gene,]))


cells_prog = tgs_matrix_tapply(t(mat@mat[names(filt_main_gset),names(mc@mc)]), filt_main_gset, sum)
mcs_prog = tgs_matrix_tapply(cells_prog, mc@mc, sum)
mc_cells_depth = colSums(mat@mat[,names(mc@mc)])
colnames(cells_prog) = names(mc@mc)
rownames(cells_prog) = names(main_gset@set_names)
cells_col = mc@colors[mc@mc]
cells_prog_n = t(cells_prog)/mc_cells_depth
plot(cells_prog_n[,8], cells_prog_n[,9], col=cells_col)
abline(0,1)
quantile(cells_prog_n[,8], (0:20)/20)
quantile(cells_prog_n[,9], (0:20)/20)
cells_prf1 = cells_prog_n[cells_prog_n[,8]>0.0062,]
cells_sell = cells_prog_n[cells_prog_n[,9]>0.0114,]
both_cells = intersect(rownames(cells_prf1), rownames(cells_sell))
joint_cells = unique(c(rownames(cells_prf1), rownames(cells_sell)))
diff_cells = joint_cells[!joint_cells %in% both_cells]
md_diff = mat@cell_metadata[diff_cells,]
md_diff$prog = ifelse(rownames(md_diff) %in% rownames(cells_prf1 ), "Prf1", "Sell")


# gp_share = sort(core_genes_ref_gene)
# write.csv(gp_share, file=paste(scfigs_dir, "/gp_share.csv", sep=""))
for(ii in 1:length(core_genes)){
  grp = core_genes_ref[ii]
  gn2 = core_genes_ref_gene[ii]
  gn = core_genes[ii]
  fn = paste("figs/main_prog/egc/", grp, "_", gn2, "_", gn, ".png",sep="")
  #plt(grp, gn, egc, mc@colors, ofn=fn)
  mcell_mc_plot_gg(mc_id, gn2, gn, use_egc=T, fig_fn=fn, md_mode="f")
}

for(ii in 1:length(core_genes)){
  grp = core_genes_ref[ii]
  gn2 = core_genes_ref_gene[ii]
  gn = core_genes[ii]
  fn = paste("figs/main_prog/lfp/", grp, "_", gn2, "_", gn, ".png",sep="")
  #plt(grp, gn, egc, mc@colors, ofn=fn)
  mcell_mc_plot_gg(mc_id, gn2, gn, use_egc=F, fig_fn=fn, md_mode="f")
}

main_prog_mc_sum = tgs_matrix_tapply(mat@mat[names(filt_main_gset),names(mc@mc)], mc@mc, sum)
td2_prog_mc_sum = tgs_matrix_tapply(mat@mat[names(filt_main_gset)[filt_main_gset==14],names(mc@mc)], mc@mc, sum)
mc_sums = tgs_matrix_tapply(mat@mat[,names(mc@mc)], mc@mc, sum)

#mc_depths = rowSums(main_prog_mc_sum)
mc_depths = rowSums(mc_sums)
td2_sums = rowSums(td2_prog_mc_sum)
td2_per_mc = td2_sums / mc_depths
td2_per_mc = sort(td2_per_mc, decreasing=T)

td2_ordered = egc["Tcf7",names(td2_per_mc)]

for (gene_ind in 1:length(filt_main_gset[names(filt_main_gset) %in% tf_list])){ #[core_genes %in% tf_list]
  gene = names(filt_main_gset)[names(filt_main_gset) %in% tf_list][gene_ind]
  prog = core_genes_ref[gene]
  dir1 = paste("figs/main_prog/dynamics/TFs/", prog, "/", sep="")
  dir.create(dir1, showWarnings=F)
  fig_fn = paste(dir1, gene, ".png", sep="")
  png(fig_fn, w=600, h=600)
  gene_ordered = 2**(egc[gene,names(td2_per_mc)])
  plot(gene_ordered, bg=mc@colors[as.integer(names(td2_ordered))], ylab=gene, xlab="mc by Tcf7 prog expression order", pch=21, cex=2)
  text(gene_ordered, names(td2_ordered), cex=0.5)
  dev.off()
}

for (gene_ind in 1:length(filt_main_gset)){ #[core_genes %in% tf_list]
  gene = names(filt_main_gset)[gene_ind]
  prog = core_genes_ref[gene]
  dir1 = paste("figs/main_prog/dynamics/all_features/", prog, "/", sep="")
  dir.create(dir1, showWarnings=F)
  fig_fn = paste(dir1, gene, ".png", sep="")
  png(fig_fn, w=600, h=600)
  gene_ordered = 2**(egc[gene,names(td2_per_mc)])
  plot(gene_ordered, bg=mc@colors[as.integer(names(td2_ordered))], ylab=gene, xlab="mc by Tcf7 prog expression order", pch=21, cex=2)
  text(gene_ordered, names(td2_ordered), cex=0.5)
  dev.off()
}



# Tsne analysis from program sums features
a = tgs_matrix_tapply(t(mat@mat[names(filt_main_gset),]), filt_main_gset, sum) #[!filt_main_gset %in% c(3,4,5,9)] [!filt_main_gset %in% c(3,4,5,9)]
colnames(a) = colnames(mat@mat)

a = a[,names(mc@mc)]
#look at the distribution of total umi per program, e.g.:
quantile(a[5,], (0:20)/20)

af = a[,colSums(a)>0]
#normalize and project - no MCs are needed when having the "right" features, at least for first approximation

an = t(af)/colSums(af)
#ts = Rtsne::Rtsne(an)
#plot(ts$Y[,1], ts$Y[,2], cex=0.2, pch=19, col=mc@colors[mc@mc[colnames(mat@mat)]])



ts_exp = Rtsne::Rtsne(an[!duplicated(an),], perplexity=length(mc@mc)/100, eta=length(mc@mc)/12)#,perplexity=length(mc@mc)/100, eta=length(mc@mc)/12, num_threads=0) #, perplexity=ncol(mat@mat)/100, eta=ncol(mat@mat)/12, num_threads=0
tsne_fn = paste("figs/main_prog/", "tsne_exp_full.png", sep="")
png(tsne_fn, w=800, h=800)
plot(ts_exp$Y[,1], ts_exp$Y[,2], cex=0.5,pch=19, col=mc@colors[as.integer(mc@mc[rownames(an[!duplicated(an),])])])
dev.off()

##### --- Facs MC plots ---- ########
md = mat@cell_metadata[names(mc@mc),]
plot(sort(log10(md[,"TIM3.BV421.A_Ab"])))
plot(sort(log10(md[,"PD1.APC.A_Ab"])))
# md[is.na(md[,"PD1.APC.A_Ab"]),"PD1.APC.A_Ab"] = 0
# md[is.na(md[,"TIM3.BV421.A_Ab"]),"TIM3.BV421.A_Ab"] = 0
# md[md[,"TIM3.BV421.A_Ab"]<200,"TIM3.BV421.A_Ab"] = 0
# md[md[,"PD1.APC.A_Ab"]<200,"PD1.APC.A_Ab"] = 0
mc_tim3_50 = tapply(md[names(mc@mc),"TIM3.BV421.A_Ab"], mc@mc, median, na.rm=T)
mc_tim3_10 = tapply(md[names(mc@mc),"TIM3.BV421.A_Ab"], mc@mc, quantile, 0.1, na.rm=T)
mc_tim3_90 = tapply(md[names(mc@mc),"TIM3.BV421.A_Ab"], mc@mc, quantile, 0.9, na.rm=T)
mc_pd1_50 = tapply(md[names(mc@mc),"PD1.APC.A_Ab"], mc@mc, median, na.rm=T)
mc_pd1_10 = tapply(md[names(mc@mc),"PD1.APC.A_Ab"], mc@mc, quantile, 0.1, na.rm=T)
mc_pd1_90 = tapply(md[names(mc@mc),"PD1.APC.A_Ab"], mc@mc, quantile, 0.9, na.rm=T)
mc_tim3_50 = tapply(md[names(mc@mc),"TIM3_Ab"], mc@mc, median, na.rm=T)
mc_tim3_10 = tapply(md[names(mc@mc),"TIM3_Ab"], mc@mc, quantile, 0.1, na.rm=T)
mc_tim3_90 = tapply(md[names(mc@mc),"TIM3_Ab"], mc@mc, quantile, 0.9, na.rm=T)
mc_pd1_50 = tapply(md[names(mc@mc),"PD1_Ab"], mc@mc, median, na.rm=T)
mc_pd1_10 = tapply(md[names(mc@mc),"PD1_Ab"], mc@mc, quantile, 0.1, na.rm=T)
mc_pd1_90 = tapply(md[names(mc@mc),"PD1_Ab"], mc@mc, quantile, 0.9, na.rm=T)
mc_ccr7_50 = tapply(md[names(mc@mc),"CCR7_Ab"], mc@mc, median, na.rm=T)
mc_ccr7_10 = tapply(md[names(mc@mc),"CCR7_Ab"], mc@mc, quantile, 0.1, na.rm=T)
mc_ccr7_90 = tapply(md[names(mc@mc),"CCR7_Ab"], mc@mc, quantile, 0.9, na.rm=T)
mc_cd25_50 = tapply(md[names(mc@mc),"CD25_Ab"], mc@mc, median, na.rm=T)
mc_tigit_50 = tapply(md[names(mc@mc),"Tigit_Ab"], mc@mc, median, na.rm=T)
mc_sell_50 = tapply(md[names(mc@mc),"Sell_Ab"], mc@mc, median, na.rm=T)
mc_il7r_50 = tapply(md[names(mc@mc),"IL7R_Ab"], mc@mc, median, na.rm=T)

# xnew = ifelse(log10(mc_tim3_50)<3, log10(mc_tim3_50), log10(mc_tim3_50)/3)
# ynew = ifelse(log10(mc_pd1_50)<3, log10(mc_pd1_50), log10(mc_pd1_50)/3)
png("figs/tim3_pd1_facs.png", w=400,h=400)
plot(log10(mc_tim3_50), log10(mc_pd1_50), pch=19, cex=2, col=mc@colors)
text(log10(mc_tim3_50), log10(mc_pd1_50), 1:length(log10(mc_pd1_50)),cex=0.5)
dev.off()

png("figs/ccr7_pd1_facs.png", w=400,h=400)
plot(log10(mc_ccr7_50), log10(mc_pd1_50), pch=19, cex=2, col=mc@colors)
text(log10(mc_ccr7_50), log10(mc_pd1_50), 1:length(log10(mc_pd1_50)),cex=0.5)
dev.off()

png("figs/tigit_cd25_facs.png", w=400,h=400)
plot(log10(mc_tigit_50), log10(mc_cd25_50), pch=19, cex=2, col=mc@colors)
text(log10(mc_tigit_50), log10(mc_cd25_50), 1:length(log10(mc_pd1_50)),cex=0.5)
dev.off()

png("figs/tim3_cd25_facs.png", w=400,h=400)
plot(log10(mc_tim3_50), log10(mc_cd25_50), pch=19, cex=2, col=mc@colors)
text(log10(mc_tim3_50), log10(mc_cd25_50), 1:length(log10(mc_pd1_50)),cex=0.5)
dev.off()

png("figs/il7r_sell_facs.png", w=400,h=400)
plot(log10(mc_il7r_50), log10(mc_sell_50), pch=19, cex=2, col=mc@colors)
text(log10(mc_il7r_50), log10(mc_sell_50), 1:length(log10(mc_pd1_50)),cex=0.5)
dev.off()

plot(log10(mc_il7r_50), log10(mc_sell_50), pch=19, cex=2, col=mc@colors)
text(log10(mc_il7r_50), log10(mc_sell_50), 1:length(log10(mc_pd1_50)),cex=0.5)

png("figs/il2ra_cd25_facs.png", w=400,h=400)
plot(lfp["Il2ra",], log10(mc_cd25_50), pch=19, cex=2, col=mc@colors)
text(log10(mc_tim3_50), log10(mc_cd25_50), 1:length(log10(mc_pd1_50)),cex=0.5)
dev.off()

plot(log10(mc_il7r_50), log10(mc_pd1_50), pch=19, cex=2, col=mc@colors)
plot(log10(mc_tim3_90), log10(mc_pd1_90), pch=19, cex=2, col=mc@colors)

tail(sort(apply(lfp, 1, cor, mc_pd1_50)),50)
tail(sort(apply(lfp, 1, cor, mc_tim3_50)),50)

png("figs/tim3_havcr2_facs.png", w=800,h=800)
plot(log10(mc_tim3_50), lfp["Havcr2",], pch=19, cex=2, col=mc@colors)
text(log10(mc_tim3_50), lfp["Havcr2",], 1:length(lfp["Havcr2",]),cex=0.5)
dev.off()

png("figs/pdcd1_pd1_facs.png", w=800,h=800)
plot(log10(mc_pd1_50), lfp["Pdcd1",], pch=19, cex=2, col=mc@colors)
text(log10(mc_pd1_50), lfp["Pdcd1",], 1:length(lfp["Pdcd1",]),cex=0.5)
dev.off()

png("figs/tim3_pdcd1_facs.png", w=400,h=400)
plot(log10(mc_tim3_50), lfp["Pdcd1",], pch=19, cex=2, col=mc@colors)
text(log10(mc_tim3_50), lfp["Pdcd1",], 1:length(log10(mc_pd1_50)),cex=0.5)
dev.off()

#####----mm10dysf NMF testing ----######
umap_obj = umap_py$UMAP(n_neighbors = as.integer(10), min_dist=0.5, spread=1) #, metric="manhattan"
mc_id = mc_id_2
mc = scdb_mc(mc_id)
rel_mcs = (1:max(mc@mc))[!mc@colors %in% c("brown", 'red')]
rel_cells = names(mc@mc)[mc@mc %in% rel_mcs]
rel_mc_cols = mc@colors[!mc@colors %in% c("brown", 'red')]
rel_col_table = table(mc@colors)[!names(table(mc@colors)) %in% c("brown", 'red')]
#rel_main_cols = main_cols[main_cols != "brown"]
egc = mc@e_gc
# knn_mat = matrix(rep(0, ncol(egc)**2), ncol=ncol(egc), nrow=ncol(egc))
# mgraph_edges = mc2d_comp_mgraph(mc_id = mc_id, graph_id = id, ignore_mismatch=F)
# for (edge in 1:nrow(mgraph_edges)){
#   knn_mat[mgraph_edges$mc1[edge], mgraph_edges$mc2[edge]] = 1
# }
# filt_egc = mc@e_gc[,rel_mcs]
# knn_mat = knn_mat[rel_mcs, rel_mcs]
# filt_lfp = log2(mc@mc_fp)[,rel_mcs]
gs = scdb_gset(mm_filt_id)
genes = names(gs@gene_set)
#genes = unique(c(genes, all_genes))
# max_lfps = apply(lfp, 1, max)
# genes = names(max_lfps)[max_lfps > 0.8] 
lat_gs = scdb_gset(mm_lateral_gset_id)
genes = genes[!genes %in% names(lat_gs@gene_set)]

#genes = unique(c(core_genes[core_genes_ref %in% non_lat_gp], genes))
#genes = unique(c(genes,  "Fos", "Jun"))

genes = genes[genes %in% rownames(egc)]
mc_2_feat_genes = mc_2$var_names[mc_2$var$top_feature_gene]
mc_2_lat_genes = mc_2$var_names[mc_2$var$forbidden_gene]

all_lat_genes = unique(c(mc_2_lat_genes, names(lat_gs@gene_set),  c("B2m", "Ldha", "Srgn", "Ms4a4b", "Itgb1", "Ly6a", "Trbc1", "Trac", "Trbc2", "Anxa2", "Sparc", "Fth1", "Fth2", "Ftl1", "Ctsd", "Lgals3", "Ccnd2", "Chd4",
                                                                    "S100a11", "Nkg7", "Epsti1", "AC153498.1", "Hsp90ab1", "Npm1", "Lgals1", "Serbp1", "Ybx1", "Tuba1b", "S100a6", "S100a4", "Ighm", "Dct", "Cd74", "Mlana", "Pmel", "Tyrobp", "Timp2", "AC117232.5", "Mir7067", "Hnrnpab", "Ran", "Ptma", "Mt1", "Mt2")))
write(all_lat_genes, file = "metacells_vignette/cd8_all_lat.txt")
write.csv(all_lat_genes, "S4-Forbidden lateral genes.csv")
write(all_lat_genes, "S4-Forbidden lateral genes.txt")

write.csv(mc@colors, "Sx-mc_colors.csv")

genes = genes[!genes %in% all_lat_genes]
gene_means_2 = apply(mat_d[mc_2_feat_genes,], 1, mean)
mc_2_feat_genes_unique = mc_2_feat_genes[!mc_2_feat_genes %in% genes]
low_mc_2_feat_genes = mc_2_feat_genes[gene_means_2 < 0.55]
genes = unique(c(genes, low_mc_2_feat_genes))

#genes = core_genes[core_genes_ref %in% non_lat_gp]
#write.csv(as.data.frame(genes), paste("figs","/nmf_plots/umap_full_progs/", "atlas_genes.csv",sep=""))
#at_genes = read.csv(paste("figs","/nmf_plots/umap_full_progs/", "atlas_genes.csv",sep=""))
gs_mgraph = gs
gs_mgraph@gene_set = gs_mgraph@gene_set[genes]
scdb_add_gset("flows_gs", gs_mgraph)


mat = scdb_mat(mm_filt_id)
mat_d = mat@mat#scm_downsamp(mat@mat, 800)
#mat_d = mat_d[,colnames(mat_d) %in% rel_cells]
# egc = 10000 * egc
# egc = log2(1+egc)
gene_means = apply(mat_d[genes,], 1, mean)
tail(sort(gene_means), 40)
scfigs_dir = .scfigs_base
main_gset = scdb_gset(g_id)
# main_cols = c("#cee6b9", "#5fb8f4", "yellow", "brown", "#e2abd6", "#1f78b4", "orange", "#975e9f", "#6a8255")
# names(main_cols) = names(table(core_genes_ref_gene))

# Preparing main gset score mat from main topics
# main_top_sum = tgs_matrix_tapply(t(mat@mat[names(main_gset@gene_set),names(mc@mc)]), main_gset@gene_set, sum)
cell_depths = colSums(mat@mat)
# # main_top_percentage = t(t(main_top_sum) / cell_depths[names(mc@mc)])
# main_top_mat = matrix(rep(0, length(main_cols)*length(genes)), nrow=length(main_cols), ncol=length(genes))
# colnames(main_top_mat) = genes
# rownames(main_top_mat) = names(main_cols)
# for (gs in 1:length(main_cols)){
#   gs_names = names(main_gset@gene_set)[main_gset@gene_set==gs]
#   gs_names = gs_names[gs_names %in% genes]
#   if (length(gs_names) > 0){
#     main_top_mat[gs,gs_names] = 1 / length(gs_names)  
#   }
#   
# }
# main_top_mat = t(main_top_mat)
# 
# mc_sums = vector()
# for (mc_ind in 1:max(mc@mc)){
#   mc_sum = sum(mat@mat[,names(mc@mc[mc@mc==mc_ind])])
#   mc_sums = c(mc_sums, mc_sum)
# }
# names(mc_sums) = 1:max(mc@mc)
# reg_const = 0.3 / quantile(mc_sums, c(0.05))

importlib = import("importlib")
umap_py = import("umap")
source("~/src/mc_local_nmf_gene_programs.R")
source("~/src/mc_local_nmf_gene_programs_by_mcs.R")
umap_obj = umap_py$UMAP(n_neighbors = as.integer(30), min_dist=0.5, spread=1) #, metric="manhattan"

rl_nmf()
graph_i = id
mc2d_i = mc_id
max_mc = max(mc@mc)
mat_i = id

lfp = log2(mc@mc_fp)
egc = log2(mc@e_gc+1e-5)

downs_mat_all = scm_downsamp(mat@mat[,names(mc@mc)], 750)


g_sd = apply(egc,1,sd)
g_mean = apply(egc,1,mean)
g_norm_sd = g_sd / g_mean
g_norm_sd = sort(-1*g_norm_sd)
g_norm_sd = g_norm_sd[!names(g_norm_sd) %in% all_lat_genes]

g_hits = names(which(g_norm_sd>0.06))
tf_norm_sd = g_norm_sd[tf_list_all_non_lat]
tf_norm_sd = tf_norm_sd[!is.na(tf_norm_sd)]
tf_norm_sd = sort(tf_norm_sd)
tf_exp = names(tf_norm_sd[tf_norm_sd > 0.03])
ind = 1
for (tf in tf_exp){
  mcell_mc_plot_gg(mc_id, "Tcf7", tf, fig_fn = paste("figs/tf_gg_plots/",ind, "_", tf, "_vs_Tcf7.png", sep=""))
  ind = ind + 1
}
interesting_tfs_mem = c("Mta1", "Ubtf", "Mbd3", "Cul1", "Gtf2h5", "Pfdn1", "Rbbp4", "Pcbp1", "Hnrnpd", "Gtf2a2", "Ash2l", "Rest", "Ddb1", "Ikzf1", "Suz12", "Sf3a2", "Tsc22d4", "Hnrnpu", "Rbbp7", "Ssrp1", "Sfpq", "Trim28", "Taf9", "Ccnh", "Supt16", "Parp1", "Mtf2", "Nrip1", "Bmi1", "Phf10", "Nfatc1", "Dnmt3b", "Zmym1", "Zfp367", "Msh6", "Ptch1", "Nmral1", "Zfpm1", "Ezh2", "Pnrc1", "Ung", "Fen1", "Eomes", "Pou2af1", "Ikzf2", "Myb", "Nr4a2")
interesting_tfs_eff_cyt = c("Zc3h15", "Rnf4", "Rela", "Polr2e", "Tbl1x", "Kras", "Snd1", "Hcls1", "Chd1", "Twistnb", "Bclaf1", "Prdx3", "Chp1", "Traf1", "Fubp1", "Top1", "Psmc5", "Cnbp", "Tardbp", "Psmd11", "Dnajc2", "Mtdh", "Eif3c", "Eif2ak2", "Maf1", "Nfatc3", "Cebpz", "Zfp330", "Hnrnpdl", "Jak2", "Usp39", "Fosl1", "Banf1", "Erh", "Cbx4", "Tet2", "E2f5","Zbtb21", "Psmd12", "Atf4", "Dnajc21", "Larp1", "Mllt3", "Creb1", "Dnmt3a", "Phb2", "Stat5a", "Utp6", "U2af1", "Btaf1", "Dnaja3", "Gtf2f2", "Gtf2e2", "Nfkb1", "Prkcq", "Tcerg1", "Farsb", "Nufip1", "Polr2l", "Phb", "Morf4l2", "Larp4", "Pum3", "Aatf", "Pwp1", "Noct", "Zfp593", "Nfkbib", "Akna", "Nfat5", "Cd40lg", "Noc4l", "Mafg", "Dot1l", "Ruvbl1", "Rpa2","Noc3l", "Lgals9", "Pml", "Ruvbl2", "Apex1", "Aes", "Pprc1", "Hp1bp3", "Polr1b", "Rbl2", "Pa2g4", "Rel", "Nop2", "Atf7ip", "Smyd5", "Tnfsf11", "Foxp1", "Atf3", "Sema4a", "Utf1", "Zeb2", "Egr2", "Tnf", "Pou2f2", "Nfkbid", "Rgcc", "Myc", "Eea1", "Irf4", "Nr4a1", "Zbtb32", "Spry1", "Hivep3", "Nr4a3", "Irf8", "Xcl1")
interesting_tfs_dysf = c("Pfdn5", "Sirt2", "Eid1", "Bloc1s1", "Ikzf3", "Sp4", "Pbx2", "Zfp318", "Irf2", "Cers4", "Cebpb", "Zfp683", "Pbxip1", "Rnf166", "Zbtb38", "Btg1", "Rora", "Spry2", "Mxd4", "Adam8", "Tox", "Id2", "Havcr2") # "Sub1"
interesting_tfs_eff_cytox = c("Thrap3", "Akt1", "Fus", "Cers2", "Smarca5", "Hdac1", "Ssbp4", "Gnptab", "Srebf2", "Smarcc1", "Gtf2h4", "Pin1", "Nfkb2", "Stat3", "Runx2", "Siva1", "Trps1", "Prdm1", "Litaf", "Tnfrsf4")
interesting_tfs_byst = c("Sp110", "Eya2", "Crebbp", "Mbnl1", "Rnf114", "Bptf", "Daxx", "Nmi", "Fli1", "Nr3c1", "Arid4b", "Rnf138", "Gata3", "Runx3", "Pnrc2", "Nlrc5", "Elf4", "Elf1", "Irf9", "Card11", "Arrb2", "Stat2", "Smad7", "Pycard", "Thap3", "Zfp36l2", "Flna", "Dtx3l", "Smad3", "Tsc22d3", "Sgk1", "Bcl11b", "Nod1", "Prex1", "Sp100", "Ssh2")
interesting_tfs_naive = c("Tcf7", "Tcf12", "Smad4", "Klf13", "Pgs1", "Pou6f1", "Hopx", "Elk4", "Txk", "Rere", "Zfp652", "Foxo1", "E2f2", "Chd3", "Rbpj", "Dtx1", "Lef1", "Aff3", "Id3", "Klf3")

ind = 1
dir.create("figs/tf_gg_plots/interesting_tfs_mem")
for (tf in interesting_tfs_mem){
  mcell_mc_plot_gg(mc_id, "Tcf7", tf, fig_fn = paste("figs/tf_gg_plots/interesting_tfs_mem/",ind, "_", tf, "_vs_Tcf7.png", sep=""))
  ind = ind + 1
}
ind = 1
dir.create("figs/tf_gg_plots/interesting_tfs_eff_cyt")
for (tf in interesting_tfs_eff_cyt){
  mcell_mc_plot_gg(mc_id, "Tcf7", tf, fig_fn = paste("figs/tf_gg_plots/interesting_tfs_eff_cyt/",ind, "_", tf, "_vs_Tcf7.png", sep=""))
  ind = ind + 1
}
ind = 1
dir.create("figs/tf_gg_plots/interesting_tfs_dysf")
for (tf in interesting_tfs_dysf){
  mcell_mc_plot_gg(mc_id, "Tcf7", tf, fig_fn = paste("figs/tf_gg_plots/interesting_tfs_dysf/",ind, "_", tf, "_vs_Tcf7.png", sep=""))
  ind = ind + 1
}
ind = 1
dir.create("figs/tf_gg_plots/interesting_tfs_eff_cytox")
for (tf in interesting_tfs_eff_cytox){
  mcell_mc_plot_gg(mc_id, "Tcf7", tf, fig_fn = paste("figs/tf_gg_plots/interesting_tfs_eff_cytox/",ind, "_", tf, "_vs_Tcf7.png", sep=""))
  ind = ind + 1
}
ind = 1
dir.create("figs/tf_gg_plots/interesting_tfs_byst")
for (tf in interesting_tfs_byst){
  mcell_mc_plot_gg(mc_id, "Tcf7", tf, fig_fn = paste("figs/tf_gg_plots/interesting_tfs_byst/",ind, "_", tf, "_vs_Tcf7.png", sep=""))
  ind = ind + 1
}
ind = 1
dir.create("figs/tf_gg_plots/interesting_tfs_naive")
for (tf in interesting_tfs_naive){
  mcell_mc_plot_gg(mc_id, "Tcf7", tf, fig_fn = paste("figs/tf_gg_plots/interesting_tfs_naive/",ind, "_", tf, "_vs_Tcf7.png", sep=""))
  ind = ind + 1
}

genes2plt_v2 = c("Nfkbid", "Myc", "Xcl1", "Ccl3", "Ccl4", "Utf1", "Irf4", "Irf8", "Zbtb32","Tnfrsf9", "Il2ra", "Prf1", "Gzmb", "Stat3", "Prdm1", "Il2rb", "Havcr2", "Pdcd1", "Ctla4", "Stat1", "Id2", "Ccl5", "Cxcr6", "Klrd1", "Tigit", "Cd244", "Lag3", "Tox", "Myb", "Zfpm1", "Ikzf1", "Nfatc1", "Tcf7", "Id3", "Sell","Klf2","S1pr1", "Pcna", "Top2a", "Hells") # "Nr4a3", "Nr4a1", "Egr2", 


col2group = names(leaders_col)
names(col2group) = leaders_col
names(col2group)[names(col2group) == "steelblue3"] = "steelblue4"

col2group = c(col2group, "cytokine_burst")
names(col2group)[length(col2group)] = "deeppink"


my_mel_plot_e_gc_barplots(mc_id, "mc_v2", genes=genes2plt_v2, egc=t(burst_genes_frac),ord_first_by_color=T, col2group = col2group, ord_by_id = c("cytokine_burst", names(leaders_col)[c(6,3,5,10,1,11,2,9,8,4,7)]))

all_interesting_tfs = c(interesting_tfs_mem, interesting_tfs_eff_cyt, interesting_tfs_dysf, interesting_tfs_eff_cytox, interesting_tfs_byst, interesting_tfs_naive)

#genes = unique(c(genes, all_interesting_tfs))
sum(!g_hits %in% genes)
sum(g_hits %in% genes)

genes = genes[!genes %in% c("Nrgn", "Igkc", "Serpinb6b", "Gas5", "Dapl1", "Jak1")]
tail(sort(gene_means[genes]), 40)

set.seed(1009)
#####----Gene topics from MCs NMFs parameters ----######
thresh_factor = 0.95 # 0.95
top_n_full = 11 # 17
split_elipse_f = 0.5 # split_elipse_f
max_genes_per_prog = 250
#filt_prog = c(12, 14 ,21)
filt_prog = c()
#####----Running MetaMan----######
main_prog_cols = NULL
main_gset = NULL
full_model_list = mcell_gene_prog_neighborhoods(mc_id, mat_i, NULL, mc_id, genes=genes, n_progs=top_n_full, cells=NULL, thresh_factor=thresh_factor, max_genes_per_prog = max_genes_per_prog, filt_prog = filt_prog, overlap_req=0, nmf_dir = "nmf_plots_mc2_burst/", glob_nmf_tol = 5e-7, scale_egc = 1e+3, manifold_reg = 1, knn_mat=umap_edges_2_sparse)

full_model_list = mcell_gene_prog_neighborhoods_plots(full_model_list, mc_id, mc2d_id = mc_id, genes, annot_colors=NULL, elipse_plots=T, nmf_dir = full_model_list$nmf_dir, n_words=20, umap_n = 15, umap_min_dist=0.5, umap_progs=T, split_elipse_f = split_elipse_f)

full_model_list$filt_prog = c()
#### Heatmaps of top down neighborhoods and gene programs ####
reg_full_h = full_model_list$reg_full_h
#save(reg_full_h, file=paste("scrna_db/reg_full_h_2.Rda"))
reg_full_w = full_model_list$reg_full_w
reg_full_w_nf = full_model_list$reg_full_w_non_filt

full_model_list$leaders_col[11] = "steelblue3"
full_model_list$leaders_col[2] = "#5fb8f4"
# full_model_list$leaders_col[3] = "#6a8255"
full_model_list$leaders_col[4] = "#cee6b9"
# full_model_list$leaders_col[5] = "#e2abd6"
# full_model_list$leaders_col[6] = "#5fb8f4"
# # full_model_list$leaders_col[7] = "#e2abd6"
# full_model_list$leaders_col[8] = "deeppink"
# full_model_list$leaders_col[9] = "yellow"
# full_model_list$leaders_col[15] = "magenta"
# full_model_list$leaders_col[13] = "steelblue2"
# full_model_list$leaders_col[10] = "cyan"
# full_model_list$leaders_col[14] = "deeppink"
# full_model_list$leaders_col[16] = "#1f78b4"

# full_model_list$leaders_col[16] = "#975e9f"
# full_model_list$leaders_col[6] = "gold"
# full_model_list$leaders_col[17] = "springgreen"
# full_model_list$leaders_col[10] = "blue"
# full_model_list$leaders_col[1] = "#cee6b9"
# full_model_list$leaders_col[11] = "#5fb8f4"
# full_model_list$leaders_col[22] = "cyan"
full_model_list$filt_prog = c(7)
filt_prog = c(7)
full_model_list$leaders_col[full_model_list$filt_prog] = "grey"

# names(full_model_list$leaders_col) = c("eff_cytokine_Xcl1", "byst_Ccl5", "eff_mem_progenitor_Myb", "eff_cytotox_Tnfrsf9", "eff_cytokine_Ccl3", "naive_Tcf7", "byst_eff_Gzma", "eff_cytokine_Nfkbid", "eff_cytotox_Gzmf", "dysf_Tigit")
names(full_model_list$leaders_col) = c("Klrk1_TME_activation", "naive_mem1_Tcf7", "eff_cytotox_Tnfrsf9", "byst_em__Id2", "eff_cytotox_Prf1","eff_cytokine_Myc", "filt1_Ahnak", "naive_byst_Ctla2a", "naive2_Il7r", "dysf_Tigit", "mem_progenitor_Myb")
leaders_col = full_model_list$leaders_col
colnames(reg_full_w) = names(leaders_col)
rownames(reg_full_h) = names(leaders_col)
save(leaders_col, file=paste("scrna_db/leaders_col.Rda"))

colnames(reg_full_w_nf) = names(leaders_col)

burst_genes =  c("Ccl3", "Ccl4", "Ccl5", "Xcl1","Ccl1", "Ccl9", "Gzmb", "Gzma", "Gzmf", "Gzmc", "Gzmd", "Il2", "Il3", "Spp1") #, "Nrgn"
write.csv(burst_genes, file="S6-Excluded bursty genes")

burst_genes_mc_sum = tgs_matrix_tapply(mat@mat[unique(c(rownames(mc@e_gc),burst_genes)),names(mc@mc)], mc@mc, sum)
burst_genes_depth = rowSums(burst_genes_mc_sum)
burst_genes_frac = burst_genes_mc_sum / burst_genes_depth
burst_genes_frac_f = burst_genes_frac[,burst_genes]

burst_cor = tgs_cor(reg_full_w_nf, burst_genes_frac_f)
pheatmap(burst_cor)

burst_gene = "Spp1"
burst_prog = 3
png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_burst_",burst_gene, "_", burst_prog, ".png",sep=""), 10 * 100, 10 * 100, res = 200)
plot(sort(reg_full_w_nf[,burst_prog]), burst_genes_frac_f[order(reg_full_w_nf[,burst_prog]),burst_gene], xlab=colnames(reg_full_w_nf)[burst_prog], ylab=burst_gene, col=mc@colors[order(reg_full_w_nf[,burst_prog])], pch=19)
dev.off()

rel_mcs = rownames(reg_full_w)
rel_cells_from_mcs = names(mc@mc)[mc@mc %in% as.integer(rel_mcs)]
downs_mat = scm_downsamp(mat@mat[,rel_cells_from_mcs], 750)
#downs_mat = cbind(downs_mat, mat@mat[,setdiff(rel_cells_from_mcs, colnames(downs_mat))])
#downs_mat = downs_mat[,rel_cells_from_mcs]
save(downs_mat, file=paste("downs_mat.Rda"))
rel_mat = as.matrix(downs_mat[colnames(reg_full_h),rel_downs_cells])
rel_mat_cor = tgs_cor(t(rel_mat))
save(rel_mat_cor, file=paste("scrna_db/reg_full_h_cells_cor.Rda"))

# rel_depths = colSums(mat@mat[,rel_cells_from_mcs])
# rel_mat = t(t(rel_mat) / rel_depths)
# quantile(rel_depths)

W_nnls_loc_downs = .fcnnls(t(reg_full_h), rel_mat)
W_loc_downs = t(W_nnls_loc_downs$coef)

W_nnls_loc_downs = nmf_py$non_negative_factorization(t(as.matrix(rel_mat)), n_components=as.integer(nrow(reg_full_h)), H = reg_full_h, update_H = F, beta_loss = 'kullback-leibler', solver='mu', alpha=0.5,l1_ratio=0.9)
W_loc_downs = W_nnls_loc_downs[[1]]
rownames(W_loc_downs) = colnames(rel_mat)

mc_prog_cors = matrix(nrow = length(mc@colors), ncol=length(burst_genes))


colnames(mc_prog_cors) = burst_genes
rownames(mc_prog_cors) = 1:length(mc@colors)

for (mc_ind in 1:length(mc@colors)){
  which_mc_cells = rel_downs_cells[mc@mc[rel_downs_cells] == as.integer(mc_ind)]
  downs_mc = downs_mat[burst_genes,which_mc_cells]
  mc_prog_cors[mc_ind,1:ncol(mc_prog_cors)] = cor(as.matrix(t(downs_mc)), W_loc_downs[which_mc_cells,8, drop=F], method = "spearman")
}

rel_downs_cells = colnames(downs_mat)
burst_gene = "Ccl5"
burst_prog = 8
which_mc = rownames(reg_full_w_nf)[which.max(burst_genes_frac_f[,burst_gene])]
which_mc = "682"
which_mc_cells = rel_downs_cells[mc@mc[rel_downs_cells] == as.integer(which_mc)]
png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/", which_mc,"_prog_vs_burst_",burst_gene, "_", burst_prog, ".png",sep=""), 10 * 100, 10 * 100, res = 200)
plot(W_loc_downs[which_mc_cells,burst_prog], downs_mat[burst_gene,which_mc_cells], xlab=colnames(reg_full_w_nf)[burst_prog], ylab=burst_gene, pch=19, main=which_mc)
dev.off()

W_loc_mcs_mean = tgs_matrix_tapply(t(W_loc_downs), mc@mc[rel_cells_from_mcs], mean)
W_loc_mcs_var = tgs_matrix_tapply(t(W_loc_downs), mc@mc[rel_cells_from_mcs], var)

# downs_mat = scm_downsamp(mat@mat[,names(mc@mc)], 750)
# rel_cells_from_mcs = colnames(downs_mat)
W_loc_mcs_mean_mc = tgs_matrix_tapply(downs_mat[colnames(local_h),rel_cells_from_mcs], mc@mc[rel_cells_from_mcs], mean)
W_loc_mcs_var_mc = tgs_matrix_tapply(downs_mat[colnames(local_h),rel_cells_from_mcs], mc@mc[rel_cells_from_mcs], var)

W_loc_mcs_mean_mc = W_loc_mcs_mean_mc[table(mc@mc[rel_cells_from_mcs]) > 40,]
W_loc_mcs_var_mc = W_loc_mcs_var_mc[table(mc@mc[rel_cells_from_mcs]) > 40,]
mean_feat_genes = colMeans(W_loc_mcs_mean_mc)
high_mean_feat_genes = names(mean_feat_genes)[mean_feat_genes > 0.1]

varmean = W_loc_mcs_var_mc/W_loc_mcs_mean_mc
varmean[!is.finite(varmean)] = 0
varmean = varmean[,high_mean_feat_genes]
mean_W_loc_mcs_varmean_mc = apply(varmean, 2, quantile,1)
tail(sort(mean_W_loc_mcs_varmean_mc), 100)

mean_W_loc_mcs_varmean_which_mc = apply(varmean, 2, which.max)
max_varmean_mcs = rownames(varmean)[mean_W_loc_mcs_varmean_which_mc]
names(max_varmean_mcs) = names(mean_W_loc_mcs_varmean_which_mc)

which_gene = "Xcl1"
which_mc = max_varmean_mcs[which_gene]
which_mc_cells = rel_cells_from_mcs[mc@mc[rel_cells_from_mcs] == as.integer(which_mc)]
which_mc_genes = names(tail(sort(egc[genes[genes %in% rownames(downs_mat)],which_mc]), 50))
which_mc_genes = setdiff(which_mc_genes, which_gene)
which_mc_score = apply(downs_mat[which_mc_genes,which_mc_cells], 2, sum)
plot(which_mc_score, downs_mat[which_gene,which_mc_cells], ylab=which_gene, xlab="Score", main=which_mc)


W_loc_mcs_mean_mc_high_mean = W_loc_mcs_mean_mc[,high_mean_feat_genes]
W_loc_mcs_var_mc_high_mean = W_loc_mcs_var_mc[,high_mean_feat_genes]
varmean = W_loc_mcs_var_mc/W_loc_mcs_mean_mc
varmean[!is.finite(varmean)] = 0
mean_W_loc_mcs_varmean_mc = apply(varmean, 2, quantile, 0.95)
tail(sort(mean_W_loc_mcs_varmean_mc), 100)
high_varmean_genes = names(mean_W_loc_mcs_varmean_mc)[mean_W_loc_mcs_varmean_mc > 1.7]


overlap_list = full_model_list$overlap_list
#full_model_list$leaders_col[12] = "yellow"

genes_not_shown = c()
#prog_ord =  c(5,17,2,14,13,10,6,12,3,4,11,18,15,8,16,1,7,9) # 1:ncol(mc_usage)#c(3,5,1,4,7,2,6) (1,12,13,14,16,18,23,25,26,27)
prog_ord =  1:top_n_full

##### UMAP testing #####
filt_mcs = (1:max(mc@mc))[mc@colors %in% c("brown", 'red')]
dend_prog=""
umap_py = import("umap")
norm_global_w = reg_full_w / rowSums(reg_full_w)
glob_reconst_egc = norm_global_w %*% reg_full_h

norm_global_mat = full_model_list$reg_full_w_non_filt
norm_global_mat = norm_global_mat / rowSums(norm_global_mat)
reconst_egc_full = norm_global_mat %*% reg_full_h

dist_w = as.matrix(dist(norm_global_mat))
if (length(filt_prog)>0){
  dist_w = as.matrix(dist(norm_global_mat[,-filt_prog]))
}

is_spec_atlas = md$cell_type == "spec"
is_tumor_atlas =md$location == "tumor"
is_day6_atlas = md$days_post_transfer == 6
is_atlas_compare = is_spec_atlas & is_tumor_atlas & is_day6_atlas & is_ctrl
norm_global_w_filt = norm_global_w[,-7]
norm_global_w_cells_spec = tgs_matrix_tapply(t(norm_global_w_filt[as.character(mc@mc)[is_atlas_compare],]), md$mouse_id[is_atlas_compare], mean)
norm_global_w_cells_spec[is.na(norm_global_w_cells_spec)] = 0
norm_global_w_cells_spec = norm_global_w_cells_spec[rowSums(norm_global_w_cells_spec) > 0,]
save(norm_global_w_cells_spec, file="norm_global_w_cells_spec")

sigma_k = 0.5
sim_w = exp((-dist_w^2)/sigma_k^2)
dist_w_clean = max(sim_w) - sim_w
diag(dist_w_clean) = 1
write.csv(x=sim_w, file=paste("scrna_db/similarity_mat_mcs_cd8.csv"))
mc_cols = mc@colors
names(mc_cols) = rownames(dist_w_clean)
write.csv(x=mc_cols, file=paste("scrna_db/colors_mcs_cd8.csv"))

# norm_global_mat[filt_mcs, -dend_prog] = 0
# norm_global_w[filt_mcs, -dend_prog] = 0
norm_full_h = t(reg_full_h) / colSums(reg_full_h)
colnames(norm_full_h) = names(leaders_col)
genes_hm = c()
for (dd in 1:ncol(norm_full_h)){
  genes_hm = c(genes_hm, names(tail(sort(norm_full_h[,dd]), 20)))
  # if (dd == dend_prog){
  #   dend_genes_norm = names(tail(sort(norm_full_h[,dd]), 20))
  # }
}
genes_hm = unique(genes_hm)

mcell_gene_prog_neigborhoods_heatmaps(full_model_list$reg_full_w, full_model_list$reg_full_h, full_model_list$leaders_col, genes_not_shown=c(), main_prog_cols=main_prog_cols, nmf_dir="nmf_plots_mc2_burst/", burn_quant = 0.995, zero_factor = 5, main_gset=main_gset, use_slanter=T, genes2show = genes_hm)# 

# dend_genes_norm = c()
# genes_hm = genes_hm[!genes_hm %in% dend_genes_norm]

#### set up metadata booleans ####
mat@cell_metadata[mat@cell_metadata$amp_batch_id == "AB10697","mouse_id"] = "m224"
rel_cells = names(mc@mc)
md = mc_md
md = mat@cell_metadata[rel_cells,]



times = as.numeric(as.character(md$days_post_transfer))
is_spec = md$cell_type == "spec"
is_endo = md$cell_type == "endo"
is_byst = md$cell_type == "byst"
is_tumor = md$location == "tumor"
is_spleen= md$location == "spleen"
is_ln= md$location == "LN"

is_mc38 = md$batch_set_id %in% c('MC38_41BB+PD1_2', 'MC38_41BB_2', 'MC38_41BB+PD1_3', 'MC38_Ctrl_2', 'MC38_PD1_3', 'MC38_PD1_2', 'MC38_41BB_3', 'MC38_Ctrl_3')
names(is_mc38) = names(mc@mc)
is_b16_endo = is_endo & !is_mc38
names(is_b16_endo) = names(mc@mc)
is_pd1 = md$batch_set_id == "PD1_5" | md$batch_set_id == "PD1_treatment" | md$batch_set_id == "MC38_PD1_3" | md$batch_set_id == "MC38_PD1_2" | md$batch_set_id == "S6_Panel_PD1_1" | md$batch_set_id == "pd1_8" | md$batch_set_id == "pd1_9"
names(is_pd1) = names(mc@mc)
is_41bb = md$batch_set_id == "MC38_41BB_3" | md$batch_set_id == "MC38_41BB_2" | md$batch_set_id == "41bb_8" | md$batch_set_id == "41bb_9"  
names(is_41bb) = names(mc@mc)
is_41bb_pd1 = md$batch_set_id == "MC38_41BB+PD1_3" | md$batch_set_id == "41bb+pd1_8" | md$batch_set_id == "MC38_41BB+PD1_2" | md$batch_set_id == "41bb+pd1_9"
names(is_41bb_pd1) = names(mc@mc)
is_filt = md$batch_set_id %in% c('BFP_pd1', '41BB_ko_pd1', 'Delay_OT1_GFPOT1_PD1', 'Delay_only_OT1_ctrl', 'Delay_only_OT1_PD1', 'Delay_OT1_GFPOT1_ctrl', 'PD1_ko_pd1', 'PD1_ko_ctrl', 'Zbdb32_ko_pd1', 'ID3_KO', 'Zbdb32_ko_Ctrl', 'BFP_Ctrl', 'ID3_control', '41BB_ko_Ctrl', "PD1_7", "naive_transfer", "cd8_in_vitro")

m_id = md$mouse_id
t_points = sort(unique(times))
names(is_filt) = names(mc@mc)
is_ctrl = (!is_pd1) & (!is_41bb) & (!is_41bb_pd1) & (!is_filt)
names(is_ctrl) = names(mc@mc)

md$condition = rep("na", nrow(md))
md$condition[is_pd1] = "pd1"
md$condition[is_ctrl] = "ctrl"
md$condition[is_41bb] = "41bb"
md$condition[is_41bb_pd1] = "41bb_pd1"

is_ko_wt = md$cell_type == "ctrl_ko_spec"
is_ko_pd1 = md$cell_type == "pd1_ko_spec"
is_ko_ctrl = md$batch_set_id == "PD1_ko_ctrl"
is_ko_treatment = md$batch_set_id == "PD1_ko_pd1"

levels(md$cell_type) = c(levels(md$cell_type), "endo_mc38")
md$cell_type[is_mc38] = "endo_mc38"

#### cal umap globally ####
umap_full = full_model_list$umap_full
png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/", "umap_prog_full2.png",sep=""), 12 * 100, 12 * 100, res = 200)
par(mar=c(0,0,0,0))
plot(umap_full[,1], umap_full[,2], bg=mc@colors, pch = 21, col="black")
dev.off()


mc2d = scdb_mc2d(mc_id)
mc2d@mc_x = umap_full[,1]
mc2d@mc_y = umap_full[,2]
mc_xy = data.frame('mc_x' = mc2d@mc_x, 'mc_y' = mc2d@mc_y)
mgraph = mc2d_comp_mgraph(mc_id = mc_id, graph_id = id, ignore_mismatch=F)
xy = mc2d_comp_cell_coord(mc_id, id, mgraph, mc_xy, symmetrize=F)
mc2d@sc_x = xy$x
mc2d@sc_y = xy$y
scdb_add_mc2d("atlas_gene_prog_umap", mc2d)

mc_age = tapply(as.integer(as.character(mat@cell_metadata[names(mc@mc),'days_post_transfer'])), mc@mc, mean)
mc_annot = data.frame('metacell' = (1:max(mc@mc)), cell_type = col2group[mc@colors], mc_col = mc@colors, mc_age = mc_age)
write.table(file = "scrna_db/mc_metadata.tsv", x = mc_annot)
cell_type_annot = data.frame('cell_type' = names(group2col), 'col' = group2col, 'ord' = 1:length(group2col))
write.table(file = "scrna_db/cell_type_annot.tsv", x = cell_type_annot)


.plot_start(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/", "umap_mc_2.svg",sep=""), 12 * 100, 12 * 100)
par(mar=c(0,0,0,0))
plot(mc2d_2@mc_x, mc2d_2@mc_y, bg=mc@colors, pch = 21, col=mc@colors)
dev.off()


colspec = get_param("mcell_mc2d_gene_shades", package = "metacell")
colspec = c("white", colspec[1], colspec[11])
cc_sig = apply(lfp[cc_genes,],2, mean)
prog_usage = cc_sig
max_usage = max(prog_usage)
min_usage = min(prog_usage)
x = pmin(pmax(prog_usage, min_usage), max_usage) - min_usage
shades = colorRampPalette(colspec)(100 * (max_usage-min_usage) + 1)
mc_cols = shades[round(100 * x) + 1]
png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/", "umap_prog_full_cc.png",sep=""), 12 * 100, 12 * 100, res = 200)
par(mar=c(0,0,0,0))
plot(mc2d_2@mc_x, mc2d_2@mc_y, bg=mc_cols, pch = 21, col="black")
dev.off()
prog_usage = cc_g1_sig
max_usage = max(prog_usage)
min_usage = min(prog_usage)
x = pmin(pmax(prog_usage, min_usage), max_usage) - min_usage
shades = colorRampPalette(colspec)(100 * (max_usage-min_usage) + 1)
mc_cols = shades[round(100 * x) + 1]
png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/", "umap_prog_full_cc_g1.png",sep=""), 12 * 100, 12 * 100, res = 200)
par(mar=c(0,0,0,0))
plot(umap_full[,1], umap_full[,2], bg=mc_cols, pch = 21, col="black")
dev.off()
prog_usage = cc_g2_sig
max_usage = max(prog_usage)
min_usage = min(prog_usage)
x = pmin(pmax(prog_usage, min_usage), max_usage) - min_usage
shades = colorRampPalette(colspec)(100 * (max_usage-min_usage) + 1)
mc_cols = shades[round(100 * x) + 1]
png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/", "umap_prog_full_cc_g2.png",sep=""), 12 * 100, 12 * 100, res = 200)
par(mar=c(0,0,0,0))
plot(umap_full[,1], umap_full[,2], bg=mc_cols, pch = 21, col="black")
dev.off()
prog_usage = reconst_egc_full[,"Pdcd1"]
max_usage = max(prog_usage)
min_usage = min(prog_usage)
x = pmin(pmax(prog_usage, min_usage), max_usage) - min_usage
shades = colorRampPalette(colspec)(100 * (max_usage-min_usage) + 1)
mc_cols = shades[round(100 * x) + 1]
png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/", "umap_prog_full_pdcd1.png",sep=""), 12 * 100, 12 * 100, res = 200)
par(mar=c(0,0,0,0))
plot(mc2d_2@mc_x, mc2d_2@mc_y, bg=mc_cols, pch = 21, col="black")
dev.off()
prog_usage = reconst_egc_full[,"Tnfrsf9"]
max_usage = max(prog_usage)
min_usage = min(prog_usage)
x = pmin(pmax(prog_usage, min_usage), max_usage) - min_usage
shades = colorRampPalette(colspec)(100 * (max_usage-min_usage) + 1)
mc_cols = shades[round(100 * x) + 1]
png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/", "umap_prog_full_tnfrsf9.png",sep=""), 12 * 100, 12 * 100, res = 200)
par(mar=c(0,0,0,0))
plot(mc2d_2@mc_x, mc2d_2@mc_y, bg=mc_cols, pch = 21, col="black")
dev.off()

mc_locations = table(mc@mc, md[names(mc@mc), "location"])
mc_locations = mc_locations / rowSums(mc_locations)
for (loc in 1:ncol(mc_locations)){
  loc_nm = colnames(mc_locations)[loc]
  prog_usage = mc_locations[,loc]
  max_usage = max(prog_usage)
  min_usage = 0
  x = pmin(pmax(prog_usage, min_usage), max_usage) - min_usage
  shades = colorRampPalette(colspec)(100 * (max_usage-min_usage) + 1)
  mc_cols = shades[round(100 * x) + 1]
  png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/", "umap_prog_full_", loc_nm, ".png",sep=""), 12 * 100, 12 * 100, res = 200)
  par(mar=c(0,0,0,0))
  plot(mc2d_2@mc_x, mc2d_2@mc_y, bg=mc_cols, pch = 21, col="black")
  dev.off()
}


mc_locations = table(mc@mc, md[names(mc@mc), "condition"])
mc_locations = mc_locations / rowSums(mc_locations)
for (loc in 1:ncol(mc_locations)){
  loc_nm = colnames(mc_locations)[loc]
  prog_usage = mc_locations[,loc]
  max_usage = max(prog_usage)
  min_usage = 0
  x = pmin(pmax(prog_usage, min_usage), max_usage) - min_usage
  shades = colorRampPalette(colspec)(100 * (max_usage-min_usage) + 1)
  mc_cols = shades[round(100 * x) + 1]
  png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/", "umap_prog_full_", loc_nm, ".png",sep=""), 12 * 100, 12 * 100, res = 200)
  par(mar=c(0,0,0,0))
  plot(mc2d_2@mc_x, mc2d_2@mc_y, bg=mc_cols, pch = 21, col="black")
  dev.off()
}


md$bigbatch = rep("None", nrow(md))
md$bigbatch[is_mc38 & is_tumor] = "mc38"
md$bigbatch[is_b16_endo & is_tumor] = "b16"
md$bigbatch[is_spec & is_tumor] = "spec"
mc_batch = table(mc@mc, md[names(mc@mc), "bigbatch"])
#mc_batch = mc_batch[,-3]
mc_batch = mc_batch / rowSums(mc_batch)
for (loc in 1:ncol(mc_batch)){
  loc_nm = colnames(mc_batch)[loc]
  prog_usage = mc_batch[,loc]
  max_usage = 0.5*max(prog_usage)
  min_usage = min(prog_usage)
  x = pmin(pmax(prog_usage, min_usage), max_usage) - min_usage
  shades = colorRampPalette(colspec)(100 * (max_usage-min_usage) + 1)
  mc_cols = shades[round(100 * x) + 1]
  png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/", "umap_prog_full_", loc_nm, ".png",sep=""), 12 * 100, 12 * 100, res = 200)
  par(mar=c(0,0,0,0))
  plot(umap_full[,1], umap_full[,2], bg=mc_cols, pch = 21, col="black")
  dev.off()
}

# mc_hc = confs[[id]]$mc_hc
# sup_clust = rep(0, nrow(reg_full_w))
# names(sup_clust) = rownames(reg_full_w)
# 
# sup_clust[sup[[722]]$mcs] = 1
# sup_clust[sup[[744]]$mcs] = 2
# sup_clust[sup[[836]]$mcs] = 3
# sup_clust[sup[[1009]]$mcs] = 4
# sup_clust[sup[[1413]]$mcs] = 5
# sup_clust[sup[[1444]]$mcs] = 6
# sup_clust[sup[[1452]]$mcs] = 6
# sup_clust[sup[[1150]]$mcs] = 9
# sup_clust[sup[[1286]]$mcs] = 8
# sup_clust[sup[[1295]]$mcs] = 7
# sup_clust[sup[[226]]$mcs] = 10
# sup_clust[sup[[64]]$mcs] = 10
# sup_clust[sup[[153]]$mcs] = 10

sup_clust = mc_clusts
clust_names = col2group[as.character(mc2_cols$color)]

non_filt_gset = scdb_gset(mm_filt_id)
cc_genes = names(non_filt_gset@gene_set)[non_filt_gset@gene_set %in% c(110,11,60,130)]
cc_sig = apply(lfp[cc_genes,],2, mean)
cc_sig_by_cell = cc_sig[mc@mc]
names(cc_sig_by_cell) = names(mc@mc)
clust_by_cell = sup_clust[mc@mc]
#clust_names = c("cycling_trans", "central_mem", "res_mem", "byst", "effector_early", "cytokine", "cytotoxic-dysf", "cytotoxic-grzm", "dysf", "naive")
for (cc in 1:length(clust_names)){
  c_name = clust_names[cc]
  act_mat = matrix(rep(0, 4*3), nrow=4, ncol=3)
  rownames(act_mat) = c("Ctrl", "aPD1", "a41BB", "aPD1+a41BB")
  colnames(act_mat) = c("OT1", "B16-endo", "MC38-endo")
  
  cond_cells = rel_cells[clust_by_cell == cc & is_ctrl & is_spec & is_tumor]
  act_mat["Ctrl", "OT1"] = mean(cc_sig_by_cell[cond_cells])
  cond_cells = rel_cells[clust_by_cell == cc & is_pd1 & is_spec & is_tumor]
  act_mat["aPD1", "OT1"] = mean(cc_sig_by_cell[cond_cells])
  cond_cells = rel_cells[clust_by_cell == cc & is_41bb & is_spec & is_tumor]
  act_mat["a41BB", "OT1"] = mean(cc_sig_by_cell[cond_cells])
  cond_cells = rel_cells[clust_by_cell == cc & is_41bb_pd1 & is_spec & is_tumor]
  act_mat["aPD1+a41BB", "OT1"] = mean(cc_sig_by_cell[cond_cells])
  cond_cells = rel_cells[clust_by_cell == cc & is_ctrl & is_b16_endo & is_tumor]
  act_mat["Ctrl", "B16-endo"] = mean(cc_sig_by_cell[cond_cells])
  cond_cells = rel_cells[clust_by_cell == cc & is_pd1 & is_b16_endo & is_tumor]
  act_mat["aPD1", "B16-endo"] = mean(cc_sig_by_cell[cond_cells])
  cond_cells = rel_cells[clust_by_cell == cc & is_41bb & is_b16_endo & is_tumor]
  act_mat["a41BB", "B16-endo"] = mean(cc_sig_by_cell[cond_cells])
  cond_cells = rel_cells[clust_by_cell == cc & is_41bb_pd1 & is_b16_endo & is_tumor]
  act_mat["aPD1+a41BB", "B16-endo"] = mean(cc_sig_by_cell[cond_cells])
  cond_cells = rel_cells[clust_by_cell == cc & is_ctrl & is_mc38 & is_tumor]
  act_mat["Ctrl", "MC38-endo"] = mean(cc_sig_by_cell[cond_cells])
  cond_cells = rel_cells[clust_by_cell == cc & is_pd1 & is_mc38 & is_tumor]
  act_mat["aPD1", "MC38-endo"] = mean(cc_sig_by_cell[cond_cells])
  cond_cells = rel_cells[clust_by_cell == cc & is_41bb & is_mc38 & is_tumor]
  act_mat["a41BB", "MC38-endo"] = mean(cc_sig_by_cell[cond_cells])
  cond_cells = rel_cells[clust_by_cell == cc & is_41bb_pd1 & is_mc38 & is_tumor]
  act_mat["aPD1+a41BB", "MC38-endo"] = mean(cc_sig_by_cell[cond_cells])
  
  dir.create(paste0(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_distributions/"))
  png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_distributions/", "barplot_cell_cycle_",c_name, ".png",sep=""), 12 * 100, 12 * 100, res = 200)
  barplot(act_mat, names.arg=colnames(act_mat), beside=T, col = c("grey", "black", "pink", "red"), main=paste("Cell cycle clust = ", c_name, sep=""), ylim=c(-0.3, 0.5))
  dev.off()
  
}

dir.create(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/", sep=""),showWarnings = F)
for (ii in 1:ncol(norm_global_w)){
  png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/", "no_arrows_umap_proj_prog_",ii,".png",sep=""), 50 + 2 * 12 * 100, 12 * 100, res = 200)
  proj_prog = ii
  colspec = get_param("mcell_mc2d_gene_shades", package = "metacell")
  colspec = c("white", colspec[1], colspec[11])
  min_usage = 0
  prog_usage = norm_global_w[,proj_prog]
  max_usage = max(prog_usage)
  x = pmin(pmax(prog_usage, min_usage), max_usage) - min_usage
  shades = colorRampPalette(colspec)(100 * (max_usage-min_usage) + 1)
  mc_cols = shades[round(100 * x) + 1]
  par(mfrow=c(1, 2))
  #par(mar=c(0,0,0,1))
  #plot(umap_16[,1], umap_16[,2], col=c_cols, pch = 19)
  plot(umap_full[,1], umap_full[,2], bg=mc_cols, pch = 21, col="black")
  prog_genes = sort(norm_full_h[,ii], decreasing = T)
  n_words = 20
  barplot(prog_genes[n_words:1], names.arg=names(prog_genes[n_words:1]), horiz = T,col=leaders_col[ii],las=2, cex.names = 0.7)
  #plot(umap_16[,1], umap_16[,2], col=c_cols, pch = 19)
  dev.off()
}
for (ii in 1:ncol(norm_global_w)){
  png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/", "non_norm_umap_proj_prog_",ii,".png",sep=""), 50 + 2 * 12 * 100, 12 * 100, res = 200)
  proj_prog = ii
  colspec = get_param("mcell_mc2d_gene_shades", package = "metacell")
  colspec = c("white", colspec[1], colspec[11])
  min_usage = 0
  zero_sc_v = 0
  one_sc_v = 1
  two_sc_v=2
  prog_usage = norm_global_mat[,proj_prog]
  max_usage = max(prog_usage)
  x = pmin(pmax(prog_usage, min_usage), max_usage) - min_usage
  shades = colorRampPalette(colspec)(100 * (max_usage-min_usage) + 1)
  mc_cols = shades[round(100 * x) + 1]
  par(mfrow=c(1, 2))
  #par(mar=c(0,0,0,1))
  #plot(umap_16[,1], umap_16[,2], col=c_cols, pch = 19)
  plot(umap_full[,1], umap_full[,2], bg=mc_cols, pch = 21, col="black")
  prog_genes = sort(reg_full_h[ii,], decreasing = T)
  n_words = 20
  barplot(prog_genes[n_words:1], names.arg=names(prog_genes[n_words:1]), horiz = T,col=leaders_col[ii],las=2, cex.names = 0.7)
  #plot(umap_16[,1], umap_16[,2], col=c_cols, pch = 19)
  dev.off()
}
for (ii in 1:ncol(norm_global_w)){
  png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/", "no_arrows_mc2d_2_proj_prog_",ii,".png",sep=""), 50 + 2 * 12 * 100, 12 * 100, res = 200)
  proj_prog = ii
  colspec = get_param("mcell_mc2d_gene_shades", package = "metacell")
  colspec = colspec = c("white", colspec[1], colspec[11])
  min_usage = 0
  zero_sc_v = 0
  one_sc_v = 1
  two_sc_v=2
  prog_usage = norm_global_w[,proj_prog]
  max_usage = max(prog_usage)
  x = pmin(pmax(prog_usage, min_usage), max_usage) - min_usage
  shades = colorRampPalette(colspec)(100 * (max_usage-min_usage) + 1)
  mc_cols = shades[round(100 * x) + 1]
  par(mfrow=c(1, 2))
  #par(mar=c(0,0,0,1))
  #plot(umap_16[,1], umap_16[,2], col=c_cols, pch = 19)
  plot(mc2d@mc_x, mc2d@mc_y, bg=mc_cols, pch = 21, col="black")
  prog_genes = sort(norm_full_h[,ii], decreasing = T)
  n_words = 20
  barplot(prog_genes[n_words:1], names.arg=names(prog_genes[n_words:1]), horiz = T,col=leaders_col[ii],las=2, cex.names = 0.7)
  dev.off()
}
for (ii in 1:ncol(norm_global_w)){
  .plot_start(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/", "no_bar_mc2d_2_proj_prog_",ii,".svg",sep=""),  6 * 100, 6 * 100)
  par(mar=c(0,0,0,0))
  proj_prog = ii
  colspec = get_param("mcell_mc2d_gene_shades", package = "metacell")
  colspec = colspec = c("white","red")
  min_usage = 0
  zero_sc_v = 0
  one_sc_v = 1
  two_sc_v=2
  prog_usage = norm_global_w[,proj_prog]
  max_usage = max(prog_usage)
  x = pmin(pmax(prog_usage, min_usage), max_usage) - min_usage
  shades = colorRampPalette(colspec)(10 * (max_usage-min_usage) + 1)
  mc_cols = shades[round(10 * x) + 1]
  #par(mfrow=c(1, 2))
  #par(mar=c(0,0,0,1))
  #plot(umap_16[,1], umap_16[,2], col=c_cols, pch = 19)
  plot(mc2d_2@mc_x, mc2d_2@mc_y, bg="white", pch = 21, col="black", cex=2)
  points(mc2d_2@mc_x[x > 0], mc2d_2@mc_y[x > 0], bg=mc_cols[x > 0], pch = 21, col="black", cex=2.5)
  # prog_genes = sort(norm_full_h[,ii], decreasing = T)
  # n_words = 20
  # barplot(prog_genes[n_words:1], names.arg=names(prog_genes[n_words:1]), horiz = T,col=leaders_col[ii],las=2, cex.names = 0.7)
  dev.off()
}

for (ii in 1:ncol(norm_global_w)){
  png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/", "non_norm_mc2d_proj_prog_",ii,".png",sep=""), 50 + 2 * 12 * 100, 12 * 100, res = 200)
  proj_prog = ii
  colspec = get_param("mcell_mc2d_gene_shades", package = "metacell")
  colspec = colspec = c("white", colspec[1], colspec[11])
  min_usage = 0
  zero_sc_v = 0
  one_sc_v = 1
  two_sc_v=2
  prog_usage = norm_global_w[,proj_prog]
  max_usage = max(prog_usage)
  x = pmin(pmax(prog_usage, min_usage), max_usage) - min_usage
  shades = colorRampPalette(colspec)(100 * (max_usage-min_usage) + 1)
  mc_cols = shades[round(100 * x) + 1]
  par(mfrow=c(1, 2))
  #par(mar=c(0,0,0,1))
  #plot(umap_16[,1], umap_16[,2], col=c_cols, pch = 19)
  plot(mc2d@mc_x, mc2d@mc_y, bg=mc_cols, pch = 21, col="black")
  prog_genes = sort(reg_full_h[ii,], decreasing = T)
  n_words = 20
  barplot(prog_genes[n_words:1], names.arg=names(prog_genes[n_words:1]), horiz = T,col=leaders_col[ii],las=2, cex.names = 0.7)
  dev.off()
}

for (ii in 1:ncol(norm_global_w)){
  .plot_start(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/", "bargraph_fig_",ii,".svg",sep=""), 6 * 100, 10 * 100)
  par(mar=c(5,8,4,2) + 0.1)
  proj_prog = ii
  colspec = get_param("mcell_mc2d_gene_shades", package = "metacell")
  colspec = colspec = c("white", colspec[1], colspec[11])
  prog_genes = sort(reg_full_h[ii,], decreasing = T)
  n_words = 10
  barplot(prog_genes[n_words:1], names.arg=names(prog_genes[n_words:1]), horiz = T,col=leaders_col[ii],las=2, cex.names = 1)
  dev.off()
}

prog_by_t_all_pd1 = c()
prog_by_t_all_ctrl = c()
prog_by_t_all_41bb = c()
prog_by_t_all_41bb_pd1 = c()
prog_by_t_all_pd1_mc38 = c()
prog_by_t_all_ctrl_mc38 = c()
prog_by_t_all_41bb_mc38 = c()
prog_by_t_all_41bb_pd1_mc38 = c()
prog_by_t_all_pd1_endo = c()
prog_by_t_all_ctrl_endo = c()
prog_by_t_all_41bb_endo = c()
prog_by_t_all_41bb_pd1_endo = c()
#rel_cells = names(mc@mc)
for (t in t_points){
  min_cells = 100
  cond_cells = rel_cells[times == t & is_pd1 & is_spec & is_tumor]
  if (length(cond_cells) > min_cells){
    prog_by_t = apply(norm_global_w[mc@mc[cond_cells],], 2, mean)
  }
  else{
    prog_by_t = rep(0, ncol(norm_global_w))
  }
  prog_by_t_all_pd1 = rbind(prog_by_t_all_pd1, prog_by_t)
  
  cond_cells = rel_cells[times == t & is_ctrl & is_spec & is_tumor]
  if (length(cond_cells) > min_cells){
    prog_by_t = apply(norm_global_w[mc@mc[cond_cells],], 2, mean)
  }
  else{
    prog_by_t = rep(0, ncol(norm_global_w))
  }
  prog_by_t_all_ctrl = rbind(prog_by_t_all_ctrl, prog_by_t)
  
  cond_cells = rel_cells[times == t & is_pd1 & is_b16_endo & is_tumor]
  if (length(cond_cells) > min_cells){
    prog_by_t = apply(norm_global_w[mc@mc[cond_cells],], 2, mean)
  }
  else{
    prog_by_t = rep(0, ncol(norm_global_w))
  }
  prog_by_t_all_pd1_endo = rbind(prog_by_t_all_pd1_endo, prog_by_t)
  
  cond_cells = rel_cells[times == t & is_ctrl & is_b16_endo & is_tumor]
  if (length(cond_cells) > min_cells){
    prog_by_t = apply(norm_global_w[mc@mc[cond_cells],], 2, mean)
  }
  else{
    prog_by_t = rep(0, ncol(norm_global_w))
  }
  prog_by_t_all_ctrl_endo = rbind(prog_by_t_all_ctrl_endo, prog_by_t)
  
  cond_cells = rel_cells[times == t & is_41bb & is_b16_endo & is_tumor]
  if (length(cond_cells) > min_cells){
    prog_by_t = apply(norm_global_w[mc@mc[cond_cells],], 2, mean)
  }
  else{
    prog_by_t = rep(0, ncol(norm_global_w))
  }
  prog_by_t_all_41bb_endo = rbind(prog_by_t_all_41bb_endo, prog_by_t)
  
  cond_cells = rel_cells[times == t & is_41bb_pd1 & is_b16_endo & is_tumor]
  if (length(cond_cells) > min_cells){
    prog_by_t = apply(norm_global_w[mc@mc[cond_cells],], 2, mean)
  }
  else{
    prog_by_t = rep(0, ncol(norm_global_w))
  }
  prog_by_t_all_41bb_pd1_endo = rbind(prog_by_t_all_41bb_pd1_endo, prog_by_t)
  
  cond_cells = rel_cells[times == t & is_41bb & is_spec & is_tumor]
  if (length(cond_cells) > min_cells){
    prog_by_t = apply(norm_global_w[mc@mc[cond_cells],], 2, mean)
  }
  else{
    prog_by_t = rep(0, ncol(norm_global_w))
  }
  prog_by_t_all_41bb = rbind(prog_by_t_all_41bb, prog_by_t)
  
  cond_cells = rel_cells[times == t & is_41bb_pd1 & is_spec & is_tumor]
  if (length(cond_cells) > min_cells){
    prog_by_t = apply(norm_global_w[mc@mc[cond_cells],], 2, mean)
  }
  else{
    prog_by_t = rep(0, ncol(norm_global_w))
  }
  prog_by_t_all_41bb_pd1 = rbind(prog_by_t_all_41bb_pd1, prog_by_t)
  
  cond_cells = rel_cells[times == t & is_pd1 & is_mc38 & is_tumor]
  if (length(cond_cells) > min_cells){
    prog_by_t = apply(norm_global_w[mc@mc[cond_cells],], 2, mean)
  }
  else{
    prog_by_t = rep(0, ncol(norm_global_w))
  }
  prog_by_t_all_pd1_mc38 = rbind(prog_by_t_all_pd1_mc38, prog_by_t)
  
  cond_cells = rel_cells[times == t & is_ctrl & is_mc38 & is_tumor]
  if (length(cond_cells) > min_cells){
    prog_by_t = apply(norm_global_w[mc@mc[cond_cells],], 2, mean)
  }
  else{
    prog_by_t = rep(0, ncol(norm_global_w))
  }
  prog_by_t_all_ctrl_mc38 = rbind(prog_by_t_all_ctrl_mc38, prog_by_t)
  
  cond_cells = rel_cells[times == t & is_41bb & is_mc38 & is_tumor]
  if (length(cond_cells) > min_cells){
    prog_by_t = apply(norm_global_w[mc@mc[cond_cells],], 2, mean)
  }
  else{
    prog_by_t = rep(0, ncol(norm_global_w))
  }
  prog_by_t_all_41bb_mc38 = rbind(prog_by_t_all_41bb_mc38, prog_by_t)
  
  cond_cells = rel_cells[times == t & is_41bb_pd1 & is_mc38 & is_tumor]
  if (length(cond_cells) > min_cells){
    prog_by_t = apply(norm_global_w[mc@mc[cond_cells],], 2, mean)
  }
  else{
    prog_by_t = rep(0, ncol(norm_global_w))
  }
  prog_by_t_all_41bb_pd1_mc38 = rbind(prog_by_t_all_41bb_pd1_mc38, prog_by_t)
}
for (prog in 1:ncol(norm_global_w)){
  t_s = c(3:8)
  t_s_41bb = c(4,7)
  png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/", "prog_on_time_",prog,"_pd1_ctrl.png",sep=""), 12 * 100, 12 * 100, res = 200)
  plot(t_points[t_s], prog_by_t_all_pd1[t_s,prog], col="darkblue", pch=19, ylim=c(0, max(c(prog_by_t_all_pd1[t_s,prog], prog_by_t_all_ctrl[t_s,prog], prog_by_t_all_41bb_pd1[t_s_41bb,prog]))), ylab=paste(prog, "_prog_norm_exp", sep=""))
  points(t_points[t_s], prog_by_t_all_ctrl[t_s,prog], col="lightblue", pch=19)
  points(t_points[t_s_41bb], prog_by_t_all_41bb_pd1[t_s_41bb,prog], col="purple", pch=19)
  points(t_points[t_s_41bb], prog_by_t_all_41bb[t_s_41bb,prog], col="plum", pch=19)
  dev.off()
  
  png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/", "prog_on_time_with_endo_",prog,"_pd1_ctrl_41bb.png",sep=""), 12 * 100, 12 * 100, res = 200)
  plot(t_points[t_s], prog_by_t_all_pd1[t_s,prog], col="darkblue", pch=19, ylim=c(0, max(c(prog_by_t_all_pd1[t_s,prog], prog_by_t_all_ctrl[t_s,prog], prog_by_t_all_pd1_endo[t_s,prog], prog_by_t_all_ctrl_endo[t_s,prog], prog_by_t_all_41bb_pd1[t_s_41bb,prog], prog_by_t_all_41bb[t_s_41bb,prog], prog_by_t_all_pd1_mc38[t_s_41bb,prog], prog_by_t_all_ctrl_mc38[t_s_41bb,prog], prog_by_t_all_41bb_pd1_mc38[t_s_41bb,prog], prog_by_t_all_41bb_mc38[t_s_41bb,prog], prog_by_t_all_41bb_endo[t_s_41bb,prog], prog_by_t_all_41bb_pd1_endo[t_s_41bb,prog]))), ylab=paste(prog, "_prog_norm_exp", sep=""))
  points(t_points[t_s], prog_by_t_all_ctrl[t_s,prog], col="lightblue", pch=19)
  points(t_points[t_s_41bb], prog_by_t_all_41bb_pd1[t_s_41bb,prog], col="purple", pch=19)
  points(t_points[t_s_41bb], prog_by_t_all_41bb[t_s_41bb,prog], col="plum", pch=19)
  points(t_points[t_s_41bb], prog_by_t_all_pd1_mc38[t_s_41bb,prog], col="orange", pch=19)
  points(t_points[t_s_41bb], prog_by_t_all_ctrl_mc38[t_s_41bb,prog], col="yellow", pch=19)
  points(t_points[t_s_41bb], prog_by_t_all_41bb_pd1_mc38[t_s_41bb,prog], col="red", pch=19)
  points(t_points[t_s_41bb], prog_by_t_all_41bb_mc38[t_s_41bb,prog], col="pink", pch=19)
  points(t_points[t_s], prog_by_t_all_pd1_endo[t_s,prog], col="darkgreen", pch=19)
  points(t_points[t_s], prog_by_t_all_ctrl_endo[t_s,prog], col="lightgreen", pch=19)
  points(t_points[t_s_41bb], prog_by_t_all_41bb_endo[t_s_41bb,prog], col="black", pch=19)
  points(t_points[t_s_41bb], prog_by_t_all_41bb_pd1_endo[t_s_41bb,prog], col="grey", pch=19)
  dev.off()
  
  png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/", "prog_on_time_only_endo_",prog,"_pd1_ctrl_41bb.png",sep=""), 12 * 100, 12 * 100, res = 200)
  plot(t_points[t_s], prog_by_t_all_pd1_endo[t_s,prog], col="darkgreen", pch=19, ylim=c(0, max(c(prog_by_t_all_pd1_endo[t_s,prog], prog_by_t_all_ctrl_endo[t_s,prog], prog_by_t_all_pd1_mc38[t_s_41bb,prog], prog_by_t_all_ctrl_mc38[t_s_41bb,prog], prog_by_t_all_41bb_pd1_mc38[t_s_41bb,prog], prog_by_t_all_41bb_mc38[t_s_41bb,prog], prog_by_t_all_41bb_endo[t_s_41bb,prog], prog_by_t_all_41bb_pd1_endo[t_s_41bb,prog]))), ylab=paste(prog, "_prog_norm_exp", sep=""))
  points(t_points[t_s], prog_by_t_all_ctrl_endo[t_s,prog], col="lightgreen", pch=19)
  points(t_points[t_s_41bb], prog_by_t_all_pd1_mc38[t_s_41bb,prog], col="orange", pch=19)
  points(t_points[t_s_41bb], prog_by_t_all_ctrl_mc38[t_s_41bb,prog], col="yellow", pch=19)
  points(t_points[t_s_41bb], prog_by_t_all_41bb_pd1_mc38[t_s_41bb,prog], col="red", pch=19)
  points(t_points[t_s_41bb], prog_by_t_all_41bb_mc38[t_s_41bb,prog], col="pink", pch=19)
  points(t_points[t_s_41bb], prog_by_t_all_41bb_endo[t_s_41bb,prog], col="black", pch=19)
  points(t_points[t_s_41bb], prog_by_t_all_41bb_pd1_endo[t_s_41bb,prog], col="grey", pch=19)
  dev.off()
  
}
png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/", "prog_on_time_legend.png",sep=""), 14 * 100, 14 * 100, res = 200)
plot.new()
legend("topleft", legend=c("spec_pd1", "spec_ctrl", "spec_41bb", "spec_41bb_pd1", "b16_endo_pd1", "b16_endo_ctrl", "b16_endo_41bb" ,"b16_endo_41bb_pd1", "mc38_pd1", "mc38_ctrl", "mc38_41bb", "mc38_41bb_pd1"), pch=19, col=c("darkblue", "lightblue", "plum", "purple", "darkgreen", "lightgreen", "black", "grey", "orange", "yellow", "pink", "red"), cex=2)
dev.off()
png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/", "cond_on_time_legend.png",sep=""), 12 * 100, 13 * 100, res = 200)
plot.new()
legend("topleft", legend=c("Ctrl", "aPD1", "a41bb", "a41bb+pd1"), pch=19, col=c("grey", "black", "pink", "red"), cex=2)
dev.off()

active_default = 0.2
strong_default = 0.6
dens_eps = 1e-7
active_thres_vec = rep(active_default, ncol(norm_global_mat))
strong_thres_vec = rep(strong_default, ncol(norm_global_mat))
active_thres_vec[c(6)] = 0.1
active_thres_vec[c(1,10)] = 0.3
# active_thres_vec[c(4,6)] = 0.3
# active_thres_vec[c(3,6,10)] = 0.5
# strong_thres_vec[c(6,9)] = 0.25
for (prog in 1:ncol(norm_global_mat)){
  active_thres = active_thres_vec[prog]
  strong_thres = strong_thres_vec[prog]
  png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_distributions/", "prog_linear_act_",prog, active_thres, ".png",sep=""), 12 * 100, 12 * 100, res = 200)
  #plot(sort(log2(norm_global_mat[rel_mcs,prog] + dens_eps)))
  plot(sort(norm_global_mat[rel_mcs,prog]))
  abline(h=active_thres, lty=2, col="navyblue")
  abline(h=strong_thres, lty=2, col="springgreen")
  dev.off()
}
for (prog in 1:ncol(norm_global_mat)){
  active_thres = active_thres_vec[prog]
  strong_thres = strong_thres_vec[prog]
  png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_distributions/", "prog_linear_dist_",prog, active_thres, ".png",sep=""), 12 * 100, 12 * 100, res = 200)
  #plot(sort(log2(norm_global_mat[rel_mcs,prog] + dens_eps)))
  plot(density(norm_global_w[rel_mcs,prog]))
  abline(v=active_thres, lty=2, col="navyblue")
  abline(v=strong_thres, lty=2, col="springgreen")
  dev.off()
}

active_thres_vec_log = log2(active_thres_vec + dens_eps)
# active_thres_vec[c(1,5,7)] = log2(0.07 + dens_eps)
strong_thres_vec_log = log2(strong_thres_vec + dens_eps)
# strong_thres_vec[c(3,11,12)] = log2(0.35 + dens_eps)
# strong_thres_vec[c(6,9)] = log2(0.25 + dens_eps)
for (prog in 1:ncol(norm_global_mat)){
  active_thres = active_thres_vec_log[prog]
  strong_thres = strong_thres_vec_log[prog]
  png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_distributions/", "prog_log_act_",prog, active_thres, ".png",sep=""), 12 * 100, 12 * 100, res = 200)
  plot(sort(log2(norm_global_mat[rel_mcs,prog] + dens_eps)))
  #plot(sort(norm_global_mat[rel_mcs,prog]))
  abline(h=active_thres, lty=2, col="navyblue")
  abline(h=strong_thres, lty=2, col="springgreen")
  dev.off()
}
for (prog in 1:ncol(norm_global_mat)){
  active_thres = active_thres_vec_log[prog]
  strong_thres = strong_thres_vec_log[prog]
  png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_distributions/", "prog_log_dist_",prog, active_thres, ".png",sep=""), 12 * 100, 12 * 100, res = 200)
  plot(density(log2(norm_global_w[rel_mcs,prog] + dens_eps)))
  #plot(sort(norm_global_mat[rel_mcs,prog]))
  abline(v=active_thres, lty=2, col="navyblue")
  abline(v=strong_thres, lty=2, col="springgreen")
  dev.off()
}


state_by_t_all_pd1 = c()
state_by_t_all_ctrl = c()
state_by_t_all_41bb = c()
state_by_t_all_41bb_pd1 = c()
state_by_t_all_pd1_mc38 = c()
state_by_t_all_ctrl_mc38 = c()
state_by_t_all_41bb_mc38 = c()
state_by_t_all_41bb_pd1_mc38 = c()
state_by_t_all_pd1_endo = c()
state_by_t_all_ctrl_endo = c()
state_by_t_all_41bb_endo = c()
state_by_t_all_41bb_pd1_endo = c()


for (prog in 1:ncol(norm_global_mat)){
  active_thres = active_thres_vec[prog]
  strong_thres = strong_thres_vec[prog]
  active_mcs = (1:max(mc@mc))[norm_global_mat[,prog] >= active_thres]
  active_cells = rel_cells[mc@mc[rel_cells] %in% active_mcs]
  
  min_cells = 30
  state_by_t = c()
  for (t in t_points){
    cond_cells = rel_cells[times == t & is_pd1 & is_spec & is_tumor]
    num_cells = length(cond_cells)
    if (num_cells > min_cells){
      state_by_t = c(state_by_t, sum(cond_cells %in% active_cells)/num_cells)
    }
    else{
      state_by_t = c(state_by_t, 0)
    }
  }
  state_by_t_all_pd1 = rbind(state_by_t_all_pd1, state_by_t)
  
  state_by_t = c()
  for (t in t_points){
    cond_cells = rel_cells[times == t & is_ctrl & is_spec & is_tumor]
    num_cells = length(cond_cells)
    if (num_cells > min_cells){
      state_by_t = c(state_by_t, sum(cond_cells %in% active_cells)/num_cells)
    }
    else{
      state_by_t = c(state_by_t, 0)
    }
  }
  state_by_t_all_ctrl = rbind(state_by_t_all_ctrl, state_by_t)
  
  state_by_t = c()
  for (t in t_points){
    cond_cells = rel_cells[times == t & is_pd1 & is_b16_endo & is_tumor]
    num_cells = length(cond_cells)
    if (num_cells > min_cells){
      state_by_t = c(state_by_t, sum(cond_cells %in% active_cells)/num_cells)
    }
    else{
      state_by_t = c(state_by_t, 0)
    }
  }
  state_by_t_all_pd1_endo = rbind(state_by_t_all_pd1_endo, state_by_t)
  
  state_by_t = c()
  for (t in t_points){
    cond_cells = rel_cells[times == t & is_ctrl & is_b16_endo & is_tumor]
    num_cells = length(cond_cells)
    if (num_cells > min_cells){
      state_by_t = c(state_by_t, sum(cond_cells %in% active_cells)/num_cells)
    }
    else{
      state_by_t = c(state_by_t, 0)
    }
  }
  state_by_t_all_ctrl_endo = rbind(state_by_t_all_ctrl_endo, state_by_t)
  
  state_by_t = c()
  for (t in t_points){
    cond_cells = rel_cells[times == t & is_41bb & is_b16_endo & is_tumor]
    num_cells = length(cond_cells)
    if (num_cells > min_cells){
      state_by_t = c(state_by_t, sum(cond_cells %in% active_cells)/num_cells)
    }
    else{
      state_by_t = c(state_by_t, 0)
    }
  }
  state_by_t_all_41bb_endo = rbind(state_by_t_all_41bb_endo, state_by_t)
  
  state_by_t = c()
  for (t in t_points){
    cond_cells = rel_cells[times == t & is_41bb_pd1 & is_b16_endo & is_tumor]
    num_cells = length(cond_cells)
    if (num_cells > min_cells){
      state_by_t = c(state_by_t, sum(cond_cells %in% active_cells)/num_cells)
    }
    else{
      state_by_t = c(state_by_t, 0)
    }
  }
  state_by_t_all_41bb_pd1_endo = rbind(state_by_t_all_41bb_pd1_endo, state_by_t)
  
  state_by_t = c()
  for (t in t_points){
    cond_cells = rel_cells[times == t & is_41bb & is_spec & is_tumor]
    num_cells = length(cond_cells)
    if (num_cells > min_cells){
      state_by_t = c(state_by_t, sum(cond_cells %in% active_cells)/num_cells)
    }
    else{
      state_by_t = c(state_by_t, 0)
    }
  }
  state_by_t_all_41bb = rbind(state_by_t_all_41bb, state_by_t)
  
  state_by_t = c()
  for (t in t_points){
    cond_cells = rel_cells[times == t & is_41bb_pd1 & is_spec & is_tumor]
    num_cells = length(cond_cells)
    if (num_cells > min_cells){
      state_by_t = c(state_by_t, sum(cond_cells %in% active_cells)/num_cells)
    }
    else{
      state_by_t = c(state_by_t, 0)
    }
  }
  state_by_t_all_41bb_pd1 = rbind(state_by_t_all_41bb_pd1, state_by_t)
  
  state_by_t = c()
  for (t in t_points){
    cond_cells = rel_cells[times == t & is_pd1 & is_mc38 & is_tumor]
    num_cells = length(cond_cells)
    if (num_cells > min_cells){
      state_by_t = c(state_by_t, sum(cond_cells %in% active_cells)/num_cells)
    }
    else{
      state_by_t = c(state_by_t, 0)
    }
  }
  state_by_t_all_pd1_mc38 = rbind(state_by_t_all_pd1_mc38, state_by_t)
  
  state_by_t = c()
  for (t in t_points){
    cond_cells = rel_cells[times == t & is_ctrl & is_mc38 & is_tumor]
    num_cells = length(cond_cells)
    if (num_cells > min_cells){
      state_by_t = c(state_by_t, sum(cond_cells %in% active_cells)/num_cells)
    }
    else{
      state_by_t = c(state_by_t, 0)
    }
  }
  state_by_t_all_ctrl_mc38 = rbind(state_by_t_all_ctrl_mc38, state_by_t)
  
  state_by_t = c()
  for (t in t_points){
    cond_cells = rel_cells[times == t & is_41bb & is_mc38 & is_tumor]
    num_cells = length(cond_cells)
    if (num_cells > min_cells){
      state_by_t = c(state_by_t, sum(cond_cells %in% active_cells)/num_cells)
    }
    else{
      state_by_t = c(state_by_t, 0)
    }
  }
  state_by_t_all_41bb_mc38 = rbind(state_by_t_all_41bb_mc38, state_by_t)
  
  state_by_t = c()
  for (t in t_points){
    cond_cells = rel_cells[times == t & is_41bb_pd1 & is_mc38 & is_tumor]
    num_cells = length(cond_cells)
    if (num_cells > min_cells){
      state_by_t = c(state_by_t, sum(cond_cells %in% active_cells)/num_cells)
    }
    else{
      state_by_t = c(state_by_t, 0)
    }
  }
  state_by_t_all_41bb_pd1_mc38 = rbind(state_by_t_all_41bb_pd1_mc38, state_by_t)
  
}

for (prog in 1:ncol(norm_global_w)){
  active_thres = active_thres_vec[prog]
  strong_thres = strong_thres_vec[prog]
  t_s = c(3:8)
  t_s_41bb = c(4,7)
  png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_distributions/", "barplot_active_time_ot1_",prog, active_thres, ".png",sep=""), 12 * 100, 12 * 100, res = 200)
  barplot(rbind(state_by_t_all_ctrl[prog,t_s], state_by_t_all_pd1[prog,t_s], state_by_t_all_41bb[prog,t_s], state_by_t_all_41bb_pd1[prog,t_s]), names.arg=t_points[t_s], beside=T, col = c("lightblue", "darkblue", "plum", "purple"))
  dev.off()
  
  png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_distributions/", "barplot_active_time_with_endo_",prog, active_thres, ".png",sep=""), 12 * 100, 12 * 100, res = 200)
  barplot(rbind(state_by_t_all_ctrl[prog,t_s_41bb], state_by_t_all_pd1[prog,t_s_41bb], state_by_t_all_41bb[prog,t_s_41bb], state_by_t_all_41bb_pd1[prog,t_s_41bb], state_by_t_all_ctrl_endo[prog,t_s_41bb], state_by_t_all_pd1_endo[prog,t_s_41bb], state_by_t_all_41bb_endo[prog,t_s_41bb], state_by_t_all_41bb_pd1_endo[prog,t_s_41bb], state_by_t_all_ctrl_mc38[prog,t_s_41bb], state_by_t_all_pd1_mc38[prog,t_s_41bb], state_by_t_all_41bb_mc38[prog,t_s_41bb], state_by_t_all_41bb_pd1_mc38[prog,t_s_41bb]), names.arg=t_points[t_s_41bb], beside=T, col = c("lightblue", "darkblue", "plum", "purple", "lightgreen", "darkgreen", "black", "grey", "yellow", "orange", "pink", "red"))
  dev.off()
  
  png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_distributions/", "barplot_active_time_only_endo_",prog, active_thres, ".png",sep=""), 12 * 100, 12 * 100, res = 200)
  barplot(rbind(state_by_t_all_ctrl_endo[prog,t_s_41bb], state_by_t_all_pd1_endo[prog,t_s_41bb], state_by_t_all_41bb_endo[prog,t_s_41bb], state_by_t_all_41bb_pd1_endo[prog,t_s_41bb], state_by_t_all_ctrl_mc38[prog,t_s_41bb], state_by_t_all_pd1_mc38[prog,t_s_41bb], state_by_t_all_41bb_mc38[prog,t_s_41bb], state_by_t_all_41bb_pd1_mc38[prog,t_s_41bb]), names.arg=t_points[t_s_41bb], beside=T, col = c("lightgreen", "darkgreen", "black", "grey", "yellow", "orange", "pink", "red"))
  dev.off()
  
}
t_s_41bb = c(4,7)
for (prog in 1:ncol(norm_global_w)){
  active_thres = active_thres_vec[prog]
  strong_thres = strong_thres_vec[prog]
  for (tt in t_s_41bb){
    act_mat = matrix(rep(0, 4*3), nrow=4, ncol=3)
    rownames(act_mat) = c("Ctrl", "aPD1", "a41BB", "aPD1+a41BB")
    colnames(act_mat) = c("OT1", "B16-endo", "MC38-endo")
    act_mat["Ctrl", "OT1"] = state_by_t_all_ctrl[prog,tt]
    act_mat["aPD1", "OT1"] = state_by_t_all_pd1[prog,tt]
    act_mat["a41BB", "OT1"] = state_by_t_all_41bb[prog,tt]
    act_mat["aPD1+a41BB", "OT1"] = state_by_t_all_41bb_pd1[prog,tt]
    act_mat["Ctrl", "B16-endo"] = state_by_t_all_ctrl_endo[prog,tt]
    act_mat["aPD1", "B16-endo"] = state_by_t_all_pd1_endo[prog,tt]
    act_mat["a41BB", "B16-endo"] = state_by_t_all_41bb_endo[prog,tt]
    act_mat["aPD1+a41BB", "B16-endo"] = state_by_t_all_41bb_pd1_endo[prog,tt]
    act_mat["Ctrl", "MC38-endo"] = state_by_t_all_ctrl_mc38[prog,tt]
    act_mat["aPD1", "MC38-endo"] = state_by_t_all_pd1_mc38[prog,tt]
    act_mat["a41BB", "MC38-endo"] = state_by_t_all_pd1_mc38[prog,tt]
    act_mat["aPD1+a41BB", "MC38-endo"] = state_by_t_all_41bb_mc38[prog,tt]
    
    png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_distributions/", "barplot_active_with_endo_",prog, "_time_", t_points[tt], "_", active_thres, ".png",sep=""), 12 * 100, 12 * 100, res = 200)
    barplot(act_mat, names.arg=colnames(act_mat), beside=T, col = c("grey", "black", "pink", "red"), main=paste("T=", t_points[tt], sep=""), ylim=c(0,0.6))
    dev.off()
  }
}


state_by_t_all_pd1 = c()
state_by_t_all_ctrl = c()
state_by_t_all_41bb = c()
state_by_t_all_41bb_pd1 = c()
state_by_t_all_pd1_mc38 = c()
state_by_t_all_ctrl_mc38 = c()
state_by_t_all_41bb_mc38 = c()
state_by_t_all_41bb_pd1_mc38 = c()
state_by_t_all_pd1_endo = c()
state_by_t_all_ctrl_endo = c()
state_by_t_all_41bb_endo = c()
state_by_t_all_41bb_pd1_endo = c()

for (prog in 1:ncol(norm_global_mat)){
  active_thres = active_thres_vec[prog]
  strong_thres = strong_thres_vec[prog]
  strong_mcs = (1:max(mc@mc))[norm_global_mat[,prog] >= strong_thres]
  strong_cells = rel_cells[mc@mc[rel_cells] %in% strong_mcs]
  min_cells = 30
  state_by_t = c()
  for (t in t_points){
    cond_cells = rel_cells[times == t & is_pd1 & is_spec & is_tumor]
    num_cells = length(cond_cells)
    if (num_cells > min_cells){
      state_by_t = c(state_by_t, sum(cond_cells %in% strong_cells)/num_cells)
    }
    else{
      state_by_t = c(state_by_t, 0)
    }
  }
  state_by_t_all_pd1 = rbind(state_by_t_all_pd1, state_by_t)
  
  state_by_t = c()
  for (t in t_points){
    cond_cells = rel_cells[times == t & is_ctrl & is_spec & is_tumor]
    num_cells = length(cond_cells)
    if (num_cells > min_cells){
      state_by_t = c(state_by_t, sum(cond_cells %in% strong_cells)/num_cells)
    }
    else{
      state_by_t = c(state_by_t, 0)
    }
  }
  state_by_t_all_ctrl = rbind(state_by_t_all_ctrl, state_by_t)
  
  state_by_t = c()
  for (t in t_points){
    cond_cells = rel_cells[times == t & is_pd1 & is_b16_endo & is_tumor]
    num_cells = length(cond_cells)
    if (num_cells > min_cells){
      state_by_t = c(state_by_t, sum(cond_cells %in% strong_cells)/num_cells)
    }
    else{
      state_by_t = c(state_by_t, 0)
    }
  }
  state_by_t_all_pd1_endo = rbind(state_by_t_all_pd1_endo, state_by_t)
  
  state_by_t = c()
  for (t in t_points){
    cond_cells = rel_cells[times == t & is_ctrl & is_b16_endo & is_tumor]
    num_cells = length(cond_cells)
    if (num_cells > min_cells){
      state_by_t = c(state_by_t, sum(cond_cells %in% strong_cells)/num_cells)
    }
    else{
      state_by_t = c(state_by_t, 0)
    }
  }
  state_by_t_all_ctrl_endo = rbind(state_by_t_all_ctrl_endo, state_by_t)
  
  state_by_t = c()
  for (t in t_points){
    cond_cells = rel_cells[times == t & is_41bb & is_b16_endo & is_tumor]
    num_cells = length(cond_cells)
    if (num_cells > min_cells){
      state_by_t = c(state_by_t, sum(cond_cells %in% strong_cells)/num_cells)
    }
    else{
      state_by_t = c(state_by_t, 0)
    }
  }
  state_by_t_all_41bb_endo = rbind(state_by_t_all_41bb_endo, state_by_t)
  
  state_by_t = c()
  for (t in t_points){
    cond_cells = rel_cells[times == t & is_41bb_pd1 & is_b16_endo & is_tumor]
    num_cells = length(cond_cells)
    if (num_cells > min_cells){
      state_by_t = c(state_by_t, sum(cond_cells %in% strong_cells)/num_cells)
    }
    else{
      state_by_t = c(state_by_t, 0)
    }
  }
  state_by_t_all_41bb_pd1_endo = rbind(state_by_t_all_41bb_pd1_endo, state_by_t)
  
  state_by_t = c()
  for (t in t_points){
    cond_cells = rel_cells[times == t & is_41bb & is_spec & is_tumor]
    num_cells = length(cond_cells)
    if (num_cells > min_cells){
      state_by_t = c(state_by_t, sum(cond_cells %in% strong_cells)/num_cells)
    }
    else{
      state_by_t = c(state_by_t, 0)
    }
  }
  state_by_t_all_41bb = rbind(state_by_t_all_41bb, state_by_t)
  
  state_by_t = c()
  for (t in t_points){
    cond_cells = rel_cells[times == t & is_41bb_pd1 & is_spec & is_tumor]
    num_cells = length(cond_cells)
    if (num_cells > min_cells){
      state_by_t = c(state_by_t, sum(cond_cells %in% strong_cells)/num_cells)
    }
    else{
      state_by_t = c(state_by_t, 0)
    }
  }
  state_by_t_all_41bb_pd1 = rbind(state_by_t_all_41bb_pd1, state_by_t)
  
  state_by_t = c()
  for (t in t_points){
    cond_cells = rel_cells[times == t & is_pd1 & is_mc38 & is_tumor]
    num_cells = length(cond_cells)
    if (num_cells > min_cells){
      state_by_t = c(state_by_t, sum(cond_cells %in% strong_cells)/num_cells)
    }
    else{
      state_by_t = c(state_by_t, 0)
    }
  }
  state_by_t_all_pd1_mc38 = rbind(state_by_t_all_pd1_mc38, state_by_t)
  
  state_by_t = c()
  for (t in t_points){
    cond_cells = rel_cells[times == t & is_ctrl & is_mc38 & is_tumor]
    num_cells = length(cond_cells)
    if (num_cells > min_cells){
      state_by_t = c(state_by_t, sum(cond_cells %in% strong_cells)/num_cells)
    }
    else{
      state_by_t = c(state_by_t, 0)
    }
  }
  state_by_t_all_ctrl_mc38 = rbind(state_by_t_all_ctrl_mc38, state_by_t)
  
  state_by_t = c()
  for (t in t_points){
    cond_cells = rel_cells[times == t & is_41bb & is_mc38 & is_tumor]
    num_cells = length(cond_cells)
    if (num_cells > min_cells){
      state_by_t = c(state_by_t, sum(cond_cells %in% strong_cells)/num_cells)
    }
    else{
      state_by_t = c(state_by_t, 0)
    }
  }
  state_by_t_all_41bb_mc38 = rbind(state_by_t_all_41bb_mc38, state_by_t)
  
  state_by_t = c()
  for (t in t_points){
    cond_cells = rel_cells[times == t & is_41bb_pd1 & is_mc38 & is_tumor]
    num_cells = length(cond_cells)
    if (num_cells > min_cells){
      state_by_t = c(state_by_t, sum(cond_cells %in% strong_cells)/num_cells)
    }
    else{
      state_by_t = c(state_by_t, 0)
    }
  }
  state_by_t_all_41bb_pd1_mc38 = rbind(state_by_t_all_41bb_pd1_mc38, state_by_t)
  
}
for (prog in 1:ncol(norm_global_w)){
  active_thres = active_thres_vec[prog]
  strong_thres = strong_thres_vec[prog]
  t_s = c(3:8)
  t_s_41bb = c(4,7)
  png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_distributions/", "barplot_strong_time_ot1_",prog, active_thres, ".png",sep=""), 12 * 100, 12 * 100, res = 200)
  barplot(rbind(state_by_t_all_ctrl[prog,t_s], state_by_t_all_pd1[prog,t_s], state_by_t_all_41bb[prog,t_s], state_by_t_all_41bb_pd1[prog,t_s]), names.arg=t_points[t_s], beside=T, col = c("lightblue", "darkblue", "plum", "purple"))
  dev.off()
  
  png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_distributions/", "barplot_strong_time_with_endo_",prog, active_thres, ".png",sep=""), 12 * 100, 12 * 100, res = 200)
  barplot(rbind(state_by_t_all_ctrl[prog,t_s_41bb], state_by_t_all_pd1[prog,t_s_41bb], state_by_t_all_41bb[prog,t_s_41bb], state_by_t_all_41bb_pd1[prog,t_s_41bb], state_by_t_all_ctrl_endo[prog,t_s_41bb], state_by_t_all_pd1_endo[prog,t_s_41bb], state_by_t_all_41bb_endo[prog,t_s_41bb], state_by_t_all_41bb_pd1_endo[prog,t_s_41bb], state_by_t_all_ctrl_mc38[prog,t_s_41bb], state_by_t_all_pd1_mc38[prog,t_s_41bb], state_by_t_all_41bb_mc38[prog,t_s_41bb], state_by_t_all_41bb_pd1_mc38[prog,t_s_41bb]), names.arg=t_points[t_s_41bb], beside=T, col = c("lightblue", "darkblue", "plum", "purple", "lightgreen", "darkgreen", "black", "grey", "yellow", "orange", "pink", "red"))
  dev.off()
  
  png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_distributions/", "barplot_strong_time_only_endo_",prog, active_thres, ".png",sep=""), 12 * 100, 12 * 100, res = 200)
  barplot(rbind(state_by_t_all_ctrl_endo[prog,t_s_41bb], state_by_t_all_pd1_endo[prog,t_s_41bb], state_by_t_all_41bb_endo[prog,t_s_41bb], state_by_t_all_41bb_pd1_endo[prog,t_s_41bb], state_by_t_all_ctrl_mc38[prog,t_s_41bb], state_by_t_all_pd1_mc38[prog,t_s_41bb], state_by_t_all_41bb_mc38[prog,t_s_41bb], state_by_t_all_41bb_pd1_mc38[prog,t_s_41bb]), names.arg=t_points[t_s_41bb], beside=T, col = c("lightgreen", "darkgreen", "black", "grey", "yellow", "orange", "pink", "red"))
  dev.off()
  
}
t_s_41bb = c(4,7)
for (prog in 1:ncol(norm_global_w)){
  active_thres = active_thres_vec[prog]
  strong_thres = strong_thres_vec[prog]
  for (tt in t_s_41bb){
    act_mat = matrix(rep(0, 4*3), nrow=4, ncol=3)
    rownames(act_mat) = c("Ctrl", "aPD1", "a41BB", "aPD1+a41BB")
    colnames(act_mat) = c("OT1", "B16-endo", "MC38-endo")
    act_mat["Ctrl", "OT1"] = state_by_t_all_ctrl[prog,tt]
    act_mat["aPD1", "OT1"] = state_by_t_all_pd1[prog,tt]
    act_mat["a41BB", "OT1"] = state_by_t_all_41bb[prog,tt]
    act_mat["aPD1+a41BB", "OT1"] = state_by_t_all_41bb_pd1[prog,tt]
    act_mat["Ctrl", "B16-endo"] = state_by_t_all_ctrl_endo[prog,tt]
    act_mat["aPD1", "B16-endo"] = state_by_t_all_pd1_endo[prog,tt]
    act_mat["a41BB", "B16-endo"] = state_by_t_all_41bb_endo[prog,tt]
    act_mat["aPD1+a41BB", "B16-endo"] = state_by_t_all_41bb_pd1_endo[prog,tt]
    act_mat["Ctrl", "MC38-endo"] = state_by_t_all_ctrl_mc38[prog,tt]
    act_mat["aPD1", "MC38-endo"] = state_by_t_all_pd1_mc38[prog,tt]
    act_mat["a41BB", "MC38-endo"] = state_by_t_all_pd1_mc38[prog,tt]
    act_mat["aPD1+a41BB", "MC38-endo"] = state_by_t_all_41bb_mc38[prog,tt]
    
    png(paste(scfigs_dir,"/nmf_plots_mc2/umap_full_progs/prog_distributions/", "barplot_strong_with_endo_",prog, "_time_", t_points[tt], "_", active_thres, ".png",sep=""), 12 * 100, 12 * 100, res = 200)
    barplot(act_mat, names.arg=colnames(act_mat), beside=T, col = c("grey", "black", "pink", "red"), main=paste("T=", t_points[tt], sep=""), ylim=c(0,0.6))
    dev.off()
  }
}


nmf_dir = "nmf_plots_mc2_burst"
tf_genes = tf_list_all_non_lat[tf_list_all_non_lat %in% rownames(mc@mc_fp)]
tf_egc = log2(mc@e_gc[tf_genes,] + 1e-5)
prog_score = norm_global_w
tf_pg_cor = tgs_cor(t(tf_egc), prog_score)
cor_no_na = rowSums(is.na(tf_pg_cor)) == 0
tf_pg_cor = tf_pg_cor[cor_no_na,]

tf_genes = setdiff(tf_genes, "Havcr2")

tf_cor_thresh = 0.3
for (gp_ind in 1:ncol(norm_global_w)){
  gp_col = leaders_col[gp_ind]
  gp_name = col2group[gp_col]
  gp_name = gsub("/", "-", gp_name)
  c_gp = sort(tf_pg_cor[,gp_ind])
  c_gp_pos = c_gp[c_gp > tf_cor_thresh]
  c_gp_neg = c_gp[c_gp < -(tf_cor_thresh - 0.1)]
  if(length(c_gp_pos > 0)){
    if (length(c_gp_pos) > 10){
      c_gp_pos = tail(c_gp_pos,10)
    }
    png(paste(scfigs_dir,"/",nmf_dir, "/tfs_gp_pos_cor_", gp_ind,"_",gp_name, "_", gp_col, ".png", sep=""), 200, 400)
    par(mar=c(5,8,4,1)+.1)
    barplot(c_gp_pos, col=gp_col, horiz = T, las = 2, cex.names = 1.5)
    dev.off()
  }
  if(length(c_gp_neg > 0)){
    if (length(c_gp_neg) > 10){
      c_gp_neg = head(c_gp_neg,10)
    }
    png(paste(scfigs_dir,"/",nmf_dir, "/tfs_gp_neg_cor_", gp_ind,"_",gp_name, "_", gp_col, ".png", sep=""), 200, 400)
    par(mar=c(5,8,4,1)+.1)
    barplot(sort(c_gp_neg, decreasing = T), col=gp_col, horiz = T, las = 2, cex.names = 1.5)
    dev.off()
  }
}

tf_genes = surfs_list[surfs_list %in% rownames(mc@mc_fp)]
tf_egc = log2(mc@e_gc[tf_genes,] + 1e-5)
prog_score = norm_global_w
tf_pg_cor = tgs_cor(t(tf_egc), prog_score)
cor_no_na = rowSums(is.na(tf_pg_cor)) == 0
tf_pg_cor = tf_pg_cor[cor_no_na,]

tf_cor_thresh = 0.3
for (gp_ind in 1:ncol(norm_global_w)){
  gp_col = leaders_col[gp_ind]
  gp_name = col2group[gp_col]
  gp_name = gsub("/", "-", gp_name)
  c_gp = sort(tf_pg_cor[,gp_ind])
  c_gp_pos = c_gp[c_gp > tf_cor_thresh]
  c_gp_neg = c_gp[c_gp < -(tf_cor_thresh - 0.1)]
  if(length(c_gp_pos > 0)){
    if (length(c_gp_pos) > 10){
      c_gp_pos = tail(c_gp_pos,10)
    }
    png(paste(scfigs_dir,"/",nmf_dir, "/surface_gp_pos_cor_", gp_ind,"_",gp_name, "_", gp_col, ".png", sep=""), 200, 400)
    par(mar=c(5,8,4,1)+.1)
    barplot(c_gp_pos, col=gp_col, horiz = T, las = 2, cex.names = 1.5)
    dev.off()
  }
  if(length(c_gp_neg > 0)){
    if (length(c_gp_neg) > 10){
      c_gp_neg = head(c_gp_neg,10)
    }
    png(paste(scfigs_dir,"/",nmf_dir, "/surface_gp_neg_cor_", gp_ind,"_",gp_name, "_", gp_col, ".png", sep=""), 200, 400)
    par(mar=c(5,8,4,1)+.1)
    barplot(sort(c_gp_neg, decreasing = T), col=gp_col, horiz = T, las = 2, cex.names = 1.5)
    dev.off()
  }
}

tf_genes = tf_list_all_non_lat
surf_genes = c(surfs_list[surfs_list %in% rownames(mc@mc_fp)], "Prf1")
tf_surf_genes = c(rep(1, length(tf_genes)),rep(2, length(surf_genes)))
names(tf_surf_genes) = c(tf_genes, surf_genes)

paletteLength <- 500
myColor <- colorRampPalette(c("white", "steelblue"))(paletteLength)#colorRampPalette(c("blue","green", "white", "orange","red"))(paletteLength)

top_ann = data.frame(row.names=names(leaders_col),main_top=names(leaders_col))#annotation_col=top_ann, annotation_colors=list(main_top=main_prog_cols), 
myBreaks <- c(seq(0, max(norm_full_h[genes_hm,]), length.out=floor(paletteLength)))
#png(paste(scfigs_dir,"/nmf_plots/_194/local_h_heatmap.png",sep=""), 12 * 100, 12 * 100, res = 200)
pheatmap(norm_full_h[genes_hm,], cluster_rows = T, cluster_cols = T, fontsize = 4, color = myColor, breaks = myBreaks, annotation_col=top_ann, annotation_colors=list(main_top=leaders_col), filename=paste("figs","/nmf_plots_mc2_burst/full_h_heatmap.png",sep=""), width=10, height=10, border_color = NA, treeheight_row = 0)

genes_hm_2 = c()
for (dd in 1:nrow(reg_full_h)){
  genes_hm_2 = c(genes_hm_2, names(tail(sort(reg_full_h[dd,]), 25)))
  # if (dd == dend_prog){
  #   dend_genes = names(tail(sort(reg_full_h[dd,]), 25))
  # }
}
genes_hm_2 = unique(genes_hm_2)
#genes_hm_2 = genes_hm_2[!genes_hm_2 %in% dend_genes]
rownames(reg_full_h) = names(leaders_col)
myBreaks <- c(seq(0, max(norm_full_h[genes_hm_2,]), length.out=floor(paletteLength)))
#png(paste(scfigs_dir,"/nmf_plots/_194/local_h_heatmap.png",sep=""), 12 * 100, 12 * 100, res = 200)
pheatmap(norm_full_h[genes_hm_2,], cluster_rows = T, cluster_cols = T, fontsize = 4, color = myColor, breaks = myBreaks, annotation_col=top_ann, annotation_colors=list(main_top=leaders_col), filename=paste("figs","/nmf_plots_mc2_burst/full_h_heatmap_no_norm.png",sep=""), width=10, height=10, border_color = NA, treeheight_row = 0)
#dev.off()

dir.create(paste0("figs","/nmf_plots_mc2_burst/gg_plots"),recursive = T)
gg_dir = paste("figs","/nmf_plots_mc2_burst/gg_plots/", sep="")
gns_1 = c("Nr4a1", "Myc", "Ifng", "Utf1", "Stat1", "Sp100", "Rbpj", "Rbpj", "Id2", "Id2", "Tox", "Prdm1", "Il2ra", "Id2", "Myb", "Nfatc1", "Ikzf1", "Pdcd1", "Tcf7")
gns_2 = c("Egr2", "Nr4a3", "Cd40lg", "Tnf", "Myc", "Rel", "Adam8", "Irf8", "Prdm1", "Tnfrsf4", "Cd244", "Il2rb", "Il2rb", "Id3", "Zfpm1", "Cd93", "Nt5e", "Tnfrsf9", "Slamf7")
for (g in 1:length(gns_1)){
  mcell_mc_plot_gg(mc_id, gns_1[g], gns_2[g], use_egc = T, fig_fn = paste(gg_dir, gns_1[g], "_", gns_2[g], ".png", sep=""))
}

 

group_mc_by = c('location', 'cell_type', 'days_post_transfer', 'condition', 'mouse_id')
mc_bd = table(mc@mc, apply(md[names(mc@mc), group_mc_by], 1, paste0, collapse=" "))
mc_bd_n = mc_bd / rowSums(mc_bd)

mc_bd_ann = as.data.frame(matrix(unlist(strsplit(colnames(mc_bd_n), split=" ")), ncol=length(group_mc_by), byrow=T))
rownames(mc_bd_ann) = colnames(mc_bd)
colnames(mc_bd_ann) = group_mc_by

mc_bd_ann$mouse_id = as.character(mc_bd_ann$mouse_id)
save(mc_bd_ann, file="mc_bd_ann")

mc_bd_ann_tum = mc_bd_ann[mc_bd_ann$location == "tumor",]
# levels(mc_bd_ann_tum$cell_type) = c(levels(mc_bd_ann_tum$cell_type), "endo_mc38")
# mc_bd_ann_tum$cell_type[mc_bd_ann_tum$mouse_id %in% mc_bd_ann_tum_endo_mc38$mouse_id] = "endo_mc38"
mc_bd_ann_tum = mc_bd_ann_tum[order(mc_bd_ann_tum$cell_type),]
mc_bd_ann_tum = mc_bd_ann_tum[!mc_bd_ann_tum$condition == "na",]

m_id_cond = table(md$mouse_id, md$condition)
m_id_cond = m_id_cond[as.character(unique(mc_bd_ann_tum$mouse_id)),]

mc_bd_ann_tum_3 = mc_bd_ann_tum[mc_bd_ann_tum$days_post_transfer == 3,]
mc_bd_ann_tum_6 = mc_bd_ann_tum[mc_bd_ann_tum$days_post_transfer == 6,]

mc_bd_ann_tum_spec = mc_bd_ann_tum[mc_bd_ann_tum$cell_type == "spec",]

mc_bd_ann_tum_endo_mc38 = mc_bd_ann_tum[mc_bd_ann_tum$cell_type == "endo_mc38",]
mc_bd_ann_tum_endo_b16 = mc_bd_ann_tum[mc_bd_ann_tum$cell_type == "endo",]



state_by_m_spec = c()
state_by_m_mc38 = c()
state_by_m_b16 = c()
state_by_m_byst = c()

cc_by_p_spec = c()
non_filt_gset = scdb_gset(mm_filt_id)
cc_genes = names(non_filt_gset@gene_set)[non_filt_gset@gene_set %in% c(110,11,60,130)]
cc_sig = apply(egc[cc_genes,],2, mean)
cc_sig_by_cell = cc_sig[mc@mc]
names(cc_sig_by_cell) = names(mc@mc)

for (prog in 1:ncol(norm_global_mat)){
  active_thres = active_thres_vec[prog]
  strong_thres = strong_thres_vec[prog]
  active_mcs = (1:max(mc@mc))[norm_global_mat[,prog] >= active_thres]
  active_cells = rel_cells[mc@mc[rel_cells] %in% active_mcs]
  spec_cells = rel_cells[is_spec & is_tumor]
  spec_cells_active = spec_cells[spec_cells %in% active_cells]
  cc_by_p_spec = c(cc_by_p_spec, mean(cc_sig_by_cell[spec_cells_active]))
  
  min_cells = 30
  state_by_t = c()
  for (m_id in unique(as.character(mc_bd_ann_tum$mouse_id))){
    cond_cells = rel_cells[as.character(md$mouse_id) == m_id & is_spec & is_tumor]
    num_cells = length(cond_cells)
    if (num_cells > min_cells){
      state_by_t = c(state_by_t, sum(cond_cells %in% active_cells)/num_cells)
    }
    else{
      state_by_t = c(state_by_t, 0)
    }
  }
  state_by_m_spec = rbind(state_by_m_spec, state_by_t)
  
  state_by_t = c()
  for (m_id in unique(as.character(mc_bd_ann_tum$mouse_id))){
    cond_cells = rel_cells[as.character(md$mouse_id) == m_id & is_b16_endo & is_tumor]
    num_cells = length(cond_cells)
    if (num_cells > min_cells){
      state_by_t = c(state_by_t, sum(cond_cells %in% active_cells)/num_cells)
    }
    else{
      state_by_t = c(state_by_t, 0)
    }
  }
  state_by_m_b16 = rbind(state_by_m_b16, state_by_t)
  
  state_by_t = c()
  for (m_id in unique(as.character(mc_bd_ann_tum$mouse_id))){
    cond_cells = rel_cells[as.character(md$mouse_id) == m_id & is_mc38 & is_tumor]
    num_cells = length(cond_cells)
    if (num_cells > min_cells){
      state_by_t = c(state_by_t, sum(cond_cells %in% active_cells)/num_cells)
    }
    else{
      state_by_t = c(state_by_t, 0)
    }
  }
  state_by_m_mc38 = rbind(state_by_m_mc38, state_by_t)
  
  state_by_t = c()
  for (m_id in unique(as.character(mc_bd_ann_tum$mouse_id))){
    cond_cells = rel_cells[as.character(md$mouse_id) == m_id & is_byst & is_tumor]
    num_cells = length(cond_cells)
    if (num_cells > min_cells){
      state_by_t = c(state_by_t, sum(cond_cells %in% active_cells)/num_cells)
    }
    else{
      state_by_t = c(state_by_t, 0)
    }
  }
  state_by_m_byst = rbind(state_by_m_byst, state_by_t)
  
}
rownames(state_by_m_spec) = names(leaders_col)
colnames(state_by_m_spec) = unique(mc_bd_ann_tum$mouse_id)
rownames(state_by_m_mc38) = names(leaders_col)
colnames(state_by_m_mc38) = unique(mc_bd_ann_tum$mouse_id)
rownames(state_by_m_b16) = names(leaders_col)
colnames(state_by_m_b16) = unique(mc_bd_ann_tum$mouse_id)
rownames(state_by_m_byst) = names(leaders_col)
colnames(state_by_m_byst) = unique(mc_bd_ann_tum$mouse_id)

state_by_m_spec = state_by_m_spec[,colSums(state_by_m_spec) > 0]
state_by_m_mc38 = state_by_m_mc38[,colSums(state_by_m_mc38) > 0]
state_by_m_b16 = state_by_m_b16[,colSums(state_by_m_b16) > 0]
state_by_m_byst = state_by_m_byst[,colSums(state_by_m_byst) > 0]

state_by_m_ct_ls = list()
state_by_m_ct_ls$state_by_m_spec = state_by_m_spec
state_by_m_ct_ls$state_by_m_mc38 = state_by_m_mc38
state_by_m_ct_ls$state_by_m_b16 = state_by_m_b16
state_by_m_ct_ls$state_by_m_byst = state_by_m_byst
save(state_by_m_ct_ls, file="state_by_m_ct_ls")

m_id_time = table(md$mouse_id, md$days_post_transfer)

mc_bd_ann_tum_ctrl = mc_bd_ann_tum[as.character(mc_bd_ann_tum$condition) == "ctrl",]
state_by_m_spec_ctrl = state_by_m_spec[,colnames(state_by_m_spec) %in% mc_bd_ann_tum_ctrl$mouse_id]
state_by_m_mc38_ctrl = state_by_m_mc38[,colnames(state_by_m_mc38) %in% mc_bd_ann_tum_ctrl$mouse_id]
state_by_m_b16_ctrl = state_by_m_b16[,colnames(state_by_m_b16) %in% mc_bd_ann_tum_ctrl$mouse_id]
state_by_m_byst_ctrl = state_by_m_byst[,colnames(state_by_m_byst) %in% mc_bd_ann_tum_ctrl$mouse_id]




#### Tumor
#### Spec Ctrl
state_by_m_spec_ctrl = state_by_m_spec_ctrl[,colSums(state_by_m_spec_ctrl) > 0]
mc_bd_ann_tum_ctrl_spec = mc_bd_ann_tum_ctrl[as.character(mc_bd_ann_tum_ctrl$cell_type) == "spec",]

t_points_scat = t_points[1:8]

t_points_scat_late = t_points[11:14]


ylims = c(0.8,0.7,0.8,0.3,0.5,0.7,0.7,0.1, 0.6,0.6,0.5)
time_genes = c("Tcf7", "Sell", "Tnfrsf9", "Prf1", "Pdcd1", "Cxcr6", "Irf8", "Zbtb2", "Id2", "Id3", "Tigit", "Myc", "Lag3", "Havcr2", "Prdm1", "Il2ra", "Il2rb", "Ccr7")
for (prog in 1:nrow(state_by_m_spec_ctrl)){
  x_vec_late = c()
  y_vec_late = c()
  for (t in t_points_scat_late){
    m_ids = rownames(m_id_time)[m_id_time[,as.character(t)] > 0]
    prog_acts = state_by_m_spec_ctrl[prog,colnames(state_by_m_spec_ctrl) %in% m_ids]
    x_vec_late = c(x_vec_late, rep(t, length(prog_acts)))
    y_vec_late = c(y_vec_late, prog_acts)
  }
  
  x_vec = c()
  y_vec = c()
  for (t in t_points_scat){
    m_ids = rownames(m_id_time)[m_id_time[,as.character(t)] > 0]
    prog_acts = state_by_m_spec_ctrl[prog,colnames(state_by_m_spec_ctrl) %in% m_ids]
    x_vec = c(x_vec, rep(t, length(prog_acts)))
    y_vec = c(y_vec, prog_acts)
  }
  early_late_df = data.frame(act=c(y_vec_late, y_vec), time=c(rep("late", length(y_vec_late)), rep("early", length(y_vec))))
  colord = c("grey", "brown")
  names(colord) = c("early", "late")
  
  p_spec = ggviolin(early_late_df, x="time", y="act", fill="time", add="jitter") + scale_fill_manual(values=colord, labels=names(colord)) + scale_x_discrete(limits=c("early", "late"))  + ylab(paste0("%_active_",names(leaders_col)[prog]))
  p_spec = p_spec + theme(legend.position="none", text = element_text(size=20)) + ylim(0, ylims[prog]) 
  ggsave(filename = paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/violin/", "early_late_ctrl_",prog, ".png",sep=""), plot=p_spec)
  
  y_vec_ord = order(y_vec)
  y.lo = loess(y_vec ~ x_vec, span=1.2)
  pred_y = predict(y.lo, se=T)
  x_vec_ord = order(x_vec)
  x_vec_r = x_vec + runif(length(x_vec), min = -0.2, max=0.2)
  
  png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_conds_m/", "scatter_all_time_",prog, "_","spec_ctrl", ".png",sep=""), 8 * 100, 8 * 100, res = 200)
  plot(x_vec_r, y_vec, pch=19, xlab="time point", ylab=paste0("%_active_",names(leaders_col)[prog]), col=leaders_col[prog], ylim=c(0, ylims[prog]), xlim=c(0,9), cex=0.7)
  polygon(x = c(x_vec[x_vec_ord], rev(x_vec[x_vec_ord])), y = c(pred_y$fit[x_vec_ord] - qt(0.975,pred_y$df)*pred_y$se[x_vec_ord], rev(pred_y$fit[x_vec_ord] + qt(0.975,pred_y$df)*pred_y$se[x_vec_ord])), col = adjustcolor(leaders_col[prog], alpha.f=0.7), border=NA)
  lines(x_vec, pred_y$fit[x_vec_ord], col="grey", lwd=2, lty=2)
  
  dev.off()
}

t_points_scat = t_points[2:8]
downs_tum_spec = downs_mat[,md[colnames(downs_mat),"location"] == "tumor" & md[colnames(downs_mat),"cell_type"] %in% c("spec", "cd8_in_vitro") & as.character(md[colnames(downs_mat),"condition"]) == "ctrl"]
downs_tum_spec_genes = downs_tum_spec[time_genes,]
for (prog in time_genes){
  x_vec = c()
  y_vec = c()
  for (t in t_points_scat){
    m_ids = rownames(m_id_time)[m_id_time[,as.character(t)] > 0]
    # print(m_ids)
    # print(length(m_ids))
    downs_tum_spec_t = downs_tum_spec_genes[prog,as.character(mat@cell_metadata[colnames(downs_tum_spec_genes),"mouse_id"]) %in% m_ids, drop=F]
    prog_acts = tgs_matrix_tapply(downs_tum_spec_t, as.character(mat@cell_metadata[colnames(downs_tum_spec_t),"mouse_id"]), mean)
    x_vec = c(x_vec, rep(t, length(prog_acts)))
    y_vec = c(y_vec, prog_acts)
  }
  y_vec_ord = order(y_vec)
  y.lo = loess(y_vec ~ x_vec, span=1.2)
  pred_y = predict(y.lo, se=T)
  x_vec_ord = order(x_vec)
  x_vec_r = x_vec + runif(length(x_vec), min = -0.2, max=0.2)
  
  png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_conds_m/", "scatter_all_time_",prog, "_","spec_ctrl", ".png",sep=""), 8 * 100, 8 * 100, res = 200)
  plot(x_vec_r, y_vec, pch=19, xlab="time point", ylab=paste0("umi_",prog), ylim=c(0, max(y_vec)), xlim=c(0,9), cex=0.7)
  polygon(x = c(x_vec[x_vec_ord], rev(x_vec[x_vec_ord])), y = c(pred_y$fit[x_vec_ord] - qt(0.975,pred_y$df)*pred_y$se[x_vec_ord], rev(pred_y$fit[x_vec_ord] + qt(0.975,pred_y$df)*pred_y$se[x_vec_ord])), col = adjustcolor("grey", alpha.f=0.7), border=NA)
  lines(x_vec, pred_y$fit[x_vec_ord], col="grey", lwd=2, lty=2)
  
  dev.off()
}

state_by_m_spec_ctrl = state_by_m_spec_ctrl[,colSums(state_by_m_spec_ctrl) > 0]
mc_bd_ann_tum_ctrl_spec = mc_bd_ann_tum_ctrl[as.character(mc_bd_ann_tum_ctrl$cell_type) == "spec",]
  
#### B16 Ctrl
state_by_m_b16_ctrl = state_by_m_b16_ctrl[,colSums(state_by_m_b16_ctrl) > 0]
mc_bd_ann_tum_ctrl_b16 = mc_bd_ann_tum_ctrl[as.character(mc_bd_ann_tum_ctrl$cell_type) == "endo",]

t_points_scat = t_points[1:8]


ylims = c(0.8,0.7,0.8,0.5,0.5,0.7,0.7,0.1, 0.6,0.6,0.5)
for (prog in 1:nrow(state_by_m_spec_ctrl)){
  x_vec = c()
  y_vec = c()
  for (t in t_points_scat){
    m_ids = rownames(m_id_time)[m_id_time[,as.character(t)] > 0]
    prog_acts = state_by_m_b16_ctrl[prog,colnames(state_by_m_b16_ctrl) %in% m_ids]
    x_vec = c(x_vec, rep(t, length(prog_acts)))
    #x_vec = x_vec + runif(length(x_vec), min = -0.2, max=0.2)
    y_vec = c(y_vec, prog_acts)
  }
  y_vec_ord = order(y_vec)
  y.lo = loess(y_vec ~ x_vec, span=1.2)
  pred_y = predict(y.lo, se=T)
  x_vec_ord = order(x_vec)
  x_vec_r = x_vec + runif(length(x_vec), min = -0.2, max=0.2)
  png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_conds_m/", "scatter_all_time_",prog, "_","endo_ctrl", ".png",sep=""), 8 * 100, 8 * 100, res = 200)
  plot(x_vec_r, y_vec, pch=19, xlab="time point", ylab=paste0("%_active_",names(leaders_col)[prog]), col=leaders_col[prog], ylim=c(0, ylims[prog]), xlim=c(0,9), cex=0.7)
  polygon(x = c(x_vec[x_vec_ord], rev(x_vec[x_vec_ord])), y = c(pred_y$fit[x_vec_ord] - qt(0.975,pred_y$df)*pred_y$se[x_vec_ord], rev(pred_y$fit[x_vec_ord] + qt(0.975,pred_y$df)*pred_y$se[x_vec_ord])), col = adjustcolor(leaders_col[prog], alpha.f=0.7), border=NA)
  lines(x_vec, pred_y$fit[x_vec_ord], col="grey", lwd=2, lty=2)
  dev.off()
}

#### MC38 Ctrl
state_by_m_mc38_ctrl = state_by_m_mc38_ctrl[,colSums(state_by_m_mc38_ctrl) > 0]

t_points_scat = t_points[1:8]

ylims = c(0.8,0.7,0.8,0.5,0.5,0.7,0.7,0.1, 0.6,0.6,0.5)
for (prog in 1:nrow(state_by_m_spec_ctrl)){
  x_vec = c()
  y_vec = c()
  for (t in t_points_scat){
    m_ids = rownames(m_id_time)[m_id_time[,as.character(t)] > 0]
    prog_acts = state_by_m_mc38_ctrl[prog,colnames(state_by_m_mc38_ctrl) %in% m_ids]
    x_vec = c(x_vec, rep(t, length(prog_acts)))
    #x_vec = x_vec + runif(length(x_vec), min = -0.2, max=0.2)
    y_vec = c(y_vec, prog_acts)
  }
  y_vec_ord = order(y_vec)
  y.lo = loess(y_vec ~ x_vec, span=1.2)
  pred_y = predict(y.lo, se=T)
  x_vec_ord = order(x_vec)
  x_vec_r = x_vec + runif(length(x_vec), min = -0.2, max=0.2)
  png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_conds_m/", "scatter_all_time_",prog, "_","mc38_ctrl", ".png",sep=""), 8 * 100, 8 * 100, res = 200)
  plot(x_vec_r, y_vec, pch=19, xlab="time point", ylab=paste0("%_active_",names(leaders_col)[prog]), col=leaders_col[prog], ylim=c(0, ylims[prog]), xlim=c(0,9), cex=0.7)
  polygon(x = c(x_vec[x_vec_ord], rev(x_vec[x_vec_ord])), y = c(pred_y$fit[x_vec_ord] - qt(0.975,pred_y$df)*pred_y$se[x_vec_ord], rev(pred_y$fit[x_vec_ord] + qt(0.975,pred_y$df)*pred_y$se[x_vec_ord])), col = adjustcolor(leaders_col[prog], alpha.f=0.7), border=NA)
  lines(x_vec, pred_y$fit[x_vec_ord], col="grey", lwd=2, lty=2)
  dev.off()
}

#### byst Ctrl
state_by_m_byst_ctrl = state_by_m_byst_ctrl[,colSums(state_by_m_byst_ctrl) > 0]
mc_bd_ann_tum_ctrl_byst = mc_bd_ann_tum_ctrl[as.character(mc_bd_ann_tum_ctrl$cell_type) == "byst",]

t_points_scat = t_points[1:8]


ylims = c(0.8,0.7,0.8,0.5,0.5,0.7,0.7,0.1, 0.6,0.6,0.5)
for (prog in 1:nrow(state_by_m_spec_ctrl)){
  x_vec = c()
  y_vec = c()
  for (t in t_points_scat){
    m_ids = rownames(m_id_time)[m_id_time[,as.character(t)] > 0]
    prog_acts = state_by_m_byst_ctrl[prog,colnames(state_by_m_byst_ctrl) %in% m_ids]
    x_vec = c(x_vec, rep(t, length(prog_acts)))
    #x_vec = x_vec + runif(length(x_vec), min = -0.2, max=0.2)
    y_vec = c(y_vec, prog_acts)
  }
  y_vec_ord = order(y_vec)
  y.lo = loess(y_vec ~ x_vec, span=1.2)
  pred_y = predict(y.lo, se=T)
  x_vec_ord = order(x_vec)
  x_vec_r = x_vec + runif(length(x_vec), min = -0.2, max=0.2)
  png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_conds_m/", "scatter_all_time_",prog, "_","byst_ctrl", ".png",sep=""), 8 * 100, 8 * 100, res = 200)
  plot(x_vec_r, y_vec, pch=19, xlab="time point", ylab=paste0("%_active_",names(leaders_col)[prog]), col=leaders_col[prog], ylim=c(0, ylims[prog]), xlim=c(0,9), cex=0.7)
  polygon(x = c(x_vec[x_vec_ord], rev(x_vec[x_vec_ord])), y = c(pred_y$fit[x_vec_ord] - qt(0.975,pred_y$df)*pred_y$se[x_vec_ord], rev(pred_y$fit[x_vec_ord] + qt(0.975,pred_y$df)*pred_y$se[x_vec_ord])), col = adjustcolor(leaders_col[prog], alpha.f=0.7), border=NA)
  lines(x_vec, pred_y$fit[x_vec_ord], col="grey", lwd=2, lty=2)
  dev.off()
}

mc_bd_ann_tum_pd1 = mc_bd_ann_tum[as.character(mc_bd_ann_tum$condition) == "pd1",]
state_by_m_spec_pd1 = state_by_m_spec[,colnames(state_by_m_spec) %in% mc_bd_ann_tum_pd1$mouse_id]
state_by_m_mc38_pd1 = state_by_m_mc38[,colnames(state_by_m_mc38) %in% mc_bd_ann_tum_pd1$mouse_id]
state_by_m_b16_pd1 = state_by_m_b16[,colnames(state_by_m_b16) %in% mc_bd_ann_tum_pd1$mouse_id]

#### Spec PD1
state_by_m_spec_pd1 = state_by_m_spec_pd1[,colSums(state_by_m_spec_pd1) > 0]
mc_bd_ann_tum_pd1_spec = mc_bd_ann_tum_pd1[as.character(mc_bd_ann_tum_pd1$cell_type) == "spec",]

for (prog in 1:nrow(state_by_m_spec_pd1)){
  x_vec = c()
  y_vec = c()
  for (t in t_points_scat){
    m_ids = rownames(m_id_time)[m_id_time[,as.character(t)] > 0]
    prog_acts = state_by_m_spec_pd1[prog,colnames(state_by_m_spec_pd1) %in% m_ids]
    x_vec = c(x_vec, rep(t, length(prog_acts)))
    #x_vec = x_vec + runif(length(x_vec), min = -0.2, max=0.2)
    y_vec = c(y_vec, prog_acts)
  }
  y_vec_ord = order(y_vec)
  y.lo = loess(y_vec ~ x_vec, span=1.2)
  pred_y = predict(y.lo, se=T)
  x_vec_ord = order(x_vec)
  x_vec_r = x_vec + runif(length(x_vec), min = -0.2, max=0.2)
  png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_conds_m/", "scatter_all_time_",prog, "_","spec_pd1", ".png",sep=""), 8 * 100, 8 * 100, res = 200)
  plot(x_vec_r, y_vec, pch=19, xlab="time point", ylab=paste0("%_active_",names(leaders_col)[prog]), col=leaders_col[prog], ylim=c(0, ylims[prog]), xlim=c(0,9), cex=0.7)
  polygon(x = c(x_vec[x_vec_ord], rev(x_vec[x_vec_ord])), y = c(pred_y$fit[x_vec_ord] - qt(0.975,pred_y$df)*pred_y$se[x_vec_ord], rev(pred_y$fit[x_vec_ord] + qt(0.975,pred_y$df)*pred_y$se[x_vec_ord])), col = adjustcolor(leaders_col[prog], alpha.f=0.7), border=NA)
  lines(x_vec, pred_y$fit[x_vec_ord], col="grey", lwd=2, lty=2)
  dev.off()
}

#### B16 PD1
state_by_m_b16_pd1 = state_by_m_b16_pd1[,colSums(state_by_m_b16_pd1) > 0]
mc_bd_ann_tum_pd1_b16 = mc_bd_ann_tum_pd1[as.character(mc_bd_ann_tum_pd1$cell_type) == "endo",]

for (prog in 1:nrow(state_by_m_spec_pd1)){
  x_vec = c()
  y_vec = c()
  for (t in t_points_scat){
    m_ids = rownames(m_id_time)[m_id_time[,as.character(t)] > 0]
    prog_acts = state_by_m_b16_pd1[prog,colnames(state_by_m_b16_pd1) %in% m_ids]
    x_vec = c(x_vec, rep(t, length(prog_acts)))
    #x_vec = x_vec + runif(length(x_vec), min = -0.2, max=0.2)
    y_vec = c(y_vec, prog_acts)
  }
  y_vec_ord = order(y_vec)
  y.lo = loess(y_vec ~ x_vec, span=1.2)
  pred_y = predict(y.lo, se=T)
  x_vec_ord = order(x_vec)
  x_vec_r = x_vec + runif(length(x_vec), min = -0.2, max=0.2)
  png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_conds_m/", "scatter_all_time_",prog, "_","endo_pd1", ".png",sep=""), 8 * 100, 8 * 100, res = 200)
  plot(x_vec_r, y_vec, pch=19, xlab="time point", ylab=paste0("%_active_",names(leaders_col)[prog]), col=leaders_col[prog], ylim=c(0, ylims[prog]), xlim=c(0,9), cex=0.7)
  polygon(x = c(x_vec[x_vec_ord], rev(x_vec[x_vec_ord])), y = c(pred_y$fit[x_vec_ord] - qt(0.975,pred_y$df)*pred_y$se[x_vec_ord], rev(pred_y$fit[x_vec_ord] + qt(0.975,pred_y$df)*pred_y$se[x_vec_ord])), col = adjustcolor(leaders_col[prog], alpha.f=0.7), border=NA)
  lines(x_vec, pred_y$fit[x_vec_ord], col="grey", lwd=2, lty=2)
  dev.off()
}

#### MC38 PD1
state_by_m_mc38_pd1 = state_by_m_mc38_pd1[,colSums(state_by_m_mc38_pd1) > 0]

for (prog in 1:nrow(state_by_m_spec_pd1)){
  x_vec = c()
  y_vec = c()
  for (t in t_points_scat){
    m_ids = rownames(m_id_time)[m_id_time[,as.character(t)] > 0]
    prog_acts = state_by_m_mc38_pd1[prog,colnames(state_by_m_mc38_pd1) %in% m_ids]
    x_vec = c(x_vec, rep(t, length(prog_acts)))
    #x_vec = x_vec + runif(length(x_vec), min = -0.2, max=0.2)
    y_vec = c(y_vec, prog_acts)
  }
  y_vec_ord = order(y_vec)
  y.lo = loess(y_vec ~ x_vec, span=1.2)
  pred_y = predict(y.lo, se=T)
  x_vec_ord = order(x_vec)
  x_vec_r = x_vec + runif(length(x_vec), min = -0.2, max=0.2)
  png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_conds_m/", "scatter_all_time_",prog, "_","mc38_pd1", ".png",sep=""), 8 * 100, 8 * 100, res = 200)
  plot(x_vec_r, y_vec, pch=19, xlab="time point", ylab=paste0("%_active_",names(leaders_col)[prog]), col=leaders_col[prog], ylim=c(0, ylims[prog]), xlim=c(0,9), cex=0.7)
  polygon(x = c(x_vec[x_vec_ord], rev(x_vec[x_vec_ord])), y = c(pred_y$fit[x_vec_ord] - qt(0.975,pred_y$df)*pred_y$se[x_vec_ord], rev(pred_y$fit[x_vec_ord] + qt(0.975,pred_y$df)*pred_y$se[x_vec_ord])), col = adjustcolor(leaders_col[prog], alpha.f=0.7), border=NA)
  lines(x_vec, pred_y$fit[x_vec_ord], col="grey", lwd=2, lty=2)
  dev.off()
}


#### Non Tumor
mc_bd_ann = as.data.frame(matrix(unlist(strsplit(colnames(mc_bd_n), split=" ")), ncol=length(group_mc_by), byrow=T))
rownames(mc_bd_ann) = colnames(mc_bd)
colnames(mc_bd_ann) = group_mc_by

mc_bd_ann$mouse_id = as.character(mc_bd_ann$mouse_id)

mc_bd_ann_peri = mc_bd_ann[mc_bd_ann$location == "LN",]
# levels(mc_bd_ann_tum$cell_type) = c(levels(mc_bd_ann_tum$cell_type), "endo_mc38")
# mc_bd_ann_tum$cell_type[mc_bd_ann_tum$mouse_id %in% mc_bd_ann_tum_endo_mc38$mouse_id] = "endo_mc38"
mc_bd_ann_peri = mc_bd_ann_peri[order(mc_bd_ann_peri$cell_type),]
mc_bd_ann_peri = mc_bd_ann_peri[!mc_bd_ann_peri$condition == "na",]

m_id_cond = table(md$mouse_id, md$condition)
m_id_cond_peri = m_id_cond[as.character(unique(mc_bd_ann_peri$mouse_id)),]

mc_bd_ann_peri_3 = mc_bd_ann_peri[mc_bd_ann_peri$days_post_transfer == 3,]
mc_bd_ann_peri_6 = mc_bd_ann_peri[mc_bd_ann_peri$days_post_transfer == 6,]

mc_bd_ann_peri_spec = mc_bd_ann_peri[mc_bd_ann_peri$cell_type == "spec",]

mc_bd_ann_peri_endo_mc38 = mc_bd_ann_peri[mc_bd_ann_peri$cell_type == "endo_mc38",]
mc_bd_ann_peri_endo_b16 = mc_bd_ann_peri[mc_bd_ann_peri$cell_type == "endo",]



state_by_m_spec_peri = c()
state_by_m_mc38_peri = c()
state_by_m_b16_peri = c()
state_by_m_byst_peri = c()

for (prog in 1:ncol(norm_global_mat)){
  active_thres = active_thres_vec[prog]
  strong_thres = strong_thres_vec[prog]
  active_mcs = (1:max(mc@mc))[norm_global_mat[,prog] >= active_thres]
  active_cells = rel_cells[mc@mc[rel_cells] %in% active_mcs]
  
  min_cells = 30
  state_by_t = c()
  for (m_id in unique(as.character(mc_bd_ann_peri$mouse_id))){
    cond_cells = rel_cells[as.character(md$mouse_id) == m_id & is_spec & is_ln]
    num_cells = length(cond_cells)
    if (num_cells > min_cells){
      state_by_t = c(state_by_t, sum(cond_cells %in% active_cells)/num_cells)
    }
    else{
      state_by_t = c(state_by_t, 0)
    }
  }
  state_by_m_spec_peri = rbind(state_by_m_spec_peri, state_by_t)
  
  state_by_t = c()
  for (m_id in unique(as.character(mc_bd_ann_peri$mouse_id))){
    cond_cells = rel_cells[as.character(md$mouse_id) == m_id & is_b16_endo & is_ln]
    num_cells = length(cond_cells)
    if (num_cells > min_cells){
      state_by_t = c(state_by_t, sum(cond_cells %in% active_cells)/num_cells)
    }
    else{
      state_by_t = c(state_by_t, 0)
    }
  }
  state_by_m_b16_peri = rbind(state_by_m_b16_peri, state_by_t)
  
  state_by_t = c()
  for (m_id in unique(as.character(mc_bd_ann_peri$mouse_id))){
    cond_cells = rel_cells[as.character(md$mouse_id) == m_id & is_mc38 & is_ln]
    num_cells = length(cond_cells)
    if (num_cells > min_cells){
      state_by_t = c(state_by_t, sum(cond_cells %in% active_cells)/num_cells)
    }
    else{
      state_by_t = c(state_by_t, 0)
    }
  }
  state_by_m_mc38_peri = rbind(state_by_m_mc38_peri, state_by_t)
  
  state_by_t = c()
  for (m_id in unique(as.character(mc_bd_ann_peri$mouse_id))){
    cond_cells = rel_cells[as.character(md$mouse_id) == m_id & is_byst & is_ln]
    num_cells = length(cond_cells)
    if (num_cells > min_cells){
      state_by_t = c(state_by_t, sum(cond_cells %in% active_cells)/num_cells)
    }
    else{
      state_by_t = c(state_by_t, 0)
    }
  }
  state_by_m_byst_peri = rbind(state_by_m_byst_peri, state_by_t)
  
}
rownames(state_by_m_spec_peri) = names(leaders_col)
colnames(state_by_m_spec_peri) = unique(mc_bd_ann_peri$mouse_id)
rownames(state_by_m_mc38_peri) = names(leaders_col)
colnames(state_by_m_mc38_peri) = unique(mc_bd_ann_peri$mouse_id)
rownames(state_by_m_b16_peri) = names(leaders_col)
colnames(state_by_m_b16_peri) = unique(mc_bd_ann_peri$mouse_id)
rownames(state_by_m_byst_peri) = names(leaders_col)
colnames(state_by_m_byst_peri) = unique(mc_bd_ann_peri$mouse_id)

state_by_m_spec_peri = state_by_m_spec_peri[,colSums(state_by_m_spec_peri) > 0]
state_by_m_mc38_peri = state_by_m_mc38_peri[,colSums(state_by_m_mc38_peri) > 0]
state_by_m_b16_peri = state_by_m_b16_peri[,colSums(state_by_m_b16_peri) > 0]
state_by_m_byst_peri = state_by_m_b16_peri[,colSums(state_by_m_b16_peri) > 0]

m_id_time = table(md$mouse_id, md$days_post_transfer)

mc_bd_ann_peri_ctrl = mc_bd_ann_peri[as.character(mc_bd_ann_peri$condition) == "ctrl",]
state_by_m_spec_peri_ctrl = state_by_m_spec_peri[,colnames(state_by_m_spec_peri) %in% mc_bd_ann_peri_ctrl$mouse_id]
state_by_m_mc38_peri_ctrl = state_by_m_mc38_peri[,colnames(state_by_m_mc38_peri) %in% mc_bd_ann_peri_ctrl$mouse_id]
state_by_m_b16_peri_ctrl = state_by_m_b16_peri[,colnames(state_by_m_b16_peri) %in% mc_bd_ann_peri_ctrl$mouse_id]
state_by_m_byst_peri_ctrl = state_by_m_byst_peri[,colnames(state_by_m_byst_peri) %in% mc_bd_ann_peri_ctrl$mouse_id]


#### Periphery
#### Spec Ctrl
state_by_m_spec_peri_ctrl = state_by_m_spec_peri_ctrl[,colSums(state_by_m_spec_peri_ctrl) > 0]
mc_bd_ann_peri_ctrl_spec = mc_bd_ann_peri_ctrl[as.character(mc_bd_ann_peri_ctrl$cell_type) == "spec",]

t_points_scat = t_points[1:8]

ylims = c(0.8,1,0.8,0.5,0.5,0.7,0.7,0.1, 0.9,0.5,0.5)
for (prog in 1:nrow(state_by_m_spec_peri_ctrl)){
  x_vec = c()
  y_vec = c()
  for (t in t_points_scat){
    m_ids = rownames(m_id_time)[m_id_time[,as.character(t)] > 0]
    prog_acts = state_by_m_spec_peri_ctrl[prog,colnames(state_by_m_spec_peri_ctrl) %in% m_ids]
    x_vec = c(x_vec, rep(t, length(prog_acts)))
    #x_vec = x_vec + runif(length(x_vec), min = -0.2, max=0.2)
    y_vec = c(y_vec, prog_acts)
  }
  y_vec_ord = order(y_vec)
  y.lo = loess(y_vec ~ x_vec, span=1.2)
  pred_y = predict(y.lo, se=T)
  x_vec_ord = order(x_vec)
  x_vec_r = x_vec + runif(length(x_vec), min = -0.2, max=0.2)
  png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_conds_m/", "scatter_all_time_",prog, "_","peri_spec_ctrl", ".png",sep=""), 8 * 100, 8 * 100, res = 200)
  plot(x_vec_r, y_vec, pch=19, xlab="time point", ylab=paste0("%_active_",names(leaders_col)[prog]), col=leaders_col[prog], ylim=c(0, ylims[prog]), xlim=c(0,9), cex=0.7)
  polygon(x = c(x_vec[x_vec_ord], rev(x_vec[x_vec_ord])), y = c(pred_y$fit[x_vec_ord] - qt(0.975,pred_y$df)*pred_y$se[x_vec_ord], rev(pred_y$fit[x_vec_ord] + qt(0.975,pred_y$df)*pred_y$se[x_vec_ord])), col = adjustcolor(leaders_col[prog], alpha.f=0.7), border=NA)
  lines(x_vec, pred_y$fit[x_vec_ord], col="grey", lwd=2, lty=2)
  dev.off()
}

#### B16 Ctrl
state_by_m_b16_peri_ctrl = state_by_m_b16_peri_ctrl[,colSums(state_by_m_b16_peri_ctrl) > 0]
mc_bd_ann_peri_ctrl_b16 = mc_bd_ann_peri_ctrl[as.character(mc_bd_ann_peri_ctrl$cell_type) == "endo",]

t_points_scat = t_points[1:8]
ylims = c(0.8,1,0.8,0.5,0.5,0.7,0.7,0.1, 0.9,0.5,0.5)
for (prog in 1:nrow(state_by_m_spec_peri_ctrl)){
  x_vec = c()
  y_vec = c()
  for (t in t_points_scat){
    m_ids = rownames(m_id_time)[m_id_time[,as.character(t)] > 0]
    prog_acts = state_by_m_b16_peri_ctrl[prog,colnames(state_by_m_b16_peri_ctrl) %in% m_ids]
    x_vec = c(x_vec, rep(t, length(prog_acts)))
    #x_vec = x_vec + runif(length(x_vec), min = -0.2, max=0.2)
    y_vec = c(y_vec, prog_acts)
  }
  y_vec_ord = order(y_vec)
  y.lo = loess(y_vec ~ x_vec, span=1.2)
  pred_y = predict(y.lo, se=T)
  x_vec_ord = order(x_vec)
  x_vec_r = x_vec + runif(length(x_vec), min = -0.2, max=0.2)
  png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_conds_m/", "scatter_all_time_",prog, "_","peri_endo_ctrl", ".png",sep=""), 8 * 100, 8 * 100, res = 200)
  plot(x_vec_r, y_vec, pch=19, xlab="time point", ylab=paste0("%_active_",names(leaders_col)[prog]), col=leaders_col[prog], ylim=c(0, ylims[prog]), xlim=c(0,9), cex=0.7)
  # polygon(x = c(x_vec[x_vec_ord], rev(x_vec[x_vec_ord])), y = c(pred_y$fit[x_vec_ord] - qt(0.975,pred_y$df)*pred_y$se[x_vec_ord], rev(pred_y$fit[x_vec_ord] + qt(0.975,pred_y$df)*pred_y$se[x_vec_ord])), col = adjustcolor(leaders_col[prog], alpha.f=0.7), border=NA)
  polygon(x = c(x_vec[x_vec_ord], rev(x_vec[x_vec_ord])), y = c(pred_y$fit[x_vec_ord] - 0.02, rev(pred_y$fit[x_vec_ord] + 0.02)), col = adjustcolor(leaders_col[prog], alpha.f=0.7), border=NA)
  lines(x_vec, pred_y$fit[x_vec_ord], col="grey", lwd=2, lty=2)
  dev.off()
}

#### byst Ctrl
state_by_m_byst_peri_ctrl = state_by_m_byst_peri_ctrl[,colSums(state_by_m_byst_peri_ctrl) > 0]
mc_bd_ann_peri_ctrl_byst = mc_bd_ann_peri_ctrl[as.character(mc_bd_ann_peri_ctrl$cell_type) == "byst",]

t_points_scat = t_points[1:8]
ylims = c(0.8,1,0.8,0.5,0.5,0.7,0.7,0.1, 0.9,0.5,0.5)
for (prog in 1:nrow(state_by_m_spec_peri_ctrl)){
  x_vec = c()
  y_vec = c()
  for (t in t_points_scat){
    m_ids = rownames(m_id_time)[m_id_time[,as.character(t)] > 0]
    prog_acts = state_by_m_byst_peri_ctrl[prog,colnames(state_by_m_byst_peri_ctrl) %in% m_ids]
    x_vec = c(x_vec, rep(t, length(prog_acts)))
    #x_vec = x_vec + runif(length(x_vec), min = -0.2, max=0.2)
    y_vec = c(y_vec, prog_acts)
  }
  y_vec_ord = order(y_vec)
  y.lo = loess(y_vec ~ x_vec, span=1.2)
  pred_y = predict(y.lo, se=T)
  x_vec_ord = order(x_vec)
  x_vec_r = x_vec + runif(length(x_vec), min = -0.2, max=0.2)
  png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_conds_m/", "scatter_all_time_",prog, "_","peri_byst_ctrl", ".png",sep=""), 8 * 100, 8 * 100, res = 200)
  plot(x_vec_r, y_vec, pch=19, xlab="time point", ylab=paste0("%_active_",names(leaders_col)[prog]), col=leaders_col[prog], ylim=c(0, ylims[prog]), xlim=c(0,9), cex=0.7)
  # polygon(x = c(x_vec[x_vec_ord], rev(x_vec[x_vec_ord])), y = c(pred_y$fit[x_vec_ord] - qt(0.975,pred_y$df)*pred_y$se[x_vec_ord], rev(pred_y$fit[x_vec_ord] + qt(0.975,pred_y$df)*pred_y$se[x_vec_ord])), col = adjustcolor(leaders_col[prog], alpha.f=0.7), border=NA)
  polygon(x = c(x_vec[x_vec_ord], rev(x_vec[x_vec_ord])), y = c(pred_y$fit[x_vec_ord] - 0.02, rev(pred_y$fit[x_vec_ord] + 0.02)), col = adjustcolor(leaders_col[prog], alpha.f=0.7), border=NA)
  lines(x_vec, pred_y$fit[x_vec_ord], col="grey", lwd=2, lty=2)
  dev.off()
}


mc_bd_ann_peri_pd1 = mc_bd_ann_peri[as.character(mc_bd_ann_peri$condition) == "pd1",]
state_by_m_spec_peri_pd1 = state_by_m_spec_peri[,colnames(state_by_m_spec_peri) %in% mc_bd_ann_peri_pd1$mouse_id]
state_by_m_b16_peri_pd1 = state_by_m_b16_peri[,colnames(state_by_m_b16_peri) %in% mc_bd_ann_peri_pd1$mouse_id]

#### Spec PD1
state_by_m_spec_peri_pd1 = state_by_m_spec_peri_pd1[,colSums(state_by_m_spec_peri_pd1) > 0]
mc_bd_ann_peri_pd1_spec = mc_bd_ann_peri_pd1[as.character(mc_bd_ann_peri_pd1$cell_type) == "spec",]

for (prog in 1:nrow(state_by_m_spec_pd1)){
  x_vec = c()
  y_vec = c()
  for (t in t_points_scat){
    m_ids = rownames(m_id_time)[m_id_time[,as.character(t)] > 0]
    prog_acts = state_by_m_spec_peri_pd1[prog,colnames(state_by_m_spec_peri_pd1) %in% m_ids]
    x_vec = c(x_vec, rep(t, length(prog_acts)))
    #x_vec = x_vec + runif(length(x_vec), min = -0.2, max=0.2)
    y_vec = c(y_vec, prog_acts)
  }
  y_vec_ord = order(y_vec)
  y.lo = loess(y_vec ~ x_vec, span=1.2)
  pred_y = predict(y.lo, se=T)
  x_vec_ord = order(x_vec)
  x_vec_r = x_vec + runif(length(x_vec), min = -0.2, max=0.2)
  png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_conds_m/", "scatter_all_time_",prog, "_","peri_spec_pd1", ".png",sep=""), 8 * 100, 8 * 100, res = 200)
  plot(x_vec_r, y_vec, pch=19, xlab="time point", ylab=paste0("%_active_",names(leaders_col)[prog]), col=leaders_col[prog], ylim=c(0, ylims[prog]), xlim=c(0,9), cex=0.7)
  polygon(x = c(x_vec[x_vec_ord], rev(x_vec[x_vec_ord])), y = c(pred_y$fit[x_vec_ord] - qt(0.975,pred_y$df)*pred_y$se[x_vec_ord], rev(pred_y$fit[x_vec_ord] + qt(0.975,pred_y$df)*pred_y$se[x_vec_ord])), col = adjustcolor(leaders_col[prog], alpha.f=0.7), border=NA)
  lines(x_vec, pred_y$fit[x_vec_ord], col="grey", lwd=2, lty=2)
  dev.off()
}

#### B16 PD1
state_by_m_b16_peri_pd1 = state_by_m_b16_peri_pd1[,colSums(state_by_m_b16_peri_pd1) > 0]
mc_bd_ann_peri_pd1_b16 = mc_bd_ann_peri_pd1[as.character(mc_bd_ann_peri_pd1$cell_type) == "endo",]

for (prog in 1:nrow(state_by_m_spec_pd1)){
  x_vec = c()
  y_vec = c()
  for (t in t_points_scat){
    m_ids = rownames(m_id_time)[m_id_time[,as.character(t)] > 0]
    prog_acts = state_by_m_b16_peri_pd1[prog,colnames(state_by_m_b16_peri_pd1) %in% m_ids]
    x_vec = c(x_vec, rep(t, length(prog_acts)))
    #x_vec = x_vec + runif(length(x_vec), min = -0.2, max=0.2)
    y_vec = c(y_vec, prog_acts)
  }
  y_vec_ord = order(y_vec)
  y.lo = loess(y_vec ~ x_vec, span=1.2)
  pred_y = predict(y.lo, se=T)
  x_vec_ord = order(x_vec)
  x_vec_r = x_vec + runif(length(x_vec), min = -0.2, max=0.2)
  png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_conds_m/", "scatter_all_time_",prog, "_","peri_endo_pd1", ".png",sep=""), 8 * 100, 8 * 100, res = 200)
  plot(x_vec_r, y_vec, pch=19, xlab="time point", ylab=paste0("%_active_",names(leaders_col)[prog]), col=leaders_col[prog], ylim=c(0, ylims[prog]), xlim=c(0,9), cex=0.7)
  polygon(x = c(x_vec[x_vec_ord], rev(x_vec[x_vec_ord])), y = c(pred_y$fit[x_vec_ord] - qt(0.975,pred_y$df)*pred_y$se[x_vec_ord], rev(pred_y$fit[x_vec_ord] + qt(0.975,pred_y$df)*pred_y$se[x_vec_ord])), col = adjustcolor(leaders_col[prog], alpha.f=0.7), border=NA)
  lines(x_vec, pred_y$fit[x_vec_ord], col="grey", lwd=2, lty=2)
  dev.off()
}

state_by_m_spec_3 = state_by_m_spec[,colnames(state_by_m_spec) %in% mc_bd_ann_tum_3$mouse_id]
state_by_m_spec_6 = state_by_m_spec[,colnames(state_by_m_spec) %in% mc_bd_ann_tum_6$mouse_id]
state_by_m_mc38_3 = state_by_m_mc38[,colnames(state_by_m_mc38) %in% mc_bd_ann_tum_3$mouse_id]
state_by_m_mc38_6 = state_by_m_mc38[,colnames(state_by_m_mc38) %in% mc_bd_ann_tum_6$mouse_id]
state_by_m_b16_3 = state_by_m_b16[,colnames(state_by_m_b16) %in% mc_bd_ann_tum_3$mouse_id]
state_by_m_b16_6 = state_by_m_b16[,colnames(state_by_m_b16) %in% mc_bd_ann_tum_6$mouse_id]

cond2col =  c("grey", "black", "pink", "red")
names(cond2col) = c("ctrl", "pd1", "41bb", "41bb_pd1")

col2cond = c("ctrl", "pd1", "41bb", "41bb_pd1")
names(col2cond) =  c("grey", "black", "pink", "red")


mouse_cols_spec = rep("white", ncol(state_by_m_spec_3))
names(mouse_cols_spec) = colnames(state_by_m_spec_3)
mouse_cols_mc38 = rep("white", ncol(state_by_m_mc38_3))
names(mouse_cols_mc38) = colnames(state_by_m_mc38_3)
mouse_cols_b16 = rep("white", ncol(state_by_m_b16_3))
names(mouse_cols_b16) = colnames(state_by_m_b16_3)
for (cond in colnames(m_id_cond)){
  if (cond != "na"){
    m_ids = rownames(m_id_cond)[m_id_cond[,cond] > 0]
    mouse_cols_spec[colnames(state_by_m_spec_3) %in% m_ids] = cond2col[cond]
    
    mouse_cols_mc38[colnames(state_by_m_mc38_3) %in% m_ids] = cond2col[cond]
    
    mouse_cols_b16[colnames(state_by_m_b16_3) %in% m_ids] = cond2col[cond]
    
  }
}


state_by_m_spec_3 = state_by_m_spec_3[,order(mouse_cols_spec)]
state_by_m_spec_3_df = as.data.frame(t(state_by_m_spec_3))
state_by_m_spec_3_df$cond =  col2cond[mouse_cols_spec[colnames(state_by_m_spec_3)]]
state_by_m_mc38_3 = state_by_m_mc38_3[,order(mouse_cols_mc38)]
state_by_m_mc38_3_df = as.data.frame(t(state_by_m_mc38_3))
state_by_m_mc38_3_df$cond =  col2cond[mouse_cols_mc38[colnames(state_by_m_mc38_3)]]
state_by_m_b13_3 = state_by_m_b16_3[,order(mouse_cols_b16)]
state_by_m_b16_3_df = as.data.frame(t(state_by_m_b16_3))
state_by_m_b16_3_df$cond =  col2cond[mouse_cols_b16[colnames(state_by_m_b16_3)]]

dir.create(paste0(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_conds_m/"))

col2ind = c(1,2,3,4)
names(col2ind) = c("grey", "black", "pink", "red")


dir.create(paste0(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/violin_endo/"))
violin_dir = paste0(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/violin_endo/")
state_by_m_spec_3_df_pd1 = state_by_m_spec_3_df[as.character(state_by_m_spec_3_df$cond) %in% c("ctrl", "pd1"),]
state_by_m_mc38_3_df_pd1 = state_by_m_mc38_3_df[as.character(state_by_m_mc38_3_df$cond) %in% c("ctrl"),]
state_by_m_b16_3_df_pd1 = state_by_m_b16_3_df[as.character(state_by_m_b16_3_df$cond) %in% c("ctrl"),]
for (prog in 1:ncol(norm_global_w)){
  ylim_max = max(c(0.8, max(state_by_m_spec_3[prog,]), max(state_by_m_mc38_3[prog,]), max(state_by_m_b16_3[prog,]), max(state_by_m_spec_6[prog,]), max(state_by_m_mc38_6[prog,]), max(state_by_m_b16_6[prog,]))) + 0.1
  colord = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF")
  names(colord) = c("pd1", "ctrl", "41bb", "41bb_pd1")
  
  p_spec = ggviolin(state_by_m_spec_3_df, x="cond", y=colnames(state_by_m_spec_3_df)[prog], fill="cond", add="jitter", palette = "jco") + scale_fill_manual(values=colord, labels=names(colord)) + scale_x_discrete(limits=c("pd1", "ctrl", "41bb", "41bb_pd1")) 
  p_spec = p_spec + theme(legend.position="none", text = element_text(size=20)) + ylim(0, ylim_max) 
  ggsave(filename = paste( violin_dir,"time_3_violin_spec_41bb_pd1_",prog, ".png",sep=""), plot=p_spec, width = 2, height = 1, scale = 3)
  
  p_mc38 = ggviolin(state_by_m_mc38_3_df, x="cond", y=colnames(state_by_m_mc38_3_df)[prog], fill="cond", add="jitter", palette = "jco") + scale_fill_manual(values=colord, labels=names(colord)) + scale_x_discrete(limits=c("pd1", "ctrl", "41bb", "41bb_pd1")) 
  p_mc38 = p_mc38 + theme(legend.position="none", text = element_text(size=20)) + ylim(0, ylim_max)
  ggsave(filename = paste(violin_dir, "time_3_violin_mc38_41bb_pd1_",prog, ".png",sep=""), plot=p_mc38, width = 2, height = 1, scale = 3)
  
  p_b16 = ggviolin(state_by_m_b16_3_df, x="cond", y=colnames(state_by_m_b16_3_df)[prog], fill="cond", add="jitter", palette = "jco") + scale_fill_manual(values=colord, labels=names(colord)) + scale_x_discrete(limits=c("pd1", "ctrl", "41bb", "41bb_pd1")) 
  p_b16 = p_b16  + theme(legend.position="none", text = element_text(size=20)) + ylim(0, ylim_max)
  ggsave(filename = paste(violin_dir, "time_3_violin_b16_41bb_pd1_",prog, ".png",sep=""), plot=p_b16, width = 2, height = 1, scale = 3)
}

# Per cell type violins
for (prog in 1:ncol(norm_global_w)){
  ylim_max = max(c(0.7, max(state_by_m_spec_3[prog,]), max(state_by_m_mc38_3[prog,]), max(state_by_m_b16_3[prog,]), max(state_by_m_spec_6[prog,]), max(state_by_m_mc38_6[prog,]), max(state_by_m_b16_6[prog,]))) + 0.1
  p_spec = ggviolin(state_by_m_spec_3_df_pd1, x="cond", y=colnames(state_by_m_spec_3_df_pd1)[prog], fill="cond", add="jitter", palette = "jco")  + scale_x_discrete(limits=c("ctrl", "pd1"))
  p_spec = p_spec + theme(legend.position="none", text = element_text(size=30)) #+ stat_compare_means(comparisons=list(c("ctrl", "pd1")),method = "wilcox.test", label = "p.signif", bracket.size = 3, ) #ylim(0, ylim_max) + 
  ggsave(filename = paste(violin_dir, "time_3_violin_spec_pd1_",prog, ".png",sep=""), plot=p_spec)
  
  p_mc38 = ggviolin(state_by_m_mc38_3_df, x="cond", y=colnames(state_by_m_mc38_3_df)[prog], fill="cond", add="jitter", palette = "jco")  + scale_x_discrete(limits=c("ctrl", "pd1"))
  p_mc38 = p_mc38 + theme(legend.position="none", text = element_text(size=30)) + ylim(0, ylim_max) #+ stat_compare_means(comparisons=list(c("ctrl", "pd1")),method = "wilcox.test", label = "p.signif", bracket.size = 3) 
  ggsave(filename = paste(violin_dir, "time_3_violin_mc38_pd1_",prog, ".png",sep=""), plot=p_mc38)
  
  colord = c("#EFC000FF", "#0073C2FF")
  names(colord) = c("ctrl", "pd1")
  p_b16 = ggviolin(state_by_m_b16_3_df, x="cond", y=colnames(state_by_m_b16_3_df)[prog], fill="cond", add="jitter", palette = "jco") + scale_fill_manual(values=colord, labels=names(colord)) + scale_x_discrete(limits=c("ctrl", "pd1")) 
  p_b16 = p_b16 + theme(legend.position="none", text = element_text(size=30)) + ylim(0, ylim_max) #+ stat_compare_means(comparisons=list(c("ctrl", "pd1")),method = "wilcox.test", label = "p.signif", bracket.size = 3)
  ggsave(filename = paste(violin_dir, "time_3_violin_b16_pd1_",prog, ".png",sep=""), plot=p_b16)
}

# B16-MC38 plot
state_by_m_mc38_3_df_pd1$ct = "MC38"
state_by_m_mc38_3_df_pd1$ct = as.character(state_by_m_mc38_3_df_pd1$ct)
state_by_m_b16_3_df_pd1$ct = "B16"
state_by_m_b16_3_df_pd1$ct = as.character(state_by_m_b16_3_df_pd1$ct)
state_by_m_endo_3_df_ctrl = rbind(state_by_m_b16_3_df_pd1, state_by_m_mc38_3_df_pd1)

#state_by_m_endo_3_df_ctrl = state_by_m_endo_3_df_pd1[state_by_m_endo_3_df_pd1$cond == "ctrl",]
for (prog in 1:ncol(norm_global_w)){
  print(prog)
  ylim_max = max(state_by_m_endo_3_df_ctrl[,prog]) + 0.1
  
  p_endo = ggviolin(state_by_m_endo_3_df_ctrl, x="ct", y=colnames(state_by_m_endo_3_df_ctrl)[prog], fill="ct", add="jitter", palette = "jama")  + scale_x_discrete(limits=c("B16", "MC38")) + xlab("Tumor model")
  p_endo = p_endo + theme(legend.position="none", text = element_text(size=30)) + ylim(0, ylim_max) #+ stat_compare_means(comparisons=list(c("B16", "MC38")),method = "wilcox.test", label = "p.signif", bracket.size = 1, size = 8) 
  ggsave(filename = paste(violin_dir, "time_3_violin_endo_ctrl_",prog, ".png",sep=""), plot=p_endo)
  
}

for (prog in 1:ncol(norm_global_w)){
  png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_conds_m/", "scatter_active_",prog, "_","spec", "_time_", "3", ".png",sep=""), 8 * 100, 8 * 100, res = 200)
  ylim_max = max(c(0.4, max(state_by_m_spec_3[prog,]), max(state_by_m_mc38_3[prog,]), max(state_by_m_b16_3[prog,]), max(state_by_m_spec_6[prog,]), max(state_by_m_mc38_6[prog,]), max(state_by_m_b16_6[prog,])))
  plot(0.2*runif(ncol(state_by_m_spec_3)) + col2ind[mouse_cols_spec[colnames(state_by_m_spec_3)]], state_by_m_spec_3[prog,], col=mouse_cols_spec[colnames(state_by_m_spec_3)], pch=19, xaxt = "n", ylim=c(0,ylim_max), main=paste("spec 3 p=", names(leaders_col)[prog], sep=""), xlab="Cond", ylab="% Active")
  axis(1, at=1:4, labels=c("ctrl", "pd1", "41bb", "41bb_pd1"))
  dev.off()
  
  png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_conds_m/", "scatter_active_",prog, "_","mc38", "_time_", "3", ".png",sep=""), 8 * 100, 8 * 100, res = 200)
  plot(0.1*runif(ncol(state_by_m_mc38_3)) + col2ind[mouse_cols_mc38[colnames(state_by_m_mc38_3)]], state_by_m_mc38_3[prog,], col=mouse_cols_mc38[colnames(state_by_m_mc38_3)], pch=19, xaxt = "n", ylim=c(0,ylim_max), main=paste("mc38 3 p=", names(leaders_col)[prog], sep=""), xlab="Cond", ylab="% Active")
  axis(1, at=1:4, labels=c("ctrl", "pd1", "41bb", "41bb_pd1"))
  dev.off()
  
  png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_conds_m/", "scatter_active_",prog, "_","b16", "_time_", "3", ".png",sep=""), 8 * 100, 8 * 100, res = 200)
  plot(0.1*runif(ncol(state_by_m_b16_3)) + col2ind[mouse_cols_b16[colnames(state_by_m_b16_3)]], state_by_m_b16_3[prog,], col=mouse_cols_b16[colnames(state_by_m_b16_3)], pch=19, xaxt = "n", ylim=c(0,ylim_max), main=paste("b16 3 p=", names(leaders_col)[prog], sep=""), xlab="Cond", ylab="% Active")
  axis(1, at=1:4, labels=c("ctrl", "pd1", "41bb", "41bb_pd1"))
  dev.off()
}

for (prog in 1:ncol(norm_global_w)){
  png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_conds_m/", "barplot_active_",prog, "_","spec", "_time_", "3", ".png",sep=""), 12 * 100, 12 * 100, res = 200)
  barplot(state_by_m_spec_3[prog,], horiz =F, col = mouse_cols_spec[colnames(state_by_m_spec_3)], main=paste("spec 3 p=", names(leaders_col)[prog], sep=""), ylim=c(0,max(c(0.6, max(state_by_m_spec_3[prog,])))), las=2)
  dev.off()
  
  png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_conds_m/", "barplot_active_",prog, "_","mc38", "_time_", "3", ".png",sep=""), 12 * 100, 12 * 100, res = 200)
  barplot(state_by_m_mc38_3[prog,], horiz =F, col = mouse_cols_mc38[colnames(state_by_m_mc38_3)], main=paste("mc38 3 p=", names(leaders_col)[prog], sep=""), ylim=c(0,max(c(0.6, max(state_by_m_mc38_3[prog,])))), las=2)
  dev.off()
  
  png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_conds_m/", "barplot_active_",prog, "_","b16", "_time_", "3", ".png",sep=""), 12 * 100, 12 * 100, res = 200)
  barplot(state_by_m_b16_3[prog,], horiz =F, col = mouse_cols_b16[colnames(state_by_m_b16_3)], main=paste("b16 3 p=", names(leaders_col)[prog], sep=""), ylim=c(0,max(c(0.6, max(state_by_m_b16_3[prog,])))), las=2)
  dev.off()
}

mouse_cols_spec = rep("white", ncol(state_by_m_spec_6))
names(mouse_cols_spec) = colnames(state_by_m_spec_6)
mouse_cols_mc38 = rep("white", ncol(state_by_m_mc38_6))
names(mouse_cols_mc38) = colnames(state_by_m_mc38_6)
mouse_cols_b16 = rep("white", ncol(state_by_m_b16_6))
names(mouse_cols_b16) = colnames(state_by_m_b16_6)
for (cond in colnames(m_id_cond)){
  if (cond != "na"){
    m_ids = rownames(m_id_cond)[m_id_cond[,cond] > 0]
    mouse_cols_spec[colnames(state_by_m_spec_6) %in% m_ids] = cond2col[cond]

    mouse_cols_mc38[colnames(state_by_m_mc38_6) %in% m_ids] = cond2col[cond]

    mouse_cols_b16[colnames(state_by_m_b16_6) %in% m_ids] = cond2col[cond]

  }
}

# state_by_m_spec_3 = state_by_m_spec_3[,order(mouse_cols_spec)]
# state_by_m_spec_3_df = as.data.frame(t(state_by_m_spec_3))
# state_by_m_spec_3_df$cond =  col2cond[mouse_cols_spec[colnames(state_by_m_spec_3)]]
# state_by_m_mc38_3 = state_by_m_mc38_3[,order(mouse_cols_mc38)]
# state_by_m_mc38_3_df = as.data.frame(t(state_by_m_mc38_3))
# state_by_m_mc38_3_df$cond =  col2cond[mouse_cols_mc38[colnames(state_by_m_mc38_3)]]
# state_by_m_b13_3 = state_by_m_b16_3[,order(mouse_cols_b16)]
# state_by_m_b16_3_df = as.data.frame(t(state_by_m_b16_3))
# state_by_m_b16_3_df$cond =  col2cond[mouse_cols_b16[colnames(state_by_m_b16_3)]]

state_by_m_spec_6 = state_by_m_spec_6[,order(mouse_cols_spec)]
state_by_m_spec_6_df = as.data.frame(t(state_by_m_spec_6))
state_by_m_spec_6_df$cond =  col2cond[mouse_cols_spec[colnames(state_by_m_spec_6)]]
state_by_m_mc38_6 = state_by_m_mc38_6[,order(mouse_cols_mc38)]
state_by_m_mc38_6_df = as.data.frame(t(state_by_m_mc38_6))
state_by_m_mc38_6_df$cond =  col2cond[mouse_cols_mc38[colnames(state_by_m_mc38_6)]]
state_by_m_b16_6 = state_by_m_b16_6[,order(mouse_cols_b16)]
state_by_m_b16_6_df = as.data.frame(t(state_by_m_b16_6))
state_by_m_b16_6_df$cond =  col2cond[mouse_cols_b16[colnames(state_by_m_b16_6)]]

library(ggplot2)
library(ggpubr)

state_by_m_spec_6_df_pd1 = state_by_m_spec_6_df[as.character(state_by_m_spec_6_df$cond) %in% c("ctrl", "pd1"),]
state_by_m_mc38_6_df_pd1 = state_by_m_mc38_6_df[as.character(state_by_m_mc38_6_df$cond) %in% c("ctrl", "pd1"),]
state_by_m_b16_6_df_pd1 = state_by_m_b16_6_df[as.character(state_by_m_b16_6_df$cond) %in% c("ctrl", "pd1"),]
for (prog in 1:ncol(norm_global_w)){
  ylim_max = max(c(0.8, max(state_by_m_spec_3[prog,]), max(state_by_m_mc38_3[prog,]), max(state_by_m_b16_3[prog,]), max(state_by_m_spec_6[prog,]), max(state_by_m_mc38_6[prog,]), max(state_by_m_b16_6[prog,]))) + 0.1
  colord = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF")
  names(colord) = c("pd1", "ctrl", "41bb", "41bb_pd1")
  
  p_spec = ggviolin(state_by_m_spec_6_df, x="cond", y=colnames(state_by_m_spec_6_df)[prog], fill="cond", add="jitter", palette = "jco") + scale_fill_manual(values=colord, labels=names(colord)) + scale_x_discrete(limits=c("pd1", "ctrl", "41bb", "41bb_pd1")) 
  p_spec = p_spec + theme(legend.position="none", text = element_text(size=20)) + ylim(0, ylim_max)
  ggsave(filename = paste(violin_dir, "time_6_violin_spec_41bb_pd1_",prog, ".png",sep=""), plot=p_spec, width = 4, height = 1, scale = 3) #+ stat_compare_means(comparisons=list(c("ctrl", "pd1"), c("ctrl", "41bb"), c("ctrl", "41bb_pd1"),c("pd1", "41bb_pd1")),method = "wilcox.test", label = "p.signif") 
  
  p_mc38 = ggviolin(state_by_m_mc38_6_df, x="cond", y=colnames(state_by_m_mc38_6_df)[prog], fill="cond", add="jitter", palette = "jco") + scale_fill_manual(values=colord, labels=names(colord)) + scale_x_discrete(limits=c("pd1", "ctrl", "41bb", "41bb_pd1")) 
  p_mc38 = p_mc38 + theme(legend.position="none", text = element_text(size=20)) + ylim(0, ylim_max) #+ stat_compare_means(comparisons=list(c("ctrl", "pd1"), c("ctrl", "41bb"), c("ctrl", "41bb_pd1"),c("pd1", "41bb_pd1")),method = "wilcox.test", label = "p.signif") 
  ggsave(filename = paste(violin_dir, "time_6_violin_mc38_41bb_pd1_",prog, ".png",sep=""), plot=p_mc38, width = 4, height = 1, scale = 3)
  
  p_b16 = ggviolin(state_by_m_b16_6_df, x="cond", y=colnames(state_by_m_b16_6_df)[prog], fill="cond", add="jitter", palette = "jco") + scale_fill_manual(values=colord, labels=names(colord)) + scale_x_discrete(limits=c("pd1", "ctrl", "41bb", "41bb_pd1")) 
  p_b16 = p_b16 + theme(legend.position="none", text = element_text(size=20)) + ylim(0, ylim_max) # + stat_compare_means(comparisons=list(c("ctrl", "pd1"), c("ctrl", "41bb"), c("ctrl", "41bb_pd1"),c("pd1", "41bb_pd1")),method = "wilcox.test", label = "p.signif") 
  ggsave(filename = paste(violin_dir, "time_6_violin_b16_41bb_pd1_",prog, ".png",sep=""), plot=p_b16, width = 4, height = 1, scale = 3)
}

state_by_m_spec_6_df_pd1$time = "6"
state_by_m_spec_6_df_pd1$cond_time = paste(state_by_m_spec_6_df_pd1$cond, state_by_m_spec_6_df_pd1$time, sep="_")
state_by_m_spec_3_df_pd1$time = "3"
state_by_m_spec_3_df_pd1$cond_time = paste(state_by_m_spec_3_df_pd1$cond, state_by_m_spec_3_df_pd1$time, sep="_")
state_by_m_spec_6_df_pd1_save = state_by_m_spec_6_df_pd1
state_by_m_spec_6_df_pd1 = rbind(state_by_m_spec_3_df_pd1, state_by_m_spec_6_df_pd1)
ylim_spec_6 = c(0.7, 0.7, 0.5, 0.4, 0.3, 0.4, 0.6, 0.4, 0.6, 0.3, 0.4)
for (prog in 1:ncol(norm_global_w)){
  ylim_max = max(c(0.8, max(state_by_m_spec_3[prog,]), max(state_by_m_mc38_3[prog,]), max(state_by_m_b16_3[prog,]), max(state_by_m_spec_6[prog,]), max(state_by_m_mc38_6[prog,]), max(state_by_m_b16_6[prog,]))) #+ 0.1
  colord = c("#EFC000FF", "#EFC000FF", "#0073C2FF", "#0073C2FF")
  names(colord) = c("ctrl_3", "ctrl_6", "pd1_3", "pd1_6")
  
  p_spec = ggviolin(state_by_m_spec_6_df_pd1, x="cond_time", y=colnames(state_by_m_spec_6_df_pd1)[prog], fill="cond_time", add="jitter", palette = "jco")
  p_spec = p_spec + theme(legend.position="none", text = element_text(size=20)) + scale_x_discrete(limits=c("ctrl_3", "ctrl_6", "pd1_3", "pd1_6")) + ylim(0, ylim_spec_6[prog]) + scale_fill_manual(values=colord, labels=names(colord))# + stat_compare_means(comparisons=list(c("ctrl", "pd1")),method = "wilcox.test", label = "p.signif", bracket.size = 3) 
  ggsave(filename = paste(violin_dir, "time_6_and_3_violin_spec_pd1_",prog, ".png",sep=""), plot=p_spec) 
  
  p_mc38 = ggviolin(state_by_m_mc38_6_df, x="cond", y=colnames(state_by_m_mc38_6_df)[prog], fill="cond", add="jitter", palette = "jco")
  p_mc38 = p_mc38 + theme(legend.position="none", text = element_text(size=30)) + ylim(0, ylim_max) + scale_x_discrete(limits=c("ctrl", "pd1")) + stat_compare_means(comparisons=list(c("ctrl", "pd1")),method = "wilcox.test", label = "p.signif", bracket.size = 3) 
  ggsave(filename = paste(violin_dir, "time_6_violin_mc38_pd1_",prog, ".png",sep=""), plot=p_mc38)
  
  p_b16 = ggviolin(state_by_m_b16_6_df, x="cond", y=colnames(state_by_m_b16_6_df)[prog], fill="cond", add="jitter", palette = "jco")
  p_b16 = p_b16 + theme(legend.position="none", text = element_text(size=30)) + ylim(0, ylim_max) + scale_x_discrete(limits=c("ctrl", "pd1")) + stat_compare_means(comparisons=list(c("ctrl", "pd1")),method = "wilcox.test", label = "p.signif", bracket.size = 3) 
  ggsave(filename=paste(violin_dir, "time_6_violin_b16_pd1_", prog, ".png", sep=""), plot=p_b16)
}


#### CC analysis
non_filt_gset = scdb_gset(mm_filt_id)
frac_mat = t(t(mat@mat[,names(mc@mc)]) / colSums(mat@mat[,names(mc@mc)]))
cc_genes = names(non_filt_gset@gene_set)[non_filt_gset@gene_set %in% c(110,11,60,130)]
write.csv(cc_genes, "")

# cc_genes = c(cc_g1_genes, cc_g2_genes)
# cc_genes = cc_g2_genes
cc_sig = apply(mc@e_gc[cc_genes,],2, mean)
cc_sig_by_cell = cc_sig[mc@mc]

cc_sig_by_cell = apply(frac_mat[cc_genes,names(mc@mc)], 2, sum)
names(cc_sig_by_cell) = names(mc@mc)

all_gene_means = rowMeans(mat@mat)
s_genes = cc_g1_genes#unique(c(names(non_filt_gset@gene_set)[non_filt_gset@gene_set %in% c(110)], cc_g1_genes))
s_genes = setdiff(s_genes, c("Dek", "Dut","Rrm2"))
m_genes = cc_g2_genes#unique(c(names(non_filt_gset@gene_set)[non_filt_gset@gene_set %in% c(11,60,130)], cc_g2_genes))
m_genes = setdiff(m_genes, c("Tuba1b", "Hmgb2", "Tubb5", "Tubb4b", "Kif2c", "Aurka", "Cdca2", "Ankrd11"))
cc1_sig_by_cell = apply(frac_mat[s_genes,names(mc@mc)], 2, sum)
names(cc1_sig_by_cell) = names(mc@mc)
cc2_sig_by_cell = apply(frac_mat[m_genes,names(mc@mc)], 2, sum)
names(cc2_sig_by_cell) = names(mc@mc)

cc1_sig = apply(mc@e_gc[s_genes,], 2, sum)
cc2_sig = apply(mc@e_gc[m_genes,], 2, sum)

hits_egc_cc = mc@e_gc[,cc1_sig > 0.001 | cc_sig > 0.0025]
hits_egc_cc1_cor = as.numeric(cor(t(hits_egc_cc), cc1_sig[cc1_sig > 0.001 | cc_sig > 0.0025]))
hits_egc_cc2_cor = as.numeric(cor(t(hits_egc_cc), cc2_sig[cc1_sig > 0.001 | cc_sig > 0.0025]))
names(hits_egc_cc1_cor) = rownames(hits_egc_cc)
names(hits_egc_cc2_cor) = rownames(hits_egc_cc)
add_s_genes = names(tail(sort(hits_egc_cc1_cor), 50))[!names(tail(sort(hits_egc_cc1_cor), 50)) %in% m_genes]
add_m_genes = names(tail(sort(hits_egc_cc2_cor), 50))[!names(tail(sort(hits_egc_cc2_cor), 50)) %in% s_genes]

add_s_genes_f = setdiff(add_s_genes, add_m_genes)
add_s_genes_f = add_s_genes_f[hits_egc_cc2_cor[add_s_genes_f] < 0.4]
add_m_genes_f = setdiff(add_m_genes, add_s_genes)
add_m_genes_f = add_m_genes_f[hits_egc_cc1_cor[add_m_genes_f] < 0.4]

s_genes = unique(c(s_genes, add_s_genes_f))
m_genes = unique(c(m_genes, add_m_genes_f))

cc1_sig_by_cell = apply(frac_mat[s_genes,names(mc@mc)], 2, sum)
names(cc1_sig_by_cell) = names(mc@mc)
cc2_sig_by_cell = apply(frac_mat[m_genes,names(mc@mc)], 2, sum)
names(cc2_sig_by_cell) = names(mc@mc)

cc1_sig = apply(mc@e_gc[s_genes,], 2, sum)
cc2_sig = apply(mc@e_gc[m_genes,], 2, sum)

plot(density(cc1_sig_by_cell))
plot(density(cc2_sig_by_cell))


cc_sig = apply(mc@e_gc[cc_genes,], 2, sum)
cc_sig_by_cell = cc_sig[mc@mc]
names(cc_sig_by_cell) = names(mc@mc)

par(mar=c(3,3,3,3))
plot(density(cc_sig_by_cell))
abline(v=0.01, col="red")

all_cc_genes = cc_genes#c(s_genes, m_genes)
cc_sig_by_cell = apply(frac_mat[all_cc_genes,names(mc@mc)], 2, sum)
names(cc_sig_by_cell) = names(mc@mc)

cc_pos_cells = names(cc_sig_by_cell)[cc_sig_by_cell > 0.01]

s_pos_cells = names(cc1_sig_by_cell)[cc1_sig_by_cell > 0.005]
m_pos_cells = names(cc2_sig_by_cell)[cc2_sig_by_cell > 0.005]

cc_cells = unique(c(s_pos_cells, m_pos_cells))
cc_cells_cols = mc@colors[mc@mc[cc_cells]]
cc_gene_cor = tgs_cor(as.matrix(t(frac_mat[all_cc_genes,cc_cells])))
diag(cc_gene_cor) = 0
pheatmap(cc_gene_cor)

cc_gene_cor_mc = tgs_cor(as.matrix(t(mc@e_gc[all_cc_genes,mc@mc[cc_cells]])))
diag(cc_gene_cor_mc) = 0
pheatmap(cc_gene_cor_mc)

plot(log2(cc2_sig_by_cell[cc_cells]), log2(cc1_sig_by_cell[cc_cells]), col=cc_cells_cols, pch=19)

cc_gene_cor_all = tgs_cor(as.matrix(t(frac_mat[all_cc_genes,names(mc@mc)])))
diag(cc_gene_cor_all) = 0
pheatmap(cc_gene_cor_all)

n_cells_mc = table(mc@mc)
n_cells_s = rep(0, length(n_cells_mc))
names(n_cells_s) = names(n_cells_mc)
n_cells_m = rep(0, length(n_cells_mc))
names(n_cells_m) = names(n_cells_mc)
n_cells_cc = rep(0, length(n_cells_mc))
names(n_cells_cc) = names(n_cells_mc)

n_cells_s[names(table(mc@mc[s_pos_cells]))] = table(mc@mc[s_pos_cells])
n_cells_m[names(table(mc@mc[m_pos_cells]))] = table(mc@mc[m_pos_cells])
n_cells_cc[names(table(mc@mc[cc_pos_cells]))] = table(mc@mc[cc_pos_cells])
p_cells_s = as.numeric(n_cells_s / n_cells_mc)
p_cells_m = as.numeric(n_cells_m / n_cells_mc)
p_cells_cc = as.numeric(n_cells_cc / n_cells_mc)
plot(p_cells_m, p_cells_s, col=mc@colors, pch=19)

for (prog in 1:ncol(norm_global_mat)){
  active_thres = active_thres_vec[prog]
  active_mcs = (1:max(mc@mc))[norm_global_mat[,prog] >= active_thres]
  png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_conds_m/", "cc_vs_prog_",prog, "_","active_mcs", ".png",sep=""), 8 * 100, 8 * 100, res = 150)
  plot(norm_global_mat[active_mcs,prog], p_cells_s[active_mcs], col=mc@colors[active_mcs], pch=19, main=colnames(reg_full_w)[prog])
  print(cor(norm_global_mat[active_mcs,prog], p_cells_s[active_mcs]))
  dev.off()
}

s_sig = apply(mc@e_gc[s_genes,],2, mean)
m_sig = apply(mc@e_gc[m_genes,],2, mean)
plot(m_sig, s_sig, col=mc@colors, pch=19)


s_sig = apply(mc@e_gc[s_genes,],2, sum)
m_sig = apply(mc@e_gc[m_genes,],2, sum)
plot(m_sig, s_sig, col=mc@colors, pch=19)

cc_by_p_spec_all = c()

cc_by_p_spec_ctrl = c()
cc_by_p_spec_ctrl_cells = c()
cc_by_p_spec_pd1 = c()
cc_by_p_spec_pd1_cells = c()
cc_min_cells = 200
for (prog in 1:ncol(norm_global_mat)){
  active_thres = active_thres_vec[prog]
  strong_thres = strong_thres_vec[prog]
  active_mcs = (1:max(mc@mc))[norm_global_mat[,prog] >= active_thres]
  active_cells = rel_cells[mc@mc[rel_cells] %in% active_mcs]
  spec_cells_ctrl = rel_cells[is_spec & is_tumor & is_ctrl & md$days_post_transfer == 6]
  spec_cells_pd1 = rel_cells[is_spec & is_tumor & is_pd1 & md$days_post_transfer == 6]
  spec_cells_ctrl_active = spec_cells_ctrl[spec_cells_ctrl %in% active_cells]
  cc_by_p_spec_ctrl_cells = c(cc_by_p_spec_ctrl_cells, length(spec_cells_ctrl_active))
  spec_cells_pd1_active = spec_cells_pd1[spec_cells_pd1 %in% active_cells]
  cc_by_p_spec_pd1_cells = c(cc_by_p_spec_pd1_cells, length(spec_cells_pd1_active))
  if (length(spec_cells_ctrl_active) > cc_min_cells){
    cc_by_p_spec_ctrl = c(cc_by_p_spec_ctrl, sum(spec_cells_ctrl_active %in% m_pos_cells) / length(spec_cells_ctrl_active))
  }
  else{
    cc_by_p_spec_ctrl = c(cc_by_p_spec_ctrl,0)
  }
  if (length(spec_cells_pd1_active) > cc_min_cells){
    cc_by_p_spec_pd1 = c(cc_by_p_spec_pd1, sum(spec_cells_pd1_active %in% m_pos_cells) / length(spec_cells_pd1_active))
  }
  else{
    cc_by_p_spec_pd1 = c(cc_by_p_spec_pd1,0)
  }
}
names(cc_by_p_spec_ctrl) = names(leaders_col)
cc_by_p_spec_ctrl_tum = cc_by_p_spec_ctrl
names(cc_by_p_spec_pd1) = names(leaders_col)
cc_by_p_spec_pd1_tum = cc_by_p_spec_pd1
cc_by_p_spec = rbind(cc_by_p_spec_ctrl, cc_by_p_spec_pd1)

cc_by_p_spec_all = rbind(cc_by_p_spec_all, cc_by_p_spec)
barplot(cc_by_p_spec, beside = T, horiz = F, las=2, col = c("#EFC000FF", "#0073C2FF"), cex.names = 0.7, cex.axis = 0.9,ylab="CC_mean", main="Tumor day 6 CC_all Ctrl vs PD1")
plot(cc_by_p_spec_ctrl, cc_by_p_spec_pd1, bg="black", pch=19) #, xlim=c(0.003, 0.015), ylim=c(0.003, 0.015)
abline(0,1,col="red")
text(cc_by_p_spec_ctrl, cc_by_p_spec_pd1, labels = names(cc_by_p_spec_ctrl), pos = 3)

cc_by_p_spec_ctrl = c()
cc_by_p_spec_ctrl_cells = c()
cc_by_p_spec_pd1 = c()
cc_by_p_spec_pd1_cells = c()
cc_min_cells = 100
for (prog in 1:ncol(norm_global_mat)){
  active_thres = active_thres_vec[prog]
  strong_thres = strong_thres_vec[prog]
  active_mcs = (1:max(mc@mc))[norm_global_mat[,prog] >= active_thres]
  active_cells = rel_cells[mc@mc[rel_cells] %in% active_mcs]
  spec_cells_ctrl = rel_cells[is_spec & (is_ln | is_spleen) & is_ctrl & md$days_post_transfer == 3]
  spec_cells_pd1 = rel_cells[is_spec & (is_ln | is_spleen) & is_pd1 & md$days_post_transfer == 3]
  spec_cells_ctrl_active = spec_cells_ctrl[spec_cells_ctrl %in% active_cells]
  cc_by_p_spec_ctrl_cells = c(cc_by_p_spec_ctrl_cells, length(spec_cells_ctrl_active))
  spec_cells_pd1_active = spec_cells_pd1[spec_cells_pd1 %in% active_cells]
  cc_by_p_spec_pd1_cells = c(cc_by_p_spec_pd1_cells, length(spec_cells_pd1_active))
  if (length(spec_cells_ctrl_active) > cc_min_cells){
    cc_by_p_spec_ctrl = c(cc_by_p_spec_ctrl, sum(spec_cells_pd1_active %in% m_pos_cells) / length(spec_cells_ctrl_active))
  }
  else{
    cc_by_p_spec_ctrl = c(cc_by_p_spec_ctrl,0)
  }
  if (length(spec_cells_pd1_active) > cc_min_cells){
    cc_by_p_spec_pd1 = c(cc_by_p_spec_pd1, sum(spec_cells_pd1_active %in% m_pos_cells) / length(spec_cells_pd1_active))
  }
  else{
    cc_by_p_spec_pd1 = c(cc_by_p_spec_pd1,0)
  }
}
names(cc_by_p_spec_ctrl) = names(leaders_col)
names(cc_by_p_spec_pd1) = names(leaders_col)
cc_by_p_spec_ctrl_peri = cc_by_p_spec_ctrl
cc_by_p_spec = rbind(cc_by_p_spec_ctrl, cc_by_p_spec_pd1)

cc_by_p_spec_all = rbind(cc_by_p_spec_all, cc_by_p_spec)
barplot(cc_by_p_spec_all, beside = T, horiz = F, las=2, cex.names = 0.7, cex.axis = 0.9,ylab="CC_mean", main="Peri day 3 CC Ctrl vs PD1") # col = c("#EFC000FF", "#0073C2FF"),


plot(cc_by_p_spec_ctrl, cc_by_p_spec_pd1, bg="black", pch=19, xlim=c(0.00, 0.018), ylim=c(0.00, 0.018))
abline(0,1,col="red")
text(cc_by_p_spec_ctrl, cc_by_p_spec_pd1, labels = names(cc_by_p_spec_ctrl), pos = 4)

plot(cc_by_p_spec_ctrl_peri, cc_by_p_spec_ctrl_tum, bg="black", pch=19)
abline(0,1,col="red")
text(cc_by_p_spec_ctrl_peri, cc_by_p_spec_ctrl_tum, labels = names(cc_by_p_spec_ctrl_peri), pos = 4)


naive_mcs = rownames(reg_full_w)[norm_global_mat[,2] >= active_thres_vec[2] | norm_global_mat[,9] >= active_thres_vec[9] | norm_global_mat[,11] >= active_thres_vec[11] | norm_global_mat[,1] >= active_thres_vec[1]]
# naive_mcs = rownames(reg_full_w)[norm_global_mat[,2] >= active_thres_vec[2] | norm_global_mat[,9] >= active_thres_vec[9] | norm_global_mat[,11] >= active_thres_vec[11]]
naive_score = norm_global_mat[,9] + norm_global_mat[,2] + norm_global_mat[,11] - norm_global_mat[,1]
# naive_score = norm_global_mat[,9] + norm_global_mat[,2] - norm_global_mat[,11] 
naive_score = naive_score / max(naive_score)
names(naive_score) = rownames(reg_full_w)

naive_mcs_ord = naive_mcs[order(naive_score[naive_mcs], decreasing = T)]
naive_mcs_div = seq(1, length(naive_mcs_ord), by=floor(length(naive_mcs)/10))
naive_mcs_div[length(naive_mcs_div)] = length(naive_mcs_ord)
naive_mcs_ord_ls = list()
for (ii in 1:(length(naive_mcs_div) - 1)){
  naive_mcs_ord_ls[[ii]] = naive_mcs_ord[naive_mcs_div[ii]:naive_mcs_div[ii+1]]
}


min_cells = 20
day_naive = c(1:3)

time_genes = c("Tcf7", "Sell", "Tnfrsf9", "Prf1", "Pdcd1", "Cxcr6", "Irf8", "Zbtb2", "Id2", "Id3", "Tigit", "Myc", "Lag3", "Havcr2", "Prdm1", "Il2ra", "Il2rb", "Ccr7")
time_cc_genes = c(time_genes, all_cc_genes)
# downs_tum_spec = downs_mat[,md[colnames(downs_mat),"location"] == "tumor" & md[colnames(downs_mat),"cell_type"] %in% c("spec", "cd8_in_vitro") & as.character(md[colnames(downs_mat),"condition"]) == "ctrl"]
# downs_tum_spec_genes = downs_tum_spec[time_cc_genes,]

downs_mat_all = scm_downsamp(mat@mat[,names(mc@mc)], 750)
downs_mat_all = cbind(downs_mat_all, mat@mat[,setdiff(names(mc@mc), colnames(downs_mat_all))])
downs_mat_all = downs_mat_all[,names(mc@mc)]

downs_mat_time_genes = downs_mat_all[time_cc_genes,]
for (gene in all_cc_genes){
  naive_by_p_spec_ctrl_tum = c()
  naive_by_p_spec_ctrl_ln = c()
  naive_by_p_spec_ctrl_spleen = c()
  
  naive_by_p_spec_pd1_tum = c()
  naive_by_p_spec_pd1_ln = c()
  naive_by_p_spec_pd1_spleen = c()
  for (quant in 1:length(naive_mcs_ord_ls)){
    active_mcs = naive_mcs_ord_ls[[quant]]
    active_cells = rel_cells[mc@mc[rel_cells] %in% active_mcs]
    
    spec_cells_ctrl = rel_cells[is_spec & (is_tumor) & is_ctrl & md$days_post_transfer %in% day_naive]
    spec_cells_pd1 = rel_cells[is_spec & (is_tumor) & is_pd1 & md$days_post_transfer %in% day_naive]
    spec_cells_ctrl_active = spec_cells_ctrl[spec_cells_ctrl %in% active_cells]
    if (min_cells < length(spec_cells_ctrl_active)){
      naive_by_p_spec_ctrl_tum = c(naive_by_p_spec_ctrl_tum, mean(downs_mat_time_genes[gene,spec_cells_ctrl_active]))
    }
    else{
      naive_by_p_spec_ctrl_tum = c(naive_by_p_spec_ctrl_tum, 0)
    }
    
    spec_cells_pd1_active = spec_cells_pd1[spec_cells_pd1 %in% active_cells]
    if (min_cells < length(spec_cells_pd1_active)){
      naive_by_p_spec_pd1_tum = c(naive_by_p_spec_pd1_tum, mean(downs_mat_time_genes[gene,spec_cells_pd1_active]))
      
    }
    else{
      naive_by_p_spec_pd1_tum = c(naive_by_p_spec_pd1_tum, 0)
    }
    
    
    spec_cells_ctrl = rel_cells[is_spec & (is_ln) & is_ctrl & md$days_post_transfer %in% day_naive]
    spec_cells_pd1 = rel_cells[is_spec & (is_ln) & is_pd1 & md$days_post_transfer %in% day_naive]
    spec_cells_ctrl_active = spec_cells_ctrl[spec_cells_ctrl %in% active_cells]
    if (min_cells < length(spec_cells_ctrl_active)){
      naive_by_p_spec_ctrl_ln = c(naive_by_p_spec_ctrl_ln, mean(downs_mat_time_genes[gene,spec_cells_ctrl_active]))
      
    }
    else{
      naive_by_p_spec_ctrl_ln = c(naive_by_p_spec_ctrl_ln, 0)
      
    }
    spec_cells_pd1_active = spec_cells_pd1[spec_cells_pd1 %in% active_cells]
    if (min_cells < length(spec_cells_pd1_active)){
      naive_by_p_spec_pd1_ln = c(naive_by_p_spec_pd1_ln, mean(downs_mat_time_genes[gene,spec_cells_pd1_active]))
      
    }
    else{
      naive_by_p_spec_pd1_ln = c(naive_by_p_spec_pd1_ln, 0)
      
    }
    
    
    spec_cells_ctrl = rel_cells[is_spec & (is_spleen) & is_ctrl & md$days_post_transfer %in% day_naive]
    spec_cells_pd1 = rel_cells[is_spec & (is_spleen) & is_pd1 & md$days_post_transfer %in% day_naive]
    spec_cells_ctrl_active = spec_cells_ctrl[spec_cells_ctrl %in% active_cells]
    if (min_cells < length(spec_cells_ctrl_active)){
      naive_by_p_spec_ctrl_spleen = c(naive_by_p_spec_ctrl_spleen, mean(downs_mat_time_genes[gene,spec_cells_ctrl_active]))
      
    }
    else{
      naive_by_p_spec_ctrl_spleen = c(naive_by_p_spec_ctrl_spleen, 0)
      c
    }
    
    spec_cells_pd1_active = spec_cells_pd1[spec_cells_pd1 %in% active_cells]
    if (min_cells < length(spec_cells_ctrl_active)){
      naive_by_p_spec_pd1_spleen = c(naive_by_p_spec_pd1_spleen, mean(downs_mat_time_genes[gene,spec_cells_pd1_active]))
    }
    else{
      naive_by_p_spec_pd1_spleen = c(naive_by_p_spec_pd1_spleen, 0)
    }
    
  }
  png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_conds_m/", "naive_traj_active_cells_early_ctrl_pd1_spleen_cc_", gene, ".png",sep=""), 8 * 100, 8 * 100, res = res_cc)
  plot(naive_by_p_spec_ctrl_spleen, pch=19, col=leaders_col[9], ylim=c(0, max(c(naive_by_p_spec_ctrl_spleen, naive_by_p_spec_pd1_spleen))), xlab="naive_quant", ylab="% in quantile", main=gene)
  lines(naive_by_p_spec_ctrl_spleen, col=leaders_col[9])
  
  points(naive_by_p_spec_pd1_spleen, pch=21, col=leaders_col[9])
  lines(naive_by_p_spec_pd1_spleen, col=leaders_col[9], lty=3)
  dev.off()
  
  png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_conds_m/", "naive_traj_active_cells_early_ctrl_pd1_tum_cc_", gene, ".png",sep=""), 8 * 100, 8 * 100, res = res_cc)
  plot(naive_by_p_spec_ctrl_tum, pch=19, col=leaders_col[11], ylim=c(0, max(c(naive_by_p_spec_ctrl_tum, naive_by_p_spec_pd1_tum))), xlab="naive_quant", ylab="% in quantile", main=gene)
  lines(naive_by_p_spec_ctrl_tum, col=leaders_col[11])
  
  points(naive_by_p_spec_pd1_tum, pch=21, col=leaders_col[11])
  lines(naive_by_p_spec_pd1_tum, col=leaders_col[11], lty=3)
  dev.off()
  
  png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_conds_m/", "naive_traj_active_cells_early_ctrl_pd1_ln_cc_", gene, ".png",sep=""), 8 * 100, 8 * 100, res = res_cc)
  plot(naive_by_p_spec_ctrl_ln, pch=19, col=leaders_col[2], ylim=c(0, max(c(naive_by_p_spec_ctrl_ln, naive_by_p_spec_pd1_ln))), xlab="naive_quant", ylab="% in quantile", main=gene)
  lines(naive_by_p_spec_ctrl_ln, col=leaders_col[2])
  
  points(naive_by_p_spec_pd1_ln, pch=21, col=leaders_col[2])
  lines(naive_by_p_spec_pd1_ln, col=leaders_col[2], lty=3)
  dev.off()
}



naive_by_p_spec_ctrl_tum = c()
naive_by_p_spec_ctrl_ln = c()
naive_by_p_spec_ctrl_spleen = c()
cc_naive_by_p_spec_ctrl_tum = c()
cc_naive_by_p_spec_ctrl_ln = c()
cc_naive_by_p_spec_ctrl_spleen = c()

naive_by_p_spec_pd1_tum = c()
naive_by_p_spec_pd1_ln = c()
naive_by_p_spec_pd1_spleen = c()
cc_naive_by_p_spec_pd1_tum = c()
cc_naive_by_p_spec_pd1_ln = c()
cc_naive_by_p_spec_pd1_spleen = c()
for (quant in 1:length(naive_mcs_ord_ls)){
  active_mcs = naive_mcs_ord_ls[[quant]]
  active_cells = rel_cells[mc@mc[rel_cells] %in% active_mcs]
  
  spec_cells_ctrl = rel_cells[is_spec & (is_tumor) & is_ctrl & md$days_post_transfer %in% day_naive]
  spec_cells_pd1 = rel_cells[is_spec & (is_tumor) & is_pd1 & md$days_post_transfer %in% day_naive]
  spec_cells_ctrl_active = spec_cells_ctrl[spec_cells_ctrl %in% active_cells]
  if (min_cells < length(spec_cells_ctrl_active)){
    naive_by_p_spec_ctrl_tum = c(naive_by_p_spec_ctrl_tum, length(spec_cells_ctrl_active)/length(spec_cells_ctrl))
    cc_naive_by_p_spec_ctrl_tum = c(cc_naive_by_p_spec_ctrl_tum, sum(spec_cells_ctrl_active %in% cc_pos_cells)/length(spec_cells_ctrl_active))
  }
  else{
    naive_by_p_spec_ctrl_tum = c(naive_by_p_spec_ctrl_tum, 0)
    cc_naive_by_p_spec_ctrl_tum = c(cc_naive_by_p_spec_ctrl_tum, 0)
  }
  
  spec_cells_pd1_active = spec_cells_pd1[spec_cells_pd1 %in% active_cells]
  if (min_cells < length(spec_cells_pd1_active)){
    naive_by_p_spec_pd1_tum = c(naive_by_p_spec_pd1_tum, length(spec_cells_pd1_active)/length(spec_cells_pd1))
    cc_naive_by_p_spec_pd1_tum = c(cc_naive_by_p_spec_pd1_tum, sum(spec_cells_pd1_active %in% cc_pos_cells)/length(spec_cells_pd1_active))
  }
  else{
    naive_by_p_spec_pd1_tum = c(naive_by_p_spec_pd1_tum, 0)
    cc_naive_by_p_spec_pd1_tum = c(cc_naive_by_p_spec_pd1_tum, 0)
  }
  
  
  spec_cells_ctrl = rel_cells[is_spec & (is_ln) & is_ctrl & md$days_post_transfer %in% day_naive]
  spec_cells_pd1 = rel_cells[is_spec & (is_ln) & is_pd1 & md$days_post_transfer %in% day_naive]
  spec_cells_ctrl_active = spec_cells_ctrl[spec_cells_ctrl %in% active_cells]
  if (min_cells < length(spec_cells_ctrl_active)){
    naive_by_p_spec_ctrl_ln = c(naive_by_p_spec_ctrl_ln, length(spec_cells_ctrl_active)/length(spec_cells_ctrl))
    cc_naive_by_p_spec_ctrl_ln = c(cc_naive_by_p_spec_ctrl_ln, sum(spec_cells_ctrl_active %in% cc_pos_cells)/length(spec_cells_ctrl_active))
  }
  else{
    naive_by_p_spec_ctrl_ln = c(naive_by_p_spec_ctrl_ln, 0)
    cc_naive_by_p_spec_ctrl_ln = c(cc_naive_by_p_spec_ctrl_ln, 0)
  }
  spec_cells_pd1_active = spec_cells_pd1[spec_cells_pd1 %in% active_cells]
  if (min_cells < length(spec_cells_pd1_active)){
    naive_by_p_spec_pd1_ln = c(naive_by_p_spec_pd1_ln, length(spec_cells_pd1_active)/length(spec_cells_pd1))
    cc_naive_by_p_spec_pd1_ln = c(cc_naive_by_p_spec_pd1_ln, sum(spec_cells_pd1_active %in% cc_pos_cells)/length(spec_cells_pd1_active))
  }
  else{
    naive_by_p_spec_pd1_ln = c(naive_by_p_spec_pd1_ln, 0)
    cc_naive_by_p_spec_pd1_ln = c(cc_naive_by_p_spec_pd1_ln, 0)
  }
  
  
  spec_cells_ctrl = rel_cells[is_spec & (is_spleen) & is_ctrl & md$days_post_transfer %in% day_naive]
  spec_cells_pd1 = rel_cells[is_spec & (is_spleen) & is_pd1 & md$days_post_transfer %in% day_naive]
  spec_cells_ctrl_active = spec_cells_ctrl[spec_cells_ctrl %in% active_cells]
  if (min_cells < length(spec_cells_ctrl_active)){
    naive_by_p_spec_ctrl_spleen = c(naive_by_p_spec_ctrl_spleen, length(spec_cells_ctrl_active)/length(spec_cells_ctrl))
    cc_naive_by_p_spec_ctrl_spleen = c(cc_naive_by_p_spec_ctrl_spleen, sum(spec_cells_ctrl_active %in% cc_pos_cells)/length(spec_cells_ctrl_active))
  }
  else{
    naive_by_p_spec_ctrl_spleen = c(naive_by_p_spec_ctrl_spleen, 0)
    cc_naive_by_p_spec_ctrl_spleen = c(cc_naive_by_p_spec_ctrl_spleen, 0)
  }
  
  spec_cells_pd1_active = spec_cells_pd1[spec_cells_pd1 %in% active_cells]
  if (min_cells < length(spec_cells_ctrl_active)){
    naive_by_p_spec_pd1_spleen = c(naive_by_p_spec_pd1_spleen, length(spec_cells_pd1_active)/length(spec_cells_pd1))
    cc_naive_by_p_spec_pd1_spleen = c(cc_naive_by_p_spec_pd1_spleen, sum(spec_cells_pd1_active %in% cc_pos_cells)/length(spec_cells_pd1_active))
  }
  else{
    naive_by_p_spec_pd1_spleen = c(naive_by_p_spec_pd1_spleen, 0)
    cc_naive_by_p_spec_pd1_spleen = c(cc_naive_by_p_spec_pd1_spleen, 0)
  }
  
}

res_cc = 250
png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_conds_m/", "naive_traj_active_cells_early_ctrl", ".png",sep=""), 8 * 100, 8 * 100, res = res_cc)
plot(naive_by_p_spec_ctrl_spleen, pch=19, col=leaders_col[9], ylim=c(0,0.25), xlab="naive_quant", ylab="% in quantile", main="early_ctrl")
lines(naive_by_p_spec_ctrl_spleen, col=leaders_col[9])
points(naive_by_p_spec_ctrl_ln, pch=19, col=leaders_col[2])
lines(naive_by_p_spec_ctrl_ln, col=leaders_col[2])
points(naive_by_p_spec_ctrl_tum, pch=19, col=leaders_col[11])
lines(naive_by_p_spec_ctrl_tum, col=leaders_col[11])
dev.off()

png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_conds_m/", "naive_traj_active_cells_early_pd1", ".png",sep=""), 8 * 100, 8 * 100, res = res_cc)
plot(naive_by_p_spec_pd1_spleen, pch=19, col=leaders_col[9], ylim=c(0,0.25), xlab="naive_quant", ylab="% in quantile", main="early_pd1")
lines(naive_by_p_spec_pd1_spleen, col=leaders_col[9])
points(naive_by_p_spec_pd1_ln, pch=19, col=leaders_col[2])
lines(naive_by_p_spec_pd1_ln, col=leaders_col[2])
points(naive_by_p_spec_pd1_tum, pch=19, col=leaders_col[11])
lines(naive_by_p_spec_pd1_tum, col=leaders_col[11])
dev.off()


png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_conds_m/", "naive_traj_active_cells_early_ctrl_pd1_spleen", ".png",sep=""), 8 * 100, 8 * 100, res = res_cc)
plot(naive_by_p_spec_ctrl_spleen, pch=19, col=leaders_col[9], ylim=c(0,0.25), xlab="naive_quant", ylab="% in quantile", main="early_ctrl_pd1")
lines(naive_by_p_spec_ctrl_spleen, col=leaders_col[9])

points(naive_by_p_spec_pd1_spleen, pch=21, col=leaders_col[9])
lines(naive_by_p_spec_pd1_spleen, col=leaders_col[9], lty=3)
dev.off()

png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_conds_m/", "naive_traj_active_cells_early_ctrl_pd1_tum", ".png",sep=""), 8 * 100, 8 * 100, res = res_cc)
plot(naive_by_p_spec_ctrl_tum, pch=19, col=leaders_col[11], ylim=c(0,0.25), xlab="naive_quant", ylab="% in quantile", main="early_ctrl_pd1")
lines(naive_by_p_spec_ctrl_tum, col=leaders_col[11])

points(naive_by_p_spec_pd1_tum, pch=21, col=leaders_col[11])
lines(naive_by_p_spec_pd1_tum, col=leaders_col[11], lty=3)
dev.off()

png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_conds_m/", "naive_traj_active_cells_early_ctrl_pd1_ln", ".png",sep=""), 8 * 100, 8 * 100, res = res_cc)
plot(naive_by_p_spec_ctrl_ln, pch=19, col=leaders_col[2], ylim=c(0,0.25), xlab="naive_quant", ylab="% in quantile", main="early_ctrl_pd1")
lines(naive_by_p_spec_ctrl_ln, col=leaders_col[2])

points(naive_by_p_spec_pd1_ln, pch=21, col=leaders_col[2])
lines(naive_by_p_spec_pd1_ln, col=leaders_col[2], lty=3)
dev.off()

png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_conds_m/", "naive_traj_cc_cells_early_ctrl", ".png",sep=""), 8 * 100, 8 * 100, res = res_cc)
plot(cc_naive_by_p_spec_ctrl_spleen, pch=19, col=leaders_col[9], ylim=c(0,0.6), xlab="naive_quant", ylab="% cycling", main="early_ctrl_cc")
lines(cc_naive_by_p_spec_ctrl_spleen, col=leaders_col[9])
points(cc_naive_by_p_spec_ctrl_ln, pch=19, col=leaders_col[2])
lines(cc_naive_by_p_spec_ctrl_ln, col=leaders_col[2])
points(cc_naive_by_p_spec_ctrl_tum, pch=19, col=leaders_col[11])
lines(cc_naive_by_p_spec_ctrl_tum, col=leaders_col[11])
dev.off()

png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_conds_m/", "naive_traj_cc_cells_early_pd1", ".png",sep=""), 8 * 100, 8 * 100, res = res_cc)
plot(cc_naive_by_p_spec_pd1_spleen, pch=19, col=leaders_col[9], ylim=c(0,0.6), xlab="naive_quant", ylab="% cycling", main="early_pd1_cc")
lines(cc_naive_by_p_spec_pd1_spleen, col=leaders_col[9])
points(cc_naive_by_p_spec_pd1_ln, pch=19, col=leaders_col[2])
lines(cc_naive_by_p_spec_pd1_ln, col=leaders_col[2])
points(cc_naive_by_p_spec_pd1_tum, pch=19, col=leaders_col[11])
lines(cc_naive_by_p_spec_pd1_tum, col=leaders_col[11])
dev.off()

png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_conds_m/", "naive_traj_cc_cells_early_ctrl_pd1_spleen", ".png",sep=""), 8 * 100, 8 * 100, res = res_cc)
plot(cc_naive_by_p_spec_ctrl_spleen, pch=19, col=leaders_col[9], ylim=c(0,0.6), xlab="naive_quant", ylab="% cycling", main="early_ctrl_pd1_cc")
lines(cc_naive_by_p_spec_ctrl_spleen, col=leaders_col[9])

points(cc_naive_by_p_spec_pd1_spleen, pch=21, col=leaders_col[9])
lines(cc_naive_by_p_spec_pd1_spleen, col=leaders_col[9], lty=3)
dev.off()

png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_conds_m/", "naive_traj_cc_cells_early_ctrl_pd1_tum", ".png",sep=""), 8 * 100, 8 * 100, res = res_cc)
plot(cc_naive_by_p_spec_ctrl_tum, pch=19, col=leaders_col[11], ylim=c(0,0.6), xlab="naive_quant", ylab="% cycling", main="early_ctrl_pd1_cc")
lines(cc_naive_by_p_spec_ctrl_tum, col=leaders_col[11])

points(cc_naive_by_p_spec_pd1_tum, pch=21, col=leaders_col[11])
lines(cc_naive_by_p_spec_pd1_tum, col=leaders_col[11], lty=3)
dev.off()

png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_conds_m/", "naive_traj_cc_cells_early_ctrl_pd1_ln", ".png",sep=""), 8 * 100, 8 * 100, res = res_cc)
plot(cc_naive_by_p_spec_ctrl_ln, pch=19, col=leaders_col[2], ylim=c(0,0.6), xlab="naive_quant", ylab="% cycling", main="early_ctrl_pd1_cc")
lines(cc_naive_by_p_spec_ctrl_ln, col=leaders_col[2])

points(cc_naive_by_p_spec_pd1_ln, pch=21, col=leaders_col[2])
lines(cc_naive_by_p_spec_pd1_ln, col=leaders_col[2], lty=3)
dev.off()


naive_by_p_spec_41bb_tum = c()
cc_naive_by_p_spec_41bb_tum = c()

naive_by_p_spec_41bb_pd1_tum = c()
cc_naive_by_p_spec_41bb_pd1_tum = c()

min_cells = 20

day_naive = c(1:3)
for (quant in 1:length(naive_mcs_ord_ls)){
  active_mcs = naive_mcs_ord_ls[[quant]]
  active_cells = rel_cells[mc@mc[rel_cells] %in% active_mcs]
  
  # spec_cells_ctrl = rel_cells[is_spec & (is_tumor) & is_41bb & md$days_post_transfer %in% day_naive]
  # spec_cells_pd1 = rel_cells[is_spec & (is_tumor) & is_41bb_pd1 & md$days_post_transfer %in% day_naive]
  spec_cells_ctrl = rel_cells[is_byst & (is_ln) & is_ctrl & md$days_post_transfer %in% day_naive]
  spec_cells_pd1 = rel_cells[is_byst & (is_ln) & is_pd1 & md$days_post_transfer %in% day_naive]
  spec_cells_ctrl_active = spec_cells_ctrl[spec_cells_ctrl %in% active_cells]
  if (min_cells < length(spec_cells_ctrl_active)){
    naive_by_p_spec_41bb_tum = c(naive_by_p_spec_41bb_tum, length(spec_cells_ctrl_active)/length(spec_cells_ctrl))
    cc_naive_by_p_spec_41bb_tum = c(cc_naive_by_p_spec_41bb_tum, sum(spec_cells_ctrl_active %in% cc_pos_cells)/length(spec_cells_ctrl_active))
  }
  else{
    naive_by_p_spec_41bb_tum = c(naive_by_p_spec_41bb_tum, 0)
    cc_naive_by_p_spec_41bb_tum = c(cc_naive_by_p_spec_41bb_tum, 0)
  }
  
  spec_cells_pd1_active = spec_cells_pd1[spec_cells_pd1 %in% active_cells]
  if (min_cells < length(spec_cells_pd1_active)){
    naive_by_p_spec_41bb_pd1_tum = c(naive_by_p_spec_41bb_pd1_tum, length(spec_cells_pd1_active)/length(spec_cells_pd1))
    cc_naive_by_p_spec_41bb_pd1_tum = c(cc_naive_by_p_spec_41bb_pd1_tum, sum(spec_cells_pd1_active %in% cc_pos_cells)/length(spec_cells_pd1_active))
  }
  else{
    naive_by_p_spec_41bb_pd1_tum = c(naive_by_p_spec_41bb_pd1_tum, 0)
    cc_naive_by_p_spec_41bb_pd1_tum = c(cc_naive_by_p_spec_41bb_pd1_tum, 0)
  }
}


png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_conds_m/", "naive_traj_active_cells_early_41bb_pd1_tum", ".png",sep=""), 8 * 100, 8 * 100, res = res_cc)
plot(naive_by_p_spec_ctrl_tum, pch=19, col=leaders_col[11], ylim=c(0,0.25), xlab="naive_quant", ylab="% in quantile", main="early_41bb_pd1")
lines(naive_by_p_spec_ctrl_tum, col=leaders_col[11])

points(naive_by_p_spec_pd1_tum, pch=21, col=leaders_col[11])
lines(naive_by_p_spec_pd1_tum, col=leaders_col[11], lty=3)

points(naive_by_p_spec_41bb_pd1_tum, pch=21, col=leaders_col[3])
lines(naive_by_p_spec_41bb_pd1_tum, col=leaders_col[3], lty=3)

points(naive_by_p_spec_41bb_tum, pch=21, col=leaders_col[6])
lines(naive_by_p_spec_41bb_tum, col=leaders_col[6], lty=3)
dev.off()

png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_conds_m/", "naive_traj_cc_cells_early_41bb_pd1_tum", ".png",sep=""), 8 * 100, 8 * 100, res = res_cc)
plot(cc_naive_by_p_spec_ctrl_tum, pch=19, col=leaders_col[11], ylim=c(0,0.8), xlab="naive_quant", ylab="% cycling", main="early_41bb_pd1_cc")
lines(cc_naive_by_p_spec_ctrl_tum, col=leaders_col[11])

points(cc_naive_by_p_spec_pd1_tum, pch=21, col=leaders_col[11])
lines(cc_naive_by_p_spec_pd1_tum, col=leaders_col[11], lty=3)

points(cc_naive_by_p_spec_41bb_pd1_tum, pch=21, col=leaders_col[3])
lines(cc_naive_by_p_spec_41bb_pd1_tum, col=leaders_col[3], lty=3)
points(cc_naive_by_p_spec_41bb_tum, pch=21, col=leaders_col[6])
lines(cc_naive_by_p_spec_41bb_tum, col=leaders_col[6], lty=3)
dev.off()

png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_conds_m/", "naive_traj_active_cells_early_pd1_byst_ln", ".png",sep=""), 8 * 100, 8 * 100, res = res_cc)
plot(naive_by_p_spec_ctrl_ln, pch=19, col=leaders_col[2], ylim=c(0,0.3), xlab="naive_quant", ylab="% in quantile", main="early_byst_pd1")
lines(naive_by_p_spec_ctrl_ln, col=leaders_col[2])

points(naive_by_p_spec_pd1_ln, pch=21, col=leaders_col[2])
lines(naive_by_p_spec_pd1_ln, col=leaders_col[2], lty=3)

points(naive_by_p_spec_41bb_pd1_tum, pch=21, col=leaders_col[4])
lines(naive_by_p_spec_41bb_pd1_tum, col=leaders_col[4], lty=3)

points(naive_by_p_spec_41bb_tum, pch=21, col=leaders_col[4])
lines(naive_by_p_spec_41bb_tum, col=leaders_col[4])
dev.off()

png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_conds_m/", "naive_traj_cc_cells_early_pd1_byst_ln", ".png",sep=""), 8 * 100, 8 * 100, res = res_cc)
plot(cc_naive_by_p_spec_ctrl_ln, pch=19, col=leaders_col[2], ylim=c(0,0.8), xlab="naive_quant", ylab="% cycling", main="early_byst_pd1_cc")
lines(cc_naive_by_p_spec_ctrl_ln, col=leaders_col[2])

points(cc_naive_by_p_spec_pd1_ln, pch=21, col=leaders_col[2])
lines(cc_naive_by_p_spec_pd1_ln, col=leaders_col[2], lty=3)

points(cc_naive_by_p_spec_41bb_pd1_tum, pch=21, col=leaders_col[4])
lines(cc_naive_by_p_spec_41bb_pd1_tum, col=leaders_col[4], lty=3)
points(cc_naive_by_p_spec_41bb_tum, pch=21, col=leaders_col[4])
lines(cc_naive_by_p_spec_41bb_tum, col=leaders_col[4])
dev.off()






naive_by_p_spec_ctrl_tum = c()
naive_by_p_spec_ctrl_ln = c()
naive_by_p_spec_ctrl_spleen = c()
cc_naive_by_p_spec_ctrl_tum = c()
cc_naive_by_p_spec_ctrl_ln = c()
cc_naive_by_p_spec_ctrl_spleen = c()

naive_by_p_spec_pd1_tum = c()
naive_by_p_spec_pd1_ln = c()
naive_by_p_spec_pd1_spleen = c()
cc_naive_by_p_spec_pd1_tum = c()
cc_naive_by_p_spec_pd1_ln = c()
cc_naive_by_p_spec_pd1_spleen = c()

time_genes = c("Tcf7", "Sell", "Tnfrsf9", "Prf1", "Pdcd1", "Cxcr6", "Irf8", "Zbtb2", "Id2", "Id3", "Tigit", "Myc", "Lag3", "Havcr2", "Prdm1", "Il2ra", "Il2rb", "Ccr7")
# downs_tum_spec = downs_mat[,md[colnames(downs_mat),"location"] == "tumor" & md[colnames(downs_mat),"cell_type"] %in% c("spec", "cd8_in_vitro") & as.character(md[colnames(downs_mat),"condition"]) == "ctrl"]
# downs_tum_spec_genes = downs_tum_spec[time_genes,]

day_naive = c(4:6)
for (quant in 1:length(naive_mcs_ord_ls)){
  active_mcs = naive_mcs_ord_ls[[quant]]
  active_cells = rel_cells[mc@mc[rel_cells] %in% active_mcs]
  
  spec_cells_ctrl = rel_cells[is_spec & (is_tumor) & is_ctrl & md$days_post_transfer %in% day_naive]
  spec_cells_pd1 = rel_cells[is_spec & (is_tumor) & is_pd1 & md$days_post_transfer %in% day_naive]
  spec_cells_ctrl_active = spec_cells_ctrl[spec_cells_ctrl %in% active_cells]
  if (min_cells < length(spec_cells_ctrl_active)){
    naive_by_p_spec_ctrl_tum = c(naive_by_p_spec_ctrl_tum, length(spec_cells_ctrl_active)/length(spec_cells_ctrl))
    cc_naive_by_p_spec_ctrl_tum = c(cc_naive_by_p_spec_ctrl_tum, sum(spec_cells_ctrl_active %in% cc_pos_cells)/length(spec_cells_ctrl_active))
  }
  else{
    naive_by_p_spec_ctrl_tum = c(naive_by_p_spec_ctrl_tum, 0)
    cc_naive_by_p_spec_ctrl_tum = c(cc_naive_by_p_spec_ctrl_tum, 0)
  }
  
  spec_cells_pd1_active = spec_cells_pd1[spec_cells_pd1 %in% active_cells]
  if (min_cells < length(spec_cells_pd1_active)){
    naive_by_p_spec_pd1_tum = c(naive_by_p_spec_pd1_tum, length(spec_cells_pd1_active)/length(spec_cells_pd1))
    cc_naive_by_p_spec_pd1_tum = c(cc_naive_by_p_spec_pd1_tum, sum(spec_cells_pd1_active %in% cc_pos_cells)/length(spec_cells_pd1_active))
  }
  else{
    naive_by_p_spec_pd1_tum = c(naive_by_p_spec_pd1_tum, 0)
    cc_naive_by_p_spec_pd1_tum = c(cc_naive_by_p_spec_pd1_tum, 0)
  }
  
  
  spec_cells_ctrl = rel_cells[is_spec & (is_ln) & is_ctrl & md$days_post_transfer %in% day_naive]
  spec_cells_pd1 = rel_cells[is_spec & (is_ln) & is_pd1 & md$days_post_transfer %in% day_naive]
  spec_cells_ctrl_active = spec_cells_ctrl[spec_cells_ctrl %in% active_cells]
  if (min_cells < length(spec_cells_ctrl_active)){
    naive_by_p_spec_ctrl_ln = c(naive_by_p_spec_ctrl_ln, length(spec_cells_ctrl_active)/length(spec_cells_ctrl))
    cc_naive_by_p_spec_ctrl_ln = c(cc_naive_by_p_spec_ctrl_ln, sum(spec_cells_ctrl_active %in% cc_pos_cells)/length(spec_cells_ctrl_active))
  }
  else{
    naive_by_p_spec_ctrl_ln = c(naive_by_p_spec_ctrl_ln, 0)
    cc_naive_by_p_spec_ctrl_ln = c(cc_naive_by_p_spec_ctrl_ln, 0)
  }
  spec_cells_pd1_active = spec_cells_pd1[spec_cells_pd1 %in% active_cells]
  if (min_cells < length(spec_cells_pd1_active)){
    naive_by_p_spec_pd1_ln = c(naive_by_p_spec_pd1_ln, length(spec_cells_pd1_active)/length(spec_cells_pd1))
    cc_naive_by_p_spec_pd1_ln = c(cc_naive_by_p_spec_pd1_ln, sum(spec_cells_pd1_active %in% cc_pos_cells)/length(spec_cells_pd1_active))
  }
  else{
    naive_by_p_spec_pd1_ln = c(naive_by_p_spec_pd1_ln, 0)
    cc_naive_by_p_spec_pd1_ln = c(cc_naive_by_p_spec_pd1_ln, 0)
  }
  
  
  spec_cells_ctrl = rel_cells[is_spec & (is_spleen) & is_ctrl & md$days_post_transfer %in% day_naive]
  spec_cells_pd1 = rel_cells[is_spec & (is_spleen) & is_pd1 & md$days_post_transfer %in% day_naive]
  spec_cells_ctrl_active = spec_cells_ctrl[spec_cells_ctrl %in% active_cells]
  if (min_cells < length(spec_cells_ctrl_active)){
    naive_by_p_spec_ctrl_spleen = c(naive_by_p_spec_ctrl_spleen, length(spec_cells_ctrl_active)/length(spec_cells_ctrl))
    cc_naive_by_p_spec_ctrl_spleen = c(cc_naive_by_p_spec_ctrl_spleen, sum(spec_cells_ctrl_active %in% cc_pos_cells)/length(spec_cells_ctrl_active))
  }
  else{
    naive_by_p_spec_ctrl_spleen = c(naive_by_p_spec_ctrl_spleen, 0)
    cc_naive_by_p_spec_ctrl_spleen = c(cc_naive_by_p_spec_ctrl_spleen, 0)
  }
  
  spec_cells_pd1_active = spec_cells_pd1[spec_cells_pd1 %in% active_cells]
  if (min_cells < length(spec_cells_ctrl_active)){
    naive_by_p_spec_pd1_spleen = c(naive_by_p_spec_pd1_spleen, length(spec_cells_pd1_active)/length(spec_cells_pd1))
    cc_naive_by_p_spec_pd1_spleen = c(cc_naive_by_p_spec_pd1_spleen, sum(spec_cells_pd1_active %in% cc_pos_cells)/length(spec_cells_pd1_active))
  }
  else{
    naive_by_p_spec_pd1_spleen = c(naive_by_p_spec_pd1_spleen, 0)
    cc_naive_by_p_spec_pd1_spleen = c(cc_naive_by_p_spec_pd1_spleen, 0)
  }
  
}

png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_conds_m/", "naive_traj_active_cells_middle_ctrl", ".png",sep=""), 8 * 100, 8 * 100, res = res_cc)
plot(naive_by_p_spec_ctrl_spleen, pch=19, col=leaders_col[9], ylim=c(0,0.25), xlab="naive_quant", ylab="% in quantile", main="middle_ctrl")
lines(naive_by_p_spec_ctrl_spleen, col=leaders_col[9])
points(naive_by_p_spec_ctrl_ln, pch=19, col=leaders_col[2])
lines(naive_by_p_spec_ctrl_ln, col=leaders_col[2])
points(naive_by_p_spec_ctrl_tum, pch=19, col=leaders_col[11])
lines(naive_by_p_spec_ctrl_tum, col=leaders_col[11])
dev.off()

png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_conds_m/", "naive_traj_active_cells_middle_pd1", ".png",sep=""), 8 * 100, 8 * 100, res = res_cc)
plot(naive_by_p_spec_pd1_spleen, pch=19, col=leaders_col[9], ylim=c(0,0.25), xlab="naive_quant", ylab="% in quantile", main="middle_pd1")
lines(naive_by_p_spec_pd1_spleen, col=leaders_col[9])
points(naive_by_p_spec_pd1_ln, pch=19, col=leaders_col[2])
lines(naive_by_p_spec_pd1_ln, col=leaders_col[2])
points(naive_by_p_spec_pd1_tum, pch=19, col=leaders_col[11])
lines(naive_by_p_spec_pd1_tum, col=leaders_col[11])
dev.off()


png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_conds_m/", "naive_traj_active_cells_middle_ctrl_pd1_spleen", ".png",sep=""), 8 * 100, 8 * 100, res = res_cc)
plot(naive_by_p_spec_ctrl_spleen, pch=19, col=leaders_col[9], ylim=c(0,0.25), xlab="naive_quant", ylab="% in quantile", main="middle_ctrl_pd1")
lines(naive_by_p_spec_ctrl_spleen, col=leaders_col[9])

points(naive_by_p_spec_pd1_spleen, pch=21, col=leaders_col[9])
lines(naive_by_p_spec_pd1_spleen, col=leaders_col[9], lty=3)
dev.off()

png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_conds_m/", "naive_traj_active_cells_middle_ctrl_pd1_tum", ".png",sep=""), 8 * 100, 8 * 100, res = res_cc)
plot(naive_by_p_spec_ctrl_tum, pch=19, col=leaders_col[11], ylim=c(0,0.25), xlab="naive_quant", ylab="% in quantile", main="middle_ctrl_pd1")
lines(naive_by_p_spec_ctrl_tum, col=leaders_col[1])

points(naive_by_p_spec_pd1_tum, pch=21, col=leaders_col[11])
lines(naive_by_p_spec_pd1_tum, col=leaders_col[11], lty=3)
dev.off()

png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_conds_m/", "naive_traj_active_cells_middle_ctrl_pd1_ln", ".png",sep=""), 8 * 100, 8 * 100, res = res_cc)
plot(naive_by_p_spec_ctrl_ln, pch=19, col=leaders_col[2], ylim=c(0,0.25), xlab="naive_quant", ylab="% in quantile", main="middle_ctrl_pd1")
lines(naive_by_p_spec_ctrl_ln, col=leaders_col[2])

points(naive_by_p_spec_pd1_ln, pch=21, col=leaders_col[2])
lines(naive_by_p_spec_pd1_ln, col=leaders_col[2], lty=3)
dev.off()



png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_conds_m/", "naive_traj_cc_cells_middle_ctrl", ".png",sep=""), 8 * 100, 8 * 100, res = res_cc)
plot(cc_naive_by_p_spec_ctrl_spleen, pch=19, col=leaders_col[9], ylim=c(0,0.6), xlab="naive_quant", ylab="% cycling", main="middle_ctrl_cc")
lines(cc_naive_by_p_spec_ctrl_spleen, col=leaders_col[9])
points(cc_naive_by_p_spec_ctrl_ln, pch=19, col=leaders_col[2])
lines(cc_naive_by_p_spec_ctrl_ln, col=leaders_col[2])
points(cc_naive_by_p_spec_ctrl_tum, pch=19, col=leaders_col[11])
lines(cc_naive_by_p_spec_ctrl_tum, col=leaders_col[11])
dev.off()

png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_conds_m/", "naive_traj_cc_cells_middle_pd1", ".png",sep=""), 8 * 100, 8 * 100, res = res_cc)
plot(cc_naive_by_p_spec_pd1_spleen, pch=19, col=leaders_col[9], ylim=c(0,0.6), xlab="naive_quant", ylab="% cycling", main="middle_pd1_cc")
lines(cc_naive_by_p_spec_pd1_spleen, col=leaders_col[9])
points(cc_naive_by_p_spec_pd1_ln, pch=19, col=leaders_col[2])
lines(cc_naive_by_p_spec_pd1_ln, col=leaders_col[2])
points(cc_naive_by_p_spec_pd1_tum, pch=19, col=leaders_col[11])
lines(cc_naive_by_p_spec_pd1_tum, col=leaders_col[11])
dev.off()


png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_conds_m/", "naive_traj_cc_cells_middle_ctrl_pd1_spleen", ".png",sep=""), 8 * 100, 8 * 100, res = res_cc)
plot(cc_naive_by_p_spec_ctrl_spleen, pch=19, col=leaders_col[9], ylim=c(0,0.6), xlab="naive_quant", ylab="% cycling", main="middle_ctrl_pd1_cc")
lines(cc_naive_by_p_spec_ctrl_spleen, col=leaders_col[9])

points(cc_naive_by_p_spec_pd1_spleen, pch=21, col=leaders_col[9])
lines(cc_naive_by_p_spec_pd1_spleen, col=leaders_col[9], lty=3)
dev.off()

png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_conds_m/", "naive_traj_cc_cells_middle_ctrl_pd1_tum", ".png",sep=""), 8 * 100, 8 * 100, res = res_cc)
plot(cc_naive_by_p_spec_ctrl_tum, pch=19, col=leaders_col[11], ylim=c(0,0.6), xlab="naive_quant", ylab="% cycling", main="middle_ctrl_pd1_cc")
lines(cc_naive_by_p_spec_ctrl_tum, col=leaders_col[1])

points(cc_naive_by_p_spec_pd1_tum, pch=21, col=leaders_col[11])
lines(cc_naive_by_p_spec_pd1_tum, col=leaders_col[11], lty=3)
dev.off()

png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_conds_m/", "naive_traj_cc_cells_middle_ctrl_pd1_ln", ".png",sep=""), 8 * 100, 8 * 100, res = res_cc)
plot(cc_naive_by_p_spec_ctrl_ln, pch=19, col=leaders_col[2], ylim=c(0,0.6), xlab="naive_quant", ylab="% cycling", main="middle_ctrl_pd1_cc")
lines(cc_naive_by_p_spec_ctrl_ln, col=leaders_col[2])

points(cc_naive_by_p_spec_pd1_ln, pch=21, col=leaders_col[2])
lines(cc_naive_by_p_spec_pd1_ln, col=leaders_col[2], lty=3)
dev.off()


naive_by_p_spec_41bb_tum = c()
cc_naive_by_p_spec_41bb_tum = c()

naive_by_p_spec_41bb_pd1_tum = c()
cc_naive_by_p_spec_41bb_pd1_tum = c()

min_cells = 20

day_naive = c(4:6)
for (quant in 1:length(naive_mcs_ord_ls)){
  active_mcs = naive_mcs_ord_ls[[quant]]
  active_cells = rel_cells[mc@mc[rel_cells] %in% active_mcs]
  
  # spec_cells_ctrl = rel_cells[is_spec & (is_tumor) & is_41bb & md$days_post_transfer %in% day_naive]
  # spec_cells_pd1 = rel_cells[is_spec & (is_tumor) & is_41bb_pd1 & md$days_post_transfer %in% day_naive]
  spec_cells_ctrl = rel_cells[is_byst & (is_ln) & is_ctrl & md$days_post_transfer %in% day_naive]
  spec_cells_pd1 = rel_cells[is_byst & (is_ln) & is_pd1 & md$days_post_transfer %in% day_naive]
  
  spec_cells_ctrl_active = spec_cells_ctrl[spec_cells_ctrl %in% active_cells]
  if (min_cells < length(spec_cells_ctrl_active)){
    # print(quant)
    # print(length(spec_cells_ctrl_active))
    naive_by_p_spec_41bb_tum = c(naive_by_p_spec_41bb_tum, length(spec_cells_ctrl_active)/length(spec_cells_ctrl))
    cc_naive_by_p_spec_41bb_tum = c(cc_naive_by_p_spec_41bb_tum, sum(spec_cells_ctrl_active %in% cc_pos_cells)/length(spec_cells_ctrl_active))
  }
  else{
    naive_by_p_spec_41bb_tum = c(naive_by_p_spec_41bb_tum, 0)
    cc_naive_by_p_spec_41bb_tum = c(cc_naive_by_p_spec_41bb_tum, 0)
  }
  
  spec_cells_pd1_active = spec_cells_pd1[spec_cells_pd1 %in% active_cells]
  if (min_cells < length(spec_cells_pd1_active)){
    print(quant)
    print(length(spec_cells_ctrl_active))
    naive_by_p_spec_41bb_pd1_tum = c(naive_by_p_spec_41bb_pd1_tum, length(spec_cells_pd1_active)/length(spec_cells_pd1))
    cc_naive_by_p_spec_41bb_pd1_tum = c(cc_naive_by_p_spec_41bb_pd1_tum, sum(spec_cells_pd1_active %in% cc_pos_cells)/length(spec_cells_pd1_active))
  }
  else{
    naive_by_p_spec_41bb_pd1_tum = c(naive_by_p_spec_41bb_pd1_tum, 0)
    cc_naive_by_p_spec_41bb_pd1_tum = c(cc_naive_by_p_spec_41bb_pd1_tum, 0)
  }
}


png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_conds_m/", "naive_traj_active_cells_middle_41bb_pd1_tum", ".png",sep=""), 8 * 100, 8 * 100, res = res_cc)
plot(naive_by_p_spec_ctrl_tum, pch=19, col=leaders_col[11], ylim=c(0,0.25), xlab="naive_quant", ylab="% in quantile", main="middle_41bb_pd1")
lines(naive_by_p_spec_ctrl_tum, col=leaders_col[1])

points(naive_by_p_spec_pd1_tum, pch=21, col=leaders_col[11])
lines(naive_by_p_spec_pd1_tum, col=leaders_col[11], lty=3)

points(naive_by_p_spec_41bb_pd1_tum, pch=21, col=leaders_col[3])
lines(naive_by_p_spec_41bb_pd1_tum, col=leaders_col[3], lty=3)

points(naive_by_p_spec_41bb_tum, pch=21, col=leaders_col[6])
lines(naive_by_p_spec_41bb_tum, col=leaders_col[6], lty=3)
dev.off()

png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_conds_m/", "naive_traj_cc_cells_middle_41bb_pd1_tum", ".png",sep=""), 8 * 100, 8 * 100, res = res_cc)
plot(cc_naive_by_p_spec_ctrl_tum, pch=19, col=leaders_col[11], ylim=c(0,0.8), xlab="naive_quant", ylab="% cycling", main="middle_41bb_pd1_cc")
lines(cc_naive_by_p_spec_ctrl_tum, col=leaders_col[1])

points(cc_naive_by_p_spec_pd1_tum, pch=21, col=leaders_col[11])
lines(cc_naive_by_p_spec_pd1_tum, col=leaders_col[11], lty=3)

points(cc_naive_by_p_spec_41bb_pd1_tum, pch=21, col=leaders_col[3])
lines(cc_naive_by_p_spec_41bb_pd1_tum, col=leaders_col[3], lty=3)
points(cc_naive_by_p_spec_41bb_tum, pch=21, col=leaders_col[6])
lines(cc_naive_by_p_spec_41bb_tum, col=leaders_col[6], lty=3)
dev.off()

png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_conds_m/", "naive_traj_active_cells_middle_pd1_byst_ln", ".png",sep=""), 8 * 100, 8 * 100, res = res_cc)
plot(naive_by_p_spec_ctrl_ln, pch=19, col=leaders_col[2], ylim=c(0,0.3), xlab="naive_quant", ylab="% in quantile", main="middle_byst_pd1")
lines(naive_by_p_spec_ctrl_ln, col=leaders_col[2])

points(naive_by_p_spec_pd1_ln, pch=21, col=leaders_col[2])
lines(naive_by_p_spec_pd1_ln, col=leaders_col[2], lty=3)

points(naive_by_p_spec_41bb_pd1_tum, pch=21, col=leaders_col[4])
lines(naive_by_p_spec_41bb_pd1_tum, col=leaders_col[4], lty=3)

points(naive_by_p_spec_41bb_tum, pch=21, col=leaders_col[4])
lines(naive_by_p_spec_41bb_tum, col=leaders_col[4])
dev.off()

png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_conds_m/", "naive_traj_cc_cells_middle_pd1_byst_ln", ".png",sep=""), 8 * 100, 8 * 100, res = res_cc)
plot(cc_naive_by_p_spec_ctrl_ln, pch=19, col=leaders_col[2], ylim=c(0,0.8), xlab="naive_quant", ylab="% cycling", main="middle_byst_pd1_cc")
lines(cc_naive_by_p_spec_ctrl_ln, col=leaders_col[2])

points(cc_naive_by_p_spec_pd1_ln, pch=21, col=leaders_col[2])
lines(cc_naive_by_p_spec_pd1_ln, col=leaders_col[2], lty=3)

points(cc_naive_by_p_spec_41bb_pd1_tum, pch=21, col=leaders_col[4])
lines(cc_naive_by_p_spec_41bb_pd1_tum, col=leaders_col[4], lty=3)
points(cc_naive_by_p_spec_41bb_tum, pch=21, col=leaders_col[4])
lines(cc_naive_by_p_spec_41bb_tum, col=leaders_col[4])
dev.off()





png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_conds_m/", "naive_legend.png",sep=""), 14 * 100, 14 * 100, res = 200)
plot.new()
legend("topleft", legend=c("Spleen", "LN", "Tumor"), pch=19, col=c(leaders_col[9], leaders_col[2], leaders_col[11]), cex=2)
dev.off()

png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_conds_m/", "naive_treat_legend.png",sep=""), 14 * 100, 14 * 100, res = 200)
plot.new()
legend("topleft", legend=c("Ctrl", "aPD1"), pch=c(19,21), col=c(leaders_col[11]), cex=2)
dev.off()


# naive_score = (norm_global_mat[,9] + norm_global_mat[,2]) / (norm_global_mat[,9] + norm_global_mat[,2] + norm_global_mat[,11] + norm_global_mat[,1]) 

prog_usage =naive_score
max_usage = max(prog_usage)
min_usage = min(prog_usage)
x = pmin(pmax(prog_usage, min_usage), max_usage) - min_usage
colspec_2d = c("navyblue", "white", "red")
shades = colorRampPalette(colspec_2d)(100 * (max_usage-min_usage) + 1)
mc_cols = shades[round(100 * x) + 1]
# mc_cols[!(norm_global_mat[,2] >= active_thres_vec[2] | norm_global_mat[,9] >= active_thres_vec[9] | norm_global_mat[,11] >= active_thres_vec[11] | norm_global_mat[,1] >= active_thres_vec[1])] = "grey"
naive_mask = (norm_global_mat[,2] >= active_thres_vec[2] | norm_global_mat[,9] >= active_thres_vec[9] | norm_global_mat[,11] >= active_thres_vec[11] | norm_global_mat[,1] >= active_thres_vec[1])
mc_cols[!naive_mask] = "grey"
# png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/", "umap_prog_full_", loc_nm, ".png",sep=""), 12 * 100, 12 * 100, res = 200)
par(mar=c(0,0,0,0))
plot(mc2d_2@mc_x[naive_mask], mc2d_2@mc_y[naive_mask], bg=mc_cols[naive_mask], pch = 21, col="black")
plot(umap_full[naive_mask,1], umap_full[naive_mask,2], bg=mc_cols[naive_mask], pch = 21, col="black")
# dev.off()


prog=4
spec_6_41bb = rbind(mouse_cols_spec[colnames(state_by_m_spec_6)], state_by_m_spec_6[prog,])
colnames(spec_6_41bb) = colnames(state_by_m_spec_6)
rownames(spec_6_41bb) = c("cond", "%Active")
write.csv(spec_6_41bb, file = "spec_6_trnfrsf9.csv")

prog=6
spec_6_41bb = rbind(mouse_cols_spec[colnames(state_by_m_spec_6)], state_by_m_spec_6[prog,])
colnames(spec_6_41bb) = colnames(state_by_m_spec_6)
rownames(spec_6_41bb) = c("cond", "%Active")
write.csv(spec_6_41bb, file = "spec_6_tcf7.csv")

prog=5
spec_6_41bb = rbind(mouse_cols_spec[colnames(state_by_m_spec_6)], state_by_m_spec_6[prog,])
colnames(spec_6_41bb) = colnames(state_by_m_spec_6)
rownames(spec_6_41bb) = c("cond", "%Active")
write.csv(spec_6_41bb, file = "spec_6_ccl3.csv")

for (prog in 1:ncol(norm_global_w)){
  ylim_max = max(c(0.4, max(state_by_m_spec_3[prog,]), max(state_by_m_mc38_3[prog,]), max(state_by_m_b16_3[prog,]), max(state_by_m_spec_6[prog,]), max(state_by_m_mc38_6[prog,]), max(state_by_m_b16_6[prog,])))
  png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_conds_m/", "scatter_active_",prog, "_","spec", "_time_", "6", ".png",sep=""), 8 * 100, 8 * 100, res = 200)
  plot(0.2*runif(ncol(state_by_m_spec_6)) + col2ind[mouse_cols_spec[colnames(state_by_m_spec_6)]], state_by_m_spec_6[prog,], col=mouse_cols_spec[colnames(state_by_m_spec_6)], pch=19, xaxt = "n", ylim=c(0,ylim_max), main=paste("spec 6 p=", names(leaders_col)[prog], sep=""), xlab="Cond", ylab="% Active")
  axis(1, at=1:4, labels=c("ctrl", "pd1", "41bb", "41bb_pd1"))
  dev.off()
  
  png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_conds_m/", "scatter_active_",prog, "_","mc38", "_time_", "6", ".png",sep=""), 8 * 100, 8 * 100, res = 200)
  plot(0.1*runif(ncol(state_by_m_mc38_6)) + col2ind[mouse_cols_mc38[colnames(state_by_m_mc38_6)]], state_by_m_mc38_6[prog,], col=mouse_cols_mc38[colnames(state_by_m_mc38_6)], pch=19, xaxt = "n", ylim=c(0,ylim_max), main=paste("mc38 6 p=", names(leaders_col)[prog], sep=""), xlab="Cond", ylab="% Active")
  axis(1, at=1:4, labels=c("ctrl", "pd1", "41bb", "41bb_pd1"))
  dev.off()
  
  png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_conds_m/", "scatter_active_",prog, "_","b16", "_time_", "6", ".png",sep=""), 8 * 100, 8 * 100, res = 200)
  plot(0.1*runif(ncol(state_by_m_b16_6)) + col2ind[mouse_cols_b16[colnames(state_by_m_b16_6)]], state_by_m_b16_6[prog,], col=mouse_cols_b16[colnames(state_by_m_b16_6)], pch=19, xaxt = "n", ylim=c(0,ylim_max), main=paste("b16 6 p=", names(leaders_col)[prog], sep=""), xlab="Cond", ylab="% Active")
  axis(1, at=1:4, labels=c("ctrl", "pd1", "41bb", "41bb_pd1"))
  dev.off()
}

for (prog in 1:ncol(norm_global_w)){
  png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_conds_m/", "barplot_active_",prog, "_","spec", "_time_", "6", ".png",sep=""), 12 * 100, 12 * 100, res = 200)
  barplot(state_by_m_spec_6[prog,], horiz =F, col = mouse_cols_spec[colnames(state_by_m_spec_6)], main=paste("spec 6 p=", names(leaders_col)[prog], sep=""), ylim=c(0,max(c(0.6, max(state_by_m_spec_6[prog,])))), las=2)
  dev.off()
  
  png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_conds_m/", "barplot_active_",prog, "_","mc38", "_time_", "6", ".png",sep=""), 12 * 100, 12 * 100, res = 200)
  barplot(state_by_m_mc38_6[prog,], horiz =F, col = mouse_cols_mc38[colnames(state_by_m_mc38_6)], main=paste("mc38 6 p=", names(leaders_col)[prog], sep=""), ylim=c(0,max(c(0.6, max(state_by_m_mc38_6[prog,])))), las=2)
  dev.off()
  
  png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/prog_conds_m/", "barplot_active_",prog, "_","b16", "_time_", "6", ".png",sep=""), 12 * 100, 12 * 100, res = 200)
  barplot(state_by_m_b16_6[prog,], horiz =F, col = mouse_cols_b16[colnames(state_by_m_b16_6)], main=paste("b16 6 p=", names(leaders_col)[prog], sep=""),ylim=c(0,max(c(0.6, max(state_by_m_b16_6[prog,])))), las=2)
  dev.off()
}
##### local modeling #####
full_model_list$reg_full_w[,3] = full_model_list$reg_full_w_non_filt[,3]
full_model_list$reg_full_w[full_model_list$reg_full_w[,3]<9441.864,3] = 0
#full_model_list$reg_full_w[mc@colors == "#cee6b9",3] = 0


for (ii in 1:length(full_model_list$overlap_list)){
  full_model_list$overlap_list[[ii]] = c(ii, full_model_list$overlap_list[[ii]])
}
full_model_list$overlap_list[[3]] = c(1,2,4,5,6,7,8,10)
start_ovlps = rep(1, top_n_full)
start_ovlps[1] = 5
start_ovlps[4] = 5
start_ovlps[3] = 9

root_env = rep(F, top_n_full)
root_env[6] = T

local_models = mcell_gene_prog_local_modeling(full_model_list, mc_id, mat_i, mc2d_id = NULL,cells=NULL, prog_ord=NULL, nmf_dir="nmf_plots_mc2_burst/", env_dir="", global_add_prog=1)
dev.off()

local_models = mcell_gene_prog_local_modeling_plots(local_models, full_model_list, mc_id, mat_i, mc2d_id=NULL, cells=NULL, prog_ord=NULL, nmf_dir=full_model_list$nmf_dir, env_dir="", plots=F, sim_plots=F, sim_progs=NULL, cd34_prog=mpp_neig, main_prog_cols=NULL, main_gset=NULL, tf_list=tf_exp, start_order_ovlps = NULL, root_env=NULL, filter_end = c(), lead_cor_thresh=0.35, lead_col=T, tfs_vel=F, tf_cor_thresh=0.1, annot_colors=NULL, choose_lead=NULL)


mpp_neig = 5
km_list_2 = list()
for (mpp_neig in 1:11){
  local_model = local_models[[mpp_neig]]
  local_w = local_models[[mpp_neig]]$local_w
  local_h = local_models[[mpp_neig]]$local_h
  norm_local_w = local_w / rowSums(local_w)
  max_norm_local_w = t(t(norm_local_w) / apply(norm_local_w, 2, max))
  norm_reconst_mat = norm_local_w %*% local_h
  
  n_clusts = ceiling(ncol(local_w) * 1.5)
  km_loc_5 = tglkmeans::TGL_kmeans(norm_local_w, n_clusts, id_column=F)
  names(km_loc_5$cluster) = rownames(norm_local_w)
  km_list_2[[mpp_neig]] = km_loc_5
  
  cell_clusts = km_loc_5$cluster
  clust_cols = table(mc@colors[mc@mc[names(cell_clusts)]], cell_clusts)
  clust_cols_chosen = rownames(clust_cols)[apply(clust_cols/rowSums(clust_cols), 2, which.max)]
  names(clust_cols_chosen) = colnames(clust_cols)
  paletteLength <- 500
  myColor <- colorRampPalette(c("white", "royalblue", "red"))(paletteLength)#colorRampPalette(c("blue","green", "white", "orange","red"))(paletteLength)
  myBreaks <- c(seq(0, 1, length.out=floor(paletteLength)))
  cells_clust_progs = norm_local_w[names(sort(cell_clusts)),]
  colnames(cells_clust_progs) = c("prog_lead", paste0("overlap_prog_", 2:ncol(cells_clust_progs)))
  
  prog_ann_loc = data.frame(row.names=rownames(cells_clust_progs),clust=names(clust_cols_chosen)[cell_clusts[rownames(cells_clust_progs)]])
  png(paste(scfigs_dir,"/nmf_plots_mc2_burst/cell_clust_double_heatmap_",mpp_neig, ".png",sep=""), 16 * 100, 16 * 100, res = 300)
  pheatmap(cells_clust_progs, cluster_rows = F, cluster_cols = T, color=myColor, breaks=myBreaks, annotation_row = prog_ann_loc, annotation_colors=list(clust=clust_cols_chosen))
  dev.off()
}

umap_list_2 = list()
for (mpp_neig in 1:11){
  local_model = local_models[[mpp_neig]]
  local_w = local_models[[mpp_neig]]$local_w
  local_h = local_models[[mpp_neig]]$local_h
  norm_local_w = local_w / rowSums(local_w)
  max_norm_local_w = t(t(norm_local_w) / apply(norm_local_w, 2, max))
  norm_reconst_mat = norm_local_w %*% local_h
  
  km_loc_5 = km_list_2[[mpp_neig]]
  
  cell_clusts = km_loc_5$cluster
  clust_cols = table(mc@colors[mc@mc[names(cell_clusts)]], cell_clusts)
  clust_cols_chosen = rownames(clust_cols)[apply(clust_cols/rowSums(clust_cols), 2, which.max)]
  names(clust_cols_chosen) = colnames(clust_cols)
  
  dist_w_loc = as.matrix(dist(norm_local_w))
  sigma_k_loc = 0.5
  sim_w_loc = exp((-dist_w_loc^2)/sigma_k_loc^2)
  dist_w_loc_clean = max(sim_w_loc) - sim_w_loc
  diag(dist_w_loc_clean) = 0
  umap_obj = umap_py$UMAP(n_neighbors = as.integer(15), min_dist=0.8, random_state = as.integer(1009), metric="precomputed") #, metric="manhattan", 
  umap_loc_env = umap_obj$fit_transform(dist_w_loc)
  rownames(umap_loc_env) = rownames(norm_local_w)
  png(paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_", mpp_neig, ".png",sep=""), 12 * 100, 12 * 100, res = 300)
  plot(umap_loc_env[,1], umap_loc_env[,2], pch=19,col=clust_cols_chosen[cell_clusts[rownames(umap_loc_env)]], cex=0.3)
  dev.off()
  umap_list_2[[mpp_neig]] = umap_loc_env
  
}


per_list_spec = vector("list",11)

for (mpp_neig in c(1,2,3,5,6,8,9,10,11)){
  umap_loc_env = umap_list_2[[mpp_neig]]
  m_id = mat@cell_metadata[rownames(umap_loc_env),"mouse_id"]
  all_mouse_id = names(table(m_id))[table(m_id) > 0]
  b_id = mat@cell_metadata[rownames(umap_loc_env),"batch_set_id"]
  all_batch_id = names(table(b_id))[table(b_id) > 0]
  km_loc_5 = km_list_2[[mpp_neig]]
  
  cell_clusts = km_loc_5$cluster
  clust_cols = table(mc@colors[mc@mc[names(cell_clusts)]], cell_clusts)
  clust_cols_chosen = rownames(clust_cols)[apply(clust_cols/rowSums(clust_cols), 2, which.max)]
  names(clust_cols_chosen) = colnames(clust_cols)
  
  rel_cells = rownames(umap_loc_env)
  rel_md = md[names(mc@mc)[names(mc@mc) %in% rel_cells],]
  rel_md = rel_md[rel_md$condition != "na",]
  rel_cells = rownames(rel_md)
  rel_mice = names(table(rel_md$mouse_id)[table(rel_md$mouse_id)>0])
  m_ct = table(rel_md$mouse_id, rel_md$cell_type)[rel_mice,]
  rel_mice = rel_mice[m_ct[,"spec"]>0]
  rel_md = rel_md[rel_md$mouse_id %in% rel_mice,]
  rel_cells = rownames(rel_md)
  
  colum_order_by = c('location', 'days_post_transfer', 'condition', 'mouse_id')
  unique_md = unique(rel_md[, colum_order_by])
  unique_md = unique_md[order(as.integer(as.character(unique_md$days_post_transfer)), unique_md$condition, unique_md$mouse_id, method='shell'), ]
  rel_mice_ord = as.character(unique(unique_md$mouse_id))
  rel_mice = rel_mice_ord
  
  dir.create(paste0("figs/nmf_plots_mc2_burst/2d_mice_traj_new_local_mcs_comb/spec_",mpp_neig), recursive = T)
  per_list_spec[[mpp_neig]] = plot_2d_mouse_traj_local_mcs(umap_loc_env, md, paste0("figs/nmf_plots_mc2_burst/2d_mice_traj_new_local_mcs_comb/spec_",mpp_neig), col_by_ct=T,col_by_clust=F, clust_cols_chosen=clust_cols_chosen, cell_clusts=cell_clusts, ct="spec")
  
}

per_list_endo = vector("list",11)
for (mpp_neig in c(1,2,3,5,6,8,9,10,11)){
  umap_loc_env = umap_list_2[[mpp_neig]]
  m_id = mat@cell_metadata[rownames(umap_loc_env),"mouse_id"]
  all_mouse_id = names(table(m_id))[table(m_id) > 0]
  b_id = mat@cell_metadata[rownames(umap_loc_env),"batch_set_id"]
  all_batch_id = names(table(b_id))[table(b_id) > 0]
  km_loc_5 = km_list_2[[mpp_neig]]
  
  cell_clusts = km_loc_5$cluster
  clust_cols = table(mc@colors[mc@mc[names(cell_clusts)]], cell_clusts)
  clust_cols_chosen = rownames(clust_cols)[apply(clust_cols/rowSums(clust_cols), 2, which.max)]
  names(clust_cols_chosen) = colnames(clust_cols)
  
  rel_cells = rownames(umap_loc_env)
  rel_md = md[names(mc@mc)[names(mc@mc) %in% rel_cells],]
  rel_md = rel_md[rel_md$condition != "na",]
  rel_cells = rownames(rel_md)
  rel_mice = names(table(rel_md$mouse_id)[table(rel_md$mouse_id)>0])
  m_ct = table(rel_md$mouse_id, rel_md$cell_type)[rel_mice,]
  rel_mice = rel_mice[m_ct[,"spec"]==0]
  rel_md = rel_md[rel_md$mouse_id %in% rel_mice,]
  rel_cells = rownames(rel_md)
  
  colum_order_by = c('location', 'days_post_transfer', 'condition', 'mouse_id')
  unique_md = unique(rel_md[, colum_order_by])
  unique_md = unique_md[order(as.integer(as.character(unique_md$days_post_transfer)), unique_md$condition, unique_md$mouse_id, method='shell'), ]
  rel_mice_ord = as.character(unique(unique_md$mouse_id))
  rel_mice = rel_mice_ord
  
  dir.create(paste0("figs/nmf_plots_mc2_burst/2d_mice_traj_new_local_mcs_comb/endo_",mpp_neig), recursive = T)
  per_list_endo[[mpp_neig]] = plot_2d_mouse_traj_local_mcs(umap_loc_env, md, paste0("figs/nmf_plots_mc2_burst/2d_mice_traj_new_local_mcs_comb/endo_",mpp_neig), col_by_ct=T,col_by_clust=F, clust_cols_chosen=clust_cols_chosen, cell_clusts=cell_clusts, ct="endo")
  
}

mpp_neig = 3
day_pt = 6
umap_loc_env = umap_list_2[[mpp_neig]]
m_id = mat@cell_metadata[rownames(umap_loc_env),"mouse_id"]
all_mouse_id = names(table(m_id))[table(m_id) > 0]
b_id = mat@cell_metadata[rownames(umap_loc_env),"batch_set_id"]
all_batch_id = names(table(b_id))[table(b_id) > 0]
km_loc_5 = km_list_2[[mpp_neig]]

cell_clusts = km_loc_5$cluster
clust_cols = table(mc@colors[mc@mc[names(cell_clusts)]], cell_clusts)
clust_cols_chosen = rownames(clust_cols)[apply(clust_cols/rowSums(clust_cols), 2, which.max)]
names(clust_cols_chosen) = colnames(clust_cols)

rel_cells = rownames(umap_loc_env)
rel_md = md[names(mc@mc)[names(mc@mc) %in% rel_cells],]
rel_md = rel_md[rel_md$condition != "na",]
rel_cells = rownames(rel_md)
rel_mice = names(table(rel_md$mouse_id)[table(rel_md$mouse_id)>0])
m_ct = table(rel_md$mouse_id, rel_md$cell_type)[rel_mice,]
rel_mice = rel_mice[m_ct[,"spec"]>0]
rel_md = rel_md[rel_md$mouse_id %in% rel_mice,]
rel_cells = rownames(rel_md)

colum_order_by = c('location', 'days_post_transfer', 'condition', 'mouse_id')
unique_md = unique(rel_md[, colum_order_by])
unique_md = unique_md[order(as.integer(as.character(unique_md$days_post_transfer)), unique_md$condition, unique_md$mouse_id, method='shell'), ]
rel_mice_ord = as.character(unique(unique_md$mouse_id))
rel_mice = rel_mice_ord
unique_md_day = unique_md[as.integer(as.character(unique_md$days_post_transfer))==day_pt,]
ctrl_mice = unique_md_day[as.character(unique_md_day$condition)=="ctrl","mouse_id"]
pd1_mice = unique_md_day[as.character(unique_md_day$condition)=="pd1","mouse_id"]
mice_41bb = unique_md_day[as.character(unique_md_day$condition)=="41bb","mouse_id"]
pd1_41bb_mice = unique_md_day[as.character(unique_md_day$condition)=="41bb_pd1","mouse_id"]

mean(per_list_spec[[mpp_neig]][as.character(ctrl_mice)])
mean(per_list_spec[[mpp_neig]][as.character(pd1_mice)])
mean(per_list_spec[[mpp_neig]][as.character(mice_41bb)])
mean(per_list_spec[[mpp_neig]][as.character(pd1_41bb_mice)])

mean(per_list_endo[[mpp_neig]][as.character(ctrl_mice)])
mean(per_list_endo[[mpp_neig]][as.character(pd1_mice)])
mean(per_list_endo[[mpp_neig]][as.character(mice_41bb)])
mean(per_list_endo[[mpp_neig]][as.character(pd1_41bb_mice)])



png(paste(scfigs_dir,"/nmf_plots_mc2_burst/2d_mice_traj_new_local_mcs_comb/", "prog_on_time_legend.png",sep=""), 14 * 100, 14 * 100, res = 200)
plot.new()
legend("topleft", legend=c("endo_ctrl", "endo_pd1", "endo_41bb", "endo_41bb_pd1", "spec_ctrl", "spec_pd1", "spec_41bb", "spec_41bb_pd1", "byst", "in_vitro"), pch=19, col=c("springgreen", "gold", "pink", "purple", "darkgreen", "steelblue4", "deeppink", "plum","#cee6b9", "steelblue3"), cex=2)
dev.off()


mouse_id = all_mouse_id[50]
plot(umap_loc_env[,1], umap_loc_env[,2], pch=19,col="grey", cex=0.3)
points(umap_loc_env[m_id == mouse_id,1], umap_loc_env[m_id == mouse_id,2], pch=19,col=clust_cols_chosen[cell_clusts[rownames(umap_loc_env)[m_id == mouse_id]]], cex=0.5)

batch_id = all_batch_id[11]
batch_id
plot(umap_loc_env[,1], umap_loc_env[,2], pch=19,col="grey", cex=0.3)
points(umap_loc_env[b_id == batch_id,1], umap_loc_env[b_id == batch_id,2], pch=19,col=clust_cols_chosen[cell_clusts[rownames(umap_loc_env)[b_id == batch_id]]], cex=0.5)


full_model_list$reg_full_w = full_model_list$reg_full_w_non_filt
# for (prog in c(1,4,5,6,7,8,9,11)){
#   full_model_list$reg_full_w[norm_global_mat[,prog] < active_thres_vec[prog],prog] = 0
# }

full_model_list$reg_full_w[(norm_global_mat[,2] < active_thres_vec[2]) & (norm_global_mat[,1] < active_thres_vec[1]) & (norm_global_mat[,9] < active_thres_vec[9]) & (norm_global_mat[,11] < active_thres_vec[11]),2] = 0

full_model_list$reg_full_w[(norm_global_mat[,1] < active_thres_vec[1]) & (norm_global_mat[,3] < active_thres_vec[3]) & (norm_global_mat[,6] < active_thres_vec[6])  & (norm_global_mat[,10] < active_thres_vec[10])  & (norm_global_mat[,5] < active_thres_vec[5]),3] = 0

full_model_list$reg_full_w[(norm_global_mat[,5] < active_thres_vec[5]) & (norm_global_mat[,10] < active_thres_vec[10]),10] = 0

n_progs = ncol(full_model_list$reg_full_w)

reg_full_w = full_model_list$reg_full_w

overlap_req = 1

# Make list of overlapping prog per env
overlap_list = vector(length = n_progs, mode="list")
for (k in 1:n_progs){
  k_mat = reg_full_w[reg_full_w[,k]>0, , drop=F]
  if (dim(k_mat)[1] > 2) { #& !k %in% c(14,17)
    overlap_list[[k]] = (1:n_progs)[colSums(k_mat > 0)>overlap_req]
  }
  else{
    overlap_list[[k]] = (1:n_progs)[colSums(k_mat > 0)>0]
  }
  
  if (length(overlap_list[[k]]) == 0 ){
    overlap_list[[k]] = k
  }
  if (length(overlap_list[[k]]) > 1){
    overlap_list[[k]] = overlap_list[[k]][overlap_list[[k]] != k]
  }
}

full_model_list$overlap_list = overlap_list
full_model_list$overlap_list[[2]] = c(1,9,11)
full_model_list$overlap_list[[3]] = c(1,6,5,10,11)
full_model_list$overlap_list[[10]] = c(5)

local_models_mcs = mcell_gene_prog_local_modeling_mcs(full_model_list, mc_id, mat_i, mc2d_id = NULL,mcs=NULL, prog_ord=NULL, nmf_dir="nmf_plots_mc2_burst/", env_dir="traj",envs_to_run = c(2,3, 10), global_add_prog=1)
dev.off()


local_models_mcs = mcell_gene_prog_local_modeling_mcs_plots(local_models_mcs, full_model_list, mc_id, mat_i, mc2d_id=NULL, mcs=NULL, prog_ord=NULL, nmf_dir="nmf_plots_mc2_burst/", env_dir="traj", plots=F, sim_plots=F, sim_progs=NULL, cd34_prog=c(2,3,10), main_prog_cols=NULL, main_gset=NULL, tf_list=tf_exp, start_order_ovlps = NULL, root_env=NULL, filter_end = c(), lead_cor_thresh=0.35, lead_col=T, tfs_vel=F, envs_to_run=c(2,3, 10), tf_cor_thresh=0.1, annot_colors=NULL, choose_lead=NULL)

mpp_neig = 2
local_model = local_models_mcs[[mpp_neig]]
local_w = local_models_mcs[[mpp_neig]]$local_w
local_h = local_models_mcs[[mpp_neig]]$local_h
norm_local_w = local_w / rowSums(local_w)
max_norm_local_w = t(t(norm_local_w) / apply(norm_local_w, 2, max))
norm_reconst_mat = norm_local_w %*% local_h

rel_mcs = rownames(local_w)
rel_cells_from_mcs = names(mc@mc)[mc@mc %in% as.integer(rel_mcs)]
# downs_mat = scm_downsamp(mat@mat[,rel_cells_from_mcs], 750)
# downs_mat = cbind(downs_mat, mat@mat[,setdiff(rel_cells_from_mcs, colnames(downs_mat))])
# downs_mat = downs_mat[,rel_cells_from_mcs]

# W_nnls_loc_cells_umis = nmf_py$non_negative_factorization(t(as.matrix(downs_mat_loc[colnames(local_h),rel_cells_from_mcs])), n_components=as.integer(nrow(local_h)), H = local_h, update_H = F, beta_loss = 'kullback-leibler', solver='mu', alpha=0.5,l1_ratio=0.9)
# W_loc_umis = W_nnls_loc_cells_umis[[1]]
# colnames(W_loc_umis) = paste0("prog_", 1:ncol(W_loc_umis))
# rownames(W_loc_umis) = rel_cells_from_mcs
# 
# reconst_umis = W_loc_umis %*% local_h
# rownames(reconst_umis) = rel_cells_from_mcs

# 
# local_models_mcs_3 = mcell_gene_prog_local_modeling_mcs(full_model_list, mc_id, mat_i, mc2d_id = NULL,mcs=NULL, prog_ord=NULL, nmf_dir="nmf_plots_mc2_by_mcs_myb/", env_dir="",envs_to_run = mpp_neig, global_add_prog=2, refine_mc_id = mc_id_2)
# dev.off()
# 
# 
# local_models_mcs_3 = mcell_gene_prog_local_modeling_mcs_plots(local_models_mcs_3, full_model_list, mc_id, mat_i, mc2d_id=NULL, mcs=NULL, prog_ord=NULL, nmf_dir="nmf_plots_mc2_by_mcs_myb/", env_dir="", plots=F, sim_plots=F, sim_progs=NULL, cd34_prog=mpp_neig, main_prog_cols=NULL, main_gset=NULL, tf_list=tf_exp, start_order_ovlps = NULL, root_env=NULL, filter_end = c(), lead_cor_thresh=0.35, lead_col=T, tfs_vel=F, envs_to_run=mpp_neig, tf_cor_thresh=0.1, annot_colors=NULL, choose_lead=NULL, refine_mc_id = mc_id_2)
# 
# local_model = local_models_mcs_3[[mpp_neig]]
# local_w = local_models_mcs_3[[mpp_neig]]$local_w
# local_h = local_models_mcs_3[[mpp_neig]]$local_h
# norm_local_w = local_w / rowSums(local_w)
# max_norm_local_w = t(t(norm_local_w) / apply(norm_local_w, 2, max))
# norm_reconst_mat = norm_local_w %*% local_h



sp_traj_ls_mcs_2 = mcell_gene_prog_shortest_path_traj_mcs(local_models_mcs, mpp_neig, full_model_list, mc_id, root_prog=3, target_progs=NULL, other_overlapps=NULL, mc2d_id=mc_id, rel_cells=NULL, sigma_k=0.5, graph_thresh_quant=0.2, max_dist_quant=1, target_usg = 0.1, n_neighbors_umap=20, min_dist_umap=1, umap_rand_state=88,nmf_dir="nmf_plots_mc2_burst/", env_dir="traj", umap_graph=F)

other_overlapps = list()
# target_progs=c(2,3,5,10)
# other_overlapps[[1]] = c(2,3:6,8:10)[-1]
# other_overlapps[[2]] = c(2,3:6,8:10)[-2]
# other_overlapps[[3]] = c(2,3:6,8:10)[-4]
# other_overlapps[[4]] = c(2,3:6,8:10)[-8]

target_progs=c(2)
other_overlapps[[1]] = c(1)
# other_overlapps[[2]] = c(1,2)
#other_overlapps[[2]] = c(1,5)


sp_traj_ls_mcs_2 = mcell_gene_prog_shortest_path_traj_mcs(local_models_mcs, mpp_neig, full_model_list, mc_id, root_prog=3, target_progs=target_progs, other_overlapps=other_overlapps, mc2d_id=mc_id, rel_cells=NULL, sigma_k=1, graph_thresh_quant=1, max_dist_quant=1, target_usg = 0.1, n_neighbors_umap=20, min_dist_umap=1, umap_rand_state=88,nmf_dir="nmf_plots_mc2_burst/", env_dir="traj", umap_graph=F, dist_mat = sp_traj_ls_mcs_2$dist_w_clean, graph_w=sp_traj_ls_mcs_2$graph_w)

mcell_gene_prog_shortest_path_traj_mcs_plots(sp_traj_ls_mcs_2, local_models_mcs, mpp_neig, full_model_list, mc_id, target_progs=target_progs, other_overlapps=other_overlapps, mc2d_id=mc_id, rel_cells=NULL, max_dist_quant=1, nmf_dir="nmf_plots_mc2_burst/", env_dir="traj", smooth_f=10, arrows_width=10, arrows_length=0.0, arrows_col=local_model$local_prog_cols[target_progs], text_2d=F, root_traj_prog_thresh=1, arrow_quant = 0.15, add_prog=1, dev=".svg")


mpp_neig = 3
local_model = local_models_mcs[[mpp_neig]]
local_w = local_models_mcs[[mpp_neig]]$local_w
local_h = local_models_mcs[[mpp_neig]]$local_h
norm_local_w = local_w / rowSums(local_w)
max_norm_local_w = t(t(norm_local_w) / apply(norm_local_w, 2, max))
norm_reconst_mat = norm_local_w %*% local_h

rel_mcs = rownames(local_w)
rel_cells_from_mcs = names(mc@mc)[mc@mc %in% as.integer(rel_mcs)]

sp_traj_ls_mcs_3 = mcell_gene_prog_shortest_path_traj_mcs(local_models_mcs, mpp_neig, full_model_list, mc_id, root_prog=6, target_progs=NULL, other_overlapps=NULL, mc2d_id=mc_id, rel_cells=NULL, sigma_k=0.5, graph_thresh_quant=0.15, max_dist_quant=0.2, target_usg = 0.1, n_neighbors_umap=20, min_dist_umap=1, umap_rand_state=88,nmf_dir="nmf_plots_mc2_burst/", env_dir="traj", umap_graph=F)

other_overlapps = list()
# target_progs=c(2,3,5,10)
# other_overlapps[[1]] = c(2,3:6,8:10)[-1]
# other_overlapps[[2]] = c(2,3:6,8:10)[-2]
# other_overlapps[[3]] = c(2,3:6,8:10)[-4]
# other_overlapps[[4]] = c(2,3:6,8:10)[-8]

target_progs=c(1,4,5)
other_overlapps[[1]] = c(3,4,5)
other_overlapps[[2]] = c(1,3,5)
other_overlapps[[3]] = c(1,3,4)
#other_overlapps[[2]] = c(1,5)


sp_traj_ls_mcs_3 = mcell_gene_prog_shortest_path_traj_mcs(local_models_mcs, mpp_neig, full_model_list, mc_id, root_prog=6, target_progs=target_progs, other_overlapps=other_overlapps, mc2d_id=mc_id, rel_cells=NULL, sigma_k=1, graph_thresh_quant=0.15, max_dist_quant=0.2, target_usg = 0.1, n_neighbors_umap=35, min_dist_umap=1.5, umap_rand_state=88,nmf_dir="nmf_plots_mc2_burst/", env_dir="traj", umap_graph=F, dist_mat = sp_traj_ls_mcs_3$dist_w_clean, graph_w=sp_traj_ls_mcs_3$graph_w)

mcell_gene_prog_shortest_path_traj_mcs_plots(sp_traj_ls_mcs_3, local_models_mcs, mpp_neig, full_model_list, mc_id, target_progs=target_progs, other_overlapps=other_overlapps, mc2d_id=mc_id, rel_cells=NULL, max_dist_quant=0.2, nmf_dir="nmf_plots_mc2_burst/", env_dir="traj", smooth_f=10, arrows_width=7, arrows_length=0.0, arrows_col=local_model$local_prog_cols[target_progs], text_2d=F, root_traj_prog_thresh=1, arrow_quant = 0.08, add_prog=1, dev=".svg")


ot1_byst = read.csv("Ot1_byst_ctrl.csv")
ot1_byst = ot1_byst[1:10,]
ot1_byst_2 = ot1_byst[1:5,]

ot1_byst$day = ot1_byst$day + 1

p_ot1 = ggboxplot(ot1_byst, x="day", y = "Tumor", fill=ot1_byst$Loc) #, add="jitter"
p_ot1 = p_ot1 + theme(legend.position="right", text = element_text(size=20)) + ylim(0, 10) + ylab("OT1/Byst Ratio") + xlab("Day") + scale_fill_manual(values = c("steelblue1", "#6a8255", "steelblue4")) + geom_hline(yintercept=1, color="red", linetype="dotted") + geom_text(aes(0.5, 1, label="ratio=1", vjust=-0.5))
ggsave(filename = paste(scfigs_dir,"/nmf_plots_mc2_burst/umap_full_progs/violin/", "ot1_byst_tumor", ".svg",sep=""), plot=p_ot1)


mpp_neig = 10
local_model = local_models_mcs[[mpp_neig]]
local_w = local_models_mcs[[mpp_neig]]$local_w
local_h = local_models_mcs[[mpp_neig]]$local_h
norm_local_w = local_w / rowSums(local_w)
max_norm_local_w = t(t(norm_local_w) / apply(norm_local_w, 2, max))
norm_reconst_mat = norm_local_w %*% local_h

rel_mcs = rownames(local_w)
rel_cells_from_mcs = names(mc@mc)[mc@mc %in% as.integer(rel_mcs)]
sp_traj_ls_mcs_10 = mcell_gene_prog_shortest_path_traj_mcs(local_models_mcs, mpp_neig, full_model_list, mc_id, root_prog=2, target_progs=NULL, other_overlapps=NULL, mc2d_id=mc_id, rel_cells=NULL, sigma_k=0.5, graph_thresh_quant=0.2, max_dist_quant=0.2, target_usg = 0.1, n_neighbors_umap=20, min_dist_umap=1, umap_rand_state=88,nmf_dir="nmf_plots_mc2_burst/", env_dir="traj", umap_graph=F)

other_overlapps = list()
# target_progs=c(2,3,5,10)
# other_overlapps[[1]] = c(2,3:6,8:10)[-1]
# other_overlapps[[2]] = c(2,3:6,8:10)[-2]
# other_overlapps[[3]] = c(2,3:6,8:10)[-4]
# other_overlapps[[4]] = c(2,3:6,8:10)[-8]

target_progs=c(1)
other_overlapps[[1]] = c()
#other_overlapps[[2]] = c(1,5)


sp_traj_ls_mcs_10 = mcell_gene_prog_shortest_path_traj_mcs(local_models_mcs, mpp_neig, full_model_list, mc_id, root_prog=2, target_progs=target_progs, other_overlapps=other_overlapps, mc2d_id=mc_id, rel_cells=NULL, sigma_k=1, graph_thresh_quant=0.2, max_dist_quant=0.2, target_usg = 0.1, n_neighbors_umap=20, min_dist_umap=1, umap_rand_state=88,nmf_dir="nmf_plots_mc2_burst/", env_dir="traj", umap_graph=F, dist_mat = sp_traj_ls_mcs_10$dist_w_clean, graph_w=sp_traj_ls_mcs_10$graph_w)

mcell_gene_prog_shortest_path_traj_mcs_plots(sp_traj_ls_mcs_10, local_models_mcs, mpp_neig, full_model_list, mc_id, target_progs=target_progs, other_overlapps=other_overlapps, mc2d_id=mc_id, rel_cells=NULL, max_dist_quant=0.2, nmf_dir="nmf_plots_mc2_burst/", env_dir="traj", smooth_f=10, arrows_width=5, arrows_length=0.15, arrows_col=local_model$local_prog_cols[target_progs], text_2d=F, root_traj_prog_thresh=1, arrow_quant = 0.1, add_prog=1)

#sp_traj_ls_mcs_2_sv = sp_traj_ls_mcs_2
#sp_traj_ls_mcs_3_sv = sp_traj_ls_mcs_2
#sp_traj_ls_mcs_10_sv = sp_traj_ls_mcs_2

#local_models[[mpp_neig]]$local_prog_cols[1] = "#cee6b9"
local_models[[mpp_neig]]$local_prog_cols[1] = "#1f78b4"


local_model = local_models[[mpp_neig]]
local_w = local_models[[mpp_neig]]$local_w
local_h = local_models[[mpp_neig]]$local_h
norm_local_w = local_w / rowSums(local_w)
max_norm_local_w = t(t(norm_local_w) / apply(norm_local_w, 2, max))
norm_reconst_mat = norm_local_w %*% local_h
norm_mat = t(mat@mat[,rownames(norm_local_w)]) / colSums(mat@mat[,rownames(norm_local_w)])

rel_cells = rownames(norm_local_w)
dist_w = as.matrix(dist(norm_local_w[rel_cells,]))
#dist_w_sparse = as(dist_w, "sparseMatrix")
sigma_k = 0.5
sim_w = exp((-dist_w^2)/sigma_k^2)

dist_w_clean = max(sim_w) - sim_w
diag(dist_w_clean) = 0

graph_thresh = quantile(dist_w_clean, 0.2)

sp_traj_ls = mcell_gene_prog_shortest_path_traj(local_models, mpp_neig, full_model_list, mc_id, root_prog=9, target_progs=NULL, other_overlapps=NULL, mc2d_id=mc_id, rel_cells=NULL, sigma_k=0.5, graph_thresh_quant=0.2, max_dist_quant=0.2, target_usg = 0.15, n_neighbors_umap=15, min_dist_umap=0.05, umap_rand_state=88,nmf_dir=full_model_list$nmf_dir, env_dir="", umap_graph=F, graph_thresh=graph_thresh, max_dist=graph_thresh, dist_mat = dist_w_clean)


# local_model$local_prog_cols[8] = "grey"
# local_models[[mpp_neig]]$local_prog_cols[8] = "grey"
# local_model$local_prog_cols[5] = "lightblue"
# local_models[[mpp_neig]]$local_prog_cols[5] = "lightblue"

other_overlapps = list()
target_progs=c(1,5,7,12)
other_overlapps[[1]] = c(1,5:8,10:12)[-1]
other_overlapps[[2]] = c(1,5:8,10:12)[-2]
other_overlapps[[3]] = c(1,5:8,10:12)[-4]
other_overlapps[[4]] = c(1,5:8,10:12)[-8]


sp_traj_ls = mcell_gene_prog_shortest_path_traj(local_models, mpp_neig, full_model_list, mc_id, root_prog=9, target_progs=target_progs, other_overlapps=other_overlapps, mc2d_id=mc_id, rel_cells=NULL, sigma_k=1, graph_thresh_quant=0.2, max_dist_quant=0.1, target_usg = 0.1, n_neighbors_umap=15, min_dist_umap=0.05, umap_rand_state=88,nmf_dir=full_model_list$nmf_dir, env_dir="", umap_graph=F, dist_mat = sp_traj_ls$dist_w_clean, graph_thresh=graph_thresh, max_dist=graph_thresh, graph_w=sp_traj_ls$graph_w)

mcell_gene_prog_shortest_path_traj_plots(sp_traj_ls, local_models, mpp_neig, full_model_list, mc_id, target_progs=target_progs, other_overlapps=other_overlapps, mc2d_id=mc_id, rel_cells=NULL, max_dist_quant=0.2, nmf_dir=full_model_list$nmf_dir, env_dir="", smooth_f=400, arrows_width=5, arrows_length=0.15, arrows_col=local_model$local_prog_cols[target_progs], text_2d=F, root_traj_prog_thresh=1, arrow_quant = 0.15, add_prog=4, max_dist=graph_thresh)

rel_md = mat@cell_metadata[rel_cells,]
rel_mice = names(table(rel_md$mouse_id))


ord_list = list()
ord_list[[1]] = sp_traj_ls_mcs$sp_princ_order_ls[[1]]
ord_list[[2]] = sp_traj_ls_mcs$sp_princ_order_ls[[2]]
ord_list[[3]] = sp_traj_ls_mcs$sp_princ_order_ls[[3]]
ord_list[[4]] = sp_traj_ls_mcs$sp_princ_order_ls[[4]]
#ord_list[[5]] = sp_traj_ls$sp_princ_order_ls[[5]]



genes_km = NULL
leader_norm_reconst_traj_cor = NULL
mc_env = mpp_neig
n_gene_clusts = 35
smooth_n = 300
smooth_n = 10
leader_norm_reconst_mat = t(norm_reconst_mat)
reconst_sds = apply(leader_norm_reconst_mat, 1, sd)
g_hits = names(tail(reconst_sds,350))
leader_norm_reconst_mat = leader_norm_reconst_mat[g_hits,]
# leader_norm_reconst_mat = leader_norm_reconst_mat - rowMeans(leader_norm_reconst_mat)
# leader_norm_reconst_mat = leader_norm_reconst_mat / reconst_sds[g_hits]
#leader_norm_reconst_cor = tgs_cor(t(leader_norm_reconst_mat))
leader_norm_reconst_mat_traj = c()
for (op in 1:length(ord_list)){
  #order_env = local_env_orders_pcs[[op]]$ord
  order_env = ord_list[[op]]
  
  traj_mat = t(apply(leader_norm_reconst_mat[,order_env], 1, mov_a, n=min(smooth_n, length(order_env)/4)))
  traj_mat[is.na(traj_mat)] = 0
  traj_mat = traj_mat[,colSums(traj_mat)!=0]
  #colnames(traj_mat) = paste(colnames(traj_mat), op, sep="_")
  leader_norm_reconst_mat_traj = cbind(leader_norm_reconst_mat_traj, traj_mat)
}
leader_norm_reconst_traj_cor = tgs_cor(t(leader_norm_reconst_mat_traj))
genes_km = tglkmeans::TGL_kmeans(leader_norm_reconst_traj_cor, n_gene_clusts, id_column=F)
names(genes_km$cluster) = g_hits

gene_clusts = genes_km$cluster
color_hm = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)
gene_cor_traj_clust = leader_norm_reconst_traj_cor[names(sort(gene_clusts)),names(sort(gene_clusts))]
png(paste(scfigs_dir,"/nmf_plots_mc2_by_mcs/gene_clusts_traj_heatmap_30_new.png",sep=""), 30 * 100, 30 * 100, res = 300)
pheatmap(gene_cor_traj_clust, cluster_rows = F, cluster_cols = F, fontsize = 2, color=color_hm)
dev.off()

top_ann_loc = data.frame(row.names=names(local_model$local_prog_col),main_top=names(local_model$local_prog_col))#annotation_col=top_ann, annotation_colors=list(main_top=main_prog_cols), 
#myBreaks <- c(seq(0, max(norm_full_h[genes_hm,]), length.out=floor(paletteLength)))
#png(paste(scfigs_dir,"/nmf_plots/_194/local_h_heatmap.png",sep=""), 12 * 100, 12 * 100, res = 200)

color_hm = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)
png(paste(scfigs_dir,"/nmf_plots_mc2_by_mcs/gene_progs_heatmap.png",sep=""), 10 * 100, 10 * 100, res = 200)
colnames(norm_local_w) = names(local_model$local_prog_cols)
pheatmap(norm_local_w, cluster_rows = T, cluster_cols = T, fontsize = 3, color=color_hm, annotation_col=top_ann_loc, annotation_row = mc_annot[rownames(norm_local_w),c("cell_type"), drop=F], annotation_colors=list(main_top=local_model$local_prog_col, cell_type=group2col))
dev.off()

plot_genes_by_order(norm_reconst_mat, norm_reconst_mat, "Tcf7", ord_list[[1]], ord2=ord_list[[2]], smooth_n=10, nmf_dir = NULL)

rel_mcs = rownames(local_w)
rel_cells_from_mcs = names(mc@mc)[mc@mc %in% as.integer(rel_mcs)]
W_nnls_loc = .fcnnls(t(local_h), as.matrix(mat@mat[colnames(local_h),rel_cells_from_mcs] ))
W_loc = t(W_nnls_loc$coef)
rownames(W_loc) = rel_cells_from_mcs
norm_W_loc = W_loc / rowSums(W_loc)
W_loc_mcs = tgs_matrix_tapply(t(W_loc), as.character(mc@mc[rel_cells_from_mcs]), sum)
W_loc_mcs = W_loc_mcs[rel_mcs,]
norm_W_loc_mcs = W_loc_mcs / rowSums(W_loc_mcs)


for (prog_comp in 1:ncol(norm_local_w)){
  png(paste(scfigs_dir,"/nmf_plots_mc2_by_mcs/cells_vs_mcs_prog", prog_comp, ".png",sep=""), 8 * 100, 8 * 100, res = 200)
  plot(norm_local_w[,prog_comp], norm_W_loc_mcs[rel_mcs,prog_comp], col=mc@colors[as.integer(rel_mcs)], pch=19)
  abline(0,1)
  dev.off()
}

p1 = c(7,2,2,2,5,7,7,7)
p2 = c(2,5,10,3,10,3,10,5)
for (prog_comp in 1:length(p1)){
  png(paste(scfigs_dir,"/nmf_plots_mc2_by_mcs/cells_vs_mcs_prog_", p1[prog_comp], "_by_", p2[prog_comp], ".png",sep=""), 12 * 100, 12 * 100, res = 200)
  plot(norm_W_loc_mcs[rel_mcs,p1[prog_comp]], norm_W_loc_mcs[rel_mcs,p2[prog_comp]], col=mc@colors[as.integer(rel_mcs)], pch=19)
  dev.off()
}

rel_mcs = rownames(local_w)
rel_cells_from_mcs = names(mc@mc)[mc@mc %in% as.integer(rel_mcs)]
downs_mat = scm_downsamp(mat@mat[,rel_cells_from_mcs], 750)
downs_mat = cbind(downs_mat, mat@mat[,setdiff(rel_cells_from_mcs, colnames(downs_mat))])
downs_mat = downs_mat[,rel_cells_from_mcs]

rel_mat = as.matrix(mat@mat[colnames(local_h),rel_cells_from_mcs])
rel_depths = colSums(mat@mat[,rel_cells_from_mcs])
rel_mat = t(t(rel_mat) / rel_depths)
quantile(rel_depths)

W_nnls_loc_downs = .fcnnls(t(local_h), rel_mat)
W_loc_downs = t(W_nnls_loc_downs$coef)

W_nnls_loc_downs = nmf_py$non_negative_factorization(t(as.matrix(rel_mat)), n_components=as.integer(ncol(local_w)), H = local_h, update_H = F, beta_loss = 'kullback-leibler', solver='mu', alpha=0.5,l1_ratio=0.9)
W_loc_downs = W_nnls_loc_downs[[1]]
rownames(W_loc_downs) = colnames(rel_mat)

W_loc_mcs_mean = tgs_matrix_tapply(t(W_loc_downs), mc@mc[rel_cells_from_mcs], mean)
W_loc_mcs_var = tgs_matrix_tapply(t(W_loc_downs), mc@mc[rel_cells_from_mcs], var)

downs_mat = scm_downsamp(mat@mat[,names(mc@mc)], 750)
rel_cells_from_mcs = colnames(downs_mat)
W_loc_mcs_mean_mc = tgs_matrix_tapply(downs_mat[colnames(local_h),rel_cells_from_mcs], mc@mc[rel_cells_from_mcs], mean)
W_loc_mcs_var_mc = tgs_matrix_tapply(downs_mat[colnames(local_h),rel_cells_from_mcs], mc@mc[rel_cells_from_mcs], var)

W_loc_mcs_mean_mc = W_loc_mcs_mean_mc[table(mc@mc[rel_cells_from_mcs]) > 40,]
W_loc_mcs_var_mc = W_loc_mcs_var_mc[table(mc@mc[rel_cells_from_mcs]) > 40,]
mean_feat_genes = colMeans(W_loc_mcs_mean_mc)
high_mean_feat_genes = names(mean_feat_genes)[mean_feat_genes > 0.1]

varmean = W_loc_mcs_var_mc/W_loc_mcs_mean_mc
varmean[!is.finite(varmean)] = 0
varmean = varmean[,high_mean_feat_genes]
mean_W_loc_mcs_varmean_mc = apply(varmean, 2, quantile,1)
tail(sort(mean_W_loc_mcs_varmean_mc), 100)

mean_W_loc_mcs_varmean_which_mc = apply(varmean, 2, which.max)
max_varmean_mcs = rownames(varmean)[mean_W_loc_mcs_varmean_which_mc]
names(max_varmean_mcs) = names(mean_W_loc_mcs_varmean_which_mc)

which_gene = "Xcl1"
which_mc = max_varmean_mcs[which_gene]
which_mc = "228"
which_mc_cells = rel_cells_from_mcs[mc@mc[rel_cells_from_mcs] == as.integer(which_mc)]
which_mc_genes = names(tail(sort(egc[genes[genes %in% rownames(downs_mat)],which_mc]), 50))
which_mc_genes = setdiff(which_mc_genes, which_gene)
which_mc_score = apply(downs_mat[which_mc_genes,which_mc_cells], 2, sum)
plot(which_mc_score, downs_mat[which_gene,which_mc_cells], ylab=which_gene, xlab="Score", main=which_mc)


W_loc_mcs_mean_mc_high_mean = W_loc_mcs_mean_mc[,high_mean_feat_genes]
W_loc_mcs_var_mc_high_mean = W_loc_mcs_var_mc[,high_mean_feat_genes]
varmean = W_loc_mcs_var_mc/W_loc_mcs_mean_mc
varmean[!is.finite(varmean)] = 0
mean_W_loc_mcs_varmean_mc = apply(varmean, 2, quantile, 0.95)
tail(sort(mean_W_loc_mcs_varmean_mc), 100)
high_varmean_genes = names(mean_W_loc_mcs_varmean_mc)[mean_W_loc_mcs_varmean_mc > 1.7]

mc_ind = rownames(W_loc_mcs_mean_mc_high_mean)[1]
varmean_downs_mat = downs_mat[high_varmean_genes,names(mc@mc)[mc@mc == as.integer(mc_ind)]]

mc1 = scdb_mc(mm_mc_endo_spec_ids[5])
downs_mat_1 = scm_downsamp(mat@mat[,names(mc1@mc)], 450)
  
rel_cells_from_mcs_1 = rel_cells_from_mcs[rel_cells_from_mcs %in% names(mc1@mc)]
W_loc_mcs_mean_mc_1 = tgs_matrix_tapply(downs_mat_1[colnames(local_h),], mc1@mc, mean)
W_loc_mcs_var_mc_1 = tgs_matrix_tapply(downs_mat_1[colnames(local_h),], mc1@mc, var)

varmean_1 = W_loc_mcs_var_mc_1/W_loc_mcs_mean_mc_1
varmean_1[!is.finite(varmean_1)] = 0
mean_W_loc_mcs_varmean_mc_1 = apply(varmean_1, 2, quantile, 0.95)
tail(sort(mean_W_loc_mcs_varmean_mc_1), 100)


one_mc_cells = rel_cells_from_mcs[mc@mc[rel_cells_from_mcs] == 1043]
one_mc_w = W_loc_downs[one_mc_cells,]
one_mc_means = colMeans(one_mc_w)
one_mc_var = apply(one_mc_w, 2, var)

reconst_downs = W_loc_downs %*% local_h
reconst_downs_gene_mean = apply(t(reconst_downs)[,one_mc_cells], 1, mean)
reconst_downs_gene_var = apply(t(reconst_downs)[,one_mc_cells], 1, var)

for (prog2plot in 1:ncol(W_loc_mcs_mean)){
  png(paste(scfigs_dir,"/nmf_plots_mc2_by_mcs/mean_varmean_prog_", prog2plot, ".png",sep=""), 12 * 100, 12 * 100, res = 200)
  active_mcs = W_loc_mcs_mean[,prog2plot] > 0
  plot(W_loc_mcs_mean[active_mcs,prog2plot], W_loc_mcs_var[active_mcs,prog2plot]/W_loc_mcs_mean[active_mcs,prog2plot], col=mc@colors[as.integer(rownames(W_loc_mcs_mean)[active_mcs])], pch=19, xlab='mean', ylab='var/mean', main=names(local_model$local_prog_col)[prog2plot])
  dev.off()

}

ccl3_varmean = varmean[,"Ccl3"]
tail(sort(ccl3_varmean))
cor_1295 = cor(t(as.matrix(downs_mat[genes,rel_cells_from_mcs[mc@mc[rel_cells_from_mcs]==1233]])))
cor_1295[is.na(cor_1295)] = 0
tail(sort(cor_1295[,"Ccl3"]),10)

emb_varmean = varmean[,"Ccl5"]
tail(sort(emb_varmean))
cor_1287 = cor(t(as.matrix(downs_mat[genes,names(mc@mc)[mc@mc==1287]])))
cor_1287[is.na(cor_1287)] = 0
tail(sort(cor_1287[,"Ccl5"]),100)

Dnmt3b_varmean = varmean[,"Dnmt3b"]
tail(sort(Dnmt3b_varmean))
cor_210 = cor(t(as.matrix(downs_mat[genes,names(mc@mc)[mc@mc==210]])))
cor_210[is.na(cor_210)] = 0
tail(sort(cor_210[,"Dnmt3b"]),30)



norm_W_loc_downs = W_loc_downs / rowSums(W_loc_downs)


dist_w_loc = as.matrix(dist(norm_W_loc))
if (length(filt_prog)>0){
  dist_w_loc = as.matrix(dist(dist_w_loc[,-filt_prog]))
}


sigma_k_loc = 0.5
sim_w_loc = exp((-dist_w_loc^2)/sigma_k_loc^2)
dist_w_loc_clean = max(sim_w_loc) - sim_w_loc
diag(dist_w_loc_clean) = 0
umap_obj = umap_py$UMAP(n_neighbors = as.integer(15), min_dist=1, spread=1, random_state = as.integer(109), metric="precomputed") #, metric="manhattan", 
umap_loc_env = umap_obj$fit_transform(dist_w_loc_clean)
rownames(umap_loc_env) = rownames(norm_W_loc)


##### Interesting genes #####
mc_env = 11
local_model = local_models[[mc_env]]
local_w = local_models[[mc_env]]$local_w
local_h = local_models[[mc_env]]$local_h
norm_local_w = local_w / rowSums(local_w)
max_norm_local_w = t(t(norm_local_w) / apply(norm_local_w, 2, max))
reconst_mat = local_w %*% local_h
norm_reconst_mat = norm_local_w %*% local_h
norm_mat = t(mat@mat) / colSums(mat@mat)
local_mat = mat@mat[colnames(local_h), rownames(local_w)]


##### UMAP testing local #####
umap_py = import("umap")
rel_cells = rownames(norm_local_w)#[!rownames(norm_local_w) %in% c("ica_bone_marrow_MantonBM4_HiSeq_8-1_156961:TGGCGCACACACAGAG:G-1267")]

umap_obj = umap_py$UMAP(n_neighbors = as.integer(10), min_dist=0.5, spread=1) #, metric="manhattan"
norm_local_w_umap = norm_local_w[rel_cells,]
umap_16 = umap_obj$fit_transform(norm_local_w_umap)
rownames(umap_16) = rownames(norm_local_w_umap)
#umap_16_sv = umap_16
#umap_16 = umap_16_sv
dist_w = as.matrix(dist(norm_local_w[rel_cells,]))
sigma_k = 1
# k_neig = 2
# dist_k = apply(dist_w, 1, function(x){sort(x)[k_neig]})
# dist_k_pair = matrix(rep(dist_k,length(dist_k)),nrow=length(dist_k), ncol=length(dist_k))
# dist_k_pair = dist_k_pair + t(dist_k_pair)
# sim_w = (1 / sqrt(2*pi*(dist_k_pair))) * exp((-0.5*dist_w^2)/dist_k_pair)
sim_w = exp((-dist_w^2)/sigma_k^2)
# dist_w_clean = sim_w + (diag(nrow(sim_w))*10)
# min_dists = apply(dist_w_clean, 1, min)
# dist_w_clean = dist_w_clean - (diag(nrow(dist_w))*10)
# min_alpha = 0.999
# for (ii in 1:nrow(dist_w_clean)){
#   min_vec = pmin(min_dists, min_dists[ii])
#   dist_w_clean[ii, -ii] = dist_w_clean[ii, -ii] - min_alpha*min_vec[-ii]
# }

dist_w_clean = max(sim_w) - sim_w
diag(dist_w_clean) = 0
graph_thresh = quantile(dist_w_clean, 0.10)#0.2397419#0.193369 0.00435#
max_dist = quantile(dist_w_clean, 0.1)
graph_weig = dist_w_clean / graph_thresh
graph_weig[graph_weig>1] = 0
# graph_edges = reshape2::melt(graph_weig)
# graph_edges = graph_edges[graph_edges$value>0,]
graph_w =  graph.adjacency(graph_weig, mode = "undirected", weighted = T)

is.connected(graph_w)

length(sp_2$vpath[[1]])
for (ii in 1:ncol(norm_local_w)){
  png(paste(scfigs_dir,"/nmf_plots/_144/", "no_arrows_umap_proj_prog_",ii,".png",sep=""), 12 * 100, 12 * 100, res = 200)
  proj_prog = ii
  colspec = get_param("mcell_mc2d_gene_shades", package = "metacell")
  min_usage = 0
  zero_sc_v = 0
  one_sc_v = 1
  two_sc_v=2
  prog_usage = norm_local_w[rel_cells,proj_prog]
  max_usage = max(prog_usage)
  x = pmin(pmax(prog_usage, min_usage), max_usage) - min_usage
  shades = colorRampPalette(colspec)(100 * (max_usage-min_usage) + 1)
  mc_cols = shades[round(100 * x) + 1]
  
  #plot(umap_16[,1], umap_16[,2], col=c_cols, pch = 19)
  plot(umap_16[,1], umap_16[,2], bg=mc_cols, pch = 21, col="black")
  dev.off()
}
md = mat@cell_metadata[rownames(norm_local_w),]
times = as.numeric(as.character(md$days_post_transfer))

t_points = sort(unique(times))
is_pd1 = md$batch_set_id == "PD1_5" | md$batch_set_id == "PD1_treatment" | md$batch_set_id == "MC38_PD1_3" | md$batch_set_id == "MC38_PD1_2" | md$batch_set_id == "S6_Panel_PD1_1" | md$batch_set_id == "pd1_8"
is_41bb = md$batch_set_id == "MC38_41BB_3" | md$batch_set_id == "MC38_41BB_2" | md$batch_set_id == "41bb_8" 
is_41bb_pd1 = md$batch_set_id == "MC38_41BB+PD1_3" | md$batch_set_id == "41bb+pd1_8" | md$batch_set_id == "MC38_41BB+PD1_2"
is_filt = md$batch_set_id %in% c('BFP_pd1', '41BB_ko_pd1', 'Delay_OT1_GFPOT1_PD1', 'Delay_only_OT1_ctrl', 'Delay_only_OT1_PD1', 'Delay_OT1_GFPOT1_ctrl', 'PD1_ko_pd1', 'PD1_ko_ctrl', 'Zbdb32_ko_pd1', 'ID3_KO', 'Zbdb32_ko_Ctrl', 'BFP_Ctrl', 'ID3_control', '41BB_ko_Ctrl')
is_ctrl = (!is_pd1) & (!is_41bb) & (!is_41bb_pd1) & (!is_filt)
is_spec = md$cell_type == "spec"
is_endo = md$cell_type == "endo"
is_byst = md$cell_type == "byst"
is_tumor = md$location == "tumor"
is_spleen= md$location == "spleen"
is_ln= md$location == "LN"

is_mc38 = md$batch_set_id %in% c('MC38_41BB+PD1_2', 'MC38_41BB_2', 'MC38_41BB+PD1_3', 'MC38_Ctrl_2', 'MC38_PD1_3', 'MC38_PD1_2', 'MC38_41BB_3', 'MC38_Ctrl_3')
is_b16_endo = is_endo & !is_mc38

col_by_mc = mc@colors[mc@mc[rownames(norm_local_w)]]
names(col_by_mc) = rownames(norm_local_w)
for (t in t_points){
  max_cells = 300
  cond_cells = rownames(norm_local_w)[times == t & is_pd1]
  cond_cells = sample(cond_cells, ifelse(length(cond_cells)>max_cells,max_cells, length(cond_cells)))
  if (length(cond_cells > 10)){
    png(paste(scfigs_dir,"/nmf_plots/_144/", "umap_pd1_time_",t,".png",sep=""), 600,600)
    par(mar=c(0,0,0,0))
    #plot(umap_16[,1], umap_16[,2], col=c_cols, pch = 19)
    plot(umap_16[,1], umap_16[,2], pch = 19, col="grey", cex=0.2)
    points(umap_16[cond_cells,1], umap_16[cond_cells,2], bg=col_by_mc[cond_cells], pch = 21, col="black")
    text(8,-5, sprintf("T=%d aPD1", t), cex=2.5)
    dev.off()
  }
  
  cond_cells = rownames(norm_local_w)[times == t & is_ctrl]
  cond_cells = sample(cond_cells, ifelse(length(cond_cells)>max_cells,max_cells, length(cond_cells)))
  if (length(cond_cells > 10)){
    png(paste(scfigs_dir,"/nmf_plots/_144/", "umap_ctrl_time_",t,".png",sep=""), 600,600)
    par(mar=c(0,0,0,0))
    #plot(umap_16[,1], umap_16[,2], col=c_cols, pch = 19)
    plot(umap_16[,1], umap_16[,2], pch = 19, col="grey", cex=0.2)
    points(umap_16[cond_cells,1], umap_16[cond_cells,2], bg=col_by_mc[cond_cells], pch = 21, col="black")
    text(8,-5, sprintf("T=%d no treat", t), cex=2.5)
    dev.off()
  }
  
  cond_cells = rownames(norm_local_w)[times == t & is_41bb]
  cond_cells = sample(cond_cells, ifelse(length(cond_cells)>max_cells,max_cells, length(cond_cells)))
  if (length(cond_cells > 10)){
    png(paste(scfigs_dir,"/nmf_plots/_144/", "umap_41bb_time_",t,".png",sep=""), 600,600)
    par(mar=c(0,0,0,0))
    #plot(umap_16[,1], umap_16[,2], col=c_cols, pch = 19)
    plot(umap_16[,1], umap_16[,2], pch = 19, col="grey", cex=0.2)
    points(umap_16[cond_cells,1], umap_16[cond_cells,2], bg=col_by_mc[cond_cells], pch = 21, col="black")
    text(8,-5, sprintf("T=%d a41BB", t), cex=2.5)
    dev.off()
  }
  
  cond_cells = rownames(norm_local_w)[times == t & is_41bb_pd1]
  cond_cells = sample(cond_cells, ifelse(length(cond_cells)>max_cells,max_cells, length(cond_cells)))
  if (length(cond_cells > 10)){
    png(paste(scfigs_dir,"/nmf_plots/_144/", "umap_41bb_pd1_time_",t,".png",sep=""), 600,600)
    par(mar=c(0,0,0,0))
    #plot(umap_16[,1], umap_16[,2], col=c_cols, pch = 19)
    plot(umap_16[,1], umap_16[,2], pch = 19, col="grey", cex=0.2)
    points(umap_16[cond_cells,1], umap_16[cond_cells,2], bg=col_by_mc[cond_cells], pch = 21, col="black")
    text(8,-5, sprintf("T=%d a41BB+aPD1", t), cex=2.5)
    dev.off()
  }
}

prog_by_t_all_pd1 = c()
prog_by_t_all_ctrl = c()
prog_by_t_all_41bb = c()
prog_by_t_all_41bb_pd1 = c()
for (t in t_points){
  
  cond_cells = rownames(norm_local_w)[times == t & is_pd1 & is_spec & is_tumor]
  if (length(cond_cells > 10)){
    prog_by_t = apply(norm_local_w[cond_cells,], 2, mean)
  }
  else{
    prog_by_t = rep(0, ncol(norm_local_w))
  }
  prog_by_t_all_pd1 = rbind(prog_by_t_all_pd1, prog_by_t)
  
  cond_cells = rownames(norm_local_w)[times == t & is_ctrl & is_spec & is_tumor]
  if (length(cond_cells > 10)){
    prog_by_t = apply(norm_local_w[cond_cells,], 2, mean)
  }
  else{
    prog_by_t = rep(0, ncol(norm_local_w))
  }
  prog_by_t_all_ctrl = rbind(prog_by_t_all_ctrl, prog_by_t)
  
  cond_cells = rownames(norm_local_w)[times == t & is_41bb & is_spec & is_tumor]
  if (length(cond_cells > 10)){
    prog_by_t = apply(norm_local_w[cond_cells,], 2, mean)
  }
  else{
    prog_by_t = rep(0, ncol(norm_local_w))
  }
  prog_by_t_all_41bb = rbind(prog_by_t_all_41bb, prog_by_t)
  
  cond_cells = rownames(norm_local_w)[times == t & is_41bb_pd1 & is_spec & is_tumor]
  if (length(cond_cells > 10)){
    prog_by_t = apply(norm_local_w[cond_cells,], 2, mean)
  }
  else{
    prog_by_t = rep(0, ncol(norm_local_w))
  }
  prog_by_t_all_41bb_pd1 = rbind(prog_by_t_all_41bb_pd1, prog_by_t)
}
for (prog in 1:ncol(norm_local_w)){
  t_s = c(2:7)
  png(paste(scfigs_dir,"/nmf_plots/_144/", "prog_on_time_",prog,"_pd1_ctrl.png",sep=""), 12 * 100, 12 * 100, res = 200)
  plot(t_points[t_s], prog_by_t_all_pd1[t_s,prog], col="darkblue", pch=19, ylim=c(0, max(c(prog_by_t_all_pd1[t_s,prog], prog_by_t_all_ctrl[t_s,prog]))), ylab=paste(prog, "_prog_norm_exp", sep=""))
  points(t_points[t_s], prog_by_t_all_ctrl[t_s,prog], col="lightblue", pch=19)
  dev.off()
}

smooth_n = 200
for (ge in genes2plt[genes2plt %in% colnames(norm_reconst_mat)]){
  gene_plt = ge
  png(paste(scfigs_dir,"/nmf_plots/_144/", "Gene_by_order_", "_", gene_plt, "_", 1, ".png",sep=""), 12 * 100, 8 * 100, res = 200)
  plot(mov_a(norm_reconst_mat[sp_order_1,gene_plt],smooth_n), xlab="kinetics", ylab=gene_plt)
  dev.off()
  png(paste(scfigs_dir,"/nmf_plots/_144/", "Gene_by_order_", "_", gene_plt, "_", 2, ".png",sep=""), 12 * 100, 8 * 100, res = 200)
  plot(mov_a(norm_reconst_mat[sp_order_2,gene_plt],smooth_n), xlab="kinetics", ylab=gene_plt)
  dev.off()
  png(paste(scfigs_dir,"/nmf_plots/_144/", "Gene_by_order_", "_", gene_plt, "_", 3, ".png",sep=""), 12 * 100, 8 * 100, res = 200)
  plot(mov_a(norm_reconst_mat[sp_order_3,gene_plt],smooth_n), xlab="kinetics", ylab=gene_plt)
  dev.off()
}

nmf_dir = "nmf_plots/"
env_dir = "/_14"
mc2d_id = mc2d_i
n_gene_clusts = 10

norm_reconst_mat = norm_local_w %*% local_h
leader_norm_reconst_mat = t(norm_reconst_mat[, colnames(reg_full_h)[reg_full_h[mc_env,]>0]])
reconst_sds = apply(leader_norm_reconst_mat, 1, sd)
g_hits = names(reconst_sds[reconst_sds > 2^-11])
leader_norm_reconst_mat = leader_norm_reconst_mat[g_hits,]
leader_norm_reconst_mat = leader_norm_reconst_mat - rowMeans(leader_norm_reconst_mat)
leader_norm_reconst_mat = leader_norm_reconst_mat / reconst_sds[g_hits]
#leader_norm_reconst_cor = tgs_cor(t(leader_norm_reconst_mat))
leader_norm_reconst_mat_traj = c()
rest_ovlp = c(11,7,1)
all_orders = list()
all_orders[[7]] = sp_order_2
all_orders[[11]] = sp_order_1
all_orders[[1]] = sp_order_3
for (op in rest_ovlp){
  #order_env = local_env_orders_pcs[[op]]$ord
  order_env = all_orders[[op]]
  my_mcell_mc2d_plot(mc2d_id, legend_pos="None", plot_single_cells=FALSE, plot_nm_add=paste("_no_overlap_clust_",mc_env,"_order_pc_",op,sep=""), clust_only=reg_full_w[,mc_env]>0, plot_edges=F, mc_cexs=6, nmf_fn = paste(nmf_dir, paste(env_dir, mc_env,sep=""), "/",sep=""), prog_usage = reg_full_w, main_prog = mc_env, overlapp_prog = overlap_list[[mc_env]], filt_prog = filt_prog, split_factor=split_elipse_f, arrow_order=order_env, arrow_quant = 0.10, clust_overlapp = apply(reg_full_w[,overlap_list[[mc_env]], drop=F]>0, 1, sum)>0, arrows_width=arrows_width)
  traj_mat = t(apply(leader_norm_reconst_mat[,order_env], 1, mov_a, n=min(smooth_n, length(order_env)/4)))
  traj_mat[is.na(traj_mat)] = 0
  traj_mat = traj_mat[,colSums(traj_mat)!=0]
  #colnames(traj_mat) = paste(colnames(traj_mat), op, sep="_")
  leader_norm_reconst_mat_traj = cbind(leader_norm_reconst_mat_traj, traj_mat)
}
leader_norm_reconst_traj_cor = tgs_cor(t(leader_norm_reconst_mat_traj))
genes_km = tglkmeans::TGL_kmeans(leader_norm_reconst_traj_cor, n_gene_clusts, id_column=F)
names(genes_km$cluster) = g_hits
for (clust_i in 1:n_gene_clusts){
  dir.create(paste(scfigs_dir,"/", nmf_dir, paste(env_dir, mc_env,sep=""),"/gene_cluster_", clust_i, sep=""), showWarnings = F)
  for (g in g_hits[genes_km$cluster == clust_i]){
    for (op in rest_ovlp){
      #order_env = local_env_orders_pcs[[op]]$ord
      order_env = all_orders[[op]]
      png(paste(scfigs_dir,"/", nmf_dir, paste(env_dir, mc_env,sep=""),"/gene_cluster_", clust_i,"/", op, "_clust_", genes_km$cluster[g],"_", g, "pc_order.png",sep=""), 16 * 100, 8 * 100, res = 100)
      plot(mov_a(norm_reconst_mat[order_env, g], n=min(smooth_n, length(order_env)/4)), xlab=paste("order cells trajectory", op), ylab=paste(g, "clust", genes_km$cluster[g], "norm reconst expression", sep=" "))
      # lo = loess(norm_reconst_mat[order_env, g]~c(1:length(order_env)))
      # lines(predict(lo), col='red', lwd=2)
      dev.off()
    }
    
  }
}

gene_clusts = genes_km$cluster
#gene_cor = tgs_cor(norm_reconst_mat[,names(sort(gene_clusts))])
gene_cor = leader_norm_reconst_traj_cor[names(sort(gene_clusts)),names(sort(gene_clusts))]

paletteLength <- 500
myColor <- colorRampPalette(c("royalblue", "white", "red"))(paletteLength)#colorRampPalette(c("blue","green", "white", "orange","red"))(paletteLength)
myBreaks <- c(seq(-1, 1, length.out=floor(paletteLength)))

png(paste(scfigs_dir,"/nmf_plots/_144/gene_clusts_heatmap.png",sep=""), 12 * 100, 12 * 100, res = 200)
pheatmap(gene_cor, cluster_rows = F, cluster_cols = F, fontsize = 3, color = myColor, breaks = myBreaks)
dev.off()  

norm_local_h = t(local_h) / colSums(local_h)
names(local_model$local_prog_cols)[1:global_add_prog] = paste(names(local_model$local_prog_cols)[1:global_add_prog], 1:global_add_prog, sep="_")
colnames(norm_local_h) = names(local_model$local_prog_cols)
gene_means_loc = rowMeans(local_mat) 
genes_hm = c()
for (dd in 1:ncol(norm_local_h)){
  genes_hm = c(genes_hm, names(tail(sort(norm_local_h[gene_means_loc>0.03,dd]), 20)))
}
genes_hm = unique(genes_hm)
for (ii in 1:ncol(norm_local_w)){
  png(paste(scfigs_dir,"/nmf_plots/_144/", "genes_norm_with_umap_proj_prog_",ii,".png",sep=""), 50 + 2 * 12 * 100, 12 * 100, res = 200)
  proj_prog = ii
  colspec = get_param("mcell_mc2d_gene_shades", package = "metacell")
  colspec = c("white", colspec[1], colspec[11])
  min_usage = 0
  zero_sc_v = 0
  one_sc_v = 1
  two_sc_v=2
  prog_usage = norm_local_w[,proj_prog]
  max_usage = max(prog_usage)
  x = pmin(pmax(prog_usage, min_usage), max_usage) - min_usage
  shades = colorRampPalette(colspec)(100 * (max_usage-min_usage) + 1)
  mc_cols = shades[round(100 * x) + 1]
  par(mfrow=c(1, 2))
  #par(mar=c(0,0,0,1))
  #plot(umap_16[,1], umap_16[,2], col=c_cols, pch = 19)
  plot(umap_16[,1], umap_16[,2], bg=mc_cols, pch = 21, col="black")
  prog_genes = sort(norm_full_h[,ii], decreasing = T)
  n_words = 20
  barplot(prog_genes[n_words:1], names.arg=names(prog_genes[n_words:1]), horiz = T,col=local_model$local_prog_cols[ii],las=2, cex.names = 0.7)
  #plot(umap_16[,1], umap_16[,2], col=c_cols, pch = 19)
  dev.off()
}
for (ii in 1:ncol(norm_local_w)){
  png(paste(scfigs_dir,"/nmf_plots/_144/", "genes_with_umap_proj_prog_",ii,".png",sep=""), 50 + 2 * 12 * 100, 12 * 100, res = 200)
  proj_prog = ii
  colspec = get_param("mcell_mc2d_gene_shades", package = "metacell")
  colspec = c("white", colspec[1], colspec[11])
  min_usage = 0
  zero_sc_v = 0
  one_sc_v = 1
  two_sc_v=2
  prog_usage = norm_local_w[,proj_prog]
  max_usage = max(prog_usage)
  x = pmin(pmax(prog_usage, min_usage), max_usage) - min_usage
  shades = colorRampPalette(colspec)(100 * (max_usage-min_usage) + 1)
  mc_cols = shades[round(100 * x) + 1]
  par(mfrow=c(1, 2))
  #par(mar=c(0,0,0,1))
  #plot(umap_16[,1], umap_16[,2], col=c_cols, pch = 19)
  plot(umap_16[,1], umap_16[,2], bg=mc_cols, pch = 21, col="black")
  prog_genes = sort(local_h[ii,], decreasing = T)
  n_words = 20
  barplot(prog_genes[n_words:1], names.arg=names(prog_genes[n_words:1]), horiz = T,col=local_model$local_prog_cols[ii],las=2, cex.names = 0.7)
  #plot(umap_16[,1], umap_16[,2], col=c_cols, pch = 19)
  dev.off()
}

paletteLength <- 500
myColor <- colorRampPalette(c("white", "steelblue"))(paletteLength)#colorRampPalette(c("blue","green", "white", "orange","red"))(paletteLength)

top_ann = data.frame(row.names=names(local_model$local_prog_cols),main_top=names(local_model$local_prog_cols))#annotation_col=top_ann, annotation_colors=list(main_top=main_prog_cols), 
myBreaks <- c(seq(0, max(norm_local_h[genes_hm,]), length.out=floor(paletteLength)))
#png(paste(scfigs_dir,"/nmf_plots/_194/local_h_heatmap.png",sep=""), 12 * 100, 12 * 100, res = 200)
pheatmap(norm_local_h[genes_hm,], cluster_rows = T, cluster_cols = T, fontsize = 4, color = myColor, breaks = myBreaks, annotation_col=top_ann, annotation_colors=list(main_top=local_model$local_prog_cols), filename=paste("figs","/nmf_plots/_144/local_norm_h_heatmap.png",sep=""), width=10, height=10, border_color = NA, treeheight_row = 0)


# Tsne of MCs based on usage, doesn't really work
source("~/src/FIt-SNE/fast_tsne.R", chdir=T)

# mean_nmf_mc_fts = fftRtsne(mean_w, perplexity=0, perplexity_list=c(3), learning_rate=nrow(mean_w)/12)
# #plot(mean_nmf_mc_fts[,1], mean_nmf_mc_fts[,2], pch=21, col=mc@colors)
# tsne_fn = paste(glob_top_dir,"mean_nmf_mc_tsne.png",sep="")
# png(tsne_fn, w=800, h=800)
# plot(mean_nmf_mc_fts[,1], mean_nmf_mc_fts[,2], pch=21, col=mc@colors)
# dev.off()

ts_exp = Rtsne::Rtsne(mean_w) #, perplexity=ncol(mat@mat)/100, eta=ncol(mat@mat)/12, num_threads=0
tsne_fn = paste(glob_top_dir,"mean_nmf_mc_tsne_oob.png",sep="")
png(tsne_fn, w=800, h=800)
plot(ts_exp$Y[,1], ts_exp$Y[,2], pch=19, col=mc@colors)
dev.off()

# Projecting on cells from MC full nmf
cells_mc_nmf = nmf_py$non_negative_factorization(t(mat@mat[genes,rel_cells]), n_components=as.integer(top_n_full), init='nndsvda', use_const_topics=F, const_topics=NULL, solver='mu',max_iter=as.integer(1000), init_const_topics=FALSE, beta_loss='kullback-leibler', alpha=0.0,l1_ratio=1.0, regularization='both', H=full_h, update_H=F)
cells_mc_w = cells_mc_nmf[[1]]
rownames(cells_mc_w)= rel_cells#colnames(mat@mat)
# mean_nmf_cells_from_mc_fts = fftRtsne(cells_mc_w, perplexity=0, perplexity_list=c(3), learning_rate=nrow(cells_mc_w)/12)
# #plot(mean_nmf_mc_fts[,1], mean_nmf_mc_fts[,2], pch=21, col=mc@colors)
# tsne_fn = paste(glob_top_dir,"mean_nmf_cells_from_mc_tsne.png",sep="")
# png(tsne_fn, w=800, h=800)
# plot(mean_nmf_cells_from_mc_fts[,1], mean_nmf_cells_from_mc_fts[,2], cex=0.2, pch=19, col=mc@colors[mc@mc[colnames(mat@mat)]])
# dev.off()

ts_exp_cells_from_exp = Rtsne::Rtsne(cells_mc_w) #, perplexity=ncol(new_mat@mat)/100, eta=ncol(new_mat@mat)/12, num_threads=0
tsne_fn = paste(glob_top_dir,"mean_nmf_cells_from_mc_tsne_oob.png",sep="")
png(tsne_fn, w=800, h=800)
plot(ts_exp_cells_from_exp$Y[,1], ts_exp_cells_from_exp$Y[,2], cex=0.5, pch=19, col=mc@colors[mc@mc[rel_cells]])
dev.off()

# Finding cells usage of topics, and comparing tsnes of mean topics to naive nmf topics
cells_mean_nmf = nmf_py$non_negative_factorization(t(mat@mat[genes,rel_cells]), n_components=as.integer(glob_topics_n), init='nndsvda', use_const_topics=F, const_topics=NULL, solver='mu',max_iter=as.integer(1000), init_const_topics=FALSE, beta_loss='kullback-leibler', alpha=0.0,l1_ratio=1.0, regularization='both', H=mean_topics, update_H=F)
cells_mean_w = cells_mean_nmf[[1]]
rownames(cells_mean_w)= rel_cells#colnames(mat@mat)

# mean_nmf_cells_fts = fftRtsne(cells_mean_w, perplexity=0, perplexity_list=c(30,300), learning_rate=nrow(cells_mean_w)/12)
# plot(mean_nmf_cells_fts[,1], mean_nmf_cells_fts[,2], cex=0.2, pch=19, col=mc@colors)
# tsne_fn = paste(glob_top_dir,"mean_nmf_cells_tsne.png",sep="")
# png(tsne_fn, w=800, h=800)
# plot(mean_nmf_cells_fts[,1], mean_nmf_cells_fts[,2], cex=0.2, pch=21, col=mc@colors[mc@mc[colnames(mat@mat)]])
# dev.off()

ts_exp = Rtsne::Rtsne(cells_mean_w) #, perplexity=ncol(new_mat@mat)/100, eta=ncol(new_mat@mat)/12, num_threads=0
tsne_fn = paste(glob_top_dir,"mean_nmf_cells_tsne_oob.png",sep="")
png(tsne_fn, w=800, h=800)
plot(ts_exp$Y[,1], ts_exp$Y[,2], cex=0.5, pch=19, col=mc@colors[mc@mc[rel_cells]])
dev.off()

# # Naive NMF with same number of topics for comparison
cells_mean_naive_nmf = nmf_py$non_negative_factorization(t(mat@mat[genes,rel_cells]), n_components=as.integer(glob_topics_n), init='nndsvd', use_const_topics=F, const_topics=NULL, solver='mu',max_iter=as.integer(1000), init_const_topics=FALSE, beta_loss='kullback-leibler', alpha=0.0,l1_ratio=1.0, regularization='both')
cells_mean_w_naive = cells_mean_naive_nmf[[1]]
rownames(cells_mean_w_naive)= rel_cells#colnames(mat@mat)
# 
# # mean_nmf_cells_fts_naive = fftRtsne(cells_mean_w_naive, perplexity=0, perplexity_list=c(30,300), learning_rate=nrow(cells_mean_w)/12)
# # #plot(mean_nmf_cells_fts[,1], mean_nmf_cells_fts[,2], cex=0.2, pch=19, col=mc@colors)
# # tsne_fn = paste(glob_top_dir,"mean_nmf_naive_cells_tsne.png",sep="")
# # png(tsne_fn, w=800, h=800)
# # plot(mean_nmf_cells_fts_naive[,1], mean_nmf_cells_fts_naive[,2], cex=0.2, pch=21, col=mc@colors[mc@mc[colnames(mat@mat)]])
# # dev.off()
# # 
ts_exp_naive = Rtsne::Rtsne(cells_mean_w_naive) #, perplexity=ncol(new_mat@mat)/100, eta=ncol(new_mat@mat)/12, num_threads=0
tsne_fn = paste(glob_top_dir,"mean_nmf_naive_cells_tsne_oob.png",sep="")
png(tsne_fn, w=800, h=800)
plot(ts_exp_naive$Y[,1], ts_exp_naive$Y[,2], cex=0.5, pch=19, col=mc@colors[mc@mc[rel_cells]])
dev.off()



#####----Create joint spec byst metacell model ----######
# Dense time series analysis
sell_genes = c("Tcf7", "Klf2", "Sell", "S1pr1", "Dapl1")
sell_mat = mat@mat[sell_genes, ]
sell_sum = colSums(sell_mat)
md = mat@cell_metadata[colnames(mat@mat),]

for (ctype in c('spec', 'byst', 'cd8_in_vitro')) {
  for (ttype in c('tumor','spleen', 'LN', 'cd8_in_vitro')){
    for (dtype in 0:6){
      sell_cells = mat@cell_metadata[mat@cells, ] %>%
        tibble::rownames_to_column("cell_id") %>%
        filter(days_post_transfer == dtype & cell_type == ctype & location == ttype)
      sell_cells$batch_set_id[sell_cells$batch_set_id == "PD1_5"] = "PD1_5"
      sell_cells$batch_set_id[sell_cells$batch_set_id == "PD1_6"] = "PD1_5"
      sell_cells$batch_set_id[sell_cells$batch_set_id == "Ctrl_5"] = "Ctrl_5"
      sell_cells$batch_set_id[sell_cells$batch_set_id == "Ctrl_6"] = "Ctrl_5"
      if (nrow(sell_cells) > 0){
        print(paste("figs/sell_plots/sell", ctype, ttype, dtype, "density.png", sep="_"))
        sell_cells$sell_sum = sell_sum[sell_cells$cell_id]
        if (dtype == 2){
          ggplot(sell_cells, aes(x = sell_sum, fill = batch_set_id)) + geom_density(alpha=0.2, adjust=3)  + xlim(0,30) + ylim(0,0.8) + ggtitle(paste(ctype, ttype, dtype, "sell density", sep="_"))
        }
        else{
          ggplot(sell_cells, aes(x = sell_sum, fill = batch_set_id)) + geom_density(alpha=0.2, adjust=1.5)  + xlim(0,30) + ylim(0,0.4) + ggtitle(paste(ctype, ttype, dtype, "sell density", sep="_"))
        }
        ggsave(filename = paste("figs/sell_plots/sell", ctype, ttype, dtype, "density.png", sep="_"))
      }
      
    }
  }
}


min_max_umi=32
n_genes_to_show=40
min_enr=2
min_cells=20
ids=c(5,15)

des = c()
mod_nms = c("Atlas", "Zbtb32_KO")
for (ii in 1:(length(ids)-1)){
  for (jj in (ii+1):length(ids)){
    print(mm_mc_endo_spec_ids[ids[ii]])
    print(mm_mc_endo_spec_ids[ids[jj]])
    mc_1 = scdb_mc(mm_mc_endo_spec_ids[ids[ii]])
    col2grp_1 = get_mc_col2group(mc_1)
    grp2col_1 = get_mc_group2col(mc_1)
    mc_2 = scdb_mc(mm_mc_endo_spec_ids[ids[jj]])
    col2grp_2 = get_mc_col2group(mc_2)
    grp2col_2 = get_mc_group2col(mc_2)
    
    mc_comb = mc_1
    mc_comb@mc = c(mc_1@mc, mc_2@mc + max(mc_1@mc))
    
    mat_1 = scdb_mat(mm_endo_spec_ids[ids[ii]])
    mat_comb = cbind(mat_1@mat, scdb_mat(mm_endo_spec_ids[ids[jj]])@mat)
    
    mat_ds = scm_downsamp(mat_comb, scm_which_downsamp_n(mat_1))
    mod_nm1 = mod_nms[ii]
    mod_nm2 = mod_nms[jj]
    mod_nm = paste(mod_nm1, "vs", mod_nm2, sep="_")
    for (col1 in names(col2grp_1)){
      print(col1)
      if(col1 %in% names(col2grp_2)){
        nms1 = names(mc_1@mc[mc_1@colors[mc_1@mc]==col1])
        nms2 = names(mc_2@mc[mc_2@colors[mc_2@mc]==col1])
        nms2 = nms2[mat@cell_metadata[nms2,"batch_set_id"] %in% c("Zbdb32_ko_Ctrl", 'Zbdb32_ko_pd1')]
        if (length(nms1) > 0 & length(nms2) > 0){
          de = diff_expr(mc_comb, mat_ds, mcs1=NULL, mcs2=NULL, min_max_umi=min_max_umi, geo_mean=T, nms1=nms1, nms2=nms2, reg=1)
          nm = sprintf("%s_%s_vs_%s_%s_de", mod_nm1, col2grp_1[col1], mod_nm2, col2grp_2[col1])
          des[[nm]] = de
          pd1_inner_diff_expr_plot_func(de, min_enr, n_genes_to_show, mod_nm, nm, xlab=paste(mm_mc_endo_spec_ids[ids[ii]],col2grp_1[col1]), ylab=paste(mm_mc_endo_spec_ids[ids[jj]],col2grp_2[col1]), min_max_umi = min_max_umi) 
        }
        
      }
    }
    
  }
}

pd1_inner_diff_expr_plot_func = function(de, min_enr, n_genes_to_show, mod_nm, ofn, xlab, ylab, min_max_umi = 16) 
{
  png(scfigs_fn(mod_nm, ofn, scfigs_dir(mod_nm, 'diff_expr')), 800, 800)
  par(mar=c(4,4,8,8))
  x1 = log2(de$tot1 + 1)
  x2 = log2(de$tot2 + 1)
  names(x1) = names(x2) = de$gene
  ind = abs(de$enr) > min_enr
  lim = range(c(x1, x2))
  plot(x1, x2, ylab=sprintf('%s, %d cells (log2)', ylab, nrow(x2)), xlab=sprintf('%s, %d cells (log2)', xlab, nrow(x1)), xlim=lim, ylim=lim, pch=19, cex=0.5, col=ifelse(ind, 'red', 'darkgrey'), main=ofn)
  abline(a=0, b=1, col='black')
  
  enr_g = de %>% filter(enr > min_enr & (tot1 >= min_max_umi | tot2 >= min_max_umi) & (tot1 >= 0 & tot2 >= 0)) %>% head(n_genes_to_show) %>% arrange(tot2)
  dpl_g = de %>% filter(enr < -min_enr & (tot1 >= min_max_umi | tot2 >= min_max_umi) & (tot1 >= 0 & tot2 >= 0)) %>% tail(n_genes_to_show) %>% arrange(tot1)
  lab_off = (lim[2] - lim[1]) / 10 
  if (nrow(enr_g) > 0) {
    mtext(enr_g$gene, side=4, line=0.2, at=seq(lim[1] + lab_off, lim[2] - lab_off, length=nrow(enr_g)), las=2, cex=0.9)
    segments(x0=lim[2], y0=seq(lim[1] + lab_off, lim[2] - lab_off, length=nrow(enr_g)), x1=x1[enr_g$gene], y1=x2[enr_g$gene], col='grey')
    
  }
  if (nrow(dpl_g) > 0) {
    text(seq(lim[1] + lab_off, lim[2] - lab_off, length=nrow(dpl_g)), par("usr")[4] + par("csi"), offset=0, dpl_g$gene, srt = 90, xpd = TRUE, pos = 4, cex=0.9)
    
    segments(x0=seq(lim[1] + lab_off, lim[2] - lab_off, length=nrow(dpl_g)), y0=lim[2], x1=x1[dpl_g$gene], y1=x2[dpl_g$gene], col='grey')
  }
  dev.off()
}

ids = 5
des = c()
for (ii in 1:(length(ids)-1)){
  
  mc_1 = scdb_mc(mm_mc_endo_spec_ids[ids[ii]])
  col2grp_1 = get_mc_col2group(mc_1)
  grp2col_1 = get_mc_group2col(mc_1)
  
  mc_comb = mc_1
  
  mat_1 = scdb_mat(mm_endo_spec_ids[ids[ii]])
  mat_comb = mat_1@mat#cbind(mat_1@mat, scdb_mat(mm_endo_spec_ids[ids[jj]])@mat)
  
  mat_ds = scm_downsamp(mat_comb, scm_which_downsamp_n(mat_1))
  
  mod_nm = "Combined"#paste(mod_nm1, "vs", mod_nm2, sep="_")
  for (ind1 in 5:length(col2grp_1)){
    for (ind2 in ind1:length(col2grp_1)){
      col1 = names(col2grp_1)[ind1]
      col2 = names(col2grp_1)[ind2]
      if (col1 != col2){
        nms1 = names(mc_1@mc[mc_1@colors[mc_1@mc]==col1])
        nms2 = names(mc_1@mc[mc_1@colors[mc_1@mc]==col2])
        de = diff_expr(mc_comb, mat_ds, mcs1=NULL, mcs2=NULL, min_max_umi=min_max_umi, geo_mean=T, nms1=nms1, nms2=nms2, reg=1)
        nm = sprintf("%s_%s_vs_%s_de", mod_nm, col2grp_1[col1], col2grp_1[col2])
        des[[nm]] = de
        pd1_inner_diff_expr_plot_func(de, min_enr, n_genes_to_show, mod_nm, nm, xlab=paste(mm_mc_endo_spec_ids[ids[ii]],col2grp_1[col1]), ylab=paste(mm_mc_endo_spec_ids[ids[ii]],col2grp_1[col2]), min_max_umi = min_max_umi) 
      }
      
    }
  }
  
  
}


min_max_umi=16
n_genes_to_show=40
min_enr=2
min_cells=20
ids=c(8,9)

des = c()
mod_nms = c("PD1_KO", "Ctrl_KO")
for (ii in 1:(length(ids)-1)){
  for (jj in (ii+1):length(ids)){
    mc_1 = scdb_mc(mm_mc_endo_spec_ids[ids[ii]])
    col2grp_1 = get_mc_col2group(mc_1)
    grp2col_1 = get_mc_group2col(mc_1)
    mc_2 = scdb_mc(mm_mc_endo_spec_ids[ids[jj]])
    col2grp_2 = get_mc_col2group(mc_2)
    grp2col_2 = get_mc_group2col(mc_2)
    
    mc_comb = mc_1
    mc_comb@mc = c(mc_1@mc, mc_2@mc + max(mc_1@mc))
    
    mat_1 = scdb_mat(mm_endo_spec_ids[ids[ii]])
    mat_comb = cbind(mat_1@mat, scdb_mat(mm_endo_spec_ids[ids[jj]])@mat)
    
    mat_ds = scm_downsamp(mat_comb, scm_which_downsamp_n(mat_1))
    mod_nm1 = mod_nms[ii]
    mod_nm2 = mod_nms[jj]
    mod_nm = paste(mod_nm1, "vs", mod_nm2, sep="_")
    nms1 = names(mc_1@mc[mc_1@colors[mc_1@mc]==col1])
    nms2 = names(mc_2@mc[mc_2@colors[mc_2@mc]==col2])
    
    for (col1 in names(col2grp_1)){
      for (col2 in names(col2grp_2)){
        nms1 = names(mc_1@mc[mc_1@colors[mc_1@mc]==col1])
        nms2 = names(mc_2@mc[mc_2@colors[mc_2@mc]==col2])
        de = diff_expr(mc_comb, mat_ds, mcs1=NULL, mcs2=NULL, min_max_umi=min_max_umi, geo_mean=T, nms1=nms1, nms2=nms2, reg=1)
        nm = sprintf("%s_%s_vs_%s_%s_de", mod_nm1, col2grp_1[col1], mod_nm2, col2grp_2[col2])
        des[[nm]] = de
        pd1_inner_diff_expr_plot_func(de, min_enr, n_genes_to_show, mod_nm, nm, xlab=paste(mm_mc_endo_spec_ids[ids[ii]],col2grp_1[col1]), ylab=paste(mm_mc_endo_spec_ids[ids[jj]],col2grp_2[col2]), min_max_umi = min_max_umi) 
      }
    }
    
  }
}


min_max_umi=16
n_genes_to_show=40
min_enr=2
min_cells=20
id=c(5)

des = c()
md = mat@cell_metadata[names(mc@mc),]
cells_pd1 = rownames(md[md$batch_set_id=="PD1_5" & md$cell_type=="spec" & md$location=="tumor",])
cells_ctrl = rownames(md[md$batch_set_id=="Ctrl_5" & md$cell_type=="spec" & md$location=="tumor",])
mat_ds = scm_downsamp(mat@mat, scm_which_downsamp_n(mat))
for (ii in 1:length(col2group)){
  
  mod_nm = paste("Tumor_spec_PD1", "vs", "Ctrl", sep="_")
  col = names(col2group)[ii]
  if (col != "brown"){
    nms1 = cells_pd1[mc@colors[mc@mc[cells_pd1]]==col]
    nms2 = cells_ctrl[mc@colors[mc@mc[cells_pd1]]==col]
    
    de = diff_expr(mc, mat_ds, mcs1=NULL, mcs2=NULL, min_max_umi=min_max_umi, geo_mean=T, nms1=nms1, nms2=nms2, reg=1)
    nm = sprintf("%s_%s_de", mod_nm, col2group[col])
    des[[nm]] = de
    pd1_inner_diff_expr_plot_func(de, min_enr, n_genes_to_show, mod_nm, nm, xlab=paste(mm_mc_endo_spec_ids[id],col2group[col], "PD1"), ylab=paste(mm_mc_endo_spec_ids[id],col2group[col], "Ctrl"), min_max_umi = min_max_umi) 
  }
  
  
  
  
}

mc2d_comp_weighted_mgraph = function(mc_id, graph_id, ignore_mismatch=F, symmetrize=F)
{
  mc2d_K = get_param("mcell_mc2d_K", package="metacell")
  mc2d_T_edge = get_param("mcell_mc2d_T_edge", package="metacell")
  mc2d_max_confu_deg = get_param("mcell_mc2d_max_confu_deg", package="metacell")
  mc2d_edge_asym = get_param("mcell_mc2d_edge_asym", package="metacell")
  mc2d_k_expand_inout_factor = get_param("mcell_mc2d_expand_inout_factor", package="metacell")
  
  if(is.null(mc2d_max_confu_deg)) {
    stop("MC-ERR:max_confu_deg must be defined - currently both are null")
  }
  
  graph = scdb_cgraph(graph_id)
  if(is.null(graph)) {
    stop("MC-ERR: cell graph id ", graph_id, " is missing when mc 2d projection")
  }
  mc = scdb_mc(mc_id)
  if(is.null(mc)) {
    stop("MC-ERR: mc id ", mc_id, " is missing when running add_mc_from_graph")
  }
  restrict_in_degree = T
  
  if(!is.null(mc2d_max_confu_deg)) {
    
    message("comp mc graph using the graph ", graph_id, " and K ", mc2d_K)
    confu = mcell_mc_confusion_mat(mc_id, graph_id, mc2d_K, 
                                   ignore_mismatch=ignore_mismatch)
    if(symmetrize) {
      confu =confu + t(confu)
    }
    # k_expand_inout_factor=k_expand_inout_factor
    
    csize = as.matrix(table(mc@mc))
    csize = pmax(csize, 20)
    csize2 = csize %*% t(csize)
    csize2 = csize2 / median(csize)**2
    confu = confu / csize2
    
    confu_p_from = confu/rowSums(confu)
    confu_p_to = t(confu)/colSums(confu)
    
    
    rank_fr = t(apply(confu_p_from, 1, rank))
    rank_to = t(apply(confu_p_to, 1, rank))
    rank2 = rank_fr * rank_to
    diag(rank2) = 1e+6
    amgraph = apply(rank2, 1, function(x) {  rank(-x) <= (1+mc2d_max_confu_deg) })
    mgraph = amgraph * ((confu_p_from + confu_p_to) > mc2d_T_edge)
    if(restrict_in_degree) {
      amgraph2 = t(apply(rank2, 2, function(x) {  rank(-x) <= (1+mc2d_max_confu_deg) }))
      mgraph = mgraph * amgraph2
    }
    
    mgraph = mgraph>0 | t(mgraph>0)
    
  }
  weighted_mgraph = (mgraph * rank2) * 1e-6
  
  return(weighted_mgraph)
}


##### MC38 scripting #############
table(col2group[mc@colors])
col2group
cd8_mcs = (1:max(mc@mc))[mc@colors %in% c("magenta", "blue", "#1f78b4", "pink", "#5fb8f4")]
cd8_cells = names(mc@mc)[mc@mc %in% cd8_mcs]
length(cd8_cells)

ctype ='endo'
ttype = 'tumor'
pd1_bsets = c('MC38_Ctrl_1')


mel_basic_mat2mc(mm_filt_id, pd1_bsets, mm_lateral_gset_id,
                 T_top3 = 2, T_lfc=3, cgraph_knn=80, min_mc_size = 40, mc_K = 40, T_vm=0.08, T_tot=50, name=paste(ctype, pd1_bsets, "cd8",sep="_", collapse = "_"), cells=cd8_cells) #"spleen_LN"
rl()


##### MC38 Graph projection scripting ########
i = 5
mat = scdb_mat(mm_filt_id)
id = mm_endo_spec_ids[i]
mc_id = mm_mc_endo_spec_ids[i]

mc = scdb_mc(mc_id)
lfp = log2(mc@mc_fp)
egc = mc@e_gc
col2group = get_mc_col2group(mc)
col2grp = col2group
group2col = get_mc_group2col(mc)
mc2d = scdb_mc2d(mc_id)
md = mat@cell_metadata[colnames(mat@mat),]
gset = scdb_gset(id)
ref_cls_color = mc@colors[mc@mc]
names(ref_cls_color) = names(mc@mc)
atlas_cells = names(mc@mc)

is_query_mc38 = md$batch_set_id %in% c("MC38_Ctrl_2", "MC38_PD1_2", "MC38_41BB_2", "MC38_41BB+PD1_2")
query_cells = rownames(md)[is_query_mc38]
query_cells = query_cells[query_cells %in% cells$cell_id]
mutual_cells = c(atlas_cells, query_cells)

mutual_graph_id = "mmdysf_atlas_mc38_filtered_proj"
mcell_mat_ignore_cells(new_mat_id = mutual_graph_id, mat_id = mm_filt_id, ig_cells = unique(c(mat@ignore_cells, colnames(mat@mat)[!colnames(mat@mat) %in% mutual_cells])))

mcell_gset_split_by_dsmat(id, mutual_graph_id, 30)
mcell_plot_gset_cor_mats(id, mutual_graph_id)

mut_mat = scdb_mat(mutual_graph_id)

# Generate mutual cgraph for atlas and query
cgraph_knn = 100
cgraph_downsamp = T
mcell_add_cgraph_from_mat_bknn(mat_id = mutual_graph_id, gset_id = id, graph_id = mutual_graph_id, K=cgraph_knn, dsamp=cgraph_downsamp)

mut_cgraph = scdb_cgraph(mutual_graph_id)

cmp_annot = color_query_by_reference(mutual_graph_id, query_cells, ref_cls_color, 
                                     max_knn_cgraph = 50,
                                     knn_color_space = 1, ref_fraction_thresh = 0.5,
                                     color_query_cls_with_few_ref_neighbors = T,
                                     fig_dir = NULL)

is_query_mc38_ctrl = md$batch_set_id_saved %in% c("MC38_Ctrl_2")
is_query_mc38_pd1 = md$batch_set_id_saved %in% c("MC38_PD1_2")
is_query_mc38_41bb = md$batch_set_id_saved %in% c("MC38_41BB_2")
is_query_mc38_41bb_pd1 = md$batch_set_id_saved %in% c("MC38_41BB+PD1_2")
query_cells_ctrl = rownames(md)[is_query_mc38_ctrl]
query_cells_pd1 = rownames(md)[is_query_mc38_pd1]
query_cells_41bb = rownames(md)[is_query_mc38_41bb]
query_cells_41bb_pd1 = rownames(md)[is_query_mc38_41bb_pd1]

annot_frac_ctrl = table(cmp_annot$query_cls_col[query_cells_ctrl[!query_cells_ctrl %in% cmp_annot$cls_with_low_fraction_of_ref_cell_neighbors]]) / sum(table(cmp_annot$query_cls_col[query_cells_ctrl[!query_cells_ctrl %in% cmp_annot$cls_with_low_fraction_of_ref_cell_neighbors]]))
annot_frac_pd1 = table(cmp_annot$query_cls_col[query_cells_pd1[!query_cells_pd1 %in% cmp_annot$cls_with_low_fraction_of_ref_cell_neighbors]]) / sum(table(cmp_annot$query_cls_col[query_cells_pd1[!query_cells_pd1 %in% cmp_annot$cls_with_low_fraction_of_ref_cell_neighbors]]))
annot_frac_41bb = table(cmp_annot$query_cls_col[query_cells_41bb[!query_cells_41bb %in% cmp_annot$cls_with_low_fraction_of_ref_cell_neighbors]]) / sum(table(cmp_annot$query_cls_col[query_cells_41bb[!query_cells_41bb %in% cmp_annot$cls_with_low_fraction_of_ref_cell_neighbors]]))
annot_frac_41bb_pd1 = table(cmp_annot$query_cls_col[query_cells_41bb_pd1[!query_cells_41bb_pd1 %in% cmp_annot$cls_with_low_fraction_of_ref_cell_neighbors]]) / sum(table(cmp_annot$query_cls_col[query_cells_41bb_pd1[!query_cells_41bb_pd1 %in% cmp_annot$cls_with_low_fraction_of_ref_cell_neighbors]]))

annot_frac_41bb = c(annot_frac_41bb, 0)
names(annot_frac_41bb)[11] = "deeppink"
annot_frac_41bb = annot_frac_41bb[names(annot_frac_pd1)]

annot_frac = rbind(annot_frac_ctrl, annot_frac_pd1, annot_frac_41bb, annot_frac_41bb_pd1)

fn = paste("figs/", mutual_graph_id,".color_query_by_reference/proj_color_treat_vs_ctrl.png",sep="")
png(fn, w=800, h=500)
barplot(t(annot_frac), col=colnames(annot_frac), horiz=T, las=2)
dev.off()

annot_by_mouse = table(paste(md[names(cmp_annot$query_cls_col),"days_post_transfer"], md[names(cmp_annot$query_cls_col),"batch_set_id"], md[names(cmp_annot$query_cls_col),"mouse_id"]),cmp_annot$query_cls_col)
annot_by_mouse_n = annot_by_mouse / rowSums(annot_by_mouse)
fn = paste("figs/", mutual_graph_id,".color_query_by_reference/col_by_n_id.png",sep="")
png(fn, w=1200, h=500)
barplot(t(annot_by_mouse), col=colnames(annot_by_mouse), horiz=T, las=2)
dev.off()
fn = paste("figs/", mutual_graph_id,".color_query_by_reference/col_by_n_id_norm.png",sep="")
png(fn, w=1200, h=500)
barplot(t(annot_by_mouse_n), col=colnames(annot_by_mouse_n), horiz=T, las=2)
dev.off()

nbars=23
.plot_start(scfigs_fn(mutual_graph_id, "color_query_by_reference/col_by_n_id_norm"), 520, 30 * nbars + 100)
layout(matrix(1:3, nrow=1), widths=c(400, 40, 40))
par(mar=c(4, 15, 4, 0.5))

barplot(t(annot_by_mouse_n), col=colnames(annot_by_mouse_n), las=2, horiz=T, cex.names=1.3, ylim=c(0, nbars))

par(mar=c(4,0.5, 4, 0.5))
barplot(rep(1, nrow(annot_by_mouse_n)), col=c(rep("lightgreen",12), rep("darkgreen",8)), horiz=T, ylim=c(0, nbars), xaxt='n', main='days')
barplot(rep(1, nrow(annot_by_mouse_n)), col=c(rep("purple",3),rep("springgreen3",3), rep("lightblue", 3), rep("darkblue", 3), rep("purple",1),rep("springgreen3",1), rep("lightblue", 3), rep("darkblue", 3)), horiz=T, ylim=c(0, nbars), xaxt='n', main='cond')

dev.off() 

nbars=23
.plot_start(scfigs_fn(mutual_graph_id, "color_query_by_reference/col_by_n_id_non_norm"), 520, 30 * nbars + 100)
layout(matrix(1:3, nrow=1), widths=c(400, 40, 40))
par(mar=c(4, 15, 4, 0.5))

barplot(t(annot_by_mouse), col=colnames(annot_by_mouse), las=2, horiz=T, cex.names=1.3, ylim=c(0, nbars))

par(mar=c(4,0.5, 4, 0.5))
barplot(rep(1, nrow(annot_by_mouse)), col=c(rep("lightgreen",12), rep("darkgreen",8)), horiz=T, ylim=c(0, nbars), xaxt='n', main='days')
barplot(rep(1, nrow(annot_by_mouse)), col=c(rep("purple",3),rep("springgreen3",3), rep("lightblue", 3), rep("darkblue", 3), rep("purple",1),rep("springgreen3",1), rep("lightblue", 3), rep("darkblue", 3)), horiz=T, ylim=c(0, nbars), xaxt='n', main='cond')


dev.off()

# supid color name
# 2   "brown" tumor-dend
# 4   "brown" tumor-dend
# 40  "#1f78b4" cycling-trans
# 8 "#6a8255" central-memory-like
# 6 "#e2abd6" short-lived-eff-cytokines
# 20  "magenta" dysfunctional
# 28  "#975e9f" cycling-effector-trans
# 24  "gold"  cytotoxic-effector-high
# 11  "cyan"  cytotoxic-effector-low

#### endo projection ####
i = 5
mat = scdb_mat(mm_filt_id)
id = mm_endo_spec_ids[i]
mc_id = mm_mc_endo_spec_ids[i]

mc = scdb_mc(mc_id)
lfp = log2(mc@mc_fp)
egc = mc@e_gc
col2group = get_mc_col2group(mc)
col2grp = col2group
group2col = get_mc_group2col(mc)
mc2d = scdb_mc2d(mc_id)
md = mat@cell_metadata[colnames(mat@mat),]
gset = scdb_gset(id)
ref_cls_color = mc@colors[mc@mc]
names(ref_cls_color) = names(mc@mc)
atlas_cells = names(mc@mc)

is_query_endo = md$cell_type == "endo" & md$days_post_transfer == "6" & md$batch_set_id_saved %in% c("PD1_6", "Ctrl_6")
query_cells = rownames(md)[is_query_endo]
mutual_cells = c(atlas_cells, query_cells)

mutual_graph_id = "mmdysf_atlas_endo_6_proj"
mcell_mat_ignore_cells(new_mat_id = mutual_graph_id, mat_id = mm_filt_id, ig_cells = unique(c(mat@ignore_cells, colnames(mat@mat)[!colnames(mat@mat) %in% mutual_cells])))

mcell_gset_split_by_dsmat(id, mutual_graph_id, 30)
mcell_plot_gset_cor_mats(id, mutual_graph_id)

mut_mat = scdb_mat(mutual_graph_id)

# Generate mutual cgraph for atlas and query
cgraph_knn = 100
cgraph_downsamp = T
mcell_add_cgraph_from_mat_bknn(mat_id = mutual_graph_id, gset_id = id, graph_id = mutual_graph_id, K=cgraph_knn, dsamp=cgraph_downsamp)

mut_cgraph = scdb_cgraph(mutual_graph_id)

cmp_annot = color_query_by_reference(mutual_graph_id, query_cells, ref_cls_color, 
                                     max_knn_cgraph = 50,
                                     knn_color_space = 1, ref_fraction_thresh = 0.5,
                                     color_query_cls_with_few_ref_neighbors = T,
                                     fig_dir = NULL)

is_query_endo_ctrl = md$cell_type == "endo" & md$days_post_transfer == "6" & md$batch_set_id_saved %in% c("Ctrl_6")
is_query_endo_treatment = md$cell_type == "endo" & md$days_post_transfer == "6" & md$batch_set_id_saved %in% c("PD1_6")
query_cells_ctrl = rownames(md)[is_query_endo_ctrl]
query_cells_treat = rownames(md)[is_query_endo_treatment]

annot_frac_ctrl = table(cmp_annot$query_cls_col[query_cells_ctrl[!query_cells_ctrl %in% cmp_annot$cls_with_low_fraction_of_ref_cell_neighbors]]) / sum(table(cmp_annot$query_cls_col[query_cells_ctrl[!query_cells_ctrl %in% cmp_annot$cls_with_low_fraction_of_ref_cell_neighbors]]))
annot_frac_treat = table(cmp_annot$query_cls_col[query_cells_treat[!query_cells_treat %in% cmp_annot$cls_with_low_fraction_of_ref_cell_neighbors]]) / sum(table(cmp_annot$query_cls_col[query_cells_treat[!query_cells_treat %in% cmp_annot$cls_with_low_fraction_of_ref_cell_neighbors]]))
annot_frac_ctrl = c(annot_frac_ctrl, 0)
names(annot_frac_ctrl)[11] = "deeppink"
annot_frac_ctrl = annot_frac_ctrl[names(annot_frac_treat)]

annot_frac = rbind(annot_frac_ctrl, annot_frac_treat)

fn = paste("figs/", mutual_graph_id,".color_query_by_reference/proj_color_treat_vs_ctrl.png",sep="")
png(fn, w=800, h=500)
barplot(t(annot_frac), col=colnames(annot_frac), horiz=T)
dev.off()


annot_by_mouse = table(md$mouse_id[names(cmp_annot$query_cls_col)[!names(cmp_annot$query_cls_col) %in% cmp_annot$cls_with_low_fraction_of_ref_cell_neighbors]],cmp_annot$query_cls_col[!names(cmp_annot$query_cls_col) %in% cmp_annot$cls_with_low_fraction_of_ref_cell_neighbors])
annot_by_mouse = annot_by_mouse[rowSums(annot_by_mouse)>0,]
annot_by_mouse_n = annot_by_mouse / rowSums(annot_by_mouse)
fn = paste("figs/", mutual_graph_id,".color_query_by_reference/col_by_n_id.png",sep="")
png(fn, w=800, h=500)
barplot(t(annot_by_mouse), col=colnames(annot_by_mouse), horiz=T)
dev.off()
fn = paste("figs/", mutual_graph_id,".color_query_by_reference/col_by_n_id_norm.png",sep="")
png(fn, w=800, h=500)
barplot(t(annot_by_mouse_n), col=colnames(annot_by_mouse_n), horiz=T)
dev.off()

is_ko_wt = md$cell_type == "ctrl_ko_spec"
is_ko_pd1 = md$cell_type == "pd1_ko_spec"
is_ko_ctrl = md$batch_set_id == "PD1_ko_ctrl"
is_ko_treatment = md$batch_set_id == "PD1_ko_pd1"
ko_cells = is_ko_ctrl | is_ko_treatment

is_query_ko = is_ko_ctrl | is_ko_treatment
query_cells = rownames(md)[is_query_ko]
mutual_cells = c(atlas_cells, query_cells)

mutual_graph_id = "mmdysf_atlas_pd1_ko_proj"
mcell_mat_ignore_cells(new_mat_id = mutual_graph_id, mat_id = mm_filt_id, ig_cells = unique(c(mat@ignore_cells, colnames(mat@mat)[!colnames(mat@mat) %in% mutual_cells])))

mcell_gset_split_by_dsmat(id, mutual_graph_id, 30)
mcell_plot_gset_cor_mats(id, mutual_graph_id)
mut_mat = scdb_mat(mutual_graph_id)

# Generate mutual cgraph for atlas and query
cgraph_knn = 100
cgraph_downsamp = T
mcell_add_cgraph_from_mat_bknn(mat_id = mutual_graph_id, gset_id = id, graph_id = mutual_graph_id, K=cgraph_knn, dsamp=cgraph_downsamp)

mut_cgraph = scdb_cgraph(mutual_graph_id)

cmp_annot = color_query_by_reference(mutual_graph_id, query_cells, ref_cls_color, 
                                     max_knn_cgraph = 50,
                                     knn_color_space = 1,
                                     color_query_cls_with_few_ref_neighbors = T,
                                     fig_dir = NULL)


query_cells_ko_ctrl = rownames(md)[is_ko_pd1 & is_ko_ctrl]
query_cells_ko_treat = rownames(md)[is_ko_pd1 & is_ko_treatment]
query_cells_wt_ctrl = rownames(md)[is_ko_wt & is_ko_ctrl]
query_cells_wt_treat = rownames(md)[is_ko_wt & is_ko_treatment]

annot_frac_ko_ctrl = table(cmp_annot$query_cls_col[query_cells_ko_ctrl]) / sum(table(cmp_annot$query_cls_col[query_cells_ko_ctrl]))
annot_frac_ko_treat = table(cmp_annot$query_cls_col[query_cells_ko_treat]) / sum(table(cmp_annot$query_cls_col[query_cells_ko_treat]))
annot_frac_wt_ctrl = table(cmp_annot$query_cls_col[query_cells_wt_ctrl]) / sum(table(cmp_annot$query_cls_col[query_cells_wt_ctrl]))
annot_frac_wt_treat = table(cmp_annot$query_cls_col[query_cells_wt_treat]) / sum(table(cmp_annot$query_cls_col[query_cells_wt_treat]))

annot_frac = rbind(annot_frac_ko_ctrl, annot_frac_ko_treat, annot_frac_wt_ctrl, annot_frac_wt_treat)

fn = paste("figs/", mutual_graph_id,".color_query_by_reference/proj_color_treat_vs_ctrl.png",sep="")
png(fn, w=800, h=500)
barplot(t(annot_frac), col=colnames(annot_frac), horiz=T, las=2)
dev.off()


annot_by_mouse = table(paste(md[names(cmp_annot$query_cls_col),"mouse_id"],md[names(cmp_annot$query_cls_col),"cell_type"]),cmp_annot$query_cls_col)
annot_by_mouse_n = annot_by_mouse / rowSums(annot_by_mouse)
fn = paste("figs/", mutual_graph_id,".color_query_by_reference/col_by_n_id.png",sep="")
png(fn, w=1200, h=500)
barplot(t(annot_by_mouse), col=colnames(annot_by_mouse), horiz=T, las=2)
dev.off()
fn = paste("figs/", mutual_graph_id,".color_query_by_reference/col_by_n_id_norm.png",sep="")
png(fn, w=1200, h=500)
barplot(t(annot_by_mouse_n), col=colnames(annot_by_mouse_n), horiz=T, las=2)
dev.off()

nbars=8
.plot_start(scfigs_fn(mutual_graph_id, "color_query_by_reference/col_by_n_id_norm"), 520, 30 * nbars + 100)
layout(matrix(1:3, nrow=1), widths=c(400, 40, 40))
par(mar=c(4, 15, 4, 0.5))

barplot(t(annot_by_mouse_n), col=colnames(annot_by_mouse_n), las=2, horiz=T, cex.names=1.3, ylim=c(0, nbars))

par(mar=c(4,0.5, 4, 0.5))
barplot(rep(1, nrow(annot_by_mouse_n)), col=c(rep("steelblue",4),rep("lightblue",4)), horiz=T, ylim=c(0, nbars), xaxt='n', main='cond')
barplot(rep(1, nrow(md)), col=rep(c("green" ,"red"),4), horiz=T, ylim=c(0, nbars), xaxt='n', main='cell_type')


dev.off()

#### 41bb_pd1 projection ####
i = 5
mat = scdb_mat(mm_filt_id)
id = mm_endo_spec_ids[i]
mc_id = mm_mc_endo_spec_ids[i]

mc = scdb_mc(mc_id)
lfp = log2(mc@mc_fp)
egc = mc@e_gc
col2group = get_mc_col2group(mc)
col2grp = col2group
group2col = get_mc_group2col(mc)
mc2d = scdb_mc2d(mc_id)
md = mat@cell_metadata[colnames(mat@mat),]
gset = scdb_gset(id)
full_gset = scdb_gset(mm_filt_id)
ref_cls_color = mc@colors[mc@mc]
names(ref_cls_color) = names(mc@mc)
atlas_cells = names(mc@mc)

is_query_endo = md$cell_type == "endo" & md$days_post_transfer == "3" & md$batch_set_id_saved %in% c("41bb+pd1_9", "41bb_9", "pd1_9", "Ctrl_9")
query_cells = rownames(md)[is_query_endo]
dend_gens = names(full_gset@gene_set[full_gset@gene_set == 72])
tum_genes = names(full_gset@gene_set[full_gset@gene_set == 1])
cells_depth = colSums(mat@mat)
dend_frac = colSums(mat@mat[dend_gens,query_cells]) / cells_depth[query_cells]
tum_frac = colSums(mat@mat[tum_genes,query_cells]) / cells_depth[query_cells]

dend_cells = query_cells[dend_frac > 0.004]
tum_cells = query_cells[tum_frac > 0.0003]
filt_query = unique(c(dend_cells, tum_cells))

query_cells = query_cells[!query_cells %in% filt_query]
md = md[!rownames(md) %in% filt_query,]
mutual_cells = c(atlas_cells, query_cells)

mutual_graph_id = "mmdysf_atlas_41bb_pd1_9_proj"
mcell_mat_ignore_cells(new_mat_id = mutual_graph_id, mat_id = mm_filt_id, ig_cells = unique(c(mat@ignore_cells, colnames(mat@mat)[!colnames(mat@mat) %in% mutual_cells])))

mcell_gset_split_by_dsmat(id, mutual_graph_id, 30)
mcell_plot_gset_cor_mats(id, mutual_graph_id)

mut_mat = scdb_mat(mutual_graph_id)

# Generate mutual cgraph for atlas and query
cgraph_knn = 100
cgraph_downsamp = T
mcell_add_cgraph_from_mat_bknn(mat_id = mutual_graph_id, gset_id = id, graph_id = mutual_graph_id, K=cgraph_knn, dsamp=cgraph_downsamp)

mut_cgraph = scdb_cgraph(mutual_graph_id)

cmp_annot = color_query_by_reference(mutual_graph_id, query_cells, ref_cls_color, 
                                     max_knn_cgraph = 50,
                                     knn_color_space = 1, ref_fraction_thresh = 0.5,
                                     color_query_cls_with_few_ref_neighbors = T,
                                     fig_dir = NULL)

# is_query_endo_ctrl = md$batch_set_id_saved %in% c("Ctrl_8")
# is_query_endo_treatment = md$batch_set_id_saved %in% c("41bb+pd1_8", "41bb_8", "pd1_8")
# query_cells_ctrl = rownames(md)[is_query_endo_ctrl]
# query_cells_treat = rownames(md)[is_query_endo_treatment]
# 
# annot_frac_ctrl = table(cmp_annot$query_cls_col[query_cells_ctrl[!query_cells_ctrl %in% cmp_annot$cls_with_low_fraction_of_ref_cell_neighbors]]) / sum(table(cmp_annot$query_cls_col[query_cells_ctrl[!query_cells_ctrl %in% cmp_annot$cls_with_low_fraction_of_ref_cell_neighbors]]))
# annot_frac_treat = table(cmp_annot$query_cls_col[query_cells_treat[!query_cells_treat %in% cmp_annot$cls_with_low_fraction_of_ref_cell_neighbors]]) / sum(table(cmp_annot$query_cls_col[query_cells_treat[!query_cells_treat %in% cmp_annot$cls_with_low_fraction_of_ref_cell_neighbors]]))
# #annot_frac_ctrl = c(annot_frac_ctrl, 0)
# #names(annot_frac_ctrl)[11] = "deeppink"
# annot_frac_ctrl = annot_frac_ctrl[names(annot_frac_treat)]
# 
# annot_frac = rbind(annot_frac_ctrl, annot_frac_treat)
# 
# fn = paste("figs/", mutual_graph_id,".color_query_by_reference/proj_color_treat_vs_ctrl.png",sep="")
# png(fn, w=800, h=500)
# barplot(t(annot_frac), col=colnames(annot_frac), horiz=T)
# dev.off()
# 
# 
# annot_by_mouse = table(md$mouse_id[names(cmp_annot$query_cls_col)[!names(cmp_annot$query_cls_col) %in% cmp_annot$cls_with_low_fraction_of_ref_cell_neighbors]],cmp_annot$query_cls_col[!names(cmp_annot$query_cls_col) %in% cmp_annot$cls_with_low_fraction_of_ref_cell_neighbors])
# annot_by_mouse = annot_by_mouse[rowSums(annot_by_mouse)>0,]
# annot_by_mouse_n = annot_by_mouse / rowSums(annot_by_mouse)
# fn = paste("figs/", mutual_graph_id,".color_query_by_reference/col_by_n_id.png",sep="")
# png(fn, w=800, h=500)
# barplot(t(annot_by_mouse), col=colnames(annot_by_mouse), horiz=T)
# dev.off()
# fn = paste("figs/", mutual_graph_id,".color_query_by_reference/col_by_n_id_norm.png",sep="")
# png(fn, w=800, h=500)
# barplot(t(annot_by_mouse_n), col=colnames(annot_by_mouse_n), horiz=T)
# dev.off()

#c("41bb+pd1_8", "41bb_8", "pd1_8", "Ctrl_8")
is_ctrl = md$batch_set_id == "Ctrl_9"
is_pd1 = md$batch_set_id == "pd1_9"
is_41bb = md$batch_set_id == "41bb_9"
is_41bb_pd1 = md$batch_set_id == "41bb+pd1_9"
# ko_cells = is_ko_ctrl | is_ko_treatment

# is_query_ko = is_ko_ctrl | is_ko_treatment
# query_cells = rownames(md)[is_query_ko]
# mutual_cells = c(atlas_cells, query_cells)
# 
# mutual_graph_id = "mmdysf_atlas_pd1_ko_proj"
# mcell_mat_ignore_cells(new_mat_id = mutual_graph_id, mat_id = mm_filt_id, ig_cells = unique(c(mat@ignore_cells, colnames(mat@mat)[!colnames(mat@mat) %in% mutual_cells])))
# 
# mcell_gset_split_by_dsmat(id, mutual_graph_id, 30)
# mcell_plot_gset_cor_mats(id, mutual_graph_id)
# mut_mat = scdb_mat(mutual_graph_id)
# 
# # Generate mutual cgraph for atlas and query
# cgraph_knn = 100
# cgraph_downsamp = T
# mcell_add_cgraph_from_mat_bknn(mat_id = mutual_graph_id, gset_id = id, graph_id = mutual_graph_id, K=cgraph_knn, dsamp=cgraph_downsamp)
# 
# mut_cgraph = scdb_cgraph(mutual_graph_id)
# 
# cmp_annot = color_query_by_reference(mutual_graph_id, query_cells, ref_cls_color, 
#                                      max_knn_cgraph = 50,
#                                      knn_color_space = 1,
#                                      color_query_cls_with_few_ref_neighbors = T,
#                                      fig_dir = NULL)


query_cells_ctrl = rownames(md)[is_ctrl]
query_cells_pd1 = rownames(md)[is_pd1]
query_cells_41bb = rownames(md)[is_41bb]
query_cells_41bb_pd1 = rownames(md)[is_41bb_pd1]

annot_frac_ctrl = table(cmp_annot$query_cls_col[query_cells_ctrl]) / sum(table(cmp_annot$query_cls_col[query_cells_ctrl]))
annot_frac_pd1 = table(cmp_annot$query_cls_col[query_cells_pd1]) / sum(table(cmp_annot$query_cls_col[query_cells_pd1]))
annot_frac_41bb = table(cmp_annot$query_cls_col[query_cells_41bb]) / sum(table(cmp_annot$query_cls_col[query_cells_41bb]))
annot_frac_41bb_pd1 = table(cmp_annot$query_cls_col[query_cells_41bb_pd1]) / sum(table(cmp_annot$query_cls_col[query_cells_41bb_pd1]))

annot_frac_41bb = c(annot_frac_41bb, 0)
names(annot_frac_41bb)[11] = "#5fb8f4"
annot_frac_41bb = annot_frac_41bb[names(annot_frac_ctrl)]
annot_frac = rbind(annot_frac_ctrl, annot_frac_pd1, annot_frac_41bb, annot_frac_41bb_pd1)

fn = paste("figs/", mutual_graph_id,".color_query_by_reference/proj_color_treat_vs_ctrl.png",sep="")
png(fn, w=800, h=500)
barplot(t(annot_frac), col=colnames(annot_frac), horiz=T, las=2)
dev.off()


annot_by_mouse = table(paste(md[names(cmp_annot$query_cls_col),"mouse_id"],md[names(cmp_annot$query_cls_col),"batch_set_id"]),cmp_annot$query_cls_col)
annot_by_mouse = annot_by_mouse[rev(c(9,5,7,4,6,3,8,2,1)),]
annot_by_mouse_n = annot_by_mouse / rowSums(annot_by_mouse)
fn = paste("figs/", mutual_graph_id,".color_query_by_reference/col_by_n_id.png",sep="")
png(fn, w=1200, h=500)
barplot(t(annot_by_mouse), col=colnames(annot_by_mouse), horiz=T, las=2)
dev.off()
fn = paste("figs/", mutual_graph_id,".color_query_by_reference/col_by_n_id_norm.png",sep="")
png(fn, w=1200, h=500)
barplot(t(annot_by_mouse_n), col=colnames(annot_by_mouse_n), horiz=T, las=2)
dev.off()

nbars=10
.plot_start(scfigs_fn(mutual_graph_id, "color_query_by_reference/col_by_n_id_norm"), 520, 30 * nbars + 100)
layout(matrix(1:3, nrow=1), widths=c(400, 40, 40))
par(mar=c(4, 15, 4, 0.5))

barplot(t(annot_by_mouse_n), col=colnames(annot_by_mouse_n), las=2, horiz=T, cex.names=1.3, ylim=c(0, nbars))

par(mar=c(4,0.5, 4, 0.5))
barplot(rep(1, nrow(annot_by_mouse_n)), col=c(rep("steelblue",4),rep("lightblue",4)), horiz=T, ylim=c(0, nbars), xaxt='n', main='cond')
barplot(rep(1, nrow(md)), col=rep(c("green" ,"red"),4), horiz=T, ylim=c(0, nbars), xaxt='n', main='cell_type')


dev.off()


#### 41bb_pd1 projection ####
dend_cells = names(mc@mc)[mc@mc %in% sup[[13]]$mcs]

i = 5
mat = scdb_mat(mm_filt_id)
id = mm_endo_spec_ids[i]
mc_id = mm_mc_endo_spec_ids[i]

mc = scdb_mc(mc_id)
lfp = log2(mc@mc_fp)
egc = mc@e_gc
col2group = get_mc_col2group(mc)
col2grp = col2group
group2col = get_mc_group2col(mc)
mc2d = scdb_mc2d(mc_id)
md = mat@cell_metadata[colnames(mat@mat),]
gset = scdb_gset(id)
ref_cls_color = mc@colors[mc@mc]
names(ref_cls_color) = names(mc@mc)
atlas_cells = names(mc@mc)

is_query_endo = md$cell_type == "spec" & md$days_post_transfer == "6" & md$batch_set_id_saved %in% c('Zbdb32_ko_Ctrl', 'Zbdb32_ko_pd1', 'BFP_Ctrl', 'BFP_pd1') #c('Zbdb32_ko_Ctrl', 'Zbdb32_ko_pd1', 'BFP_Ctrl', 'BFP_pd1', '41bb_ko_Ctrl', '41bb_ko_pd1')
query_cells = rownames(md)[is_query_endo]
query_cells = query_cells[!query_cells %in% dend_cells]
mutual_cells = c(atlas_cells, query_cells)

mutual_graph_id = "mmdysf_atlas_zbdb32_ko_proj"
mcell_mat_ignore_cells(new_mat_id = mutual_graph_id, mat_id = mm_filt_id, ig_cells = unique(c(mat@ignore_cells, colnames(mat@mat)[!colnames(mat@mat) %in% mutual_cells])))

mcell_gset_split_by_dsmat(id, mutual_graph_id, 30)
mcell_plot_gset_cor_mats(id, mutual_graph_id)

mut_mat = scdb_mat(mutual_graph_id)

# Generate mutual cgraph for atlas and query
cgraph_knn = 100
cgraph_downsamp = T
mcell_add_cgraph_from_mat_bknn(mat_id = mutual_graph_id, gset_id = id, graph_id = mutual_graph_id, K=cgraph_knn, dsamp=cgraph_downsamp)

mut_cgraph = scdb_cgraph(mutual_graph_id)

cmp_annot = color_query_by_reference(mutual_graph_id, query_cells, ref_cls_color, 
                                     max_knn_cgraph = 50,
                                     knn_color_space = 1, ref_fraction_thresh = 0.5,
                                     color_query_cls_with_few_ref_neighbors = T,
                                     fig_dir = NULL)

# is_query_endo_ctrl = md$batch_set_id_saved %in% c("Ctrl_8")
# is_query_endo_treatment = md$batch_set_id_saved %in% c("41bb+pd1_8", "41bb_8", "pd1_8")
# query_cells_ctrl = rownames(md)[is_query_endo_ctrl]
# query_cells_treat = rownames(md)[is_query_endo_treatment]
# 
# annot_frac_ctrl = table(cmp_annot$query_cls_col[query_cells_ctrl[!query_cells_ctrl %in% cmp_annot$cls_with_low_fraction_of_ref_cell_neighbors]]) / sum(table(cmp_annot$query_cls_col[query_cells_ctrl[!query_cells_ctrl %in% cmp_annot$cls_with_low_fraction_of_ref_cell_neighbors]]))
# annot_frac_treat = table(cmp_annot$query_cls_col[query_cells_treat[!query_cells_treat %in% cmp_annot$cls_with_low_fraction_of_ref_cell_neighbors]]) / sum(table(cmp_annot$query_cls_col[query_cells_treat[!query_cells_treat %in% cmp_annot$cls_with_low_fraction_of_ref_cell_neighbors]]))
# #annot_frac_ctrl = c(annot_frac_ctrl, 0)
# #names(annot_frac_ctrl)[11] = "deeppink"
# annot_frac_ctrl = annot_frac_ctrl[names(annot_frac_treat)]
# 
# annot_frac = rbind(annot_frac_ctrl, annot_frac_treat)
# 
# fn = paste("figs/", mutual_graph_id,".color_query_by_reference/proj_color_treat_vs_ctrl.png",sep="")
# png(fn, w=800, h=500)
# barplot(t(annot_frac), col=colnames(annot_frac), horiz=T)
# dev.off()
# 
# 
# annot_by_mouse = table(md$mouse_id[names(cmp_annot$query_cls_col)[!names(cmp_annot$query_cls_col) %in% cmp_annot$cls_with_low_fraction_of_ref_cell_neighbors]],cmp_annot$query_cls_col[!names(cmp_annot$query_cls_col) %in% cmp_annot$cls_with_low_fraction_of_ref_cell_neighbors])
# annot_by_mouse = annot_by_mouse[rowSums(annot_by_mouse)>0,]
# annot_by_mouse_n = annot_by_mouse / rowSums(annot_by_mouse)
# fn = paste("figs/", mutual_graph_id,".color_query_by_reference/col_by_n_id.png",sep="")
# png(fn, w=800, h=500)
# barplot(t(annot_by_mouse), col=colnames(annot_by_mouse), horiz=T)
# dev.off()
# fn = paste("figs/", mutual_graph_id,".color_query_by_reference/col_by_n_id_norm.png",sep="")
# png(fn, w=800, h=500)
# barplot(t(annot_by_mouse_n), col=colnames(annot_by_mouse_n), horiz=T)
# dev.off()

#c("41bb+pd1_8", "41bb_8", "pd1_8", "Ctrl_8")
is_ctrl = md$batch_set_id == "BFP_Ctrl"
is_pd1 = md$batch_set_id == "BFP_pd1"
is_41bb = md$batch_set_id == "Zbdb32_ko_Ctrl"
is_41bb_pd1 = md$batch_set_id == "Zbdb32_ko_pd1"
# ko_cells = is_ko_ctrl | is_ko_treatment

# is_query_ko = is_ko_ctrl | is_ko_treatment
# query_cells = rownames(md)[is_query_ko]
# mutual_cells = c(atlas_cells, query_cells)
# 
# mutual_graph_id = "mmdysf_atlas_pd1_ko_proj"
# mcell_mat_ignore_cells(new_mat_id = mutual_graph_id, mat_id = mm_filt_id, ig_cells = unique(c(mat@ignore_cells, colnames(mat@mat)[!colnames(mat@mat) %in% mutual_cells])))
# 
# mcell_gset_split_by_dsmat(id, mutual_graph_id, 30)
# mcell_plot_gset_cor_mats(id, mutual_graph_id)
# mut_mat = scdb_mat(mutual_graph_id)
# 
# # Generate mutual cgraph for atlas and query
# cgraph_knn = 100
# cgraph_downsamp = T
# mcell_add_cgraph_from_mat_bknn(mat_id = mutual_graph_id, gset_id = id, graph_id = mutual_graph_id, K=cgraph_knn, dsamp=cgraph_downsamp)
# 
# mut_cgraph = scdb_cgraph(mutual_graph_id)
# 
# cmp_annot = color_query_by_reference(mutual_graph_id, query_cells, ref_cls_color, 
#                                      max_knn_cgraph = 50,
#                                      knn_color_space = 1,
#                                      color_query_cls_with_few_ref_neighbors = T,
#                                      fig_dir = NULL)


query_cells_ctrl = rownames(md)[is_ctrl]
query_cells_pd1 = rownames(md)[is_pd1]
query_cells_41bb = rownames(md)[is_41bb]
query_cells_41bb_pd1 = rownames(md)[is_41bb_pd1]

annot_frac_ctrl = table(cmp_annot$query_cls_col[query_cells_ctrl]) / sum(table(cmp_annot$query_cls_col[query_cells_ctrl]))
annot_frac_pd1 = table(cmp_annot$query_cls_col[query_cells_pd1]) / sum(table(cmp_annot$query_cls_col[query_cells_pd1]))
annot_frac_41bb = table(cmp_annot$query_cls_col[query_cells_41bb]) / sum(table(cmp_annot$query_cls_col[query_cells_41bb]))
annot_frac_41bb_pd1 = table(cmp_annot$query_cls_col[query_cells_41bb_pd1]) / sum(table(cmp_annot$query_cls_col[query_cells_41bb_pd1]))

annot_frac_41bb = c(annot_frac_41bb, 0)
names(annot_frac_41bb)[11] = "#5fb8f4"
annot_frac_41bb = annot_frac_41bb[names(annot_frac_ctrl)]
annot_frac = rbind(annot_frac_ctrl, annot_frac_pd1, annot_frac_41bb, annot_frac_41bb_pd1)

fn = paste("figs/", mutual_graph_id,".color_query_by_reference/proj_color_treat_vs_ctrl.png",sep="")
png(fn, w=800, h=500)
barplot(t(annot_frac), col=colnames(annot_frac), horiz=T, las=2)
dev.off()


annot_by_mouse = table(paste(md[names(cmp_annot$query_cls_col),"mouse_id"],md[names(cmp_annot$query_cls_col),"batch_set_id"]),cmp_annot$query_cls_col)
annot_by_mouse = annot_by_mouse[c(1,5,8,2,3,4,6,9,7),]
annot_by_mouse_n = annot_by_mouse / rowSums(annot_by_mouse)
fn = paste("figs/", mutual_graph_id,".color_query_by_reference/col_by_n_id.png",sep="")
png(fn, w=1200, h=500)
barplot(t(annot_by_mouse), col=colnames(annot_by_mouse), horiz=T, las=2)
dev.off()
fn = paste("figs/", mutual_graph_id,".color_query_by_reference/col_by_n_id_norm.png",sep="")
png(fn, w=1200, h=500)
barplot(t(annot_by_mouse_n), col=colnames(annot_by_mouse_n), horiz=T, las=2)
dev.off()

nbars=10
.plot_start(scfigs_fn(mutual_graph_id, "color_query_by_reference/col_by_n_id_norm"), 520, 30 * nbars + 100)
layout(matrix(1:3, nrow=1), widths=c(400, 40, 40))
par(mar=c(4, 15, 4, 0.5))

barplot(t(annot_by_mouse_n), col=colnames(annot_by_mouse_n), las=2, horiz=T, cex.names=1.3, ylim=c(0, nbars))

par(mar=c(4,0.5, 4, 0.5))
barplot(rep(1, nrow(annot_by_mouse_n)), col=c(rep("steelblue",4),rep("lightblue",4)), horiz=T, ylim=c(0, nbars), xaxt='n', main='cond')
barplot(rep(1, nrow(md)), col=rep(c("green" ,"red"),4), horiz=T, ylim=c(0, nbars), xaxt='n', main='cell_type')


dev.off()

### cd25 markers projection ####
i = 5
mat = scdb_mat(mm_filt_id)
id = mm_endo_spec_ids[i]
mc_id = mm_mc_endo_spec_ids[i]

mc = scdb_mc(mc_id)
lfp = log2(mc@mc_fp)
egc = mc@e_gc
col2group = get_mc_col2group(mc)
col2grp = col2group
group2col = get_mc_group2col(mc)
mc2d = scdb_mc2d(mc_id)
md = mat@cell_metadata[colnames(mat@mat),]
gset = scdb_gset(id)
ref_cls_color = mc@colors[mc@mc]
names(ref_cls_color) = names(mc@mc)
atlas_cells = names(mc@mc)

is_query_endo = md$amp_batch_id %in% c("AB10120", "AB10121", "AB10122")
query_cells = rownames(md)[is_query_endo]
#query_cells = query_cells[!query_cells %in% dend_cells]
mutual_cells = c(atlas_cells, query_cells)

mutual_graph_id = "mmdysf_atlas_cd25_marks_proj"
mcell_mat_ignore_cells(new_mat_id = mutual_graph_id, mat_id = mm_filt_id, ig_cells = unique(c(mat@ignore_cells, colnames(mat@mat)[!colnames(mat@mat) %in% mutual_cells])))

mcell_gset_split_by_dsmat(id, mutual_graph_id, 30)
mcell_plot_gset_cor_mats(id, mutual_graph_id)

mut_mat = scdb_mat(mutual_graph_id)

# Generate mutual cgraph for atlas and query
cgraph_knn = 100
cgraph_downsamp = T
mcell_add_cgraph_from_mat_bknn(mat_id = mutual_graph_id, gset_id = id, graph_id = mutual_graph_id, K=cgraph_knn, dsamp=cgraph_downsamp)

mut_cgraph = scdb_cgraph(mutual_graph_id)

cmp_annot = color_query_by_reference(mutual_graph_id, query_cells, ref_cls_color, 
                                     max_knn_cgraph = 50,
                                     knn_color_space = 1, ref_fraction_thresh = 0.5,
                                     color_query_cls_with_few_ref_neighbors = T,
                                     fig_dir = NULL)

# is_query_endo_ctrl = md$batch_set_id_saved %in% c("Ctrl_8")
# is_query_endo_treatment = md$batch_set_id_saved %in% c("41bb+pd1_8", "41bb_8", "pd1_8")
# query_cells_ctrl = rownames(md)[is_query_endo_ctrl]
# query_cells_treat = rownames(md)[is_query_endo_treatment]
# 
# annot_frac_ctrl = table(cmp_annot$query_cls_col[query_cells_ctrl[!query_cells_ctrl %in% cmp_annot$cls_with_low_fraction_of_ref_cell_neighbors]]) / sum(table(cmp_annot$query_cls_col[query_cells_ctrl[!query_cells_ctrl %in% cmp_annot$cls_with_low_fraction_of_ref_cell_neighbors]]))
# annot_frac_treat = table(cmp_annot$query_cls_col[query_cells_treat[!query_cells_treat %in% cmp_annot$cls_with_low_fraction_of_ref_cell_neighbors]]) / sum(table(cmp_annot$query_cls_col[query_cells_treat[!query_cells_treat %in% cmp_annot$cls_with_low_fraction_of_ref_cell_neighbors]]))
# #annot_frac_ctrl = c(annot_frac_ctrl, 0)
# #names(annot_frac_ctrl)[11] = "deeppink"
# annot_frac_ctrl = annot_frac_ctrl[names(annot_frac_treat)]
# 
# annot_frac = rbind(annot_frac_ctrl, annot_frac_treat)
# 
# fn = paste("figs/", mutual_graph_id,".color_query_by_reference/proj_color_treat_vs_ctrl.png",sep="")
# png(fn, w=800, h=500)
# barplot(t(annot_frac), col=colnames(annot_frac), horiz=T)
# dev.off()
# 
# 
# annot_by_mouse = table(md$mouse_id[names(cmp_annot$query_cls_col)[!names(cmp_annot$query_cls_col) %in% cmp_annot$cls_with_low_fraction_of_ref_cell_neighbors]],cmp_annot$query_cls_col[!names(cmp_annot$query_cls_col) %in% cmp_annot$cls_with_low_fraction_of_ref_cell_neighbors])
# annot_by_mouse = annot_by_mouse[rowSums(annot_by_mouse)>0,]
# annot_by_mouse_n = annot_by_mouse / rowSums(annot_by_mouse)
# fn = paste("figs/", mutual_graph_id,".color_query_by_reference/col_by_n_id.png",sep="")
# png(fn, w=800, h=500)
# barplot(t(annot_by_mouse), col=colnames(annot_by_mouse), horiz=T)
# dev.off()
# fn = paste("figs/", mutual_graph_id,".color_query_by_reference/col_by_n_id_norm.png",sep="")
# png(fn, w=800, h=500)
# barplot(t(annot_by_mouse_n), col=colnames(annot_by_mouse_n), horiz=T)
# dev.off()

#c("41bb+pd1_8", "41bb_8", "pd1_8", "Ctrl_8")
is_cd25 = md$amp_batch_id == "AB10120"
is_cd25_neg = md$amp_batch_id == "AB10121"
is_cd25_all = md$amp_batch_id == "AB10122"
# ko_cells = is_ko_ctrl | is_ko_treatment

# is_query_ko = is_ko_ctrl | is_ko_treatment
# query_cells = rownames(md)[is_query_ko]
# mutual_cells = c(atlas_cells, query_cells)
# 
# mutual_graph_id = "mmdysf_atlas_pd1_ko_proj"
# mcell_mat_ignore_cells(new_mat_id = mutual_graph_id, mat_id = mm_filt_id, ig_cells = unique(c(mat@ignore_cells, colnames(mat@mat)[!colnames(mat@mat) %in% mutual_cells])))
# 
# mcell_gset_split_by_dsmat(id, mutual_graph_id, 30)
# mcell_plot_gset_cor_mats(id, mutual_graph_id)
# mut_mat = scdb_mat(mutual_graph_id)
# 
# # Generate mutual cgraph for atlas and query
# cgraph_knn = 100
# cgraph_downsamp = T
# mcell_add_cgraph_from_mat_bknn(mat_id = mutual_graph_id, gset_id = id, graph_id = mutual_graph_id, K=cgraph_knn, dsamp=cgraph_downsamp)
# 
# mut_cgraph = scdb_cgraph(mutual_graph_id)
# 
# cmp_annot = color_query_by_reference(mutual_graph_id, query_cells, ref_cls_color, 
#                                      max_knn_cgraph = 50,
#                                      knn_color_space = 1,
#                                      color_query_cls_with_few_ref_neighbors = T,
#                                      fig_dir = NULL)


query_cells_ctrl = rownames(md)[is_cd25]
query_cells_pd1 = rownames(md)[is_cd25_neg]
query_cells_41bb = rownames(md)[is_cd25_all]

annot_frac_ctrl = table(cmp_annot$query_cls_col[query_cells_ctrl]) / sum(table(cmp_annot$query_cls_col[query_cells_ctrl]))
annot_frac_pd1 = table(cmp_annot$query_cls_col[query_cells_pd1]) / sum(table(cmp_annot$query_cls_col[query_cells_pd1]))
annot_frac_41bb = table(cmp_annot$query_cls_col[query_cells_41bb]) / sum(table(cmp_annot$query_cls_col[query_cells_41bb]))
annot_frac_41bb_n = annot_frac_ctrl
annot_frac_41bb_n[names(annot_frac_41bb)] = annot_frac_41bb
annot_frac_41bb_n["#5fb8f4"] = 0
annot_frac_41bb_n["deeppink"] = 0
annot_frac_41bb = annot_frac_41bb_n
# annot_frac_41bb = c(annot_frac_41bb, 0)
# names(annot_frac_41bb)[11] = "#5fb8f4"
annot_frac_41bb = annot_frac_41bb[names(annot_frac_ctrl)]
annot_frac = rbind(annot_frac_ctrl, annot_frac_pd1, annot_frac_41bb)

fn = paste("figs/", mutual_graph_id,".color_query_by_reference/proj_color_treat_vs_ctrl.png",sep="")
png(fn, w=800, h=500)
barplot(t(annot_frac), col=colnames(annot_frac), horiz=T, las=2)
dev.off()


annot_by_mouse = table(paste(md[names(cmp_annot$query_cls_col),"amp_batch_id"],md[names(cmp_annot$query_cls_col),"batch_set_id"]),cmp_annot$query_cls_col)
#annot_by_mouse = annot_by_mouse[c(1,5,8,2,3,4,6,9,7),]
annot_by_mouse_n = annot_by_mouse / rowSums(annot_by_mouse)
fn = paste("figs/", mutual_graph_id,".color_query_by_reference/col_by_n_id.png",sep="")
png(fn, w=1200, h=500)
barplot(t(annot_by_mouse), col=colnames(annot_by_mouse), horiz=T, las=2)
dev.off()
fn = paste("figs/", mutual_graph_id,".color_query_by_reference/col_by_n_id_norm.png",sep="")
png(fn, w=1200, h=500)
barplot(t(annot_by_mouse_n), col=colnames(annot_by_mouse_n), horiz=T, las=2)
dev.off()

nbars=4
.plot_start(scfigs_fn(mutual_graph_id, "color_query_by_reference/col_by_n_id_norm"), 520, 30 * nbars + 100)
layout(matrix(1:3, nrow=1), widths=c(400, 40, 40))
par(mar=c(4, 15, 4, 0.5))

barplot(t(annot_by_mouse_n), col=colnames(annot_by_mouse_n), las=2, horiz=T, cex.names=1.3, ylim=c(0, nbars))

par(mar=c(4,0.5, 4, 0.5))
barplot(rep(1, nrow(annot_by_mouse_n)), col=c(rep("pink",1),rep("lightblue",1), "blue"), horiz=T, ylim=c(0, nbars), xaxt='n', main='cond')
#barplot(rep(1, nrow(md)), col=rep(c("green" ,"red"),4), horiz=T, ylim=c(0, nbars), xaxt='n', main='cell_type')


dev.off()
