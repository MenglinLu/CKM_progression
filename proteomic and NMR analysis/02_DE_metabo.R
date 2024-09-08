


library(dplyr)

# data
sustain_res <- readr::read_csv("./data/result_subtype_stage_allcohort.csv") %>% 
    dplyr::rename(eid=Eid, Group=ml_subtype)

sustain_res[sustain_res$Flag_disease==0,"Group"] <- 0
sample_meta <- dplyr::select(sustain_res, eid, Group)

prot_tib_tidy <- readRDS("/home/zjia/ukb2/aging_ukb/tidy_data/metabol_tib_tidy.RDS")


meta_prot <- dplyr::inner_join(sample_meta, prot_tib_tidy )

table(meta_prot$Group)
################################################################################
g3_mean <- dplyr::group_by(meta_prot, Group) %>% 
    dplyr::summarise_all(mean, na.rm=TRUE ) 
g3_mean_mat <- t(g3_mean[,-c(1:2)])
colnames(g3_mean_mat) <- c("C_average", "S1_average", "S2_average", "S3_average", "S4_average")
g3_mean_tib <- tibble::rownames_to_column(as.data.frame(g3_mean_mat), "gene")



################################################################################
# T-test
S1_DEA <- matrixTests::col_t_equalvar(subset(meta_prot[meta_prot$Group==1,], select= -c(1:2) ) , 
                                      subset(meta_prot[meta_prot$Group !=1,], select= -c(1:2) ) ) %>% 
    tibble::rownames_to_column("gene") %>% 
    dplyr::mutate(S1_fdr=p.adjust(pvalue, "fdr"), S1_logFC=log2(mean.x/mean.y) ) %>% 
    dplyr::arrange(pvalue) %>% 
    dplyr::left_join(g3_mean_tib)

S2_DEA <- matrixTests::col_t_equalvar(subset(meta_prot[meta_prot$Group==2,], select= -c(1:2) ) , 
                                      subset(meta_prot[meta_prot$Group !=2,], select= -c(1:2) ) ) %>% 
    tibble::rownames_to_column("gene") %>% 
    dplyr::mutate(S2_fdr=p.adjust(pvalue, "fdr"), S2_logFC=log2(mean.x/mean.y)  ) %>% 
    dplyr::arrange(pvalue) %>% 
    dplyr::left_join(g3_mean_tib)
S3_DEA <- matrixTests::col_t_equalvar(subset(meta_prot[meta_prot$Group==3,], select= -c(1:2) ) , 
                                      subset(meta_prot[meta_prot$Group !=3,], select= -c(1:2) ) ) %>% 
    tibble::rownames_to_column("gene") %>% 
    dplyr::mutate(S3_fdr=p.adjust(pvalue, "fdr"), S3_logFC=log2(mean.x/mean.y)  ) %>% 
    dplyr::arrange(pvalue) %>% 
    dplyr::left_join(g3_mean_tib)
S4_DEA <- matrixTests::col_t_equalvar(subset(meta_prot[meta_prot$Group==4,], select= -c(1:2) ) , 
                                      subset(meta_prot[meta_prot$Group !=4,], select= -c(1:2) ) ) %>% 
    tibble::rownames_to_column("gene") %>% 
    dplyr::mutate(S4_fdr=p.adjust(pvalue, "fdr"), S4_logFC=log2(mean.x/mean.y)  ) %>% 
    dplyr::arrange(pvalue) %>% 
    dplyr::left_join(g3_mean_tib)

#
meta_prot_3Group <- as.data.frame(meta_prot[,-c(1:2)])
rownames(meta_prot_3Group) <- meta_prot$eid

meta_prot$Group <- factor(meta_prot$Group, levels=c(0,1,2,3,4) )

all_subtype_DEA <- matrixTests::col_oneway_equalvar(meta_prot_3Group, meta_prot$Group ) %>% 
    tibble::rownames_to_column("gene") %>% 
    dplyr::mutate(fdr=p.adjust(pvalue, "fdr") ) %>% 
    dplyr::arrange(pvalue) %>% 
    dplyr::left_join(g3_mean_tib) 




S1_logFC <- S1_DEA[,c("gene", "S1_fdr", "S1_logFC")] 
S2_logFC <- S2_DEA[,c("gene", "S2_fdr", "S2_logFC")] 
S3_logFC <- S3_DEA[,c("gene", "S3_fdr", "S3_logFC")] 
S4_logFC <- S4_DEA[,c("gene", "S4_fdr", "S4_logFC", "C_average", "S1_average",
                      "S2_average", "S3_average", "S4_average")] 

all_subtype_average <- all_subtype_DEA[,c("gene", "C_average", "S1_average",
                                      "S2_average", "S3_average", "S4_average")]

metabolite_lib <- readxl::read_xlsx("./data/biobank.ndph.ox.ac.uk_showcase_ukb_docs_Nightingale_biomarker_groups.xlsx")

logFC_fdr <- dplyr::left_join(S1_logFC, S2_logFC) %>% 
    dplyr::left_join(S3_logFC) %>% 
    dplyr::left_join(S4_logFC) %>% 
    tidyr::separate_wider_delim(gene, delim="_f2", names=c("metabolite","field_id") ) %>% 
    dplyr::mutate(field_id=gsub("_.*","",field_id), field_id=as.numeric(paste0(2, field_id)) ) %>% 
    dplyr::left_join(metabolite_lib)


cutoff <- 0.3
logFC_fdr_sig <- dplyr::filter(logFC_fdr, abs(S1_logFC) > cutoff | abs(S2_logFC) > cutoff | 
                                   abs(S3_logFC) > cutoff | abs(S4_logFC) > cutoff ) %>% 
    dplyr::arrange(Group)
library(ComplexHeatmap)
library(circlize)
mat1 <- logFC_fdr_sig[,c("S1_logFC", "S2_logFC", "S3_logFC", "S4_logFC")] %>% as.matrix()
mat0 <- logFC_fdr_sig[,c("C_average", "S1_average", 
                         "S2_average", "S3_average", "S4_average")] %>% as.matrix()
mat0 = t(scale(t(mat0)))

mat2 <- logFC_fdr_sig[,c("S1_fdr", "S2_fdr", "S3_fdr", "S4_fdr")] %>% as.matrix()
rownames(mat0) <- rownames(mat1) <- rownames(mat2) <- logFC_fdr_sig$title

col_fun = colorRamp2(c(-1.2, 0, 1.16), c("#F0E442", "white", "#E69F00"))

ht0 <- Heatmap(mat0, name="abundance", left_annotation=rowAnnotation(Group=logFC_fdr_sig$Group,
                                                                     Subgroup=logFC_fdr_sig$Subgroup),
               row_names_side="left", col = col_fun, 
               cluster_rows=FALSE, show_row_dend=FALSE, cluster_columns=FALSE)
ht <- Heatmap(mat1, name="logFC", show_row_names=FALSE,
              cluster_rows=FALSE, show_row_dend=FALSE, cluster_columns=FALSE,
              cell_fun = function(j, i, x, y, w, h, fill) {
                  if(mat2[i, j] < 0.05) {
                      grid.text("*", x, y)
                  }
              })
draw(ht0+ht, padding = unit(c(2, 100, 2, 2), "mm"))

res_xlsx <- list(S1_DEA=S1_DEA, 
                 S2_DEA=S2_DEA,
                 S3_DEA=S3_DEA,
                 S4_DEA=S4_DEA,
                 all_subtype_DEA=all_subtype_DEA,
                 logFC_fdr=logFC_fdr,
                 logFC_fdr_sig=logFC_fdr_sig,
                 sampleSize=as.data.frame(table(meta_prot$Group)) )

openxlsx::write.xlsx(res_xlsx, file="./results/DEA_metabo.xlsx", asTable=TRUE )

save.image("./results/02_DE_metabo.RData")
