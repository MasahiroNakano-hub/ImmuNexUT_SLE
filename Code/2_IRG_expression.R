#---------------------------------------#
# 211124 Nakano
# IRG expression
#---------------------------------------#

##############################################
# IRG expressions of SLE and HC
##############################################
source("my.source.R")
options(stringsAsFactors=F)
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(RColorBrewer))

list1 = make_list("data/SLEinc_andHC/each_subset/LowExp_TMM_logCPM_Combat","_logCPM_Combat.txt",tag="count")
list1$subset = take_factor(list1$FILE_count,5:6,"_")
list2 = make_list("data/SLEinc_andHC/each_subset/each_subset_cond","_cond.txt",tag="cond")
list2$subset = take_factor(list2$FILE_cond,5:6,"_")
list = left_join(list1,list2,by="subset")

# Nat Immunol 2020 SLE PBMC scRNAseq ISG list
ISG_2020 = read.table_FT("data_ref/210608_ISG_pascual_sc2020_v2.txt")
celltype_reordered_v2 = c("B05_NCD4","B06_MCD4","B01_Th1","B02_Th2","B03_TH17","B04_Tfh","B09_Fra1","B07_aTreg","B10_Fra3","A01_NaiB","A03_UnswMB","A02_SwiMB","A04_DNB","C01_NCD8","C04_CmCD8","C05_EmCD8","C03_EffectorCD8","G01_NK","E02_CD16nMo","E01_CD16pMo","E04_Intermediate","E03_NonClassical","D01_mDC","D02_pDC","A05_PB","F01_Neu","F02_LDG")

for(iii in 1:length(celltype_reordered_v2)){
  subset_tmp = celltype_reordered_v2[iii]
  list_tmp   = list %>% filter(subset==subset_tmp)
  cond_tmp   = fread_FT(list_tmp$PATH_cond) 
  count_tmp  = fread_FT(list_tmp$PATH_count)
  count_ISG  = count_tmp %>% filter(Gene%in%ISG_2020$genes) %>% column_to_rownames("Gene") # Extract 100 ISGs

  cond_SLE   = cond_tmp %>% filter(disease=="1SLE")
  cond_HC    = cond_tmp %>% filter(disease=="0HC")

  ISG_SLE    = count_ISG %>% select(cond_SLE$name) %>% apply(.,1,mean) %>% as.data.frame() %>% rownames_to_column("Gene")
  ISG_HC     = count_ISG %>% select(cond_HC$name)  %>% apply(.,1,mean) %>% as.data.frame() %>% rownames_to_column("Gene")

  colnames(ISG_SLE)[2] = paste0(subset_tmp,"_SLE")
  colnames(ISG_HC)[2]  = paste0(subset_tmp,"_HC")
  
  ISG_2G = left_join(ISG_HC,ISG_SLE,by="Gene")
  
  if(iii==1){res_sum_2G=ISG_2G}else{res_sum_2G=full_join(res_sum_2G,ISG_2G,by="Gene")}
}

# For visualization, we prioritize CCL4L1, not CCL4L2
# Original article: CCL4
res_sum_2G_100 = res_sum_2G %>% filter(!Gene=="CCL4L2")
res_sum_2G_100 = res_sum_2G_100 %>% column_to_rownames("Gene") 

# Scale expression
min(res_sum_2G_100,na.rm=T)
# [1] 0.04151807 All expressed genes logCPM >0
res_sum_2G_100[is.na(res_sum_2G_100)]=0 # LowExp filtered genes are regarded as zero 
res_sum_2G_Z = res_sum_2G_100[intersect(ISG_2020$genes,rownames(res_sum_2G_100)),] %>% apply(.,1,scale) %>% t() %>% as.data.frame()
colnames(res_sum_2G_Z) = colnames(res_sum_2G_100)

write.table_n_2(res_sum_2G_Z,"Gene",paste0("ISGs/",today(),"_COI_HCSLE_ISG2020_100genes_Eachsubsets_meanexp.txt"))

##############################
colnames(res_sum_2G_Z) = rep(c("HC","SLE"),27)

col_fun = colorRamp2(c(-3,0,3), c("cyan","black","yellow"))

label_col = rep(c("Naive CD4","Mem CD4","Th1","Th2","Th17","Tfh",
              "Fr. I nTreg","Fr. II eTreg","Fr. III T",
              "Naive B","USM B","SM B","DN B",
              "Naive CD8","CM CD8","EM CD8","TEMRA CD8","NK",
              "CL Mono","CD16p Mono","Int Mono","NC Mono",
              "mDC","pDC","Plasmablast",
              "Neu","LDG"),each=2)
label_col = factor(label_col,levels=unique(label_col))

label_row = c(rep("G1",21),rep("G2",18),rep("G3",16),rep("G4",1),rep("G5",20),rep("G6",6),rep("G7",3),rep("G8",5),rep("G9",10))

p=Heatmap(as.matrix(res_sum_2G_Z),cluster_rows=FALSE,cluster_columns=FALSE,col=col_fun,rect_gp=gpar(col="black"),width=unit(30,"cm"),height=unit(30,"cm"),
                                  column_split=label_col,column_title_gp=gpar(fontsize=20),column_title_rot=90,column_title_side="bottom",
                                  row_split=label_row,row_title_gp=gpar(fontsize=20,col="blue"),row_title_rot=0,row_title_side="right",
                                  row_names_gp=gpar(fontsize=11,fontface="italic"),
                                  column_names_gp=gpar(fontsize=15,col=rep(c("#F8766D","#00BFC4"),27)),
                                  heatmap_legend_param=list(title="Scaled expression",at=c(-3,0,3),direction="horizontal",
                                                            title_position="topcenter",title_gp=gpar(fontsize=16,fontface="bold"),labels_gp=gpar(fontsize=16),
                                                            grid_height=unit(0.6,"cm"),legend_width=unit(4,"cm")))

pdf_3(paste0("ISGs/",today(),"_COI_HCSLE_ISG2020_100genes_Eachsubsets_meanexp.pdf"),h=16,w=18)
 p
dev.off()


# Retry row longer
p=Heatmap(as.matrix(res_sum_2G_Z),cluster_rows=FALSE,cluster_columns=FALSE,col=col_fun,rect_gp=gpar(col="black"),width=unit(30,"cm"),height=unit(48,"cm"),
                                  column_split=label_col,column_title_gp=gpar(fontsize=24),column_title_rot=90,column_title_side="bottom",
                                  row_split=label_row,row_title_gp=gpar(fontsize=24,col="blue"),row_title_rot=0,row_title_side="right",
                                  row_names_gp=gpar(fontsize=17.5,fontface="italic"),
                                  column_names_gp=gpar(fontsize=15,col=rep(c("#F8766D","#00BFC4"),27)),
                                  heatmap_legend_param=list(title="Scaled expression",at=c(-3,0,3),direction="horizontal",
                                                            title_position="topcenter",title_gp=gpar(fontsize=20,fontface="bold"),labels_gp=gpar(fontsize=20),
                                                            grid_height=unit(0.6,"cm"),legend_width=unit(4,"cm")))

pdf_3(paste0("Figure/",today(),"_COI_HCSLE_ISG2020_100genes_Eachsubsets_meanexp.pdf"),h=22,w=20)
 p
dev.off()




##############################
# Number of ISGs with highest in SLE-Neu or SLE-LDG
maxcol = apply(res_sum_2G_Z,1,function(x){which.max(x)})

length(which(maxcol%in%c(52,54)))
# [1] 54
table(maxcol)


sessionInfo()
R version 4.0.2 (2020-06-22)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)

Matrix products: default
BLAS:   /usr/local/package/r/4.0.2/lib64/R/lib/libRblas.so
LAPACK: /usr/local/package/r/4.0.2/lib64/R/lib/libRlapack.so

locale:
 [1] LC_CTYPE=ja_JP.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=ja_JP.UTF-8        LC_COLLATE=ja_JP.UTF-8    
 [5] LC_MONETARY=ja_JP.UTF-8    LC_MESSAGES=ja_JP.UTF-8   
 [7] LC_PAPER=ja_JP.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=ja_JP.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] grid      stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] circlize_0.4.13      ComplexHeatmap_2.6.2 RColorBrewer_1.1-2  
 [4] data.table_1.14.0    forcats_0.5.1        stringr_1.4.0       
 [7] dplyr_1.0.7          purrr_0.3.4          readr_1.4.0         
[10] tidyr_1.1.3          tibble_3.1.2         ggplot2_3.3.5       
[13] tidyverse_1.3.1     

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.6          lubridate_1.7.10    png_0.1-7          
 [4] ps_1.6.0            assertthat_0.2.1    digest_0.6.27      
 [7] utf8_1.2.1          R6_2.5.0            cellranger_1.1.0   
[10] backports_1.2.1     reprex_2.0.0        stats4_4.0.2       
[13] httr_1.4.2          pillar_1.6.1        GlobalOptions_0.1.2
[16] rlang_0.4.11        readxl_1.3.1        rstudioapi_0.13    
[19] S4Vectors_0.28.1    GetoptLong_1.0.5    labeling_0.4.2     
[22] munsell_0.5.0       broom_0.7.8         compiler_4.0.2     
[25] modelr_0.1.8        pkgconfig_2.0.3     BiocGenerics_0.36.1
[28] shape_1.4.6         tidyselect_1.1.1    IRanges_2.24.1     
[31] matrixStats_0.59.0  fansi_0.5.0         crayon_1.4.1       
[34] dbplyr_2.1.1        withr_2.4.2         jsonlite_1.7.2     
[37] gtable_0.3.0        lifecycle_1.0.0     DBI_1.1.1          
[40] magrittr_2.0.1      scales_1.1.1        cli_2.5.0          
[43] stringi_1.6.2       farver_2.1.0        fs_1.5.0           
[46] xml2_1.3.2          ellipsis_0.3.2      generics_0.1.0     
[49] vctrs_0.3.8         rjson_0.2.20        tools_4.0.2        
[52] Cairo_1.5-12.2      glue_1.4.2          hms_1.1.0          
[55] parallel_4.0.2      clue_0.3-59         colorspace_2.0-2   
[58] cluster_2.1.2       rvest_1.0.0         haven_2.4.1   














