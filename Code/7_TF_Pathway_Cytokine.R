#---------------------------------------#
# 211202 Nakano
# DEG biology 
#---------------------------------------#

source("my.source.R")
options(stringsAsFactors=F)
suppressPackageStartupMessages(library(msigdbr))
suppressPackageStartupMessages(library(ggsci))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(circlize))

celltype_reordered_27 = c("B05_NCD4","B06_MCD4","B01_Th1","B02_Th2","B03_TH17","B04_Tfh","B09_Fra1","B07_aTreg","B10_Fra3","C01_NCD8","C04_CmCD8","C05_EmCD8","C03_EffectorCD8","G01_NK","A01_NaiB","A03_UnswMB","A02_SwiMB","A04_DNB","A05_PB","E02_CD16nMo","E01_CD16pMo","E04_Intermediate","E03_NonClassical","D01_mDC","D02_pDC","F01_Neu","F02_LDG")
label = c(`B05_NCD4`="Naive CD4",`B06_MCD4`="Mem CD4",`B01_Th1`="Th1",`B02_Th2`="Th2",`B03_TH17`="Th17",`B04_Tfh`="Tfh",
          `B09_Fra1`="Fr. I nTreg",`B07_aTreg`="Fr. II eTreg",`B10_Fra3`="Fr. III T",
          `C01_NCD8`="Naive CD8",`C04_CmCD8`="CM CD8",`C05_EmCD8`="EM CD8",`C03_EffectorCD8`="TEMRA CD8",`G01_NK`="NK",
          `A01_NaiB`="Naive B",`A03_UnswMB`="USM B",`A02_SwiMB`="SM B",`A04_DNB`="DN B",`A05_PB`="Plasmablast",
          `E02_CD16nMo`="CL Mono",`E01_CD16pMo`="CD16p Mono",`E04_Intermediate`="Int Mono",`E03_NonClassical`="NC Mono",
          `D01_mDC`="mDC",`D02_pDC`="pDC",
          `F01_Neu`="Neu",`F02_LDG`="LDG")
label_df = as.data.frame(label) %>% rownames_to_column("subset")
labeller = as_labeller(label)
celltype_corresp=fread_FT("data_ref/COI_27subset_color_list.txt")

col1=c("Disease-state"="#0072B5FF","Disease-activity"="#BC3C29FF","Both significant"="#E18727FF","Not significant"="#cccccc")
col2=c("CD4"="#cc0010","CD8"="#8b7355","NK"="#e0d51a","B"="#c4f20b","Monocyte"="#1d2088","DC"="#a4fbfa","Neutrophil"="#cccccc")

col_full=c("#921b1b","#ff9994","#cc0010","#ff0000","#e6a88f","#e66557","#feaf53","#f77308","#fde0bd",
          "#995522","#cc9955","#eebb77","#bb9988",
          "#e0d51a",
          "#00561f","#228b22","#8fc31f","#33dd33","#c4f20b",
          "#1d2088","#1e90ff","#a4bdfc","#4d79f7",
          "#0bdaf2","#a4fbfa",
          "#cccccc","#999999")

column_labels = structure(label,names=celltype_reordered_27)

############################################################################################################################################################
# Enrichment analysis
############################################################################################################################################################

list1         = make_list("res/stateactivity","_DEGlist_FDR005.txt")
list1$subset  = take_factor(list1$FILE,3:4,"_")
list1$type    = take_factor(list1$FILE,5,"_")
list1 = list1 %>% filter(type%in%c("InactivevsHC","HDAvsInactive"))
list1$type    = ifelse(list1$type=="InactivevsHC","state","activity")

for(kkk in 1:nrow(list1)){
    subset_tmp    = list1$subset[kkk]
    res_tmp    = fread_FT(list1$PATH[kkk])
    DEG_tmp     = res_tmp$genes
    if(kkk==1){union_sum=DEG_tmp}else{union_sum=union(union_sum,DEG_tmp)}
   }

union_sum = data.frame(genes=union_sum)
write.table_FT_2(union_sum,paste0("res/",today(),"_COI_27subsets_stateactivity_DEGunionlist.txt"))

list2         = make_list("res/treatment","_DEGlist_FDR005.txt")
list2$subset  = take_factor(list2$FILE,3:4,"_")
list2$type    = take_factor(list2$FILE,5,"_")

list4 = list2 %>% filter(type=="TAC")%>% filter(subset%in%c("B10_Fra3"))
list5 = list2 %>% filter(type=="MMF")%>% filter(subset%in%c("B01_Th1","C04_CmCD8","A05_PB"))

list = bind_rows(list(list1,list3,list4,list5))
list$type_subset = paste0(list$type,"_",list$subset)
list$back_PATH = "/home/mnakano/ws/2021/211129_COI_SLE_HC_run27_edgeR/res/211202_COI_27subsets_stateactivity_DEGunionlist.txt"

write.table_FT_2(list,paste0("tmp_job_list/",today(),"_COI_stateactivity_Tx_DEGlist.txt"))
######################################################################################################################

############################################################################################################################################################
# TFs
# Example: One cell type disease-state and activity
############################################################################################################################################################

LIST = "tmp_job_list/211202_COI_state_activity_DEGlist.txt"
## task_id=1

list       = fread_FT(LIST)
TARGET_PATH = list[task_id,1]
subset_tmp = list[task_id,3]
type_tmp = list[task_id,4]
back_PATH = list[task_id,6]

suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(AnnotationDbi))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(GO.db))
suppressPackageStartupMessages(library(enrichplot))
suppressPackageStartupMessages(library(topGO))
suppressPackageStartupMessages(library(Rgraphviz))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(msigdbr))

########################## 
#TARGET_NAME        = TARGET_NAME
target_gene_list   = read.table_FT(TARGET_PATH)[,1]
back_gene_list     = read.table_FT(back_PATH)[,1]
########################## 

gene_full_list   = read.table_FT("data_ref/200303_HTseq_UCSC_hg38_wo_chrY_gene_w_ENTREZ_ID_list.txt")
gene_full_list_2 = gene_full_list[!is.na(gene_full_list$ENTREZID),]

target_gene_list_2 = data.frame(SYMBOL = target_gene_list) %>% inner_join(.,gene_full_list_2,by="SYMBOL")
back_gene_list_2 = data.frame(SYMBOL = back_gene_list) %>% inner_join(.,gene_full_list_2,by="SYMBOL")

paste0("TARGET GENE num is ",length(target_gene_list))          %>% print()
paste0("Analyzed TARGET GENE num is ",nrow(target_gene_list_2)) %>% print()

paste0("BACK GENE num is ",length(back_gene_list))          %>% print()
paste0("Analyzed BACK GENE num is ",nrow(back_gene_list_2)) %>% print()


# HALLMARK
HALLMARK = msigdbr(species="Homo sapiens",category="H") %>% 
           dplyr::select(gs_name,entrez_gene)

res1 = enricher(gene=target_gene_list_2$ENTREZID%>%as.character(),universe=back_gene_list_2$ENTREZID%>%as.character(),TERM2GENE=HALLMARK,maxGSSize=1000,pvalueCutoff=1,qvalueCutoff=1)
write.table_FT_2(res1,paste0("ORA/HALLMARK/",type_tmp,"/",today(),"_",subset_tmp,"_",type_tmp,"_ORA_HALLMARK.txt"))

res1_dot = dotplot(res1, showCategory=20)
pdf_3(paste0("ORA/HALLMARK/",type_tmp,"/",today(),"_",subset_tmp,"_",type_tmp,"_ORA_HALLMARK.pdf"),h=8,w=12)
 plot(res1_dot)
dev.off()

# C3_GTRD
C3_GTRD = msigdbr(species="Homo sapiens",category="C3",subcategory = "TFT:GTRD") %>% 
           dplyr::select(gs_name,entrez_gene)

res2 = enricher(gene=target_gene_list_2$ENTREZID%>%as.character(),universe=back_gene_list_2$ENTREZID%>%as.character(),TERM2GENE=C3_GTRD,maxGSSize=3000,pvalueCutoff=1,qvalueCutoff=1)
write.table_FT_2(res2,paste0("ORA/C3GTRD/",type_tmp,"/",today(),"_",subset_tmp,"_",type_tmp,"_ORA_C3GTRD.txt"))

res2_dot = dotplot(res2, showCategory=20)
pdf_3(paste0("ORA/C3GTRD/",type_tmp,"/",today(),"_",subset_tmp,"_",type_tmp,"_ORA_C3GTRD.pdf"),h=8,w=12)
 plot(res2_dot)
dev.off()

# C3_Legacy
C3_Legacy = msigdbr(species="Homo sapiens",category="C3",subcategory = "TFT:TFT_Legacy") %>% 
           dplyr::select(gs_name,entrez_gene)

res3 = enricher(gene=target_gene_list_2$ENTREZID%>%as.character(),universe=back_gene_list_2$ENTREZID%>%as.character(),TERM2GENE=C3_Legacy,maxGSSize=3000,pvalueCutoff=1,qvalueCutoff=1)
write.table_FT_2(res3,paste0("ORA/C3Legacy/",type_tmp,"/",today(),"_",subset_tmp,"_",type_tmp,"_ORA_C3Legacy.txt"))


# KEGG
enrichKEGG  = enrichKEGG(gene     = target_gene_list_2$ENTREZID %>% as.character(), 
                         universe = back_gene_list_2$ENTREZID   %>% as.character(), 
                         organism     = 'hsa',
                         maxGSSize    = 2000, 
                         pvalueCutoff = 1,
                         qvalueCutoff = 1)

write.table_FT_2(enrichKEGG ,paste0("ORA/KEGG/",type_tmp,"/",today(),"_",subset_tmp,"_",type_tmp,"_ORA_KEGG.txt"))

KEGG_dot = dotplot(enrichKEGG, showCategory=20)
pdf_3(paste0("ORA/KEGG/",type_tmp,"/",today(),"_",subset_tmp,"_",type_tmp,"_ORA_KEGG.pdf"),h=8,w=12)
 plot(KEGG_dot)
dev.off()

######################################################################################################################

######################################################################################################################
# C3 data: Combine and recalculate FDR
list = make_list("ORA",".txt")
list$subset = take_factor(list$FILE,2:3,"_")
list$type   = take_factor(list$FILE,4,"_")
list$anno   = take_factor(list$FILE,6,"_") %>% gsub(".txt","",.) 

list1=list%>%filter(type%in%c("state","activity"))
list2 = list1 %>% filter(anno%in%c("C3GTRD","C3Legacy"))

for(kkk in 1:length(celltype_reordered_27)){
   subset_tmp = celltype_reordered_27[kkk]
   list2_tmp1 = list2 %>% filter(subset==subset_tmp) %>% filter(anno=="C3GTRD")
   list2_tmp2 = list2 %>% filter(subset==subset_tmp) %>% filter(anno=="C3Legacy")

   state_tmp1 = fread_FT(list2_tmp1$PATH[2])
   paste0(subset_tmp,"; state; C3GTRD; ",nrow(state_tmp1))
   state_tmp2   = fread_FT(list2_tmp2$PATH[2])
   paste0(subset_tmp,"; state; C3Legacy; ",nrow(state_tmp2))

   act_tmp1 = fread_FT(list2_tmp1$PATH[1])
   paste0(subset_tmp,"; activity; C3GTRD; ",nrow(act_tmp1))
   act_tmp2   = fread_FT(list2_tmp2$PATH[1])
   paste0(subset_tmp,"; activity; C3Legacy; ",nrow(act_tmp2))

   state_sum=rbind(state_tmp1,state_tmp2)
   state_sum=state_sum[order(state_sum$pvalue),]
   state_sum$p.adjust_sum=p.adjust(state_sum$pvalue,method="BH")

   act_sum=rbind(act_tmp1,act_tmp2)
   act_sum=act_sum[order(act_sum$pvalue),]
   act_sum$p.adjust_sum=p.adjust(act_sum$pvalue,method="BH")

   write.table_FT_2(state_sum,paste0("ORA/C3TFsum/state/",today(),"_",subset_tmp,"_state_ORA_C3TFsum.txt"))
   write.table_FT_2(act_sum,paste0("ORA/C3TFsum/activity/",today(),"_",subset_tmp,"_activity_ORA_C3TFsum.txt"))
}

#### TFsum significant TF list
list_tmp = make_list("ORA/C3TFsum","_C3TFsum.txt")
list_tmp$subset = take_factor(list_tmp$FILE,2:3,"_")
list_tmp$type   = take_factor(list_tmp$FILE,4,"_")
list_tmp$anno   = take_factor(list_tmp$FILE,6,"_") %>% gsub(".txt","",.) 

for(kkk in 1:nrow(list_tmp)){
   sig_tmp   = fread_FT(list_tmp$PATH[kkk]) %>% filter(p.adjust_sum<0.05)
   sig_tmp   = sig_tmp$Description
   if(kkk==1){sig_sum=sig_tmp}else{sig_sum=union(sig_sum,sig_tmp)}
   } 
write.table_FT_2(data.frame(TF=sig_sum),paste0("ORA/",today(),"_ORA_C3TFsum_sig_TFlist.txt"))
length(sig_sum)
# [1] 299

# manual curation to remove redundant TFs
sig_sum_curated=fread_FT("ORA/211204_ORA_C3TFsum_sig_TFlist_curated.txt")

# Each signature Top 3
 for(kkk in 1:length(celltype_reordered_27)){
    subset_tmp = celltype_reordered_27[kkk]
    list_tmp2  = list_tmp %>% filter(subset==subset_tmp)
    list_tmp3  = list_tmp2 %>% filter(type=="state")
    list_tmp4  = list_tmp2 %>% filter(type=="activity")
       
    data_tmp3  = fread_FT(list_tmp3$PATH) %>% filter(Description%in%sig_sum_curated$TF)
    data_tmp3  = data_tmp3[1:3,] %>% filter(p.adjust_sum<0.05) 

    data_tmp4  = fread_FT(list_tmp4$PATH) %>% filter(Description%in%sig_sum_curated$TF)
    data_tmp4  = data_tmp4[1:3,] %>% filter(p.adjust_sum<0.05) 

    union_tmp=union(data_tmp3$Description,data_tmp4$Description)
    if(kkk==1){union_sum=union_tmp}else{union_sum=union(union_sum,union_tmp)}
}

write.table_FT_2(data.frame(TF=union_sum),paste0("ORA/",today(),"_ORA_C3TFsum_sig_Top3TFs.txt"))

# Clustering limited to Top3
for(kkk in 1:length(celltype_reordered_27)){
    subset_tmp = celltype_reordered_27[kkk]
    list_tmp2  = list_tmp %>% filter(subset==subset_tmp)
    list_tmp3  = list_tmp2 %>% filter(type=="state")
    list_tmp4  = list_tmp2 %>% filter(type=="activity")
       
    data_tmp3  = fread_FT(list_tmp3$PATH) %>% select(Description,pvalue)
    colnames(data_tmp3)[2] = "state_P"
    if(nrow(data_tmp3)==0){
    data_tmp3 = data.frame(Description=NA,state_P=NA)
    }else{data_tmp3 = data_tmp3}

    data_tmp4  = fread_FT(list_tmp4$PATH) %>% select(Description,pvalue)
    colnames(data_tmp4)[2] = "activity_P"
    if(nrow(data_tmp4)==0){
    data_tmp4 = data.frame(Description=NA,activity_P=NA)
    }else{data_tmp4 = data_tmp4}

    data_tmp5 = full_join(data_tmp3,data_tmp4,by="Description") %>% filter(Description%in%union_sum)
    if(nrow(data_tmp5)==0){
    data_tmp5 = data.frame(Description=NA,state_P=NA,activity_P=NA)
    }else{data_tmp5 = data_tmp5}
    data_tmp5$state_P[is.na(data_tmp5$state_P)]=1
    data_tmp5$activity_P[is.na(data_tmp5$activity_P)]=1
    data_tmp5$subset      = subset_tmp
    data_tmp5=data_tmp5%>%mutate(state_logP=-log10(state_P))%>%
                          mutate(activity_logP=-log10(activity_P))
   if(kkk==1){data_sum2=data_tmp5}else{data_sum2=rbind(data_sum2,data_tmp5)}
  }

    data_sum2 = data_sum2[!is.na(data_sum2$Description),]

    data_sum3 = data_sum2 %>% select(Description,subset,state_P,activity_P,state_logP,activity_logP)
    write.table_FT_2(data_sum3,paste0("ORA/heatmap/",today(),"_ORA_C3TFsum_eachTop3_logP.txt"))

    data_sum4 = data_sum3 %>% select(Description,subset,activity_logP) %>%
                              pivot_wider(names_from="subset",values_from="activity_logP") %>% column_to_rownames("Description")

    hc_row = hclust(dist(data_sum4), method="ward.D2")
    pdf_3(paste0("ORA/heatmap/",today(),"_ORA_C3TFsum_eachTop3_logP_hrow.pdf"),h=5,w=10)
     plot(hc_row)
    dev.off()

    hc_col = hclust(dist(t(data_sum4)), method="ward.D2")
    pdf_3(paste0("ORA/heatmap/",today(),"_ORA_C3TFsum_eachTop3_logP_hcol.pdf"),h=5,w=10)
     plot(hc_col)
    dev.off()

# rotate orders maintaining cluster structures
    label_row=hc_row$labels[hc_row$order]
    label_col=hc_col$labels[hc_col$order]
    label_row = label_row[c(1:6,10:14,7:9,42:15)]   
    label_col = label_col[c(6:1,16:17,21:18,15:7,22:27)]   

# Heatmap
 for(kkk in 1:length(celltype_reordered_27)){
    subset_tmp = celltype_reordered_27[kkk]
    list_tmp2  = list_tmp %>% filter(subset==subset_tmp)
    list_tmp3  = list_tmp2 %>% filter(type=="state")
    list_tmp4  = list_tmp2 %>% filter(type=="activity")
       
    data_tmp3  = fread_FT(list_tmp3$PATH) %>% select(Description,p.adjust_sum)
    colnames(data_tmp3)[2] = "state_FDR"
    if(nrow(data_tmp3)==0){
    data_tmp3 = data.frame(Description=NA,state_FDR=NA)
    }else{data_tmp3 = data_tmp3}

    data_tmp4  = fread_FT(list_tmp4$PATH) %>% select(Description,p.adjust_sum)
    colnames(data_tmp4)[2] = "activity_FDR"
    if(nrow(data_tmp4)==0){
    data_tmp4 = data.frame(Description=NA,activity_FDR=NA)
    }else{data_tmp4 = data_tmp4}

    data_tmp5 = full_join(data_tmp3,data_tmp4,by="Description") %>% filter(Description%in%union_sum)
    if(nrow(data_tmp5)==0){
    data_tmp5 = data.frame(Description=NA,state_FDR=NA,activity_FDR=NA)
    }else{data_tmp5 = data_tmp5}
    data_tmp5$state_FDR[is.na(data_tmp5$state_FDR)]=1
    data_tmp5$activity_FDR[is.na(data_tmp5$activity_FDR)]=1
    data_tmp5$subset      = subset_tmp
    data_tmp5$siggroup005 = ifelse(data_tmp5$state_FDR<0.05&data_tmp5$activity_FDR<0.05,"Both significant",
                                ifelse(data_tmp5$state_FDR<0.05,"Disease-state only",
                                ifelse(data_tmp5$activity_FDR<0.05,"Disease-activity only",NA)))
   if(kkk==1){data_sum2=data_tmp5}else{data_sum2=rbind(data_sum2,data_tmp5)}
  }
data_sum2 = data_sum2[!is.na(data_sum2$Description),]

data_sum3 = data_sum2 %>% select(Description,subset,siggroup005) %>%
                          pivot_wider(names_from="subset",values_from="siggroup005") %>% column_to_rownames("Description")

write.table_n_2(data_sum3,"anno",paste0("ORA/heatmap/",today(),"_ORA_C3TFsum_4class_eachTop3.txt"))

# reorder and rename
data_sum3_ordered=data_sum3[label_row,label_col]
rownames(data_sum3_ordered)=gsub("_TARGET_GENES","",rownames(data_sum3_ordered))
rownames(data_sum3_ordered)=gsub("_01","",rownames(data_sum3_ordered))
rownames(data_sum3_ordered)=gsub("_Q6","",rownames(data_sum3_ordered))
rownames(data_sum3_ordered)[26]="NFAT"
rownames(data_sum3_ordered)[27]="ETS2"
rownames(data_sum3_ordered)[30]="SRF"
rownames(data_sum3_ordered)[31]="NFY"

cm1 = ColorMapping(colors=c("Disease-state only"="#0072B5FF","Disease-activity only"="#BC3C29FF","Both significant"="#E18727FF"))
    
label_col2 = left_join(data.frame(subset=label_col),celltype_corresp%>%select(subset,lineage),by="subset")
label_col2$lineage = factor(label_col2$lineage,levels=c("CD4","CD8","NK","B","Monocyte","DC","Neutrophil"))

ha = HeatmapAnnotation(lineage=label_col2$lineage,col=list(lineage=col2),show_annotation_name=F,
                   height=unit(0.15,"cm"),simple_anno_size_adjust=T,show_legend=F,
                   annotation_legend_param=list(title_gp=gpar(fontsize=12),labels_gp=gpar(fontsize=12)))

# Figure: row left
p=Heatmap(as.matrix(data_sum3_ordered),cluster_rows=FALSE,cluster_columns=FALSE,col=cm1,na_col="white",
                              rect_gp=gpar(col="#808080",lwd=0.3),width=unit(11,"cm"),height=unit(12,"cm"),
                              row_names_gp=gpar(fontsize=10,fontface="italic"),row_names_side=c("left"),
                              column_names_gp=gpar(fontsize=12),column_labels=column_labels[colnames(data_sum3_ordered)],
                              bottom_annotation = ha,
                              heatmap_legend_param=list(title="Enrichment",title_gp=gpar(fontsize=12),labels_gp=gpar(fontsize=12)))
pdf_3(paste0("Figure/",today(),"_ORA_C3TFsum_4class_eachTop3_fig.pdf"),h=6,w=7.2)
 draw(p)
dev.off()

######################################################################################################################
# Count the significant enrichments
######################################################################################################################

### C3
list_tmp = make_list("ORA/C3TFsum",".txt")
list_tmp$subset = take_factor(list_tmp$FILE,2:3,"_")
list_tmp$type   = take_factor(list_tmp$FILE,4,"_")
list_tmp$anno   = take_factor(list_tmp$FILE,6,"_") %>% gsub(".txt","",.) 

for(kkk in 1:nrow(list_tmp)){
   subset_tmp=list_tmp$subset[kkk]
   type_tmp  =list_tmp$type[kkk]
   sig_tmp   = fread_FT(list_tmp$PATH[kkk]) %>% filter(p.adjust_sum<0.05)
   signum_tmp=data.frame(subset=subset_tmp,DEGtype=type_tmp,signum=nrow(sig_tmp))
   if(kkk==1){signum_sum=signum_tmp}else{signum_sum=union(signum_sum,signum_tmp)}
   } 
write.table_FT_2(signum_sum,paste0("ORA/C3TFsum/",today(),"_ORA_C3TFsum_sig_sum.txt"))

signum_sum_act=signum_sum%>%filter(DEGtype=="activity")
sum(signum_sum_act$signum)
# [1] 862
signum_sum_state=signum_sum%>%filter(DEGtype=="state")
sum(signum_sum_state$signum)
# [1] 366

# reorder by activity signum
signum_sum_act=signum_sum_act[order(signum_sum_act$signum),]
signum_sum$subset  =factor(signum_sum$subset,levels=signum_sum_act$subset)
signum_sum$DEGtype = factor(signum_sum$DEGtype,levels=c("state","activity"))
levels(signum_sum$DEGtype)=list(`Disease-state`="state",`Disease-activity`="activity")
signum_sum=signum_sum%>%left_join(.,celltype_corresp%>%select(subset,lineage),by="subset")
signum_sum$lineage = factor(signum_sum$lineage,levels=c("CD4","CD8","NK","B","Monocyte","DC","Neutrophil"))

# Figure
# Include both sig
for(kkk in 1:length(celltype_reordered_27)){
   subset_tmp=celltype_reordered_27[kkk]
   list_tmp_2=list_tmp%>%filter(subset==subset_tmp)
   state_tmp=fread_FT(list_tmp_2$PATH[2])%>%select(Description,pvalue,p.adjust_sum)
   colnames(state_tmp)[2:3]=paste0("state_",colnames(state_tmp)[2:3])
   act_tmp=fread_FT(list_tmp_2$PATH[1])%>%select(Description,pvalue,p.adjust_sum)
   colnames(act_tmp)[2:3]=paste0("activity_",colnames(act_tmp)[2:3])
   res_tmp=left_join(state_tmp,act_tmp,by="Description")%>%
           mutate(siggroup005=ifelse(state_p.adjust_sum<0.05&activity_p.adjust_sum<0.05,"Both significant",
                              ifelse(state_p.adjust_sum<0.05,"Disease-state only",
                              ifelse(activity_p.adjust_sum<0.05,"Disease-activity only","Not significant"))))
   table_tmp=table_freq(res_tmp$siggroup005) %>%
             mutate(subset=subset_tmp)
   if(kkk==1){table_sum=table_tmp}else{table_sum=union(table_sum,table_tmp)}
   } 
write.table_FT_2(table_sum,paste0("ORA/C3TFsum/",today(),"_ORA_C3TFsum_sig_sum_wbothinfo.txt"))

# reorder by clustering
colnames(table_sum)[1]="Enrichment"
table_sum=table_sum%>%filter(!Enrichment=="Not significant")
table_sum$subset  =factor(table_sum$subset,levels=label_col)
table_sum=table_sum[order(table_sum$subset),]
table_sum$Enrichment = factor(table_sum$Enrichment,levels=c("Both significant","Disease-state only","Disease-activity only"))
table_sum=table_sum%>%left_join(.,celltype_corresp%>%select(subset,lineage),by="subset")
table_sum$lineage = factor(table_sum$lineage,levels=c("CD4","CD8","NK","B","Monocyte","DC","Neutrophil"))
table_sum$subset  =factor(table_sum$subset,levels=label_col)

 p = ggplot()+
     geom_bar(data=table_sum, aes(x=subset,y=Freq,fill=Enrichment),stat="identity",position="stack")+
     theme_classic()+
     scale_fill_manual(values=c("#E18727FF","#0072B5FF","#BC3C29FF"))+
     theme(axis.text.x=element_blank(),
           axis.text.y=element_text(colour="black",size=14),
           axis.title.x=element_blank(),
           axis.title.y=element_text(colour="black",size=14),
           axis.ticks.x=element_blank(),
           plot.title=element_blank(),
      　　　legend.position="right",
           legend.title=element_text(colour="black",size=14),
           legend.text=element_text(colour="black",size=14))+
     scale_x_discrete(labels= label)+
　　　labs(y="Number of\nsignificant\nTF annotations")

    pdf_3(paste0("Figure/",today(),"_ORA_C3TFsum_sig_sum_wboth.pdf"),h=1.5,w=8.5)
     plot(p)
    dev.off()

############# Supple
list1=make_list("ORA/C3TFsum/state","_ORA_C3TFsum.txt")
list1$subset=take_factor(list1$FILE,2:3,"_")

for(iii in 1:nrow(list1)){
    subset_tmp=list1$subset[iii]
    res_tmp=fread_FT(list1$PATH[iii])
    res_tmp$subset=subset_tmp
    if(iii==1){res_sum1=res_tmp}else{res_sum1=rbind(res_sum1,res_tmp)}
}
res_sum1$DEGtype="Disease-state"

list2=make_list("ORA/C3TFsum/activity","_ORA_C3TFsum.txt")
list2$subset=take_factor(list2$FILE,2:3,"_")

for(iii in 1:nrow(list2)){
    subset_tmp=list2$subset[iii]
    res_tmp=fread_FT(list2$PATH[iii])
    res_tmp$subset=subset_tmp
    if(iii==1){res_sum2=res_tmp}else{res_sum2=rbind(res_sum2,res_tmp)}
}
res_sum2$DEGtype="Disease-activity"

res_sum_all=rbind(res_sum1,res_sum2)

res_sum_all$subset=factor(res_sum_all$subset,levels=celltype_reordered_27)
res_sum_all$DEGtype=factor(res_sum_all$DEGtype,levels=c("Disease-state","Disease-activity"))
res_sum_all=res_sum_all[order(res_sum_all$subset),]
res_sum_all=res_sum_all[order(res_sum_all$DEGtype),]

res_sum_lim=res_sum_all%>%left_join(.,celltype_corresp[,c(1,2)],by="subset")%>%
            select(DEGtype,label,Description,pvalue,p.adjust_sum)%>%filter(pvalue<0.05)
res_sum_lim$pvalue = formatC(res_sum_lim$pvalue,digits=2)
res_sum_lim$p.adjust_sum  = formatC(res_sum_lim$p.adjust_sum,digits=2)

colnames(res_sum_lim)=c("Signature","Cell type","Annotation","Pvalue","FDR")
write.table_FT_2(res_sum_lim,paste0(today(),"_ORA_C3TF_27subsets_stateactivity_Supple.txt"))

################################################################
### HALLMARK
list_tmp = make_list("ORA/HALLMARK",".txt")
list_tmp$subset = take_factor(list_tmp$FILE,2:3,"_")
list_tmp$type   = take_factor(list_tmp$FILE,4,"_")
list_tmp$anno   = take_factor(list_tmp$FILE,6,"_") %>% gsub(".txt","",.) 
list_tmp=list_tmp %>% filter(type%in%c("state","activity"))

for(kkk in 1:nrow(list_tmp)){
   subset_tmp=list_tmp$subset[kkk]
   type_tmp  =list_tmp$type[kkk]
   sig_tmp   = fread_FT(list_tmp$PATH[kkk]) %>% filter(p.adjust<0.05)
   signum_tmp=data.frame(subset=subset_tmp,DEGtype=type_tmp,signum=nrow(sig_tmp))
   if(kkk==1){signum_sum=signum_tmp}else{signum_sum=union(signum_sum,signum_tmp)}
   } 
write.table_FT_2(signum_sum,paste0("ORA/HALLMARK/",today(),"_ORA_HALLMARK_sig_sum.txt"))

signum_sum_act=signum_sum%>%filter(DEGtype=="activity")
sum(signum_sum_act$signum)
# [1] 95
signum_sum_state=signum_sum%>%filter(DEGtype=="state")
sum(signum_sum_state$signum)
# [1] 181

# Unique pathways
for(kkk in 1:nrow(list_tmp)){
   sig_tmp   = fread_FT(list_tmp$PATH[kkk]) %>% filter(p.adjust<0.05)
   sig_tmp   = sig_tmp$Description
   if(kkk==1){sig_sum=sig_tmp}else{sig_sum=union(sig_sum,sig_tmp)}
   } 
length(sig_sum)
# 24

### KEGG
list_tmp = make_list("ORA/KEGG",".txt")
list_tmp$subset = take_factor(list_tmp$FILE,2:3,"_")
list_tmp$type   = take_factor(list_tmp$FILE,4,"_")
list_tmp$anno   = take_factor(list_tmp$FILE,6,"_") %>% gsub(".txt","",.) 
list_tmp=list_tmp %>% filter(type%in%c("state","activity"))

for(kkk in 1:nrow(list_tmp)){
   subset_tmp=list_tmp$subset[kkk]
   type_tmp  =list_tmp$type[kkk]
   sig_tmp   = fread_FT(list_tmp$PATH[kkk]) %>% filter(p.adjust<0.05)
   signum_tmp=data.frame(subset=subset_tmp,DEGtype=type_tmp,signum=nrow(sig_tmp))
   if(kkk==1){signum_sum=signum_tmp}else{signum_sum=union(signum_sum,signum_tmp)}
   } 
write.table_FT_2(signum_sum,paste0("ORA/KEGG/",today(),"_ORA_KEGG_sig_sum.txt"))

signum_sum_act=signum_sum%>%filter(DEGtype=="activity")
sum(signum_sum_act$signum)
# [1] 220
signum_sum_state=signum_sum%>%filter(DEGtype=="state")
sum(signum_sum_state$signum)
# [1] 554

# Unique pathways
for(kkk in 1:nrow(list_tmp)){
   sig_tmp   = fread_FT(list_tmp$PATH[kkk]) %>% filter(p.adjust<0.05)
   sig_tmp   = sig_tmp$Description
   if(kkk==1){sig_sum=sig_tmp}else{sig_sum=union(sig_sum,sig_tmp)}
   } 
length(sig_sum)
# 134

####################
# Include both sig
### HALLMARK
list_tmp = make_list("ORA/HALLMARK",".txt")
list_tmp$subset = take_factor(list_tmp$FILE,2:3,"_")
list_tmp$type   = take_factor(list_tmp$FILE,4,"_")
list_tmp$anno   = take_factor(list_tmp$FILE,6,"_") %>% gsub(".txt","",.) 
list_tmp=list_tmp %>% filter(type%in%c("state","activity"))

for(kkk in 1:length(celltype_reordered_27)){
   subset_tmp=celltype_reordered_27[kkk]
   list_tmp_2=list_tmp%>%filter(subset==subset_tmp)
   state_tmp=fread_FT(list_tmp_2$PATH[2])%>%select(Description,pvalue,p.adjust)
   colnames(state_tmp)[2:3]=paste0("state_",colnames(state_tmp)[2:3])
   act_tmp=fread_FT(list_tmp_2$PATH[1])%>%select(Description,pvalue,p.adjust)
   colnames(act_tmp)[2:3]=paste0("activity_",colnames(act_tmp)[2:3])
   res_tmp=left_join(state_tmp,act_tmp,by="Description")%>%
           mutate(siggroup005=ifelse(state_p.adjust<0.05&activity_p.adjust<0.05,"Both significant",
                              ifelse(state_p.adjust<0.05,"Disease-state only",
                              ifelse(activity_p.adjust<0.05,"Disease-activity only","Not significant"))))
   table_tmp=table_freq(res_tmp$siggroup005) %>%
             mutate(subset=subset_tmp)
   if(kkk==1){table_sum1=table_tmp}else{table_sum1=union(table_sum1,table_tmp)}
   } 
write.table_FT_2(table_sum1,paste0("ORA/HALLMARK/",today(),"_ORA_HALLMARK_sig_sum_wbothinfo.txt"))

### KEGG
list_tmp = make_list("ORA/KEGG",".txt")
list_tmp$subset = take_factor(list_tmp$FILE,2:3,"_")
list_tmp$type   = take_factor(list_tmp$FILE,4,"_")
list_tmp$anno   = take_factor(list_tmp$FILE,6,"_") %>% gsub(".txt","",.) 
list_tmp=list_tmp %>% filter(type%in%c("state","activity"))

for(kkk in 1:length(celltype_reordered_27)){
   subset_tmp=celltype_reordered_27[kkk]
   list_tmp_2=list_tmp%>%filter(subset==subset_tmp)
   state_tmp=fread_FT(list_tmp_2$PATH[2])%>%select(Description,pvalue,p.adjust)
   colnames(state_tmp)[2:3]=paste0("state_",colnames(state_tmp)[2:3])
   act_tmp=fread_FT(list_tmp_2$PATH[1])%>%select(Description,pvalue,p.adjust)
   colnames(act_tmp)[2:3]=paste0("activity_",colnames(act_tmp)[2:3])
   res_tmp=left_join(state_tmp,act_tmp,by="Description")%>%
           mutate(siggroup005=ifelse(state_p.adjust<0.05&activity_p.adjust<0.05,"Both significant",
                              ifelse(state_p.adjust<0.05,"Disease-state only",
                              ifelse(activity_p.adjust<0.05,"Disease-activity only","Not significant"))))
   table_tmp=table_freq(res_tmp$siggroup005) %>%
             mutate(subset=subset_tmp)
   if(kkk==1){table_sum2=table_tmp}else{table_sum2=union(table_sum2,table_tmp)}
   } 
write.table_FT_2(table_sum2,paste0("ORA/KEGG/",today(),"_ORA_KEGG_sig_sum_wbothinfo.txt"))

table_sum1$anno="HALLMARK"
table_sum2$anno="KEGG"
table_sum=rbind(table_sum1,table_sum2)

colnames(table_sum)[1]="Enrichment"
table_sum=table_sum%>%filter(!Enrichment=="Not significant")
table_sum$subset  =factor(table_sum$subset,levels=celltype_reordered_27)
table_sum$Enrichment = factor(table_sum$Enrichment,levels=c("Both significant","Disease-state only","Disease-activity only"))

 p = ggplot()+
     geom_bar(data=table_sum, aes(x=subset,y=Freq,fill=Enrichment),stat="identity",position="stack")+
     theme_classic()+
     scale_fill_manual(values=c("#E18727FF","#0072B5FF","#BC3C29FF"))+
     facet_wrap(~anno,ncol=2,scales="free") + 
     theme(axis.text.x=element_text(colour="black",angle=45,hjust=1,vjust=1,size=16),
           axis.text.y=element_text(colour="black",size=16),
           axis.title.x=element_blank(),
           axis.title.y=element_text(colour="black",size=16),
           plot.title=element_blank(),
           strip.background=element_blank(),
           strip.text=element_text(colour="black",size=16),
           panel.background =element_rect(fill="transparent",colour="black",size=1),
           panel.grid = element_blank(),
      　　　legend.position="right",
           legend.title=element_text(colour="black",size=16),
           legend.text=element_text(colour="black",size=16))+
     scale_x_discrete(labels= label)+
　　　labs(y="Number of significant\npathways")

    pdf_3(paste0("Figure/",today(),"_ORA_HALLMARKKEGG_sig_sum_wyaxislegend.pdf"),h=4,w=18)
     plot(p)
    dev.off()

############# Supple
# HALLMARK
list1=make_list("ORA/HALLMARK/state","_ORA_HALLMARK.txt")
list1$subset=take_factor(list1$FILE,2:3,"_")
for(iii in 1:nrow(list1)){
    subset_tmp=list1$subset[iii]
    res_tmp=fread_FT(list1$PATH[iii])
    res_tmp$subset=subset_tmp
    if(iii==1){res_sum1=res_tmp}else{res_sum1=rbind(res_sum1,res_tmp)}
}
res_sum1$DEGtype="Disease-state"

list2=make_list("ORA/HALLMARK/activity","_ORA_HALLMARK.txt")
list2$subset=take_factor(list2$FILE,2:3,"_")
for(iii in 1:nrow(list2)){
    subset_tmp=list2$subset[iii]
    res_tmp=fread_FT(list2$PATH[iii])
    res_tmp$subset=subset_tmp
    if(iii==1){res_sum2=res_tmp}else{res_sum2=rbind(res_sum2,res_tmp)}
}
res_sum2$DEGtype="Disease-activity"

res_sum_all1=rbind(res_sum1,res_sum2)

res_sum_all1$subset=factor(res_sum_all1$subset,levels=celltype_reordered_27)
res_sum_all1$DEGtype=factor(res_sum_all1$DEGtype,levels=c("Disease-state","Disease-activity"))
res_sum_all1=res_sum_all1[order(res_sum_all1$subset),]
res_sum_all1=res_sum_all1[order(res_sum_all1$DEGtype),]

res_sum_lim1=res_sum_all1%>%left_join(.,celltype_corresp[,c(1,2)],by="subset")%>%
            select(DEGtype,label,Description,pvalue,p.adjust)%>%filter(pvalue<0.05)
res_sum_lim1$pvalue = formatC(res_sum_lim1$pvalue,digits=2)
res_sum_lim1$p.adjust  = formatC(res_sum_lim1$p.adjust,digits=2)
colnames(res_sum_lim1)=c("Signature","Cell type","Annotation","Pvalue","FDR")

write.table_FT_2(res_sum_lim1,paste0(today(),"_ORA_HALLMARK_27subsets_stateactivity_Supple.txt"))

# KEGG
list1=make_list("ORA/KEGG/state","_ORA_KEGG.txt")
list1$subset=take_factor(list1$FILE,2:3,"_")
for(iii in 1:nrow(list1)){
    subset_tmp=list1$subset[iii]
    res_tmp=fread_FT(list1$PATH[iii])
    res_tmp$subset=subset_tmp
    if(iii==1){res_sum1=res_tmp}else{res_sum1=rbind(res_sum1,res_tmp)}
}
res_sum1$DEGtype="Disease-state"

list2=make_list("ORA/KEGG/activity","_ORA_KEGG.txt")
list2$subset=take_factor(list2$FILE,2:3,"_")
for(iii in 1:nrow(list2)){
    subset_tmp=list2$subset[iii]
    res_tmp=fread_FT(list2$PATH[iii])
    res_tmp$subset=subset_tmp
    if(iii==1){res_sum2=res_tmp}else{res_sum2=rbind(res_sum2,res_tmp)}
}
res_sum2$DEGtype="Disease-activity"


res_sum_all2=rbind(res_sum1,res_sum2)
res_sum_all2$subset=factor(res_sum_all2$subset,levels=celltype_reordered_27)
res_sum_all2$DEGtype=factor(res_sum_all2$DEGtype,levels=c("Disease-state","Disease-activity"))
res_sum_all2=res_sum_all2[order(res_sum_all2$subset),]
res_sum_all2=res_sum_all2[order(res_sum_all2$DEGtype),]

res_sum_lim2=res_sum_all2%>%left_join(.,celltype_corresp[,c(1,2)],by="subset")%>%
            select(DEGtype,label,Description,pvalue,p.adjust)%>%filter(pvalue<0.05)
res_sum_lim2$pvalue = formatC(res_sum_lim2$pvalue,digits=2)
res_sum_lim2$p.adjust  = formatC(res_sum_lim2$p.adjust,digits=2)
colnames(res_sum_lim2)=c("Signature","Cell type","Annotation","Pvalue","FDR")
write.table_FT_2(res_sum_lim2,paste0(today(),"_ORA_KEGG_27subsets_stateactivity_Supple.txt"))

######################################################################################################################
### HALLMARK representative pathways
list_tmp = make_list("ORA/HALLMARK",".txt")
list_tmp$subset = take_factor(list_tmp$FILE,2:3,"_")
list_tmp$type   = take_factor(list_tmp$FILE,4,"_")
list_tmp$anno   = take_factor(list_tmp$FILE,6,"_") %>% gsub(".txt","",.) 
list_tmp=list_tmp %>% filter(type%in%c("state","activity"))

focus=c("HALLMARK_INTERFERON_ALPHA_RESPONSE","HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_COMPLEMENT","HALLMARK_MTORC1_SIGNALING")

 for(kkk in 1:nrow(list_tmp)){
   subset_tmp=list_tmp$subset[kkk]
   type_tmp  =list_tmp$type[kkk]
   res_tmp   = fread_FT(list_tmp$PATH[kkk]) %>% filter(Description%in%focus) %>% select(Description,pvalue)
   if(nrow(res_tmp)==0){res_tmp=data.frame(Description=focus,pvalue=1)}else{res_tmp=res_tmp} 
   res_tmp$subset = subset_tmp
   res_tmp$type   = type_tmp
   if(kkk==1){res_sum1=res_tmp}else{res_sum1=union(res_sum1,res_tmp)}
   } 
   res_sum1$subset=factor(res_sum1$subset,levels=celltype_reordered_27)
   res_sum1$type  =factor(res_sum1$type,levels=c("state","activity"))
   res_sum1$Description  =factor(res_sum1$Description,levels=focus)

   p = ggplot()+
     geom_bar(data=res_sum1, aes(x=subset,y=-log10(pvalue),fill=type),stat="identity",position="dodge")+
     theme_classic()+
     scale_fill_manual(values=c("#0072B5FF","#BC3C29FF"))+
     facet_wrap(~Description,ncol=2,scales="free") + 
     theme(axis.text.x=element_text(colour="black",angle=45,hjust=1,vjust=1,size=16),
           axis.text.y=element_text(colour="black",size=16),
           axis.title.x=element_blank(),
           axis.title.y=element_text(colour="black",size=16),
           strip.background=element_blank(),
           strip.text=element_text(colour="black",size=16),
           panel.background =element_rect(fill="transparent",colour="black",size=1),
           panel.grid = element_blank(),
           plot.title=element_blank(),
      　　　legend.position="none")+
     scale_x_discrete(labels= label)+
　　　labs(y="")

    pdf_3(paste0("ORA/",today(),"_ORA_HALLMARK_rep_logP.pdf"),h=7,w=16)
     plot(p)
    dev.off()

######################################################################################################################
### KEGG annotation class
suppressPackageStartupMessages(library(clusterProfiler))
KEGG_list_new=download_KEGG(species="hsa",keggType="KEGG",keyType="kegg")$KEGGPATHID2NAME
colnames(KEGG_list_new)=c("id","pathway")
write.table_FT_2(KEGG_list_new,paste0("KEGG_list/",today(),"_clusterProfiler_KEGG_list.txt"))
dim(KEGG_list_new)
# [1] 548   2

list_tmp = make_list("ORA/KEGG",".txt")
list_tmp$subset = take_factor(list_tmp$FILE,2:3,"_")
list_tmp$type   = take_factor(list_tmp$FILE,4,"_")
list_tmp$anno   = take_factor(list_tmp$FILE,6,"_") %>% gsub(".txt","",.) 
list_tmp=list_tmp %>% filter(type%in%c("state","activity"))

for(kkk in 1:nrow(list_tmp)){
   sig_tmp   = fread_FT(list_tmp$PATH[kkk]) %>% filter(p.adjust<0.05)
   sig_tmp   = sig_tmp$Description
   if(kkk==1){sig_sum=sig_tmp}else{sig_sum=union(sig_sum,sig_tmp)}
   } 
sig_sum=data.frame(pathway=sig_sum)

KEGG_list_new_lim=right_join(KEGG_list_new,sig_sum,by=c("pathway"))

write.table_FT_2(KEGG_list_new_lim,paste0("ORA/",today(),"_ORA_KEGG_sig_pathwaylist.txt"))
# Class info was manually added 

KEGG_list_class=fread_FT("ORA/211204_ORA_KEGG_sig_pathwaylist_wclass.txt")

KEGG_list_class_2=table_freq(KEGG_list_class$Group)

KEGG_list_class_2=KEGG_list_class_2%>%select(factor)
colnames(KEGG_list_class_2)="Group"
write.table_FT_2(KEGG_list_class_2,paste0("ORA/",today(),"_ORA_KEGG_sig_classlist.txt"))
# Add large class(Group2)
KEGG_list_class_2=fread_FT("ORA/211206_ORA_KEGG_sig_classlist_2.txt")

# Annotation count group1 and group2
for(kkk in 1:nrow(list_tmp)){
   subset_tmp=list_tmp$subset[kkk]
   type_tmp  =list_tmp$type[kkk]
   sig_tmp   = fread_FT(list_tmp$PATH[kkk]) %>% filter(p.adjust<0.05)
   sig_tmp2   = sig_tmp%>%select(Description,pvalue,p.adjust) %>%
                          left_join(.,KEGG_list_class,by=c("Description"="pathway")) %>%
                          left_join(.,KEGG_list_class_2,by="Group")
   sig_tmp2$subset=subset_tmp
   sig_tmp2$type=type_tmp

   table_tmp=table_freq(sig_tmp2$Group)
   table_tmp$subset=subset_tmp
   table_tmp$type=type_tmp

   table_tmp2=table_freq(sig_tmp2$Group2)
   table_tmp2$subset=subset_tmp
   table_tmp2$type=type_tmp

   if(kkk==1){sig_sum=sig_tmp2}else{sig_sum=union(sig_sum,sig_tmp2)}
   if(kkk==1){table_sum=table_tmp}else{table_sum=union(table_sum,table_tmp)}
   if(kkk==1){table_sum2=table_tmp2}else{table_sum2=union(table_sum2,table_tmp2)}
   } 
write.table_FT_2(sig_sum,paste0("ORA/",today(),"_ORA_KEGG_sig_allres_wclass.txt"))
write.table_FT_2(table_sum,paste0("ORA/",today(),"_ORA_KEGG_sig_allres_wGroup_count.txt"))
write.table_FT_2(table_sum2,paste0("ORA/",today(),"_ORA_KEGG_sig_allres_wGroup2_count.txt"))


###########
# pathways class summerize: pie chart
all_table_sum=table_freq(sig_sum$Group) %>% left_join(.,KEGG_list_class_2,by=c("factor"="Group"))
colnames(all_table_sum)=c("Subcategory","Freq","Category")
all_table_sum$Category=factor(all_table_sum$Category,levels=c("Immune disease/system","Non-immune disease/system","Cellular process/Cancer","Metabolism"))
all_table_sum=all_table_sum[order(all_table_sum$Freq,decreasing=T),]
all_table_sum=all_table_sum[order(all_table_sum$Category),]
write.table_FT_2(all_table_sum,paste0("ORA/",today(),"_ORA_KEGG_sig_allres_wGroup_allcount.txt"))

all_table_sum$Subcategory=factor(all_table_sum$Subcategory,levels=unique(all_table_sum$Subcategory))

colors=c(rev(brewer.pal(7,"Greens")),
         rev(colorRampPalette(brewer.pal(9,"Oranges"))(11)),
         rev(colorRampPalette(brewer.pal(9,"Purples"))(10)),
         rev(brewer.pal(6,"Reds")))

p = ggplot()+
     geom_bar(data=all_table_sum, aes(x="",y=Freq,fill=Subcategory),stat="identity")+
     coord_polar("y",start=0,direction=-1) +
     theme_void()+
     scale_fill_manual(values=colors)

pdf_3(paste0("ORA/",today(),"_ORA_KEGG_sig_allres_wGroup_allcount.pdf"),h=10,w=10)
 plot(p)
dev.off()

nrow(sig_sum)
#[1] 774

table_freq(sig_sum$Group2)
#                     factor Freq
#1   Cellular process/Cancer  169
#2     Immune disease/system  370
#3                Metabolism   36
#4 Non-immune disease/system  199

###########
######## 3categories stacked barplot
list_tmp = make_list("ORA/KEGG",".txt")
list_tmp$subset = take_factor(list_tmp$FILE,2:3,"_")
list_tmp$type   = take_factor(list_tmp$FILE,4,"_")
list_tmp$anno   = take_factor(list_tmp$FILE,6,"_") %>% gsub(".txt","",.) 
list_tmp=list_tmp %>% filter(type%in%c("state","activity"))

for(kkk in 1:length(celltype_reordered_27)){
   subset_tmp=celltype_reordered_27[kkk]
   list_tmp_2=list_tmp%>%filter(subset==subset_tmp)
   state_tmp=fread_FT(list_tmp_2$PATH[2])%>%select(Description,pvalue,p.adjust)
   colnames(state_tmp)[2:3]=paste0("state_",colnames(state_tmp)[2:3])
   act_tmp=fread_FT(list_tmp_2$PATH[1])%>%select(Description,pvalue,p.adjust)
   colnames(act_tmp)[2:3]=paste0("activity_",colnames(act_tmp)[2:3])
   res_tmp=left_join(state_tmp,act_tmp,by="Description")%>%
           mutate(siggroup005=ifelse(state_p.adjust<0.05&activity_p.adjust<0.05,"Both significant",
                              ifelse(state_p.adjust<0.05,"Disease-state only",
                              ifelse(activity_p.adjust<0.05,"Disease-activity only","Not significant"))))%>%
                          left_join(.,KEGG_list_class,by=c("Description"="pathway")) %>%
                          left_join(.,KEGG_list_class_2,by="Group")
   table_tmp=table(res_tmp$siggroup005,res_tmp$Group2) %>%as.data.frame() %>% mutate(subset=subset_tmp)
   colnames(table_tmp)[1:2]=c("Enrichment","Group2")
   if(kkk==1){table_sum=table_tmp}else{table_sum=rbind(table_sum,table_tmp)}
   } 
write.table_FT_2(table_sum,paste0("ORA/KEGG/",today(),"_ORA_KEGG_sig_sum_wbothinfo_wGroup.txt"))

table_sum=table_sum%>%filter(!Enrichment=="Not significant")
table_sum$subset  =factor(table_sum$subset,levels=celltype_reordered_27)
table_sum$Enrichment = factor(table_sum$Enrichment,levels=c("Both significant","Disease-state only","Disease-activity only"))

# Cellular process/ cancer
table_sum1=table_sum%>%filter(Group2=="Cellular process/Cancer")
 p = ggplot()+
     geom_bar(data=table_sum1, aes(x=subset,y=Freq,fill=Enrichment),stat="identity",position="stack")+
     theme_classic()+
     scale_fill_manual(values=c("#E18727FF","#0072B5FF","#BC3C29FF"))+
     #facet_wrap(~factor,ncol=1,scales="free") + 
     theme(axis.text.x=element_text(colour="black",angle=45,hjust=1,vjust=1,size=16),
           axis.text.y=element_text(colour="black",size=16),
           axis.title.x=element_blank(),
           axis.title.y=element_text(colour="black",size=16),
           axis.ticks.x = element_blank(),
           plot.title=element_blank(),
      　　　legend.position="none")+
     scale_x_discrete(labels= label)+
     scale_y_continuous(breaks=c(0,5,10))+
　　　labs(y="")

    pdf_3(paste0("Figure/",today(),"_ORA_KEGG_sig_allres_cellular_count.pdf"),h=3.5,w=7.5)
     plot(p)
    dev.off()

# Metabolism
table_sum2=table_sum%>%filter(Group2=="Metabolism")
 p = ggplot()+
     geom_bar(data=table_sum2, aes(x=subset,y=Freq,fill=Enrichment),stat="identity",position="stack")+
     theme_classic()+
     scale_fill_manual(values=c("#E18727FF","#0072B5FF","#BC3C29FF"))+
     #facet_wrap(~factor,ncol=1,scales="free") + 
     theme(axis.text.x=element_text(colour="black",angle=45,hjust=1,vjust=1,size=16),
           axis.text.y=element_text(colour="black",size=16),
           axis.title.x=element_blank(),
           axis.title.y=element_text(colour="black",size=16),
           axis.ticks.x = element_blank(),
           plot.title=element_blank(),
      　　　legend.position="none")+
     scale_x_discrete(labels= label)+
     scale_y_continuous(breaks=c(0,3,6))+
　　　labs(y="")

    pdf_3(paste0("Figure/",today(),"_ORA_KEGG_sig_allres_metab_count.pdf"),h=3.5,w=7.5)
     plot(p)
    dev.off()

###########
# Representaive log P barplot
focus=c("Cell cycle","Ribosome","Oxidative phosphorylation","Citrate cycle (TCA cycle)")

# Figure
for(iii in 1:length(focus)){
   pathway_tmp=focus[iii]

 for(kkk in 1:nrow(list_tmp)){
   subset_tmp=list_tmp$subset[kkk]
   type_tmp  =list_tmp$type[kkk]
   res_tmp   = fread_FT(list_tmp$PATH[kkk]) %>% filter(Description==pathway_tmp) %>% select(Description,pvalue)
   if(nrow(res_tmp)==0){res_tmp=data.frame(Description=pathway_tmp,pvalue=1)}else{res_tmp=res_tmp}
   res_tmp$subset = subset_tmp
   res_tmp$type   = type_tmp
   if(kkk==1){res_sum1=res_tmp}else{res_sum1=union(res_sum1,res_tmp)}
   } 
   res_sum1$subset=factor(res_sum1$subset,levels=celltype_reordered_27)
   res_sum1$type  =factor(res_sum1$type,levels=c("state","activity"))

   p = ggplot()+
     geom_bar(data=res_sum1, aes(x=subset,y=-log10(pvalue),fill=type),stat="identity",position="dodge")+
     theme_classic()+
     scale_fill_manual(values=c("#0072B5FF","#BC3C29FF"))+
     #facet_wrap(~factor,ncol=1,scales="free") + 
     theme(axis.text.x=element_blank(),
           axis.text.y=element_text(colour="black",size=16),
           axis.title.x=element_blank(),
           axis.title.y=element_blank(),
           axis.ticks.x = element_blank(),
           plot.title=element_blank(),
      　　　legend.position="none")+
     scale_x_discrete(labels= label)

    pdf_3(paste0("Figure/",today(),"_ORA_KEGG_",pathway_tmp,"_logP.pdf"),h=1,w=7.5)
     plot(p)
    dev.off()
 }

######################################################################################################################
# Cytokine DEGs heatmap
######################################################################################################################

res_sum = fread_FT("res/stateactivity/211201_COI_27subsets_state_activity_logFC_logCPM_P_FDR.txt")
res_sum$subset   = factor(res_sum$subset,levels=celltype_reordered_27)
res_sum$siggroup005   = factor(res_sum$siggroup005,levels=c("only_state","only_activity","significant_in_both","discordant","neither"))

chemo_L = read.table_FT("gene_set/210606_chemokine_ligand.txt")  %>%filter(!genes=="CCL4L2")
IL_L    = read.table_FT("gene_set/210606_interleukin_ligand.txt")%>%filter(!genes%in%c("IL21-AS1","IL12-AS1"))
IFN_L   = read.table_FT("gene_set/210606_interferon_ligand.txt") %>%filter(!genes=="IFNG-AS1")
TNF_L   = read.table_FT("gene_set/210606_TNF_ligand.txt")
cytokines = c(IFN_L$genes,IL_L$genes,chemo_L$genes,TNF_L$genes) 
# [1] 136

res_sum_ligand = res_sum %>% filter(genes%in%cytokines)
res_sum_ligand　= res_sum_ligand[order(res_sum_ligand$subset),]

# Only upregulated genes
res_sum_ligand_state_up = res_sum_ligand %>% filter(state_FDR<0.05) %>% filter(state_logFC>0)
dim(res_sum_ligand_state_up)
# [1] 245   9
length(unique(res_sum_ligand_state_up$genes))
# [1] 48

res_sum_ligand_activity_up = res_sum_ligand %>% filter(activity_FDR<0.05) %>% filter(activity_logFC>0)
dim(res_sum_ligand_activity_up)
# [1] 152   9
length(unique(res_sum_ligand_activity_up$genes))
# [1] 37

up_ligand_genes=union(unique(res_sum_ligand_state_up$genes),unique(res_sum_ligand_activity_up$genes))
length(up_ligand_genes)
# [1] 51

res_sum_ligand_up=res_sum_ligand %>% filter(genes%in%up_ligand_genes) %>% 
                  mutate(Up_DEGs=ifelse(siggroup005=="only_state"&state_logFC>0,"only_state_up",
                                    ifelse(siggroup005=="only_activity"&activity_logFC>0,"only_activity_up",
                                    ifelse(siggroup005=="significant_in_both"&activity_logFC>0,"both_up",
                                    ifelse(siggroup005=="discordant"&state_logFC>0,"only_state_up",
                                    ifelse(siggroup005=="discordant"&activity_logFC>0,"only_activity_up","not_up")))))) %>%
                  mutate(state_Z=ifelse(state_logFC<0,-qnorm(state_PValue/2,lower.tail=FALSE),
                                                 qnorm(state_PValue/2,lower.tail=FALSE))) %>%
                  mutate(activity_Z=ifelse(activity_logFC<0,-qnorm(activity_PValue/2,lower.tail=FALSE),
                                                 qnorm(activity_PValue/2,lower.tail=FALSE))) %>%
                  mutate(subset_genes=paste0(subset,"_",genes))


### logCPM scaling 
res_sum_ligand_logCPM=res_sum_ligand_up%>%select(subset,genes,state_logCPM)%>%
                                          pivot_wider(names_from="subset", values_from="state_logCPM")%>%
                                          column_to_rownames("genes")
min=min(res_sum_ligand_logCPM,na.rm=T)
res_sum_ligand_logCPM[is.na(res_sum_ligand_logCPM)]=min  # LowExp filtered genes are regarded as the minimum value of the matrix                                         
scaled_logCPM=apply(res_sum_ligand_logCPM,1,scale) %>% t() %>% as.data.frame()
colnames(scaled_logCPM)=colnames(res_sum_ligand_logCPM)

#apply(scaled_logCPM,1,sum)%>%max()
#apply(scaled_logCPM,1,sd)%>%max()
#apply(scaled_logCPM,1,sd)%>%min()

scaled_logCPM_2=scaled_logCPM%>%rownames_to_column("genes")%>%
                                pivot_longer(cols=-genes,names_to="subset", values_to='scaled_exp')%>%
                                mutate(subset_genes=paste0(subset,"_",genes)) %>%
                                select(-c(genes,subset))


### bind to original data and save
res_sum_ligand_up_2=res_sum_ligand_up%>%left_join(.,scaled_logCPM_2,by="subset_genes") %>%
                                        select(-subset_genes)

write.table_FT_2(res_sum_ligand_up_2,paste0("cytokine/",today(),"_cytokine_upDEG_sum.txt"))
table(res_sum_ligand_up_2$Up_DEGs)
# both_up           not_up only_activity_up    only_state_up 
#      91              376               61              154 

### Draw figure
# res_sum_ligand_up_2=fread_FT("cytokine/211205_cytokine_upDEG_sum.txt")

# Clustering based on activity Z score
res_sum_ligand_Z=res_sum_ligand_up_2%>%select(subset,genes,activity_Z)%>%
                                       pivot_wider(names_from="subset", values_from="activity_Z") %>%
                                       column_to_rownames("genes")
res_sum_ligand_Z[is.na(res_sum_ligand_Z)] = 0

hc_row = hclust(dist(res_sum_ligand_Z), method="ward.D2")
pdf_3(paste0("cytokine/",today(),"_cytokine_upDEG_hrow.pdf"),h=5,w=10)
 plot(hc_row)
dev.off()

hc_col = hclust(dist(t(res_sum_ligand_Z)), method="ward.D2")
pdf_3(paste0("cytokine/",today(),"_cytokine_upDEG_hcol.pdf"),h=5,w=10)
 plot(hc_col)
dev.off()

label_row=hc_row$labels[hc_row$order]
label_col=hc_col$labels[hc_col$order]
label_row = label_row[c(1:8,10:11,9,12:33,37:34,38:47,49:51,48)]   
label_col = label_col[c(1:6,12:13,16,20:22,19:17,14:15,27:23,7:11)]   

# reorder
res_sum_ligand_up_2$Up_DEGs=factor(res_sum_ligand_up_2$Up_DEGs,levels=c("not_up","only_state_up","only_activity_up","both_up"))
levels(res_sum_ligand_up_2$Up_DEGs)=list(`Not significant`="not_up",`Disease-state only`="only_state_up",`Disease-activity only`="only_activity_up",`Both significant`="both_up")
res_sum_ligand_up_2$subset=factor(res_sum_ligand_up_2$subset,levels=label_col)
res_sum_ligand_up_2$genes =factor(res_sum_ligand_up_2$genes,levels=label_row)

### Draw heatmap
p=ggplot(res_sum_ligand_up_2,aes(x=subset,y=genes,color=Up_DEGs,fill=Up_DEGs,size=scaled_exp))+
 geom_point(shape=19,stroke=0) +
 scale_radius(range=c(0.6,6)) +
 theme_bw()+
 scale_color_manual(values=c("#cccccc","#0072B5FF","#BC3C29FF","#E18727FF"))+
 scale_fill_manual(values=c("#cccccc","#0072B5FF","#BC3C29FF","#E18727FF"))+
 theme(axis.text.x=element_text(colour="black",angle=90,hjust=1,vjust=1,size=16),
      axis.text.y=element_text(colour="black",face="italic",size=13),
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title=element_blank(),
      panel.grid.minor = element_blank(),
      panel.background =element_rect(fill="transparent",colour="black",size=1),
 　　　legend.position="right",
      legend.title=element_text(colour="black",size=16),
      legend.text=element_text(colour="black",size=16))+
guides(colour=guide_legend(title="Upregulation",override.aes=list(size=4.5)))+
scale_x_discrete(labels= label)+
scale_y_discrete(position="right")

pdf_3(paste0("cytokine/",today(),"_cytokine_upDEG_sum.pdf"),h=9,w=9.5)
 plot(p)
dev.off()

#############################
# Count lineage-specific cytokines 
# res_sum_ligand_up_2=fread_FT("cytokine/211205_cytokine_upDEG_sum.txt")

res_sum_ligand_up_3 = res_sum_ligand_up_2 %>% left_join(.,celltype_corresp%>%select(subset,lineage),by="subset")

for(iii in 1:length(unique(res_sum_ligand_up_3$genes))){
    gene_tmp = unique(res_sum_ligand_up_3$genes)[iii]
    data_tmp = res_sum_ligand_up_3 %>% filter(genes==gene_tmp)
    state = data_tmp %>% filter(Up_DEGs%in%c("only_state_up","both_up"))
    activity = data_tmp %>% filter(Up_DEGs%in%c("only_activity_up","both_up"))
    state_lineage = length(unique(state$lineage))
    activity_lineage = length(unique(activity$lineage))
    res_tmp = data.frame(gene=gene_tmp,state=state_lineage,activity=activity_lineage)
    if(iii==1){res_sum=res_tmp}else{res_sum=rbind(res_sum,res_tmp)}
}

table_freq(res_sum$state)
#   factor Freq
# 1      0    3
# 2      1   21
# 3      2    8
# 4      3    9
# 5      4    6
# 6      5    1
# 7      6    1
# 8      7    2

table_freq(res_sum$activity)
#   factor Freq
# 1      0   14
# 2      1   17
# 3      2   11
# 4      3    6
# 5      5    1
# 6      6    1
# 7      7    1

#############################
# Plot IL21 and CXCL13 
# IL21 in Tfh was filtered at lowExp filter
# Need no filtered logCPM
suppressPackageStartupMessages(library(edgeR))
list         = make_list("~/ws/2021/211124_COI_SLE_HC_run27_finaldata_Combat_PCA/data/SLEunique_andHC/each_subset/each_subset_count","_count.txt")
list$subset  = take_factor(list$FILE,5:6,"_")
list=list%>%filter(subset%in%c("B01_Th1","B04_Tfh"))

clinical= fread_FT("data_ref/211124_SLEunique_andHC_COI_clinicaldata_final_lim.txt")%>%select(id,Activity)
col3=c("HC"="#F8766D","Inactive"="#A3A500","HDA"="#E76BF3")

# wo LowExp
for(iii in 1:nrow(list)){
  subset_tmp = list$subset[iii]
  count_tmp = fread_n(list$PATH[iii])
  count_tmp_2  = count_tmp
  dge     = DGEList(counts=count_tmp_2, genes=rownames(count_tmp_2))
  dge_n   = calcNormFactors( dge ,method="TMM")
  logCPM2  = log2(cpm(dge_n)+1) %>% as.data.frame()
  write.table_n_2(logCPM2,"Gene",paste0("cytokine/logCPM_wofilter/",today(),"_COI_SLEuniqueHC_",subset_tmp,"_TMM_logCPM.txt"))
}

# Plot CXCL13, IL21
list         = make_list("cytokine/logCPM_wofilter/","_TMM_logCPM.txt")
list$subset  = take_factor(list$FILE,4:5,"_")

clinical= fread_FT("data_ref/211124_SLEunique_andHC_COI_clinicaldata_final_lim.txt")%>%select(id,Activity)
col3=c("HC"="#F8766D","Inactive"="#A3A500","HDA"="#E76BF3")

for(iii in 1:nrow(list)){
    subset_tmp=list$subset[iii]
    data_tmp=fread_n(list$PATH[iii])["CXCL13",] %>% t() %>% as.data.frame()
    data_tmp2 = data_tmp %>% rownames_to_column("name") %>% mutate(id=take_factor(name,1,"_")) %>% select(-name) %>%
                            left_join(.,clinical,by="id") %>%filter(Activity%in%c("0HC","1Inactive","4HDA")) %>%
                            mutate(subset=subset_tmp)
    if(iii==1){data_sum=data_tmp2}else{data_sum=union(data_sum,data_tmp2)}
    }
 
data_sum$Activity=factor(data_sum$Activity,levels=c("0HC","1Inactive","4HDA"))
levels(data_sum$Activity)=list(`HC`="0HC",`Inactive`="1Inactive",`HDA`="4HDA")
data_sum$subset=factor(data_sum$subset,levels=c("B01_Th1","B04_Tfh"))

p =  ggplot(data_sum,aes(x=Activity,y=CXCL13,fill=Activity))+
 geom_boxplot(outlier.shape=NA)+
 #geom_point(size=0.5,color="black",position=position_jitterdodge())+
 scale_fill_manual(values=col3)+
 facet_wrap(~subset,ncol=2) + 
 theme_classic()+
 theme(axis.text.x=element_text(colour="black",angle=45,hjust=1,vjust=1,size=14),
       axis.text.y=element_text(colour="black",size=14),
       axis.title.x=element_blank(),
       axis.title.y=element_text(colour="black",size=14),
       plot.title=element_blank(),
       strip.background=element_blank(),
       strip.text=element_blank(),
       panel.background =element_rect(fill="transparent",colour="black",size=1),
       panel.grid = element_blank(),
       legend.position="none")+
 expand_limits(y=0)+
 labs(y="log(CPM+1)")
 pdf_3(paste0("cytokine/",today(),"_COI_27subset_repDEGs_exp_Th1Tfh_CXCL13.pdf"),h=2.5,w=3)
  plot(p)
 dev.off()

for(iii in 1:nrow(list)){
    subset_tmp=list$subset[iii]
    data_tmp=fread_n(list$PATH[iii])["IL21",] %>% t() %>% as.data.frame()
    data_tmp2 = data_tmp %>% rownames_to_column("name") %>% mutate(id=take_factor(name,1,"_")) %>% select(-name) %>%
                            left_join(.,clinical,by="id") %>%filter(Activity%in%c("0HC","1Inactive","4HDA")) %>%
                            mutate(subset=subset_tmp)
    if(iii==1){data_sum=data_tmp2}else{data_sum=union(data_sum,data_tmp2)}
    }
 
data_sum$Activity=factor(data_sum$Activity,levels=c("0HC","1Inactive","4HDA"))
levels(data_sum$Activity)=list(`HC`="0HC",`Inactive`="1Inactive",`HDA`="4HDA")
data_sum$subset=factor(data_sum$subset,levels=c("B01_Th1","B04_Tfh"))

p =  ggplot(data_sum,aes(x=Activity,y=IL21,fill=Activity))+
 geom_boxplot(outlier.shape=NA)+
 #geom_point(size=0.5,color="black",position=position_jitterdodge())+
 scale_fill_manual(values=col3)+
 facet_wrap(~subset,ncol=2) + 
 theme_classic()+
 theme(axis.text.x=element_text(colour="black",angle=45,hjust=1,vjust=1,size=14),
       axis.text.y=element_text(colour="black",size=14),
       axis.title.x=element_blank(),
       axis.title.y=element_text(colour="black",size=14),
       plot.title=element_blank(),
       strip.background=element_blank(),
       strip.text=element_blank(),
       panel.background =element_rect(fill="transparent",colour="black",size=1),
       panel.grid = element_blank(),
       legend.position="none")+
 expand_limits(y=0)+
 labs(y="log(CPM+1)")
 pdf_3(paste0("cytokine/",today(),"_COI_27subset_repDEGs_exp_Th1Tfh_IL21.pdf"),h=2.5,w=3)
  plot(p)
 dev.off()

###################################################################################################
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
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] clusterProfiler_3.18.1

loaded via a namespace (and not attached):
 [1] Biobase_2.50.0       viridis_0.6.1        tidyr_1.1.3         
 [4] tidygraph_1.2.0      bit64_4.0.5          viridisLite_0.4.0   
 [7] splines_4.0.2        ggraph_2.0.5         assertthat_0.2.1    
[10] shadowtext_0.0.8     DO.db_2.9            BiocManager_1.30.16 
[13] rvcheck_0.1.8        stats4_4.0.2         blob_1.2.1          
[16] ggrepel_0.9.1        pillar_1.6.1         RSQLite_2.2.7       
[19] lattice_0.20-41      glue_1.4.2           downloader_0.4      
[22] digest_0.6.27        RColorBrewer_1.1-2   polyclip_1.10-0     
[25] qvalue_2.22.0        colorspace_2.0-2     cowplot_1.1.1       
[28] Matrix_1.3-4         plyr_1.8.6           pkgconfig_2.0.3     
[31] purrr_0.3.4          GO.db_3.12.1         scales_1.1.1        
[34] tweenr_1.0.2         enrichplot_1.10.2    BiocParallel_1.24.1 
[37] ggforce_0.3.3        tibble_3.1.2         generics_0.1.0      
[40] farver_2.1.0         IRanges_2.24.1       ggplot2_3.3.5       
[43] ellipsis_0.3.2       cachem_1.0.5         BiocGenerics_0.36.1 
[46] magrittr_2.0.1       crayon_1.4.1         memoise_2.0.0       
[49] DOSE_3.16.0          fansi_0.5.0          MASS_7.3-54         
[52] tools_4.0.2          data.table_1.14.0    lifecycle_1.0.0     
[55] stringr_1.4.0        S4Vectors_0.28.1     munsell_0.5.0       
[58] AnnotationDbi_1.52.0 compiler_4.0.2       rlang_0.4.11        
[61] grid_4.0.2           igraph_1.2.6         gtable_0.3.0        
[64] DBI_1.1.1            reshape2_1.4.4       graphlayouts_0.7.1  
[67] R6_2.5.0             gridExtra_2.3        dplyr_1.0.7         
[70] fastmap_1.1.0        bit_4.0.4            utf8_1.2.1          
[73] fastmatch_1.1-0      fgsea_1.16.0         GOSemSim_2.16.1     
[76] stringi_1.6.2        parallel_4.0.2       Rcpp_1.0.6          
[79] vctrs_0.3.8          tidyselect_1.1.1     scatterpie_0.1.6   




























