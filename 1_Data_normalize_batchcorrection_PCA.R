#---------------------------------------#
# 211124 Nakano
# QC -> LowExp -> TMM -> logCPM -> Combat
#---------------------------------------#

##############################################
source("my.source.R")
options(stringsAsFactors=F)
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(sva))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ggsci))
suppressPackageStartupMessages(library(uwot))

##############################################
# Count and Metadata
##############################################

# condition, count data (open at NBDC)
condition=fread_FT("data/SLEinc_andHC/211109_COI_SLEinc_andHC_27subsets_cond.txt")
dim(condition)
# [1] 6386   15
length(unique(condition$id))
# [1] 248

count = fread_n("data/SLEinc_andHC/211109_COI_SLEinc_andHC_27subsets_count.txt")
dim(count)
# [1] 26353  6386

# Clinical info (Private)
clinical_lim = fread_FT("data_ref/clinicaldata_lim.txt") 
clinical_lim_unique = fread_FT("data_ref/clinicaldata_unique_lim.txt") 

condition_2 = condition %>% select("name","id","subset_num","subset","run","lot","sequencer","lotseq") %>%
                            left_join(.,clinical_lim,by="id")
write.table_FT_2(condition_2,paste0("data/SLEinc_andHC/",today(),"_COI_SLEinc_andHC_27subsets_cond_wclinical.txt"))

# Limit to unique pts (Discovery cohort)
condition_2_unique = condition_2　%>%　filter(id%in%SLEunique_and_HC$id)
dim(condition_2_unique)
# [1] 5797  212
write.table_FT_2(condition_2_unique,paste0("data/SLEunique_andHC/",today(),"_COI_SLEunique_andHC_27subsets_cond_wclinical.txt"))

count_unique = count[,is.element(colnames(count),condition_2_unique$name)]
dim(count_unique)
# [1] 26353  5797
write.table_n_2(count_unique,"Gene",paste0("data/SLEunique_andHC/",today(),"_COI_SLEunique_andHC_27subsets_count.txt"))

################################################################
# Subsetwise data save 
################################################################

subset_list = unique(condition_2$subset) %>% sort()

# All until Combat
for(iii in 1:length(subset_list)){
    subset_tmp = subset_list[iii]
    cond_tmp   = condition_2 %>% filter(subset==subset_tmp)
    count_tmp  = count[,cond_tmp$name]
    write.table_n_2(count_tmp,"Gene",paste0("data/SLEinc_andHC/each_subset/each_subset_count/",today(),"_COI_SLEinc_andHC_",subset_tmp,"_count.txt"))
    write.table_FT_2(cond_tmp,paste0("data/SLEinc_andHC/each_subset/each_subset_cond/",today(),"_COI_SLEinc_andHC_",subset_tmp,"_cond.txt"))
    cpm_tmp    = cpm(count_tmp)
    filter1    = rownames(count_tmp)[exp_filter(count_tmp,exp=10,p=0.15)]
    filter2    = rownames(cpm_tmp)[exp_filter(cpm_tmp,exp=1,p=0.15)]
    count_tmp_2= count_tmp[intersect(filter1,filter2),]
    dge     = DGEList(counts=count_tmp_2, genes=rownames(count_tmp_2))
    dge_n   = calcNormFactors( dge ,method="TMM")
    logCPM  = log2(cpm(dge_n)+1) %>% as.data.frame()
    write.table_n_2(logCPM,"Gene",paste0("data/SLEinc_andHC/each_subset/LowExp_TMM_logCPM/",today(),"_COI_SLEinc_andHC_",subset_tmp,"_logCPM.txt"))
    model_status   = model.matrix(~disease, data=cond_tmp)
    logCPM_Combat = ComBat(dat=logCPM, batch=cond_tmp$lotseq, mod=model_status, par.prior=TRUE, prior.plots=FALSE) %>% as.data.frame()
    write.table_n_2(logCPM_Combat,"Gene",paste0("data/SLEinc_andHC/each_subset/LowExp_TMM_logCPM_Combat/",today(),"_COI_SLEinc_andHC_",subset_tmp,"_logCPM_Combat.txt"))
}

# unique until raw count
for(iii in 1:length(subset_list)){
    subset_tmp = subset_list[iii]
    cond_tmp   = condition_2_unique %>% filter(subset==subset_tmp)
    count_tmp  = count_unique[,cond_tmp$name]
    write.table_n_2(count_tmp,"Gene",paste0("data/SLEunique_andHC/each_subset/each_subset_count/",today(),"_COI_SLEunique_andHC_",subset_tmp,"_count.txt"))
    write.table_FT_2(cond_tmp,paste0("data/SLEunique_andHC/each_subset/each_subset_cond/",today(),"_COI_SLEunique_andHC_",subset_tmp,"_cond.txt"))
}

#############################################
# celltype-color list 
#############################################
celltype_reordered_27 = c("B05_NCD4","B06_MCD4","B01_Th1","B02_Th2","B03_TH17","B04_Tfh","B09_Fra1","B07_aTreg","B10_Fra3","C01_NCD8","C04_CmCD8","C05_EmCD8","C03_EffectorCD8","G01_NK","A01_NaiB","A03_UnswMB","A02_SwiMB","A04_DNB","A05_PB","E02_CD16nMo","E01_CD16pMo","E04_Intermediate","E03_NonClassical","D01_mDC","D02_pDC","F01_Neu","F02_LDG")

label = c(`B05_NCD4`="Naive CD4",`B06_MCD4`="Mem CD4",`B01_Th1`="Th1",`B02_Th2`="Th2",`B03_TH17`="Th17",`B04_Tfh`="Tfh",
          `B09_Fra1`="Fr. I nTreg",`B07_aTreg`="Fr. II eTreg",`B10_Fra3`="Fr. III T",
          `C01_NCD8`="Naive CD8",`C04_CmCD8`="CM CD8",`C05_EmCD8`="EM CD8",`C03_EffectorCD8`="TEMRA CD8",`G01_NK`="NK",
          `A01_NaiB`="Naive B",`A03_UnswMB`="USM B",`A02_SwiMB`="SM B",`A04_DNB`="DN B",`A05_PB`="Plasmablast",
          `E02_CD16nMo`="CL Mono",`E01_CD16pMo`="CD16p Mono",`E04_Intermediate`="Int Mono",`E03_NonClassical`="NC Mono",
          `D01_mDC`="mDC",`D02_pDC`="pDC",
          `F01_Neu`="Neu",`F02_LDG`="LDG")

subset_color = fread_FT("data_ref/SUBSET_name_for_figures_210414_v5.txt") %>% 
               select(SUBSET2,parents,lineage,COLOR,COLOR_RGB) %>% 
               filter(!is.na(COLOR))

parents_col = c(`CD4`="#cc0010",`CD8`="#8b7355",`NK`="#e0d51a",`B`="#c4f20b",`Monocyte`="#1d2088",`DC`="#a4fbfa",`Neutrophil`="#cccccc")
parents_col = data.frame(parents_col) %>% rownames_to_column("parents")

celltype_corresp=data.frame(label) %>% rownames_to_column("subset") %>% 
                 left_join(.,subset_color,by=c("label"="SUBSET2")) %>%
                 left_join(.,parents_col,by="parents") %>%
                 left_join(.,subset_color[,4:5],by=c("parents_col"="COLOR"))
colnames(celltype_corresp)[3:8]=c("lineage","lineage2","COLOR1","COLOR1_RGB","COLOR2","COLOR2_RGB")

write.table_FT_2(celltype_corresp,paste0("data_ref/COI_27subset_color_list.txt"))

#############################################
# Count sample num that passed QC
#############################################

list=make_list("data/SLEinc_andHC/each_subset/each_subset_cond","_cond.txt")
list$subset = take_factor(list$FILE,5:6,"_")

for(iii in 1:nrow(list)){
  subset_tmp   = list$subset[iii]
  df_1         = read.table_FT(list$PATH[iii]) 
  df_tmp = table_freq(df_1$disease) %>% mutate(subset=subset_tmp)
  if(iii==1){res_sum=df_tmp}else{res_sum=rbind(res_sum,df_tmp)}
  }
res_sum=res_sum%>%left_join(.,celltype_corresp,by="subset")

res_sum$subset = factor(res_sum$subset,levels=celltype_reordered_27)
res_sum$factor = factor(res_sum$factor ,levels=c("0HC","1SLE"))

write.table_FT_2(res_sum,paste0("data/SLEinc_andHC/",today(),"_COI_SLEinc_andHC_samplenum.txt"))

 p = ggplot()+
    geom_bar(data = res_sum, aes(x=subset,y=Freq,fill=factor(factor,labels=c("HC","SLE"))),
           stat="identity",position="stack")+
    theme_classic()+
    #scale_fill_manual(values = col_list)+
    theme(axis.text.x=element_text(colour="black",angle=45,hjust=1,vjust=0.95,size=15),
           axis.text.y=element_text(colour="black",size=15),
           axis.title.x=element_blank(),
           axis.title.y=element_text(colour="black",size=15),
           plot.title=element_blank(),
           legend.position="right",
           legend.title=element_blank(),
           legend.text=element_text(colour="black",size=15))+
     scale_x_discrete(labels= label)+
    labs(y="Number of samples")
    
 pdf_3(paste0("data/SLEinc_andHC/",today(),"_COI_SLEinc_andHC_samplenum.pdf"),h=5,w=10)
  plot(p)
 dev.off()

#############################################
# Count genenum that passed filter
#############################################
# SLEinc+HC AfterCombat for QC
list = make_list("data/SLEinc_andHC/each_subset/LowExp_TMM_logCPM_Combat","_logCPM_Combat.txt",tag="count")
list$subset = take_factor(list$FILE_count,5:6,"_")
list$status = "AfterCombat"

list_tmp=list

 for(iii in 1:nrow(list_tmp)){
  subset_tmp  = list_tmp$subset[iii]
  count_tmp   = fread_n(list_tmp$PATH_count[iii])
  genenum_tmp = data.frame(subset=subset_tmp,genenum=nrow(count_tmp))
  if(iii==1){genenum_sum=genenum_tmp}else{genenum_sum=rbind(genenum_sum,genenum_tmp)}
  if(iii==1){gene_union=rownames(count_tmp)}else{gene_union=union(gene_union,rownames(count_tmp))}
 }
 total_tmp=data.frame(subset="AllUnion",genenum=length(gene_union))
 genenum_sum=rbind(genenum_sum,total_tmp)
 # 16000
 write.table_FT_2(genenum_sum,paste0("data/SLEinc_andHC/",today(),"_COI_SLEinc_andHC_filterpass_genenum.txt"))

 gene_union_sum=data.frame(Gene=sort(gene_union))
 write.table_FT_2(gene_union_sum,paste0("data/SLEinc_andHC/",today(),"_COI_SLEinc_andHC_filterpass_geneunion.txt"))

#################################################################
# PCA for 3 patterns
# 1. Before Combat, SLEinc+HC (All pts)
# 2. After Combat, SLEinc+HC (All pts)
# 3. After Combat, SLEunique+HC; discovery dataset
#################################################################

# SLEinc+HC BeforeCombat for QC
list1 = make_list("data/SLEinc_andHC/each_subset/LowExp_TMM_logCPM","_logCPM.txt",tag="count")
list1$subset = take_factor(list1$FILE_count,5:6,"_")
list2 = make_list("data/SLEinc_andHC/each_subset/each_subset_cond","_cond.txt",tag="cond")
list2$subset = take_factor(list2$FILE_cond,5:6,"_")
list3 = full_join(list1,list2,by="subset")
list3$status = "BeforeCombat"

# SLEinc+HC AfterCombat for QC
list4 = make_list("data/SLEinc_andHC/each_subset/LowExp_TMM_logCPM_Combat","_logCPM_Combat.txt",tag="count")
list4$subset = take_factor(list4$FILE_count,5:6,"_")
list5 = make_list("data/SLEinc_andHC/each_subset/each_subset_cond","_cond.txt",tag="cond")
list5$subset = take_factor(list5$FILE_cond,5:6,"_")
list6 = full_join(list4,list5,by="subset")
list6$status = "AfterCombat"

list7 = rbind(list3,list6)
list7$pop = "SLEinc_andHC"

# Here we substract unique patients from after Combat data
for(iii in 1:nrow(list6)){
  subset_tmp = list6$subset[iii]
  cond_tmp   = fread_FT(list6$PATH_cond[iii])
  count_tmp = fread_n(list6$PATH_count[iii])
  cond_tmp_u = cond_tmp %>% filter(id%in%SLEunique_and_HC$id)
  count_tmp_u = count_tmp[,is.element(colnames(count_tmp),cond_tmp_u$name)]
  write.table_n_2(count_tmp_u,"Gene",paste0("data/SLEunique_andHC/each_subset/LowExp_TMM_logCPM_Combat/",today(),"_COI_SLEunique_andHC_",subset_tmp,"_logCPM_Combat.txt"))  
}

# Similarly, we substract duplicated patients from after Combat data
for(iii in 1:nrow(list6)){
  subset_tmp = list6$subset[iii]
  cond_tmp   = fread_FT(list6$PATH_cond[iii])
  count_tmp = fread_n(list6$PATH_count[iii])
  cond_tmp_d = cond_tmp %>% filter(!id%in%SLEunique_and_HC$id)
  count_tmp_d = count_tmp[,is.element(colnames(count_tmp),cond_tmp_d$name)]
  write.table_n_2(count_tmp_d,"Gene",paste0("data/SLE_duplicated/each_subset/LowExp_TMM_logCPM_Combat/",today(),"_COI_SLE_duplicated_",subset_tmp,"_logCPM_Combat.txt"))  
}

# SLEunique+HC AfterCombat for subsequent analysis
list8 = make_list("data/SLEunique_andHC/each_subset/LowExp_TMM_logCPM_Combat","_logCPM_Combat.txt",tag="count")
list8$subset = take_factor(list8$FILE_count,5:6,"_")
list9 = make_list("data/SLEunique_andHC/each_subset/each_subset_cond","_cond.txt",tag="cond")
list9$subset = take_factor(list9$FILE_cond,5:6,"_")
list10 = full_join(list8,list9,by="subset")
list10$status = "AfterCombat"
list10$pop = "SLEunique_andHC"

list=rbind(list7,list10)

write.table_FT_2(list, paste0("tmp_job_list/",today(),"_COI_Eachpopcond_Eachsubsets_logCPM_forPCA_list.txt"))

############################################################################################
# PCA for one cell type (Example)
############################################################################################
LIST="tmp_job_list/211124_COI_Eachpopcond_Eachsubsets_logCPM_forPCA_list.txt"
# task_id=1

list          = fread_FT(LIST)

pop_tmp       = list[task_id,7]
status_tmp    = list[task_id,6]
subset_tmp    = list[task_id,3]
count_path    = list[task_id,1]
cond_path     = list[task_id,4]

paste0("pop_tmp : ",pop_tmp) %>% print()
paste0("status_tmp : ",status_tmp) %>% print()
paste0("subset_tmp : ",subset_tmp) %>% print()
paste0("count_path : ",count_path) %>% print()
paste0("cond_path : ",cond_path) %>% print()

logCPM_tmp   = fread_n(count_path)
condition    = fread_FT(cond_path)

name_list = data.frame(name=colnames(logCPM_tmp))
name_list = left_join(name_list,condition%>%select(name,id,run,disease,subset,lot,sequencer,lotseq),by="name") 

# PCA 
out_f  = paste0("PCA_res/",pop_tmp,"/",status_tmp,"/")
out_ff = paste0(today(),"_COI_",pop_tmp,"_",status_tmp,"_")
  
# filter 10000 genes based on variance (all cell-type: 10000 genes)
logCPM_tmp   = fread_n(count_path) 
logCPM_tmp_2   = logCPM_tmp %>% mutate(var=apply(.,1,var))
logCPM_tmp_f   = logCPM_tmp_2[order(logCPM_tmp_2$var,decreasing=T),] %>% .[1:10000,] %>% select(-var)
  
# PCA
result = prcomp(t(logCPM_tmp_f), scale = T)
df1_tmp        = as.data.frame(result$x) %>% rownames_to_column("name") %>%
                  right_join(name_list,.,by="name")
Importance_tmp = round(summary(result)$importance[2,]*100,digits=1)
rotation_tmp = summary(result)$rotation %>% as.data.frame()

write.table_FT_2(df1_tmp,paste0(out_f,"PCdata/",out_ff,subset_tmp,"_PCdata.txt"))
write.table_FT_2(data.frame(PC=names(Importance_tmp),importance=Importance_tmp),paste0(out_f,"PCimportance/",out_ff,subset_tmp,"_","PCimportance.txt"))
write.table_n_2(rotation_tmp,"Gene",paste0(out_f,"rotation/",out_ff,subset_tmp,"_rotation.txt"))
dir.create_p(paste0(out_f,"tmp_Rdata/"))
save(result, file= paste0(out_f,"tmp_Rdata/",out_ff,subset_tmp,"_PCAres_tmp1.Rdata"))

############################################################################################
# PCA plot each cell type
############################################################################################

label_df = as.data.frame(label) %>% rownames_to_column("subset")
labeller = as_labeller(label)
celltype_corresp=fread_FT("data_ref/COI_27subset_color_list.txt")

col2 = c("#ff2800","#faf500","#7f878f","#35a16b","#0041ff","#9a0079")
col3 = c("#ff2800","#faf500","NA","#35a16b","#0041ff","#9a0079")

# SLEincandHC
# Before and after Combat
# Color: Batch
list = make_list("PCA_res/SLEinc_andHC","_PCdata.txt")
list$subset = take_factor(list$FILE,7:8,"_")
list$status = take_factor(list$FILE,5,"_")
list$filter = take_factor(list$FILE,6,"_")
list$object = take_factor(list$FILE,3:4,"_")
list$object_filter_status = paste0(list$object,"_",list$filter,"_",list$status)



for(kkk in 1:length(unique(list$object_filter_status))){
   cond_tmp = unique(list$object_filter_status)[kkk]
   list_tmp   = list %>% filter(object_filter_status==cond_tmp)
   object_tmp = unique(list_tmp$object)
   filter_tmp = unique(list_tmp$filter)
   status_tmp = unique(list_tmp$status)

   for(iii in 1:nrow(list_tmp)){
     subset_tmp = list_tmp$subset[iii]
     df_tmp = fread_FT(list_tmp$PATH[iii]) %>% select(name,id,run,disease,subset,lot,sequencer,lotseq,PC1,PC2,PC3,PC4)
     if(iii==1){res_sum=df_tmp}else{res_sum=rbind(res_sum,df_tmp)}
   }
   res_sum$disease   = ifelse(res_sum$disease=="0HC","HC","SLE")
   res_sum$disease   = factor(res_sum$disease,levels=c("HC","SLE"))
   res_sum$lotseq    = factor(res_sum$lotseq,levels=sort(unique(res_sum$lotseq)))
   res_sum$subset    = factor(res_sum$subset,levels=celltype_reordered_27)

   ### color: lotseq
   #res_sum = res_sum[order(res_sum$disease),]
   p =  ggplot(res_sum,aes(x=PC1,y=PC2,fill=lotseq))+
          geom_point(size = 0.5,aes(shape=disease,color=lotseq))+theme_classic()+
          scale_shape_manual(values = c(3,21))+
          scale_fill_manual(values = col3)+
          scale_color_manual(values = col2)+
          facet_wrap(~ subset, labeller=labeller,ncol = 7) + 
          theme(axis.text.x=element_text(colour="black",size=14),
                axis.text.y=element_text(colour="black",size=14),
                axis.title.x=element_text(colour="black",size=14),
                axis.title.y=element_text(colour="black",size=14),
                plot.title=element_blank(),
                strip.background=element_blank(),
                strip.text=element_text(colour="black",size=14),
                panel.background =element_rect(fill="transparent",colour="black",size=1),
                panel.grid = element_blank(),
                legend.position="right",
                legend.title=element_blank(),
                legend.text=element_text(colour="black",size=14))+
          guides(shape=guide_legend(override.aes=list(size=3)))+
          guides(colour=guide_legend(override.aes=list(size=3)))+
          scale_x_continuous(breaks=c(-100,0,100))+
          scale_y_continuous(breaks=c(-100,0,100))+
          labs(x=paste0("PC1"),y=paste0("PC2"))

   pdf_3(paste0("PCA_plot/",cond_tmp,"/lotseq/",today(),"_COI_",cond_tmp,"_27subsets_Eachsubset_PC12plot_lotseq.pdf"),h=7,w=12)
    plot(p)
   dev.off()

   p =  ggplot(res_sum,aes(x=PC3,y=PC4,fill=lotseq))+
          geom_point(size = 0.5,aes(shape=disease,color=lotseq))+theme_classic()+
          scale_shape_manual(values = c(3,21))+
          scale_fill_manual(values = col3)+
          scale_color_manual(values = col2)+
          facet_wrap(~ subset, labeller=labeller,ncol = 7) + 
          theme(axis.text.x=element_text(colour="black",size=14),
                axis.text.y=element_text(colour="black",size=14),
                axis.title.x=element_text(colour="black",size=14),
                axis.title.y=element_text(colour="black",size=14),
                plot.title=element_blank(),
                strip.background=element_blank(),
                strip.text=element_text(colour="black",size=14),
                panel.background =element_rect(fill="transparent",colour="black",size=1),
                panel.grid = element_blank(),
                legend.position="right",
                legend.title=element_blank(),
                legend.text=element_text(colour="black",size=14))+
          guides(shape=guide_legend(override.aes=list(size=3)))+
          guides(colour=guide_legend(override.aes=list(size=3)))+
          scale_x_continuous(breaks=c(-100,0,100))+
          scale_y_continuous(breaks=c(-100,0,100))+
          labs(x=paste0("PC3"),y=paste0("PC4"))

   pdf_3(paste0("PCA_plot/",cond_tmp,"/lotseq/",today(),"_COI_",cond_tmp,"_27subsets_Eachsubset_PC34plot_lotseq.pdf"),h=7,w=12)
    plot(p)
   dev.off()
}

# SLEuniqueandHC
# after Combat
# Color: Disease
list = make_list("PCA_res/SLEunique_andHC","_PCdata.txt")
list$subset = take_factor(list$FILE,7:8,"_")
list$status = take_factor(list$FILE,5,"_")
list$filter = take_factor(list$FILE,6,"_")
list$object = take_factor(list$FILE,3:4,"_")
list$object_filter_status = paste0(list$object,"_",list$filter,"_",list$status)

for(kkk in 1:length(unique(list$object_filter_status))){
   cond_tmp = unique(list$object_filter_status)[kkk]
   list_tmp   = list %>% filter(object_filter_status==cond_tmp)
   object_tmp = unique(list_tmp$object)
   filter_tmp = unique(list_tmp$filter)
   status_tmp = unique(list_tmp$status)

   for(iii in 1:nrow(list_tmp)){
     subset_tmp = list_tmp$subset[iii]
     df_tmp = fread_FT(list_tmp$PATH[iii]) %>% select(name,id,run,disease,subset,lot,sequencer,lotseq,PC1,PC2,PC3,PC4)
     if(iii==1){res_sum=df_tmp}else{res_sum=rbind(res_sum,df_tmp)}
   }
   res_sum$disease   = ifelse(res_sum$disease=="0HC","HC","SLE")
   res_sum$disease   = factor(res_sum$disease,levels=c("HC","SLE"))
   res_sum$lotseq    = factor(res_sum$lotseq,levels=sort(unique(res_sum$lotseq)))
   res_sum$subset    = factor(res_sum$subset,levels=celltype_reordered_27)

   ### color: disease
   p =  ggplot(res_sum,aes(x=PC1,y=PC2))+
          geom_point(size = 0.5,aes(color=disease))+theme_classic()+
          facet_wrap(~ subset, labeller=labeller,ncol = 7) + 
          theme(axis.text.x=element_text(colour="black",size=14),
                axis.text.y=element_text(colour="black",size=14),
                axis.title.x=element_text(colour="black",size=14),
                axis.title.y=element_text(colour="black",size=14),
                plot.title=element_blank(),
                strip.background=element_blank(),
                strip.text=element_text(colour="black",size=14),
                panel.background =element_rect(fill="transparent",colour="black",size=1),
                panel.grid = element_blank(),
                legend.position="right",
                legend.title=element_blank(),
                legend.text=element_text(colour="black",size=14))+
          guides(colour=guide_legend(override.aes=list(size=3)))+
          scale_x_continuous(breaks=c(-100,0,100))+
          scale_y_continuous(breaks=c(-100,0,100))+
          labs(x="PC1",y="PC2")

   pdf_3(paste0("PCA_plot/",cond_tmp,"/disease/",today(),"_COI_",cond_tmp,"_27subsets_Eachsubset_PC12plot_disease.pdf"),h=7,w=12)
    plot(p)
   dev.off()

   p =  ggplot(res_sum,aes(x=PC3,y=PC4))+
          geom_point(size = 0.5,aes(color=disease))+theme_classic()+
          facet_wrap(~ subset, labeller=labeller,ncol = 7) + 
          theme(axis.text.x=element_text(colour="black",size=14),
                axis.text.y=element_text(colour="black",size=14),
                axis.title.x=element_text(colour="black",size=14),
                axis.title.y=element_text(colour="black",size=14),
                plot.title=element_blank(),
                strip.background=element_blank(),
                strip.text=element_text(colour="black",size=14),
                panel.background =element_rect(fill="transparent",colour="black",size=1),
                panel.grid = element_blank(),
                legend.position="right",
                legend.title=element_blank(),
                legend.text=element_text(colour="black",size=14))+
          guides(colour=guide_legend(override.aes=list(size=3)))+
          scale_x_continuous(breaks=c(-100,0,100))+
          scale_y_continuous(breaks=c(-100,0,100))+
          labs(x="PC3",y="PC4")
          
   pdf_3(paste0("PCA_plot/",cond_tmp,"/disease/",today(),"_COI_",cond_tmp,"_27subsets_Eachsubset_PC34plot_disease.pdf"),h=7,w=12)
    plot(p)
   dev.off()
}
############################################################################################
# representative: Th1, DNB, Neu
list_tmp = list

   for(iii in 1:nrow(list_tmp)){
     subset_tmp = list_tmp$subset[iii]
     df_tmp = fread_FT(list_tmp$PATH[iii]) %>% select(name,id,run,disease,subset,lot,sequencer,lotseq,PC1,PC2,PC3,PC4)
     if(iii==1){res_sum=df_tmp}else{res_sum=rbind(res_sum,df_tmp)}
   }
   res_sum$disease   = ifelse(res_sum$disease=="0HC","HC","SLE")
   res_sum$disease   = factor(res_sum$disease,levels=c("HC","SLE"))
   res_sum$lotseq    = factor(res_sum$lotseq,levels=sort(unique(res_sum$lotseq)))
   
   res_sum_lim = res_sum %>%filter(subset%in%c("E02_CD16nMo","B01_Th1","A04_DNB"))
   res_sum_lim$subset    = factor(res_sum_lim$subset,levels=c("E02_CD16nMo","B01_Th1","A04_DNB"))
   
   ### color: disease
   p =  ggplot(res_sum_lim,aes(x=PC1,y=PC2))+
          geom_point(size = 1,aes(color=disease))+theme_classic()+
          facet_wrap(~ subset, labeller=labeller,ncol = 3) + 
          theme(axis.text.x=element_text(colour="black",size=14),
                axis.text.y=element_text(colour="black",size=14),
                axis.title.x=element_text(colour="black",size=14),
                axis.title.y=element_text(colour="black",size=14),
                plot.title=element_blank(),
                strip.background=element_blank(),
                strip.text=element_text(colour="black",size=14),
                panel.background =element_rect(fill="transparent",colour="black",size=1),
                panel.grid = element_blank(),
                legend.position="right",
                legend.title=element_blank(),
                legend.text=element_text(colour="black",size=14))+
          guides(colour=guide_legend(override.aes=list(size=2)))+
          scale_x_continuous(limits=c(-150,150),breaks=c(-100,0,100))+
          scale_y_continuous(limits=c(-150,150),breaks=c(-100,0,100))+
          labs(x="PC1",y="PC2")

   pdf_3(paste0("PCA_plot/",cond_tmp,"/disease/",today(),"_COI_",cond_tmp,"_27subsets_Eachsubset_PC12plot_disease_rep.pdf"),h=2.5,w=7)
    plot(p)
   dev.off()

##############################################################################################################################
# Supple Table
# 27 celltype PC 1-7 loading
##############################################################################################################################
list1 = make_list("PCA_res/SLEunique_andHC/AfterCombat/rotation","_rotation.txt")
list1$subset = take_factor(list1$FILE,7:8,"_")

for (kkk in 1:nrow(list1)){
  subset_tmp = list1$subset[kkk]
  rotation= fread_FT(list1$PATH[kkk]) %>% as.data.frame() %>% 
            mutate(subset=subset_tmp) %>%
            select(Gene,subset,PC1,PC2,PC3,PC4,PC5,PC6,PC7)
  rotation$PC1 = formatC(rotation$PC1,digits=2)
  rotation$PC2 = formatC(rotation$PC2,digits=2)
  rotation$PC3 = formatC(rotation$PC3,digits=2)
  rotation$PC4 = formatC(rotation$PC4,digits=2)
  rotation$PC5 = formatC(rotation$PC5,digits=2)
  rotation$PC6 = formatC(rotation$PC6,digits=2)
  rotation$PC7 = formatC(rotation$PC7,digits=2)
  if(kkk==1){rotation_sum=rotation}else{rotation_sum=rbind(rotation_sum,rotation)}
}
  
rotation_sum$subset = factor(rotation_sum$subset,levels=celltype_reordered_27) 
rotation_sum = rotation_sum[order(rotation_sum$subset),]

label_df = as.data.frame(label) %>% rownames_to_column("subset")
rotation_sum2 = rotation_sum %>% left_join(.,label_df,by="subset")　%>%
                            select(Gene,label,PC1,PC2,PC3,PC4,PC5,PC6,PC7)

rotation_sum2$Gene   = factor(rotation_sum2$Gene,levels=sort(unique(rotation_sum2$Gene)))  
rotation_sum2$label  = factor(rotation_sum2$label,levels=unique(rotation_sum2$label))
rotation_sum2 = rotation_sum2[order(rotation_sum2$Gene),]
rotation_sum2 = rotation_sum2[order(rotation_sum2$label),]
colnames(rotation_sum2)[2] = "Cell type"

write.table_FT_2(rotation_sum2,paste0(today(),"_COI_27subset_PC1to7_loadings.txt"))

unique(rotation_sum2$Gene) %>% length()
##############################################################################################################################

############################################################################################
# PCA and UMAP plot: all 6386 samples
############################################################################################
# Each cell type Combat logCPM
# use intersection genes of each cell type
cond = fread_FT("data/SLEinc_andHC/211124_COI_SLEinc_andHC_27subsets_cond_wclinical.txt")
dim(cond)
# [1] 6386  212

list = make_list("data/SLEinc_andHC/each_subset/LowExp_TMM_logCPM_Combat","_logCPM_Combat.txt")
list$subset = take_factor(list$FILE,5:6,"_")

for(iii in 1:nrow(list)){
  subset_tmp = list$subset[iii]
  df_tmp   = fread_FT(list$PATH[iii])
  if(iii==1){res_sum=df_tmp}else{res_sum=full_join(res_sum,df_tmp,by="Gene")}
}  
res_sum = res_sum %>% column_to_rownames("Gene")
dim(res_sum)
# [1] 16000  6386
res_sum_woNA = res_sum[!apply(res_sum,1,anyNA),]
dim(res_sum_woNA)
# [1] 8397 6386

####### PCA
PCA = prcomp(t(res_sum_woNA),scale=T)
df1_tmp       = as.data.frame(PCA$x) %>% rownames_to_column("name") %>% select(name,PC1,PC2,PC3,PC4,PC5,PC6) %>%        
                left_join(.,cond%>%select(name,subset,disease),by="name") %>%
                left_join(.,celltype_corresp%>%select(subset,COLOR1),by="subset")

write.table_FT_2(df1_tmp,paste0("Allsample_PCAUMAP/",today(),"_COI_SLEinc_andHC_27subsets_PCA.txt")) 
save(PCA,file=paste0("Allsample_PCAUMAP/",today(),"_COI_SLEinc_andHC_27subsets_PCA.Rdata"))

df1_tmp$subset = factor(df1_tmp$subset,levels=celltype_reordered_27)
df1_tmp$disease = factor(df1_tmp$disease,levels=c("0HC","1SLE"))
levels(df1_tmp$disease) = list('HC'="0HC",'SLE'="1SLE")
df1_tmp2 = df1_tmp[order(df1_tmp$subset),]
cols=unique(df1_tmp2$COLOR1)

p =  ggplot(df1_tmp,aes(x=PC1,y=PC2,shape=disease,color=factor(subset,labels=label)))+
       geom_point(size = 0.6)+theme_classic()+
       scale_shape_manual(values = c(3,20))+
       scale_color_manual(values = cols)+
       theme(axis.text.x=element_text(colour="black",size=14),
             axis.text.y=element_text(colour="black",size=14),
             axis.title.x=element_text(colour="black",size=14),
             axis.title.y=element_text(colour="black",size=14),
             plot.title=element_blank(),
             legend.position="right",
             legend.title=element_blank(),
             legend.text=element_text(colour="black",size=14))+
       guides(colour=guide_legend(override.aes=list(size=3)))+
       guides(shape=guide_legend(override.aes=list(size=3)))+
       labs(x=paste0("PC1"),y=paste0("PC2"))
pdf_3(paste0("Allsample_PCAUMAP/",today(),"_COI_SLEinc_andHC_27subsets_PCA.pdf"),h=5,w=8.5)
 plot(p)
dev.off()

####### UMAP
set.seed(1)
umap = uwot::umap(t(res_sum_woNA),verbose=T) %>% as.data.frame()
rownames(umap) = colnames(res_sum_woNA)

df_1 = umap %>% rownames_to_column("name") %>%
                left_join(.,cond%>%select(name,subset,disease),by="name") %>%
                left_join(.,celltype_corresp%>%select(subset,COLOR1),by="subset")
write.table_FT_2(df_1,paste0("Allsample_PCAUMAP/",today(),"_COI_SLEinc_andHC_27subsets_UMAP.txt")) 

df_1$subset  = factor(df_1$subset,levels=celltype_reordered_27)
df_1$disease = factor(df_1$disease,levels=c("0HC","1SLE"))
levels(df_1$disease) = list('HC'="0HC",'SLE'="1SLE")

p =  ggplot(df_1,aes(x=V1,y=V2,shape=disease,color=factor(subset,labels=label)))+
       geom_point(size = 0.6)+theme_classic()+
       scale_shape_manual(values = c(3,20))+
       scale_color_manual(values = cols)+
       theme(axis.text.x=element_text(colour="black",size=14),
             axis.text.y=element_text(colour="black",size=14),
             axis.title.x=element_text(colour="black",size=14),
             axis.title.y=element_text(colour="black",size=14),
             plot.title=element_blank(),
             legend.position="right",
             legend.title=element_blank(),
             legend.text=element_text(colour="black",size=14))+
       guides(colour=guide_legend(override.aes=list(size=3)))+
       guides(shape=guide_legend(override.aes=list(size=3)))+
       labs(x=paste0("UMAP1"),y=paste0("UMAP2"))
pdf_3(paste0("Allsample_PCAUMAP/",today(),"_COI_SLEinc_andHC_27subsets_UMAP.pdf"),h=5,w=8.5)
 plot(p)
dev.off()

############################################################################################
############################################################################################

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
 [1] uwot_0.1.10         Matrix_1.3-4        ggsci_2.9          
 [4] RColorBrewer_1.1-2  sva_3.38.0          BiocParallel_1.24.1
 [7] genefilter_1.72.1   mgcv_1.8-36         nlme_3.1-152       
[10] edgeR_3.32.1        limma_3.46.0        data.table_1.14.0  
[13] forcats_0.5.1       stringr_1.4.0       dplyr_1.0.7        
[16] purrr_0.3.4         readr_1.4.0         tidyr_1.1.3        
[19] tibble_3.1.2        ggplot2_3.3.5       tidyverse_1.3.1    

loaded via a namespace (and not attached):
 [1] Biobase_2.50.0       httr_1.4.2           bit64_4.0.5         
 [4] jsonlite_1.7.2       splines_4.0.2        modelr_0.1.8        
 [7] assertthat_0.2.1     stats4_4.0.2         blob_1.2.1          
[10] cellranger_1.1.0     pillar_1.6.1         RSQLite_2.2.7       
[13] backports_1.2.1      lattice_0.20-41      glue_1.4.2          
[16] rvest_1.0.0          colorspace_2.0-2     XML_3.99-0.6        
[19] pkgconfig_2.0.3      broom_0.7.8          haven_2.4.1         
[22] xtable_1.8-4         scales_1.1.1         annotate_1.68.0     
[25] generics_0.1.0       IRanges_2.24.1       ellipsis_0.3.2      
[28] cachem_1.0.5         withr_2.4.2          BiocGenerics_0.36.1 
[31] cli_2.5.0            survival_3.2-11      magrittr_2.0.1      
[34] crayon_1.4.1         readxl_1.3.1         memoise_2.0.0       
[37] ps_1.6.0             fs_1.5.0             fansi_0.5.0         
[40] xml2_1.3.2           tools_4.0.2          hms_1.1.0           
[43] matrixStats_0.59.0   lifecycle_1.0.0      S4Vectors_0.28.1    
[46] munsell_0.5.0        reprex_2.0.0         locfit_1.5-9.4      
[49] AnnotationDbi_1.52.0 compiler_4.0.2       rlang_0.4.11        
[52] grid_4.0.2           rstudioapi_0.13      gtable_0.3.0        
[55] DBI_1.1.1            R6_2.5.0             lubridate_1.7.10    
[58] fastmap_1.1.0        bit_4.0.4            utf8_1.2.1          
[61] stringi_1.6.2        parallel_4.0.2       Rcpp_1.0.6          
[64] vctrs_0.3.8          dbplyr_2.1.1         tidyselect_1.1.1    


























