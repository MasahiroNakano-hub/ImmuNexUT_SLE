#---------------------------------------#
# 211129 Nakano
# edgeR 
#---------------------------------------#

source("my.source.R")
options(stringsAsFactors=F)

# unique pts data list
list1 = make_list("~/ws/2021/211124_COI_SLE_HC_run27_finaldata_Combat_PCA/data/SLEunique_andHC/each_subset/each_subset_count","_count.txt",tag="count")
list1$subset = take_factor(list1$FILE_count,5:6,"_")
list2 = make_list("~/ws/2021/211124_COI_SLE_HC_run27_finaldata_Combat_PCA/data/SLEunique_andHC/each_subset/each_subset_cond","_cond.txt",tag="batch")
list2$subset = take_factor(list2$FILE_batch,5:6,"_")
list = full_join(list1,list2,by="subset")
write.table_FT_2(list,paste0("tmp_job_list/",today(),"_COI_SLEandHC_edgeR_list.txt"))

q()
############################################################################################################################################################
# edgeR DEGs
# Example: One cell type disease-state and activity
############################################################################################################################################################
LIST="tmp_job_list/211129_COI_SLEandHC_edgeR_list.txt"
##  task_id=1

list       = fread_FT(LIST)
subset_tmp = list[task_id,3]
count_path = list[task_id,1]
batch_path = list[task_id,4]

paste0("subset_tmp : ",subset_tmp) %>% print()
paste0("count_path : ",count_path) %>% print()
paste0("batch_path : ",batch_path) %>% print()

###########################################################################

out_f  = paste0("~/ws/2021/211129_COI_SLE_HC_run27_edgeR/res/stateactivity/")
out_ff = paste0(today(),"_COI_",subset_tmp,"_")

###################
################################

count_tmp  = fread_n(count_path)
cond_tmp   = fread_FT(batch_path)

chrY_gene = fread_FF("data_ref/hg38_chrY_gene_list.txt")[,1]

count_tmp = count_tmp[!is.element(rownames(count_tmp),chrY_gene),]

## Check
if(ncol(count_tmp)==nrow(cond_tmp)){"sample and cond were same num"}else{"please check condition list"}

cond_tmp_2  = cond_tmp[is.element(cond_tmp$name,colnames(count_tmp)),]
count_tmp_2 = count_tmp[,cond_tmp_2$name]

## Check
print("## final use data number ##")
table_freq(cond_tmp$Activity)
print("## final use data number ##")

cond_tmp_2$Activity  = factor(cond_tmp_2$Activity, levels=sort(unique(cond_tmp_2$Activity)))
cond_tmp_2$lotseq    = paste0(cond_tmp_2$lot,"_",cond_tmp_2$sequencer) 

## low expression data filtering
cpm_tmp    = cpm(count_tmp_2)
filter1    = rownames(count_tmp_2)[exp_filter(count_tmp_2,exp=10,p=0.15)]
filter2    = rownames(cpm_tmp)[exp_filter(cpm_tmp,exp=1,p=0.15)]
count_tmp_3= count_tmp_2[intersect(filter1,filter2),]

## Check
print("## Low Exp Gene filtering ##")
paste0("original gene num : ",nrow(count_tmp_2))     %>% print()
paste0("after filter gene num : ",nrow(count_tmp_3)) %>% print()
print("## Low Exp Gene filtering ##")

### edgeR TMM norm
dge     = DGEList(counts=count_tmp_3, genes=rownames(count_tmp_3))
dge_n   = calcNormFactors( dge ,method="TMM")

#####################################
#HC vs Inactive
#####################################
design_1     = model.matrix(~ 0+Activity+lotseq+age+gender , data= cond_tmp_2 )

## estimate dispersion
y = estimateGLMCommonDisp ( dge_n , design_1 )
y = estimateGLMTrendedDisp( y     , design_1 )
y = estimateGLMTagwiseDisp( y     , design_1 )

fit = glmQLFit( y, design_1 )
con = makeContrasts(Activity1Inactive - Activity0HC, levels=design_1)
lrt = glmQLFTest( fit , contrast = con) 
result_edgeR = topTags(lrt,n=nrow(count_tmp_3))$table

## plotSmear
isDE    = decideTestsDGE( lrt ) %>% as.logical()
DEnames = row.names(y)[isDE]

df_tmp  = data.frame(Gene=result_edgeR$genes,logFC=result_edgeR$logFC,logCPM=result_edgeR$logCPM)
df_tmp$isDEG = "nonDEG"
df_tmp$isDEG[result_edgeR$FDR<0.05] = "DEG_FDR0.05"
df_tmp$isDEG[result_edgeR$FDR<0.01] = "DEG_FDR0.01"
df_tmp$isDEG = factor(df_tmp$isDEG,levels=c("nonDEG","DEG_FDR0.05","DEG_FDR0.01"))
df_tmp$up_down = NA
df_tmp$up_down[df_tmp$logFC>=0] = "Up"
df_tmp$up_down[df_tmp$logFC<0]  = "Down"
df_tmp$isDEG_2 = paste0(df_tmp$up_down,"_",df_tmp$isDEG) %>% gsub("Up_nonDEG","nonDEG",.) %>% gsub("Down_nonDEG","nonDEG",.)

DEG_num_tmp = table_freq(df_tmp$isDEG_2)
colnames(DEG_num_tmp)[2] = paste0(subset_tmp,"_InactivevsHC")
write.table_FT_2(DEG_num_tmp,paste0(out_f,"InactivevsHC/DEG_table/",out_ff,"InactivevsHC_DEG_table_list.txt"))


## DEG list
result_DEG_0.05   = result_edgeR[result_edgeR$FDR<0.05,]
result_DEG_f_0.05 = result_DEG_0.05 %>% arrange(desc(logFC))  

result_DEG_0.01   = result_edgeR[result_edgeR$FDR<0.01,]
result_DEG_f_0.01 = result_DEG_0.01 %>% arrange(desc(logFC))  

write.table_FT_2(result_DEG_f_0.05,paste0(out_f,"InactivevsHC/DEGlist005/",out_ff,"InactivevsHC_DEGlist_FDR005.txt"))
write.table_FT_2(result_DEG_f_0.01,paste0(out_f,"InactivevsHC/DEGlist001/",out_ff,"InactivevsHC_DEGlist_FDR001.txt"))
write.table_FT_2(result_edgeR,paste0(out_f,"InactivevsHC/GLMresult/",out_ff,"InactivevsHC_GLMresult.txt"))


#####################################
#HDA vs Inactive
#####################################
design_2     = model.matrix(~ 0+Activity+lotseq+age+gender+PSLmg+HCQ+MMF+TAC , data= cond_tmp_2 )

## estimate dispersion
y2 = estimateGLMCommonDisp ( dge_n , design_2 )
y2 = estimateGLMTrendedDisp( y2     , design_2 )
y2 = estimateGLMTagwiseDisp( y2     , design_2 )

fit = glmQLFit( y2, design_2 )
con = makeContrasts(Activity4HDA - Activity1Inactive, levels=design_2)
lrt = glmQLFTest( fit , contrast = con) 
result_edgeR = topTags(lrt,n=nrow(count_tmp_3))$table

## plotSmear
isDE    = decideTestsDGE( lrt ) %>% as.logical()
DEnames = row.names(y2)[isDE]

df_tmp  = data.frame(Gene=result_edgeR$genes,logFC=result_edgeR$logFC,logCPM=result_edgeR$logCPM)
df_tmp$isDEG = "nonDEG"
df_tmp$isDEG[result_edgeR$FDR<0.05] = "DEG_FDR0.05"
df_tmp$isDEG[result_edgeR$FDR<0.01] = "DEG_FDR0.01"
df_tmp$isDEG = factor(df_tmp$isDEG,levels=c("nonDEG","DEG_FDR0.05","DEG_FDR0.01"))
df_tmp$up_down = NA
df_tmp$up_down[df_tmp$logFC>=0] = "Up"
df_tmp$up_down[df_tmp$logFC<0]  = "Down"
df_tmp$isDEG_2 = paste0(df_tmp$up_down,"_",df_tmp$isDEG) %>% gsub("Up_nonDEG","nonDEG",.) %>% gsub("Down_nonDEG","nonDEG",.)

DEG_num_tmp = table_freq(df_tmp$isDEG_2)
colnames(DEG_num_tmp)[2] = paste0(subset_tmp,"_HDAvsInactive")
write.table_FT_2(DEG_num_tmp,paste0(out_f,"HDAvsInactive/DEG_table/",out_ff,"HDAvsInactive_DEG_table_list.txt"))

## DEG list
result_DEG_0.05   = result_edgeR[result_edgeR$FDR<0.05,]
result_DEG_f_0.05 = result_DEG_0.05 %>% arrange(desc(logFC))  

result_DEG_0.01   = result_edgeR[result_edgeR$FDR<0.01,]
result_DEG_f_0.01 = result_DEG_0.01 %>% arrange(desc(logFC))  

write.table_FT_2(result_DEG_f_0.05,paste0(out_f,"HDAvsInactive/DEGlist005/",out_ff,"HDAvsInactive_DEGlist_FDR005.txt"))
write.table_FT_2(result_DEG_f_0.01,paste0(out_f,"HDAvsInactive/DEGlist001/",out_ff,"HDAvsInactive_DEGlist_FDR001.txt"))
write.table_FT_2(result_edgeR,paste0(out_f,"HDAvsInactive/GLMresult/",out_ff,"HDAvsInactive_GLMresult.txt"))

######################################################################################################################################################################################################

##############################################################################
# DEGnum barplot
##############################################################################
suppressPackageStartupMessages(library(ggsci))
suppressPackageStartupMessages(library(RColorBrewer))

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
###################################################################################################

################################
# stateactivity DEGs
################################
list1         = make_list("res/stateactivity","_DEGlist_FDR005.txt")
list1$subset  = take_factor(list1$FILE,3:4,"_")
list1$type    = take_factor(list1$FILE,5,"_")

for(kkk in 1:length(unique(list1$type))){
   type_tmp  = unique(list1$type)[kkk]
   list_tmp  = list1 %>% filter(type==type_tmp)

 for(iii in 1:nrow(list_tmp)){
  subset_tmp   = list_tmp$subset[iii]
  df_1         = read.table_FT(list_tmp$PATH[iii])
  df_tmp = data.frame(DEGtype=type_tmp,subset=subset_tmp,DEGnum=nrow(df_1))
  if(iii==1){res_sum=df_tmp}else{res_sum=rbind(res_sum,df_tmp)}
  }
if(kkk==1){res_sum2=res_sum}else{res_sum2=rbind(res_sum2,res_sum)}
}

res_sum3 = left_join(res_sum2,celltype_corresp%>%select(subset,lineage),by="subset") 
write.table_FT_2(res_sum3,paste0("res/stateactivity/DEGnum/",today(),"_COI_stateactivity4type_DEGnum.txt"))

res_sum3$lineage = factor(res_sum3$lineage,levels=c("CD4","CD8","NK","B","Monocyte","DC","Neutrophil"))
res_sum3$subset  = factor(res_sum3$subset,levels=celltype_reordered_27)
res_sum3$DEGnum  = as.numeric(res_sum3$DEGnum)

# State activity DEGnum group
res_sum_pattern1 = res_sum3 %>% filter(DEGtype%in%c("InactivevsHC","HDAvsInactive"))
res_sum_pattern1$DEGtype = factor(res_sum_pattern1$DEGtype,levels=c("InactivevsHC","HDAvsInactive"))
col = pal_nejm("default")(2)[c(2,1)]

res_sum_pattern11=res_sum_pattern1 %>% filter(DEGtype=="InactivevsHC")
mean(res_sum_pattern11$DEGnum)
# [1] 2097.852
res_sum_pattern12=res_sum_pattern1 %>% filter(DEGtype=="HDAvsInactive")
mean(res_sum_pattern12$DEGnum)
# [1] 2114.444

################################
# treatment DEGs
################################
list2         = make_list("res/treatment","_DEGlist_FDR005.txt")
list2$subset  = take_factor(list2$FILE,3:4,"_")
list2$type    = take_factor(list2$FILE,5,"_")

for(kkk in 1:length(unique(list2$type))){
   type_tmp  = unique(list2$type)[kkk]
   list_tmp  = list2 %>% filter(type==type_tmp)

 for(iii in 1:nrow(list_tmp)){
  subset_tmp   = list_tmp$subset[iii]
  df_1         = read.table_FT(list_tmp$PATH[iii])
  df_tmp = data.frame(DEGtype=type_tmp,subset=subset_tmp,DEGnum=nrow(df_1))
  if(iii==1){res_sum=df_tmp}else{res_sum=rbind(res_sum,df_tmp)}
  }
if(kkk==1){res_sum2=res_sum}else{res_sum2=rbind(res_sum2,res_sum)}
}

res_sum3 = left_join(res_sum2,celltype_corresp%>%select(subset,lineage),by="subset") 
write.table_FT_2(res_sum3,paste0("res/treatment/",today(),"_COI_Tx_DEGnum.txt"))
# res_sum3=fread_FT("res/treatment/211130_COI_Tx_DEGnum.txt")
res_sum3$lineage = factor(res_sum3$lineage,levels=c("CD4","CD8","NK","B","Monocyte","DC","Neutrophil"))
res_sum3$subset  = factor(res_sum3$subset,levels=celltype_reordered_27)
res_sum3$DEGnum  = as.numeric(res_sum3$DEGnum)

# MMF only
res_sum_pattern3 = res_sum3 %>% filter(DEGtype=="MMF")

p = ggplot()+
     geom_bar(data=res_sum_pattern3, aes(x=subset,y=DEGnum,fill=lineage),stat="identity")+
     theme_classic()+
     #facet_wrap(~ DEGtype, labeller=facet5,ncol = 3) + 
     scale_fill_manual(values = col2)+
     theme(axis.text.x=element_text(colour="black",angle=45,hjust=1,vjust=1,size=16),
           axis.text.y=element_text(colour="black",size=16),
           axis.title.x=element_blank(),
           axis.title.y=element_text(colour="black",size=16),
           plot.title=element_blank(),
           strip.text=element_text(colour="black",size=16),
           legend.position="right",
           legend.title=element_text(colour="black",size=16),
           legend.text=element_text(colour="black",size=16))+
     scale_x_discrete(labels= label)+
     labs(y="Number of DEGs")

 pdf_3(paste0("res/treatment/",today(),"_COI_MMF_DEGnum_barplot_fig.pdf"),h=3.5,w=9.5)
  plot(p)
 dev.off()

###################################################################################################
# Disease-state, activity logFC
###################################################################################################

list1         = make_list("res/stateactivity/InactivevsHC/GLMresult","_GLMresult.txt",tag="state")
list1$subset  = take_factor(list1$FILE,3:4,"_")
list2         = make_list("res/stateactivity/HDAvsInactive/GLMresult","_GLMresult.txt",tag="activity")
list2$subset  = take_factor(list2$FILE,3:4,"_")

list = left_join(list1,list2,by="subset")

for(iii in 1:nrow(list)){
    subset_tmp         = list$subset[iii]
    df_tmp_s           = fread_FT(list$PATH_state[iii]) %>% select(genes,logFC,logCPM,PValue,FDR)
    colnames(df_tmp_s)[2:5] = paste0("state_",colnames(df_tmp_s)[2:5])

    df_tmp_a           = fread_FT(list$PATH_activity[iii]) %>% select(genes,logFC,logCPM,PValue,FDR)
    colnames(df_tmp_a)[2:5] = paste0("activity_",colnames(df_tmp_a)[2:5])

    df_tmp = full_join(df_tmp_s,df_tmp_a,by="genes")
    df_tmp$subset =  subset_tmp
    if(iii==1){res_sum=df_tmp}else{res_sum=rbind(res_sum,df_tmp)}
  }

# Here we color even for discordant genes (tentative)
res_sum$siggroup005 = ifelse(res_sum$state_FDR<0.05&res_sum$activity_FDR<0.05&(res_sum$state_logFC*res_sum$activity_logFC)>0,"significant_in_both",
                       ifelse(res_sum$state_FDR<0.05,"only_state",
                        ifelse(res_sum$activity_FDR<0.05,"only_activity","neither")))

res_sum$subset   = factor(res_sum$subset,levels=celltype_reordered_27)
res_sum$siggroup005   = factor(res_sum$siggroup005,levels=c("neither","only_state","only_activity","significant_in_both"))

# for visualization logFC lim =5,-5
res_sum_max5 = res_sum
res_sum_max5$state_logFC[res_sum_max5$state_logFC>5] = 5 
res_sum_max5$state_logFC[res_sum_max5$state_logFC< -5] = -5 
res_sum_max5$activity_logFC[res_sum_max5$activity_logFC>5] = 5 
res_sum_max5$activity_logFC[res_sum_max5$activity_logFC< -5] = -5 

res_sum_max5 = res_sum_max5[order(res_sum$siggroup005),]

label_group = c("Not significant","Disease-state only","Disease-activity only","Significant in both")

# 9x 3 ver
p =  ggplot(res_sum_max5,aes(x=state_logFC,y=activity_logFC))+
          geom_point(size=0.2,aes(color=factor(siggroup005,labels=label_group)))+
          theme_classic()+
          scale_color_manual(values = c("#cccccc","#0072B5FF","#BC3C29FF","#E18727FF"))+
          facet_wrap(~ subset, labeller=labeller,ncol = 9) + 
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
          scale_x_continuous(breaks=c(-4,-2,0,2,4),limits=c(-5,5))+
          scale_y_continuous(breaks=c(-4,-2,0,2,4),limits=c(-5,5))+
          labs(x=paste0("Inactive vs HC (logFC)"),y=paste0("HDA vs Inactive (logFC)"))

png(paste0("res/stateactivity/similarity/",today(),"_COI_27subsets_state_and_activity_FDR005_scatter_v2.png"),height=5,width=14,units="in",res=300)
 plot(p)
dev.off()

# representative 3 celltypes
# NO character
res_sum_max5_rep = res_sum_max5　%>% filter(subset%in%c("B05_NCD4","A05_PB","F01_Neu"))
res_sum_max5_rep$subset = factor(res_sum_max5_rep$subset,levels=c("B05_NCD4","A05_PB","F01_Neu"))

p =  ggplot(res_sum_max5_rep,aes(x=state_logFC,y=activity_logFC))+
          geom_point(size=0.2,aes(color=factor(siggroup005,labels=label_group)))+
          theme_classic()+
          scale_color_manual(values = c("#cccccc","#0072B5FF","#BC3C29FF","#E18727FF"))+
          facet_wrap(~ subset, labeller=labeller,ncol = 7) + 
          theme(axis.text.x=element_blank(),
                axis.text.y=element_blank(),
                axis.title.x=element_blank(),
                axis.title.y=element_blank(),
                plot.title=element_blank(),
                strip.background=element_blank(),
           		strip.text=element_blank(),
           		panel.background =element_rect(fill="transparent",colour="black",size=1),
           		panel.grid = element_blank(),
                legend.position="right",
                legend.title=element_blank(),
                legend.text=element_blank(),)+
          guides(colour=guide_legend(override.aes=list(size=1.5)))+
          scale_x_continuous(breaks=c(-4,-2,0,2,4),limits=c(-5,5))+
          scale_y_continuous(breaks=c(-4,-2,0,2,4),limits=c(-5,5))+
          labs(x=paste0("Inactive vs HC (logFC)"),y=paste0("HDA vs Inactive (logFC)"))

png(paste0("res/stateactivity/similarity/",today(),"_COI_27subsets_state_and_activity_FDR005_scatter_rep.png"),height=1.9,width=6,units="in",res=300)
 plot(p)
dev.off()

# Figure
res_sum_max5_rep$shape=ifelse(res_sum_max5_rep$subset=="B05_NCD4"&res_sum_max5_rep$genes=="PRDM1","Diamond",
                        ifelse(res_sum_max5_rep$subset=="A05_PB"&res_sum_max5_rep$genes=="PTP4A3","Diamond",
                         ifelse(res_sum_max5_rep$subset=="F01_Neu"&res_sum_max5_rep$genes=="ISG15","Diamond","Point")))
res_sum_max5_rep$shape=factor(res_sum_max5_rep$shape,levels=c("Diamond","Point"))

res_sum_max5_rep$siggroup005_2=ifelse(res_sum_max5_rep$shape=="Diamond","Diamond",as.character(res_sum_max5_rep$siggroup005))
res_sum_max5_rep$siggroup005_2   = factor(res_sum_max5_rep$siggroup005_2,levels=c("neither","only_state","only_activity","significant_in_both","Diamond"))
label_group_2 = c("Not significant","Disease-state only","Disease-activity only","Significant in both","")

res_sum_max5_rep=res_sum_max5_rep[order(res_sum_max5_rep$shape,decreasing=T),]

p =  ggplot(res_sum_max5_rep,aes(x=state_logFC,y=activity_logFC))+
          geom_point(aes(size=shape,shape=shape,stroke=1,fill=factor(siggroup005,labels=label_group),color=factor(siggroup005_2,labels=label_group_2)))+
          theme_classic()+
          scale_fill_manual(values = c("#cccccc","#0072B5FF","#BC3C29FF","#E18727FF"))+
          scale_color_manual(values = c("#cccccc","#0072B5FF","#BC3C29FF","#E18727FF","#2b2b2b"))+
          scale_shape_manual(values = c(23,16))+
          scale_size_manual(values = c(3,0.25))+
          facet_wrap(~ subset, labeller=labeller,ncol = 7) + 
          theme(axis.text.x=element_blank(),
                axis.text.y=element_blank(),
                axis.title.x=element_blank(),
                axis.title.y=element_blank(),
                plot.title=element_blank(),
                strip.background=element_blank(),
                strip.text=element_blank(),
                panel.background =element_rect(fill="transparent",colour="black",size=1),
                panel.grid = element_blank(),
                legend.position="right",
                legend.title=element_blank(),
                legend.text=element_blank(),)+
          guides(colour=guide_legend(override.aes=list(size=1.5)))+
          scale_x_continuous(breaks=c(-4,-2,0,2,4),limits=c(-5,5))+
          scale_y_continuous(breaks=c(-4,-2,0,2,4),limits=c(-5,5))+
          labs(x=paste0("Inactive vs HC (logFC)"),y=paste0("HDA vs Inactive (logFC)"))

png(paste0("Figure/",today(),"_COI_27subsets_state_and_activity_FDR005_scatter_rep.png"),height=1.9,width=6,units="in",res=300)
 plot(p)
dev.off()

##############################################################################
# Jaccard similarity
##############################################################################

res_sum$siggroup005 = ifelse(res_sum$state_FDR<0.05&res_sum$activity_FDR<0.05&(res_sum$state_logFC*res_sum$activity_logFC)>0,"significant_in_both",
					   ifelse(res_sum$state_FDR<0.05&res_sum$activity_FDR<0.05&(res_sum$state_logFC*res_sum$activity_logFC)<0,"discordant",
                        ifelse(res_sum$state_FDR<0.05,"only_state",
                         ifelse(res_sum$activity_FDR<0.05,"only_activity","neither"))))

write.table_FT_2(res_sum,paste0("res/stateactivity/",today(),"_COI_27subsets_state_activity_logFC_logCPM_P_FDR.txt"))

for(iii in 1:length(celltype_reordered_27)){
    subset_tmp = celltype_reordered_27[iii]
    res_tmp = res_sum %>% filter(subset==subset_tmp)
    state_tmp = res_tmp %>% filter(siggroup005=="only_state")
    activity_tmp = res_tmp %>% filter(siggroup005=="only_activity")
    both_tmp  = res_tmp %>% filter(siggroup005=="significant_in_both")
    dis_tmp   = res_tmp %>% filter(siggroup005=="discordant")
    jaccard_tmp = data.frame(subset=subset_tmp,state=nrow(state_tmp),activity=nrow(activity_tmp),both=nrow(both_tmp),
    						 discordant=nrow(dis_tmp),union=nrow(state_tmp)+nrow(activity_tmp)+nrow(both_tmp),
    						 all=nrow(state_tmp)+nrow(activity_tmp)+nrow(both_tmp)+nrow(dis_tmp),
    						 state_P=nrow(state_tmp)/(nrow(state_tmp)+nrow(activity_tmp)+nrow(both_tmp)),
    						 activity_P=nrow(activity_tmp)/(nrow(state_tmp)+nrow(activity_tmp)+nrow(both_tmp)),
    						 jaccard=nrow(both_tmp)/(nrow(state_tmp)+nrow(activity_tmp)+nrow(both_tmp)))
    if(iii==1){jaccard_sum=jaccard_tmp}else{jaccard_sum=rbind(jaccard_sum,jaccard_tmp)}
}

jaccard_sum=jaccard_sum%>%mutate(group=ifelse(jaccard>0.15,"Shared",
										ifelse(state>activity,"Disease-state dominant","Disease-activity dominant")))

jaccard_sum_state=jaccard_sum%>%filter(group=="Disease-state dominant")
jaccard_sum_state=jaccard_sum_state[order(jaccard_sum_state$state_P),]
jaccard_sum_activity=jaccard_sum%>%filter(group=="Disease-activity dominant")
jaccard_sum_activity=jaccard_sum_activity[order(jaccard_sum_activity$activity_P),]
jaccard_sum_shared=jaccard_sum%>%filter(group=="Shared")
jaccard_sum_shared=jaccard_sum_shared[order(jaccard_sum_shared$jaccard),]

jaccard_sum_ordered=bind_rows(list(jaccard_sum_state,jaccard_sum_activity,jaccard_sum_shared)) %>%
                    left_join(.,celltype_corresp[,c(1,3)],by="subset")

write.table_FT_2(jaccard_sum_ordered,paste0("res/stateactivity/similarity/",today(),"_COI_27subset_state_and_activity_Jaccard.txt"))

jaccard_sum_ordered$subset = factor(jaccard_sum_ordered$subset,levels=unique(jaccard_sum_ordered$subset))

p1 = ggplot()+
     geom_bar(data=jaccard_sum_ordered, aes(x=subset,y=union,fill=lineage),stat="identity")+
     theme_classic()+
     scale_fill_manual(values=col2)+
     theme(axis.text.x=element_blank(),
           axis.text.y=element_blank(),
           axis.title.x=element_blank(),
           axis.title.y=element_blank(),
           plot.title=element_blank(),
      　　　legend.position="none")+
     #scale_x_discrete(labels= label)+
　　　labs(y="Number of DEGs")

pdf_3(paste0("res/stateactivity/similarity/",today(),"_COI_27subset_union_DEGnum.pdf"),h=1.8,w=7.5)
 plot(p1)
dev.off()

p2 = ggplot()+
     geom_bar(data=jaccard_sum_ordered, aes(x=subset,y=union,fill=lineage),stat="identity")+
     theme_classic()+
     scale_fill_manual(values=col2)+
     theme(axis.text.x=element_blank(),
           axis.text.y=element_text(colour="black",size=16),
           axis.title.x=element_blank(),
           axis.title.y=element_text(colour="black",size=16),
           plot.title=element_blank(),
      　　　legend.position="right",
           legend.title=element_text(colour="black",size=16),
           legend.text=element_text(colour="black",size=16))+
     #scale_x_discrete(labels= label)+
　　　labs(y="Number of DEGs")

pdf_3(paste0("res/stateactivity/similarity/",today(),"_COI_27subset_union_DEGnum_wyaxislegend.pdf"),h=1.8,w=10)
 plot(p2)
dev.off()

jaccard_sum2=jaccard_sum_ordered%>%select(subset,jaccard,activity_P,state_P)
colnames(jaccard_sum2)[2:4]=c("Both significant","Disease-activity only","Disease-state only")
jaccard_sum2=jaccard_sum2 %>% pivot_longer(-subset,names_to="DEGtype",values_to="value")
jaccard_sum2$DEGtype=factor(jaccard_sum2$DEGtype,levels=unique(jaccard_sum2$DEGtype))

p3 = ggplot()+
    geom_bar(data=jaccard_sum2,aes(x=subset,y=value,fill=DEGtype),stat="identity",position="stack")+
    theme_classic()+
    scale_fill_manual(values=c("#E18727FF","#BC3C29FF","#0072B5FF"))+
    theme(axis.text.x=element_text(colour="black",angle=45,hjust=1,vjust=1,size=16),
           axis.text.y=element_blank(),
           axis.title.x=element_blank(),
           axis.title.y=element_blank(),
           plot.title=element_blank(),
           legend.position="none")+
     scale_x_discrete(labels= label)+
    labs(y="Proportion of DEGs")

pdf_3(paste0("res/stateactivity/similarity/",today(),"_COI_27subset_state_and_activity_Jaccard.pdf"),h=3.5,w=7.5)
 plot(p3)
dev.off()

p4 = ggplot()+
    geom_bar(data=jaccard_sum2,aes(x=subset,y=value,fill=DEGtype),stat="identity",position="stack")+
    theme_classic()+
    scale_fill_manual(values=c("#E18727FF","#BC3C29FF","#0072B5FF"))+
    theme(axis.text.x=element_text(colour="black",angle=45,hjust=1,vjust=1,size=16),
           axis.text.y=element_text(colour="black",size=16),
           axis.title.x=element_blank(),
           axis.title.y=element_text(colour="black",size=16),
           plot.title=element_blank(),
           legend.position="right",
           legend.title=element_text(colour="black",size=16),
           legend.text=element_text(colour="black",size=16))+
     scale_x_discrete(labels= label)+
    labs(y="Proportion of DEGs")

pdf_3(paste0("res/stateactivity/similarity/",today(),"_COI_27subset_state_and_activity_Jaccard_wyaxislegend.pdf"),h=3.5,w=10)
 plot(p4)
dev.off()

##############################################################################
# Extended: state and activity DEGs reordered by celltype group
##############################################################################
res_sum3=fread_FT("res/stateactivity/DEGnum/211130_COI_stateactivity4type_DEGnum.txt")

res_sum3$lineage = factor(res_sum3$lineage,levels=c("CD4","CD8","NK","B","Monocyte","DC","Neutrophil"))
res_sum3$subset  = factor(res_sum3$subset,levels=levels(jaccard_sum_ordered$subset))
res_sum3$DEGnum  = as.numeric(res_sum3$DEGnum)

# State activity DEGnum groupe
res_sum_pattern1 = res_sum3 %>% filter(DEGtype%in%c("InactivevsHC","HDAvsInactive"))
res_sum_pattern1$DEGtype = factor(res_sum_pattern1$DEGtype,levels=c("InactivevsHC","HDAvsInactive"))
levels(res_sum_pattern1$DEGtype)=list(`Disease-state`="InactivevsHC",`Disease-activity`="HDAvsInactive")
col = pal_nejm("default")(2)[c(2,1)]

# w xaxis short
p3 = ggplot()+
     geom_bar(data=res_sum_pattern1, aes(x=subset,y=DEGnum,fill=DEGtype),
              stat="identity",position="dodge")+
     theme_classic()+
     scale_fill_manual(values = c("#0072B5FF","#BC3C29FF"))+
     theme(axis.text.x=element_text(colour="black",angle=45,hjust=1,vjust=1,size=16),
           axis.text.y=element_text(colour="black",size=16),
           axis.title.x=element_blank(),
           axis.title.y=element_text(colour="black",size=16),
           plot.title=element_blank(),
      　　　legend.position="right",
           legend.title=element_text(colour="black",size=16),
           legend.text=element_text(colour="black",size=16))+
     scale_x_discrete(labels= label)+
     scale_y_continuous(limits=c(0,6000))+
　　　labs(y="Number of DEGs")

pdf_3(paste0("res/stateactivity/DEGnum/",today(),"_COI_stateactivity_DEGnum_simpledodge_ordered_wxaxisyaxislegend.pdf"),h=2.9,w=10)
 plot(p3)
dev.off()

##############################################################################
# Extended: Spearman cor
##############################################################################
res_sum = fread_FT("res/stateactivity/211201_COI_27subsets_state_activity_logFC_logCPM_P_FDR.txt")

for(iii in 1:length(unique(res_sum$subset))){
 subset_tmp   = celltype_reordered_27[iii]
 df_tmp       = res_sum %>% filter(subset==subset_tmp)
 spearman     = cor.test(df_tmp$state_logFC,df_tmp$activity_logFC,method="spearman")
 res_tmp = data.frame(subset=subset_tmp,spe_cor=spearman$estimate,spe_P=spearman$p.value)
 if(iii==1){corres_sum=res_tmp}else{corres_sum=rbind(corres_sum,res_tmp)}
 }

corres_sum$spe_Bonf = p.adjust(corres_sum$spe_P,method="bonferroni")
corres_sum$subset = factor(corres_sum$subset,levels=levels(jaccard_sum_ordered$subset)) # ordered by jaccard group
corres_sum = corres_sum[order(corres_sum$subset),]

corres_sum = left_join(corres_sum,celltype_corresp%>%select(subset,lineage),by="subset")
write.table_FT_2(corres_sum,paste0("res/stateactivity/similarity/",today(),"_COI_27subset_state_and_activity_spearman.txt"))

corres_sum$lineage=factor(corres_sum$lineage,levels=c("CD4","CD8","NK","B","Monocyte","DC","Neutrophil"))
corres_sum$subset =factor(corres_sum$subset,levels=unique(corres_sum$subset))

# 220102 resize
p4 = ggplot()+
    geom_bar(data=corres_sum,aes(x=subset,y=spe_cor,fill=lineage),stat="identity")+
    theme_classic()+
    scale_fill_manual(values=col2)+
    theme(axis.text.x=element_text(colour="black",angle=45,hjust=1,vjust=1,size=16),
           axis.text.y=element_text(colour="black",size=16),
           axis.title.x=element_blank(),
           axis.title.y=element_text(colour="black",size=16),
           plot.title=element_blank(),
           legend.position="right",
           legend.title=element_text(colour="black",size=16),
           legend.text=element_text(colour="black",size=16))+
     scale_x_discrete(labels= label)+
    labs(y="Correlation\ncoefficient")

pdf_3(paste0("res/stateactivity/similarity/",today(),"_COI_27subset_state_and_activity_Spearman_wyaxislegend.pdf"),h=2.9,w=9.5)
 plot(p4)
dev.off()

##############################################################################
# Across celltype shareness and specificity
##############################################################################
list         = make_list("res/stateactivity","_DEGlist_FDR005.txt")
list$subset  = take_factor(list$FILE,3:4,"_")
list$type    = take_factor(list$FILE,5,"_")
list = list %>% filter(type%in%c("InactivevsHC","HDAvsInactive"))

for(kkk in 1:length(unique(list$type))){
   type_tmp  = unique(list$type)[kkk]
   list_tmp  = list %>% filter(type==type_tmp)

 for(iii in 1:nrow(list_tmp)){
  subset_tmp   = list_tmp$subset[iii]
  df_tmp         = read.table_FT(list_tmp$PATH[iii])
  
  df_tmp = df_tmp %>% select(1)
  df_tmp$subset  = 1
  colnames(df_tmp)[2] = subset_tmp
  if(iii==1){res_sum=df_tmp}else{res_sum=full_join(res_sum,df_tmp,by="genes")}
  }
  res_sum = res_sum %>% column_to_rownames("genes")
  
  subset_sum = apply(res_sum,1,function(x){sum(x,na.rm=T)}) %>% as.data.frame()
  colnames(subset_sum) = "celltype_num"
  

  res_sum_B   = res_sum %>% select(1:5)
  res_sum_CD4 = res_sum %>% select(6:14)
  res_sum_CD8 = res_sum %>% select(15:18)
  res_sum_DC  = res_sum %>% select(19:20)
  res_sum_Mono= res_sum %>% select(21:24)
  res_sum_Neu = res_sum %>% select(25:26)
  res_sum_NK  = res_sum %>% select(27)

  B_sum     = apply(res_sum_B,1,function(x){ifelse(sum(x,na.rm=T)>0,1,0)}) %>% as.data.frame()
  colnames(B_sum) ="B"
  CD4_sum   = apply(res_sum_CD4,1,function(x){ifelse(sum(x,na.rm=T)>0,1,0)}) %>% as.data.frame()
  colnames(CD4_sum) ="CD4"
  CD8_sum   = apply(res_sum_CD8,1,function(x){ifelse(sum(x,na.rm=T)>0,1,0)}) %>% as.data.frame()
  colnames(CD8_sum) ="CD8"
  DC_sum    = apply(res_sum_DC,1,function(x){ifelse(sum(x,na.rm=T)>0,1,0)}) %>% as.data.frame()
  colnames(DC_sum) ="DC"
  Mono_sum  = apply(res_sum_Mono,1,function(x){ifelse(sum(x,na.rm=T)>0,1,0)}) %>% as.data.frame()
  colnames(Mono_sum) ="Mono"
  Neu_sum   = apply(res_sum_Neu,1,function(x){ifelse(sum(x,na.rm=T)>0,1,0)}) %>% as.data.frame()
  colnames(Neu_sum) ="Neu"
  NK_sum   = apply(res_sum_NK,1,function(x){ifelse(sum(x,na.rm=T)>0,1,0)}) %>% as.data.frame()
  colnames(NK_sum) ="NK"

  lineage_sum = cbind(B_sum,CD4_sum) %>% cbind(.,CD8_sum) %>% cbind(.,DC_sum) %>% cbind(.,Mono_sum) %>% cbind(.,Neu_sum) %>% cbind(.,NK_sum)
  lineage_num = apply(lineage_sum,1,function(x){sum(x,na.rm=T)}) %>% as.data.frame()
  colnames(lineage_num) = "lineage_num"

  data_sum = cbind(subset_sum,lineage_num) %>% mutate(lineage_num2=ifelse(lineage_num>4,5,lineage_num))
  write.table_n_2(data_sum,"Gene",paste0("res/stateactivity/across_celltype/",today(),"_COI_27subsets_",type_tmp,"_DEGinfosum.txt"))
  
  data_sum$celltype_num =factor(data_sum$celltype_num,levels=sort(unique(data_sum$celltype_num)))
  data_sum$lineage_num2 =factor(data_sum$lineage_num2,levels=sort(unique(data_sum$lineage_num2)))

  prop_table=prop.table(table(data_sum$celltype_num,data_sum$lineage_num2)) %>% as.matrix()
  prop_table2=prop_table%>%as.data.frame()
  colnames(prop_table2)[1:2]=c("celltype_num","lineage_num")
  prop_table2$DEGtype=type_tmp

  if(kkk==1){prop_table_sum=prop_table2}else{prop_table_sum=rbind(prop_table_sum,prop_table2)}
}

write.table_FT_2(prop_table_sum,paste0("res/stateactivity/across_celltype/",today(),"_COI_27subsets_specificity_prop.txt"))

prop_table_sum$celltype_num =factor(prop_table_sum$celltype_num,levels=sort(unique(prop_table_sum$celltype_num)))
prop_table_sum$lineage_num =factor(prop_table_sum$lineage_num,levels=sort(unique(prop_table_sum$lineage_num)))

prop_table_sum$DEGtype = factor(prop_table_sum$DEGtype,levels=c("InactivevsHC","HDAvsInactive"))
facet1= as_labeller(c(`InactivevsHC`="Disease-state",`HDAvsInactive`="Disease-activity"))

col4=rev(brewer.pal(n=6,name="Greens"))[1:5]


p = ggplot()+
    geom_bar(data=prop_table_sum,aes(x=celltype_num,y=Freq,fill=factor(lineage_num,labels=c("1","2","3","4",">4"))),
             stat="identity",position="stack")+
    theme_classic()+
    facet_wrap(~DEGtype,labeller=facet1,ncol=2)+ 
    scale_fill_manual(values=col4) +
    theme(axis.text.x=element_text(colour="black",size=14),
           axis.text.y=element_text(colour="black",size=14),
           axis.title.x=element_text(colour="black",size=14),
           axis.title.y=element_text(colour="black",size=14),
           plot.title=element_blank(),
           strip.background=element_blank(),
           strip.text=element_text(colour="black",size=15),
           panel.background =element_rect(fill="transparent",colour="black",size=1),
           panel.grid = element_blank(),
           legend.position="right",
           legend.title=element_text(colour="black",size=14),
           legend.text=element_text(colour="black",size=14))+
     guides(fill=guide_legend(title="lineage"))+
     scale_x_discrete(breaks=c(1,5,10,15,20,25))+
     labs(x="Number of cell types shared by differential expression",y="Proportion of DEGs")

pdf_3(paste0("Figure/",today(),"_COI_27subsets_specificity_prop_wlegend.pdf"),h=3.5,w=7)
 plot(p)
dev.off()


###################################################################################################
# Cell type clustering jaccard + spearman
###################################################################################################
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(circlize))


res_sum = fread_FT("res/stateactivity/211201_COI_27subsets_state_activity_logFC_logCPM_P_FDR.txt")
res_sum$subset   = factor(res_sum$subset,levels=celltype_reordered_27)

# Jaccard
for(kkk in 1:length(celltype_reordered_27)){
    subset_1 = celltype_reordered_27[kkk]
    ref_tmp = res_sum %>% filter(subset==subset_1)
    ref_state = ref_tmp %>% filter(siggroup005%in%c("only_state","significant_in_both"))
    ref_act   = ref_tmp %>% filter(siggroup005%in%c("only_activity","significant_in_both"))

    for(mmm in 1:length(celltype_reordered_27)){
    	subset_2 = celltype_reordered_27[mmm]
    	data_tmp = res_sum %>% filter(subset==subset_2)
        data_state = data_tmp %>% filter(siggroup005%in%c("only_state","significant_in_both"))
        data_act   = data_tmp %>% filter(siggroup005%in%c("only_activity","significant_in_both"))
        concord_state = inner_join(ref_state,data_state,by="genes") %>% filter((state_logFC.x*state_logFC.y)>0)
        concord_act = inner_join(ref_act,data_act,by="genes") %>% filter((activity_logFC.x*activity_logFC.y)>0)

    	jaccard_state = data.frame(subset1=subset_1,subset2=subset_2,jaccard=(nrow(concord_state)/length(union(ref_state$genes,data_state$genes))))
        jaccard_act = data.frame(subset1=subset_1,subset2=subset_2,jaccard=(nrow(concord_act)/length(union(ref_act$genes,data_act$genes))))

    	if(mmm==1){state_sum=jaccard_state}else{state_sum=rbind(state_sum,jaccard_state)}
        if(mmm==1){act_sum=jaccard_act}else{act_sum=rbind(act_sum,jaccard_act)}
    }
    if(kkk==1){state_sum2=state_sum}else{state_sum2=rbind(state_sum2,state_sum)}
    if(kkk==1){act_sum2=act_sum}else{act_sum2=rbind(act_sum2,act_sum)}
}
    
state_sum3 = state_sum2 %>% pivot_wider(names_from="subset2",values_from="jaccard")
act_sum3 = act_sum2 %>% pivot_wider(names_from="subset2",values_from="jaccard")

write.table_FT_2(state_sum3,paste0("res/stateactivity/across_celltype/",today(),"_COI_state_27subset_Jaccard.txt"))
write.table_FT_2(act_sum3,paste0("res/stateactivity/across_celltype/",today(),"_COI_activity_27subset_Jaccard.txt"))

state_sum = fread_n("res/stateactivity/across_celltype/211201_COI_state_27subset_Jaccard.txt")
hc_state = hclust(as.dist(1-state_sum), method="ward.D2")
pdf_3(paste0("res/stateactivity/across_celltype/",today(),"_COI_state_27subset_Jaccard_cluster.pdf"),h=5,w=10)
 plot(hc_state)
dev.off()

act_sum = fread_n("res/stateactivity/across_celltype/211201_COI_activity_27subset_Jaccard.txt")
hc_act = hclust(as.dist(1-act_sum), method="ward.D2")
pdf_3(paste0("res/stateactivity/across_celltype/",today(),"_COI_activity_27subset_Jaccard_cluster.pdf"),h=5,w=10)
 plot(hc_act)
dev.off()

# rotate state_sum order maintaining cluster structures
labels_state = rev(hc_state$labels[hc_state$order])
labels_state = labels_state[c(4:6,2:1,3,14:13,11:12,7:10,18:15,22:21,19:20,23:27)]   
state_sum = state_sum[labels_state,labels_state]
act_sum = act_sum[labels_state,labels_state]

# Both integration
jac_both = state_sum
for(iii in 1:nrow(jac_both)){
  for(kkk in 1:ncol(jac_both)){
     if(kkk>iii){jac_both[iii,kkk]=act_sum[iii,kkk]}else{jac_both[iii,kkk]=state_sum[iii,kkk]}
 }
}
jac_both = jac_both %>% as.data.frame()
max(jac_both[jac_both!=max(jac_both)])
# [1] 0.5232558

# Heatmap
column_labels = structure(label,names=celltype_reordered_27)
col_fun = colorRamp2(seq(0,0.6,length=8), rev(brewer.pal(11,"RdYlBu")[1:8]))

labels_state2 = left_join(data.frame(subset=labels_state),celltype_corresp%>%select(subset,lineage),by="subset")
labels_state2$lineage = factor(labels_state2$lineage,levels=c("CD4","CD8","NK","B","Monocyte","DC","Neutrophil"))

ha = HeatmapAnnotation(lineage=labels_state2$lineage,col=list(lineage=col2),show_annotation_name=F,
					   height=unit(0.2,"cm"),simple_anno_size_adjust=T,show_legend=F,
					   annotation_legend_param=list(title_gp=gpar(fontsize=12),labels_gp=gpar(fontsize=12)))
ra =     rowAnnotation(lineage=labels_state2$lineage,col=list(lineage=col2),show_annotation_name=F,
					   width=unit(0.2,"cm"),simple_anno_size_adjust=T,show_legend=F,
					   annotation_legend_param=list(title_gp=gpar(fontsize=12),labels_gp=gpar(fontsize=12)))

p=Heatmap(jac_both,cluster_rows=FALSE,cluster_columns=FALSE,,col=col_fun,width=unit(12,"cm"),height=unit(12,"cm"),
                       row_names_gp=gpar(fontsize=12),row_labels=column_labels[rownames(jac_both)],
                       column_names_gp=gpar(fontsize=12),column_labels=column_labels[colnames(jac_both)],
                       bottom_annotation = ha,
                       right_annotation = ra,
                       heatmap_legend_param=list(title="Jaccard similarity index",at=c(0,0.2,0.4,0.6),
                       							 direction="horizontal",title_position="topcenter",
                       							 title_gp=gpar(fontsize=12,fontface="bold"),labels_gp=gpar(fontsize=12),
                       							 grid_height=unit(0.6,"cm"),legend_width=unit(4,"cm"))
                       )

pdf_3(paste0("res/stateactivity/across_celltype/",today(),"_COI_27subset_Jaccard_cluster_fig.pdf"),h=7,w=6)
 draw(p,heatmap_legend_side="bottom")
dev.off()

# Extended: Spearman
list1         = make_list("res/stateactivity/InactivevsHC/GLMresult","_GLMresult.txt",tag="state")
list1$subset  = take_factor(list1$FILE,3:4,"_")
list2         = make_list("res/stateactivity/HDAvsInactive/GLMresult","_GLMresult.txt",tag="activity")
list2$subset  = take_factor(list2$FILE,3:4,"_")

for(iii in 1:nrow(list1)){
    subset_tmp = list1$subset[iii]
　　 data_tmp   = read.table_FT(list1$PATH_state[iii]) %>% select(genes,logFC)
    colnames(data_tmp)[2] = subset_tmp
    if(iii==1){state_sum=data_tmp}else{state_sum=full_join(state_sum,data_tmp,by="genes")}
}
state_sum = state_sum %>% column_to_rownames("genes")
dim(state_sum)
# [1] 15967    27
spearman_state = cor(state_sum,use="pairwise.complete.obs",method="spearman") %>% as.data.frame()

for(iii in 1:nrow(list2)){
    subset_tmp = list2$subset[iii]
　　 data_tmp   = read.table_FT(list2$PATH_activity[iii]) %>% select(genes,logFC)
    colnames(data_tmp)[2] = subset_tmp
    if(iii==1){act_sum=data_tmp}else{act_sum=full_join(act_sum,data_tmp,by="genes")}
}
act_sum = act_sum %>% column_to_rownames("genes")
dim(act_sum)
# [1] 15967    27
spearman_act = cor(act_sum,use="pairwise.complete.obs",method="spearman") %>% as.data.frame()

write.table_n_2(spearman_state,"subset1",paste0("res/stateactivity/across_celltype/",today(),"_COI_state_27subset_spearman.txt"))
write.table_n_2(spearman_act,"subset1",paste0("res/stateactivity/across_celltype/",today(),"_COI_activity_27subset_spearman.txt"))

hc_state_spearman = hclust(as.dist(1-spearman_state), method="ward.D2")
pdf_3(paste0("res/stateactivity/across_celltype/",today(),"_COI_state_27subset_spearman_cluster.pdf"),h=5,w=10)
 plot(hc_state_spearman)
dev.off()

hc_act_spearman = hclust(as.dist(1-spearman_act), method="ward.D2")
pdf_3(paste0("res/stateactivity/across_celltype/",today(),"_COI_activity_27subset_spearman_cluster.pdf"),h=5,w=10)
 plot(hc_act_spearman)
dev.off()

# rotate act_sum order maintaining cluster structures
labels_act = hc_act_spearman$labels[hc_act_spearman$order]
labels_act = labels_act[c(2:3,1,8:9,4:7,10:20,26:27,23:25,21:22)]   
spearman_state = spearman_state[labels_act,labels_act]
spearman_act = spearman_act[labels_act,labels_act]

# Both integration
spearman_both = spearman_state
for(iii in 1:nrow(spearman_both)){
  for(kkk in 1:ncol(spearman_both)){
     if(kkk>iii){spearman_both[iii,kkk]=spearman_act[iii,kkk]}else{spearman_both[iii,kkk]=spearman_state[iii,kkk]}
 }
}
sort(spearman_both[spearman_both!=1],decreasing=T)
# 0.8261

col_fun = colorRamp2(seq(0.3,0.9,length=8), rev(brewer.pal(11,"RdYlBu")[1:8]))

labels_act2 = left_join(data.frame(subset=labels_act),celltype_corresp%>%select(subset,lineage),by="subset")
labels_act2$lineage = factor(labels_act2$lineage,levels=c("CD4","CD8","NK","B","Monocyte","DC","Neutrophil"))

ha = HeatmapAnnotation(lineage=labels_act2$lineage,col=list(lineage=col2),show_annotation_name=F,
					   height=unit(0.2,"cm"),simple_anno_size_adjust=T,show_legend=F,
					   annotation_legend_param=list(title_gp=gpar(fontsize=12),labels_gp=gpar(fontsize=12)))
ra =     rowAnnotation(lineage=labels_act2$lineage,col=list(lineage=col2),show_annotation_name=F,
					   width=unit(0.2,"cm"),simple_anno_size_adjust=T,show_legend=F,
					   annotation_legend_param=list(title_gp=gpar(fontsize=12),labels_gp=gpar(fontsize=12)))

p=Heatmap(spearman_both,cluster_rows=FALSE,cluster_columns=FALSE,,col=col_fun,width=unit(12,"cm"),height=unit(12,"cm"),
                       row_names_gp=gpar(fontsize=12),row_labels=column_labels[rownames(spearman_both)],
                       column_names_gp=gpar(fontsize=12),column_labels=column_labels[colnames(spearman_both)],
                       bottom_annotation = ha,
                       right_annotation = ra,
                       heatmap_legend_param=list(title="Correlation coefficient",at=c(0.3,0.5,0.7,0.9),
                       							 direction="horizontal",title_position="topcenter",
                       							 title_gp=gpar(fontsize=12,fontface="bold"),labels_gp=gpar(fontsize=12),
                       							 grid_height=unit(0.6,"cm"),legend_width=unit(4,"cm"))
                       )

pdf_3(paste0("res/stateactivity/across_celltype/",today(),"_COI_27subset_spearman_cluster_fig.pdf"),h=7,w=6)
 draw(p,heatmap_legend_side="bottom")
dev.off()


################################
# MMF DEGs Jaccard
################################
# DEGs
res_sum3=fread_FT("res/treatment/211130_COI_Tx_DEGnum.txt")
res_sum3$lineage = factor(res_sum3$lineage,levels=c("CD4","CD8","NK","B","Monocyte","DC","Neutrophil"))
res_sum3$subset  = factor(res_sum3$subset,levels=celltype_reordered_27)
res_sum3$DEGnum  = as.numeric(res_sum3$DEGnum)
# MMF only
DEG_sum = res_sum3 %>% filter(DEGtype=="MMF")

# Jaccard with activity
list1 = make_list("res/stateactivity/HDAvsInactive/GLMresult","_HDAvsInactive_GLMresult.txt",tag="activity")
list1$subset = take_factor(list1$FILE_activity,3:4,"_")
list2 = make_list("res/treatment/MMF/GLMresult","_MMF_GLMresult.txt",tag="MMF")
list2$subset = take_factor(list2$FILE_MMF,3:4,"_")

list_tmp_2 = left_join(list1,list2,by="subset") 

for(iii in 1:nrow(list_tmp_2)){
    subset_tmp         = list_tmp_2$subset[iii]
    df_tmp_m           = fread_FT(list_tmp_2$PATH_MMF[iii]) %>% select(genes,logCPM,logFC,PValue,FDR)
    colnames(df_tmp_m)[2:5] = paste0("MMF_",colnames(df_tmp_m)[2:5])

    df_tmp_a           = read.table_FT(list_tmp_2$PATH_activity[iii]) %>% select(genes,logCPM,logFC,PValue,FDR)
    colnames(df_tmp_a)[2:5] = paste0("activity_",colnames(df_tmp_a)[2:5])
       
    df_tmp = inner_join(df_tmp_a,df_tmp_m,by="genes")
    df_tmp$subset =  subset_tmp

    if(iii==1){res_sum=df_tmp}else{res_sum=rbind(res_sum,df_tmp)}
  }

  res_sum$siggroup005 = ifelse(res_sum$activity_FDR<0.05&res_sum$MMF_FDR<0.05&(res_sum$activity_logFC*res_sum$MMF_logFC)<0,"significant_in_both",
                         ifelse(res_sum$activity_FDR<0.05&res_sum$MMF_FDR<0.05&(res_sum$activity_logFC*res_sum$MMF_logFC)>0,"discordant",
                          ifelse(res_sum$activity_FDR<0.05,"only_activity",
                          ifelse(res_sum$MMF_FDR<0.05,"only_MMF","neither"))))
  write.table_FT_2(res_sum,paste0("res/treatment/",today(),"_COI_SLE_activity_MMF_logFC_logCPM_P_FDR.txt"))

for(iii in 1:length(celltype_reordered_27)){
    subset_tmp = celltype_reordered_27[iii]
    res_tmp = res_sum %>% filter(subset==subset_tmp)
    MMF_tmp = res_tmp %>% filter(siggroup005=="only_MMF")
    activity_tmp = res_tmp %>% filter(siggroup005=="only_activity")
    both_tmp  = res_tmp %>% filter(siggroup005=="significant_in_both")
    dis_tmp   = res_tmp %>% filter(siggroup005=="discordant")
    jaccard_tmp = data.frame(subset=subset_tmp,MMF=nrow(MMF_tmp),activity=nrow(activity_tmp),both=nrow(both_tmp),
                             discordant=nrow(dis_tmp),union=nrow(MMF_tmp)+nrow(activity_tmp)+nrow(both_tmp),
                             all=nrow(MMF_tmp)+nrow(activity_tmp)+nrow(both_tmp)+nrow(dis_tmp),
                             MMF_P=nrow(MMF_tmp)/(nrow(MMF_tmp)+nrow(activity_tmp)+nrow(both_tmp)),
                             activity_P=nrow(activity_tmp)/(nrow(MMF_tmp)+nrow(activity_tmp)+nrow(both_tmp)),
                             jaccard=nrow(both_tmp)/(nrow(MMF_tmp)+nrow(activity_tmp)+nrow(both_tmp)))
    if(iii==1){jaccard_sum=jaccard_tmp}else{jaccard_sum=rbind(jaccard_sum,jaccard_tmp)}
    }
write.table_FT_2(jaccard_sum,paste0("res/treatment/",today(),"_COI_MMF_activity_jaccard.txt"))

all_sum=left_join(DEG_sum,jaccard_sum,by="subset")
all_sum$jaccard2=all_sum$jaccard*10000
all_sum=all_sum %>% select(subset,DEGnum,jaccard2) %>% 
                    pivot_longer(col=-subset,names_to="parameter",values_to="value")
all_sum$subset  = factor(all_sum$subset,levels=celltype_reordered_27)
all_sum$parameter  = factor(all_sum$parameter,levels=c("DEGnum","jaccard2"))
levels(all_sum$parameter)=list(`Number of DEGs`="DEGnum",`Jaccard similarity index`="jaccard2")

pal_nejm("default")(6)
#[1] "#BC3C29FF" "#0072B5FF" "#E18727FF" "#20854EFF" "#7876B1FF" "#6F99ADFF"

p = ggplot()+
     geom_bar(data=all_sum, aes(x=subset,y=value,fill=parameter),stat="identity",position="dodge")+
     theme_classic()+
     scale_fill_manual(values=c("#7876B1FF","#E18727FF"))+
     theme(axis.text.x=element_text(colour="black",angle=45,hjust=1,vjust=1,size=16),
           axis.text.y=element_text(colour="black",size=16),
           axis.title.x=element_blank(),
           axis.title.y=element_text(colour="black",size=16),
           plot.title=element_blank(),
           legend.position="right",
           legend.title=element_blank(),
           legend.text=element_text(colour="black",size=16))+
     scale_x_discrete(labels= label)+
     scale_y_continuous(name="Number of DEGs",
                        sec.axis=sec_axis(trans=~.*(1/10000), name="",breaks=c(0,0.1,0.2)))

 pdf_3(paste0("res/treatment/",today(),"_COI_MMF_DEGnum_Jaccard.pdf"),h=3,w=11.5)
  plot(p)
 dev.off()


######################
# For Supple
######################
res_sum = fread_FT("res/stateactivity/211201_COI_27subsets_state_activity_logFC_logCPM_P_FDR.txt") %>% 
          left_join(.,celltype_corresp[,c(1,2)],by="subset") 
res_sum$subset=factor(res_sum$subset,levels=celltype_reordered_27)
res_sum = res_sum[order(res_sum$genes),]
res_sum = res_sum[order(res_sum$subset),]

# state
res_sum_state = res_sum %>% select(label,genes,state_logCPM,state_logFC,state_PValue,state_FDR) 
res_sum_state$state_logCPM = formatC(res_sum_state$state_logCPM,digits=2)
res_sum_state$state_logFC  = formatC(res_sum_state$state_logFC,digits=2)
res_sum_state$state_PValue = formatC(res_sum_state$state_PValue,digits=2)
res_sum_state$state_FDR    = formatC(res_sum_state$state_FDR,digits=2)
colnames(res_sum_state) =c("Cell type","Gene","logCPM","logFC","Pvalue","FDR")
write.table_FT_2(res_sum_state,paste0(today(),"_COI_state_logFC_logCPM_P_FDR_Supple.txt"))

# activity
res_sum_activity = res_sum %>% select(label,genes,activity_logCPM,activity_logFC,activity_PValue,activity_FDR) 
res_sum_activity$activity_logCPM = formatC(res_sum_activity$activity_logCPM,digits=2)
res_sum_activity$activity_logFC  = formatC(res_sum_activity$activity_logFC,digits=2)
res_sum_activity$activity_PValue = formatC(res_sum_activity$activity_PValue,digits=2)
res_sum_activity$activity_FDR    = formatC(res_sum_activity$activity_FDR,digits=2)
colnames(res_sum_activity) =c("Cell type","Gene","logCPM","logFC","Pvalue","FDR")
write.table_FT_2(res_sum_activity,paste0(today(),"_COI_activity_logFC_logCPM_P_FDR_Supple.txt"))

# Tx
list         = make_list("res/treatment","_GLMresult.txt")
list$subset  = take_factor(list$FILE,3:4,"_")
list$type    = take_factor(list$FILE,5,"_")

for(iii in 1:nrow(list)){
    subset_tmp    = list$subset[iii]
    type_tmp      = list$type[iii]
    df_tmp        = fread_FT(list$PATH[iii]) %>% select(genes,logFC,PValue,FDR)
    df_tmp$type   =  type_tmp
    df_tmp$subset =  subset_tmp
    if(iii==1){res_sum=df_tmp}else{res_sum=rbind(res_sum,df_tmp)}
  }

res_sum = res_sum %>% left_join(.,celltype_corresp[,c(1,2)],by="subset") 
res_sum$subset=factor(res_sum$subset,levels=celltype_reordered_27)
res_sum$type=factor(res_sum$type,levels=c("MMF","HCQ","TAC"))
res_sum = res_sum[order(res_sum$genes),]
res_sum = res_sum[order(res_sum$subset),]
res_sum = res_sum[order(res_sum$type),]

res_sum = res_sum %>% select(type,label,genes,logFC,PValue,FDR) 
res_sum$logFC  = formatC(res_sum$logFC,digits=2)
res_sum$PValue = formatC(res_sum$PValue,digits=2)
res_sum$FDR    = formatC(res_sum$FDR,digits=2)
colnames(res_sum) =c("Agent","Cell type","Gene","logFC","Pvalue","FDR")
write.table_FT_2(res_sum,paste0(today(),"_COI_Tx_logFC_P_FDR_Supple.txt"))


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
[1] grid      stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] msigdbr_7.4.1        wesanderson_0.3.6    circlize_0.4.13     
 [4] ComplexHeatmap_2.6.2 RColorBrewer_1.1-2   ggsci_2.9           
 [7] data.table_1.14.0    forcats_0.5.1        stringr_1.4.0       
[10] dplyr_1.0.7          purrr_0.3.4          readr_1.4.0         
[13] tidyr_1.1.3          tibble_3.1.2         ggplot2_3.3.5       
[16] tidyverse_1.3.1     

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
[55] parallel_4.0.2      babelgene_21.4      clue_0.3-59        
[58] colorspace_2.0-2    cluster_2.1.2       rvest_1.0.0        
[61] haven_2.4.1 































