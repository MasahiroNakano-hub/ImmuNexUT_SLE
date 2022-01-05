#---------------------------------------#
# External dataset validation
# edgeR
# 211223 Nakano
#---------------------------------------#

source("my.source.R")
options(stringsAsFactors=F)

#
################################################################################################
############### edgeR
#### Example: One cell type in cohort 1
################################################################################################
LIST="tmp_job_list/211224_GenomeRes2021_edgeR_list.txt"
##  task_id=1

list       = fread_FT(LIST)
subset_tmp = list[task_id,3]
count_path = list[task_id,1]
batch_path = list[task_id,4]
cohort_tmp = list[task_id,6]

paste0("cohort_tmp : ",cohort_tmp) %>% print()
paste0("subset_tmp : ",subset_tmp) %>% print()
paste0("count_path : ",count_path) %>% print()
paste0("batch_path : ",batch_path) %>% print()

###########################################################################

out_f  = paste0("~/ws/2021/211223_COI_SLE_HC_replication/res/",cohort_tmp,"/stateactivity/")
out_ff = paste0(today(),"_",cohort_tmp,"_",subset_tmp,"_")

###################
################################

count_tmp  = fread_n(count_path)
cond_tmp   = fread_FT(batch_path)

chrY_gene = fread_FF("data_ref/hg38_chrY_gene_list.txt")[,1]

count_tmp = count_tmp[!is.element(rownames(count_tmp),chrY_gene),]

## Check
if(ncol(count_tmp)==nrow(cond_tmp)){"sample and cond were same num"}else{"please check condition list"}

cond_tmp_2  = cond_tmp[is.element(cond_tmp$Sample,colnames(count_tmp)),]
cond_tmp_2=cond_tmp_2[!is.na(cond_tmp_2$Activity),]
count_tmp_2 = count_tmp[,cond_tmp_2$Sample]

## Check
print("## final use data number ##")
table_freq(cond_tmp_2$Activity)
print("## final use data number ##")

#REFERNCE = unique(cond_tmp_2$batch) %>% grep(TARGET,.,invert=T,value=T)
cond_tmp_2$Activity  = factor(cond_tmp_2$Activity, levels=sort(unique(cond_tmp_2$Activity)))

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
design_1     = model.matrix(~ 0+Activity+Age+Gender , data= cond_tmp_2 )

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

print(paste0(subset_tmp," InactivevsHC was done"))

#####################################
#HDA vs Inactive
#####################################

con = makeContrasts(Activity4HDA - Activity1Inactive, levels=design_1)
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


################################################################################################
# Differential expression sign test
################################################################################################

###############
# State sig
###############
# AMED
list1 = make_list("res/AMED/stateactivity/InactivevsHC/GLMresult","_GLMresult.txt")
list1$subset = take_factor(list1$FILE,3:4,"_")
list1$cohort = take_factor(list1$FILE,2,"_")
# GenomeRes2021
list2 = make_list("res/GenomeRes2021/stateactivity/InactivevsHC/GLMresult","_GLMresult.txt")
list2$subset = take_factor(list2$FILE,3,"_")
list2$cohort = take_factor(list2$FILE,2,"_")
list=rbind(list1,list2)

# Discovery
listd=make_list("~/ws/2021/211129_COI_SLE_HC_run27_edgeR/res/stateactivity/InactivevsHC/DEGlist005","_DEGlist_FDR005.txt")
listd$subset = take_factor(listd$FILE,3:4,"_")

for(iii in 1:length(unique(list$cohort))){
    cohort_tmp=unique(list$cohort)[iii]
    list_tmp=list%>%filter(cohort==cohort_tmp)

    for(kkk in 1:nrow(list_tmp)){
        subset_tmp_r = list_tmp$subset[kkk]
        data_tmp_r   = read.table_FT(list_tmp$PATH[kkk]) %>% select(genes,logFC)
        colnames(data_tmp_r)[2] = "r_logFC"
        
      for(mmm in 1:nrow(listd)){
        subset_tmp_d = listd$subset[mmm]
        data_tmp_d   = read.table_FT(listd$PATH[mmm]) %>%filter(FDR<0.05)　%>% select(genes,logFC)
        colnames(data_tmp_d)[2] = "d_logFC"
        data_tmp_d$d_state_sign = sign_pick(data_tmp_d$d_logFC)
        
        data_tmp = inner_join(data_tmp_d,data_tmp_r,by="genes")
        data_tmp$r_logFC2 = data_tmp$d_state_sign * data_tmp$r_logFC

        dir_same = length(which(data_tmp$r_logFC2>0))
        dir_same_prop = dir_same/nrow(data_tmp)
        binom = binom.test(c(dir_same,nrow(data_tmp)-dir_same),alternative="greater")
        binom_logP = -log10(binom$p.value)

        res_tmp= data.frame(cohort=cohort_tmp,discovery=subset_tmp_d,replication=subset_tmp_r,binom_prop=dir_same_prop,binom_p=binom$p.value,binom_logP=binom_logP)
        if(mmm==1){res_sum=res_tmp}else{res_sum=rbind(res_sum,res_tmp)}
    }
    if(kkk==1){res_sum2=res_sum}else{res_sum2=rbind(res_sum2,res_sum)}
  }

  write.table_FT_2(res_sum2,paste0("res/",today(),"_",cohort_tmp,"_state_onebinom_res.txt"))
}

### representative plot
# ImmuNexUT USMB vs AMED USMB, pDC 
for(kkk in c(3,16)){
    subset_tmp_r = list1$subset[kkk]
    data_tmp_r   = read.table_FT(list1$PATH[kkk]) %>% select(genes,logFC)
    colnames(data_tmp_r)[2] = "r_logFC"
        
    mmm=3
    subset_tmp_d = listd$subset[mmm]
    data_tmp_d   = read.table_FT(listd$PATH[mmm]) %>%filter(FDR<0.05)　%>% select(genes,logFC)
    colnames(data_tmp_d)[2] = "d_logFC"
        
    data_tmp = inner_join(data_tmp_d,data_tmp_r,by="genes")
    data_tmp$subset_d = subset_tmp_d
    data_tmp$subset_r = subset_tmp_r

    if(kkk==3){res_sum=data_tmp}else{res_sum=rbind(res_sum,data_tmp)}
}
    
res_sum_max5 = res_sum
res_sum_max5$r_logFC[res_sum_max5$r_logFC>5] = 5 
res_sum_max5$r_logFC[res_sum_max5$r_logFC< -5] = -5 
res_sum_max5$d_logFC[res_sum_max5$d_logFC>5] = 5 
res_sum_max5$d_logFC[res_sum_max5$d_logFC< -5] = -5 

res_sum_max5$subset_r = factor(res_sum_max5$subset_r,levels=c("A03_UnswMB","D02_pDC"))
label = as_labeller(c(`A03_UnswMB`="Cohort2 USM B",`D02_pDC`="Cohort2 pDC"))

p =  ggplot(res_sum_max5,aes(x=d_logFC,y=r_logFC))+
          geom_point(size=0.2,color="#0072B5FF")+
          theme_classic()+
          facet_wrap(~ subset_r, labeller=label,nrow = 1) + 
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
                legend.text=element_blank())+
          guides(colour=guide_legend(override.aes=list(size=1.5)))+
          scale_x_continuous(breaks=c(-4,-2,0,2,4),limits=c(-5,5))+
          scale_y_continuous(breaks=c(-4,-2,0,2,4),limits=c(-5,5))

png(paste0("res/",today(),"_AMED_state_logFC_rep.png"),height=1.9,width=3.8,units="in",res=300)
 plot(p)
dev.off()

# ImmuNexUT CL Mono vs GenomeRes2021 cMo, T 
for(kkk in c(3,6)){
    subset_tmp_r = list2$subset[kkk]
    data_tmp_r   = read.table_FT(list2$PATH[kkk]) %>% select(genes,logFC)
    colnames(data_tmp_r)[2] = "r_logFC"
        
    mmm=22
    subset_tmp_d = listd$subset[mmm]
    data_tmp_d   = read.table_FT(listd$PATH[mmm]) %>%filter(FDR<0.05)　%>% select(genes,logFC)
    colnames(data_tmp_d)[2] = "d_logFC"
        
    data_tmp = inner_join(data_tmp_d,data_tmp_r,by="genes")
    data_tmp$subset_d = subset_tmp_d
    data_tmp$subset_r = subset_tmp_r

    if(kkk==3){res_sum=data_tmp}else{res_sum=rbind(res_sum,data_tmp)}
}
    
res_sum_max5 = res_sum
res_sum_max5$r_logFC[res_sum_max5$r_logFC>5] = 5 
res_sum_max5$r_logFC[res_sum_max5$r_logFC< -5] = -5 
res_sum_max5$d_logFC[res_sum_max5$d_logFC>5] = 5 
res_sum_max5$d_logFC[res_sum_max5$d_logFC< -5] = -5 

res_sum_max5$subset_r = factor(res_sum_max5$subset_r,levels=c("cMo","T"))
label = as_labeller(c(`cMo`="Cohort1 CL Mono",`T`="Cohort1 T"))

p =  ggplot(res_sum_max5,aes(x=d_logFC,y=r_logFC))+
          geom_point(size=0.2,color="#0072B5FF")+
          theme_classic()+
          facet_wrap(~ subset_r, labeller=label,nrow = 1) + 
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
                legend.text=element_blank())+
          guides(colour=guide_legend(override.aes=list(size=1.5)))+
          scale_x_continuous(breaks=c(-4,-2,0,2,4),limits=c(-5,5))+
          scale_y_continuous(breaks=c(-4,-2,0,2,4),limits=c(-5,5))

png(paste0("res/",today(),"_GenomeRes2021_state_logFC_rep_v2.png"),height=1.9,width=3.8,units="in",res=300)
 plot(p)
dev.off()

###############
# Activity sig
###############
# GenomeRes2021
list2 = make_list("res/GenomeRes2021/stateactivity/HDAvsInactive/GLMresult","_GLMresult.txt")
list2$subset = take_factor(list2$FILE,3,"_")
list2$cohort = take_factor(list2$FILE,2,"_")
# CommuBiol2021
list3 = make_list("res/CommuBiol2021/stateactivity/HDAvsInactive_woRace/GLMresult","_GLMresult.txt")
list3$subset = take_factor(list3$FILE,3,"_")
list3$cohort = take_factor(list3$FILE,2,"_")

list=rbind(list2,list3)

# Discovery
listd=make_list("~/ws/2021/211129_COI_SLE_HC_run27_edgeR/res/stateactivity/HDAvsInactive/DEGlist005","_DEGlist_FDR005.txt")
listd$subset = take_factor(listd$FILE,3:4,"_")

### representative plot
# ImmuNexUT CL Mono vs GenomeRes2021 cMo, T 
for(kkk in c(3,6)){
    subset_tmp_r = list2$subset[kkk]
    data_tmp_r   = read.table_FT(list2$PATH[kkk]) %>% select(genes,logFC)
    colnames(data_tmp_r)[2] = "r_logFC"
        
    mmm=22
    subset_tmp_d = listd$subset[mmm]
    data_tmp_d   = read.table_FT(listd$PATH[mmm]) %>%filter(FDR<0.05)　%>% select(genes,logFC)
    colnames(data_tmp_d)[2] = "d_logFC"
        
    data_tmp = inner_join(data_tmp_d,data_tmp_r,by="genes")
    data_tmp$subset_d = subset_tmp_d
    data_tmp$subset_r = subset_tmp_r

    if(kkk==3){res_sum=data_tmp}else{res_sum=rbind(res_sum,data_tmp)}
}
    
res_sum_max5 = res_sum
res_sum_max5$r_logFC[res_sum_max5$r_logFC>5] = 5 
res_sum_max5$r_logFC[res_sum_max5$r_logFC< -5] = -5 
res_sum_max5$d_logFC[res_sum_max5$d_logFC>5] = 5 
res_sum_max5$d_logFC[res_sum_max5$d_logFC< -5] = -5 

res_sum_max5$subset_r = factor(res_sum_max5$subset_r,levels=c("cMo","T"))
label = as_labeller(c(`cMo`="Cohort1 CL Mono",`T`="Cohort1 T"))

p =  ggplot(res_sum_max5,aes(x=d_logFC,y=r_logFC))+
          geom_point(size=0.2,color="#BC3C29FF")+
          theme_classic()+
          facet_wrap(~ subset_r, labeller=label,nrow = 1) + 
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
                legend.text=element_blank())+
          guides(colour=guide_legend(override.aes=list(size=1.5)))+
          scale_x_continuous(breaks=c(-4,-2,0,2,4),limits=c(-5,5))+
          scale_y_continuous(breaks=c(-4,-2,0,2,4),limits=c(-5,5))

png(paste0("res/",today(),"_GenomeRes2021_activity_logFC_rep_v2.png"),height=1.9,width=3.8,units="in",res=300)
 plot(p)
dev.off()

# ImmuNexUT NK vs CommuBiol2021 NK, Mono 
for(kkk in c(3,4)){
    subset_tmp_r = list3$subset[kkk]
    data_tmp_r   = read.table_FT(list3$PATH[kkk]) %>% select(genes,logFC)
    colnames(data_tmp_r)[2] = "r_logFC"
        
    mmm=27
    subset_tmp_d = listd$subset[mmm]
    data_tmp_d   = read.table_FT(listd$PATH[mmm]) %>%filter(FDR<0.05)　%>% select(genes,logFC)
    colnames(data_tmp_d)[2] = "d_logFC"
        
    data_tmp = inner_join(data_tmp_d,data_tmp_r,by="genes")
    data_tmp$subset_d = subset_tmp_d
    data_tmp$subset_r = subset_tmp_r

    if(kkk==3){res_sum=data_tmp}else{res_sum=rbind(res_sum,data_tmp)}
}
    
res_sum_max5 = res_sum
res_sum_max5$r_logFC[res_sum_max5$r_logFC>5] = 5 
res_sum_max5$r_logFC[res_sum_max5$r_logFC< -5] = -5 
res_sum_max5$d_logFC[res_sum_max5$d_logFC>5] = 5 
res_sum_max5$d_logFC[res_sum_max5$d_logFC< -5] = -5 

res_sum_max5$subset_r = factor(res_sum_max5$subset_r,levels=c("NK","Mono"))
label = as_labeller(c(`NK`="Cohort3 NK",`Mono`="Cohort3 Mono"))

p =  ggplot(res_sum_max5,aes(x=d_logFC,y=r_logFC))+
          geom_point(size=0.2,color="#BC3C29FF")+
          theme_classic()+
          facet_wrap(~ subset_r, labeller=label,nrow = 1) + 
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
                legend.text=element_blank())+
          guides(colour=guide_legend(override.aes=list(size=1.5)))+
          scale_x_continuous(breaks=c(-4,-2,0,2,4),limits=c(-5,5))+
          scale_y_continuous(breaks=c(-4,-2,0,2,4),limits=c(-5,5))

png(paste0("res/",today(),"_CommuBiol2021_activity_logFC_rep.png"),height=1.9,width=3.8,units="in",res=300)
 plot(p)
dev.off()

###############
# -logP Heatmap
###############
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(RColorBrewer))
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
celltype_corresp$lineage=factor(celltype_corresp$lineage,levels=c("CD4","CD8","NK","B","Monocyte","DC","Neutrophil")) 

dis_rep_corresp=fread_FT("data_ref/subset_replication_corresp.txt")
#####

##### State, Cohort1
celltype_reordered_6=c("T","B","cMo","cDC","pDC","PMN")
colnames6 = c("T","B","CL Mono","mDC","pDC","Neu")
row_labels6 = structure(colnames6,names=celltype_reordered_6)
celltype_corresp6 = dis_rep_corresp %>% select(I_subset,G_subset,lineage)%>%filter(!is.na(G_subset))
celltype_corresp6$lineage=factor(celltype_corresp6$lineage,levels=c("CD4","B","Monocyte","DC","Neutrophil"))
                 
state1 = fread_FT("res/211224_GenomeRes2021_state_onebinom_res.txt") %>% filter(discovery%in%celltype_corresp6$I_subset) %>%
      select(discovery,replication,binom_logP) %>% pivot_wider(names_from="discovery",values_from="binom_logP") %>% column_to_rownames("replication") 

state1 = state1[celltype_reordered_6,c("B05_NCD4","A01_NaiB","E02_CD16nMo","D01_mDC","D02_pDC","F01_Neu")]
column_labels6=column_labels[c("B05_NCD4","A01_NaiB","E02_CD16nMo","D01_mDC","D02_pDC","F01_Neu")]
# Corresponding p<10e-34

bonf_threshold = -log10(0.05/(ncol(state1)*nrow(state1)))

col_fun = colorRamp2(c(bonf_threshold,50), c("#FFFFFF","#0072B5FF"))

ha = HeatmapAnnotation(lineage=celltype_corresp6$lineage,col=list(lineage=col2),show_annotation_name=F,
                    height=unit(0.2,"cm"),simple_anno_size_adjust=T,show_legend=T,
                    annotation_legend_param=list(title_gp=gpar(fontsize=12),labels_gp=gpar(fontsize=12)))

ra = rowAnnotation(lineage=celltype_corresp6$lineage,col=list(lineage=col2),show_annotation_name=F,
                    width=unit(0.2,"cm"),simple_anno_size_adjust=T,show_legend=F,
                    annotation_legend_param=list(title_gp=gpar(fontsize=12),labels_gp=gpar(fontsize=12)))

p=Heatmap(state1,cluster_rows=FALSE,cluster_columns=FALSE,,col=col_fun,width=unit(2.5,"cm"),height=unit(2.5,"cm"),
                       row_names_gp=gpar(fontsize=12),row_labels=row_labels6,
                       column_names_gp=gpar(fontsize=12),column_labels=column_labels6,
                       right_annotation = ra,
                       bottom_annotation = ha,
                       heatmap_legend_param=list(title="-log10(P)",at=c(2.5,25,50),
                                            direction="horizontal",title_position="topcenter",
                                            grid_width=unit(2,"cm"),legend_height=unit(0.4,"cm"),
                                            title_gp=gpar(fontsize=12,fontface="bold"),labels_gp=gpar(fontsize=12)))

pdf_3(paste0("res/",today(),"_GenomeRes2021_state_onebinom_res.pdf"),h=3,w=3.5)
 draw(p,heatmap_legend_side="bottom", annotation_legend_side="right")
dev.off()

##### State, Cohort2
celltype_reordered_19=c("B05_NCD4","B06_MCD4","B01_Th1","B02_Th2","B03_TH17","B04_Tfh","B07_aTreg","C01_NCD8","C02_MCD8","G01_NK","A01_NaiB","A03_UnswMB","A02_SwiMB","A04_DNB","A05_PB","E02_CD16nMo","E01_CD16pMo","D01_mDC","D02_pDC")
colnames19 = c("Naive CD4","Mem CD4","Th1","Th2","Th17","Tfh","Fr. II eTreg","Naive CD8","Mem CD8","NK","Naive B","USM B","SM B","DN B","Plasmablast","CL Mono","CD16p Mono","mDC","pDC")
row_labels19 = structure(colnames19,names=celltype_reordered_19)
celltype_corresp19 = dis_rep_corresp %>% select(I_subset,A_subset,lineage)%>%filter(!is.na(A_subset))
celltype_corresp19$lineage=factor(celltype_corresp19$lineage,levels=c("CD4","CD8","NK","B","Monocyte","DC"))       

state2 = fread_FT("res/211224_AMED_state_onebinom_res.txt") %>% 
      select(discovery,replication,binom_logP) %>% pivot_wider(names_from="discovery",values_from="binom_logP") %>% column_to_rownames("replication") 

state2 = state2[celltype_reordered_19,c("B05_NCD4","B06_MCD4","B01_Th1","B02_Th2","B03_TH17","B04_Tfh","B07_aTreg","C01_NCD8","C05_EmCD8","G01_NK","A01_NaiB","A03_UnswMB","A02_SwiMB","A04_DNB","A05_PB","E02_CD16nMo","E01_CD16pMo","D01_mDC","D02_pDC")]
state2[state2=="Inf"] = 300
column_labels19=column_labels[c("B05_NCD4","B06_MCD4","B01_Th1","B02_Th2","B03_TH17","B04_Tfh","B07_aTreg","C01_NCD8","C05_EmCD8","G01_NK","A01_NaiB","A03_UnswMB","A02_SwiMB","A04_DNB","A05_PB","E02_CD16nMo","E01_CD16pMo","D01_mDC","D02_pDC")]
# Corresponding p<10e-76

bonf_threshold = -log10(0.05/(ncol(state2)*nrow(state2)))

col_fun = colorRamp2(c(bonf_threshold,200), c("#FFFFFF","#0072B5FF"))

ha = HeatmapAnnotation(lineage=celltype_corresp19$lineage,col=list(lineage=col2),show_annotation_name=F,
                    height=unit(0.2,"cm"),simple_anno_size_adjust=T,show_legend=T,
                    annotation_legend_param=list(title_gp=gpar(fontsize=12),labels_gp=gpar(fontsize=12)))

ra = rowAnnotation(lineage=celltype_corresp19$lineage,col=list(lineage=col2),show_annotation_name=F,
                    width=unit(0.2,"cm"),simple_anno_size_adjust=T,show_legend=F,
                    annotation_legend_param=list(title_gp=gpar(fontsize=12),labels_gp=gpar(fontsize=12)))

p=Heatmap(state2,cluster_rows=FALSE,cluster_columns=FALSE,,col=col_fun,width=unit(7.5,"cm"),height=unit(7.5,"cm"),
                       row_names_gp=gpar(fontsize=12),row_labels=row_labels19,
                       column_names_gp=gpar(fontsize=12),column_labels=column_labels19,
                       right_annotation = ra,
                       bottom_annotation = ha,
                       heatmap_legend_param=list(title="-log10(P)",at=c(3.5,100,200),
                                            direction="horizontal",title_position="topcenter",
                                            grid_width=unit(2,"cm"),legend_height=unit(0.4,"cm"),
                                            title_gp=gpar(fontsize=12,fontface="bold"),labels_gp=gpar(fontsize=12)))

pdf_3(paste0("res/",today(),"_AMED_state_onebinom_res.pdf"),h=5,w=5.5)
 draw(p,heatmap_legend_side="bottom")
dev.off()

##### Activity, Cohort1
act1 = fread_FT("res/211224_GenomeRes2021_activity_onebinom_res.txt") %>% 
      select(discovery,replication,binom_logP) %>% pivot_wider(names_from="discovery",values_from="binom_logP") %>% column_to_rownames("replication") 

act1 = act1[celltype_reordered_6,c("B05_NCD4","A01_NaiB","E02_CD16nMo","D01_mDC","D02_pDC","F01_Neu")]
# Corresponding p<10e-8

bonf_threshold = -log10(0.05/(ncol(act1)*nrow(act1)))

col_fun = colorRamp2(c(bonf_threshold,50), c("#FFFFFF","#BC3C29FF"))

ha = HeatmapAnnotation(lineage=celltype_corresp6$lineage,col=list(lineage=col2),show_annotation_name=F,
                    height=unit(0.2,"cm"),simple_anno_size_adjust=T,show_legend=T,
                    annotation_legend_param=list(title_gp=gpar(fontsize=12),labels_gp=gpar(fontsize=12)))

ra = rowAnnotation(lineage=celltype_corresp6$lineage,col=list(lineage=col2),show_annotation_name=F,
                    width=unit(0.2,"cm"),simple_anno_size_adjust=T,show_legend=F,
                    annotation_legend_param=list(title_gp=gpar(fontsize=12),labels_gp=gpar(fontsize=12)))

p=Heatmap(act1,cluster_rows=FALSE,cluster_columns=FALSE,,col=col_fun,width=unit(2.5,"cm"),height=unit(2.5,"cm"),
                       row_names_gp=gpar(fontsize=12),row_labels=row_labels6,
                       column_names_gp=gpar(fontsize=12),column_labels=column_labels6,
                       right_annotation = ra,
                       bottom_annotation = ha,
                       heatmap_legend_param=list(title="-log10(P)",at=c(2.5,25,50),
                                            direction="horizontal",title_position="topcenter",
                                            grid_width=unit(2,"cm"),legend_height=unit(0.4,"cm"),
                                            title_gp=gpar(fontsize=12,fontface="bold"),labels_gp=gpar(fontsize=12)))

pdf_3(paste0("res/",today(),"_GenomeRes2021_activity_onebinom_res.pdf"),h=3,w=3.5)
 draw(p,heatmap_legend_side="bottom", annotation_legend_side="right")
dev.off()

##### Activity, Cohort3
celltype_reordered_4=c("CD4","NK","B","Mono")
colnames4 = c("CD4","NK","B","Mono")
row_labels4 = structure(colnames4,names=celltype_reordered_4)
celltype_corresp4 = dis_rep_corresp %>% select(I_subset,C_subset,lineage)%>%filter(!is.na(C_subset))
celltype_corresp4$lineage=factor(celltype_corresp4$lineage,levels=c("CD4","NK","B","Monocyte"))

act2 = fread_FT("res/211224_CommuBiol2021_activity_onebinom_res.txt") %>% 
      select(discovery,replication,binom_logP) %>% pivot_wider(names_from="discovery",values_from="binom_logP") %>% column_to_rownames("replication") 

act2 = act2[celltype_reordered_4,c("B05_NCD4","G01_NK","A01_NaiB","E02_CD16nMo")]
column_labels4=column_labels[c("B05_NCD4","G01_NK","A01_NaiB","E02_CD16nMo")]
# Corresponding p<10e-3

bonf_threshold = -log10(0.05/(ncol(act2)*nrow(act2)))

col_fun = colorRamp2(c(bonf_threshold,40), c("#FFFFFF","#BC3C29FF"))

ha = HeatmapAnnotation(lineage=celltype_corresp4$lineage,col=list(lineage=col2),show_annotation_name=F,
                    height=unit(0.2,"cm"),simple_anno_size_adjust=T,show_legend=T,
                    annotation_legend_param=list(title_gp=gpar(fontsize=12),labels_gp=gpar(fontsize=12)))

ra = rowAnnotation(lineage=celltype_corresp4$lineage,col=list(lineage=col2),show_annotation_name=F,
                    width=unit(0.2,"cm"),simple_anno_size_adjust=T,show_legend=F,
                    annotation_legend_param=list(title_gp=gpar(fontsize=12),labels_gp=gpar(fontsize=12)))

p=Heatmap(act2,cluster_rows=FALSE,cluster_columns=FALSE,,col=col_fun,width=unit(1.7,"cm"),height=unit(1.7,"cm"),
                       row_names_gp=gpar(fontsize=12),row_labels=row_labels4,
                       column_names_gp=gpar(fontsize=12),column_labels=column_labels4,
                       right_annotation = ra,
                       bottom_annotation = ha,
                       heatmap_legend_param=list(title="-log10(P)",at=c(2.5,20,40),
                                            direction="horizontal",title_position="topcenter",
                                            grid_width=unit(2,"cm"),legend_height=unit(0.4,"cm"),
                                            title_gp=gpar(fontsize=12,fontface="bold"),labels_gp=gpar(fontsize=12)))

pdf_3(paste0("res/",today(),"_CommuBiol2021_activity_onebinom_res.pdf"),h=3,w=3.5)
 draw(p,heatmap_legend_side="bottom", annotation_legend_side="right")
dev.off()

######### Supple Table
# Cell type corresp table
dis_rep_corresp2=dis_rep_corresp%>%select(I_label,G_subset,A_subset,C_subset)
dis_rep_corresp2$A_subset=ifelse(is.na(dis_rep_corresp2$A_subset)==T,NA,
                          ifelse(dis_rep_corresp2$A_subset=="C02_MCD8","Mem CD8",dis_rep_corresp2$I_label))
dis_rep_corresp2$G_subset=ifelse(dis_rep_corresp2$G_subset=="cMo","CL Mono",
                          ifelse(dis_rep_corresp2$G_subset=="cDC","mDC",
                          ifelse(dis_rep_corresp2$G_subset=="PMN","Neu",dis_rep_corresp2$G_subset)))
colnames(dis_rep_corresp2)=c("ImmuNexUT","Cohort 1","Cohort 2","Cohort 3")

write.table_FT_2(dis_rep_corresp2,paste0(today(),"_rep_celltype_corresp_forSupple.txt"))


