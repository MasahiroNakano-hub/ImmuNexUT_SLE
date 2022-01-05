#---------------------------------------#
# 211204 Nakano
# Belimumab
# lme4-GLMM
#---------------------------------------#

source("my.source.R")
options(stringsAsFactors=F)

# BLM  DEGs
############################################################################################
BLM_samples = fread_FT("data_ref/211031_BEL_RNAseq_sample.txt")[,c(1:2,4,8)]
BLM_samples_2 = BLM_samples %>% pivot_longer(cols=-c(indiv,response),names_to="time",values_to="id") %>% select(id,time,indiv,response)
write.table_FT_2(BLM_samples_2, paste0("data/",today(),"_COI_BEL_pair_list.txt"))

# count cond SLEinc+HC list
list_1 = make_list("~/ws/2021/211124_COI_SLE_HC_run27_finaldata_Combat_PCA/data/SLEinc_andHC/each_subset/each_subset_count","_count.txt",tag="count")
list_1$subset = take_factor(list_1$FILE_count,5:6,"_")
list_2 = make_list("~/ws/2021/211124_COI_SLE_HC_run27_finaldata_Combat_PCA/data/SLEinc_andHC/each_subset/each_subset_cond","_cond.txt",tag="cond")
list_2$subset = take_factor(list_2$FILE_cond,5:6,"_")
list_3 = make_list("~/ws/2021/211129_COI_SLE_HC_run27_edgeR/res/stateactivity/HDAvsInactive/GLMresult","_GLMresult.txt",tag="uniqueres")
list_3$subset = take_factor(list_3$FILE_uniqueres,3:4,"_")
list = full_join(list_1,list_2,by="subset") %>% full_join(.,list_3,by="subset")


# Devide by 500genes
for(iii in 1:nrow(list_3)){
   subset_tmp=list_3$subset[iii]
   data_tmp=fread_FT(list_3$PATH_uniqueres[iii])
   genenum_tmp=nrow(data_tmp)
   paste0(subset_tmp,"; ",genenum_tmp) %>% print()

   chunknum_tmp=floor(genenum_tmp/500)+1

   list_tmp=list%>%filter(subset==subset_tmp)
   list_tmp_2=data.frame(matrix(rep(NA,ncol(list_tmp)*chunknum_tmp), nrow=chunknum_tmp))
   colnames(list_tmp_2)=colnames(list_tmp)

   for(ppp in 1:ncol(list_tmp)){
    list_tmp_2[,ppp]=list_tmp[,ppp]   
   }
   list_tmp_2$chunk=(1:chunknum_tmp)

   if(iii==1){list_sum=list_tmp_2}else{list_sum=rbind(list_sum,list_tmp_2)}
}

[1] "A01_NaiB; 12462"
[1] "A02_SwiMB; 12723"
[1] "A03_UnswMB; 12701"
[1] "A04_DNB; 12796"
[1] "A05_PB; 11447"
[1] "B01_Th1; 12729"
[1] "B02_Th2; 12697"
[1] "B03_TH17; 12640"
[1] "B04_Tfh; 12670"
[1] "B05_NCD4; 12626"
[1] "B06_MCD4; 12791"
[1] "B07_aTreg; 12385"
[1] "B09_Fra1; 12890"
[1] "B10_Fra3; 12749"
[1] "C01_NCD8; 12685"
[1] "C03_EffectorCD8; 12457"
[1] "C04_CmCD8; 12844"
[1] "C05_EmCD8; 12643"
[1] "D01_mDC; 12447"
[1] "D02_pDC; 12351"
[1] "E01_CD16pMo; 11810"
[1] "E02_CD16nMo; 11976"
[1] "E03_NonClassical; 11758"
[1] "E04_Intermediate; 11827"
[1] "F01_Neu; 10109"
[1] "F02_LDG; 10956"
[1] "G01_NK; 12590"

write.table_FT_2(list_sum, paste0("tmp_job_list/",today(),"_COI_SLEinc_andHC_count_chunkby500_list.txt"))


######################################################################
######################################################################
# GLMM
# Example: All pts, One chunk
######################################################################
######################################################################

LIST="tmp_job_list/211204_COI_SLEinc_andHC_count_chunkby500_list.txt"
##  task_id=1

list       = fread_FT(LIST)
subset_tmp = list[task_id,3]
count_path = list[task_id,1]
batch_path = list[task_id,4]
unique_path= list[task_id,6]
chunk_tmp  = list[task_id,8]

paste0("subset_tmp : ",subset_tmp) %>% print()
paste0("count_path  : ",  count_path )       %>% print()
paste0("batch_path : ",batch_path) %>% print()
paste0("unique_path : ",unique_path) %>% print()
paste0("chunk_tmp : ",chunk_tmp) %>% print()

###########################################################################

pts_list=fread_FT("data/211204_COI_BEL_pair_list.txt")
###################
################################

count_tmp  = fread_n(count_path)
cond_tmp   = fread_FT(batch_path)
unique_tmp   = fread_FT(unique_path)

cond_tmp = inner_join(pts_list,cond_tmp,by="id")
count_tmp = count_tmp[,cond_tmp$name]

## Check
if(ncol(count_tmp)==nrow(cond_tmp)){"sample and cond were same num"}else{"please check condition list"}

cond_tmp_2  = cond_tmp[is.element(cond_tmp$name,colnames(count_tmp)),]
count_tmp_2 = count_tmp[is.element(rownames(count_tmp),unique_tmp$genes),cond_tmp_2$name] # in line with discovery cohort


print("## final use data number ##")
table_freq(cond_tmp_2$time)
print("## final use data number ##")

#REFERNCE = unique(cond_tmp_2$batch) %>% grep(TARGET,.,invert=T,value=T)
cond_tmp_2$time  = factor(cond_tmp_2$time, levels=sort(unique(cond_tmp_2$time)))
count_tmp_3= count_tmp_2

## Check
print("## Low Exp Gene filtering ##")
paste0("original gene num : ",nrow(count_tmp))     %>% print()
paste0("after filter gene num : ",nrow(count_tmp_3)) %>% print()
print("## Low Exp Gene filtering ##")

## edgeR norm
dge     = DGEList(counts=count_tmp_3, genes=rownames(count_tmp_3))
dge_n   = calcNormFactors( dge ,method="TMM")

design_1     = model.matrix(~ 0+time+indiv+lotseq , data= cond_tmp_2 )
################################################################################

# GLMM
offset_tmp = dge_n$samples %>% mutate(norm_libsize=lib.size*norm.factors) %>%
             rownames_to_column("name") %>% .[,c(1,5)]

chunk_num=floor(nrow(count_tmp_3)/500)+1

if(chunk_tmp==chunk_num){data_tmp=dge_n$counts[(1+(chunk_tmp-1)*500):nrow(count_tmp_3),]}else{
                         data_tmp=dge_n$counts[(1+(chunk_tmp-1)*500):(chunk_tmp*500),]}

data_tmp = data_tmp %>% t() %>% as.data.frame() %>% 
           rownames_to_column("name") %>% left_join(.,cond_tmp_2 %>% dplyr::select(name,time,indiv,response,lotseq,Activity),by="name") %>%
           left_join(.,offset_tmp,by="name") %>% column_to_rownames("name")

res_sum=data.frame(matrix(rep(NA,(ncol(data_tmp)-6)*7),ncol=7))
colnames(res_sum)=c("subset","gene","full_beta","full_P","full_theta","null_theta","anova_P")

for(iii in 1:nrow(res_sum)){
   data_gene_tmp = data_tmp[,c(iii,(nrow(res_sum)+1):ncol(data_tmp))]
   gene_tmp = colnames(data_gene_tmp)[1]
   colnames(data_gene_tmp)[1] ="exp"
   print(paste0(subset_tmp,"; ",gene_tmp," start"))
   a = try(glmer.nb(exp~time+(1|indiv)+lotseq,offset=log(norm_libsize),data_gene_tmp))
   b = try(glmer.nb(exp~     (1|indiv)+lotseq,offset=log(norm_libsize),data_gene_tmp))
   if(class(a)=="try-error"|class(b)=="try-error"){
      res_sum[iii,] = c(subset_tmp,gene_tmp,NA,NA,NA,NA,NA)
   }else{
   theta_a=getME(a, "glmer.nb.theta")
   theta_b=getME(b, "glmer.nb.theta")
   beta_a=summary(a)$coef[2,1]
   pval_a=summary(a)$coef[2,4]
   pval_anova=anova(a,b)[2,8]
      res_sum[iii,] = c(subset_tmp,gene_tmp,beta_a,pval_a,theta_a,theta_b,pval_anova)
   }
}

write.table_FT_2(res_sum,paste0("~/ws/2021/211204_COI_SLE_HC_run27_BLM/GLMM/All/",today(),"_COI_BLM_",subset_tmp,"_All_Chunk",chunk_tmp,"_GLMMres.txt"))

######################################################################
######################################################################

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

##########################
# Summarize results
##########################
# All
list1        = make_list("GLMM/All","_GLMMres.txt")
list1$subset = take_factor(list1$FILE,4:5,"_")
list1$chunk  = take_factor(list1$FILE,7,"_") %>% gsub("Chunk","",.)
dim(list1)
# [1] 677   4   OK!!

for(iii in 1:length(unique(list1$subset))){
    subset_tmp=unique(list1$subset)[iii]
    list_tmp=list1%>%filter(subset==subset_tmp)
    list_tmp$chunk=as.numeric(list_tmp$chunk)
    list_tmp=list_tmp[order(list_tmp$chunk),]

  for(kkk in 1:nrow(list_tmp)){
    data_tmp=fread_FT(list_tmp$PATH[kkk])
    if(kkk==1){data_sum=data_tmp}else{data_sum=rbind(data_sum,data_tmp)}
  }
  data_sum$full_FDR=p.adjust(data_sum$full_P,method="BH")
  data_sum$anova_FDR=p.adjust(data_sum$anova_P,method="BH")
  #paste0(subset_tmp,"; ",nrow(data_sum)) %>% print()
  write.table_FT_2(data_sum, paste0("GLMM/All_sum/",today(),"_COI_BLM_",subset_tmp,"_All_GLMMres.txt"))
}


### DEGnum barplot
list        = make_list("GLMM/All_sum","_GLMMres.txt")
list$subset = take_factor(list$FILE,4:5,"_")

for(iii in 1:nrow(list)){
  subset_tmp   = list$subset[iii]
  df_1         = fread_FT(list$PATH[iii]) %>% filter(anova_FDR<0.05)
  df_tmp = data.frame(DEGtype="All",subset=subset_tmp,DEGnum=nrow(df_1))
  if(iii==1){res_sum=df_tmp}else{res_sum=rbind(res_sum,df_tmp)}
}
res_sum = left_join(res_sum,celltype_corresp%>%select(subset,lineage),by="subset") 
write.table_FT_2(res_sum,paste0("GLMM/",today(),"_COI_BLM_All_DEGnum.txt"))

res_sum$lineage = factor(res_sum$lineage,levels=c("CD4","CD8","NK","B","Monocyte","DC","Neutrophil"))
res_sum$DEGnum  = as.numeric(res_sum$DEGnum)
res_sum$subset  = factor(res_sum$subset,levels=celltype_reordered_27)

p = ggplot()+
     geom_bar(data=res_sum, aes(x=subset,y=DEGnum,fill=lineage),stat="identity")+
     theme_classic()+
     #facet_wrap(~ response, ncol = 1) + 
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

 pdf_3(paste0("GLMM/",today(),"_COI_BLM_All_DEGnum.pdf"),h=3.5,w=9.2)
  plot(p)
 dev.off()
############

# Good, Poor
list1        = make_list("GLMM/Good","_GLMMres.txt")
list1$subset = take_factor(list1$FILE,4:5,"_")
list1$chunk  = take_factor(list1$FILE,7,"_") %>% gsub("Chunk","",.)
list1$type   ="Good"
dim(list1)
# [1] 677   5   OK!!

list2        = make_list("GLMM/Poor","_GLMMres.txt")
list2$subset = take_factor(list2$FILE,4:5,"_")
list2$chunk  = take_factor(list2$FILE,7,"_") %>% gsub("Chunk","",.)
list2$type   ="Poor"
dim(list2)
# [1] 677   5   OK!!

list3=rbind(list1,list2)
list3$subset_type=paste0(list3$subset,"_",list3$type)

for(iii in 1:length(unique(list3$subset_type))){
    subset_type_tmp=unique(list3$subset_type)[iii]
    list_tmp=list3%>%filter(subset_type==subset_type_tmp)
    subset_tmp=unique(list_tmp$subset)
    type_tmp=unique(list_tmp$type)
    
    list_tmp$chunk=as.numeric(list_tmp$chunk)
    list_tmp=list_tmp[order(list_tmp$chunk),]

  for(kkk in 1:nrow(list_tmp)){
    data_tmp=fread_FT(list_tmp$PATH[kkk])
    if(kkk==1){data_sum=data_tmp}else{data_sum=rbind(data_sum,data_tmp)}
  }
  data_sum$full_FDR=p.adjust(data_sum$full_P,method="BH")
  data_sum$anova_FDR=p.adjust(data_sum$anova_P,method="BH")
  #paste0(subset_type_tmp,"; ",nrow(data_sum)) %>% print()
  write.table_FT_2(data_sum, paste0("GLMM/",type_tmp,"_sum/",today(),"_COI_BLM_",subset_type_tmp,"_GLMMres.txt"))
}

### DEGnum barplot
list1        = make_list("GLMM/Good_sum","_GLMMres.txt")
list1$subset = take_factor(list1$FILE,4:5,"_")
list1$type   = "Good"
list2        = make_list("GLMM/Poor_sum","_GLMMres.txt")
list2$subset = take_factor(list2$FILE,4:5,"_")
list2$type   = "Poor"
list=rbind(list1,list2)

for(iii in 1:nrow(list)){
  subset_tmp   = list$subset[iii]
  type_tmp     = list$type[iii]
  df_1         = fread_FT(list$PATH[iii]) %>% filter(anova_FDR<0.05)
  df_tmp = data.frame(DEGtype=type_tmp,subset=subset_tmp,DEGnum=nrow(df_1))
  if(iii==1){res_sum=df_tmp}else{res_sum=rbind(res_sum,df_tmp)}
}
res_sum = left_join(res_sum,celltype_corresp%>%select(subset,lineage),by="subset") 
write.table_FT_2(res_sum,paste0("GLMM/",today(),"_COI_BLM_GoodPoor_DEGnum.txt"))

###################################################################################################
# Disease-activity and BLM-DEGs logFC
###################################################################################################

list5 = make_list("~/ws/2021/211129_COI_SLE_HC_run27_edgeR/res/stateactivity/HDAvsInactive/GLMresult","_HDAvsInactive_GLMresult.txt",tag="activity")
list5$subset = take_factor(list5$FILE_activity,3:4,"_")

# Good, Poor logFC plot
for(kkk in 1:length(unique(list$type))){
   type_tmp  = unique(list$type)[kkk]
   list_tmp  = list %>% filter(type==type_tmp)
　　list_tmp_2 = left_join(list_tmp,list5,by="subset") 

   for(iii in 1:nrow(list_tmp_2)){
       subset_tmp         = list_tmp_2$subset[iii]
       df_tmp_b           = fread_FT(list_tmp_2$PATH[iii]) %>% select(gene,full_beta,anova_P,anova_FDR) 
       colnames(df_tmp_b) = c("genes","belimumab_logFC","belimumab_P","belimumab_FDR")

       df_tmp_a           = read.table_FT(list_tmp_2$PATH_activity[iii]) %>% select(genes,logCPM,logFC,PValue,FDR)
       colnames(df_tmp_a)[3:5] = c("activity_logFC","activity_P","activity_FDR")

       all.equal(df_tmp_b$genes%>%sort(),df_tmp_a$genes%>%sort()) %>% paste0(.,";",type_tmp,";",subset_tmp) %>% print()
       
       df_tmp_b = df_tmp_b %>% .[!apply(.,1,anyNA),]
       df_tmp = inner_join(df_tmp_a,df_tmp_b,by="genes")
       df_tmp$subset =  subset_tmp

       df_tmp_lim=df_tmp%>%filter(activity_FDR<0.05)

       lm_lim       = lm(belimumab_logFC~activity_logFC,data=df_tmp_lim) 
       spe_lim=cor.test(df_tmp_lim$activity_logFC,df_tmp_lim$belimumab_logFC,method="spearman")
       lm_tmp       = data.frame(subset=subset_tmp,lm_beta=summary(lm_lim)$coef[2,1],lm_P=summary(lm_lim)$coef[2,4],
                                    spe_cor=spe_lim$estimate,spe_P=spe_lim$p.value)

       if(iii==1){res_sum=df_tmp_lim}else{res_sum=rbind(res_sum,df_tmp_lim)}
       if(iii==1){lm_sum=lm_tmp}else{lm_sum=rbind(lm_sum,lm_tmp)}
     }

  write.table_FT_2(res_sum,paste0("GLMM/similarity/",today(),"_COI_",type_tmp,"_activity_BEL_logFC_logCPM_P_FDR_actvititysig.txt"))
  write.table_FT_2(lm_sum,paste0("GLMM/similarity/",today(),"_COI_",type_tmp,"_activity_BEL_lmsperes_activitysig.txt"))

  res_sum$subset   = factor(res_sum$subset,levels=celltype_reordered_27)
  
  # logFC>5,-5 → 5
  res_sum_max5 = res_sum
  res_sum_max5$activity_logFC[res_sum_max5$activity_logFC>5] = 5 
  res_sum_max5$activity_logFC[res_sum_max5$activity_logFC< -5] = -5 

  res_sum_max5$belimumab_logFC[res_sum_max5$belimumab_logFC>5] = 5 
  res_sum_max5$belimumab_logFC[res_sum_max5$belimumab_logFC< -5] = -5 

  p = ggplot(res_sum_max5,aes(x=activity_logFC,y=belimumab_logFC))+
          geom_point(size=0.2,color="#BC3C29FF")+
          theme_classic()+
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
          scale_x_continuous(breaks=c(-4,-2,0,2,4),limits=c(-5,5))+
          scale_y_continuous(breaks=c(-4,-2,0,2,4),limits=c(-5,5))+
          labs(x=paste0("HDA vs Inactive (logFC)"),y=paste0("Post- vs Pre-BLM (logFC)"))

  png(paste0("GLMM/similarity/",today(),"_COI_27subsets_activity_BLM_",type_tmp,"_FDR005_scatter_activitysig.png"),height=7,width=10,units="in",res=300)
   plot(p)
  dev.off()
}

###################################################################################################
# DEGnum and Jaccard barplot
###################################################################################################

list=make_list("GLMM/similarity","_activity_BEL_logFC_logCPM_P_FDR.txt")
list$type=take_factor(list$FILE,3,"_")

for(kkk in 1:length(unique(list$type))){
   type_tmp  = unique(list$type)[kkk]
   list_tmp  = list %>% filter(type==type_tmp)
   res_sum=fread_FT(list_tmp$PATH)

   for(iii in 1:length(celltype_reordered_27)){
    subset_tmp = celltype_reordered_27[iii]
    res_tmp = res_sum %>% filter(subset==subset_tmp)
    BLM_tmp = res_tmp %>% filter(siggroup005=="only_belimumab")
    activity_tmp = res_tmp %>% filter(siggroup005=="only_activity")
    both_tmp  = res_tmp %>% filter(siggroup005=="significant_in_both")
    dis_tmp   = res_tmp %>% filter(siggroup005=="discordant")
    jaccard_tmp = data.frame(type=type_tmp,subset=subset_tmp,BLM=nrow(BLM_tmp),activity=nrow(activity_tmp),both=nrow(both_tmp),
                             discordant=nrow(dis_tmp),union=nrow(BLM_tmp)+nrow(activity_tmp)+nrow(both_tmp),
                             all=nrow(BLM_tmp)+nrow(activity_tmp)+nrow(both_tmp)+nrow(dis_tmp),
                             BLM_P=nrow(BLM_tmp)/(nrow(BLM_tmp)+nrow(activity_tmp)+nrow(both_tmp)),
                             activity_P=nrow(activity_tmp)/(nrow(BLM_tmp)+nrow(activity_tmp)+nrow(both_tmp)),
                             jaccard=nrow(both_tmp)/(nrow(BLM_tmp)+nrow(activity_tmp)+nrow(both_tmp)))
    if(iii==1){jaccard_sum=jaccard_tmp}else{jaccard_sum=rbind(jaccard_sum,jaccard_tmp)}
    }
    if(kkk==1){jaccard_sum2=jaccard_sum}else{jaccard_sum2=rbind(jaccard_sum2,jaccard_sum)}
}
write.table_FT_2(jaccard_sum2,paste0("GLMM/similarity/",today(),"_COI_BLM_GoodPoor_jaccard.txt"))

DEG_sum=fread_FT("GLMM/211207_COI_BLM_GoodPoor_DEGnum.txt") %>% mutate(type_subset=paste0(DEGtype,"_",subset))

jaccard_sum2=jaccard_sum2 %>% mutate(type_subset=paste0(type,"_",subset)) %>% select(-c(type,subset))

all_sum=left_join(DEG_sum,jaccard_sum2,by="type_subset")
all_sum$jaccard2=all_sum$jaccard*30000
all_sum$jaccard2=ifelse(all_sum$BLM>300,all_sum$jaccard2,0) # Jaccard only for DEGs>300
all_sum=all_sum %>% select(DEGtype,subset,DEGnum,jaccard2) %>% 
                    pivot_longer(col=-c(DEGtype,subset),names_to="parameter",values_to="value")
colnames(all_sum)[1]="Response"
all_sum$subset  = factor(all_sum$subset,levels=celltype_reordered_27)
all_sum$Response  = factor(all_sum$Response,levels=c("Good","Poor"))
levels(all_sum$Response)=list(`BLM (Good)`="Good",`BLM (Poor)`="Poor")
all_sum$parameter  = factor(all_sum$parameter,levels=c("DEGnum","jaccard2"))
levels(all_sum$parameter)=list(`Number of DEGs`="DEGnum",`Jaccard similarity index`="jaccard2")

p = ggplot()+
     geom_bar(data=all_sum, aes(x=subset,y=value,fill=parameter),stat="identity",position="dodge")+
     theme_classic()+
     facet_wrap(~ Response, ncol=1) + 
     scale_fill_manual(values=c("#20854EFF","#E18727FF"))+
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
           legend.title=element_blank(),
           legend.text=element_text(colour="black",size=16))+
     scale_x_discrete(labels= label)+
     scale_y_continuous(name="Number of DEGs",
                        sec.axis=sec_axis(trans=~.*(1/30000), name="",breaks=c(0,0.025,0.05)))

 pdf_3(paste0("GLMM/similarity/",today(),"_COI_BLM_GoodPoor_DEGnum_Jaccard.pdf"),h=5,w=11.5)
  plot(p)
 dev.off()

# jaccard_sum2=fread_FT("GLMM/similarity/211207_COI_BLM_GoodPoor_jaccard.txt")
jaccard_sum_lim=jaccard_sum2%>%filter(subset%in%c("A01_NaiB","A02_SwiMB","A03_UnswMB")) %>%
                               select(type,subset,jaccard)%>%
                               pivot_wider(names_from=type,values_from=jaccard)
jaccard_sum_lim$ratio=jaccard_sum_lim$Good/jaccard_sum_lim$Poor  
mean(jaccard_sum_lim$ratio)                        

###################################################################################################
# lm beta
###################################################################################################

list=make_list("GLMM/similarity","_activity_BEL_lmres.txt")
list$type=take_factor(list$FILE,3,"_")

for(kkk in 1:length(unique(list$type))){
   type_tmp  = unique(list$type)[kkk]
   list_tmp  = list %>% filter(type==type_tmp)
   res_sum=fread_FT(list_tmp$PATH)%>%left_join(.,celltype_corresp%>%select(subset,lineage),by="subset")
   res_sum$Response=type_tmp
   if(kkk==1){res_sum2=res_sum}else{res_sum2=rbind(res_sum2,res_sum)}
}
res_sum2$FDR=p.adjust(res_sum2$P,method="BH")
res_sum2$BonfP=p.adjust(res_sum2$P,method="bonferroni")
write.table_FT_2(res_sum2,paste0("GLMM/similarity/",today(),"_COI_BLM_GoodPoor_activity_lmres.txt"))

res_sum2$lineage = factor(res_sum2$lineage,levels=c("CD4","CD8","NK","B","Monocyte","DC","Neutrophil"))
res_sum2$subset  = factor(res_sum2$subset,levels=celltype_reordered_27)
res_sum2$Response  = factor(res_sum2$Response,levels=c("Good","Poor"))
levels(res_sum2$Response)=list(`BLM (Good)`="Good",`BLM (Poor)`="Poor")

p = ggplot()+
  geom_bar(data=res_sum2, aes(x=subset,y=beta,fill=lineage),stat="identity")+
  theme_classic()+
  facet_wrap(~ Response, ncol = 1) + 
  scale_fill_manual(values = col2)+
  theme(axis.text.x=element_text(colour="black",angle=45,hjust=1,vjust=1,size=16),
        axis.text.y=element_text(colour="black",size=16),
        axis.title.x=element_blank(),
        axis.title.y=element_text(colour="black",size=16),
        plot.title=element_blank(),
        strip.background=element_blank(),
        strip.text=element_text(colour="black",size=18),
        panel.background =element_rect(fill="transparent",colour="black",size=1),
        panel.grid = element_blank(),
        legend.position="right",
        legend.title=element_text(colour="black",size=16),
        legend.text=element_text(colour="black",size=16))+
  scale_x_discrete(labels= label)+
  scale_y_continuous(breaks=c(-0.4,-0.2,0,0.2,0.4),limits=c(-0.4,0.4))+
  labs(y="Effect size")

    pdf_3(paste0("GLMM/similarity/",today(),"_COI_BLM_GoodPoor_activity_lmres.pdf"),h=6,w=9.2)
     plot(p)
    dev.off()

# Supple

list1        = make_list("GLMM/Good_sum","_GLMMres.txt")
list1$subset = take_factor(list1$FILE,4:5,"_")
list1$type   = "Good"
list2        = make_list("GLMM/Poor_sum","_GLMMres.txt")
list2$subset = take_factor(list2$FILE,4:5,"_")
list2$type   = "Poor"
list=rbind(list1,list2)

for(iii in 1:nrow(list)){
  subset_tmp   = list$subset[iii]
  type_tmp     = list$type[iii]
  res_tmp      = fread_FT(list$PATH[iii]) %>%mutate(type=type_tmp)
  if(iii==1){res_sum=res_tmp}else{res_sum=rbind(res_sum,res_tmp)}
}

res_sum = res_sum %>% left_join(.,celltype_corresp[,c(1,2)],by="subset") 
res_sum$subset=factor(res_sum$subset,levels=celltype_reordered_27)
res_sum$type=factor(res_sum$type,levels=c("Good","Poor"))
res_sum = res_sum[order(res_sum$gene),]
res_sum = res_sum[order(res_sum$subset),]
res_sum = res_sum[order(res_sum$type),]

res_sum = res_sum %>% select(type,label,gene,full_beta,anova_P,anova_FDR) 
res_sum$full_beta  = formatC(res_sum$full_beta,digits=2)
res_sum$anova_P = formatC(res_sum$anova_P,digits=2)
res_sum$anova_FDR    = formatC(res_sum$anova_FDR,digits=2)
colnames(res_sum) =c("Reponse","Cell type","Gene","logFC","Pvalue","FDR")

write.table_FT_2(res_sum,paste0(today(),"_COI_BLM_logFC_P_FDR_Supple.txt"))

###################################################################################################
# Pathway TF
###################################################################################################
# DEG list 4celltype
list1         = make_list("GLMM/All_sum","_GLMMres.txt")[1:4,]
list1$subset  = take_factor(list1$FILE,4:5,"_")
list1$type   = "All"
list2         = make_list("GLMM/Good_sum","_GLMMres.txt")[1:4,]
list2$subset  = take_factor(list2$FILE,4:5,"_")
list2$type   = "Good"
list3         = make_list("GLMM/Poor_sum","_GLMMres.txt")[1:4,]
list3$subset  = take_factor(list3$FILE,4:5,"_")
list3$type   = "Poor"

list = rbind(list1,list2) %>% rbind(.,list3)
list$type_subset = paste0(list$type,"_",list$subset)

list$back_PATH = "/home/mnakano/ws/2021/211129_COI_SLE_HC_run27_edgeR/res/211202_COI_27subsets_stateactivity_DEGunionlist.txt"
write.table_FT_2(list,paste0("tmp_job_list/",today(),"_COI_BEL_GLMMres_list_backDEGunion.txt"))

###########################################################
# Example: one signature
###########################################################

LIST = "tmp_job_list/211214_COI_BEL_GLMMres_list_backDEGunion.txt"
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
target_gene_list   =  read.table_FT(TARGET_PATH) %>% filter(anova_FDR<0.05) %>%.[,2]
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

### HALLMARK, KEGG representative pathways
list_tmp = make_list("ORA/HALLMARK",".txt")
list_tmp$subset = take_factor(list_tmp$FILE,2:3,"_")
list_tmp$type   = take_factor(list_tmp$FILE,4,"_")
list_tmp$anno   = take_factor(list_tmp$FILE,6,"_") %>% gsub(".txt","",.) 
list_tmp=list_tmp %>% filter(type%in%c("Good","Poor"))

focus=c("HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_TNFA_SIGNALING_VIA_NFKB")

 for(kkk in 1:nrow(list_tmp)){
   subset_tmp=list_tmp$subset[kkk]
   type_tmp  =list_tmp$type[kkk]
   res_tmp   = fread_FT(list_tmp$PATH[kkk]) %>% filter(Description%in%focus) %>% select(Description,pvalue)
   if(nrow(res_tmp)==0){res_tmp=data.frame(Description=focus,pvalue=1)}else{res_tmp=res_tmp}
   res_tmp$subset = subset_tmp
   res_tmp$type   = type_tmp
   if(kkk==1){res_sum1=res_tmp}else{res_sum1=union(res_sum1,res_tmp)}
   } 

list_tmp = make_list("ORA/KEGG",".txt")
list_tmp$subset = take_factor(list_tmp$FILE,2:3,"_")
list_tmp$type   = take_factor(list_tmp$FILE,4,"_")
list_tmp$anno   = take_factor(list_tmp$FILE,6,"_") %>% gsub(".txt","",.) 
list_tmp=list_tmp %>% filter(type%in%c("Good","Poor"))

focus=c("Glycolysis / Gluconeogenesis")

 for(kkk in 1:nrow(list_tmp)){
   subset_tmp=list_tmp$subset[kkk]
   type_tmp  =list_tmp$type[kkk]
   res_tmp   = fread_FT(list_tmp$PATH[kkk]) %>% filter(Description%in%focus) %>% select(Description,pvalue)
   if(nrow(res_tmp)==0){res_tmp=data.frame(Description=focus,pvalue=1)}else{res_tmp=res_tmp}
   res_tmp$subset = subset_tmp
   res_tmp$type   = type_tmp
   if(kkk==1){res_sum2=res_tmp}else{res_sum2=union(res_sum2,res_tmp)}
   } 

res_sum=rbind(res_sum1,res_sum2)
write.table_FT_2(res_sum,paste0("ORA/",today(),"_BLM_ORA_HALLMARK_KEGG_rep_logP.txt"))


res_sum$subset=factor(res_sum$subset,levels=rev(celltype_reordered_27))
res_sum$type  =factor(res_sum$type,levels=c("Good","Poor"))
levels(res_sum$type) =c("BLM (Good)","BLM (Poor)")
res_sum$Description  =factor(res_sum$Description,levels=c("HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_TNFA_SIGNALING_VIA_NFKB","Glycolysis / Gluconeogenesis"))
levels(res_sum$Description) =c("HALLMARK_IFN_GAMMA","HALLMARK_TNFA_NFKB","KEGG_Glycolysis_Gluconeogenesis")

pal_npg("nrc")(5)
# [1] "#E64B35FF" "#4DBBD5FF" "#00A087FF" "#3C5488FF" "#F39B7FFF"

IFNG=res_sum_tmp=res_sum%>%filter(Description=="HALLMARK_IFN_GAMMA")

p = ggplot()+
     geom_bar(data=IFNG, aes(x=subset,y=-log10(pvalue)),fill="#E64B35FF",stat="identity")+
     theme_classic()+
     coord_flip()+
     facet_wrap(~type,ncol=1) + 
     theme(axis.text.x=element_text(colour="black",size=16),
           axis.text.y=element_text(colour="black",size=16),
           axis.title.x=element_text(colour="black",size=16),
           axis.title.y=element_blank(),
           strip.background=element_blank(),
           strip.text=element_text(colour="black",size=16),
           panel.background =element_rect(fill="transparent",colour="black",size=1),
           panel.grid = element_blank(),
           plot.title=element_blank(),
      　　　legend.position="none")+
     scale_x_discrete(labels= label)+
　　　labs(y="-log10(P)")

    pdf_3(paste0("ORA/",today(),"_ORA_HALLMARK_IFN_GAMMA_rep_logP.pdf"),h=3,w=3)
     plot(p)
    dev.off()

NFKB=res_sum_tmp=res_sum%>%filter(Description=="HALLMARK_TNFA_NFKB")

p = ggplot()+
     geom_bar(data=NFKB, aes(x=subset,y=-log10(pvalue)),fill="#4DBBD5FF",stat="identity")+
     theme_classic()+
     coord_flip()+
     facet_wrap(~type,ncol=1) + 
     theme(axis.text.x=element_text(colour="black",size=16),
           axis.text.y=element_blank(),
           axis.title.x=element_text(colour="black",size=16),
           axis.title.y=element_blank(),
           axis.ticks.y=element_blank(),
           strip.background=element_blank(),
           strip.text=element_text(colour="black",size=16),
           panel.background =element_rect(fill="transparent",colour="black",size=1),
           panel.grid = element_blank(),
           plot.title=element_blank(),
      　　　legend.position="none")+
     scale_x_discrete(labels= label)+
　　　labs(y="-log10(P)")

    pdf_3(paste0("ORA/",today(),"_ORA_HALLMARK_TNFA_NFKB_rep_logP.pdf"),h=3,w=2.2)
     plot(p)
    dev.off()

Glyco=res_sum_tmp=res_sum%>%filter(Description=="KEGG_Glycolysis_Gluconeogenesis")

p = ggplot()+
     geom_bar(data=Glyco, aes(x=subset,y=-log10(pvalue)),fill="#3C5488FF",stat="identity")+
     theme_classic()+
     coord_flip()+
     facet_wrap(~type,ncol=1) + 
     theme(axis.text.x=element_text(colour="black",size=16),
           axis.text.y=element_blank(),
           axis.title.x=element_text(colour="black",size=16),
           axis.title.y=element_blank(),
           axis.ticks.y=element_blank(),
           strip.background=element_blank(),
           strip.text=element_text(colour="black",size=16),
           panel.background =element_rect(fill="transparent",colour="black",size=1),
           panel.grid = element_blank(),
           plot.title=element_blank(),
      　　　legend.position="none")+
     scale_x_discrete(labels= label)+
　　　labs(y="-log10(P)")

    pdf_3(paste0("ORA/",today(),"_ORA_KEGG_Glycolysis_Gluconeogenesis_rep_logP.pdf"),h=3,w=2.2)
     plot(p)
    dev.off()

####### Supple
list1=make_list("ORA/HALLMARK/Good","_ORA_HALLMARK.txt")
list1$subset=take_factor(list1$FILE,2:3,"_")
for(iii in 1:nrow(list1)){
    subset_tmp=list1$subset[iii]
    res_tmp=fread_FT(list1$PATH[iii])
    res_tmp$subset=subset_tmp
    if(iii==1){res_sum1=res_tmp}else{res_sum1=rbind(res_sum1,res_tmp)}
}
res_sum1$DEGtype="Good"
res_sum1$Database="HALLMARK"


list2=make_list("ORA/HALLMARK/Poor","_ORA_HALLMARK.txt")
list2$subset=take_factor(list2$FILE,2:3,"_")
for(iii in 1:nrow(list2)){
    subset_tmp=list2$subset[iii]
    res_tmp=fread_FT(list2$PATH[iii])
    res_tmp$subset=subset_tmp
    if(iii==1){res_sum2=res_tmp}else{res_sum2=rbind(res_sum2,res_tmp)}
}
res_sum2$DEGtype="Poor"
res_sum2$Database="HALLMARK"

list3=make_list("ORA/KEGG/Good","_ORA_KEGG.txt")
list3$subset=take_factor(list3$FILE,2:3,"_")
for(iii in 1:nrow(list3)){
    subset_tmp=list3$subset[iii]
    res_tmp=fread_FT(list3$PATH[iii])
    res_tmp$subset=subset_tmp
    if(iii==1){res_sum3=res_tmp}else{res_sum3=rbind(res_sum3,res_tmp)}
}
res_sum3$DEGtype="Good"
res_sum3$Database="KEGG"

list4=make_list("ORA/KEGG/Poor","_ORA_KEGG.txt")[3,]
list4$subset=take_factor(list4$FILE,2:3,"_")
for(iii in 1:nrow(list4)){
    subset_tmp=list4$subset[iii]
    res_tmp=fread_FT(list4$PATH[iii])
    res_tmp$subset=subset_tmp
    if(iii==1){res_sum4=res_tmp}else{res_sum4=rbind(res_sum4,res_tmp)}
}
res_sum4$DEGtype="Poor"
res_sum4$Database="KEGG"

res_sum_all3=bind_rows(res_sum1,res_sum2,res_sum3,res_sum4)

res_sum_all3$subset=factor(res_sum_all3$subset,levels=celltype_reordered_27)
res_sum_all3$DEGtype=factor(res_sum_all3$DEGtype,levels=c("Good","Poor"))
res_sum_all3$Database=factor(res_sum_all3$Database,levels=c("HALLMARK","KEGG"))
res_sum_all3=res_sum_all3[order(res_sum_all3$subset),]
res_sum_all3=res_sum_all3[order(res_sum_all3$Database),]
res_sum_all3=res_sum_all3[order(res_sum_all3$DEGtype),]

res_sum_lim3=res_sum_all3%>%left_join(.,celltype_corresp[,c(1,2)],by="subset")%>%
            select(DEGtype,Database,label,Description,pvalue,p.adjust)%>%filter(pvalue<0.05)
res_sum_lim3$pvalue = formatC(res_sum_lim3$pvalue,digits=2)
res_sum_lim3$p.adjust  = formatC(res_sum_lim3$p.adjust,digits=2)
colnames(res_sum_lim3)=c("Response","Database","Cell type","Annotation","Pvalue","FDR")
write.table_FT_2(res_sum_lim3,paste0(today(),"_ORA_HALLMARKKEGG_27subsets_BLM_Supple.txt"))


###################################################################################################
# PC projection Belimumab
###################################################################################################

# STEP 1: we calculate mean SD of all genes from expression matrix
# Each condition
list = make_list("data/SLEunique_andHC/each_subset/LowExp_TMM_logCPM_Combat","_logCPM_Combat.txt",tag="count")
list$subset = take_factor(list$FILE_count,5:6,"_")
list$status = "AfterCombat"

 for(iii in 1:nrow(list)){
  subset_tmp      = list$subset[iii]
  combat_tmp      = fread_n(list$PATH[iii])
  mean = apply(combat_tmp,1,mean)
  var  = apply(combat_tmp,1,var)
  sd   = apply(combat_tmp,1,sd)
  res_tmp         = rbind(mean,var) %>% rbind(.,sd) %>% t() %>% as.data.frame()
  res_tmp_f       = res_tmp[order(res_tmp$var,decreasing=T),]
  write.table_n_2(res_tmp_f,"Gene",paste0("data/SLEunique_andHC/each_subset/LowExp_TMM_logCPM_meanSD/",today(),"_COI_SLEunique_andHC_",subset_tmp,"_logCPM_meanSD.txt"))
 }

# Step 2-1: To calculate expression Z scores of all genes in dupulicated samples, here we use unique patients' mean and SD 
# This is important to project dupulicated samples onto PCA spaces of unique patients 
list1 = make_list("data/SLEunique_andHC/each_subset/LowExp_TMM_logCPM_meanSD","_logCPM_meanSD.txt",tag="meansd")
list1$subset = take_factor(list1$FILE_meansd,5:6,"_")
list2 = make_list("data/SLE_duplicated/each_subset/LowExp_TMM_logCPM_Combat","_logCPM_Combat.txt",tag="combat")
list2$subset = take_factor(list2$FILE_combat,5:6,"_")
list3  = full_join(list1,list2,by="subset")
list3$status = "AfterCombat"

 for(iii in 1:nrow(list3)){
  subset_tmp      = list3$subset[iii]
  meansd_tmp      = fread_FT(list3$PATH_meansd[iii])[1:10000,] # Var Top10000 
  combat_tmp      = fread_n(list3$PATH_combat[iii])[meansd_tmp$Gene,] %>%rownames_to_column("Gene") # Var Top10000 and reorder
  minus_mean      = left_join(meansd_tmp %>% select(Gene,mean),combat_tmp,by="Gene") %>% column_to_rownames("Gene")
  minus_mean      = apply(minus_mean,1,function(x){x-x[1]}) %>% t() %>% as.data.frame() %>% select(-mean) %>% rownames_to_column("Gene")
  split_sd        = left_join(meansd_tmp %>% select(Gene,sd),minus_mean,by="Gene") %>% column_to_rownames("Gene")
  split_sd        = apply(split_sd,1,function(x){x/x[1]}) %>% t() %>% as.data.frame() %>% select(-sd)
  write.table_n_2(split_sd,"Gene",paste0("data/SLE_duplicated/each_subset/LowExp_TMM_logCPM_Zscore10000/",today(),"_COI_SLE_duplicated_",subset_tmp,"_logCPM_Zscore.txt"))
 }

# Step 2-2: For verification of our approach as descrivbed below, we also calculate unique samples' Z score
list4  = full_join(list1,list,by="subset")
list4$status = "AfterCombat"

 for(iii in 1:nrow(list4)){
  subset_tmp      = list4$subset[iii]
  meansd_tmp      = fread_FT(list4$PATH_meansd[iii])[1:10000,] # Var Top10000 
  combat_tmp      = fread_n(list4$PATH_count[iii])[meansd_tmp$Gene,] %>%rownames_to_column("Gene") # Var Top10000 and reorder
  minus_mean      = left_join(meansd_tmp %>% select(Gene,mean),combat_tmp,by="Gene") %>% column_to_rownames("Gene")
  minus_mean      = apply(minus_mean,1,function(x){x-x[1]}) %>% t() %>% as.data.frame() %>% select(-mean) %>% rownames_to_column("Gene")
  split_sd        = left_join(meansd_tmp %>% select(Gene,sd),minus_mean,by="Gene") %>% column_to_rownames("Gene")
  split_sd        = apply(split_sd,1,function(x){x/x[1]}) %>% t() %>% as.data.frame() %>% select(-sd)
  write.table_n_2(split_sd,"Gene",paste0("data/SLEunique_andHC/each_subset/LowExp_TMM_logCPM_Zscore10000/",today(),"_COI_SLEunique_andHC_",subset_tmp,"_logCPM_Zscore.txt"))
 }


# Step 3: We infer PC scores of duplicated patients in unique patients' PC spaces as the inner products of gene expression Z scores and gene PC loadings
# Step 3-1: To verify our caluculation methods, we first infer PC scores of unique patients and compare their original PC scpres
list1 = make_list("data/SLEunique_andHC/each_subset/LowExp_TMM_logCPM_Zscore10000","_Zscore.txt",tag="Zscore")
list1$subset = take_factor(list1$FILE_Zscore,5:6,"_")
list5 = make_list("PCA_res/SLEunique_andHC/AfterCombat/tmp_Rdata","_PCAres_tmp1.Rdata",tag="Rdata")
list5$subset = take_factor(list5$FILE_Rdata,7:8,"_")
list9  = make_list("PCA_res/SLEunique_andHC/AfterCombat/PCdata","_PCdata.txt",tag="PCdata")
list9$subset = take_factor(list9$FILE_PCdata,7:8,"_")
list13 = full_join(list1,list5,by="subset") %>% full_join(.,list9,by="subset")
list13$status = "AfterCombat"

list_tmp = list13

 for(iii in 1:nrow(list_tmp)){
  subset_tmp  = list_tmp$subset[iii]
  Z_tmp       = fread_n(list_tmp$PATH_Zscore[iii]) %>% t() %>% as.matrix()

  load(list_tmp$PATH_Rdata[iii])
  rotation = summary(result)$rotation %>% as.data.frame()

　Z_tmp_lim = Z_tmp[,is.element(colnames(Z_tmp),rownames(rotation))]
  all.equal(colnames(Z_tmp_lim),rownames(rotation)) %>% print()
  # [1] TRUE
  rotation_2 = rotation[colnames(Z_tmp_lim),] %>% as.matrix()
  PC_tmp    =  Z_tmp_lim %*% rotation_2 %>% as.data.frame()
  PC_ori    = fread_FT(list_tmp$PATH_PCdata[iii]) %>% 
              select(-c(2:8)) %>% column_to_rownames("name")
  
  all.equal(PC_tmp,PC_ori) %>% print()

  write.table_n_2(PC_tmp,"name",paste0("PCA_res/SLEunique_andHC/AfterCombat/PCrecal/",today(),"_COI_SLEuniqueandHC_",subset_tmp,"_PCrecal.txt"))
 }

# Calculated PC scores were completely consistent with original PC scores

# Step 3-2: Finally, we infer PC scores of duplicated patients
list1 = make_list("data/SLE_duplicated/each_subset/LowExp_TMM_logCPM_Zscore10000","_Zscore.txt",tag="Zscore")
list1$subset = take_factor(list1$FILE_Zscore,5:6,"_")

list5 = make_list("PCA_res/SLEunique_andHC/AfterCombat/tmp_Rdata","_PCAres_tmp1.Rdata",tag="Rdata")
list5$subset = take_factor(list5$FILE_Rdata,7:8,"_")

list9 = full_join(list1,list5,by="subset")
list9$status = "AfterCombat"

list_tmp = list9

 for(iii in 1:nrow(list_tmp)){
  subset_tmp  = list_tmp$subset[iii]

  Z_tmp       = fread_n(list_tmp$PATH_Zscore[iii]) %>% t() %>% as.matrix()

  load(list_tmp$PATH_Rdata[iii])
  rotation = summary(result)$rotation %>% as.data.frame()

　Z_tmp_lim = Z_tmp[,is.element(colnames(Z_tmp),rownames(rotation))]
  all.equal(colnames(Z_tmp_lim),rownames(rotation)) %>% print()
  # [1] TRUE
  rotation_2 = rotation[colnames(Z_tmp_lim),] %>% as.matrix()

  PC_tmp_2    =  Z_tmp_lim %*% rotation_2 %>% as.data.frame()

  write.table_n_2(PC_tmp_2,"name",paste0("PCA_res/SLE_duplicated/AfterCombat/PCproj/",today(),"_COI_SLEduplicated_",subset_tmp,"_PCproj.txt"))
 }
# PC projection finished

###################################################################################################
# PC score compare
###################################################################################################
suppressPackageStartupMessages(library(lme4))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ggsci))

BLM_samples_2=fread_FT("data/211204_COI_BEL_pair_list.txt")

# Unique
list1 = make_list("~/ws/2021/211124_COI_SLE_HC_run27_finaldata_Combat_PCA/PCA_res/SLEunique_andHC/AfterCombat/PCdata","_PCdata.txt",tag="unique")
list1$subset = take_factor(list1$FILE_unique,7:8,"_")
# Duplicated
list2 = make_list("~/ws/2021/211124_COI_SLE_HC_run27_finaldata_Combat_PCA/PCA_res/SLE_duplicated/AfterCombat/PCproj","_PCproj.txt",tag="dupli")
list2$subset = take_factor(list2$FILE_dupli,5:6,"_")
list = left_join(list1,list2,by="subset")

sign_info=fread_FT("~/ws/2021/211124_COI_SLE_HC_run27_finaldata_Combat_PCA/PCAres_Clinical/211208_COI_PC30andClinical_PCA_PC7_sign_info.txt")

 for(iii in 1:nrow(list)){
  subset_tmp      = list$subset[iii]
  unique_tmp      = fread_FT(list$PATH_unique[iii]) %>% select(2,9:15)
  dupli_tmp       = fread_FT(list$PATH_dupli[iii]) %>% select(1:8) %>% mutate(id=take_factor(name,1,"_")) %>% select(9,2:8)
  data_tmp        = rbind(unique_tmp,dupli_tmp)
  data_tmp        = inner_join(BLM_samples_2,data_tmp,by="id")
  data_tmp$subset = subset_tmp
  if(iii==1){data_sum=data_tmp}else{data_sum=rbind(data_sum,data_tmp)}
 }
 data_sum2 = data_sum %>% pivot_longer(cols=-c(subset,id,indiv,response,time),names_to="PC",values_to="score") %>%
                          mutate(subset_PC=paste0(subset,"_",PC)) %>% left_join(.,sign_info,by="subset_PC") %>%
                          mutate(signed_score=ifelse(sign=="neg",-score,score)) # Add sign info
write.table_FT_2(data_sum2,paste0("PC_BEL/",today(),"_COI_SLE_BEL_PC1to7_sum.txt"))

####################
# Scale the score
####################
# All PC1-7
PC_signed = fread_FT("~/ws/2021/211124_COI_SLE_HC_run27_finaldata_Combat_PCA/PCA_res/SLEunique_andHC/AfterCombat/211126_COI_SLEuniqueandHC_AfterCombat_PC7_signed.txt")%>%
            select(id,subset_PC,value)%>%
            pivot_wider(names_from="subset_PC",values_from="value")  %>% column_to_rownames("id")

# Infer mean and sd
# Necessary to calculate scaled PC scores for duplicated pts
PC_mean=apply(PC_signed,2,mean,na.rm=T) %>% as.data.frame() # all zero
PC_sd=apply(PC_signed,2,sd,na.rm=T) %>% as.data.frame() %>% rownames_to_column("subset_PC")
colnames(PC_sd)[2]="sd"

# Scale PC scores of BLM samples
data_sum3=data_sum2%>%left_join(.,PC_sd,by="subset_PC")%>%
                      mutate(scaled_score=signed_score/sd)
write.table_FT_2(data_sum3,paste0("PC_BEL/",today(),"_COI_SLE_BEL_PC1to7_sum_scaled.txt"))

####################
# Test all PC1-7
####################

all_tmp=data_sum3%>%select(-c(score,sign,signed_score,sd))
colnames(all_tmp)[8]="score"

all_tmp$time = factor(all_tmp$time,levels=c("0m","6m"))
all_tmp$indiv = factor(all_tmp$indiv,levels=unique(all_tmp$indiv))

for(iii in 1:length(unique(all_tmp$subset_PC))){
    PC_tmp  = unique(all_tmp$subset_PC)[iii]
    res_tmp = all_tmp %>% filter(subset_PC==PC_tmp)
    good_tmp = res_tmp %>% filter(response=="Good")
    poor_tmp = res_tmp %>% filter(response=="Poor")
    
    all_full  = lmer(score~time+(1|indiv),data=res_tmp) 
    all_null  = lmer(score~     (1|indiv),data=res_tmp) 
    all_anova = anova(all_full,all_null)
    all_beta  = summary(all_full)$coef[2,1]
    all_tvalue= summary(all_full)$coef[2,3]
    all_P     = all_anova[2,8]

    good_full  = lmer(score~time+(1|indiv),data=good_tmp) 
    good_null  = lmer(score~     (1|indiv),data=good_tmp) 
    good_anova = anova(good_full,good_null)
    good_beta  = summary(good_full)$coef[2,1]
    good_tvalue= summary(good_full)$coef[2,3]
    good_P     = good_anova[2,8]

    poor_full  = lmer(score~time+(1|indiv),data=poor_tmp) 
    poor_null  = lmer(score~     (1|indiv),data=poor_tmp) 
    poor_anova = anova(poor_full,poor_null)
    poor_beta  = summary(poor_full)$coef[2,1]
    poor_tvalue= summary(poor_full)$coef[2,3]
    poor_P     = poor_anova[2,8]

    all_dup_tmp  = res_tmp[duplicated(res_tmp$indiv),]
    all_t_tmp    = res_tmp %>% filter(indiv%in%all_dup_tmp$indiv)
    all_t_0m     = all_t_tmp %>% filter(time=="0m")
    all_t_6m     = all_t_tmp %>% filter(time=="6m")
    all_t_res    = t.test(x=all_t_6m$score,y=all_t_0m$score,paired=T)
    all_t2       = all_t_res$statistic
    all_P2       = all_t_res$p.value

    good_dup_tmp  = good_tmp[duplicated(good_tmp$indiv),]
    good_t_tmp    = good_tmp %>% filter(indiv%in%good_dup_tmp$indiv)
    good_t_0m     = good_t_tmp %>% filter(time=="0m")
    good_t_6m     = good_t_tmp %>% filter(time=="6m")
    good_t_res    = t.test(x=good_t_6m$score,y=good_t_0m$score,paired=T)
    good_t2       = good_t_res$statistic
    good_P2       = good_t_res$p.value

    poor_dup_tmp  = poor_tmp[duplicated(poor_tmp$indiv),]
    poor_t_tmp    = poor_tmp %>% filter(indiv%in%poor_dup_tmp$indiv)
    poor_t_0m     = poor_t_tmp %>% filter(time=="0m")
    poor_t_6m     = poor_t_tmp %>% filter(time=="6m")
    poor_t_res    = t.test(x=poor_t_6m$score,y=poor_t_0m$score,paired=T)
    poor_t2       = poor_t_res$statistic
    poor_P2       = poor_t_res$p.value
　　　　　 

   all_w_res    = wilcox.test(x=all_t_6m$score,y=all_t_0m$score,paired=T)
   all_w        = all_w_res$statistic
   all_P3       = all_w_res$p.value

   good_w_res    = wilcox.test(x=good_t_6m$score,y=good_t_0m$score,paired=T)
   good_w        = good_w_res$statistic
   good_P3       = good_w_res$p.value

   poor_w_res    = wilcox.test(x=poor_t_6m$score,y=poor_t_0m$score,paired=T)
   poor_w        = poor_w_res$statistic
   poor_P3       = poor_w_res$p.value

  res_sum_tmp = data.frame(subset_PC=PC_tmp,response=c("all","good","poor"),
                         lmm_beta=c(all_beta,good_beta,poor_beta),lmm_t=c(all_tvalue,good_tvalue,poor_tvalue),lmm_P=c(all_P,good_P,poor_P),
                         pairT_t=c(all_t2,good_t2,poor_t2),pairT_P=c(all_P2,good_P2,poor_P2),
                         pairW_t=c(all_w,good_w,poor_w),pairW_P=c(all_P3,good_P3,poor_P3))
  
  if(iii==1){res_sum_sum=res_sum_tmp}else{res_sum_sum=rbind(res_sum_sum,res_sum_tmp)}
}
res_sum_sum$lmm_Q   = p.adjust(res_sum_sum$lmm_P,method="BH")
res_sum_sum$pairT_Q = p.adjust(res_sum_sum$pairT_P,method="BH")
res_sum_sum$pairW_Q = p.adjust(res_sum_sum$pairW_P,method="BH")
write.table_FT_2(res_sum_sum,paste0("PC_BEL/",today(),"_COI_SLE_BEL_PC1to7_testres_all.txt"))


#############################################################
# USM B PC4 boxplot
data_sum2 = fread_FT("PC_BEL/211208_COI_SLE_BEL_PC1to7_sum_scaled.txt")
PC_scaled_USMB = data_sum2 %>% filter(subset_PC=="A03_UnswMB_PC4") 

PC_scaled_USMB$time=factor(PC_scaled_USMB$time,levels=c("0m","6m"))
levels(PC_scaled_USMB$time)=list(`Pre-BLM`="0m",`Post-BLM`="6m")
PC_scaled_USMB$response=factor(PC_scaled_USMB$response,levels=c("Good","Poor"))
levels(PC_scaled_USMB$response)=list(`BLM (Good)`="Good",`BLM (Poor)`="Poor")

brewer.pal(3,"Greens")
# [1] "#E5F5E0" "#A1D99B" "#31A354"

 p = ggplot(data=PC_scaled_USMB,aes(x=response,y=scaled_score,fill=time))+
         geom_boxplot(position="dodge",outlier.shape=NA)+
         #geom_point(size=0.5,color="black",position=position_jitterdodge())+
         scale_fill_manual(values=c("#E5F5E0","#31A354"))+
         theme_classic()+
         theme(axis.text.x=element_text(colour="black",angle=45,hjust=1,vjust=1,size=14),
               axis.text.y=element_text(colour="black",size=14),
               axis.title.x=element_blank(),
               axis.title.y=element_text(colour="black",size=14),
               plot.title=element_blank(),
               legend.position="right",
               legend.title=element_blank(),
               legend.text=element_text(colour="black",size=14))+
         labs(y="Scaled PC score")
    
 pdf_3(paste0("Figures/",today(),"_COI_SLE_BEL_PC1to7_USMBPC4_box.pdf"),h=3.5,w=3.8)
  plot(p)
 dev.off()

##########################################################################################################################




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
 [1] ggsci_2.9          RColorBrewer_1.1-2 lme4_1.1-27.1      Matrix_1.3-4      
 [5] data.table_1.14.0  forcats_0.5.1      stringr_1.4.0      dplyr_1.0.7       
 [9] purrr_0.3.4        readr_1.4.0        tidyr_1.1.3        tibble_3.1.2      
[13] ggplot2_3.3.5      tidyverse_1.3.1   

loaded via a namespace (and not attached):
 [1] tidyselect_1.1.1 splines_4.0.2    haven_2.4.1      lattice_0.20-41 
 [5] colorspace_2.0-2 vctrs_0.3.8      generics_0.1.0   utf8_1.2.1      
 [9] rlang_0.4.11     nloptr_1.2.2.2   pillar_1.6.1     glue_1.4.2      
[13] withr_2.4.2      DBI_1.1.1        dbplyr_2.1.1     modelr_0.1.8    
[17] readxl_1.3.1     lifecycle_1.0.0  munsell_0.5.0    gtable_0.3.0    
[21] cellranger_1.1.0 rvest_1.0.0      labeling_0.4.2   ps_1.6.0        
[25] fansi_0.5.0      broom_0.7.8      Rcpp_1.0.6       scales_1.1.1    
[29] backports_1.2.1  jsonlite_1.7.2   farver_2.1.0     fs_1.5.0        
[33] digest_0.6.27    hms_1.1.0        stringi_1.6.2    grid_4.0.2      
[37] cli_2.5.0        tools_4.0.2      magrittr_2.0.1   crayon_1.4.1    
[41] pkgconfig_2.0.3  ellipsis_0.3.2   MASS_7.3-54      xml2_1.3.2      
[45] reprex_2.0.0     lubridate_1.7.10 minqa_1.2.4      assertthat_0.2.1
[49] httr_1.4.2       rstudioapi_0.13  R6_2.5.0         boot_1.3-28     
[53] nlme_3.1-152     compiler_4.0.2 












