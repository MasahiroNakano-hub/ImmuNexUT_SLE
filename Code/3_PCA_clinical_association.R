#---------------------------------------#
# 211124 Nakano
# PCA clinical assoc
#---------------------------------------#

source("my.source.R")
options(stringsAsFactors=F)
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ggsci))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(uwot))
suppressPackageStartupMessages(library(viridis))

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

##############################################
# Pie chart of disease activity
##############################################
clinical= fread_FT("data_ref/clinicaldata_lim.txt")  # Private

clinical_SLE=clinical%>%filter(disease=="1SLE")
Act_freq=table_freq(clinical_SLE$Activity)
colnames(Act_freq)[1]="Activity"
Act_freq$Activity=factor(Act_freq$Activity,levels=c("1Inactive","2LDA","3MDA","4HDA"))
levels(Act_freq$Activity)=list(`Inactive`="1Inactive",`LDA`="2LDA",`MDA`="3MDA",`HDA`="4HDA")

col=ggColorHue(n=5)[2:5]

p = ggplot()+
     geom_bar(data=Act_freq, aes(x="",y=Freq,fill=Activity),stat="identity")+
     coord_polar("y",start=0,direction=-1) +
     scale_fill_manual(values=col)+
     theme_void()+
     theme(legend.position="right",
           legend.title=element_text(colour="black",size=16),
           legend.text=element_text(colour="black",size=16))+
     guides(fill=guide_legend(title="Disease activity"))

pdf_3(paste0("PCAres_Clinical/",today(),"_COI_SLE_DA_pie.pdf"),h=5,w=5)
 plot(p)
dev.off()

# For figure1
col=c("#e0d51a","#e66557","#f77308","#cc0010")

p = ggplot()+
     geom_bar(data=Act_freq, aes(x="",y=Freq,fill=Activity),stat="identity")+
     coord_polar("y",start=0,direction=-1) +
     scale_fill_manual(values=col)+
     theme_void()+
     theme(legend.position="right",
           legend.title=element_text(colour="black",size=16),
           legend.text=element_text(colour="black",size=16))+
     guides(fill=guide_legend(title="Disease activity"))

pdf_3(paste0("Figure/",today(),"_COI_SLE_DA_pie_Fig1.pdf"),h=5,w=5)
 plot(p)
dev.off()


##############################################
# PCA-clinical
##############################################
# Each filter conditions
list = make_list("PCA_res/SLEunique_andHC/AfterCombat","_PCdata.txt")
list$subset = take_factor(list$FILE,7:8,"_")
list$filter = take_factor(list$FILE,6,"_")
list$status = take_factor(list$FILE,5,"_")

# aggregate 27 cell types, PC 1-30 data (explained variance 50~60%)
 for(kkk in 1:length(unique(list$filter))){
   filter_tmp = unique(list$filter)[kkk]
   list_tmp   = list %>% filter(filter==filter_tmp)
   
   for(iii in 1:length(unique(list_tmp$subset))){
     subset_tmp = celltype_reordered_27[iii]
     list_tmp_2 = list_tmp %>% filter(subset==subset_tmp)
     df_tmp = fread_FT(list_tmp_2$PATH) %>% select(2,5,9:38)
     if(iii==1){res_sum=df_tmp}else{res_sum=rbind(res_sum,df_tmp)}
     }
    write.table_FT_2(res_sum,paste0("PCA_res/SLEunique_andHC/AfterCombat/",today(),"_COI_SLEuniqueandHC_AfterCombat_PC30_ressum.txt"))
 }

##############################################
####### PCA-clinical
clinical= fread_FT("data_ref/clinicaldata_lim.txt")
PCA_res = res_sum %>%
          pivot_longer(col=-c(id,subset),names_to="PC",values_to="value") %>%
          mutate(subset_PC = paste0(subset,"_",PC)) %>% select(-c(subset,PC)) %>%
          pivot_wider(names_from="subset_PC",values_from="value") %>%
          column_to_rownames("id") %>% t() %>% as.data.frame()

# Scale PC score
PCA_res_Z = apply(PCA_res,1,scale) %>% as.data.frame()
rownames(PCA_res_Z) = colnames(PCA_res)
dim(PCA_res_Z)
# [1] 225 810

# Disease state and activity adjusted for covariates
  for(kkk in 1:ncol(PCA_res_Z)){
     PC_tmp  = colnames(PCA_res_Z)[kkk]
     data_1 = PCA_res_Z  %>% select(all_of(kkk)) %>% rownames_to_column("id")

    for(mmm in c(13:16)){
      if(mmm==13){
       clinical_tmp = colnames(clinical)[mmm]
       data_2 = clinical %>% select(1:10,all_of(mmm))
       data_tmp = left_join(data_1,data_2,by="id")
       colnames(data_tmp)[2] = "y"
       colnames(data_tmp)[12] = "x"
       data_tmp = data_tmp %>% filter(!is.na(y)) %>% filter(!is.na(x)) 
       res_tmp = summary(lm(y~x+age+gender,data=data_tmp))$coefficients %>% as.data.frame() %>% .[2,]
       res_tmp$clinical = clinical_tmp
       res_tmp$subset_PC = PC_tmp
       colnames(res_tmp)[1:4] = c("beta","se","t","P")
       res_tmp    = res_tmp %>% select(6,5,1:4)
       res_sum_tmp=res_tmp
      }else{
       clinical_tmp = colnames(clinical)[mmm]
       data_2 = clinical %>% select(1:10,all_of(mmm))
       data_tmp = left_join(data_1,data_2,by="id")
       colnames(data_tmp)[2] = "y"
       colnames(data_tmp)[12] = "x"
       data_tmp = data_tmp %>% filter(!is.na(y)) %>% filter(!is.na(x)) 
       res_tmp = summary(lm(y~x+age+gender+PSLmg+HCQ+MMF+TAC,data=data_tmp))$coefficients %>% as.data.frame() %>% .[2,]
       res_tmp$clinical = clinical_tmp
       res_tmp$subset_PC = PC_tmp
       colnames(res_tmp)[1:4] = c("beta","se","t","P")
       res_tmp    = res_tmp %>% select(6,5,1:4)
       res_sum_tmp=rbind(res_sum_tmp,res_tmp)
       }
     }
   if(kkk==1){res_sum=res_sum_tmp}else{res_sum=rbind(res_sum,res_sum_tmp)}
   }

   res_sum$subset = take_factor(res_sum$subset_PC,1:2,"_")
   res_sum$PC 　　 = take_factor(res_sum$subset_PC,3,"_")
   res_sum_2 = res_sum %>% select(subset_PC,subset,PC,clinical,beta,se,t,P)
   write.table_FT_2(res_sum_2,paste0("PCAres_Clinical/",today(),"_COI_PC30andClinical_PCA_DA_lm.txt"))


clinical= clinical %>% filter(disease=="1SLE") %>%
          mutate(extrarenal=ifelse(neurological=="1Yes"|cardiorespiratory=="1Yes"|gastrointestinal=="1Yes","1Yes","0No"))

# Organ multiple linear regression model adjusted for covariates
  for(kkk in 1:ncol(PCA_res_Z)){
     PC_tmp  = colnames(PCA_res_Z)[kkk]
     data_1 = PCA_res_Z  %>% select(all_of(kkk)) %>% rownames_to_column("id")

       data_2 = clinical %>% select(1:10,17:26)
       data_tmp = inner_join(data_1,data_2,by="id") %>% .[!apply(.,1,anyNA),]
       colnames(data_tmp)[2] = "y"
       res_tmp = summary(lm(y~constitutional+mucocutaneous+musculoskeletal+renal+extrarenal+hematological+serological+HCQ+MMF+TAC+PSLmg+age+gender,data=data_tmp))$coefficients %>% as.data.frame() %>% .[2:8,]
       res_tmp$subset_PC = PC_tmp
       colnames(res_tmp)[1:4] = c("beta","se","t","P")
       res_tmp    = res_tmp %>% select(5,1:4) %>% rownames_to_column("clinical") 
       res_sum_tmp=res_tmp     
   if(kkk==1){res_sum=res_sum_tmp}else{res_sum=rbind(res_sum,res_sum_tmp)}
   }

   res_sum$subset = take_factor(res_sum$subset_PC,1:2,"_")
   res_sum$PC 　　 = take_factor(res_sum$subset_PC,3,"_")
   res_sum_2 = res_sum %>% select(subset_PC,subset,PC,clinical,beta,se,t,P)
   res_sum_2$clinical = gsub("1Yes","",res_sum_2$clinical)
   write.table_FT_2(res_sum_2,paste0("PCAres_Clinical/",today(),"_COI_PC30andClinical_PCA_Organ_lm.txt"))

# Tx adjusted for covariates
for(kkk in 1:ncol(PCA_res_Z)){
     PC_tmp  = colnames(PCA_res_Z)[kkk]
     data_1 = PCA_res_Z  %>% select(all_of(kkk)) %>% rownames_to_column("id")

       data_2 = clinical %>% select(1:10,12)
       data_tmp = inner_join(data_1,data_2,by="id") %>% .[!apply(.,1,anyNA),]
       colnames(data_tmp)[2] = "y"
       res_tmp = summary(lm(y~HCQ+MMF+TAC+PSLmg+Activity+age+gender,data=data_tmp))$coefficients %>% as.data.frame() %>% .[2:4,]
       res_tmp$subset_PC = PC_tmp
       colnames(res_tmp)[1:4] = c("beta","se","t","P")
       res_tmp    = res_tmp %>% select(5,1:4) %>% rownames_to_column("clinical") 
       res_sum_tmp=res_tmp     
   if(kkk==1){res_sum=res_sum_tmp}else{res_sum=rbind(res_sum,res_sum_tmp)}
   }

   res_sum$subset = take_factor(res_sum$subset_PC,1:2,"_")
   res_sum$PC 　　 = take_factor(res_sum$subset_PC,3,"_")
   res_sum_2 = res_sum %>% select(subset_PC,subset,PC,clinical,beta,se,t,P)
   res_sum_2$clinical = gsub("1Yes","",res_sum_2$clinical)
   write.table_FT_2(res_sum_2,paste0("PCAres_Clinical/",today(),"_COI_PC30andClinical_PCA_Tx_lm.txt"))

##############################################
####### Aggregate and check significant associations
##############################################

Ac_res     = fread_FT("PCAres_Clinical/211126_COI_PC30andClinical_PCA_DA_lm.txt")
Organ_res  = fread_FT("PCAres_Clinical/211126_COI_PC30andClinical_PCA_Organ_lm.txt")
Tx_res     = fread_FT("PCAres_Clinical/211126_COI_PC30andClinical_PCA_Tx_lm.txt") 

res_sum = bind_rows(list(Ac_res,Organ_res,Tx_res))
res_sum$PC = factor(res_sum$PC,levels=unique(res_sum$PC))
res_sum$subset_clinical = paste0(res_sum$subset,"_",res_sum$clinical)
res_sum$subset_clinical = factor(res_sum$subset_clinical,levels=unique(res_sum$subset_clinical))
res_sum$Q   =p.adjust(res_sum$P,method="BH")

write.table_FT_2(res_sum,paste0("PCAres_Clinical/",today(),"_COI_PC30andClinical_PCA_All_lm.txt"))

# significant results
res_sum_sig005  = res_sum %>% filter(Q<0.05)   
data_tmp = table(res_sum_sig005$PC) %>% as.data.frame()
cutoff=105/30

p =  ggplot()+
  geom_bar(data=data_tmp,aes(x=Var1,y=Freq),stat="identity")+theme_classic()+
  geom_hline(yintercept=cutoff,col="red")+
  theme(axis.text.x=element_text(colour="black",angle=45,hjust=1,vjust=1,size=15),
        axis.text.y=element_text(colour="black",size=15),
        axis.title.x=element_blank(),
        axis.title.y=element_text(colour="black",size=15))+
  labs(y="Number of significant associations")
 
pdf_3(paste0("PCAres_Clinical/",today(),"_COI_PC30andClinical_PCA_All_lm_allsignum.pdf"),h=5,w=10)
  plot(p)
dev.off()  

##############################################################################################################################
# Define the sign of each PC based on lm result
##############################################################################################################################
#####################
# Limit to PC1-7
res_sum=fread_FT("PCAres_Clinical/211126_COI_PC30andClinical_PCA_All_lm.txt")
res_sum$PC=factor(res_sum$PC,levels=unique(res_sum$PC))
include_tmp=levels(res_sum$PC)[1:7]
res_sum_2 = res_sum %>% filter(PC%in%include_tmp)

# Define disease-state and activity PCs
data_all = res_sum_2 %>% 
           filter(clinical%in%c("InactivevsHC","HDAvsInactive")) %>% 
           select(subset_PC,subset,PC,clinical,beta,P,Q) %>% 
           pivot_wider(names_from="clinical", values_from=c("beta","P","Q")) %>% column_to_rownames("subset_PC") %>% 
           mutate(siggroup005=ifelse(Q_InactivevsHC<0.05&Q_HDAvsInactive<0.05&(beta_InactivevsHC*beta_HDAvsInactive>0),"Both_concordant",
                              ifelse(Q_InactivevsHC<0.05&Q_HDAvsInactive<0.05&(beta_InactivevsHC*beta_HDAvsInactive<0),"Both_discordant",
                              ifelse(Q_HDAvsInactive<0.05,"Activity",
                              ifelse(Q_InactivevsHC<0.05,"State","Others")))))
dim(data_all)
# [1] 189   9
table(data_all$siggroup005)
#       Activity Both_concordant Both_discordant          Others           State 
#             16               9               3             133              28 

#### Signed beta
data_1 = data_all %>% filter(P_InactivevsHC<P_HDAvsInactive)
data_11 = data_1 %>% filter(beta_InactivevsHC<0)
data_12 = data_1 %>% filter(beta_InactivevsHC>0)

statedominant_neg=rownames(data_11) # 54 These PC scores should be inversed
statedominant_pos=rownames(data_12) # 46 These PC scores should be kept
data_11[,3:6] = -data_11[,3:6]

# Define activity dominant PCs
data_2 = data_all %>% filter(P_InactivevsHC>P_HDAvsInactive)
data_21 = data_2 %>% filter(beta_HDAvsInactive<0)
data_22 = data_2 %>% filter(beta_HDAvsInactive>0)

activitydominant_neg=rownames(data_21) # 45 These PC scores should be inversed
activitydominant_pos=rownames(data_22) # 44 These PC scores should be kept
data_21[,3:6] = -data_21[,3:6]

##### Save sign info for all PCs
neg1=data.frame(subset_PC=statedominant_neg,sign="neg")
neg2=data.frame(subset_PC=activitydominant_neg,sign="neg")
pos1=data.frame(subset_PC=statedominant_pos,sign="pos")
pos2=data.frame(subset_PC=activitydominant_pos,sign="pos")
sign_list=bind_rows(list(neg1,neg2,pos1,pos2))
write.table_FT_2(sign_list,paste0("PCAres_Clinical/",today(),"_COI_PC30andClinical_PCA_PC7_sign_info.txt"))

# bind
data_signed=bind_rows(list(data_11,data_12,data_21,data_22))
data_signed$PC=factor(data_signed$PC,levels=include_tmp)
data_signed$subset=factor(data_signed$subset,levels=celltype_reordered_27)
data_signed=data_signed[order(data_signed$PC),]
data_signed=data_signed[order(data_signed$subset),]
data_signed=data_signed%>%rownames_to_column("subset_PC")
dim(data_signed)
# [1] 189  10
write.table_FT_2(data_signed,paste0("PCAres_Clinical/",today(),"_COI_PC30andClinical_PCA_PC7_lm_signed_betaPQ.txt"))

# State, activity PCs
data_signed_sig=data_signed%>%filter(siggroup005%in%c("State","Activity","Both_concordant"))
data_signed_act=data_signed%>%filter(siggroup005%in%c("Activity","Both_concordant"))

##############################################################################################################################
# Supple Table
# 27 celltype PC 1-7 clinical lmres signed
##############################################################################################################################
res_sum_inverse=res_sum_2%>%filter(subset_PC%in%c(statedominant_neg,activitydominant_neg))    # These PC scores should be inversed
res_sum_keep   =res_sum_2%>%filter(subset_PC%in%c(statedominant_pos,activitydominant_pos))    # These PC scores should be kept
res_sum_inverse$beta=-res_sum_inverse$beta
res_sum_signed=rbind(res_sum_keep,res_sum_inverse)

res_sum_signed$subset=factor(res_sum_signed$subset,levels=celltype_reordered_27)
res_sum_signed=res_sum_signed[order(res_sum_signed$PC),]
res_sum_signed=res_sum_signed[order(res_sum_signed$subset),]
write.table_FT_2(res_sum_signed,paste0("PCAres_Clinical/",today(),"_COI_27subset_PC1to7_lmres_signed_all.txt"))

res_sum_signed2 = res_sum_signed %>%
                 left_join(.,label_df,by="subset") %>% 
                 select(label,PC,clinical,beta,se,P,Q)

res_sum_signed2$beta = formatC(res_sum_signed2$beta,digits=2)
res_sum_signed2$se   = formatC(res_sum_signed2$se,digits=2)
res_sum_signed2$P    = formatC(res_sum_signed2$P,digits=2)
res_sum_signed2$Q    = formatC(res_sum_signed2$Q,digits=2)

colnames(res_sum_signed2)=c("Cell type","PC","Clinical parameter","Effect size","Standard error","Pvalue","FDR")

write.table_FT_2(res_sum_signed2,paste0(today(),"_COI_27subset_PC1to7_lmres_signed.txt"))

##############################
# signed lmres for State/Activity, organ(activity), Tx(activity)
##############################
res_sum_signed_activity = res_sum_signed %>% filter(clinical%in%c("InactivevsHC","LDAvsInactive","MDAvsInactive","HDAvsInactive")) %>%
                                             filter(subset_PC%in%data_signed_sig$subset_PC)
write.table_FT_2(res_sum_signed_activity,paste0("PCAres_Clinical/",today(),"_COI_27subset_PC1to7_lmres_signed_stateactivity.txt"))

res_sum_signed_organ = res_sum_signed %>% filter(clinical%in%c("constitutional","mucocutaneous","musculoskeletal","renal","extrarenal","hematological","serological")) %>%
                                          filter(subset_PC%in%data_signed_act$subset_PC)
write.table_FT_2(res_sum_signed_organ,paste0("PCAres_Clinical/",today(),"_COI_27subset_PC1to7_lmres_signed_organactivity.txt"))

res_sum_signed_Tx = res_sum_signed %>% filter(clinical%in%c("HCQ","MMF","TAC")) %>%
                                       filter(subset_PC%in%data_signed_act$subset_PC)
write.table_FT_2(res_sum_signed_Tx,paste0("PCAres_Clinical/",today(),"_COI_27subset_PC1to7_lmres_signed_Txactivity.txt"))

##############################################################################################################################
# PC 1-7 score signed
##############################################################################################################################
PC_sum = fread_FT("PCA_res/SLEunique_andHC/AfterCombat/211125_COI_SLEuniqueandHC_AfterCombat_PC30_ressum.txt")[,1:9] %>% 
         pivot_longer(col=-c(id,subset),names_to="PC",values_to="value") %>%
         mutate(subset_PC=paste0(subset,"_",PC))

PC_inverse = PC_sum %>% filter(subset_PC%in%c(statedominant_neg,activitydominant_neg))
PC_keep    = PC_sum %>% filter(subset_PC%in%c(statedominant_pos,activitydominant_pos))
PC_inverse$value = -(PC_inverse$value)
PC_signed = rbind(PC_keep,PC_inverse)

PC_signed$PC    =factor(PC_signed$PC,levels=include_tmp)
PC_signed$subset=factor(PC_signed$subset,levels=celltype_reordered_27)
PC_signed=PC_signed[order(PC_signed$PC),]
PC_signed=PC_signed[order(PC_signed$subset),]

write.table_FT_2(PC_signed,paste0("PCA_res/SLEunique_andHC/AfterCombat/",today(),"_COI_SLEuniqueandHC_AfterCombat_PC7_signed.txt"))

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
[1] grid      stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] viridis_0.6.1        viridisLite_0.4.0    circlize_0.4.13     
 [4] ComplexHeatmap_2.6.2 pheatmap_1.0.12      uwot_0.1.10         
 [7] Matrix_1.3-4         ggsci_2.9            RColorBrewer_1.1-2  
[10] sva_3.38.0           BiocParallel_1.24.1  genefilter_1.72.1   
[13] mgcv_1.8-36          nlme_3.1-152         edgeR_3.32.1        
[16] limma_3.46.0         data.table_1.14.0    forcats_0.5.1       
[19] stringr_1.4.0        dplyr_1.0.7          purrr_0.3.4         
[22] readr_1.4.0          tidyr_1.1.3          tibble_3.1.2        
[25] ggplot2_3.3.5        tidyverse_1.3.1     

loaded via a namespace (and not attached):
 [1] matrixStats_0.59.0   fs_1.5.0             lubridate_1.7.10    
 [4] bit64_4.0.5          RcppAnnoy_0.0.18     httr_1.4.2          
 [7] tools_4.0.2          backports_1.2.1      utf8_1.2.1          
[10] R6_2.5.0             irlba_2.3.3          DBI_1.1.1           
[13] BiocGenerics_0.36.1  colorspace_2.0-2     GetoptLong_1.0.5    
[16] withr_2.4.2          gridExtra_2.3        tidyselect_1.1.1    
[19] bit_4.0.4            compiler_4.0.2       cli_2.5.0           
[22] rvest_1.0.0          Biobase_2.50.0       Cairo_1.5-12.2      
[25] xml2_1.3.2           labeling_0.4.2       scales_1.1.1        
[28] digest_0.6.27        pkgconfig_2.0.3      dbplyr_2.1.1        
[31] fastmap_1.1.0        GlobalOptions_0.1.2  rlang_0.4.11        
[34] readxl_1.3.1         rstudioapi_0.13      RSQLite_2.2.7       
[37] shape_1.4.6          generics_0.1.0       farver_2.1.0        
[40] jsonlite_1.7.2       magrittr_2.0.1       Rcpp_1.0.6          
[43] munsell_0.5.0        S4Vectors_0.28.1     fansi_0.5.0         
[46] lifecycle_1.0.0      stringi_1.6.2        MASS_7.3-54         
[49] blob_1.2.1           parallel_4.0.2       crayon_1.4.1        
[52] lattice_0.20-41      haven_2.4.1          splines_4.0.2       
[55] annotate_1.68.0      hms_1.1.0            locfit_1.5-9.4      
[58] ps_1.6.0             pillar_1.6.1         rjson_0.2.20        
[61] codetools_0.2-18     stats4_4.0.2         reprex_2.0.0        
[64] XML_3.99-0.6         glue_1.4.2           modelr_0.1.8        
[67] png_0.1-7            vctrs_0.3.8          cellranger_1.1.0    
[70] gtable_0.3.0         clue_0.3-59          assertthat_0.2.1    
[73] cachem_1.0.5         xtable_1.8-4         broom_0.7.8         
[76] survival_3.2-11      AnnotationDbi_1.52.0 memoise_2.0.0       
[79] IRanges_2.24.1       cluster_2.1.2        ellipsis_0.3.2     





































