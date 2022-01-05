#---------------------------------------#
# 211222 Nakano
# variance partitioning analysis + jack
#---------------------------------------#


source("my.source.R")
options(stringsAsFactors=F)
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ggsci))
suppressPackageStartupMessages(library(lme4))
suppressPackageStartupMessages(library(variancePartition))

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
# Cell type PC 1-7 list
##############################################
# All BeforeCombat
list1 = make_list("PCA_res/SLEinc_andHC/BeforeCombat/tmp_Rdata","_PCAres_tmp1.Rdata")
list1$subset = take_factor(list1$FILE,7:8,"_")
list1$type   = "BeforeCombat"

# All AfterCombat
list2 = make_list("PCA_res/SLEinc_andHC/AfterCombat/tmp_Rdata","_PCAres_tmp1.Rdata")
list2$subset = take_factor(list2$FILE,7:8,"_")
list2$type   = "AfterCombat"

list3=rbind(list1,list2)

condition_clinical= read.table_FT("data/SLEinc_andHC/211124_COI_SLEinc_andHC_27subsets_cond_wclinical.txt") %>% 
                    mutate(Immunosuppressant=ifelse(HCQ=="0No"&MMF=="0No"&TAC=="0No","0No",
                                             ifelse(HCQ=="1Yes"&MMF=="0No"&TAC=="0No","1HCQ",
                                             ifelse(HCQ=="0No"&MMF=="1Yes"&TAC=="0No","2MMF",
                                             ifelse(HCQ=="0No"&MMF=="0No"&TAC=="1Yes","3TAC",
                                             ifelse(HCQ=="1Yes"&MMF=="1Yes"&TAC=="0No","4HCQMMF",
                                             ifelse(HCQ=="1Yes"&MMF=="0No"&TAC=="1Yes","5HCQTAC",
                                             ifelse(HCQ=="0No"&MMF=="1Yes"&TAC=="1Yes","6MMFTAC","7HCQMMFTAC")))))))) 

name_list = condition_clinical %>%select(name,id,disease,subset,lotseq,age,gender,Activity,constitutional,mucocutaneous,musculoskeletal,renal,extrarenal,hematological,serological,PSLmg,Immunosuppressant,HCQ,MMF,TAC) 
name_list$disease  = factor(name_list$disease,levels=c("0HC","1SLE"))
name_list$lotseq 　= factor(name_list$lotseq,levels=c("Lot3_HiSeq","Lot1_HiSeq","Lot2_HiSeq","Lot4_HiSeq","Lot5_HiSeq","Lot5_NovaSeq"))
name_list$Activity = factor(name_list$Activity,levels=c("0HC","1Inactive","2LDA","3MDA","4HDA"))
name_list$Immunosuppressant = factor(name_list$Immunosuppressant,levels=sort(unique(name_list$Immunosuppressant)))

name_list_SLE=name_list%>%filter(disease=="1SLE")

#################################
# PVCA batch All (Before/After)
# Random effects only: REML
#################################
for (nnn in 1:length(unique(list3$type))){
    type_tmp=unique(list3$type)[nnn]
    list_tmp=list3%>%filter(type==type_tmp)

 for (kkk in 1:nrow(list_tmp)){
    subset_tmp=list_tmp$subset[kkk]
    load(list_tmp$PATH[kkk])
    cumvar = summary(result)$importance[,1:7] %>% as.data.frame()
    cumvar_sum = cumvar[3,7]

    df1_tmp = as.data.frame(result$x) %>% select(1:7) %>%
              rownames_to_column("name") %>%
              right_join(name_list,.,by="name")

  for (iii in 21:27){
     PC_tmp   = colnames(df1_tmp)[iii]
     Prop_raw = cumvar %>% select(all_of(PC_tmp)) %>%.[2,]
     Prop_std = Prop_raw/cumvar_sum

     data_tmp = df1_tmp %>% select(1:20,all_of(PC_tmp))
     colnames(data_tmp)[ncol(data_tmp)] = "value"
     data_tmp_SLE = data_tmp %>% filter(disease=="1SLE")
     data_tmp_HC  = data_tmp %>% filter(disease=="0HC")
     info_table_tmp = data.frame(subset=subset_tmp,PC=PC_tmp,Totalnum=nrow(data_tmp),SLEnum=nrow(data_tmp_SLE),HCnum=nrow(data_tmp_HC))
     
     # lmm disease batch
     Rm1ML  = lmer(value ~ (1|disease)+(1|lotseq), data_tmp, REML = TRUE, verbose = FALSE, na.action = na.omit)

     Var_Random_effect = as.numeric(VarCorr(Rm1ML))
     Var_Residual = attr(VarCorr(Rm1ML), "sc")^2
     Var_Ramdom   = c(Var_Random_effect,Var_Residual)
     Sum_Random   = sum(Var_Ramdom)
     Var_Ramdom_std = Var_Ramdom/Sum_Random

     Var_std = Var_Ramdom_std %>% as.data.frame() %>%t() %>% as.data.frame()
     colnames(Var_std) = c("Batch","Disease","residual")
     rownames(Var_std) = PC_tmp
     Var_std_weighed = Var_std*Prop_std
     
     if(iii==21){Var_std_weighed_sum=Var_std_weighed}else{Var_std_weighed_sum=rbind(Var_std_weighed_sum,Var_std_weighed)}
    }

     Var_std_weighed_sum2=apply(Var_std_weighed_sum,2,sum) %>% as.data.frame() %>% t() %>% as.data.frame()
     Var_std_weighed_sum2$subset=subset_tmp
     if(kkk==1){Var_std_weighed_sum3=Var_std_weighed_sum2}else{Var_std_weighed_sum3=rbind(Var_std_weighed_sum3,Var_std_weighed_sum2)}
  }

 Batch_sum=Var_std_weighed_sum3 %>% select(subset,Batch,Disease)
 write.table_FT_2(Batch_sum,paste0("PCAres_Clinical/PVCA/",today(),"_COI_PC30andClinical_",type_tmp,"_PVCAbatch.txt"))

 Batch_sum$subset = factor(Batch_sum$subset,levels=celltype_reordered_27)

 Batch_sum2 = Batch_sum %>% pivot_longer(col=-subset,names_to="parameter",values_to="value")
 Batch_sum2$parameter = factor(Batch_sum2$parameter,levels=unique(Batch_sum2$parameter))
 Batch_sum2$subset = factor(Batch_sum2$subset,levels=celltype_reordered_27)

 col_list = c("#C77CFF","#7CAE00")

 p = ggplot()+
    geom_bar(data=Batch_sum2, aes(x=subset,y=value,fill=parameter),
           stat="identity",position="stack")+
    theme_classic()+
    scale_fill_manual(values = col_list)+
    theme(axis.text.x=element_text(colour="black",angle=45,hjust=1,vjust=1,size=16),
           axis.text.y=element_text(colour="black",size=16),
           axis.title.x=element_blank(),
           axis.title.y=element_text(colour="black",size=16),
           plot.title=element_blank(),
           legend.position="right",
           legend.title=element_blank(),
           legend.text=element_text(colour="black",size=16))+
     scale_x_discrete(labels= label)+
     ylim(0,0.65)+
     labs(y="Variance explained")
    
 pdf_3(paste0("PCAres_Clinical/PVCA/",today(),"_COI_PC30andClinical_",type_tmp,"_PVCAbatch.pdf"),h=5,w=10)
  plot(p)
 dev.off()
}

#################################
# Sum aquare: Discovery cohort
#################################
for (kkk in 1:nrow(list)){
    subset_tmp=list$subset[kkk]
    load(list$PATH[kkk])
    cumvar = summary(result)$importance[,1:7] %>% as.data.frame()
    cumvar_sum = cumvar[3,7]

    df1_tmp = as.data.frame(result$x) %>% select(1:7) %>%
              rownames_to_column("name") %>%
              right_join(name_list,.,by="name")

  for (iii in 21:27){
     PC_tmp   = colnames(df1_tmp)[iii]
     Prop_raw = cumvar %>% select(all_of(PC_tmp)) %>%.[2,]
     Prop_std = Prop_raw/cumvar_sum

     data_tmp = df1_tmp %>% select(1:20,all_of(PC_tmp))
     colnames(data_tmp)[ncol(data_tmp)] = "value"
     data_tmp_SLE = data_tmp %>% filter(disease=="1SLE")
     data_tmp_HC  = data_tmp %>% filter(disease=="0HC")
     info_table_tmp = data.frame(subset=subset_tmp,PC=PC_tmp,Totalnum=nrow(data_tmp),SLEnum=nrow(data_tmp_SLE),HCnum=nrow(data_tmp_HC))

     Total_value = data_tmp[,ncol(data_tmp)]
     SLE_value   = data_tmp_SLE[,ncol(data_tmp_SLE)]
     HC_value    = data_tmp_HC[,ncol(data_tmp_HC)]
     
     Var_Total = var(Total_value)
     Var_SLE   = var(SLE_value)
     Var_HC    = var(HC_value)

     SST     = Var_Total*(nrow(data_tmp)-1)
     SSW_SLE = Var_SLE*(nrow(data_tmp_SLE)-1)
     SSW_HC  = Var_HC*(nrow(data_tmp_HC)-1) 

     SSB = (((mean(SLE_value)-mean(Total_value))^2)*nrow(data_tmp_SLE)) + (((mean(HC_value)-mean(Total_value))^2)*nrow(data_tmp_HC))
     paste0("SST=SSW_SLE+SSW_HC+SSB; ",all.equal(SSW_SLE+SSW_HC+SSB,SST))
     aov_table = data.frame(SST=SST,SSB=SSB,SSW_SLE=SSW_SLE,SSW_HC=SSW_HC)
     aov_table_std = aov_table/SST
     #var_table = data.frame(Var_Total=Var_Total,Var_SLE=Var_SLE,Var_HC=Var_HC)
     aov_table_weighed = aov_table_std*Prop_std
     rownames(aov_table_weighed) = PC_tmp

     if(iii==21){aov_table_weighed_sum=aov_table_weighed}else{aov_table_weighed_sum=rbind(aov_table_weighed_sum,aov_table_weighed)}
  } 
  aov_table_weighed_sum2=apply(aov_table_weighed_sum,2,sum) %>% as.data.frame() %>% t() %>% as.data.frame()
  aov_table_weighed_sum2$subset=subset_tmp
  if(kkk==1){aov_table_weighed_sum3=aov_table_weighed_sum2}else{aov_table_weighed_sum3=rbind(aov_table_weighed_sum3,aov_table_weighed_sum2)}
}

aov_sum=aov_table_weighed_sum3 %>% select(subset,SSW_SLE,SSW_HC,SSB)
aov_sum_ordered=aov_sum[order(aov_sum$SSW_SLE),]
write.table_FT_2(aov_sum_ordered,paste0("PCAres_Clinical/PVCA/",today(),"_COI_PC30andClinical_anova.txt"))

aov_sum_ordered$subset = factor(aov_sum_ordered$subset,levels=aov_sum_ordered$subset)

aov_sum2 = aov_sum_ordered %>% pivot_longer(col=-subset,names_to="parameter",values_to="value")
aov_sum2$parameter = factor(aov_sum2$parameter,levels=rev(unique(aov_sum2$parameter)))
aov_sum2$subset = factor(aov_sum2$subset,levels=levels(aov_sum_ordered$subset))

# Figure
col = c("#7CAE00","#F8766D","#00BFC4")
p1 = ggplot()+
    geom_bar(data=aov_sum2,aes(x=subset,y=value,fill=factor(parameter,labels=c("Between HC and SLE","Within HC","Within SLE"))),
           stat="identity",position="stack")+
    theme_classic()+
    scale_fill_manual(values=col)+
    theme(axis.text.x=element_text(colour="black",angle=45,hjust=1,vjust=1,size=16),
           axis.text.y=element_text(colour="black",size=16),
           axis.title.x=element_blank(),
           axis.title.y=element_text(colour="black",size=16),
           plot.title=element_blank(),
           legend.position="right",
           legend.title=element_blank(),
           legend.text=element_text(colour="black",size=16))+
     scale_x_discrete(labels= label)+
    labs(y="Proportion of sum squared")
    
pdf_3(paste0("Figure/",today(),"_COI_PC30andClinical_anova.pdf"),h=4.7,w=10.5)
 plot(p1)
dev.off()

#################################
# VarPart Total Discovery dataset
# Mixed model: ML
#################################
# Unique AfterCombat: read the list
list = fread_FT("tmp_job_list/211128_COI_SLEuniqueandHC_27subsets_AfterCombat_PCres_list.txt")

for (kkk in 1:nrow(list)){
    subset_tmp=list$subset[kkk]
    load(list$PATH[kkk])
    cumvar = summary(result)$importance[,1:7] %>% as.data.frame()
    cumvar_sum = cumvar[3,7]

    Matrix=as.data.frame(result$x) %>% select(1:7) %>% t() %>% as.data.frame()

    Metadata=name_list_SLE%>%filter(name%in%colnames(Matrix))%>%select(-id) # SLE only
    Matrix=Matrix[,Metadata$name]

    Metadata$name=factor(Metadata$name,levels=colnames(Matrix))
    Metadata=Metadata[order(Metadata$name),]
    Metadata$name=as.character(Metadata$name)

    all.equal(colnames(Matrix),Metadata$name) %>% print()
   
    form = ~ (1|Activity)+PSLmg+(1|Immunosuppressant)+age+(1|gender)
    varPart = fitExtractVarPartModel( Matrix, form, Metadata ) %>% as.data.frame()
    varPart_weighed=varPart

  for (iii in 1:7){
     PC_tmp   = rownames(varPart)[iii]
     Prop_raw = cumvar %>% select(all_of(PC_tmp)) %>%.[2,]
     Prop_std = Prop_raw/cumvar_sum
     varPart_weighed[iii,]=varPart[iii,]*Prop_std
   }
     colnames(varPart_weighed) = c("Activity","Sex","Immunosuppressant","PSL","Age","residual")
     varPart_weighed_sum=apply(varPart_weighed,2,sum) %>% as.data.frame() %>% t() %>% as.data.frame()
     varPart_weighed_sum$subset=subset_tmp
     if(kkk==1){varPart_weighed_sum2=varPart_weighed_sum}else{varPart_weighed_sum2=rbind(varPart_weighed_sum2,varPart_weighed_sum)}
  }

Total_sum=varPart_weighed_sum2 %>% select(subset,Activity,PSL,Immunosuppressant,Age,Sex)
Total_sum_ordered=Total_sum[order(Total_sum$Activity),]
write.table_FT_2(Total_sum_ordered,paste0("PCAres_Clinical/variancePartition/",today(),"_COI_PC30andClinical_varPart_Total_auto.txt"))

Total_sum_ordered$subset = factor(Total_sum_ordered$subset,levels=Total_sum_ordered$subset)
Total_sum2 = Total_sum_ordered %>% pivot_longer(col=-subset,names_to="parameter",values_to="value")
Total_sum2$parameter = factor(Total_sum2$parameter,levels=rev(unique(Total_sum2$parameter)))
levels(Total_sum2$parameter)=list(`Sex`="Sex",`Age`="Age",`Immunosuppressant`="Immunosuppressant",`PSL`="PSL",`Disease activity`="Activity")
Total_sum2$subset = factor(Total_sum2$subset,levels=levels(Total_sum_ordered$subset))

col2 = rev(pal_npg("nrc")(5))

# Figure
p2 = ggplot()+
    geom_bar(data=Total_sum2, aes(x=subset,y=value,fill=parameter),
           stat="identity",position="stack")+
    theme_classic()+
    scale_fill_manual(values = col2)+
    theme(axis.text.x=element_text(colour="black",angle=45,hjust=1,vjust=1,size=16),
           axis.text.y=element_text(colour="black",size=16),
           axis.title.x=element_blank(),
           axis.title.y=element_text(colour="black",size=16),
           plot.title=element_blank(),
           legend.position="right",
           legend.title=element_blank(),
           legend.text=element_text(colour="black",size=16))+
     scale_x_discrete(labels= label)+
     scale_y_continuous(breaks= c(0,0.1,0.2))+
    labs(y="Variance explained")
    
pdf_3(paste0("Figure/",today(),"_COI_PC30andClinical_varPart_Total_auto.pdf"),h=4.7,w=10.5)
 plot(p2)
dev.off()

mean(Total_sum_ordered$Activity)
# [1] 0.07554302

Total_sum_ordered$Tx_total=Total_sum_ordered$PSL + Total_sum_ordered$Immunosuppressant
mean(Total_sum_ordered$Tx_total)
# [1] 0.04312755

Total_sum_ordered$ratio=Total_sum_ordered$Activity/Total_sum_ordered$Tx_total
mean(Total_sum_ordered$ratio)
# [1] 2.862982

#################################
# VarPart Organ Discovery dataset
# Mixed model: ML
#################################

# variancePartition version
for (kkk in 1:nrow(list)){
    subset_tmp=list$subset[kkk]
    load(list$PATH[kkk])
    cumvar = summary(result)$importance[,1:7] %>% as.data.frame()
    cumvar_sum = cumvar[3,7]

    Matrix=as.data.frame(result$x) %>% select(1:7) %>% t() %>% as.data.frame()

    Metadata=name_list_SLE%>%filter(name%in%colnames(Matrix))%>%select(-id) # SLE only
    Matrix=Matrix[,Metadata$name]

    Metadata$name=factor(Metadata$name,levels=colnames(Matrix))
    Metadata=Metadata[order(Metadata$name),]
    Metadata$name=as.character(Metadata$name)

    all.equal(colnames(Matrix),Metadata$name) %>% print()
   
    form = ~ (1|constitutional)+(1|mucocutaneous)+(1|musculoskeletal)+(1|renal)+(1|extrarenal)+(1|hematological)+(1|serological)+PSLmg+(1|HCQ)+(1|MMF)+(1|TAC)+age+(1|gender)
    varPart = fitExtractVarPartModel( Matrix, form, Metadata ) %>% as.data.frame()
    varPart_weighed=varPart

  for (iii in 1:7){
     PC_tmp   = rownames(varPart)[iii]
     Prop_raw = cumvar %>% select(all_of(PC_tmp)) %>%.[2,]
     Prop_std = Prop_raw/cumvar_sum
     varPart_weighed[iii,]=varPart[iii,]*Prop_std
   }
     colnames(varPart_weighed) = c("Constitutional","Extrarenal","Sex","HCQ","Hematological","MMF","Mucocutaneous","Musculoskeletal","Renal","Serological","TAC","PSL","Age","residual")
     varPart_weighed_sum=apply(varPart_weighed,2,sum) %>% as.data.frame() %>% t() %>% as.data.frame()
     varPart_weighed_sum$subset=subset_tmp
     if(kkk==1){varPart_weighed_sum2=varPart_weighed_sum}else{varPart_weighed_sum2=rbind(varPart_weighed_sum2,varPart_weighed_sum)}
  }
Organ_sum=varPart_weighed_sum2 %>% select("subset","Constitutional","Mucocutaneous","Musculoskeletal","Renal","Extrarenal","Hematological","Serological","PSL","HCQ","MMF","TAC","Sex","Age")
write.table_FT_2(Organ_sum,paste0("PCAres_Clinical/variancePartition/",today(),"_COI_PC30andClinical_varPart_Organ_auto.txt"))


#################################
# VarPart jack-knife resampling method
# Example: Total, one cell type
#################################
LIST="tmp_job_list/211128_COI_SLEuniqueandHC_27subsets_AfterCombat_PCres_list.txt"
#  task_id=1

list          = fread_FT(LIST)
subset_tmp    = list[task_id,3]
Rdata_path    = list[task_id,1]

paste0("subset_tmp : ",subset_tmp) %>% print()
paste0("Rdata_path : ",Rdata_path) %>% print()

out_f  = paste0("PCAres_Clinical/variancePartition/jackknife_total/")
out_ff = paste0(today(),"_COI_PC30andClinical_")

load(Rdata_path)
cumvar = summary(result)$importance[,1:7] %>% as.data.frame()
cumvar_sum = cumvar[3,7]

Matrix_full=as.data.frame(result$x) %>% select(1:7) %>% t() %>% as.data.frame()
Matrix_full=Matrix_full[,is.element(colnames(Matrix_full),name_list_SLE$name)]

for (mmm in 1:ncol(Matrix_full)){
     Matrix=Matrix_full[,-mmm] # jackknife -1

    Metadata=name_list%>%filter(name%in%colnames(Matrix))%>%select(-id)
    Metadata$name=factor(Metadata$name,levels=colnames(Matrix))
    Metadata=Metadata[order(Metadata$name),]
    Metadata$name=as.character(Metadata$name)

    all.equal(colnames(Matrix),Metadata$name) %>% print()
   
    form = ~ (1|Activity)+PSLmg+(1|Immunosuppressant)+age+(1|gender)
    varPart = fitExtractVarPartModel( Matrix, form, Metadata )
    varPart_weighed=varPart

  for (iii in 1:7){
     PC_tmp   = rownames(varPart)[iii]
     Prop_raw = cumvar %>% select(all_of(PC_tmp)) %>%.[2,]
     Prop_std = Prop_raw/cumvar_sum
     varPart_weighed[iii,]=varPart[iii,]*Prop_std
   }
     colnames(varPart_weighed) = c("Activity","Sex","Immunosuppressant","PSL","Age","residual")
     varPart_weighed_sum=apply(varPart_weighed,2,sum) 
     varPart_weighed_sum=data.frame(varPart_weighed_sum)%>% t()
     rownames(varPart_weighed_sum) = paste0("test",mmm)

     if(mmm==1){varPart_weighed_jack=varPart_weighed_sum}else{varPart_weighed_jack=rbind(varPart_weighed_jack,varPart_weighed_sum)}
}

dir.create_p(paste0(out_f,"tmp_Rdata/"))
save(varPart_weighed_jack, file= paste0(out_f,"tmp_Rdata/",out_ff,subset_tmp,"_varPartjack_total.Rdata"))

write.table_n_2(data.frame(varPart_weighed_jack),"test",paste0(out_f,out_ff,subset_tmp,"_varPartjack_total.txt"))
###################################################################################################

#################################
# VarPart jack-knife res summary
#################################

# Total activity
list = make_list("PCAres_Clinical/variancePartition/jackknife_total","_varPartjack_total.txt")
list$subset = take_factor(list$FILE,5:6,"_")
Total_sum_ordered=fread_FT("PCAres_Clinical/variancePartition/211222_COI_PC30andClinical_varPart_Total_auto.txt")
median =median(Total_sum_ordered$Activity)

for (kkk in 1:nrow(list)){
    subset_tmp= list$subset[kkk]
    perc_tmp  = fread_FT(list$PATH[kkk]) %>% select(test,Activity)
    percentile=quantile(perc_tmp$Activity,probs=c(0.025,0.975),na.rm=T) %>% as.data.frame() %>% rownames_to_column("test")
    colnames(percentile)[2]="Activity"
    perc_tmp=rbind(perc_tmp,percentile)
    perc_tmp$subset=subset_tmp
    perc_tmp=perc_tmp %>% select(subset,test,Activity)
    orig_tmp  = Total_sum_ordered %>% filter(subset==subset_tmp) %>% 
                                      mutate(test="orig") %>% select(subset,test,Activity)

    res_tmp=rbind(orig_tmp,perc_tmp) %>% pivot_longer(col=-c(subset,test),names_to="parameter",values_to="value") %>%
                                         pivot_wider(names_from="test",values_from="value") 

    jack_tmp=as.numeric(res_tmp[1,4:(ncol(res_tmp)-2)])
    res_tmp$test_num=(ncol(res_tmp)-5)
    res_tmp$smaller_num=length(which(jack_tmp<median))
    res_tmp$median=median
    res_tmp$Pjack=(res_tmp$smaller_num/res_tmp$test_num)
    res_tmp2=res_tmp%>%select("subset","parameter","orig","2.5%","97.5%","median","test_num","smaller_num","Pjack")
    colnames(res_tmp2)[4:5]=c("lower","upper")

    if(kkk==1){res_sum=res_tmp2}else{res_sum=rbind(res_sum,res_tmp2)}
}

res_sum$BonfP = p.adjust(res_sum$Pjack,method=("bonferroni"))
write.table_FT_2(res_sum,paste0("PCAres_Clinical/variancePartition/",today(),"_COI_PC30andClinical_varPartactivity_jackres.txt"))

res_sum_act=res_sum %>% left_join(.,celltype_corresp[,1:3],by="subset")
res_sum_act$subset    = factor(res_sum_act$subset,levels=rev(celltype_reordered_27))
res_sum_act$lineage   = factor(res_sum_act$lineage,levels=c("CD4","CD8","NK","B","Monocyte","DC","Neutrophil"))
res_sum_act$sig=ifelse(res_sum_act$BonfP<0.05,"*","")

p = ggplot(data=res_sum_act, aes(x=subset,y=orig,ymin=lower,ymax=upper,fill=lineage))+
  geom_bar(stat="identity")+
  geom_hline(yintercept=median,col="black",linetype="longdash")+
  geom_errorbar(width=0.2)+
  geom_text(aes(label=sig,y=0.2),size=7,vjust=0.75)+
  theme_classic()+
  coord_flip()+
  #facet_wrap(~ parameter,ncol = 5)+
  scale_fill_manual(values=col2)+
  theme(axis.text.x=element_text(colour="black",size=16),
         axis.text.y=element_text(colour="black",size=16),
         axis.title.x=element_text(colour="black",size=16),
         axis.title.y=element_blank(),
         plot.title=element_blank(),
         strip.text=element_text(colour="black",size=16),
         legend.position="none")+
  scale_x_discrete(labels= label)+
  scale_y_continuous(breaks= c(0,0.1,0.2))+
  labs(y="Variance explained")

pdf_3(paste0("Figure/",today(),"_COI_PC30andClinical_varPartactivity_jackres.pdf"),h=5.5,w=4)
 plot(p)
dev.off()


# Organ Main fig
list = make_list("PCAres_Clinical/variancePartition/jackknife_organ","_varPartjack_organ.txt")
list$subset = take_factor(list$FILE,5:6,"_")
Organ_sum_ordered = fread_FT("PCAres_Clinical/variancePartition/211222_COI_PC30andClinical_varPart_Organ_auto.txt")

for(iii in 2:4){
    parameter_tmp = colnames(Organ_sum_ordered)[iii+1]
    median = median(Organ_sum_ordered[,(iii+1)])

  for (kkk in 1:nrow(list)){
    subset_tmp= list$subset[kkk]
    perc_tmp  = fread_FT(list$PATH[kkk]) %>% select(test,all_of(parameter_tmp))
    percentile=quantile(perc_tmp[,2],probs=c(0.025,0.975),na.rm=T) %>% as.data.frame() %>% rownames_to_column("test")
    colnames(percentile)[2]=parameter_tmp
    perc_tmp=rbind(perc_tmp,percentile)
    perc_tmp$subset=subset_tmp
    perc_tmp=perc_tmp %>% select(subset,test,all_of(parameter_tmp))
    orig_tmp  = Organ_sum_ordered %>% filter(subset==subset_tmp) %>% 
                                      mutate(test="orig") %>% select(subset,test,all_of(parameter_tmp))

    res_tmp=rbind(orig_tmp,perc_tmp) %>% pivot_longer(col=-c(subset,test),names_to="parameter",values_to="value") %>%
                                         pivot_wider(names_from="test",values_from="value") 

    jack_tmp=as.numeric(res_tmp[1,4:(ncol(res_tmp)-2)])
    res_tmp$test_num=(ncol(res_tmp)-5)
    res_tmp$smaller_num=length(which(jack_tmp<median))
    res_tmp$median=median
    res_tmp$Pjack=(res_tmp$smaller_num/res_tmp$test_num)
    res_tmp2=res_tmp%>%select("subset","parameter","orig","2.5%","97.5%","median","test_num","smaller_num","Pjack")
    colnames(res_tmp2)[4:5]=c("lower","upper")

    if(kkk==1){res_sum=res_tmp2}else{res_sum=rbind(res_sum,res_tmp2)}
  }

res_sum$BonfP = p.adjust(res_sum$Pjack,method=("bonferroni"))
write.table_FT_2(res_sum,paste0("PCAres_Clinical/variancePartition/",today(),"_COI_PC30andClinical_varPart",parameter_tmp,"_jackres.txt"))

res_sum_tmp=res_sum %>% left_join(.,celltype_corresp[,1:3],by="subset")
res_sum_tmp$subset    = factor(res_sum_tmp$subset,levels=rev(celltype_reordered_27))
res_sum_tmp$lineage   = factor(res_sum_tmp$lineage,levels=c("CD4","CD8","NK","B","Monocyte","DC","Neutrophil"))
res_sum_tmp$sig=ifelse(res_sum_tmp$BonfP<0.05,"*","")

p = ggplot(data=res_sum_tmp, aes(x=subset,y=orig,ymin=lower,ymax=upper,fill=lineage))+
  geom_bar(stat="identity")+
  geom_hline(yintercept=median,col="black",linetype="longdash")+
  geom_errorbar(width=0.2)+
  geom_text(aes(label=sig,y=0.1),size=7,vjust=0.75)+
  theme_classic()+
  coord_flip()+
  #facet_wrap(~ parameter,ncol = 5)+
  scale_fill_manual(values=col2)+
  theme(axis.text.x=element_text(colour="black",size=16),
         axis.text.y=element_text(colour="black",size=16),
         axis.title.x=element_text(colour="black",size=16),
         axis.title.y=element_blank(),
         plot.title=element_blank(),
         strip.text=element_text(colour="black",size=16),
         legend.position="none")+
  scale_x_discrete(labels= label)+
  scale_y_continuous(limits=c(0,0.102),breaks= c(0,0.05,0.1))+
  labs(y="Variance explained")

  pdf_3(paste0("Figure/",today(),"_COI_PC30andClinical_varPart",parameter_tmp,"_jackres.pdf"),h=5.5,w=4)
   plot(p)
  dev.off()
}

# Organ Extended fig
for(iii in c(1,5:7)){
    parameter_tmp = colnames(Organ_sum_ordered)[iii+1]
    median = median(Organ_sum_ordered[,(iii+1)])

  for (kkk in 1:nrow(list)){
    subset_tmp= list$subset[kkk]
    perc_tmp  = fread_FT(list$PATH[kkk]) %>% select(test,all_of(parameter_tmp))
    percentile=quantile(perc_tmp[,2],probs=c(0.025,0.975),na.rm=T) %>% as.data.frame() %>% rownames_to_column("test")
    colnames(percentile)[2]=parameter_tmp
    perc_tmp=rbind(perc_tmp,percentile)
    perc_tmp$subset=subset_tmp
    perc_tmp=perc_tmp %>% select(subset,test,all_of(parameter_tmp))
    orig_tmp  = Organ_sum_ordered %>% filter(subset==subset_tmp) %>% 
                                      mutate(test="orig") %>% select(subset,test,all_of(parameter_tmp))

    res_tmp=rbind(orig_tmp,perc_tmp) %>% pivot_longer(col=-c(subset,test),names_to="parameter",values_to="value") %>%
                                         pivot_wider(names_from="test",values_from="value") 

    jack_tmp=as.numeric(res_tmp[1,4:(ncol(res_tmp)-2)])
    res_tmp$test_num=(ncol(res_tmp)-5)
    res_tmp$smaller_num=length(which(jack_tmp<median))
    res_tmp$median=median
    res_tmp$Pjack=(res_tmp$smaller_num/res_tmp$test_num)
    res_tmp2=res_tmp%>%select("subset","parameter","orig","2.5%","97.5%","median","test_num","smaller_num","Pjack")
    colnames(res_tmp2)[4:5]=c("lower","upper")

    if(kkk==1){res_sum=res_tmp2}else{res_sum=rbind(res_sum,res_tmp2)}
  }

res_sum$BonfP = p.adjust(res_sum$Pjack,method=("bonferroni"))
write.table_FT_2(res_sum,paste0("PCAres_Clinical/variancePartition/",today(),"_COI_PC30andClinical_varPart",parameter_tmp,"_jackres.txt"))

res_sum_tmp=res_sum %>% left_join(.,celltype_corresp[,1:3],by="subset")
res_sum_tmp$subset    = factor(res_sum_tmp$subset,levels=rev(celltype_reordered_27))
res_sum_tmp$lineage   = factor(res_sum_tmp$lineage,levels=c("CD4","CD8","NK","B","Monocyte","DC","Neutrophil"))
res_sum_tmp$sig=ifelse(res_sum_tmp$BonfP<0.05,"*","")

p = ggplot(data=res_sum_tmp, aes(x=subset,y=orig,ymin=lower,ymax=upper,fill=lineage))+
  geom_bar(stat="identity")+
  geom_hline(yintercept=median,col="black",linetype="longdash")+
  geom_errorbar(width=0.2)+
  geom_text(aes(label=sig,y=0.1),size=7,vjust=0.75)+
  theme_classic()+
  coord_flip()+
  #facet_wrap(~ parameter,ncol = 5)+
  scale_fill_manual(values=col2)+
  theme(axis.text.x=element_text(colour="black",size=16),
         axis.text.y=element_text(colour="black",size=16),
         axis.title.x=element_text(colour="black",size=16),
         axis.title.y=element_blank(),
         plot.title=element_blank(),
         strip.text=element_text(colour="black",size=16),
         legend.position="none")+
  scale_x_discrete(labels= label)+
  scale_y_continuous(limits=c(0,0.102),breaks= c(0,0.05,0.1))+
  labs(y="Variance explained")

  pdf_3(paste0("Figure/",today(),"_COI_PC30andClinical_varPart",parameter_tmp,"_jackres.pdf"),h=5.5,w=4)
   plot(p)
  dev.off()
}

# Tx
for(iii in 8:11){
    parameter_tmp = colnames(Organ_sum_ordered)[iii+1]
    median = median(Organ_sum_ordered[,(iii+1)])

  for (kkk in 1:nrow(list)){
    subset_tmp= list$subset[kkk]
    perc_tmp  = fread_FT(list$PATH[kkk]) %>% select(test,all_of(parameter_tmp))
    percentile=quantile(perc_tmp[,2],probs=c(0.025,0.975),na.rm=T) %>% as.data.frame() %>% rownames_to_column("test")
    colnames(percentile)[2]=parameter_tmp
    perc_tmp=rbind(perc_tmp,percentile)
    perc_tmp$subset=subset_tmp
    perc_tmp=perc_tmp %>% select(subset,test,all_of(parameter_tmp))
    orig_tmp  = Organ_sum_ordered %>% filter(subset==subset_tmp) %>% 
                                      mutate(test="orig") %>% select(subset,test,all_of(parameter_tmp))

    res_tmp=rbind(orig_tmp,perc_tmp) %>% pivot_longer(col=-c(subset,test),names_to="parameter",values_to="value") %>%
                                         pivot_wider(names_from="test",values_from="value") 

    jack_tmp=as.numeric(res_tmp[1,4:(ncol(res_tmp)-2)])
    res_tmp$test_num=(ncol(res_tmp)-5)
    res_tmp$smaller_num=length(which(jack_tmp<median))
    res_tmp$median=median
    res_tmp$Pjack=(res_tmp$smaller_num/res_tmp$test_num)
    res_tmp2=res_tmp%>%select("subset","parameter","orig","2.5%","97.5%","median","test_num","smaller_num","Pjack")
    colnames(res_tmp2)[4:5]=c("lower","upper")

    if(kkk==1){res_sum=res_tmp2}else{res_sum=rbind(res_sum,res_tmp2)}
  }

res_sum$BonfP = p.adjust(res_sum$Pjack,method=("bonferroni"))
write.table_FT_2(res_sum,paste0("PCAres_Clinical/variancePartition/",today(),"_COI_PC30andClinical_varPart",parameter_tmp,"_jackres.txt"))

res_sum_tmp=res_sum %>% left_join(.,celltype_corresp[,1:3],by="subset")
res_sum_tmp$subset    = factor(res_sum_tmp$subset,levels=rev(celltype_reordered_27))
res_sum_tmp$lineage   = factor(res_sum_tmp$lineage,levels=c("CD4","CD8","NK","B","Monocyte","DC","Neutrophil"))
res_sum_tmp$sig=ifelse(res_sum_tmp$BonfP<0.05,"*","")

p = ggplot(data=res_sum_tmp, aes(x=subset,y=orig,ymin=lower,ymax=upper,fill=lineage))+
  geom_bar(stat="identity")+
  geom_hline(yintercept=median,col="black",linetype="longdash")+
  geom_errorbar(width=0.2)+
  geom_text(aes(label=sig,y=0.2),size=7,vjust=0.75)+
  theme_classic()+
  coord_flip()+
  #facet_wrap(~ parameter,ncol = 5)+
  scale_fill_manual(values=col2)+
  theme(axis.text.x=element_text(colour="black",size=16),
         axis.text.y=element_text(colour="black",size=16),
         axis.title.x=element_text(colour="black",size=16),
         axis.title.y=element_blank(),
         plot.title=element_blank(),
         strip.text=element_text(colour="black",size=16),
         legend.position="none")+
  scale_x_discrete(labels= label)+
  scale_y_continuous(limits=c(0,0.2),breaks= c(0,0.1,0.2))+
  labs(y="Variance explained")

  pdf_3(paste0("Figure/",today(),"_COI_PC30andClinical_varPart",parameter_tmp,"_jackres.pdf"),h=5.5,w=4)
   plot(p)
  dev.off()
}

######### Supple Table
list=make_list("PCAres_Clinical/variancePartition","_jackres.txt")

for(iii in 1:nrow(list)){
    res_tmp=fread_FT(list$PATH[iii])%>%select(2,1,3:5,9:10)
    if(iii==1){res_sum=res_tmp}else{res_sum=rbind(res_sum,res_tmp)}
}
res_sum$subset=factor(res_sum$subset,levels=celltype_reordered_27)
res_sum$parameter=factor(res_sum$parameter,levels=c("Activity",
                                                    "Constitutional","Mucocutaneous","Musculoskeletal","Renal","Extrarenal","Hematological","Serological",
                                                    "MMF","TAC","HCQ","PSL"))
levels(res_sum$parameter)[c(1,6)]=c("Overall disease activity","Extrarenal severe")

res_sum=res_sum[order(res_sum$subset),]
res_sum=res_sum[order(res_sum$parameter),]
res_sum=res_sum%>%left_join(.,celltype_corresp[,c(1,2)],by="subset")%>%select(parameter,label,orig,lower,upper,Pjack,BonfP)

res_sum$orig = formatC(res_sum$orig,digits=2)
res_sum$lower   = formatC(res_sum$lower,digits=2)
res_sum$upper    = formatC(res_sum$upper,digits=2)
res_sum$Pjack    = formatC(res_sum$Pjack,digits=2)
res_sum$BonfP    = formatC(res_sum$BonfP,digits=2)

colnames(res_sum)=c("Parameter","Cell type","Variance explained","Lower 95percentile","Upper 95percentile","Pjk","Bonferonni-corrected Pjk")

write.table_FT_2(res_sum,paste0(today(),"_COI_varPart_Alljackres_forSupple.txt"))

##########################################################################################################################
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
[1] parallel  stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] variancePartition_1.20.0 Biobase_2.50.0           BiocGenerics_0.36.1     
 [4] scales_1.1.1             BiocParallel_1.24.1      limma_3.46.0            
 [7] lme4_1.1-27.1            Matrix_1.3-4             ggsci_2.9               
[10] RColorBrewer_1.1-2       data.table_1.14.0        forcats_0.5.1           
[13] stringr_1.4.0            dplyr_1.0.7              purrr_0.3.4             
[16] readr_1.4.0              tidyr_1.1.3              tibble_3.1.2            
[19] ggplot2_3.3.5            tidyverse_1.3.1         

loaded via a namespace (and not attached):
 [1] httr_1.4.2         jsonlite_1.7.2     splines_4.0.2      foreach_1.5.1     
 [5] modelr_0.1.8       gtools_3.9.2       assertthat_0.2.1   cellranger_1.1.0  
 [9] progress_1.2.2     pillar_1.6.1       backports_1.2.1    lattice_0.20-41   
[13] glue_1.4.2         digest_0.6.27      rvest_1.0.0        minqa_1.2.4       
[17] colorspace_2.0-2   plyr_1.8.6         pkgconfig_2.0.3    broom_0.7.8       
[21] haven_2.4.1        farver_2.1.0       generics_0.1.0     ellipsis_0.3.2    
[25] withr_2.4.2        pbkrtest_0.5.1     cli_2.5.0          magrittr_2.0.1    
[29] crayon_1.4.1       readxl_1.3.1       ps_1.6.0           fs_1.5.0          
[33] fansi_0.5.0        doParallel_1.0.16  nlme_3.1-152       MASS_7.3-54       
[37] gplots_3.1.1       xml2_1.3.2         prettyunits_1.1.1  tools_4.0.2       
[41] hms_1.1.0          lifecycle_1.0.0    munsell_0.5.0      reprex_2.0.0      
[45] colorRamps_2.3     compiler_4.0.2     caTools_1.18.2     rlang_0.4.11      
[49] grid_4.0.2         nloptr_1.2.2.2     iterators_1.0.13   rstudioapi_0.13   
[53] bitops_1.0-7       boot_1.3-28        gtable_0.3.0       codetools_0.2-18  
[57] DBI_1.1.1          reshape2_1.4.4     R6_2.5.0           lubridate_1.7.10  
[61] utf8_1.2.1         KernSmooth_2.23-20 stringi_1.6.2      Rcpp_1.0.6        
[65] vctrs_0.3.8        dbplyr_2.1.1       tidyselect_1.1.1 












