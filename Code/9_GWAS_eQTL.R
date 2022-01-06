#---------------------------------------#
# 211208 Nakano
# SLE GWAS eQTL DEGs
#---------------------------------------#

source("my.source.R")
options(stringsAsFactors=F)
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ggsci))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(coin))
suppressPackageStartupMessages(library(exact2x2))

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

######################################################################################################################
# LDSC Meta-analysis
######################################################################################################################
#EAS 
#EAS
d1.deg = fread_FT("LDSC_res/Lupus_ard2021.eas.deg.dist100000.ldsc_seg") %>%
         filter(PC=="state"|PC=="activity")
d1.seg = fread_FT("LDSC_res/Lupus_ard2021.eas.seg.dist100000.ldsc_seg")
d1 = rbind(d1.deg, d1.seg)
d1$tag = paste0(d1$Cell,":",d1$PC)
d1 = d1[,c("tag","Cell","PC","Coefficient","Coefficient_std_error")]

#EUR
d2.deg = fread_FT("LDSC_res/Lupus_langefeld.eur.deg.dist100000.ldsc_seg") %>%
         filter(PC=="state"|PC=="activity")
d2.seg = fread_FT("LDSC_res/Lupus_langefeld.eur.seg.dist100000.ldsc_seg")
d2 = rbind(d2.deg, d2.seg)
d2$tag = paste0(d2$Cell,":",d2$PC)
d2 = d2[,c("tag","Coefficient","Coefficient_std_error")]

df = left_join(d1,d2,by="tag")

#adjust by h2
xh2 = 0.0814/5633366 #EAS
yh2 = 0.4142/5926598 #EUR

df$Coefficient.x = df$Coefficient.x / xh2
df$Coefficient_std_error.x = df$Coefficient_std_error.x/ xh2
df$Coefficient.y = df$Coefficient.y / yh2
df$Coefficient_std_error.y = df$Coefficient_std_error.y/ yh2

df$metabeta = apply(df[,4:7],1,function(x){
    x = as.numeric(x);
    beta = x[c(1,3)];
    se = x[c(2,4)];
    we = 1 / ( se )^2;
    ntest = length( beta );
    
    metabeta = sum( beta * we ) / sum( we );
    metase = sqrt( 1/ sum( we ) );
    metaz  = metabeta / metase;
    metap = pnorm(- metaz)
    return(metabeta)})

df$metase = apply(df[,4:7],1,function(x){
    x = as.numeric(x);
    beta = x[c(1,3)];
    se = x[c(2,4)];
    we = 1 / ( se )^2;
    ntest = length( beta );
    
    metabeta = sum( beta * we ) / sum( we );
    metase = sqrt( 1/ sum( we ) );
    metaz  = metabeta / metase;
    metap = pnorm(- metaz)
    return(metase)})

df$metap = apply(df[,4:7],1,function(x){
    x = as.numeric(x);
    beta = x[c(1,3)];
    se = x[c(2,4)];
    we = 1 / ( se )^2;
    ntest = length( beta );
    
    metabeta = sum( beta * we ) / sum( we );
    metase = sqrt( 1/ sum( we ) );
    metaz  = metabeta / metase;
    metap = pnorm(- metaz)
    return(metap)})

df$phet = apply(df[,4:7],1,function(x){
    x = as.numeric(x);
    beta = x[c(1,3)];
    se = x[c(2,4)];
    we = 1 / ( se )^2;
    ntest = length( beta );
    
    metabeta = sum( beta * we ) / sum( we );
    Q = sum( we * ( beta - metabeta )^2 )
    Qpvalue = pchisq( Q, ntest - 1,lower.tail=FALSE)
    return(Qpvalue)})


df$FDR = p.adjust(df$metap,method="BH")
df = df[order(df$metap),]

df2 = left_join(df,celltype_corresp%>%select(subset,label,lineage),by=c("Cell"="subset"))
df2$Cell    = factor(df2$Cell,levels=celltype_reordered_27)
df2$lineage = factor(df2$lineage,levels=c("CD4","CD8","NK","B","Monocyte","DC","Neutrophil"))
df2$PC = factor(df2$PC,levels=c("seg","state","activity"))
levels(df2$PC) = c("HC-SEG","Disease-state","Disease-activity")

####### Compare 2 cohorts
filter(df2,metap < 0.05/nrow(df2)) %>% nrow()
filter(df2,metap < 0.05/nrow(df2))[,"tag"] %>% sort()
filter(df2,metap < 0.05/nrow(df2))[,"phet"] %>% summary()

df2$is_sig = ifelse(df2$metap<(0.05/nrow(df2)),"Sig","Nonsig")
df2$sig_type=ifelse(df2$is_sig=="Sig",as.character(df2$PC),"Not significant")
df2$siglabel=ifelse(df2$is_sig=="Sig",paste0(df2$PC,"_",df2$label),"")

write.table_FT_2(df2,paste0("LDSC_res/DrIshigaki_2/",today(),"_COI_state_activity_HCSEG_LDSCmeta_highZver_logP_res.txt"))

df2$siglabel=ifelse(df2$is_sig=="Sig",df2$label,"")
b_thresh = 0.05/nrow(df2)
df2$is_sig=factor(df2$is_sig,levels=c("Sig","Nonsig"))
df2$sig_type=factor(df2$sig_type,levels=c("Disease-state","Not significant"))

p2 = ggplot(data=df2,aes(x=Coefficient.x, y=Coefficient.y, color=sig_type)) + 
    geom_point() +
    geom_errorbarh(aes(xmin= Coefficient.x - 2*Coefficient_std_error.x, 
                       xmax= Coefficient.x + 2*Coefficient_std_error.x ), size=0.1) +
    geom_errorbar( aes(ymin= Coefficient.y - 2*Coefficient_std_error.y, 
                       ymax= Coefficient.y + 2*Coefficient_std_error.y), size=0.1) +
    geom_text_repel(aes(label=siglabel),max.overlaps=Inf,size=4.5)+
    theme_classic()+
    geom_abline(slope=1,intercept=0,size=0.2,col="black",linetype="dashed") +
    geom_vline(xintercept=0,size=0.2,linetype="longdash",col="black")  +
    geom_hline(yintercept=0,size=0.2,linetype="longdash",col="black")  +
    scale_color_manual(values=c("#0072B5FF","#cccccc")) +
    theme(axis.text.x=element_text(colour="black",size=16),
           axis.text.y=element_text(colour="black",size=16),
           axis.title.x=element_text(colour="black",size=16),
           axis.title.y=element_text(colour="black",size=16),
           plot.title=element_blank(),
           legend.position="right",
           legend.title=element_blank(),
           legend.text=element_text(colour="black",size=16))+
    guides(colour=guide_legend(override.aes=list(size=2)))+
    labs(x="Normalized coefficient (EAS)",y="Normalized coefficient (EUR)")

pdf_3(paste0("LDSC_res/DrIshigaki_2/",today(),"_COI_state_activity_HCSEG_LDSC_EAS_EUR_coef.pdf"),h=5,w=7)
 plot(p2)
dev.off()

cor.test(df2$Coefficient.x ,df2$Coefficient.y,method="spearman")$p.value
# [1] 0
cor(df2$Coefficient.x ,df2$Coefficient.y,method="spearman")
# [1] 0.6371048

cor.test(df2$Coefficient.x ,df2$Coefficient.y,method="pearson")$p.value
# [1] 1.467846e-12
cor(df2$Coefficient.x ,df2$Coefficient.y,method="pearson")
# [1] 0.6867203

####### Meta-res 3 gene set
# Barplot
df2$Cell    = factor(df2$Cell,levels=rev(celltype_reordered_27))
p = ggplot()+
     geom_bar(data=df2, aes(x=Cell,y=-log10(metap),fill=lineage),stat="identity")+
     geom_hline(yintercept=-log10(b_thresh),linetype="longdash",col="black")  +
     theme_classic()+
     coord_flip()+
     facet_wrap(~PC,ncol=3) +
     scale_fill_manual(values = col2)+
     theme(axis.text.x=element_text(colour="black",size=16),
           axis.text.y=element_text(colour="black",size=16),
           axis.title.x=element_text(colour="black",size=16),
           axis.title.y=element_blank(),
           strip.background=element_blank(),
           strip.text=element_text(colour="black",size=16),
           panel.background =element_rect(fill="transparent",colour="black",size=1),
           panel.grid = element_blank(),
           legend.position="right",
           legend.title=element_text(colour="black",size=16),
           legend.text=element_text(colour="black",size=16))+
     scale_x_discrete(labels= label)+
　　　labs(y="-log10(P)")

pdf_3(paste0("LDSC_res/DrIshigaki_2/",today(),"_COI_state_activity_SEG_LDSCmeta_highZver_logP_bar.pdf"),h=6,w=11)
 plot(p)
dev.off()

df3=df2%>%select(Cell,label,lineage,PC,metap)%>%pivot_wider(names_from=PC,values_from=metap) 
colnames(df3)[4:5]=c("state","activity")

df3$check = ifelse(df3$state<df3$activity,"state","activity")
table(df3$check)
# activity    state 
#        2       25 

wilcox.test(df3$state, df3$activity, paired = TRUE)
# p-value = 5.528e-06

##### Supple
df4=df2%>%mutate(Bonf_P=p.adjust(metap,method="bonferroni"))

df4$Cell=factor(df4$Cell,levels=celltype_reordered_27)
df4$PC=factor(df4$PC,levels=c("HC-SEG","Disease-state","Disease-activity"))
df4=df4[order(df4$Cell),]
df4=df4[order(df4$PC),]
df4=df4%>%select(PC,label,Coefficient.x,Coefficient_std_error.x,Coefficient.y,Coefficient_std_error.y,metabeta,metase,metap,Bonf_P)
df4$Coefficient.x            = formatC(df4$Coefficient.x,digits=2)
df4$Coefficient_std_error.x  = formatC(df4$Coefficient_std_error.x,digits=2)
df4$Coefficient.y            = formatC(df4$Coefficient.y,digits=2)
df4$Coefficient_std_error.y  = formatC(df4$Coefficient_std_error.y,digits=2)
df4$metabeta                 = formatC(df4$metabeta,digits=2)
df4$metase                   = formatC(df4$metase,digits=2)
df4$metap                    = formatC(df4$metap,digits=2)
df4$Bonf_P                   = formatC(df4$Bonf_P,digits=2)
colnames(df4)=c("Signature","Cell type","EAS coefficient","EAS SE","EUR coefficient","EUR SE","Meta effectsize","Meta SE","Meta Pvalue","Bonferonni-corrected P")

write.table_FT_2(df4,paste0(today(),"_LDSC_meta_Supple.txt"))

######################################################################################################################
# GWAS candidate gene enrichment 
######################################################################################################################
# GWAS catalog
# As of 2021-08-16, the GWAS Catalog contains 5273 publications and 276696 associations.
# GWAS Catalog data is currently mapped to Genome Assembly GRCh38.p13 and dbSNP Build 154.
# limited to p<5x10-8, genes were collapsed
# Gene symbols were aligned with UCSC38 list

whole_list = fread_FT("data/210911_SLE_GWAScatalog_cisgenes_new_UCSC38.txt")[,1]
res_sum = fread_FT("~/ws/2021/211129_COI_SLE_HC_run27_edgeR/res/stateactivity/211201_COI_27subsets_state_activity_logFC_logCPM_P_FDR.txt")

############
# Enrichment test of GWAS gene
res_sum_DEGs = res_sum %>% filter(state_FDR<0.05|activity_FDR<0.05)
backgound   = length(unique(res_sum$genes))

for(iii in 1:length(celltype_reordered_27)){
　　 subset_tmp = celltype_reordered_27[iii]
    res_tmp = res_sum %>% filter(subset==subset_tmp)
    GWAS_genes    = length(whole_list)   

    state_signum        = length(which(res_tmp$state_FDR<0.05))
    activity_signum     = length(which(res_tmp$activity_FDR<0.05))
    state_signum_GWAS   = length(which(res_tmp$state_FDR<0.05&res_tmp$genes%in%whole_list))
    activity_signum_GWAS= length(which(res_tmp$activity_FDR<0.05&res_tmp$genes%in%whole_list))
    
    mx1=matrix(c(state_signum_GWAS,GWAS_genes-state_signum_GWAS,state_signum-state_signum_GWAS,backgound_1-state_signum-(GWAS_genes-state_signum_GWAS)),nrow=2,byrow=T)
    fisher1=fisher.exact(mx1,alternative="greater")
    res1=data.frame(subset=subset_tmp,OR=fisher1$estimate,pvalue=fisher1$p.value,back="Allgenes",type="state")
    
    mx2=matrix(c(activity_signum_GWAS,GWAS_genes-activity_signum_GWAS,activity_signum-activity_signum_GWAS,backgound_1-activity_signum-(GWAS_genes-activity_signum_GWAS)),nrow=2,byrow=T)
    fisher2=fisher.exact(mx2,alternative="greater")
    res2=data.frame(subset=subset_tmp,OR=fisher3$estimate,pvalue=fisher3$p.value,back="Allgenes",type="activity")
    
    
    mx3=matrix(c(state_signum_GWAS,activity_signum_GWAS,state_signum-state_signum_GWAS,activity_signum-activity_signum_GWAS),nrow=2,byrow=T)
    fisher3=fisher.exact(mx3,alternative="greater")
    res3=data.frame(subset=subset_tmp,OR=fisher5$estimate,pvalue=fisher5$p.value,back="stateactivity",type="compare")

    fisher_tmp = bind_rows(list(res1,res2,res3))

    if(iii==1){fisher_sum=fisher_tmp}else{fisher_sum=rbind(fisher_sum,fisher_tmp)}
  }

fisher_sum = left_join(fisher_sum,celltype_corresp%>%select(subset,label,lineage),by="subset")
  
fisher_allgenes = fisher_sum %>% filter(back=="Allgenes")
fisher_allgenes$FDR = p.adjust(fisher_allgenes$pvalue,method="BH")
fisher_allgenes$Bonf = p.adjust(fisher_allgenes$pvalue,method="bonferroni")
write.table_FT_2(fisher_allgenes,paste0("GWAScatalog/",today(),"_COI_27subsets_state_activity_GWAS_fisheronesided_Allgenes.txt"))
 
fisher_compare = fisher_sum %>% filter(back=="stateactivity")
fisher_compare$FDR = p.adjust(fisher_compare$pvalue,method="BH")
fisher_compare$Bonf = p.adjust(fisher_compare$pvalue,method="bonferroni")
write.table_FT_2(fisher_compare,paste0("GWAScatalog/",today(),"_COI_27subsets_state_activity_GWAS_fisheronesided_directcompare.txt"))

fisher_sum$subset = factor(fisher_sum$subset,levels=celltype_reordered_27)
fisher_sum$lineage = factor(fisher_sum$lineage,levels=c("CD4","CD8","NK","B","Monocyte","DC","Neutrophil"))

# Barplot
fisher_allgenes$subset = factor(fisher_allgenes$subset,levels=celltype_reordered_27)
fisher_allgenes$lineage = factor(fisher_allgenes$lineage,levels=c("CD4","CD8","NK","B","Monocyte","DC","Neutrophil"))
fisher_allgenes$type = factor(fisher_allgenes$type,levels=c("state","activity"))
levels(fisher_allgenes$type) = c("Disease-state","Disease-activity")
b_thresh = 0.05/nrow(fisher_allgenes)

  p = ggplot()+
     geom_bar(data=fisher_allgenes, aes(x=subset,y=-log10(pvalue),fill=lineage),stat="identity")+
     geom_hline(yintercept=-log10(b_thresh),linetype="longdash",col="black")  +
     theme_classic()+
     scale_fill_manual(values = col2)+
     facet_wrap(~type,ncol=1) + 
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
　　　labs(y="-log10(Pvalue)")

pdf_3(paste0("GWAScatalog/",today(),"_COI_27subsets_state_activity_GWAS_fisheronesided_Allgenes_bar.pdf"),h=6,w=10)
 plot(p)
dev.off()

fisher_compare$subset = factor(fisher_compare$subset,levels=celltype_reordered_27)
fisher_compare$lineage = factor(fisher_compare$lineage,levels=c("CD4","CD8","NK","B","Monocyte","DC","Neutrophil"))

   p = ggplot()+
        geom_bar(data=fisher_compare, aes(x=subset,y=OR,fill=lineage),stat="identity")+
        geom_hline(yintercept=1,linetype="longdash",col="black")+
        theme_classic()+
        scale_fill_manual(values = col2)+
        #facet_wrap(~type, labeller=facet,ncol = 1) + 
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
   　　　labs(y="OR (Disease-state vs activity)")
   
pdf_3(paste0("GWAScatalog/",today(),"_COI_27subsets_state_activity_GWAS_fisheronesided_directcompare_bar.pdf"),h=4.5,w=10)
 plot(p)
dev.off()

fisher_allgenes2=fisher_allgenes%>%select(subset,label,lineage,type,pvalue)%>%pivot_wider(names_from=type,values_from=pvalue) 
colnames(fisher_allgenes2)[4:5]=c("state","activity")

wilcox.test(fisher_allgenes2$state, fisher_allgenes2$activity, paired = TRUE)
# p-value = 0.002469

fisher_allgenes3=fisher_allgenes2%>%
                 mutate(is_sig=ifelse(state<(0.05/54),"Sig","Nonsig"))
fisher_allgenes3$siglabel=ifelse(fisher_allgenes3$is_sig=="Sig",fisher_allgenes3$label,"")
fisher_allgenes3$subset = factor(fisher_allgenes3$subset,levels=celltype_reordered_27)
fisher_allgenes3$lineage= factor(fisher_allgenes3$lineage,levels=c("CD4","CD8","NK","B","Monocyte","DC","Neutrophil"))
b_thresh = 0.05/54

p = ggplot(data=fisher_allgenes3, aes(x=-log10(state),y=-log10(activity)))+
     geom_point(aes(color=lineage),stat="identity")+
     geom_text_repel(aes(label=siglabel),size=4.5)+
     geom_abline(slope=1,intercept=0,size=0.2,col="black",linetype="dashed") +
     geom_vline(xintercept=-log10(b_thresh),size=0.2,linetype="longdash",col="black")  +
     geom_hline(yintercept=-log10(b_thresh),size=0.2,linetype="longdash",col="black")  +
     theme_classic()+
     scale_color_manual(values = col2)+
     theme(axis.text.x=element_text(colour="black",size=16),
           axis.text.y=element_text(colour="black",size=16),
           axis.title.x=element_text(colour="black",size=16),
           axis.title.y=element_text(colour="black",size=16),
           plot.title=element_blank(),
           strip.text=element_text(colour="black",size=16),
           legend.position="right",
           legend.title=element_text(colour="black",size=16),
           legend.text=element_text(colour="black",size=16))+
     scale_x_continuous(breaks=c(0,1,2,3,4),limits=c(0,4.5))+
     scale_y_continuous(breaks=c(0,1,2,3,4),limits=c(0,4.5))+
     guides(colour=guide_legend(override.aes=list(size=2)))+
　　　labs(x="Disease-state  -log10(P)", y="Disease-activity  -log10(P)")

pdf_3(paste0("GWAScatalog/",today(),"_COI_state_activity_GWAS_fisheronesided_logP_point_lim.pdf"),h=5,w=6.5)
 plot(p)
dev.off()

########Supple
fisher_allgenes4=fisher_allgenes
fisher_allgenes4$subset=factor(fisher_allgenes4$subset,levels=celltype_reordered_27)
fisher_allgenes4$type=factor(fisher_allgenes4$type,levels=c("state","activity"))
levels(fisher_allgenes4$type)=c("Disease-state","Disease-activity")
fisher_allgenes4=fisher_allgenes4[order(fisher_allgenes4$subset),]
fisher_allgenes4=fisher_allgenes4[order(fisher_allgenes4$type),]
fisher_allgenes4=fisher_allgenes4%>%select(type,label,OR,pvalue,Bonf)

fisher_allgenes4$OR      = formatC(fisher_allgenes4$OR,digits=2)
fisher_allgenes4$pvalue  = formatC(fisher_allgenes4$pvalue,digits=2)
fisher_allgenes4$Bonf    = formatC(fisher_allgenes4$Bonf,digits=2)
colnames(fisher_allgenes4)=c("Signature","Cell type","Odds ratio","Pvalue","Bonferonni-corrected P")
write.table_FT_2(fisher_allgenes4,paste0(today(),"_GWAS_enrich_Supple.txt"))


######################################################################################################################
# Coherent incoherent genes
######################################################################################################################
# stringent colocalization list
Ota_coloc_list = read.table_FT("data/colocalization_for_heatmap_200310.txt") %>% filter(!SUBSET=="MCD8")
subset_corresp2 = read.table_FT("data_ref/210604_COI_subset_name_w_Abb_wo_b_Aria_only_new_DN12combine.txt") %>% dplyr::select(subset_abb,subset_name)
Ota_coloc_list2 = Ota_coloc_list %>% 
                  left_join(.,subset_corresp2,by=c("SUBSET"="subset_name")) %>%
                  mutate(subset=paste0(subset_abb,"_",SUBSET)) %>% dplyr::select(-c(subset_abb,SUBSET)) %>%
                  select(subset,GENE,GWAS)
colnames(Ota_coloc_list2) = c("subset","GENE","rs")
Ota_coloc_list2$rs = gsub("rs36059542","3_119435676_TA_T;rs36059542",Ota_coloc_list2$rs)
Ota_coloc_list2$rs = gsub("rs35552623","rs76724843;rs35552623",Ota_coloc_list2$rs)

write.table_FT_2(Ota_coloc_list2,paste0("data/",today(),"_colocalizedSNP_stringent.txt"))
# Ota_coloc_list2=fread_FT("data/211208_colocalizedSNP_stringent.txt")
##############
# Dr. Ota Figure order
Coloc_gene     =　c("ARHGAP31","IRF8","NCF1","TCF7","PTGER4","PAQR3","CDH1","GPX3","AHNAK2","SLC15A4","LBH","IRF5","RASGRP1","PTPRC","LRRC25","KAT8","UBE2L3","APOLD1","BLK","ELF1")

Ota_coloc_list2$GENE   = factor(Ota_coloc_list2$GENE,levels=Coloc_gene)
Ota_coloc_list2$subset = factor(Ota_coloc_list2$subset,levels=celltype_reordered_27) 
Ota_coloc_list2 = Ota_coloc_list2[order(Ota_coloc_list2$subset),]       
Ota_coloc_list2 = Ota_coloc_list2[order(Ota_coloc_list2$GENE),]
Ota_coloc_list2$subset_gene=paste0(Ota_coloc_list2$subset,"_",Ota_coloc_list2$GENE)

Ota_coloc_list2$direction = c(-1,-1,-1,-1,1,rep(1,2),rep(-1,2),rep(-1,2),rep(1,3),rep(1,4),c(1,rep(-1,8)),rep(1,9),c(rep(-1,7),rep(1,3)),c(rep(-1,8),rep(1,3)),rep(1,12),rep(1,14),rep(1,15),rep(1,18),rep(-1,19),rep(-1,25))

res_sum = fread_FT("~/ws/2021/211129_COI_SLE_HC_run27_edgeR/res/stateactivity/211201_COI_27subsets_state_activity_logFC_logCPM_P_FDR.txt")
res_sum2 = res_sum %>% mutate(subset_gene=paste0(subset,"_",genes)) %>% dplyr::select(subset_gene,state_logFC,state_PValue,state_FDR,activity_logFC,activity_PValue,activity_FDR)

res_sum3 = left_join(Ota_coloc_list2,res_sum2,by="subset_gene")
res_sum3$state_logFC_2 = (res_sum3$direction)*(res_sum3$state_logFC)
res_sum3$activity_logFC_2 = (res_sum3$direction)*(res_sum3$activity_logFC)

res_sum3_state = res_sum3 %>% filter(state_FDR<0.05) %>% dplyr::select(GENE,subset,state_logFC_2)
colnames(res_sum3_state)[3] = "adjusted_logFC"
res_sum3_state$type="State"

res_sum3_act = res_sum3 %>% filter(activity_FDR<0.05) %>% dplyr::select(GENE,subset,activity_logFC_2)
colnames(res_sum3_act)[3] = "adjusted_logFC"
res_sum3_act$type="Activity"

data= rbind(res_sum3_state,res_sum3_act)
data$type = factor(data$type,levels=c("State","Activity"))
data$GENE = as.character(data$GENE)
data = data[order(data$GENE),]

write.table_FT_2(data,paste0("eQTL_coherent/",today(),"_colocalizedGene_State_Activity_onlyDEGs_logFC.txt"))

table(data$type)
#   State Activity 
#      42       40

data$GENE = factor(data$GENE,levels=unique(data$GENE))
data$subset = factor(data$subset,levels=celltype_reordered_27)

data_x = data %>% filter(adjusted_logFC>0)
table(data_x$type)
#   State Activity 
#      28       10 

unique(data$GENE)
#  [1] AHNAK2  APOLD1  BLK     ELF1    IRF5    IRF8    KAT8    LBH     LRRC25 
# [10] NCF1    PTPRC   RASGRP1 SLC15A4 UBE2L3 

col_full=c("#921b1b","#ff9994","#cc0010","#ff0000","#e6a88f","#e66557","#feaf53","#f77308","#fde0bd",
          "#995522","#cc9955","#eebb77","#bb9988",
          "#e0d51a",
          "#00561f","#228b22","#8fc31f","#33dd33","#c4f20b",
          "#1d2088","#1e90ff","#a4bdfc","#4d79f7",
          "#0bdaf2","#a4fbfa",
          "#cccccc","#999999")


colnames(data)[4]="DEGtype"
data$DEGtype = factor(data$DEGtype,levels=c("State","Activity"))
levels(data$DEGtype)=c("Disease-state","Disease-activity")

data_state = data %>% filter(DEGtype=="Disease-state")
data_act = data %>% filter(DEGtype=="Disease-activity")
median_state=median(data_state$adjusted_logFC)
median_act=median(data_act$adjusted_logFC)

p = ggplot(data=data,aes(x=adjusted_logFC,fill=DEGtype))+
         geom_histogram(position="identity",binwidth=0.2,boundary=0)+
         geom_vline(xintercept=0,linetype="longdash",col="black")+
         theme_classic()+
         facet_wrap(~DEGtype,ncol=1) + 
         scale_fill_manual(values=c("#0072B5FF","#BC3C29FF"))+
         theme(axis.text.x=element_text(colour="black",size=16),
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
         scale_x_continuous(breaks=c(-2,-1,0,1,2))+
         labs(x="Adjusted logFC",y="Number of eGenes")

pdf_3(paste0("eQTL_coherent/",today(),"_colocalizedGene_State_Activity_onlyDEGs_logFC_hist_facet.pdf"),h=2.5,w=8)
 plot(p)
dev.off()

# one sided binom test
binom.test(c(length(which(data_state$adjusted_logFC>0)),length(which(data_state$adjusted_logFC<0))),alternative="greater")
# p-value = 0.02178
binom.test(c(length(which(data_act$adjusted_logFC>0)),length(which(data_act$adjusted_logFC<0))),alternative="greater")
# p-value = 0.9997

############ Supple
Ota_coloc_list_unique = Ota_coloc_list[!duplicated(Ota_coloc_list$GENE),] 
data2=data%>%left_join(.,label_df,by="subset") %>% left_join(.,Ota_coloc_list_unique,by="GENE")
data2$subset=factor(data2$subset,levels=celltype_reordered_27)
data2$type=factor(data2$type,levels=c("State","Activity"))
levels(data2$type)=c("Disease-state","Disease-activity")
data2=data2[order(data2$subset),]
data2=data2[order(data2$type),]
data2=data2%>%select(eQTL,GENE,label,type,adjusted_logFC)

data2$adjusted_logFC      = formatC(data2$adjusted_logFC,digits=2)

colnames(data2)=c("GWAS-eQTL","Gene","Cell type","Signature","Adjusted logFC")
write.table_FT_2(data2,paste0(today(),"_eQTL_DEGs_Supple.txt"))

######################################################################################################################
# Context dependent eQTL
######################################################################################################################
list         = make_list("~/ws/2021/211124_COI_SLE_HC_run27_finaldata_Combat_PCA/data/SLEunique_andHC/each_subset/LowExp_TMM_logCPM_Combat","_logCPM_Combat.txt")
list$subset  = take_factor(list$FILE,5:6,"_")

write.table_FT_2(list,paste0("tmp_job_list/",today(),"_logCPM_list.txt"))

######################################################################################################################
# Example: One cell type
######################################################################################################################

LIST="tmp_job_list/211208_logCPM_list.txt"
##  task_id=1

list       = fread_FT(LIST)
subset_tmp = list[task_id,3]
count_path = list[task_id,1]

paste0("subset_tmp : ",subset_tmp) %>% print()
paste0("count_path  : ",  count_path )       %>% print()

res_sum = fread_FT("~/ws/2021/211129_COI_SLE_HC_run27_edgeR/res/stateactivity/211201_COI_27subsets_state_activity_logFC_logCPM_P_FDR.txt")
Ota_coloc_list2 = fread_FT("data/211208_colocalizedSNP_stringent.txt")
coloc_genotype_SLE = fread_FT("data/210925_SLE_colocalizedSNP_genotype_info.txt")

Ota_coloc_tmp = Ota_coloc_list2 %>% filter(subset==subset_tmp)
    
data_all = fread_n(count_path)
data_tmp =  data_all[intersect(Ota_coloc_tmp$GENE,rownames(data_all)),]

res_tmp = res_sum %>% filter(subset==subset_tmp)

DEGs_tmp = res_tmp %>% filter(state_FDR<0.05|activity_FDR<0.05)
DEGs_list = intersect(DEGs_tmp$genes, rownames(data_all))

for(kkk in 1:nrow(data_tmp)){
    data_tmp2 = data_tmp[kkk,] %>% t() %>% data.frame() %>% rownames_to_column("id") %>%
                mutate(id=take_factor(id,1,"_"))
    gene_tmp = colnames(data_tmp2)[2]
    print(paste0(subset_tmp, "; eGene: ",gene_tmp))
    
    Ota_coloc_tmp2 = Ota_coloc_tmp %>% filter(GENE==gene_tmp) %>% .[1,] # avoid redundance
    rs_tmp = Ota_coloc_tmp2$rs　# only one rsID

    genotype_tmp = coloc_genotype_SLE %>% select(id,age,gender,PSLmg,HCQ,MMF,TAC,all_of(rs_tmp))

    gene_list_2 = setdiff(DEGs_list,gene_tmp)

        for(jjj in 1:length(gene_list_2)){
            pGene_tmp = gene_list_2[jjj]
            print(paste0("pGene: ",pGene_tmp))
            pGene_exp = data_all[pGene_tmp,] %>% t() %>% data.frame() %>% rownames_to_column("id") %>%
                    mutate(id=take_factor(id,1,"_"))
            data_tmp3 = right_join(data_tmp2,genotype_tmp,by="id") %>% left_join(.,pGene_exp,by="id")
            colnames(data_tmp3)[2]="eGene"
            colnames(data_tmp3)[9]="eQTL"
            colnames(data_tmp3)[10]="pGene"
            data_tmp3$eQTL = as.numeric(data_tmp3$eQTL)

            full_tmp = lm(eGene~eQTL+age+gender+PSLmg+HCQ+MMF+TAC+pGene+eQTL:pGene,data=data_tmp3)
            null_tmp = lm(eGene~eQTL+age+gender+PSLmg+HCQ+MMF+TAC+pGene,data=data_tmp3)
            interaction_beta=summary(full_tmp)$coef[10,1]
            anova_p = anova(full_tmp,null_tmp)[2,6]
            anova_tmp = data.frame(subset=subset_tmp,rsID=rs_tmp,eGene=gene_tmp,pGene=pGene_tmp,int_beta=interaction_beta,anovaP=anova_p)
            if(jjj==1){anova_sum=anova_tmp}else{anova_sum=rbind(anova_sum,anova_tmp)}
        }

        anova_sum2=anova_sum%>%left_join(.,res_tmp%>%select(genes,siggroup005),by=c("pGene"="genes"))

        if(kkk==1){anova_sum3=anova_sum2}else{anova_sum3=rbind(anova_sum3,anova_sum2)}
        }

write.table_FT_2(anova_sum3,paste0("context_eQTL/stringent/",today(),"_COI_SLE_",subset_tmp,"_contexteQTL_anova_stringent.txt"))
######################################################################################################################
######################################################################################################################

list1 = make_list("context_eQTL/stringent","_contexteQTL_anova_stringent.txt")
list1$subset = take_factor(list1$FILE,4:5,"_")

for(iii in 1:nrow(list1)){
    subset_tmp = list1$subset[iii]
    data_tmp1 = fread_FT(list1$PATH[iii])
    if(iii==1){data_sum1=data_tmp1}else{data_sum1=rbind(data_sum1,data_tmp1)}
   }
data_sum1$FDR = p.adjust(data_sum1$anovaP,method="BH")
dim(data_sum1)
# [1] 583015      8
write.table_FT_2(data_sum1,paste0("context_eQTL/",today(),"_contexteQTL_stringent_anovasum.txt"))


# Representative
iii=3
    subset_tmp = data_sum1_sig$subset[iii]
    eGene_tmp = data_sum1_sig$eGene[iii]
    pGene_tmp = data_sum1_sig$pGene[iii]
    eQTL_tmp = data_sum1_sig$rsID[iii]

    list_tmp = list %>% filter(subset==subset_tmp)
    data_tmp = fread_n(list_tmp$PATH)[c(eGene_tmp,pGene_tmp),] %>% t() %>% as.data.frame()
    data_tmp$id = take_factor(rownames(data_tmp),1,"_") 
    data_tmp2 = data_tmp %>% left_join(.,coloc_genotype_SLE,by="id") %>% 
                select(id,all_of(eQTL_tmp),all_of(eGene_tmp),all_of(pGene_tmp))

    colnames(data_tmp2) = c("id","rsID","eGene","pGene")
    data_tmp2 = data_tmp2[!is.na(data_tmp2$rsID),]
    data_tmp2 = data_tmp2[order(data_tmp2$pGene),]
    half = floor(nrow(data_tmp2)/2) +1

    data_tmp2$group = c(rep("low",half),rep("high",nrow(data_tmp2)-half))

    data_tmp2$rsID =factor(data_tmp2$rsID,levels=c("0","1","2"))
    data_tmp2$group =factor(data_tmp2$group,levels=c("low","high"))
    max = floor(max(data_tmp2$eGene))+1
    labels = c(paste0(pGene_tmp," low"),paste0(pGene_tmp," high"))

p = ggplot(data=data_tmp2, aes(x=rsID,y=eGene,fill=factor(group,labels=labels)))+
     geom_boxplot(outlier.shape = NA)+
     #geom_point(size=0.5,color="black",position=position_jitterdodge())+
     theme_classic()+
     facet_wrap(~group,ncol=2) +
     scale_fill_manual(values=c("#BC3C2980","#BC3C29FF"))+
     theme(axis.text.x=element_text(colour="black",angle=45,hjust=1,vjust=1,size=14),
           axis.text.y=element_text(colour="black",size=14),
           axis.title.x=element_text(colour="black",size=14),
           axis.title.y=element_text(colour="black",size=14),
           plot.title=element_blank(),
           strip.background=element_blank(),
           strip.text=element_blank(),
           panel.background =element_rect(fill="transparent",colour="black",size=1),
           panel.grid = element_blank(),
           legend.position="right",
           legend.title=element_blank(),
           legend.text=element_text(colour="black",size=14))+
     ylim(0,max)+
     scale_x_discrete(labels= c("TA/TA","TA/T","T/T"))+
     labs(x="rs36059542",y=paste0("log(CPM+1)"))

    pdf_3(paste0("context_eQTL/plot_PB/",today(),"_COI_SLE_",subset_tmp,"_",eQTL_tmp,"_",eGene_tmp,"_",pGene_tmp,"_facet.pdf"),h=3.6,w=4.5)
     plot(p)
    dev.off()

######## Supple
data_sum1=fread_FT("context_eQTL/211209_contexteQTL_stringent_anovasum.txt")
data_sum1$subset = factor(data_sum1$subset,levels=celltype_reordered_27)
data_sum1= data_sum1[order(data_sum1$eGene),]
data_sum1= data_sum1[order(data_sum1$subset),]

data_sum1 =data_sum1 %>% left_join(.,label_df,by="subset") %>%　select(rsID,label,eGene,pGene,siggroup005,int_beta,anovaP,FDR)
data_sum1$siggroup005=ifelse(data_sum1$siggroup005=="discordant","only_state",data_sum1$siggroup005)
data_sum1$siggroup005=factor(data_sum1$siggroup005,levels=c("only_activity","significant_in_both","only_state"))
levels(data_sum1$siggroup005)=c("Disease-activity","Both significant","Disease-state")
data_sum1= data_sum1[order(data_sum1$siggroup005),]

data_sum1$int_beta = formatC(data_sum1$int_beta,digits=2)
data_sum1$anovaP = formatC(data_sum1$anovaP,digits=2)
data_sum1$FDR = formatC(data_sum1$FDR,digits=2)
colnames(data_sum1)=c("GWAS-eQTL","Cell type","eGene","pGene","Signature","Interaction effect size","Pvalue","FDR")

write.table_FT_2(data_sum1,paste0(today(),"_contexteQTL_forsupple.txt"))


