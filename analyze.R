library('openxlsx')
library('sva')
library('org.Hs.eg.db')
library('clusterProfiler')
library('dplyr')
library('tidyr')
library('reshape2')
library('GSEABase')
library('GSVA')
library('ConsensusClusterPlus')
library('timeROC')
library('RCircos')
library('stringr')
library('ggsignif')
library('ggplot2')
library('pca3d')
library('pheatmap')
library('limma')
library('survival')
library('survminer')
library('Hmisc')
library('ggstance')
library('ggpubr')
library('forcats')
library('maftools')
library('car')
library('ggstatsplot')
library('forestplot')
library('ggalluvial')
library('corrplot')
library('rms')
library('regplot')
library('glmnet')
library('hdf5r')
library('estimate')
library('pRRophetic')
library('WGCNA')
library('pROC')
library('xCell')
library('RColorBrewer')
library('circlize')
library('ComplexHeatmap')
library('rmda')
library('ggvenn')
library('ggdist')
library('gghalves')
library('plyr')
library('cgdsr')

rm(list=ls())

##1、公共数据库下载乳腺癌转录组、基因组、细胞系数据；下载铁死亡因子(Ferroptosis-genes,FerrDb数据库、PMID：32760210、doi: 10.3389/fonc.2020.590861取并集）
#读取整个项目所有数据集
setwd('D:/bioinfo_softwares')
#*TCGA skcm
pheno = read.delim('原始数据/黑色素瘤/TCGA-SKCM/TCGA-SKCM.GDC_phenotypeCelltype_Annotation.tsv', header=T, sep='\t')
rna_mean = read.delim(gzfile('原始数据/黑色素瘤/TCGA-SKCM/TCGA-SKCM.htseq_fpkmCelltype_Annotation.tsv.gz', 'rt'), header=T, sep='\t', stringsAsFactors=F, check.names = F)
cnv_mean = read.delim(gzfile('原始数据/黑色素瘤/TCGA-SKCM/TCGA-SKCM.gisticCelltype_Annotation.tsv.gz', 'rt'), header=T, sep='\t', stringsAsFactors=F, check.names = F)
snv = read.maf('原始数据/黑色素瘤/TCGA-SKCM/TCGA.SKCM.mutect.4b7a5729-b83e-4837-9b61-a6002dce1c0a.DR-10.0.somatic.maf.gz')
#rna_count_skcm = read.delim(gzfile('原始数据/黑色素瘤/TCGA-skcm/TCGA-skcm.htseq_counts.tsv.gz', 'rt'), header=T, sep='\t', stringsAsFactors=F, check.names = F)
tmb_raw = read.delim('原始数据/泛癌/TCGA-PANCAN/mutation-load_updated.txt', header=T, sep='\t', stringsAsFactors=F, check.names = F)
pheno3 = openxlsx::read.xlsx(xlsxFile = '基础库文件/TCGAsubtype_PMID29628290_NIHMS958212-supplement-2.xlsx', sheet=1, check.names=F)
#*GEO1-GSE78220
geo1_rna = read.delim(gzfile('原始数据/黑色素瘤/GSE78220/GSE78220_PatientFPKM.txt.gz', 'rt'), header=T, sep='\t', stringsAsFactors=F, check.names = F)
geo1_info = read.delim('原始数据/黑色素瘤/GSE78220/SampleInfo.xls', header=T, sep='\t', stringsAsFactors=F, check.names = F)
#*GEO2-GSE91061
geo2_rna = read.delim(gzfile('原始数据/黑色素瘤/GSE91061/GSE91061_BMS038109Sample.hg19KnownGene.fpkmCelltype_Annotation.txt.gz', 'rt'), header=T, sep='\t', stringsAsFactors=F, check.names = F)
geo2_info = read.delim('原始数据/黑色素瘤/GSE91061/GSE91061_infoCelltype_Annotation.txt', header=T, sep='\t', stringsAsFactors=F, check.names = F)
#*GEO3-GSE65904
geo3_rna = read.delim(gzfile('原始数据/黑色素瘤/GSE65904/Merge_GPL10558.expoCelltype_Annotation.txt.gz', 'rt'), header=T, sep='\t', stringsAsFactors=F, check.names = F)
geo3_info = read.delim('原始数据/黑色素瘤/GSE65904/SampleInfo.xls', header=T, sep='\t', stringsAsFactors=F, check.names = F)
#*GEO4-GSE54467
geo4_rna = read.delim(gzfile('原始数据/黑色素瘤/GSE54467/Merge_GPL6884.expro.txt.gz', 'rt'), header=T, sep='\t', stringsAsFactors=F, check.names = F)
geo4_info = read.delim('原始数据/黑色素瘤/GSE54467/SampleInfo.xls', header=T, sep='\t', stringsAsFactors=F, check.names = F)
#*GEO5-GSE59334
geo5_rna = read.delim(gzfile('原始数据/黑色素瘤/GSE59334/Merge_GPL15018.expro.txt.gz', 'rt'), header=T, sep='\t', stringsAsFactors=F, check.names = F)
geo5_info = read.delim('原始数据/黑色素瘤/GSE59334/SampleInfo.xls', header=T, sep='\t', stringsAsFactors=F, check.names = F)
#*GEO6-GSE19234
geo6_rna = read.delim(gzfile('原始数据/黑色素瘤/GSE19234/Merge_GPL570.expro.txt.gz', 'rt'), header=T, sep='\t', stringsAsFactors=F, check.names = F)
geo6_info = read.delim('原始数据/黑色素瘤/GSE19234/SampleInfo.xls', header=T, sep='\t', stringsAsFactors=F, check.names = F)
#*GEO7-GSE22153
geo7_rna = read.delim(gzfile('原始数据/黑色素瘤/GSE22153/Merge_GPL6102.exproCelltype_Annotation.txt.gz', 'rt'), header=T, sep='\t', stringsAsFactors=F, check.names = F)
geo7_info = read.delim('原始数据/黑色素瘤/GSE22153/SampleInfo.xls', header=T, sep='\t', stringsAsFactors=F, check.names = F)
#*GEO8-GSE66839
geo8_rna = read.delim(gzfile('原始数据/黑色素瘤/GSE66839/Merge_GPL4091.expro.txt.gz', 'rt'), header=T, sep='\t', stringsAsFactors=F, check.names = F)
geo8_info = read.delim('原始数据/黑色素瘤/GSE66839/SampleInfo.xls', header=T, sep='\t', stringsAsFactors=F, check.names = F)
#*GEO9-GSE93157
geo9_rna = read.delim(gzfile('原始数据/黑色素瘤/GSE93157/Merge_GPL19965.expro.txt.gz', 'rt'), header=T, sep='\t', stringsAsFactors=F, check.names = F)
geo9_info = read.delim('原始数据/黑色素瘤/GSE93157/SampleInfo.xls', header=T, sep='\t', stringsAsFactors=F, check.names = F)
#*GEO10-GSE51547
geo10_rna = read.delim(gzfile('原始数据/黑色素瘤/GSE51547/Merge_GPL13534.exproCelltype_Annotation.txt.gz', 'rt'), header=T, sep='\t', stringsAsFactors=F, check.names = F)
geo10_info = read.delim('原始数据/黑色素瘤/GSE51547/SampleInfo.xls', header=T, sep='\t', stringsAsFactors=F, check.names = F)
#*GEO11-GSE100797
geo11_rna = read.delim(gzfile('原始数据/黑色素瘤/GSE100797/GSE100797_ProcessedData.txt.gz', 'rt'), header=T, sep='\t', stringsAsFactors=F, check.names = F)
geo11_info = read.delim('原始数据/黑色素瘤/GSE100797/SampleInfo.xls', header=T, sep='\t', stringsAsFactors=F, check.names = F)
#*GEO12-GSE59455
geo12_rna = read.delim(gzfile('原始数据/黑色素瘤/GSE59455/Merge_GPL8432.expro.txt.gz', 'rt'), header=T, sep='\t', stringsAsFactors=F, check.names = F)
geo12_info = read.delim('原始数据/黑色素瘤/GSE59455/SampleInfo.xls', header=T, sep='\t', stringsAsFactors=F, check.names = F)

#*hg19
#hg19_1 = read.delim('基础库文件/protein_coding.hg19.position.txt', sep='\t', stringsAsFactors=F)
#*ICG
icgs = read.delim('基础库文件/Immune checkpoint genes (ICGs).txt', header=F, sep='\t')$V1
#*
 = read.delim('基础库文件/.txt', header=F, sep='\t')$V1
#*CD8 T cell-related IFN-gamma pathway markers
ifn_genes = c('IFNGR1', 'PIAS1', 'JAK2', 'IFNG', 'STAT1', 'SUMO1', 'IFNGR2', 'JAK1', 'PTPN2', 'PTPN11', 'SOCS3', 'PTPN1')
#*immune cells recruitment
immunecell_genes = c('CCL19', 'CCL18', 'CCL13', 'CXCR3', 'CCL14', 'CCL5', 'CCL4', 'CCL8', 'CCL7', 'CXCL10', 'CCL2', 'CCR2', 'CXCL16', 'CXCL11', 'CX3CL1', 'CCL20', 'CCL21', 'CXCL13', 'CCL22', 'CCR4', 'CCR5', 'CXCR4', 'CXCL9', 'CXCL12', 'CXCR6', 'CXCL14', 'CSF3R', 'CXCL2', 'CXCL1', 'CXCL5', 'CX3CR1', 'CXCR2', 'CXCL6', 'CXCL3')
#*immune suppression genes
immunesuppression_genes = c('CD209', 'CTLA4', 'CX3CR1', 'ENTPD1', 'FAS', 'TNFRSF8', 'FGFR1', 'FOXP3', 'FUT4', 'IDO1', 'IL10', 'RELA', 'TGFB1', 'TNF', 'TNFSF9', 'VCAM1', 'PDCD1', 'PDCD1LG2', 'CD274', 'LAG3', 'HAVCR2', 'LGALS9')
#*innate immunity
innateimmunity_genes = c('BCL2L1', 'FADD', 'GNLY', 'CCR2', 'IFIH1', 'CASP1', 'IFNA1', 'GZMM', 'CASP9', 'CASP8', 'IL1B', 'TNFRSF1B', 'PSMB10', 'ISG15', 'PSMB9', 'TNFAIP3', 'IRF7', 'IRF3', 'ICAM2', 'ICAM3', 'OAS1', 'OAS2', 'OAS3', 'OASL', 'TLR3', 'TLR4', 'TLR7', 'TLR8')
#*adaptive immunity
adaptiveimmunity_genes = c('TNFRSF14', 'ICOS', 'CD3E', 'IFITM2', 'IFITM1', 'CD8A', 'PTEN', 'ICAM1', 'CD40LG', 'IFNG', 'LTB', 'CXCL13', 'CD80', 'IL4', 'IL2', 'IL10', 'CD40', 'CD28', 'CD19')
#*antigen presentation and processing
app_genes = c('B2M', 'CD74', 'CTSS', 'HERC6', 'HLA-A', 'HLA-B', 'HLA-C', 'HLA-DMA', 'HLA-DMB', 'HLA-DOA', 'HLA-DOB', 'HLA-DPA1', 'HLA-DPB1', 'HLA-DQA1', 'HLA-DQA2', 'HLA-DQB1', 'HLA-DRA', 'HLA-DRB1', 'HLA-DRB5', 'HLA-E', 'HLA-F', 'ICAM1', 'ICAM2', 'ICAM3', 'ICOSLG', 'LILRA1', 'LILRA3', 'MICA', 'PSMB10', 'PSMB5', 'PSMB8', 'PSMB9', 'TAP1', 'TAP2', 'TAPBP')
#*cytotoxicity/killing of cancer cells
cytotox_genes = c('GNLY', 'GZMA', 'GZMB', 'GZMH', 'GZMK', 'IFIT2', 'IFITM2', 'KIR2DL3', 'KIR3DL1', 'KIR3DL2', 'KLRB1', 'KLRD1', 'KLRK1', 'NCAM1', 'NCR1', 'NKG7', 'PRF1')
#*inflammation
inflammation_genes = c('IL17A', 'IL6', 'IL6R', 'IL1A', 'IL1B', 'CD40LG', 'BMP2', 'CD6', 'CEBPB', 'STAT1', 'CSF1R', 'CSF1', 'STAT4', 'C5AR1', 'TNF', 'C5', 'TREM1', 'IL15', 'IL18', 'NLRP3', 'IL33')
#*趋化因子
#chemokines = read.xlsx('基础库文件/趋化因子/PMID_35479097_chemokines.xlsx')
#*编码gene
coding_genes = read.delim('基础库文件/Homo_sapiens.GRCh37.87.protein_codingCelltype_Annotation.gtf', header=T, sep='\t', stringsAsFactors=F, check.names = F)$gene_name
#*TIP genes
colds = c('CXCL1', 'CXCL2', 'CCL20')
hots = c('CXCR3', 'CXCR4', 'CXCL9', 'CXCL10', 'CXCL11', 'CCL5', 'CD3E', 'CD4', 'CD8A', 'CD8B', 'CD274', 'PDCD1')
tips = c(colds, hots)

#数据预处理
#*TCGA-skcm + TCGA-READ合并
#**统一、整理TCGA各数据集ID
tcga_tumor_ids = as.character(pheno$sample[pheno$sample_type.samples=='Primary Tumor' | pheno$sample_type.samples=='Metastatic' | pheno$sample_type.samples=='Recurrent Tumor' | pheno$sample_type.samples=='Additional Metastatic'])
tcga_normal_ids = as.character(pheno$sample[pheno$sample_type.samples=='Solid Tissue Normal'])
tcga_all_ids = union(tcga_tumor_ids, tcga_normal_ids)
#**ID过滤的样本集
tcga_rna_tumor_mean = rna_mean[,(colnames(rna_mean) %in% tcga_tumor_ids)]
#**cnv
cnv_tumor_mean = cnv_mean[,(colnames(cnv_mean) %in% tcga_tumor_ids)]
cnv_tumor = cbind(rownames(cnv_tumor_mean), cnv_tumor_mean)
colnames(cnv_tumor)[1] = 'SYMBOL'
#*TMB
tmb = tmb_raw[(tmb_raw$Cohort %in% c('SKCM')),]
tmb = merge(tmb,
            data.frame('sample'=pheno$sample, 'Tumor_Sample_ID'=substring(pheno$sample,1,15)),
            by='Tumor_Sample_ID')
#*pheno
mycgds = CGDS('http://www.cbioportal.org/')
test(mysgds)
all_tcga_studies = getCancerStudies(mycgds)
pheno_cgds=getClinicalData(mycgds, getCaseLists(mycgds, 'skcm_tcga')[1,1])
rownames(pheno_cgds) = gsub('\\.', '-', rownames(pheno_cgds))
pheno_cgds$'sample2' = rownames(pheno_cgds)
pheno_cgds2 = merge(pheno_cgds, pheno %>% mutate('sample2'=substr(sample, 1, 15)) %>% select('sample', 'sample2'), by='sample2')
pheno_cgds2$'DFS' = mapvalues(pheno_cgds2$DFS_STATUS,c("0:DiseaseFree", "1:Recurred/Progressed", ""), c(0,1,NA)) %>% as.numeric()
pheno_cgds2$'DFS.time' = pheno_cgds2$DFS_MONTHS
#**pheno3分子分型
pheno31 = pheno3[(pheno3$TCGA.Study=='SKCM'),]
pheno31 = merge(pheno31,
                data.frame('sample'=pheno$sample, 'TCGA.Participant.Barcode'=substring(pheno$sample,1,12)),
                by='TCGA.Participant.Barcode', all.y=T)

#*GEO1-GSE78220
geo1_rna_mean = geo1_rna
rownames(geo1_rna_mean) = geo1_rna_mean$Gene
colnames(geo1_rna_mean) = strsplit2(colnames(geo1_rna_mean), split = '\\.')[,1]
geo1_rna_mean = geo1_rna_mean[,-1]
geo1_rna_tumor_mean = geo1_rna_mean
#*GEO2-GSE91061
geo2_rna_mean = geo2_rna
geo2_rna_tumor_mean = geo2_rna_mean
#*GEO3-GSE65904
geo3_rna_mean = geo3_rna
geo3_rna_tumor_mean = geo3_rna_mean
rownames(geo3_rna_tumor_mean) = str_replace_all(rownames(geo3_rna_tumor_mean), 'TUBB4', 'TUBB4A')
#*GEO4-GSE54467
geo4_rna_mean = geo4_rna[,-1]
rownames(geo4_rna_mean) = geo4_rna$Samples
geo4_rna_tumor_mean = geo4_rna_mean
#*GEO5-GSE59334
geo5_rna_mean = geo5_rna[,-1]
rownames(geo5_rna_mean) = geo5_rna$Samples
geo5_rna_tumor_mean = geo5_rna_mean
#*GEO6-GSE19234
geo6_rna_mean = geo6_rna[,-1]
rownames(geo6_rna_mean) = geo6_rna$Samples
geo6_rna_tumor_mean = geo6_rna_mean
#*GEO7-GSE22153
geo7_rna_mean = geo7_rna
geo7_rna_tumor_mean = geo7_rna_mean
rownames(geo7_rna_tumor_mean) = str_replace_all(rownames(geo7_rna_tumor_mean), 'TUBB4', 'TUBB4A')
#*GEO8-GSE66839
geo8_rna_mean = geo8_rna[,-1]
rownames(geo8_rna_mean) = geo8_rna$Samples
geo8_rna_tumor_mean = geo8_rna_mean
#*GEO9-GSE93157
geo9_rna_mean = geo9_rna[,-1]
rownames(geo9_rna_mean) = geo9_rna$Samples
geo9_rna_tumor_mean = geo9_rna_mean
#*GEO10-GSE51547
geo10_rna_mean = geo10_rna
geo10_rna_tumor_mean = geo10_rna_mean
#*GEO11-GSE100797
geo11_rna_mean = geo11_rna[,-1]
rownames(geo11_rna_mean) = geo11_rna[,1]
geo11_rna_tumor_mean = geo11_rna_mean
#*GEO12-GSE59455
geo12_rna_mean = geo12_rna[,-1]
rownames(geo12_rna_mean) = geo12_rna[,1]
geo12_rna_tumor_mean = geo12_rna_mean


#**整理数据集
combs_genes = intersect(combs_genes, coding_genes)
tcga_rna_tumor_mean = tcga_rna_tumor_mean[intersect(rownames(tcga_rna_tumor_mean), coding_genes),]


##主干调试
##Section1：肿瘤免疫表型得分关联临床
##1、表达谱计算TIPscore(方法参考末尾文献)；
##K-M曲线分析TIPscore与研究癌症的预后关系（分析不同的生存期）
tcga_rna_tumor_tips = tcga_rna_tumor_mean[tips,]
tcga_rna_tumor_tips_z = apply(tcga_rna_tumor_tips, 1, function(x){(x-mean(x))/sd(x)})
tcga_rna_tumor_tips_z = as.data.frame(tcga_rna_tumor_tips_z)
tips_tcga = apply(tcga_rna_tumor_tips_z, 1, sum) %>% data.frame(stringsAsFactors = F) %>%
  rename('TIPscore'='.') %>%
  mutate(sample=rownames(.))
tips_tcga$'TIP_cluster' = ifelse(tips_tcga$TIPscore>=median(tips_tcga$TIPscore), 'High TIPscore', 'Low TIPscore')
tips_tcga2 = merge(tips_tcga, pheno[,c('sample', 'OS', 'OS.time')], by='sample')
write.table(tips_tcga, 'output/1.1_tips_tcga.txt', sep = '\t', quote = F, row.names = F, col.names = T)

plot_surv = function(data2, labels, i, legend_title) {
  data2 <<- data2
  surv_fit2 <<- survfit(Surv(OS.time, OS)~cluster, data=data2)
  print(ggsurvplot(surv_fit2, pval=T, 
                   ggtheme = theme(axis.line = element_line(linetype = "solid"), 
                                   #axis.text = element_text(hjust = 1), 
                                   #axis.text.x = element_text(size = 11, vjust = 0.5, angle = 90),
                                   panel.background = element_rect(fill = NA)
                   ),
                   risk.table=T, 
                   legend.labs=labels,
                   #palette = c('red', 'blue'),
                   palette = ann_colors[[i]],
                   legend.title=legend_title,
                   #legend.position = c(0.05,0.15),
                   xlab='Time in months'
  ))
}
ann_colors = list(cluster=c('High TIPscore'='red', 'Low TIPscore'='blue'))
plot_surv(tips_tcga2 %>% dplyr::rename('cluster'='TIP_cluster'), c('High TIPscore', 'Low TIPscore'), 1, '')

##Section2：肿瘤表型相关预后分层系统（TIPRGPI）的构建
##1、基于表达谱和加权基因共表达网络分析（WGCNA）筛选TIPscore性状相关的关键模块（可多选；可适当加入其它性状）；
##3、基于关键模块基因单因素cox回归分析筛选预后因子，lasso或者多因素cox回归分析进一步筛选构建TIPRGPI的最优基因组合，进行K-M曲线展示（2、3图可以组图；K-M曲线展示模型基因）。
#**WGCNA
tcga_rna_all_mean_wgcna2 = as.data.frame(t(tcga_rna_tumor_mean))
## 检测缺失值
gsg = goodSamplesGenes(tcga_rna_all_mean_wgcna2, verbose = 3)
if (!gsg$allOK) {
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", 
                     paste(names(tcga_rna_all_mean_wgcna2)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", 
                     paste(rownames(tcga_rna_all_mean_wgcna2)[!gsg$goodSamples], collapse = ",")));
  # Remove the offending genes and samples from the data:
  tcga_rna_all_mean_wgcna2 = tcga_rna_all_mean_wgcna2[gsg$goodSamples, gsg$goodGenes]
}  ##1104样本*1311基因
## 查看是否有离群样品
sampleTree = hclust(dist(tcga_rna_all_mean_wgcna2), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
##软阈值筛选
powers = c(c(1:20), seq(from = 20, to=30, by=2))
#sft = pickSoftThreshold(tcga_rna_tumor_mean_deg_all2, powerVector=powers, networkType="signed", verbose=5, RsquaredCut = 0.85)
sft02 = pickSoftThreshold(tcga_rna_all_mean_wgcna2, powerVector=powers, networkType="signed", verbose=5, RsquaredCut = 0.8, blockSize = 10000)
par(mfrow = c(1,2))
cex1 = 0.9
# 横轴是Soft threshold (power)，纵轴是无标度网络的评估参数，数值越高，
# 网络越符合无标度特征 (non-scale)
plot(sft02$fitIndices[,1], -sign(sft02$fitIndices[,3])*sft02$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft02$fitIndices[,1], -sign(sft02$fitIndices[,3])*sft02$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# 筛选标准。R-square=0.95
abline(h=0.8,col="red")
# Soft threshold与平均连通性
plot(sft02$fitIndices[,1], sft02$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft02$fitIndices[,1], sft02$fitIndices[,5], labels=powers, 
     cex=cex1, col="red")
power02 = sft02$powerEstimate
#查看是否符合无尺度分布
k <- softConnectivity(tcga_rna_all_mean_deg_all202, power=power02)
hist(k)
scaleFreePlot(k,main="Check Scale free topology\n")
##一步法网络构建：One-step network construction and module detection##
# power: 上一步计算的软阈值
# maxBlockSize: 计算机能处理的最大模块的基因数量 (默认5000)；
#  4G内存电脑可处理8000-10000个，16G内存电脑可以处理2万个，32G内存电脑可以处理3万个
#  计算资源允许的情况下最好放在一个block里面。
# corType: pearson or bicor
# numericLabels: 返回数字而不是颜色作为模块的名字，后面可以再转换为颜色
# saveTOMs：最耗费时间的计算，存储起来，供后续使用
# mergeCutHeight: 合并模块的阈值，越大模块越少
corType = "pearson"
corFnc = ifelse(corType=="pearson", cor, bicor)
# 对二元变量，如样本性状信息计算相关性时，或基因表达严重依赖于疾病状态时，需设置下面参数
maxPOutliers = ifelse(corType=="pearson",1,0.05)
# 关联样品性状的二元变量时，设置
robustY = ifelse(corType==corType,T,F)
net02 = blockwiseModules(tcga_rna_all_mean_wgcna2, power = power02, maxBlockSize = 6000,
                         TOMType = 'unsigned', minModuleSize = 100,
                         reassignThreshold = 0, mergeCutHeight = 0.25,
                         numericLabels = TRUE, pamRespectsDendro = FALSE,
                         saveTOMs=TRUE, corType = corType, 
                         maxPOutliers=maxPOutliers, loadTOMs=TRUE,
                         saveTOMFileBase = paste0('output/2.1_wgcna_net02', ".tom"),
                         verbose = 3)
# 根据模块中基因数目的多少，降序排列，依次编号为 `1-最大模块数`。
# **0 (grey)**表示**未**分入任何模块的基因。 
table(net02$colors)
##层级聚类树展示各个模块
## 灰色的为**未分类**到模块的基因。
# Convert labels to colors for plotting
moduleLabels02 = net02$colors
moduleColors02 = labels2colors(moduleLabels02)
# Plot the dendrogram and the module colors underneath
# 如果对结果不满意，还可以recutBlockwiseTrees，节省计算时间
plotDendroAndColors(net02$dendrograms[[1]], moduleColors02[net02$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
##绘制模块之间相关性图
# module eigengene, 可以绘制线图，作为每个模块的基因表达趋势的展示
MEs02 = net02$MEs
### 不需要重新计算，改下列名字就好
### 官方教程是重新计算的，起始可以不用这么麻烦
MEs_col02 = MEs02
colnames(MEs_col02) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEs02),"ME",""))))
MEs_col02 = orderMEs(MEs_col02)
# 根据基因间表达量进行聚类所得到的各模块间的相关性图
# marDendro/marHeatmap 设置下、左、上、右的边距
plotEigengeneNetworks(MEs_col02, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90)
## 模块与表型数据关联
pheno2 = pheno[,c('sample', 'age_at_index.demographic', 'gender.demographic', 'tumor_stage.diagnoses', 'pathologic_T', 'pathologic_N', 'pathologic_M', 'breslow_depth_value', 'melanoma_clark_level_value')]
pheno2$tumor_stage.diagnoses = Recode(pheno2$tumor_stage.diagnoses, "c('stage ia','stage ib')='stage i';c('stage iia','stage iib','stage iic')='stage ii';c('stage iiia','stage iiib','stage iiic')='stage iii';c('i/ii nos', 'not reported', '')=NA")
pheno2$pathologic_T = Recode(pheno2$pathologic_T, "c('T1a', 'T1b')='T1';c('T2a', 'T2b')='T2';c('T3a', 'T3b')='T3';c('T4a', 'T4b')='T4';c('')=NA")
pheno2$pathologic_N = Recode(pheno2$pathologic_N, "c('N1a', 'N1b')='N1';c('N2a', 'N2b', 'N2c')='N2';c('')=NA")
pheno2$pathologic_M = Recode(pheno2$pathologic_M, "c('M1a', 'M1b', 'M1c')='M1';c('')=NA")
pheno2$melanoma_clark_level_value = Recode(pheno2$melanoma_clark_level_value, "''=NA")
pheno2$'Age' = ifelse(pheno2$age_at_index.demographic>=median(pheno2$age_at_index.demographic), '>=60', '<60')
pheno2$'Gender' = ifelse(pheno2$gender.demographic=='male', 1, 0)
annotation_col_wgcna = data.frame('sample'=tips_tcga2$sample,
                                  'TIPscore'=ifelse(tips_tcga2$TIP_cluster=='High TIPscore', 1, 0))
annotation_col_wgcna = merge(annotation_col_wgcna, pheno2[,c('sample', 'age_at_index.demographic', 'Gender', 'breslow_depth_value')], by='sample')
rownames(annotation_col_wgcna) = tips_tcga2$sample
annotation_col_wgcna = annotation_col_wgcna %>% arrange(rownames(.))
annotation_col_wgcna = select(annotation_col_wgcna, -c('sample'))
write.table(annotation_col_wgcna, file='output/2.1_annotation_col_wgcna.txt', sep='\t', quote=F)
corType="pearson"
if (corType=="pearson") {
  aa = MEs_col02[(rownames(MEs_col02) %in% rownames(annotation_col_wgcna)),]
  aa = aa[order(rownames(aa)),]
  modTraitCor02 = cor(aa, annotation_col_wgcna, use = "p")#, use = "p"
  modTraitP02 = corPvalueStudent(modTraitCor02, nrow(annotation_col_wgcna))
} else {
  aa = MEs_col[(rownames(MEs_col) %in% rownames(annotation_col_wgcna)),]
  aa = aa[order(rownames(aa)),]
  modTraitCorP = bicorAndPvalue(aa, annotation_col_wgcna, robustY=robustY)
  modTraitCor = modTraitCorP$bicor
  modTraitP   = modTraitCorP$p
}
write.table(aa, file='output/2.1_MEs_col.txt', sep='\t', quote=F)
write.table(modTraitCor02, file='output/2.1_modTraitCor.txt', sep='\t', quote=F)
write.table(tcga_rna_all_mean_wgcna2, file='output/2.1_expr.txt', sep='\t', quote=F)
## Warning in bicor(x, y, use = use, ...): bicor: zero MAD in variable 'y'.
## Pearson correlation was used for individual columns with zero (or missing)
## MAD.
# signif表示保留几位小数
par(mfrow = c(1,1))
textMatrix02 = paste(signif(modTraitCor02, 2), "\n(", signif(modTraitP02, 1), ")", sep = "")
dim(textMatrix02) = dim(modTraitCor02)
pdf("output/2.1_heatmap.pdf", width=10, height=10, onefile = F)
tiff("output/2.1_heatmap.tiff", res=300, width=10, height=10, compression="lzw", units="in")
labeledHeatmap(Matrix = modTraitCor02, xLabels = colnames(annotation_col_wgcna), 
               yLabels = colnames(MEs_col02), 
               cex.lab = 0.5, 
               ySymbols = colnames(MEs_col02), colorLabels = FALSE, 
               colors = blueWhiteRed(50), 
               textMatrix = textMatrix02, setStdMargins = TRUE, 
               cex.text = 0.5, zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()
##筛选与相关模块的hub基因
## 从上图可以看到MEmagenta与Insulin_ug_l相关
## 模块内基因与表型数据关联
# 性状跟模块虽然求出了相关性，可以挑选最相关的那些模块来分析，但是模块本身仍然包含非常多的基因，还需进一步的寻找最重要的基因。
# 所有的模块都可以跟基因算出相关系数，所有的连续型性状也可以跟基因的表达值算出相关系数。
# 如果跟性状显著相关基因也跟某个模块显著相关，那么这些基因可能就非常重要。
### 计算模块与基因的相关性矩阵
if (corType=="pearson") {
  geneModuleMembership02 = as.data.frame(cor(tcga_rna_all_mean_wgcna2, MEs_col02, use = "p"))
  MMPvalue02 = as.data.frame(corPvalueStudent(
    as.matrix(geneModuleMembership02), nrow(tcga_rna_all_mean_wgcna2)))
} else {
  geneModuleMembershipA = bicorAndPvalue(dataExpr, MEs_col, robustY=robustY)
  geneModuleMembership = geneModuleMembershipA$bicor
  MMPvalue   = geneModuleMembershipA$p
}
plotModuleSignificance(abs(as.numeric(cor(tcga_rna_all_mean_wgcna2, MEs_col02, use = "p"))), colors = unlist(lapply(colnames(geneModuleMembership02), function(x) {rep(substring(x,3),nrow(geneModuleMembership02))})))
# 计算性状与基因的相关性矩阵
## 只有连续型性状才能进行计算，如果是离散变量，在构建样品表时就转为0-1矩阵。
if (corType=="pearson") {
  geneTraitCor02 = as.data.frame(cor(tcga_rna_all_mean_wgcna2[(rownames(tcga_rna_all_mean_wgcna2) %in% rownames(annotation_col_wgcna)),], annotation_col_wgcna, use = "p"))
  geneTraitP02 = as.data.frame(corPvalueStudent(
    as.matrix(geneTraitCor02), nrow(annotation_col_wgcna)))
} else {
  geneTraitCorA = bicorAndPvalue(dataExpr, traitData, robustY=robustY)
  geneTraitCor = as.data.frame(geneTraitCorA$bicor)
  geneTraitP   = as.data.frame(geneTraitCorA$p)
}
## Warning in bicor(x, y, use = use, ...): bicor: zero MAD in variable 'y'.
## Pearson correlation was used for individual columns with zero (or missing)
## MAD.
# 最后把两个相关性矩阵联合起来,指定感兴趣模块进行分析
#modules = c("red", "brown", "yellow", "tan")
modules = substring(setdiff(rownames(modTraitCor02[abs(modTraitCor02[,1])>0.5,]), c('MEgrey', 'MEdarkgrey')), 3)
modules2 = c("blue", 'pink', 'brown')
'pink', 'brown', 
hubs02 = c()
modules_genes02 = c()
#module=modules[9]
for (module in modules2) {
  print(module)
  pheno_wgcna = "TIPscore"
  modNames02 = substring(colnames(MEs_col02), 3)
  # 获取关注的列
  module_column02 = match(module, modNames02)
  pheno_column02 = match(pheno_wgcna,colnames(annotation_col_wgcna))
  # 获取模块内的基因
  moduleGenes02 = moduleColors02 == module
  print(sum(moduleGenes02))
  sizeGrWindow(7, 7)
  par(mfrow = c(1,1))
  # 与性状高度相关的基因，也是与性状相关的模型的关键基因
  verboseScatterplot(abs(geneModuleMembership02[moduleGenes02, module_column02]),
                     abs(geneTraitCor02[moduleGenes02, pheno_column02]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = paste("Gene significance for", pheno_wgcna),
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module, abline = T)
  x = 0.5
  y = 0.05
  abline(v=x,lwd=2,col="red")
  abline(h=y,lwd=2,col="red")
  hubs_x = rownames(geneModuleMembership02[(moduleGenes02) & (abs(geneModuleMembership02[,module_column02])>=x),])
  hubs_y = rownames(geneTraitCor02[(moduleGenes02) & (abs(geneTraitCor02[,pheno_wgcna])>=y),])
  hubs_i = intersect(hubs_x, hubs_y)
  hubs02 = union(hubs02, hubs_i)
  modules_genes02 = union(modules_genes02, rownames(geneModuleMembership02[moduleGenes02,]))
}  ##lightyellow-0.36-192个,skyblue3-0.21-111个,saddlebrown--0.17-138个,yellow--0.14-1390个,darkorange-0.14-156个。取lightyellow+skyblue3 + saddlebrown。
write.table(modules_genes02, 'output/2.1_modules_genes02.txt', quote = F, sep = '\t', row.names = F, col.names = F)
write.table(hubs02, 'output/2.1_hubs02.txt', quote = F, sep = '\t', row.names = F, col.names = F)

##全部基因
#单因素回归分析
tcga_rna_tumor_mean2 = as.data.frame(t(tcga_rna_tumor_mean))
tcga_rna_tumor_mean2$'sample' = rownames(tcga_rna_tumor_mean2)
tcga_rna_tumor_mean2 = merge(tcga_rna_tumor_mean2, pheno[,c('sample', 'OS', 'OS.time')], by='sample')  #1362样本*(sample, 256基因, OS, OS.time)
cox_result = data.frame()
for (i in 2:(ncol(tcga_rna_tumor_mean2)-2)) {
  print(i)
  tcga_rna_tumor_mean2$'x' = tcga_rna_tumor_mean2[,i]
  model = coxph(Surv(OS.time, OS) ~ x, data=tcga_rna_tumor_mean2)
  coef = cbind(summary(model)$coefficients, summary(model)$conf.int)
  rownames(coef)[1] = colnames(tcga_rna_tumor_mean2)[i]
  cox_result = rbind(cox_result, coef) #254基因*9列cox属性
}
core_cox_result_01 = cox_result[which(cox_result$`Pr(>|z|)`<0.1),] #80个核心基因*9列cox属性
core_genes01 = rownames(core_cox_result_01)
core_cox_result_005 = cox_result[which(cox_result$`Pr(>|z|)`<0.05),] #80个核心基因*9列cox属性
core_genes005 = rownames(core_cox_result_005)
core_cox_result_001 = cox_result[which(cox_result$`Pr(>|z|)`<0.01),] #80个核心基因*9列cox属性
core_genes001 = rownames(core_cox_result_001)
core_cox_result_0005 = cox_result[which(cox_result$`Pr(>|z|)`<0.005),] #80个核心基因*9列cox属性
core_genes0005 = rownames(core_cox_result_0005)
core_cox_result_0001 = cox_result[which(cox_result$`Pr(>|z|)`<0.001),] #80个核心基因*9列cox属性
core_genes0001 = rownames(core_cox_result_0001)
module_core_genes005 = intersect(core_genes0001, hubs02)

tcga_rna_tumor_mean_lasso = as.data.frame(t(tcga_rna_tumor_mean))
tcga_rna_tumor_mean_lasso$'sample' = rownames(tcga_rna_tumor_mean_lasso)
tcga_rna_tumor_mean_lasso = merge(tcga_rna_tumor_mean_lasso, pheno[,c('sample', 'OS', 'OS.time')], by='sample')
tcga_rna_tumor_mean_lasso = tcga_rna_tumor_mean_lasso[complete.cases(tcga_rna_tumor_mean_lasso),]
geo1_rna_tumor_mean_lasso = as.data.frame(t(geo1_rna_tumor_mean))
geo1_rna_tumor_mean_lasso$'sample' = rownames(geo1_rna_tumor_mean_lasso)
geo1_rna_tumor_mean_lasso = merge(geo1_rna_tumor_mean_lasso, geo1_info[,c('SampleName', 'OS', 'OS.time')], by.x='sample', by.y='SampleName')
geo1_rna_tumor_mean_lasso = geo1_rna_tumor_mean_lasso[complete.cases(geo1_rna_tumor_mean_lasso),]
geo2_rna_tumor_mean_lasso = as.data.frame(t(geo2_rna_tumor_mean))
geo2_rna_tumor_mean_lasso$'sample' = rownames(geo2_rna_tumor_mean_lasso)
geo2_rna_tumor_mean_lasso = merge(geo2_rna_tumor_mean_lasso, geo2_info[,c('sample', 'OS', 'OS.time')], by='sample')
geo2_rna_tumor_mean_lasso = geo2_rna_tumor_mean_lasso[complete.cases(geo2_rna_tumor_mean_lasso),]
geo3_rna_tumor_mean_lasso = as.data.frame(t(geo3_rna_tumor_mean))
geo3_rna_tumor_mean_lasso$'sample' = rownames(geo3_rna_tumor_mean_lasso)
geo3_rna_tumor_mean_lasso = merge(geo3_rna_tumor_mean_lasso, geo3_info[,c('sample', 'OS', 'OS.time')], by='sample')
geo3_rna_tumor_mean_lasso = geo3_rna_tumor_mean_lasso[complete.cases(geo3_rna_tumor_mean_lasso),]
geo4_rna_tumor_mean_lasso = as.data.frame(t(geo4_rna_tumor_mean))
geo4_rna_tumor_mean_lasso$'sample' = rownames(geo4_rna_tumor_mean_lasso)
geo4_rna_tumor_mean_lasso = merge(geo4_rna_tumor_mean_lasso, geo4_info[,c('sample', 'OS', 'OS.time')], by='sample')
geo4_rna_tumor_mean_lasso = geo4_rna_tumor_mean_lasso[complete.cases(geo4_rna_tumor_mean_lasso),]
geo5_rna_tumor_mean_lasso = as.data.frame(t(geo5_rna_tumor_mean))
geo5_rna_tumor_mean_lasso$'sample' = rownames(geo5_rna_tumor_mean_lasso)
geo5_rna_tumor_mean_lasso = merge(geo5_rna_tumor_mean_lasso, geo5_info[,c('sample', 'OS', 'OS.time')], by='sample')
geo5_rna_tumor_mean_lasso = geo5_rna_tumor_mean_lasso[complete.cases(geo5_rna_tumor_mean_lasso),]
geo6_rna_tumor_mean_lasso = as.data.frame(t(geo6_rna_tumor_mean))
geo6_rna_tumor_mean_lasso$'sample' = rownames(geo6_rna_tumor_mean_lasso)
geo6_rna_tumor_mean_lasso = merge(geo6_rna_tumor_mean_lasso, geo6_info[,c('sample', 'OS', 'OS.time')], by='sample')
geo6_rna_tumor_mean_lasso = geo6_rna_tumor_mean_lasso[complete.cases(geo6_rna_tumor_mean_lasso),]
geo7_rna_tumor_mean_lasso = as.data.frame(t(geo7_rna_tumor_mean))
geo7_rna_tumor_mean_lasso$'sample' = rownames(geo7_rna_tumor_mean_lasso)
geo7_rna_tumor_mean_lasso = merge(geo7_rna_tumor_mean_lasso, geo7_info[,c('sample', 'OS', 'OS.time')], by='sample')
geo7_rna_tumor_mean_lasso = geo7_rna_tumor_mean_lasso[complete.cases(geo7_rna_tumor_mean_lasso),]
geo8_rna_tumor_mean_lasso = as.data.frame(t(geo8_rna_tumor_mean))
geo8_rna_tumor_mean_lasso$'sample' = rownames(geo8_rna_tumor_mean_lasso)
geo8_rna_tumor_mean_lasso = merge(geo8_rna_tumor_mean_lasso, geo8_info[,c('sample', 'OS', 'OS.time')], by='sample')
geo8_rna_tumor_mean_lasso = geo8_rna_tumor_mean_lasso[complete.cases(geo8_rna_tumor_mean_lasso),]
geo10_rna_tumor_mean_lasso = as.data.frame(t(geo10_rna_tumor_mean))
geo10_rna_tumor_mean_lasso$'sample' = rownames(geo10_rna_tumor_mean_lasso)
geo10_rna_tumor_mean_lasso = merge(geo10_rna_tumor_mean_lasso, geo10_info[,c('sample', 'OS', 'OS.time')], by='sample')
geo10_rna_tumor_mean_lasso = geo10_rna_tumor_mean_lasso[complete.cases(geo10_rna_tumor_mean_lasso),]
geo11_rna_tumor_mean_lasso = as.data.frame(t(geo11_rna_tumor_mean))
geo11_rna_tumor_mean_lasso$'sample' = rownames(geo11_rna_tumor_mean_lasso)
geo11_rna_tumor_mean_lasso = merge(geo11_rna_tumor_mean_lasso, geo11_info[,c('SampleName', 'OS', 'OS.time')], by.x='sample', by.y='SampleName')
geo11_rna_tumor_mean_lasso = geo11_rna_tumor_mean_lasso[complete.cases(geo11_rna_tumor_mean_lasso),]
geo12_rna_tumor_mean_lasso = as.data.frame(t(geo12_rna_tumor_mean))
geo12_rna_tumor_mean_lasso$'sample' = rownames(geo12_rna_tumor_mean_lasso)
geo12_rna_tumor_mean_lasso = merge(geo12_rna_tumor_mean_lasso, geo12_info[,c('sample', 'OS', 'OS.time')], by='sample')
geo12_rna_tumor_mean_lasso = geo12_rna_tumor_mean_lasso[complete.cases(geo12_rna_tumor_mean_lasso),]

wgcna_x = data.frame('x'=seq(0.5,0.8,0.1))
wgcna_y = data.frame('y'=seq(0.02,0.07,0.01))
uni_cox = data.frame('uni_cox'=c(0.001,0.005,0.01,0.05,0.1))
seeds = data.frame('seed'=seq(1,1))

params_comb0 = merge(wgcna_x, wgcna_y)
params_comb0 = merge(params_comb0, uni_cox)
params_comb0 = merge(params_comb0, seeds)
params_comb0$'n_hubs_lasso' = NA
params_comb0$'tcga_p' = NA
params_comb0$'tcga_if_high' = NA
params_comb0$'tcga_auc' = NA
params_comb0$'geo1_p' = NA
params_comb0$'geo1_if_high' = NA
params_comb0$'geo1_auc' = NA
params_comb0$'geo2_p' = NA
params_comb0$'geo2_if_high' = NA
params_comb0$'geo2_auc' = NA
params_comb0$'geo3_p' = NA
params_comb0$'geo3_if_high' = NA
params_comb0$'geo3_auc' = NA
params_comb0$'geo4_p' = NA
params_comb0$'geo4_if_high' = NA
params_comb0$'geo4_auc' = NA
params_comb0$'geo5_p' = NA
params_comb0$'geo5_if_high' = NA
params_comb0$'geo5_auc' = NA
params_comb0$'geo6_p' = NA
params_comb0$'geo6_if_high' = NA
params_comb0$'geo6_auc' = NA
params_comb0$'geo7_p' = NA
params_comb0$'geo7_if_high' = NA
params_comb0$'geo7_auc' = NA
params_comb0$'geo8_p' = NA
params_comb0$'geo8_if_high' = NA
params_comb0$'geo8_auc' = NA
params_comb0$'geo10_p' = NA
params_comb0$'geo10_if_high' = NA
params_comb0$'geo10_auc' = NA
params_comb0$'geo11_p' = NA
params_comb0$'geo11_if_high' = NA
params_comb0$'geo11_auc' = NA
params_comb0$'geo12_p' = NA
params_comb0$'geo12_if_high' = NA
params_comb0$'geo12_auc' = NA
#params_comb0 = params_comb0[order(params_comb0$cluster_method, params_comb0$cluster_seed, params_comb0$uni_cox),]
#rownames(params_comb0) = seq(1,nrow(params_comb0))

for (param_i in seq(nrow(params_comb0))) {
  #func_main = function(param_i) {
  tryCatch({
    t1 = proc.time()
    print(params_comb0[param_i,])
    #WGCNA卡x,y确定marker
    hubs02 = c()
    modules_genes02 = c()
    for (module in modules2) {
      #print(module)
      pheno_wgcna = "TIPscore"
      modNames02 = substring(colnames(MEs_col02), 3)
      # 获取关注的列
      module_column02 = match(module, modNames02)
      pheno_column02 = match(pheno_wgcna,colnames(annotation_col_wgcna))
      # 获取模块内的基因
      moduleGenes02 = moduleColors02 == module
      #print(sum(moduleGenes02))
      x_i = params_comb0[param_i,'x']
      y_i = params_comb0[param_i,'y']
      hubs_x = rownames(geneModuleMembership02[(moduleGenes02) & (abs(geneModuleMembership02[,module_column02])>=x_i),])
      hubs_y = rownames(geneTraitCor02[(moduleGenes02) & (abs(geneTraitCor02[,pheno_wgcna])>=y_i),])
      hubs_i = intersect(hubs_x, hubs_y)
      hubs02 = union(hubs02, hubs_i)
      modules_genes02 = union(modules_genes02, rownames(geneModuleMembership02[moduleGenes02,]))
    }
    print(length(hubs02))
    #过单因素cox
    if (params_comb0[param_i,'uni_cox']==0.1) {
      module_core_genes005 = intersect(core_genes01, hubs02)
    } else if (params_comb0[param_i,'uni_cox']==0.05) {
      module_core_genes005 = intersect(core_genes005, hubs02)
    } else if (params_comb0[param_i,'uni_cox']==0.01) {
      module_core_genes005 = intersect(core_genes001, hubs02)
    } else if (params_comb0[param_i,'uni_cox']==0.005) {
      module_core_genes005 = intersect(core_genes0005, hubs02)
    } else if (params_comb0[param_i,'uni_cox']==0.001) {
      module_core_genes005 = intersect(core_genes0001, hubs02)
    }
    print(length(module_core_genes005))
    #lasso
    train_lasso = tcga_rna_tumor_mean_lasso
    x = as.matrix(train_lasso[,module_core_genes005])
    y = apply(train_lasso[,(ncol(train_lasso)-1):ncol(train_lasso)], 2, as.numeric)
    colnames(y) = c('status', 'time')
    
    test_lasso_geo1 = geo1_rna_tumor_mean_lasso
    test_lasso_geo1[,setdiff(module_core_genes005, colnames(test_lasso_geo1))] = 0
    x_geo1 = as.matrix(test_lasso_geo1[,module_core_genes005])
    y_geo1 = apply(test_lasso_geo1[,c('OS', 'OS.time')], 2, as.numeric)
    colnames(y_geo1) = c('status', 'time')
    
    test_lasso_geo2 = geo2_rna_tumor_mean_lasso
    test_lasso_geo2[,setdiff(module_core_genes005, colnames(test_lasso_geo2))] = 0
    x_geo2 = as.matrix(test_lasso_geo2[,module_core_genes005])
    y_geo2 = apply(test_lasso_geo2[,c('OS', 'OS.time')], 2, as.numeric)
    colnames(y_geo2) = c('status', 'time')
    
    test_lasso_geo3 = geo3_rna_tumor_mean_lasso
    test_lasso_geo3[,setdiff(module_core_genes005, colnames(test_lasso_geo3))] = 0
    x_geo3 = as.matrix(test_lasso_geo3[,module_core_genes005])
    y_geo3 = apply(test_lasso_geo3[,c('OS', 'OS.time')], 2, as.numeric)
    colnames(y_geo3) = c('status', 'time')
    
    test_lasso_geo4 = geo4_rna_tumor_mean_lasso
    test_lasso_geo4[,setdiff(module_core_genes005, colnames(test_lasso_geo4))] = 0
    x_geo4 = as.matrix(test_lasso_geo4[,module_core_genes005])
    y_geo4 = apply(test_lasso_geo4[,c('OS', 'OS.time')], 2, as.numeric)
    colnames(y_geo4) = c('status', 'time')
    
    test_lasso_geo5 = geo5_rna_tumor_mean_lasso
    test_lasso_geo5[,setdiff(module_core_genes005, colnames(test_lasso_geo5))] = 0
    x_geo5 = as.matrix(test_lasso_geo5[,module_core_genes005])
    y_geo5 = apply(test_lasso_geo5[,c('OS', 'OS.time')], 2, as.numeric)
    colnames(y_geo5) = c('status', 'time')
    
    test_lasso_geo6 = geo6_rna_tumor_mean_lasso
    test_lasso_geo6[,setdiff(module_core_genes005, colnames(test_lasso_geo6))] = 0
    x_geo6 = as.matrix(test_lasso_geo6[,module_core_genes005])
    y_geo6 = apply(test_lasso_geo6[,c('OS', 'OS.time')], 2, as.numeric)
    colnames(y_geo6) = c('status', 'time')
    
    test_lasso_geo7 = geo7_rna_tumor_mean_lasso
    test_lasso_geo7[,setdiff(module_core_genes005, colnames(test_lasso_geo7))] = 0
    x_geo7 = as.matrix(test_lasso_geo7[,module_core_genes005])
    y_geo7 = apply(test_lasso_geo7[,c('OS', 'OS.time')], 2, as.numeric)
    colnames(y_geo7) = c('status', 'time')
    
    test_lasso_geo8 = geo8_rna_tumor_mean_lasso
    test_lasso_geo8[,setdiff(module_core_genes005, colnames(test_lasso_geo8))] = 0
    x_geo8 = as.matrix(test_lasso_geo8[,module_core_genes005])
    x_geo8 = apply(x_geo8, 2, as.numeric)
    y_geo8 = apply(test_lasso_geo8[,c('OS', 'OS.time')], 2, as.numeric)
    colnames(y_geo8) = c('status', 'time')
    
    test_lasso_geo10 = geo10_rna_tumor_mean_lasso
    test_lasso_geo10[,setdiff(module_core_genes005, colnames(test_lasso_geo10))] = 0
    x_geo10 = as.matrix(test_lasso_geo10[,module_core_genes005])
    y_geo10 = apply(test_lasso_geo10[,c('OS', 'OS.time')], 2, as.numeric)
    colnames(y_geo10) = c('status', 'time')
    
    test_lasso_geo11 = geo11_rna_tumor_mean_lasso
    test_lasso_geo11[,setdiff(module_core_genes005, colnames(test_lasso_geo11))] = 0
    x_geo11 = as.matrix(test_lasso_geo11[,module_core_genes005])
    y_geo11 = apply(test_lasso_geo11[,c('OS', 'OS.time')], 2, as.numeric)
    colnames(y_geo11) = c('status', 'time')
    
    test_lasso_geo12 = geo12_rna_tumor_mean_lasso
    test_lasso_geo12[,setdiff(module_core_genes005, colnames(test_lasso_geo12))] = 0
    x_geo12 = as.matrix(test_lasso_geo12[,module_core_genes005])
    y_geo12 = apply(test_lasso_geo12[,c('OS', 'OS.time')], 2, as.numeric)
    colnames(y_geo12) = c('status', 'time')
    t2 = proc.time()
    print(t2-t1)
    
    #拟合模型
    set.seed(params_comb0[param_i,'seed'])
    fit_lasso = glmnet(x, y, family = "cox")
    cvfit = cv.glmnet(x, y, family = "cox")
    coef_min = coef(cvfit, s=cvfit$lambda.min)
    colnames(coef_min) = 'coef_min'
    coef_lasso = as.data.frame(coef_min[which(coef_min!=0),])
    colnames(coef_lasso) = c('coef')
    hubs_lasso = rownames(coef_lasso)
    #params_comb0[param_i, 'n_hubs_lasso'] = length(hubs_lasso)
    t3 = proc.time()
    print(t3-t2)
    if (length(hubs_lasso)>=5 & length(hubs_lasso)<40) {
      riskscore = as.numeric(as.matrix(x[,hubs_lasso]) %*% as.matrix(coef_lasso))
      train_lasso$'TIPRGPI' = riskscore
      train_lasso$'cluster' = ifelse(train_lasso$TIPRGPI>median(train_lasso$TIPRGPI), 'High TIPRGPI', 'Low TIPRGPI')
      surv_diff = survdiff(Surv(OS.time, OS)~cluster, data=train_lasso)
      cluster_p = 1 - pchisq(surv_diff$chisq, length(surv_diff$n) -1)
      roc_riskscore12 = timeROC(T=train_lasso$OS.time,delta=train_lasso$OS,
                                marker=train_lasso$TIPRGPI,cause=1,
                                weighting="marginal",
                                times=c(1*12),ROC=TRUE)
      roc_riskscore24 = timeROC(T=train_lasso$OS.time,delta=train_lasso$OS,
                                marker=train_lasso$TIPRGPI,cause=1,
                                weighting="marginal",
                                times=c(2*12),ROC=TRUE)
      roc_riskscore36 = timeROC(T=train_lasso$OS.time,delta=train_lasso$OS,
                                marker=train_lasso$TIPRGPI,cause=1,
                                weighting="marginal",
                                times=c(3*12),ROC=TRUE)
      #print(roc(train_lasso$OS, train_lasso$TIPRGPI, ci=T)$auc)
      #params_comb0[param_i, 'tcga_p'] = cluster_p
      #params_comb0[param_i, 'tcga_if_high'] = ifelse(median(train_lasso$OS.time[train_lasso$cluster=='High TIPRGPI']) > median(train_lasso$OS.time[train_lasso$cluster=='Low TIPRGPI']), T, F)
      #params_comb0[param_i, 'tcga_auc'] = paste0(c(round(roc_riskscore12$AUC[2],3),
      #                                             round(roc_riskscore24$AUC[2],3),
      #                                             round(roc_riskscore36$AUC[2],3)),
      #                                           collapse = ' ')
      
      run_geo = function(test_lasso_geo, x_geo, geo_p, geo_if_high, geo_auc) {
        riskscore_geo <<- as.numeric(as.matrix(x_geo[,hubs_lasso]) %*% as.matrix(coef_lasso))
        test_lasso_geo_in <<- test_lasso_geo
        test_lasso_geo_in$'TIPRGPI' <<- riskscore_geo
        test_lasso_geo_in$'cluster' <<- ifelse(test_lasso_geo_in$TIPRGPI>median(test_lasso_geo_in$TIPRGPI), 'High TIPRGPI', 'Low TIPRGPI')
        surv_diff <<- survdiff(Surv(OS.time, OS)~cluster, data=test_lasso_geo_in)
        cluster_p <<- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) -1)
        roc_riskscore12_geo <<- timeROC(T=test_lasso_geo_in$OS.time,delta=test_lasso_geo_in$OS,
                                      marker=test_lasso_geo_in$TIPRGPI,cause=1,
                                      weighting="marginal",
                                      times=c(1*12),ROC=TRUE)
        roc_riskscore24_geo <<- timeROC(T=test_lasso_geo_in$OS.time,delta=test_lasso_geo_in$OS,
                                      marker=test_lasso_geo_in$TIPRGPI,cause=1,
                                      weighting="marginal",
                                      times=c(2*12),ROC=TRUE)
        roc_riskscore36_geo <<- timeROC(T=test_lasso_geo_in$OS.time,delta=test_lasso_geo_in$OS,
                                      marker=test_lasso_geo_in$TIPRGPI,cause=1,
                                      weighting="marginal",
                                      times=c(3*12),ROC=TRUE)
        #params_comb0[param_i, geo_p] <<- cluster_p
        #params_comb0[param_i, geo_if_high] <<- ifelse(median(test_lasso_geo_in$OS.time[test_lasso_geo_in$cluster=='High TIPRGPI']) > median(test_lasso_geo_in$OS.time[test_lasso_geo_in$cluster=='Low TIPRGPI']), T, F)
        #params_comb0[param_i, geo_auc] <<- paste0(c(round(roc_riskscore12_geo$AUC[2],3),
        #                                            round(roc_riskscore24_geo$AUC[2],3),
        #                                            round(roc_riskscore36_geo$AUC[2],3)),
        #                                          collapse = ' ')
      }
      #*geo1
      try(run_geo(test_lasso_geo1, x_geo1, 'geo1_p', 'geo1_if_high', 'geo1_auc'), silent = T)
      #*geo2
      try(run_geo(test_lasso_geo2, x_geo2, 'geo2_p', 'geo2_if_high', 'geo2_auc'), silent = T)
      #*geo3
      try(run_geo(test_lasso_geo3, x_geo3, 'geo3_p', 'geo3_if_high', 'geo3_auc'), silent = T)
      #*geo4
      try(run_geo(test_lasso_geo4, x_geo4, 'geo4_p', 'geo4_if_high', 'geo4_auc'), silent = T)
      #*geo5
      try(run_geo(test_lasso_geo5, x_geo5, 'geo5_p', 'geo5_if_high', 'geo5_auc'), silent = T)
      #*geo6
      try(run_geo(test_lasso_geo6, x_geo6, 'geo6_p', 'geo6_if_high', 'geo6_auc'), silent = T)
      #*geo7
      try(run_geo(test_lasso_geo7, x_geo7, 'geo7_p', 'geo7_if_high', 'geo7_auc'), silent = T)
      #*geo8
      try(run_geo(test_lasso_geo8, x_geo8, 'geo8_p', 'geo8_if_high', 'geo8_auc'), silent = T)
      #*geo10
      try(run_geo(test_lasso_geo10, x_geo10, 'geo10_p', 'geo10_if_high', 'geo10_auc'), silent = T)
      #*geo11
      try(run_geo(test_lasso_geo11, x_geo11, 'geo11_p', 'geo11_if_high', 'geo11_auc'), silent = T)
      #*geo12
      try(run_geo(test_lasso_geo12, x_geo12, 'geo12_p', 'geo12_if_high', 'geo12_auc'), silent = T)
      print(params_comb0[param_i,])
      t4 = proc.time()
      print(t4-t3)
    }  
  }, error = function(e){cat(param_i, conditionMessage(e),"\n")})
}
params_comb02 = params_comb0[!is.na(params_comb0$tcga_p) & (params_comb0$tcga_p<0.05),]
params_comb02[is.na(params_comb02)] = 99
params_comb02 = params_comb02[((params_comb02$geo1_p<0.05) + 
                                 (params_comb02$geo2_p<0.05) + 
                                 (params_comb02$geo3_p<0.05) + 
                                 (params_comb02$geo4_p<0.05) + 
                                 (params_comb02$geo5_p<0.05) +
                                 (params_comb02$geo6_p<0.05) + 
                                 (params_comb02$geo7_p<0.05) + 
                                 (params_comb02$geo8_p<0.05) + 
                                 (params_comb02$geo10_p<0.05) + 
                                 (params_comb02$geo11_p<0.05) + 
                                 (params_comb02$geo12_p<0.05))>=2 &
                                ((params_comb02$geo1_p<0.05) + 
                                   (params_comb02$geo2_p<0.05))>=1,]
#params_comb0_blue = params_comb0 #取param_i=98,geo3 & geo7
#params_comb02_blue = params_comb02
params_comb0_3modules = params_comb0 #取param_i=85,geo3
params_comb02_3modules = params_comb02
write.table(hubs02, 'output/2.1_hubs02.txt', sep = '\t', quote = F, row.names = F, col.names = F)
write.table(module_core_genes005, 'output/2.3_module_core_genes005.txt', sep = '\t', quote = F, row.names = F, col.names = F)
write.table(hubs_lasso, 'output/3.2_hubs_lasso.txt', sep = '\t', quote = F, row.names = F, col.names = F)

#取param_i=85,geo2 & geo3
ann_colors = c(ann_colors, list(cluster=c('High TIPRGPI'='red', 'Low TIPRGPI'='blue')))
plot_surv(train_lasso, c('High TIPRGPI','Low TIPRGPI'), 2, '')
run_geo(test_lasso_geo2, x_geo2, 'geo2_p', 'geo2_if_high', 'geo2_auc')
plot_surv(test_lasso_geo_in, c('High TIPRGPI','Low TIPRGPI'), 2, '')
run_geo(test_lasso_geo3, x_geo3, 'geo3_p', 'geo3_if_high', 'geo3_auc')
plot_surv(test_lasso_geo_in, c('High TIPRGPI','Low TIPRGPI'), 2, '')

##Section7：TIPRGPI可作为免疫治疗敏感预测因子
##2、基于免疫治疗数据集TIPRGPI评估免疫治疗预后（图2：若有癌型本身免疫治疗数据集则做图（I-K没有随访信息省略）若无则作图L-O（L-N用其它癌型免疫治疗数据集分析；O为TIDE）；若J、M不显著则分三组：CR/PR vs PD vs SD 或者CR/PR/SD vs PD尝试）
#*geo2--不缺hubs_lasso基因--用这个
run_geo(test_lasso_geo2, x_geo2, 'geo2_p', 'geo2_if_high', 'geo2_auc')
test_lasso_geo2_in = test_lasso_geo_in
test_lasso_geo2_in = merge(test_lasso_geo2_in, geo2_info %>% select(., -c('OS', 'OS.time')), by='sample')
test_lasso_geo2_in2 = test_lasso_geo2_in[,c('TIPRGPI', 'Response')] %>% filter(Response!='NE')
kruskal.test(test_lasso_geo2_in2$TIPRGPI, test_lasso_geo2_in2$Response)  #p=0.09859
wilcox.test(test_lasso_geo2_in2$TIPRGPI[test_lasso_geo2_in2$Response=='CR'],
            test_lasso_geo2_in2$TIPRGPI[test_lasso_geo2_in2$Response=='PD'])  #p=0.6727
wilcox.test(test_lasso_geo2_in2$TIPRGPI[test_lasso_geo2_in2$Response=='PR'],
            test_lasso_geo2_in2$TIPRGPI[test_lasso_geo2_in2$Response=='PD'])  #p=0.01407
wilcox.test(test_lasso_geo2_in2$TIPRGPI[test_lasso_geo2_in2$Response%in%c('PR', 'CR', 'SD')],
            test_lasso_geo2_in2$TIPRGPI[test_lasso_geo2_in2$Response=='PD'])  #p=0.03605
test_lasso_geo2_in2$response = Recode(test_lasso_geo2_in2$Response, "c('PR', 'CR', 'SD')='CR/PR/SD'")
cps = t(combn(unique(test_lasso_geo2_in2$response), 2))
cps = split(cps,list(cps[,1],cps[,2]), drop = T)
ggplot(test_lasso_geo2_in2, aes(x=response,y=TIPRGPI, fill=response))+
  stat_boxplot(position = position_dodge()) +
  geom_boxplot(size=0.5, outlier.color = 'white')+
  theme_bw()+
  theme(#legend.position="none",
    axis.text.x=element_text(colour="black",family="Times",size=10),
    axis.text.y=element_text(family="Times",size=10,face="plain"),
    axis.title.y=element_text(family="Times",size = 14,face="plain"),
    axis.title.x=element_text(family="Times",size = 14,face="plain"),
    axis.line.x = element_line(),
    axis.line.y = element_line(),
    plot.title = element_text(family="Times",size=15,face="bold",hjust = 0.5),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+
  xlab("")+ylab("TIPRGPI") +
  scale_fill_manual(values = c("firebrick", "forestgreen", 'brown', 'blue')) +
  stat_compare_means(comparisons = cps, label='p.format', method = 'wilcox.test')
#roc_geo2 = pROC::roc(test_lasso_geo2_in$OS, test_lasso_geo2_in$TIPRGPI, ci=T)
roc_riskscore12_geo2 = roc_riskscore12_geo
roc_riskscore24_geo2 = roc_riskscore24_geo
roc_riskscore36_geo2 = roc_riskscore36_geo
plot(roc_riskscore12_geo2,time=1*12,col="red")
plot(roc_riskscore24_geo2,time=2*12,add=T,col="blue")
plot(roc_riskscore36_geo2,time=3*12,add=T,col="orange")
legend("bottomright",c(paste0("1-year survival:",round(roc_riskscore12_geo2$AUC[2],3)),paste0("2-year survival:",round(roc_riskscore24_geo2$AUC[2],3)),paste0("3-year survival:",round(roc_riskscore36_geo2$AUC[2],3))),col=c("red","blue","orange"),lty=1,lwd=2)


##其余步骤
##Section1：肿瘤免疫表型得分关联临床
##2.K-M曲线分析TIPscore与研究癌症的预后关系
tips_tcga2_dfs = merge(tips_tcga, pheno_cgds2[,c('sample', 'DFS', 'DFS.time')], by='sample')
surv_fit2 <<- survfit(Surv(DFS.time, DFS)~cluster, data=tips_tcga2_dfs %>% dplyr::rename('cluster'='TIP_cluster'))
ggsurvplot(surv_fit2, pval=T, 
           ggtheme = theme(axis.line = element_line(linetype = "solid"), 
                           #axis.text = element_text(hjust = 1), 
                           #axis.text.x = element_text(size = 11, vjust = 0.5, angle = 90),
                           panel.background = element_rect(fill = NA)
           ),
           risk.table=T, 
           legend.labs=c('High TIPscore', 'Low TIPscore'),
           #palette = c('red', 'blue'),
           palette = ann_colors[[1]],
           legend.title=NULL,
           #legend.position = c(0.05,0.15),
           xlab='Time in months'
)

##3.分析TIPscore与免疫相关指标显著（ESTIMATE、免疫治marker: activated CD4/CD8 cell 、PD-1/CTLA-4检查点）
#estimate
write.table(tcga_rna_tumor_mean, 'output/1.3_tcga_rna_tumor_mean.txt', sep = '\t', col.names = T, row.names = T, quote = F)
filterCommonGenes(input.f = 'output/1.3_tcga_rna_tumor_mean.txt', output.f = 'output/1.3_tcga_rna_tumor_mean.gct', id = 'GeneSymbol')
estimateScore('output/1.3_tcga_rna_tumor_mean.gct', 'output/1.3_tcga_rna_tumor_score.gct')
estimate_score = read.table('output/1.3_tcga_rna_tumor_score.gct', skip = 2, header = T, row.names = 1)
estimate_score = estimate_score[,-1]
colnames(estimate_score) = str_replace_all(colnames(estimate_score), '\\.', '-')
estimate_score = as.data.frame(t(estimate_score))
estimate_score$'sample' = rownames(estimate_score)
#xCell评分
xcell_tcga = t(xCellAnalysis(tcga_rna_tumor_mean))
tip_tcga3 = merge(tips_tcga2, estimate_score, by='sample')
tip_tcga3 = merge(tip_tcga3, as.data.frame(xcell_tcga) %>% mutate('sample'=rownames(.)) %>% select(c('sample', "CD4+ T-cells", "CD8+ T-cells")), by='sample')
tip_tcga3 = merge(tip_tcga3, as.data.frame(t(tcga_rna_tumor_mean[c('PDCD1', 'CTLA4'),])) %>% mutate('sample'=rownames(.)), by='sample')
ggscatterhist(tip_tcga3, x = "TIPscore", y = 'StromalScore',
              color = "black", size = 3, # 点的颜色与大小
              add = "reg.line",  # 添加回归线
              add.params = list(color = "red", fill = "lightgray"), # 回归线的颜色设置为红色，区间颜色设置为灰色
              conf.int = TRUE, # 添加回归线的置信区间
              cor.coef = TRUE, # 添加相关系数
              cor.coeff.args = list(method = "spearman", label.sep = "\t", label.y=2200),#选择Pearson相关,
              margin.params = list(fill = "blue3"),
              ggtheme = theme_bw(), xlab = 'TIPscore'
)
ggscatterhist(tip_tcga3, x = "TIPscore", y = 'ImmuneScore',
              color = "black", size = 3, # 点的颜色与大小
              add = "reg.line",  # 添加回归线
              add.params = list(color = "red", fill = "lightgray"), # 回归线的颜色设置为红色，区间颜色设置为灰色
              conf.int = TRUE, # 添加回归线的置信区间
              cor.coef = TRUE, # 添加相关系数
              cor.coeff.args = list(method = "spearman", label.sep = "\t", label.y=4200),#选择Pearson相关,
              margin.params = list(fill = "blue3"),
              ggtheme = theme_bw(), xlab = 'TIPscore'
)
ggscatterhist(tip_tcga3, x = "TIPscore", y = 'ESTIMATEScore',
              color = "black", size = 3, # 点的颜色与大小
              add = "reg.line",  # 添加回归线
              add.params = list(color = "red", fill = "lightgray"), # 回归线的颜色设置为红色，区间颜色设置为灰色
              conf.int = TRUE, # 添加回归线的置信区间
              cor.coef = TRUE, # 添加相关系数
              cor.coeff.args = list(method = "spearman", label.sep = "\t", label.y=5200),#选择Pearson相关,
              margin.params = list(fill = "blue3"),
              ggtheme = theme_bw(), xlab = 'TIPscore'
)
ggscatterhist(tip_tcga3, x = "TIPscore", y = 'TumorPurity',
              color = "black", size = 3, # 点的颜色与大小
              add = "reg.line",  # 添加回归线
              add.params = list(color = "red", fill = "lightgray"), # 回归线的颜色设置为红色，区间颜色设置为灰色
              conf.int = TRUE, # 添加回归线的置信区间
              cor.coef = TRUE, # 添加相关系数
              cor.coeff.args = list(method = "spearman", label.sep = "\t", label.y=1.02),#选择Pearson相关,
              margin.params = list(fill = "blue3"),
              ggtheme = theme_bw(), xlab = 'TIPscore'
)
ggscatterhist(tip_tcga3, x = "TIPscore", y = 'CD4+ T-cells',
              color = "black", size = 3, # 点的颜色与大小
              add = "reg.line",  # 添加回归线
              add.params = list(color = "red", fill = "lightgray"), # 回归线的颜色设置为红色，区间颜色设置为灰色
              conf.int = TRUE, # 添加回归线的置信区间
              cor.coef = TRUE, # 添加相关系数
              cor.coeff.args = list(method = "spearman", label.sep = "\t"),#选择Pearson相关,
              margin.params = list(fill = "blue3"),
              ggtheme = theme_bw(), xlab = 'TIPscore'
)
ggscatterhist(tip_tcga3, x = "TIPscore", y = 'CD8+ T-cells',
              color = "black", size = 3, # 点的颜色与大小
              add = "reg.line",  # 添加回归线
              add.params = list(color = "red", fill = "lightgray"), # 回归线的颜色设置为红色，区间颜色设置为灰色
              conf.int = TRUE, # 添加回归线的置信区间
              cor.coef = TRUE, # 添加相关系数
              cor.coeff.args = list(method = "spearman", label.sep = "\t"),#选择Pearson相关,
              margin.params = list(fill = "blue3"),
              ggtheme = theme_bw(), xlab = 'TIPscore'
)
ggscatterhist(tip_tcga3, x = "TIPscore", y = 'PDCD1',
              color = "black", size = 3, # 点的颜色与大小
              add = "reg.line",  # 添加回归线
              add.params = list(color = "red", fill = "lightgray"), # 回归线的颜色设置为红色，区间颜色设置为灰色
              conf.int = TRUE, # 添加回归线的置信区间
              cor.coef = TRUE, # 添加相关系数
              cor.coeff.args = list(method = "spearman", label.sep = "\t"),#选择Pearson相关,
              margin.params = list(fill = "blue3"),
              ggtheme = theme_bw(), xlab = 'TIPscore', ylab = 'PD-1'
)
ggscatterhist(tip_tcga3, x = "TIPscore", y = 'CTLA4',
              color = "black", size = 3, # 点的颜色与大小
              add = "reg.line",  # 添加回归线
              add.params = list(color = "red", fill = "lightgray"), # 回归线的颜色设置为红色，区间颜色设置为灰色
              conf.int = TRUE, # 添加回归线的置信区间
              cor.coef = TRUE, # 添加相关系数
              cor.coeff.args = list(method = "spearman", label.sep = "\t", label.y=6.6),#选择Pearson相关,
              margin.params = list(fill = "blue3"),
              ggtheme = theme_bw(), xlab = 'TIPscore'
)

##4、分析TIPscore与其它免疫治疗signature相关性（TMB、TIDE、IPS若不显著则不展示）。
tip_tcga4 = merge(tip_tcga3, tmb[,c('sample', 'Non-silent per Mb')], by='sample')
ggscatterhist(tip_tcga4, x = "TIPscore", y = 'Non-silent per Mb',
              color = "black", size = 3, # 点的颜色与大小
              add = "reg.line",  # 添加回归线
              add.params = list(color = "red", fill = "lightgray"), # 回归线的颜色设置为红色，区间颜色设置为灰色
              conf.int = TRUE, # 添加回归线的置信区间
              cor.coef = TRUE, # 添加相关系数
              cor.coeff.args = list(method = "spearman", label.sep = "\t"),#选择Pearson相关,
              margin.params = list(fill = "blue3"),
              ggtheme = theme_bw(), xlab = 'TIPscore'
)
#*TIDE
write.table(sweep(tcga_rna_tumor_mean, 1, apply(tcga_rna_tumor_mean, 1, median)), file='output/1.4_tcga_rna_tumor_mean.txt', row.names=T, col.names=T, quote=F, sep='\t')
tide = read.csv('output/1.4_tide_result.csv')
tip_tcga4 = merge(tip_tcga4,tide[,c('Patient', 'TIDE')], by.y='Patient', by.x='sample')
ggscatterhist(tip_tcga4, x = "TIPscore", y = 'TIDE',
              color = "black", size = 3, # 点的颜色与大小
              add = "reg.line",  # 添加回归线
              add.params = list(color = "red", fill = "lightgray"), # 回归线的颜色设置为红色，区间颜色设置为灰色
              conf.int = TRUE, # 添加回归线的置信区间
              cor.coef = TRUE, # 添加相关系数
              cor.coeff.args = list(method = "spearman", label.sep = "\t", label.y=3.5),#选择Pearson相关,
              margin.params = list(fill = "blue3"),
              ggtheme = theme_bw(), xlab = 'TIPscore'
)

##Section2：肿瘤表型相关预后分层系统（TIPRGPI）的构建
##2、对关键模块进行GO、KEGG富集分析（若多个模块先合并再分析）；
eg = bitr(hubs02, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL",'SYMBOL'), OrgDb="org.Hs.eg.db")
ego_bp = enrichGO(gene=eg$ENTREZID, OrgDb=org.Hs.eg.db, ont='BP', pAdjustMethod="BH", pvalueCutoff=0.05, qvalueCutoff = 0.1, readable=T)
ego_bp2 = data.frame(ego_bp)
ego_bp2 = ego_bp2[order(ego_bp2$p.adjust),][1:10,]
ego_cc = enrichGO(gene=eg$ENTREZID, OrgDb=org.Hs.eg.db, ont='CC', pAdjustMethod="BH", pvalueCutoff=0.05, qvalueCutoff = 0.1, readable=T)
ego_cc2 = data.frame(ego_cc)
ego_cc2 = ego_cc2[order(ego_cc2$p.adjust),][1:10,]
ego_mf = enrichGO(gene=eg$ENTREZID, OrgDb=org.Hs.eg.db, ont='MF', pAdjustMethod="BH", pvalueCutoff=0.05, qvalueCutoff = 0.1, readable=T)
ego_mf2 = data.frame(ego_mf)
ego_mf2 = ego_mf2[order(ego_mf2$p.adjust),][1:10,]
kegg = enrichKEGG(gene=eg$ENTREZID, organism = 'hsa', pAdjustMethod="BH", pvalueCutoff=0.05, qvalueCutoff = 0.1, use_internal_data = T)
kegg2 = data.frame(kegg)
kegg2 = kegg2[order(kegg2$p.adjust),][1:10,]
df_enrich = rbind(ego_bp2, ego_cc2)
df_enrich = rbind(df_enrich, ego_mf2)
df_enrich$'type' = c(rep('BP', 10), rep('CC', 10), rep('MF', 10))
df_enrich = df_enrich %>% arrange(type, p.adjust) %>% mutate(order=seq(nrow(.)))
ggplot(df_enrich, aes(x=fct_reorder(Description, order), y=-log10(p.adjust), fill=type)) +
  geom_bar(stat='identity') +
  #coord_flip() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        plot.title = element_text(hjust=0.5)) +
  labs(x=NULL, y='-log10(FDR)', title = 'The Most Enriched GO Terms') +
  guides(fill=guide_legend(title = 'Category'))
enrichplot::dotplot(kegg, showCategory=20, title=NULL, color = 'p.adjust', label_format=60)

##3、基于关键模块基因单因素cox回归分析筛选预后因子，lasso或者多因素cox回归分析进一步筛选构建TIPRGPI的最优基因组合，进行K-M曲线展示（2、3图可以组图；K-M曲线展示模型基因）。
#lasso图
plot(cvfit)
plot(fit_lasso, xvar='lambda')
##森林图
ggplot(coef_lasso %>% mutate('sample'=rownames(.)) %>% mutate(color=ifelse(coef>0, 'pos', 'neg')),
       aes(x=coef, y=reorder(sample, coef, .desc=T), fill=color)) +
  geom_bar(stat='identity', color='black') +
  theme(panel.background = element_blank(),
        axis.line.x = element_line(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'none') +
  labs(x=NULL, y=NULL, title = 'Coefficients') +
  scale_fill_manual(values = c('pos'='orange', 'neg'='dodgerblue4')) +
  scale_color_manual(values = 'black') +
  scale_x_continuous(expand = c(0,0),
                     limits = c(-0.3000001,0.2))
for (gene in hubs_lasso) {
  print(gene)
  cell_lasso = tcga_rna_tumor_mean2[,c(gene, 'OS', 'OS.time')]
  cell_lasso$'cluster' = ifelse(cell_lasso[,gene]>median(cell_lasso[,gene]), 'High', 'Low')
  
  pdf(paste0("output/2.3/2.3_survival_", gene, ".pdf"), width=8, height=5, onefile = F)
  print(plot_surv(cell_lasso, c('highExp', 'lowExp'), 1, gene))
  dev.off()
  
  tiff(paste0("output/2.3/2.3_survival_", gene, ".tiff"), res=300, width=8, height=5, compression="lzw", units="in")
  print(plot_surv(cell_lasso, c('highExp', 'lowExp'), 1, gene))
  dev.off()
}

##Section3：肿瘤表型相关预后分层系统（TIPRGPI）效能评估和独立数据集验证
##1、训练集及独立外部数据集分析TIPRGPI的预测效能（图1）；
#*TCGA
#散点图
score_point = train_lasso
score_point = score_point[order(score_point$TIPRGPI),]
score_point$number = seq(dim(score_point)[1])
ggplot(data = score_point) +
  geom_point(aes(x = number, y = TIPRGPI, color = cluster)) +
  theme(axis.line = element_line(linetype = "solid"), 
        panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"), 
        #legend.position = c(0.08, 0.85),
        panel.background = element_blank()) +
  labs(x = "Patients(incresing risk score)", y = "Risk score", color = NULL) +
  geom_vline(xintercept = max(score_point$number[score_point$cluster=='Low TIPRGPI']), linetype='dashed') +
  geom_hline(yintercept = max(score_point$TIPRGPI[score_point$cluster=='Low TIPRGPI']), linetype='dashed') +
  scale_color_manual(values = c("red", "blue")) +
  geom_vline(xintercept = sum(score_point$cluster=='Low TIPRGPI'), lty=2) +
  geom_hline(yintercept = median(score_point$TIPRGPI), lty=2)
ggplot(data = score_point) +
  geom_point(aes(x = number, y = OS.time, color = Recode(score_point$OS, "0='Alive';1='Dead'"))) +
  theme(axis.line = element_line(linetype = "solid"), 
        panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"), 
        #legend.position = c(0.06, 0.85),
        panel.background = element_blank()) +
  labs(x = "Patients(incresing risk score)", y = "Survival time (month)", color = NULL) +
  geom_vline(xintercept = max(score_point$number[score_point$cluster=='Low TIPRGPI']), linetype='dashed') +
  scale_color_manual(values = c("red", "blue")) +
  geom_vline(xintercept = sum(score_point$cluster=='Low TIPRGPI'), lty=2)
#热图
ann_colors = c(ann_colors, list(TIPRGPI=c('High TIPRGPI'='red', 'Low TIPRGPI'='blue')))
expr_core_genes = t(tcga_rna_tumor_mean[(rownames(tcga_rna_tumor_mean) %in% hubs_lasso),])
annotation_col_tcga = data.frame('TIPRGPI'=train_lasso$cluster,
                                 'sample'=train_lasso$sample)
rownames(annotation_col_tcga) = train_lasso$sample
annotation_col_tcga = annotation_col_tcga[order(annotation_col_tcga$TIPRGPI),]
pheatmap(t(expr_core_genes[annotation_col_tcga$sample,]), scale = 'row', cluster_row = T, cluster_col = F, show_colnames = F, annotation_col=annotation_col_tcga[c('TIPRGPI')], color = colorRampPalette(colors = c("royalblue2","white","red2"))(100), show_rownames=T, annotation_colors = ann_colors[3], border_color = NA, breaks = unique(c(seq(-2, 2, length=100))))
wilcox_results_tcga = data.frame()
for (gene in colnames(expr_core_genes)) {
  print(gene)
  wilcox_result = wilcox.test(as.numeric(expr_core_genes[annotation_col_tcga$sample[annotation_col_tcga$TIPRGPI=='High TIPRGPI'], gene]), 
                              as.numeric(expr_core_genes[annotation_col_tcga$sample[annotation_col_tcga$TIPRGPI=='Low TIPRGPI'], gene]))
  wilcox_results_tcga = rbind(wilcox_results_tcga,
                              data.frame('statistic'=wilcox_result$statistic,
                                         'pvalue'=wilcox_result$p.value,
                                         'gene'=gene))
}
rownames(wilcox_results_tcga) = wilcox_results_tcga$gene
write.table(wilcox_results_tcga[,c(3,1,2)], 'output/3.1_wilcox_results_tcga.txt', col.names = T, row.names = F, sep = '\t', quote = F)
#ROC
plot(roc_riskscore12,time=1*12,col="red")
plot(roc_riskscore24,time=2*12,add=T,col="blue")
plot(roc_riskscore36,time=3*12,add=T,col="orange")
legend("bottomright",c(paste0("1-year survival:",round(roc_riskscore12$AUC[2],3)),paste0("2-year survival:",round(roc_riskscore24$AUC[2],3)),paste0("3-year survival:",round(roc_riskscore36$AUC[2],3))),col=c("red","blue","orange"),lty=1,lwd=2)
#小提琴图
cps2 = t(combn(unique(score_point$cluster), 2))
cps2 = split(cps2,list(cps2[,1],cps2[,2]), drop = T)
ggviolin(score_point %>% rename('TIPRGPI subtype'=cluster), 'TIPRGPI subtype', 'TIPRGPI', color = 'black', fill = 'TIPRGPI subtype',
         palette = c("blue", "red"),
         add = c('boxplot'), add.params = list(fill='white')) +
  theme(axis.title.x = element_blank(),
        legend.position = 'none') +
  labs(y='TIPRGPI') +
  stat_compare_means(comparisons = cps2, label = 'p.format', label.y = 1.2)
#堆叠条形图
score31_bar = score_point[, c('cluster', 'OS')] %>% mutate(OS = replace(OS, OS==0, 'Alive')) %>% mutate(OS = replace(OS, OS==1, 'Dead'))
score31_bar = score31_bar[complete.cases(score31_bar),]
score31_bar$num = 1
score31_bar2 = aggregate(score31_bar$num,by=list(Variety=score31_bar$cluster, score31_bar$OS), length)
score31_bar2 = score31_bar2 %>% group_by(Variety) %>% mutate(freq=round(x/sum(x)*100, 1))
chisq_bar = chisq.test(matrix(score31_bar2$x, ncol = 2))
ggplot(score31_bar2, aes(x=Variety, y=freq, fill=Group.2)) + 
  geom_histogram(position="fill", stat='identity', width=0.8)+ 
  theme(panel.background = element_blank(),
        axis.line = element_line(),
        plot.title = element_text(hjust = 0.5)) + 
  geom_text(aes(label = paste0(as.character(freq), '%')), size = 5, position = position_fill(vjust = 0.5)) +
  scale_fill_manual(values = c('dodgerblue', 'firebrick1')) +
  ggtitle('TCGA')+
  labs(x = NULL, y = "Relative Percent(x100%)") +
  guides(fill=guide_legend(title = NULL)) +
  annotate(geom = 'text', label='P<0.001', x = 1.5, y = 1.05)

#*geo3
#散点图
score_geo_point = test_lasso_geo_in
score_geo_point = score_geo_point[order(score_geo_point$TIPRGPI),]
score_geo_point$number = seq(dim(score_geo_point)[1])
ggplot(data = score_geo_point) +
  geom_point(aes(x = number, y = TIPRGPI, color = cluster)) +
  theme(axis.line = element_line(linetype = "solid"), 
        panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"), 
        #legend.position = c(0.08, 0.85),
        panel.background = element_blank()) +
  labs(x = "Patients(incresing risk score)", y = "Risk score", color = NULL) +
  geom_vline(xintercept = max(score_geo_point$number[score_geo_point$cluster=='Low TIPRGPI']), linetype='dashed') +
  geom_hline(yintercept = max(score_geo_point$TIPRGPI[score_geo_point$cluster=='Low TIPRGPI']), linetype='dashed') +
  scale_color_manual(values = c("red", "blue")) +
  geom_vline(xintercept = sum(score_geo_point$cluster=='Low TIPRGPI'), lty=2) +
  geom_hline(yintercept = median(score_geo_point$TIPRGPI), lty=2)
ggplot(data = score_geo_point) +
  geom_point(aes(x = number, y = OS.time, color = Recode(score_geo_point$OS, "0='Alive';1='Dead'"))) +
  theme(axis.line = element_line(linetype = "solid"), 
        panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"), 
        #legend.position = c(0.06, 0.85),
        panel.background = element_blank()) +
  labs(x = "Patients(incresing risk score)", y = "Survival time (month)", color = NULL) +
  geom_vline(xintercept = max(score_geo_point$number[score_geo_point$cluster=='Low TIPRGPI']), linetype='dashed') +
  scale_color_manual(values = c("red", "blue")) +
  geom_vline(xintercept = sum(score_geo_point$cluster=='Low TIPRGPI'), lty=2)
#热图
expr_core_genes_geo = test_lasso_geo_in[,(colnames(test_lasso_geo_in) %in% hubs_lasso)]
rownames(expr_core_genes_geo) = test_lasso_geo_in$sample
annotation_col_geo = data.frame('TIPRGPI'=test_lasso_geo_in$cluster,
                                 'sample'=test_lasso_geo_in$sample)
rownames(annotation_col_geo) = test_lasso_geo_in$sample
annotation_col_geo = annotation_col_geo[order(annotation_col_geo$TIPRGPI),]
pheatmap(t(expr_core_genes_geo[annotation_col_geo$sample,]), scale = 'row', cluster_row = T, cluster_col = F, show_colnames = F, annotation_col=annotation_col_geo[c('TIPRGPI')], color = colorRampPalette(colors = c("royalblue2","white","red2"))(100), show_rownames=T, annotation_colors = ann_colors[3], border_color = NA, breaks = unique(c(seq(-2,2, length=100))))
wilcox_results_geo = data.frame()
for (gene in colnames(expr_core_genes_geo)) {
  print(gene)
  wilcox_result = wilcox.test(as.numeric(expr_core_genes_geo[annotation_col_geo$sample[annotation_col_geo$TIPRGPI=='High TIPRGPI'], gene]), 
                              as.numeric(expr_core_genes_geo[annotation_col_geo$sample[annotation_col_geo$TIPRGPI=='Low TIPRGPI'], gene]))
  wilcox_results_geo = rbind(wilcox_results_geo,
                              data.frame('statistic'=wilcox_result$statistic,
                                         'pvalue'=wilcox_result$p.value,
                                         'gene'=gene))
}
rownames(wilcox_results_geo) = wilcox_results_geo$gene
write.table(wilcox_results_geo[,c(3,1,2)], 'output/4.2_wilcox_results_GSE65904.txt', col.names = T, row.names = F, sep = '\t', quote = F)
#ROC
plot(roc_riskscore12_geo,time=1*12,col="red")
plot(roc_riskscore24_geo,time=2*12,add=T,col="blue")
plot(roc_riskscore36_geo,time=3*12,add=T,col="orange")
legend("bottomright",c(paste0("1-year survival:",round(roc_riskscore12_geo$AUC[2],3)),paste0("2-year survival:",round(roc_riskscore24_geo$AUC[2],3)),paste0("3-year survival:",round(roc_riskscore36_geo$AUC[2],3))),col=c("red","blue","orange"),lty=1,lwd=2)
#小提琴图
cps2_geo = t(combn(unique(score_geo_point$cluster), 2))
cps2_geo = split(cps2_geo,list(cps2_geo[,1],cps2_geo[,2]), drop = T)
ggviolin(score_geo_point %>% rename('TIPRGPI subtype'=cluster), 'TIPRGPI subtype', 'TIPRGPI', color = 'black', fill = 'TIPRGPI subtype',
         palette = c("blue", "red"),
         add = c('boxplot'), add.params = list(fill='white')) +
  theme(axis.title.x = element_blank(),
        legend.position = 'none') +
  labs(y='TIPRGPI') +
  stat_compare_means(comparisons = cps2_geo, label = 'p.format', label.y = -8)
#堆叠条形图
score31_geo_bar = score_geo_point[, c('cluster', 'OS')] %>% mutate(OS = replace(OS, OS==0, 'Alive')) %>% mutate(OS = replace(OS, OS==1, 'Dead'))
score31_geo_bar = score31_geo_bar[complete.cases(score31_geo_bar),]
score31_geo_bar$num = 1
score31_geo_bar2 = aggregate(score31_geo_bar$num,by=list(Variety=score31_geo_bar$cluster, score31_geo_bar$OS), length)
score31_geo_bar2 = score31_geo_bar2 %>% group_by(Variety) %>% mutate(freq=round(x/sum(x)*100, 1))
chisq_bar_geo = chisq.test(matrix(score31_geo_bar2$x, ncol = 2))
ggplot(score31_geo_bar2, aes(x=Variety, y=freq, fill=Group.2)) + 
  geom_histogram(position="fill", stat='identity', width=0.8)+ 
  theme(panel.background = element_blank(),
        axis.line = element_line(),
        plot.title = element_text(hjust = 0.5)) + 
  geom_text(aes(label = as.character(freq)), size = 5, position = position_fill(vjust = 0.5)) +
  scale_fill_manual(values = c('dodgerblue', 'firebrick1')) +
  ggtitle('GSE65904')+
  labs(x = NULL, y = "Relative Percent(x100%)") +
  guides(fill=guide_legend(title = NULL)) +
  annotate(geom = 'text', label='P<0.001', x = 1.5, y = 1.05)

##2、基于不同临床病理亚组分析 TIPRGPI 特征预后价值
pheno2_score = merge(pheno2, train_lasso[,c('sample', 'cluster', 'TIPRGPI', 'OS', 'OS.time')], by='sample')
#Age
plot_surv(pheno2_score %>% filter(Age=='>=60'), c('High TIPRGPI','Low TIPRGPI'), 3, '')
plot_surv(pheno2_score %>% filter(Age=='<60'), c('High TIPRGPI','Low TIPRGPI'), 3, '')
#gender
plot_surv(pheno2_score %>% filter(Gender==1), c('High TIPRGPI','Low TIPRGPI'), 3, '')
plot_surv(pheno2_score %>% filter(Gender==0), c('High TIPRGPI','Low TIPRGPI'), 3, '')
#stage
plot_surv(pheno2_score %>% filter(tumor_stage.diagnoses %in% c('stage 0', 'stage i', 'stage ii')), c('High TIPRGPI','Low TIPRGPI'), 3, '')
plot_surv(pheno2_score %>% filter(tumor_stage.diagnoses %in% c('stage iii', 'stage iv')), c('High TIPRGPI','Low TIPRGPI'), 3, '')

##Section4：TIPRGPI可作为癌型的独立预后因子及指导临床应用的列线图构建（基于训练集分析）
##1、分析不同的临床因子在TIPRGPI高低风险组的分布差异，标注显著性P值
pheno2_score2 = merge(pheno2_score,
                         pheno31[,c('sample', 'TCGA.Subtype')], by='sample')
pheno2_score2 = merge(pheno2_score2,
                         pheno[,c('sample', 'sample_type.samples', 'bmi.exposures')], by='sample')
pheno2_score2$pathologic_T = mapvalues(pheno2_score2$pathologic_T, from = c('TX'), to = c(NA))
pheno2_score2$pathologic_N = mapvalues(pheno2_score2$pathologic_N, from = c('NX'), to = c(NA))
pheno2_score2$TCGA.Subtype = mapvalues(pheno2_score2$TCGA.Subtype, from = c('SKCM.-'), to = c(NA))
pheno2_score2$sample_type.samples = mapvalues(pheno2_score2$sample_type.samples, from = c('Additional Metastatic'), to = c('Metastatic'))

dif_cluster_clin = data.frame()
for (i in c('Age', 'gender.demographic', 'tumor_stage.diagnoses', 'pathologic_T', 'pathologic_N', 'pathologic_M', 'TCGA.Subtype', 'sample_type.samples', 'melanoma_clark_level_value')) {
  aa = setdiff(unique(pheno2_score2[,i]), c('', NA))
  aa2 = as.data.frame(matrix(rep(0, length(aa)*2), nrow=2))
  rownames(aa2) = c("High TIPRGPI", "Low TIPRGPI")
  colnames(aa2) = aa
  for (j in aa) {
    aa2['High TIPRGPI', j] = dim(pheno2_score2[(pheno2_score2$cluster=='High TIPRGPI')&(pheno2_score2[,i]==j),])[1]
    aa2['Low TIPRGPI', j] = dim(pheno2_score2[(pheno2_score2$cluster=='Low TIPRGPI')&(pheno2_score2[,i]==j),])[1]
  }
  aa_result = chisq.test(aa2)
  dif_cluster_clin = rbind(dif_cluster_clin, data.frame('clin'=i, 'statistics'=aa_result$statistic, 'pvalue'=aa_result$p.value))
}
for (i in c('breslow_depth_value', 'bmi.exposures')) {
  aa = wilcox.test(as.numeric(pheno2_score2[(pheno2_score2$cluster=='High TIPRGPI'),i]),
                   as.numeric(pheno2_score2[(pheno2_score2$cluster=='Low TIPRGPI'),i]))
  dif_cluster_clin = rbind(dif_cluster_clin, data.frame('clin'=i, 'statistics'=aa$statistic, 'pvalue'=aa$p.value))
}
rownames(dif_cluster_clin) = dif_cluster_clin$clin
write.table(dif_cluster_clin, 'output/4.1_dif_cluster_clin.txt', sep = '\t', quote = F, row.names = F, col.names = T)
#热图
set1 = brewer.pal(n = 9, name = "Set1")
set2 = brewer.pal(n = 8, name = "Set2")
set3 = brewer.pal(n = 12, name = "Set3")
col1 = list(Age=c('>=60'=set1[1], '<60'=set1[2]),
            gender=c('male'=set1[3], 'female'=set1[4]),
            stage=c('stage i'=set1[5], 'stage ii'=set1[6], 'stage iii'=set1[7], 'stage iv'=set1[8], 'stage 0'=set1[9]),
            pathologic_T=c('T0'=set2[1], 'T1'=set2[2], 'T2'=set2[3], 'T3'=set2[4], 'T4'=set2[5], 'Tis'=set2[6]),
            pathologic_N=c('N0'=set2[8], 'N1'=set3[1], 'N2'=set3[2], 'N3'=set3[3]),
            pathologic_M=c('M0'=set2[7], 'M1'=set3[4]),
            TCGA.Subtype=c('SKCM.BRAF_Hotspot_Mutants'=set3[5], 'SKCM.NF1_Any_Mutants'=set3[6], 'SKCM.RAS_Hotspot_Mutants'=set3[7], 'SKCM.Triple_WT'=set3[8]),
            sample_type.samples=c('Primary Tumor'=set3[9], 'Metastatic'=set3[10]),
            melanoma_clark_level_value=c('I'=set3[11], 'II'=set3[12], 'III'=set1[1], 'IV'=set1[2], 'V'=set1[3]),
            breslow_depth_value = colorRamp2(c(min(pheno2_score2$breslow_depth_value, na.rm = T), max(pheno2_score2$breslow_depth_value, na.rm = T)), c("white", set3[5])),
            bmi.exposures = colorRamp2(c(min(pheno2_score2$bmi.exposures, na.rm = T), max(pheno2_score2$bmi.exposures, na.rm = T)), c("white", set3[5])),
            cluster=ann_colors[[3]])
ha_clin<-HeatmapAnnotation(cluster=pheno2_score2$cluster,
                           Age=pheno2_score2$Age,
                           gender=pheno2_score2$gender.demographic,
                           stage=pheno2_score2$tumor_stage.diagnoses,
                           pathologic_T=pheno2_score2$pathologic_T,
                           pathologic_N=pheno2_score2$pathologic_N,
                           pathologic_M=pheno2_score2$pathologic_M,
                           TCGA_Subtype=pheno2_score2$TCGA.Subtype,
                           Sample_type=pheno2_score2$sample_type.samples,
                           melanoma_clark_level_value=pheno2_score2$melanoma_clark_level_value,
                           BMI=pheno2_score2$bmi.exposures,
                           breslow_depth_value=pheno2_score2$breslow_depth_value,
                           show_annotation_name = T,
                           #annotation_label = 'Pvalue',
                           annotation_name_side = 'left',
                           which = 'column',
                           annotation_name_gp = gpar(fontsize = 10),
                           height = unit(20, "mm"),
                           col=col1
)
tcga_tumor_ids_by_score = pheno2_score2 %>% arrange(cluster) %>% pull(sample)
pdf("output/4.1_heatmap.pdf", width=20, height=7, onefile = F)
Heatmap(tcga_rna_tumor_mean[2,tcga_tumor_ids_by_score], show_column_names = F, column_order = tcga_tumor_ids_by_score, top_annotation = ha_clin, cluster_rows = F, column_split = pheno2_score2$cluster, cluster_columns = T)
dev.off()

##2、分析不同的临床因子分析TIPRGPI得分差异
#Age
cps2_clin = t(combn(unique(pheno2_score$Age), 2))
cps2_clin = split(cps2_clin,list(cps2_clin[,1],cps2_clin[,2]), drop = T)
ggviolin(pheno2_score, 'Age', 'TIPRGPI', color = 'black', fill = 'Age',
         palette = c("red", "blue"),
         add = c('boxplot'), add.params = list(fill='white')) +
  theme(axis.title.x = element_blank(),
        legend.position = 'none') +
  labs(y='TIPRGPI') +
  stat_compare_means(comparisons = cps2_clin, label = 'p.signif', label.y = 1.3)
#gender.demographic
cps2_clin = t(combn(unique(pheno2_score$gender.demographic), 2))
cps2_clin = split(cps2_clin,list(cps2_clin[,1],cps2_clin[,2]), drop = T)
ggviolin(pheno2_score, 'gender.demographic', 'TIPRGPI', color = 'black', fill = 'gender.demographic',
         palette = c("red", "blue"),
         add = c('boxplot'), add.params = list(fill='white')) +
  theme(axis.title.x = element_blank(),
        legend.position = 'none') +
  labs(y='TIPRGPI') +
  stat_compare_means(comparisons = cps2_clin, label = 'p.signif', label.y = 1.3)
#tumor_stage.diagnoses
cps2_clin = t(combn(unique(Recode(pheno2_score2$tumor_stage.diagnoses, "c('stage 0', 'stage i')='stage 0,i';c('stage ii', 'stage iii', 'stage iv')='stage ii,iii,iv'")), 2))
cps2_clin = split(cps2_clin,list(cps2_clin[,1],cps2_clin[,2]), drop = T)
ggviolin(pheno2_score2 %>% mutate(stage=Recode(pheno2_score2$tumor_stage.diagnoses, "c('stage 0', 'stage i')='stage 0,i';c('stage ii', 'stage iii', 'stage iv')='stage ii,iii,iv'")) %>% filter(!is.na(stage)), 'stage', 'TIPRGPI', color = 'black', fill = 'stage',
         palette = c("red", "blue", "orange", "green", 'purple'),
         add = c('boxplot'), add.params = list(fill='white')) +
  theme(axis.title.x = element_blank(),
        legend.position = 'none') +
  labs(y='TIPRGPI') +
  stat_compare_means(comparisons = cps2_clin, label = 'p.signif', label.y = 1.2)
#pathologic_T
cps2_clin = t(combn(unique(Recode(pheno2_score2$pathologic_T, "c('Tis', 'T0', 'T1')='Tis,T0,T1';c('T2', 'T3', 'T4')='T2,T3,T4'")), 2))
cps2_clin = split(cps2_clin,list(cps2_clin[,1],cps2_clin[,2]), drop = T)
ggviolin(pheno2_score2 %>% mutate(pathologic_T=Recode(pheno2_score2$pathologic_T, "c('Tis', 'T0', 'T1')='Tis,T0,T1';c('T2', 'T3', 'T4')='T2,T3,T4'")) %>% filter(!is.na(pathologic_T)), 'pathologic_T', 'TIPRGPI', color = 'black', fill = 'pathologic_T',
         palette = c("red", "blue", "orange", "green", 'purple', 'yellow', 'pink'),
         add = c('boxplot'), add.params = list(fill='white')) +
  theme(axis.title.x = element_blank(),
        legend.position = 'none') +
  labs(y='TIPRGPI') +
  stat_compare_means(comparisons = cps2_clin, label = 'p.signif', label.y = 1.2)
#pathologic_N
cps2_clin = t(combn(unique(Recode(pheno2_score2$pathologic_N, "c('N0', 'N1')='N0,N1';c('N2', 'N3')='N2,N3'")), 2))
cps2_clin = split(cps2_clin,list(cps2_clin[,1],cps2_clin[,2]), drop = T)
ggviolin(pheno2_score2 %>% mutate(pathologic_N=Recode(pheno2_score2$pathologic_N, "c('N0', 'N1')='N0,N1';c('N2', 'N3')='N2,N3'")) %>% filter(!is.na(pathologic_N)), 'pathologic_N', 'TIPRGPI', color = 'black', fill = 'pathologic_N',
         palette = c("red", "blue", "orange", "green", 'purple', 'yellow', 'pink'),
         add = c('boxplot'), add.params = list(fill='white')) +
  theme(axis.title.x = element_blank(),
        legend.position = 'none') +
  labs(y='TIPRGPI') +
  stat_compare_means(comparisons = cps2_clin, label = 'p.signif', label.y = 1.2)
#pathologic_M
cps2_clin = t(combn(unique(pheno2_score$pathologic_M), 2))
cps2_clin = split(cps2_clin,list(cps2_clin[,1],cps2_clin[,2]), drop = T)
ggviolin(pheno2_score %>% filter(!is.na(pathologic_M)), 'pathologic_M', 'TIPRGPI', color = 'black', fill = 'pathologic_M',
         palette = c("red", "blue", "orange", "green", 'purple', 'yellow', 'pink'),
         add = c('boxplot'), add.params = list(fill='white')) +
  theme(axis.title.x = element_blank(),
        legend.position = 'none') +
  labs(y='TIPRGPI') +
  stat_compare_means(comparisons = cps2_clin, label = 'p.signif', label.y=1.3)
#breslow_depth_value
pheno2_score2$breslow_depth_value_cluster = ifelse(pheno2_score2$breslow_depth_value<=4.5, '<=4.5', '>4.5')
cps2_clin = t(combn(unique(pheno2_score2$breslow_depth_value_cluster), 2))
cps2_clin = split(cps2_clin,list(cps2_clin[,1],cps2_clin[,2]), drop = T)
ggviolin(pheno2_score2 %>% filter(!is.na(breslow_depth_value_cluster)), 'breslow_depth_value_cluster', 'TIPRGPI', color = 'black', fill = 'breslow_depth_value_cluster',
         palette = c("red", "blue", "orange", "green", 'purple', 'yellow', 'pink'),
         add = c('boxplot'), add.params = list(fill='white')) +
  theme(axis.title.x = element_blank(),
        legend.position = 'none') +
  labs(y='TIPRGPI') +
  stat_compare_means(comparisons = cps2_clin, label = 'p.signif', label.y = 1.2)
#TCGA.Subtype
cps2_clin = t(combn(unique(pheno2_score2$TCGA.Subtype), 2))
cps2_clin = split(cps2_clin,list(cps2_clin[,1],cps2_clin[,2]), drop = T)
ggviolin(pheno2_score2 %>% filter(!is.na(TCGA.Subtype)), 'TCGA.Subtype', 'TIPRGPI', color = 'black', fill = 'TCGA.Subtype',
         palette = c("red", "blue", "orange", "green", 'purple', 'yellow', 'pink'),
         add = c('boxplot'), add.params = list(fill='white')) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 10, hjust=1),
        legend.position = 'none') +
  labs(y='TIPRGPI') +
  stat_compare_means(comparisons = cps2_clin, label = 'p.signif')
#sample_type.samples
cps2_clin = t(combn(unique(pheno2_score2$sample_type.samples), 2))
cps2_clin = split(cps2_clin,list(cps2_clin[,1],cps2_clin[,2]), drop = T)
ggviolin(pheno2_score2 %>% filter(!is.na(sample_type.samples)), 'sample_type.samples', 'TIPRGPI', color = 'black', fill = 'sample_type.samples',
         palette = c("red", "blue", "orange", "green", 'purple', 'yellow', 'pink'),
         add = c('boxplot'), add.params = list(fill='white')) +
  theme(axis.title.x = element_blank(),
        legend.position = 'none') +
  labs(y='TIPRGPI') +
  stat_compare_means(comparisons = cps2_clin, label = 'p.signif', label.y = 1.2)
#melanoma_clark_level_value
cps2_clin = t(combn(unique(Recode(pheno2_score2$melanoma_clark_level_value, "c('I', 'II', 'III')='I,II,III';c('IV', 'V')='IV,V'")), 2))
cps2_clin = split(cps2_clin,list(cps2_clin[,1],cps2_clin[,2]), drop = T)
ggviolin(pheno2_score2 %>% mutate(melanoma_clark_level_value=Recode(pheno2_score2$melanoma_clark_level_value, "c('I', 'II', 'III')='I,II,III';c('IV', 'V')='IV,V'")) %>% filter(!is.na(melanoma_clark_level_value)), 'melanoma_clark_level_value', 'TIPRGPI', color = 'black', fill = 'melanoma_clark_level_value',
         palette = c("red", "blue", "orange", "green", 'purple', 'yellow', 'pink'),
         add = c('boxplot'), add.params = list(fill='white')) +
  theme(axis.title.x = element_blank(),
        legend.position = 'none') +
  labs(y='TIPRGPI') +
  stat_compare_means(comparisons = cps2_clin, label = 'p.signif', label.y = 1.2)
#bmi.exposures
pheno2_score2$bmi_cluster = ifelse(pheno2_score2$bmi.exposures<24, '<24', '>=24')
cps2_clin = t(combn(unique(pheno2_score2$bmi_cluster), 2))
cps2_clin = split(cps2_clin,list(cps2_clin[,1],cps2_clin[,2]), drop = T)
ggviolin(pheno2_score2 %>% filter(!is.na(bmi_cluster)), 'bmi_cluster', 'TIPRGPI', color = 'black', fill = 'bmi_cluster',
         palette = c("red", "blue", "orange", "green", 'purple', 'yellow', 'pink'),
         add = c('boxplot'), add.params = list(fill='white')) +
  theme(axis.title.x = element_blank(),
        legend.position = 'none') +
  labs(y='TIPRGPI') +
  stat_compare_means(comparisons = cps2_clin, label = 'p.signif', label.y = 1.2)
##3、单因素、多因素cox回归分析TIPRGPI的独立预后效能
#*tcga
#单因素
score_pheno = pheno2_score[,c('sample', 'Age', 'gender.demographic', 'pathologic_N', 'pathologic_M', 'cluster', 'OS', 'OS.time')]
score_pheno[score_pheno=='']=NA
ind_cox_result = data.frame()
for (i in c(2:6)) {
  print(i)
  score_pheno$'x' = score_pheno[,i]
  model = coxph(Surv(OS.time, OS) ~ x, data=score_pheno[!is.na(score_pheno[,i]),])
  print(c(colnames(score_pheno)[i], cox.zph(model)$table))
  coef = cbind(summary(model)$coefficients, summary(model)$conf.int)
  #rownames(coef)[1] = colnames(pheno15)[i]
  ind_cox_result = rbind(ind_cox_result, coef) #254基因*9列cox属性
}
write.table(ind_cox_result, file='output/4.3_uni_cox_result.txt', row.names=T, col.names=T, quote=F, sep='\t')
#画森林图
ind_cox_result24 = ind_cox_result
ind_cox_result24$`upper .95` = ifelse(ind_cox_result24$`upper .95`==Inf, max(ind_cox_result24$`upper .95`[ind_cox_result24$`upper .95`!=Inf]), ind_cox_result24$`upper .95`)
ind_table_text24 = rbind(c('Risk factors', '', 'HR (95% CI)', 'Pvalue'),
                         cbind(c('Age\n(vs <60)', 'gender.demographic\n(vs Female)', 'pathologic_N\n(vs N0)', '', '', '', 'pathologic_M\n(vs M0)', 'TIPRGPI cluster\n(vs High TIPRGPI)'),
                               c('>=60', 'Male', 'N1', 'N2', 'N3', 'NX', 'M1', 'Low TIPRGPI'),
                               as.character(t(paste0(round(ind_cox_result24$`exp(coef)`,2), ' (', round(ind_cox_result24$`lower .95`,2), '-', round(ind_cox_result24$`upper .95`,2),')\n'))),
                               unlist(lapply(ind_cox_result24$`Pr(>|z|)`, function(x) {c(round(x, 3))})))
)
ind_table_text24[ind_table_text24=='NA (NA-NA)\n'] = '1\n'
ind_table_text24[ind_table_text24==0] = '<0.001'
pdf("output/4.3_forestplot_uni.pdf", width=10, height=7, onefile = F)
forestplot(labeltext=ind_table_text24, graph.pos=3,
           mean=c(NA, round(ind_cox_result24$`exp(coef)`,2)),
           lower=c(NA, round(ind_cox_result24$`lower .95`,2)), 
           upper=c(NA, round(ind_cox_result24$`upper .95`,2)),
           title="Univariate Cox Regression",
           #           txt_gp=fpTxtGp(label=gpar(cex=1.2),
           #                          ticks=gpar(cex=0.8),
           #                          xlab=gpar(cex = 1.2),
           #                          title=gpar(cex = 1.2)),
           col=fpColors(box='red', lines='black', zero = "gray50"),
           zero=1,
           cex=0.9, lineheight = "auto",
           #colgap=unit(8,"mm"),
           lwd.ci=2, 
           boxsize=0.2,
           xlab="<---Lesser hazard---    ---Greater hazard--->",
           hrzl_lines=list("2" = gpar(lwd=1, col="black"))
           #ci.vertices=TRUE, 
           #ci.vertices.height = 0.4
)
dev.off()
#*多因素
model_cox = coxph(formula(paste0('Surv(OS.time, OS) ~ ', paste0('`', colnames(score_pheno)[c(2,6)], '`', collapse = '+'))), data=score_pheno)
print(cox.zph(model_cox))
ind_cox_result2 = cbind(summary(model_cox)$coefficients, summary(model_cox)$conf.int)
write.table(ind_cox_result2, file='output/4.3_multi_cox_result.txt', row.names=T, col.names=T, quote=F, sep='\t')
#*画森林图
ind_cox_result242 = as.data.frame(ind_cox_result2)
ind_cox_result242$`upper .95`[ind_cox_result242$`upper .95`==Inf] = max(ind_cox_result242$`upper .95`[ind_cox_result242$`upper .95`!=Inf])
ind_table_text242 = rbind(c('Risk factors', '', 'HR (95% CI)', 'Pvalue'),
                          cbind(c('Age\n(vs <60)', 'TIPRGPI cluster\n(vs High TIPRGPI)'),
                                c('>=60', 'Low TIPRGPI'),
                                as.character(t(paste0(signif(ind_cox_result242$`exp(coef)`,2), ' (', signif(ind_cox_result242$`lower .95`,2), '-', signif(ind_cox_result242$`upper .95`,2),')\n'))),
                                unlist(lapply(ind_cox_result242$`Pr(>|z|)`, function(x) {c(round(x, 3))})))
)
ind_table_text242[ind_table_text242=='NA (NA-NA)\n'] = '1\n'
ind_table_text242[ind_table_text242==0] = '<0.001'
pdf("output/4.3_forestplot_multi.pdf", width=10, height=5, onefile = F)
forestplot(labeltext=ind_table_text242, graph.pos=3,
           mean=c(NA, signif(ind_cox_result242$`exp(coef)`,2)),
           lower=c(NA, signif(ind_cox_result242$`lower .95`,2)), 
           upper=c(NA, signif(ind_cox_result242$`upper .95`,2)),
           title="Multivariate Cox Regression",
           #           txt_gp=fpTxtGp(label=gpar(cex=1.2),
           #                          ticks=gpar(cex=0.8),
           #                          xlab=gpar(cex = 1.2),
           #                          title=gpar(cex = 1.2)),
           col=fpColors(box='red', lines='black', zero = "gray50"),
           zero=1,
           cex=0.9, lineheight = "auto",
           #colgap=unit(8,"mm"),
           lwd.ci=2, 
           boxsize=0.2,
           xlab="<---Lesser hazard---    ---Greater hazard--->",
           hrzl_lines=list("2" = gpar(lwd=1, col="black"))
           #ci.vertices=TRUE, 
           #ci.vertices.height = 0.4
)
dev.off()

##4、结合临床因子整合TIPRGPI构建列线图nomogram指导临床分析预测显著性
#nomogram列线图
d <- datadist(score_pheno); options(datadist='d')
f = cph(formula(paste0('Surv(OS.time, OS) ~ ', paste0('`', colnames(score_pheno)[c(2,6)], '`', collapse = '+'))), data=score_pheno, x=T, y=T, surv = T)
regplot(f, failtime = c(12,24,36), prfail = T, droplines=T, title = 'coxph regression')
#校准曲线
p1<- calibrate(f,#模型名称
               cmethod='KM',
               method='boot',#检测方法
               u=12*1,#评估的时间，注：一定要与模型的时间一致
               m=nrow(score_pheno)/5, #每次抽样的样本量，
               B=1000)#抽样次数
#注，m值的确定：m=数据总数/3-4,即你想让最终的校准曲线有3个点，那就是m=数据总数/3
#B值一般1000，电脑配置不好可以选500,300,100等
plot(p1,
     add=F,#增加第二条线
     conf.int=T,#95%CI
     subtitles = F,#副标题
     cex.subtitles=0.8, #副标题大小
     lwd=2,#95%CI粗细
     lty=1,#95%CI实线，2=虚线
     errbar.col="black",#95%CI颜色
     xlim=c(0.0,1),#x轴范围
     ylim=c(0.0,1),
     xlab="Nomogram-predicted OS(%)",
     ylab="Observed OS(%)",
     col="red")#曲线颜色
p2<- calibrate(f,#模型名称
               cmethod='KM',
               method='boot',#检测方法
               u=12*2,#评估的时间，注：一定要与模型的时间一致
               m=nrow(score_pheno)/5, #每次抽样的样本量，
               B=1000)#抽样次数
#注，m值的确定：m=数据总数/3-4,即你想让最终的校准曲线有3个点，那就是m=数据总数/3
#B值一般1000，电脑配置不好可以选500,300,100等
plot(p2,
     add=T,#增加第二条线
     conf.int=T,#95%CI
     subtitles = F,#副标题
     cex.subtitles=0.8, #副标题大小
     lwd=2,#95%CI粗细
     lty=1,#95%CI实线，2=虚线
     errbar.col="black",#95%CI颜色
     xlim=c(0.0,1),#x轴范围
     ylim=c(0.0,1),
     xlab="Nomogram-predicted 3 year OS(%)",
     ylab="Observed 3 year OS(%)",
     col="orange")#曲线颜色
p3<- calibrate(f,#模型名称
               cmethod='KM',
               method='boot',#检测方法
               u=12*3,#评估的时间，注：一定要与模型的时间一致
               m=nrow(score_pheno)/5, #每次抽样的样本量，
               B=1000)#抽样次数
#注，m值的确定：m=数据总数/3-4,即你想让最终的校准曲线有3个点，那就是m=数据总数/3
#B值一般1000，电脑配置不好可以选500,300,100等
plot(p3,
     add=T,#增加第二条线
     conf.int=T,#95%CI
     subtitles = F,#副标题
     cex.subtitles=0.8, #副标题大小
     lwd=2,#95%CI粗细
     lty=1,#95%CI实线，2=虚线
     errbar.col="black",#95%CI颜色
     xlim=c(0.0,1),#x轴范围
     ylim=c(0.0,1),
     xlab="Nomogram-predicted 3 year OS(%)",
     ylab="Observed 3 year OS(%)",
     col="blue")#曲线颜色
legend("bottomright",c('1-year', '2-year', '3-year'),col=c("red","orange","blue"),lty=1,lwd=2)
#决策曲线
dca_age <- decision_curve(OS~Age,data = pheno2_score, 
                      family = binomial(link ='logit'),#模型类型，这里是二分类
                      thresholds= seq(0,1, by = 0.01),
                      confidence.intervals =0.95,#95可信区间
                      study.design = 'cohort')
dca_cluster <- decision_curve(OS~`TIPRGPI subtype`,data = pheno2_score %>% rename('TIPRGPI subtype'=cluster), 
                          family = binomial(link ='logit'),#模型类型，这里是二分类
                          thresholds= seq(0,1, by = 0.01),
                          confidence.intervals =0.95,#95可信区间
                          study.design = 'cohort')
dca_all <- decision_curve(OS~`TIPRGPI subtype`+Age,data = pheno2_score %>% rename('TIPRGPI subtype'=cluster), 
                      family = binomial(link ='logit'),#模型类型，这里是二分类
                      thresholds= seq(0,1, by = 0.01),
                      confidence.intervals =0.95,#95可信区间
                      study.design = 'cohort')
dca_list = list(dca_age, dca_cluster, dca_all)

plot_decision_curve(dca_list,curve.names= c('Age','TIPRGPI', 'Age+TIPRGPI'),
                    cost.benefit.axis =FALSE,col = c('red','blue','green'),
                    confidence.intervals =FALSE,standardize = FALSE)
#KM生存曲线
pheno2_score2 = cbind(pheno2_score, data.frame('monogram'=f$linear.predictors))
pheno2_score2$monogram_cluster = ifelse(pheno2_score2$monogram>median(pheno2_score2$monogram), 'High Risk', 'Low Risk')
ann_colors = c(ann_colors, list(cluster=c('High Risk'='red', 'Low Risk'='blue')))
plot_surv(pheno2_score2 %>% select(-cluster) %>% rename('cluster'='monogram_cluster'), c('High Risk', 'Low Risk'), 4, '')

##Section5：TIPRGPI的潜在分子机制
##1.分析高低TIPRGPI风险组突变基因展示（附表给出差异突变基因显著性P值）、共突变分析（可取top）
tcga_tumor_ids_high = pheno2_score$sample[pheno2_score$cluster=='High TIPRGPI']
tcga_tumor_ids_low = pheno2_score$sample[pheno2_score$cluster=='Low TIPRGPI']
snv_tumor_high = subsetMaf(snv, tsb = as.character(snv@clinical.data$Tumor_Sample_Barcode)[unlist(lapply(as.character(snv@clinical.data$Tumor_Sample_Barcode), FUN = function(x) {paste(strsplit(x, '-')[[1]][1:4], collapse = '-') %in% tcga_tumor_ids_high}))])
snv_tumor_low = subsetMaf(snv, tsb = as.character(snv@clinical.data$Tumor_Sample_Barcode)[unlist(lapply(as.character(snv@clinical.data$Tumor_Sample_Barcode), FUN = function(x) {paste(strsplit(x, '-')[[1]][1:4], collapse = '-') %in% tcga_tumor_ids_low}))])
pdf('output/5.1_snv_high.pdf', width = 6, height = 6)
oncoplot(maf = snv_tumor_high, draw_titv = T, writeMatrix = T, top = 20)  #writeMatrix会在工作目录下生成onco_matrix.txt文件
dev.off()
snv_tumor_high_matrix = read.delim('output/5.2_onco_matrix_high.txt', sep = '\t', header = T, check.names = F)
snv_genes_high = rownames(snv_tumor_high_matrix)
pdf('output/5.1_snv_low.pdf', width = 6, height = 6)
oncoplot(maf = snv_tumor_low, draw_titv = T, writeMatrix = T, top = 20)
dev.off()
snv_tumor_low_matrix = read.delim('output/5.2_onco_matrix_low.txt', sep = '\t', header = T, check.names = F)
snv_genes_low = rownames(snv_tumor_low_matrix)
snv_p = data.frame()
for (gene in intersect(snv_genes_high, snv_genes_low)) {
  aa = data.frame('Mutation'=c(sum(snv_tumor_high_matrix[gene,]!=''), sum(snv_tumor_low_matrix[gene,]!='')),
                  'Non-Mutation'=c(sum(snv_tumor_high_matrix[gene,]==''), sum(snv_tumor_low_matrix[gene,]=='')),
                  check.names = F)
  aa_result = fisher.test(aa)
  snv_p = rbind(snv_p,
                data.frame('log2_OR'=log2(aa_result$estimate),
                           'lower.95'=log2(aa_result$conf.int[1]), 'upper.95'=log2(aa_result$conf.int[2]),
                           'pvalue'=aa_result$p.value,
                           'high'=sum(snv_tumor_high_matrix[gene,]!=''),
                           'low'=sum(snv_tumor_low_matrix[gene,]!='')))
}
rownames(snv_p) = intersect(snv_genes_high, snv_genes_low)
snv_p_signif = snv_p %>% arrange(pvalue) %>% filter(pvalue<0.05) %>% .[1:10,]
snv_p_signif = snv_p_signif %>% mutate(`p-value` = case_when(pvalue<=0.0001~'****', pvalue<=0.001~'***', pvalue<=0.01~'**', pvalue<=0.05~'*'))  ##log2(OR)值均<0，结合low risk中这些基因的突变数比high risk中多，说明这些基因为保护因素
write.table(snv_p_signif, 'output/5.1_snv_p.txt', sep = '\t', quote = F, row.names = T, col.names = T)
#森林图
snv_table_text24 = rbind(c('', 'low TIPRGPI', 'high TIPRGPI', 'Pvalue'),
                         cbind(rownames(snv_p_signif),
                               snv_p_signif$low,
                               snv_p_signif$high,
                               snv_p_signif$`p-value`)
)
pdf("output/5.1_forestplot_snv.pdf", width=7, height=4, onefile = F)
forestplot(labeltext=snv_table_text24, graph.pos=1,
           mean=c(NA, round(snv_p_signif$log2_OR,2)),
           lower=c(NA, round(snv_p_signif$lower.95,2)), 
           upper=c(NA, round(snv_p_signif$upper.95,2)),
           title='High TIPRGPI v/s high TIPRGPI',
           #           txt_gp=fpTxtGp(label=gpar(cex=1.2),
           #                          ticks=gpar(cex=0.8),
           #                          xlab=gpar(cex = 1.2),
           #                          title=gpar(cex = 1.2)),
           col=fpColors(box='red', lines='black', zero = "gray50"),
           zero=0,
           cex=0.9, lineheight = "auto",
           #colgap=unit(8,"mm"),
           lwd.ci=2, 
           boxsize=0.2,
           xlab="Log odds ratio",
           hrzl_lines=list("2" = gpar(lwd=1, col="black"))
           #ci.vertices=TRUE, 
           #ci.vertices.height = 0.4
)
dev.off()
#棒棒糖图
lollipopPlot2(m1 = snv_tumor_high, m2 = snv_tumor_low, gene = 'TENM1', m1_name = 'High TIPRGPI', m2_name = 'Low TIPRGPI')
#相关性图
snv_high = snv_tumor_high_matrix
snv_high[snv_high!=''] = 1
snv_high[snv_high==''] = 0
snv_high = snv_high[1:20,]
snv_cor_high = rcorr(t(snv_high))
snv_cor_high2 = snv_cor_high
snv_cor_high2$r = snv_cor_high2$r[order(rownames(snv_cor_high2$r)), order(colnames(snv_cor_high2$r))]
snv_cor_high2$P = snv_cor_high2$P[order(rownames(snv_cor_high2$P)), order(colnames(snv_cor_high2$P))]
snv_cor_high2$P[snv_cor_high2$P==0] = 1e-99
snv_cor_high2$P2 = -log10(snv_cor_high2$P)
snv_cor_high2$P2[is.na(snv_cor_high2$P2)] = 0
snv_cor_high2$P2 = snv_cor_high2$P2 * (snv_cor_high2$r/abs(snv_cor_high2$r))
snv_cor_high2$P = ifelse(snv_cor_high2$P<=0.0001, '****', ifelse(snv_cor_high2$P<=0.001, '***', ifelse(snv_cor_high2$P<0.01, '**', ifelse(snv_cor_high2$P<0.05, '*', ''))))
snv_cor_high2$P[is.na(snv_cor_high2$P)] = ''
snv_genes_high_top20 = snv_genes_high[1:20]
pdf("output/5.1_cibersort_high.pdf", width=13/2, height=10/2, onefile = F)
pheatmap(snv_cor_high2$P2[snv_genes_high_top20, snv_genes_high_top20], color = colorRampPalette(colors = c("blue","white","red"))(100), display_numbers = snv_cor_high2$P[snv_genes_high_top20, snv_genes_high_top20], cluster_rows = F, cluster_cols = F, angle_col = '90', breaks = unique(c(seq(-3,3, length=100))), fontsize = 10)
dev.off()

snv_low = snv_tumor_low_matrix
snv_low[snv_low!=''] = 1
snv_low[snv_low==''] = 0
snv_low = snv_low[1:20,]
snv_cor_low = rcorr(t(snv_low))
snv_cor_low2 = snv_cor_low
snv_cor_low2$r = snv_cor_low2$r[order(rownames(snv_cor_low2$r)), order(colnames(snv_cor_low2$r))]
snv_cor_low2$P = snv_cor_low2$P[order(rownames(snv_cor_low2$P)), order(colnames(snv_cor_low2$P))]
snv_cor_low2$P[snv_cor_low2$P==0] = 1e-99
snv_cor_low2$P2 = -log10(snv_cor_low2$P)
snv_cor_low2$P2[is.na(snv_cor_low2$P2)] = 0
snv_cor_low2$P2 = snv_cor_low2$P2 * (snv_cor_low2$r/abs(snv_cor_low2$r))
snv_cor_low2$P = ifelse(snv_cor_low2$P<=0.0001, '****', ifelse(snv_cor_low2$P<=0.001, '***', ifelse(snv_cor_low2$P<0.01, '**', ifelse(snv_cor_low2$P<0.05, '*', ''))))
snv_cor_low2$P[is.na(snv_cor_low2$P)] = ''
snv_genes_low_top20 = snv_genes_low[1:20]
pdf("output/5.1_cibersort_low.pdf", width=13/2, height=10/2, onefile = F)
pheatmap(snv_cor_low2$P2[snv_genes_low_top20, snv_genes_low_top20], color = colorRampPalette(colors = c("blue","white","red"))(100), display_numbers = snv_cor_low2$P[snv_genes_low_top20, snv_genes_low_top20], cluster_rows = F, cluster_cols = F, angle_col = '90', breaks = unique(c(seq(-3,3, length=100))), fontsize = 10)
dev.off()

snv_cor_P2 = snv_cor_low2$P2[snv_genes_low_top20, snv_genes_low_top20]
for (i in 1:20) {
  for (j in 1:i) {
    print(c(i,j))
    snv_cor_P2[i,j] = snv_cor_high2$P2[snv_genes_high_top20, snv_genes_high_top20][i,j]
  }
}
snv_cor_P = snv_cor_low2$P[snv_genes_low_top20, snv_genes_low_top20]
for (i in 1:20) {
  for (j in 1:i) {
    print(c(i,j))
    snv_cor_P[i,j] = snv_cor_high2$P[snv_genes_high_top20, snv_genes_high_top20][i,j]
  }
}
annotation_col_snv = data.frame(row.names = snv_genes_low_top20, 'colname'=snv_genes_low_top20)


pdf("output/5.1_cibersort_all.pdf", width=13/2, height=10/2, onefile = F)
pheatmap(snv_cor_P2, color = colorRampPalette(colors = c("blue","white","red"))(100), display_numbers = snv_cor_P, cluster_rows = F, cluster_cols = F, angle_col = '90', breaks = unique(c(seq(-3,3, length=100))), fontsize = 10)
dev.off()




##2.高低TIPRGPI风险组CNV差异、模型包含基因在高低分组的CNV展示以及表达在CNV分组见的展示（取最显著的，若没有这一小步取研究癌型驱动基因cosmic数据库）
#*GISTIC-BRCA2.R下载文件做到marker_file前一步，到网上做gistic，下载返回的scores.gistic文件进行cnv.R下边的代码步骤
#提取01A结尾的样本
seqlengths=c(248956422,242193529,198295559,190214555,181538259,170805979,159345973,145138636,138394717,133797422,135086622,133275309,114364328,107043718,101991189,90338345,83257441,80373285,58617616,64444167,46709983,50818468,156040895,57227415,3088269832)  #hg38
idx_x<-vector()
for (i in 1:25) {
  idx_x[i]=sum(seqlengths[1:i])
}
cnv_h = read.delim('output/5.2_scores_high.gistic',header = T,sep = '\t', check.names=F)
cnv_l = read.delim('output/5.2_scores_low.gistic',header = T,sep = '\t', check.names=F)
cnv_h[cnv_h$Type=='Del',]$frequency=-cnv_h[cnv_h$Type=='Del',]$frequency
cnv_h = cnv_h[order(cnv_h[,'Chromosome'], cnv_h[,'Start']), ]
cnv_h$pos<-cnv_h$Start
for (i in 2:22) {
  cnv_h$pos[cnv_h$Chromosome==i]<-cnv_h$pos[cnv_h$Chromosome==i]+idx_x[i-1]
}
dash_loci_h = data.frame()
for (i in 1:22) {
  dash_loci_h = rbind(dash_loci_h,
                      data.frame('Var1'=i, 'Freq'=length(unique(cnv_h$pos[cnv_h$Chromosome==i]))))
}
cnv_h$dashLine<-0
for (i in 1:22) {
  cnv_h$dashLine[cnv_h$Chromosome==i]<-sum(dash_loci_h$Freq[1:i])
}
cnv_h$'gistic_score' = cnv_h$`G-score`
cnv_h[cnv_h$Type=='Del',]$gistic_score=-cnv_h[cnv_h$Type=='Del',]$gistic_score
cnv_h$'label_x' = cnv_h[1,'dashLine']/2
for (i in 2:22) {
  cnv_h$label_x[cnv_h$Chromosome==i] = (cnv_h$dashLine[cnv_h$Chromosome==(i-1)] + cnv_h$dashLine[cnv_h$Chromosome==i])/2
}
cnv_h$label_y = ifelse(cnv_h$Chromosome%%2==1, max(cnv_h$frequency)*1.2, min(cnv_h$frequency)*1.2)
cnv_h$'label' = paste0('chr', cnv_h$Chromosome)
write.table(cnv_h, 'output/5.2_cnv_h.txt', sep = '\t', quote = F, row.names = F, col.names = T)
pdf("output/5.2_high.pdf", width=14, height=3, onefile = F)
tiff("output/5.2_high.tiff", res=300, width=14, height=3, compression="lzw", units="in")
ggbarplot(cnv_h,x='pos',y='gistic_score',col='Type',fill='Type',palette = c('darkred','darkblue'))+
  theme_bw()+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  geom_vline(aes(xintercept=dashLine), colour="grey", linetype="dashed") +
  annotate('text', x=cnv_h$label_x, y=cnv_h$label_y, label=cnv_h$label) +
  labs(x=NULL, y='gistic score')
dev.off()

cnv_l = read.delim('output/5.2_scores_low.gistic',header = T,sep = '\t', check.names=F)
cnv_l[cnv_l$Type=='Del',]$frequency=-cnv_l[cnv_l$Type=='Del',]$frequency
cnv_l = cnv_l[order(cnv_l[,'Chromosome'], cnv_l[,'Start']), ]
cnv_l$pos<-cnv_l$Start
for (i in 2:22) {
  cnv_l$pos[cnv_l$Chromosome==i]<-cnv_l$pos[cnv_l$Chromosome==i]+idx_x[i-1]
}
dash_loci_l = data.frame()
for (i in 1:22) {
  dash_loci_l = rbind(dash_loci_l,
                      data.frame('Var1'=i, 'Freq'=length(unique(cnv_l$pos[cnv_l$Chromosome==i]))))
}
cnv_l$dashLine<-0
for (i in 1:22) {
  cnv_l$dashLine[cnv_l$Chromosome==i]<-sum(dash_loci_l$Freq[1:i])
}
cnv_l$'gistic_score' = cnv_l$`G-score`
cnv_l[cnv_l$Type=='Del',]$gistic_score=-cnv_l[cnv_l$Type=='Del',]$gistic_score
cnv_l$'label_x' = cnv_l[1,'dashLine']/2
for (i in 2:22) {
  cnv_l$label_x[cnv_l$Chromosome==i] = (cnv_l$dashLine[cnv_l$Chromosome==(i-1)] + cnv_l$dashLine[cnv_l$Chromosome==i])/2
}
cnv_l$label_y = ifelse(cnv_l$Chromosome%%2==1, max(cnv_l$frequency)*1.2, min(cnv_l$frequency)*1.2)
cnv_l$'label' = paste0('chr', cnv_l$Chromosome)
write.table(cnv_l, 'output/5.2_cnv_l.txt', sep = '\t', quote = F, row.names = F, col.names = T)
pdf("output/5.2_low.pdf", width=14, height=3, onefile = F)
tiff("output/5.2_low.tiff", res=300, width=14, height=3, compression="lzw", units="in")
ggbarplot(cnv_l,x='pos',y='gistic_score',col='Type',fill='Type',palette = c('darkred','darkblue'))+
  theme_bw()+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  geom_vline(aes(xintercept=dashLine), colour="grey", linetype="dashed") +
  annotate('text', x=cnv_l$label_x, y=cnv_l$label_y, label=cnv_l$label) +
  labs(x=NULL, y='gistic score')
dev.off()
#堆叠条形图
cnv_high = cnv_tumor_mean[,which(colnames(cnv_tumor_mean) %in% tcga_tumor_ids_high)]
cnv_low = cnv_tumor_mean[,which(colnames(cnv_tumor_mean) %in% tcga_tumor_ids_low)]
cnv_p = data.frame()
for (gene in hubs_lasso) {
  print(gene)
  aa = data.frame('amp'=c(sum(cnv_high[gene,]==1), sum(cnv_low[gene,]==1)),
                  'non-amp'=c(sum(cnv_high[gene,]!=1), sum(cnv_low[gene,]!=1)),
                  check.names = F)
  aa_result = fisher.test(aa)
  cnv_p = rbind(cnv_p,
                data.frame('log2_OR'=log2(aa_result$estimate),
                           'lower.95'=log2(aa_result$conf.int[1]), 'upper.95'=log2(aa_result$conf.int[2]),
                           'pvalue'=aa_result$p.value,
                           'high'=sum(cnv_high[gene,]==1),
                           'low'=sum(cnv_low[gene,]==1)))
}
rownames(cnv_p) = hubs_lasso
cnv_p_signif = cnv_p %>% arrange(pvalue) %>% filter(pvalue<0.05)
cnv_high2 = cnv_high[rownames(cnv_p_signif),]
cnv_low2 = cnv_low[rownames(cnv_p_signif),]
cnv_high_gain = as.data.frame(apply(cnv_high2, 1, function(x) {sum(x==1)/ncol(cnv_high2)}))
cnv_high_gain$'gene' = rownames(cnv_high_gain)
cnv_high_gain$'Type' = 'High TIPRGPI'
colnames(cnv_high_gain)[1] = 'percentage'
cnv_low_gain = as.data.frame(apply(cnv_low2, 1, function(x) {sum(x==1)/ncol(cnv_low2)}))
cnv_low_gain$'gene' = rownames(cnv_low_gain)
cnv_low_gain$'Type' = 'low TIPRGPI'
colnames(cnv_low_gain)[1] = 'percentage'
cnv_plot = rbind(cnv_high_gain, cnv_low_gain)
cnv_plot = cnv_plot[order(cnv_plot$gene),]
cnv_max = cnv_plot %>% arrange(desc(percentage)) %>% group_by(gene) %>% mutate(sort = 1:n()) %>% filter(sort==1) %>% dplyr::select(-sort) %>% data.frame(stringsAsFactors = F)
cnv_min = cnv_plot %>% arrange(gene, percentage) %>% group_by(gene) %>% mutate(sort = 1:n()) %>% filter(sort==1) %>% dplyr::select(-sort) %>% data.frame(stringsAsFactors = F)
cnv_plot2 = data.frame()
for (gene in cnv_max[order(cnv_max$percentage),]$gene) {
  cnv_plot2 = rbind(cnv_plot2, cnv_plot[which(cnv_plot$gene==gene),])
}
ggplot(cnv_plot2, aes(x=reorder(gene, percentage), y=percentage*100)) +
  geom_bar(data = cnv_max, stat='identity', width=0.1, fill='grey') +
  geom_bar(data = cnv_min, stat='identity', width=0.1, fill='white') +
  geom_point(size=7, aes(color=Type)) +
  scale_color_manual(values = c("red", "blue")) +
  theme(axis.line = element_line(linetype = "solid"), 
        axis.text = element_text(hjust = 1), 
        axis.text.x = element_text(size = 10, vjust = 0.5, angle = 90),
        panel.background = element_rect(fill = NA),
        legend.text = element_text(size=12),
        legend.position = 'top') +
  labs(x = NULL, y = "AMP.frequency(*100%)", fill = "Type") +
  guides(color=guide_legend(title = NULL))
#小提琴图
cnv_high_expr = data.frame('sample'=colnames(cnv_high))
for (gene in rownames(cnv_p_signif)) {
  cnv_high_expr_i = merge(as.data.frame(t(Recode(cnv_high[gene,], "1='Amplication';0='Diploid';-1='Deletion'"))) %>% tibble::rownames_to_column('sample'),
                      as.data.frame(t(tcga_rna_tumor_mean[gene,])) %>% tibble::rownames_to_column('sample'),
                      by='sample')
  cnv_high_expr = merge(cnv_high_expr, cnv_high_expr_i, by='sample')
}
#*gene=CFB
cps_cnv = t(combn(unique(cnv_high_expr$CFB.x), 2))
cps_cnv = split(cps_cnv,list(cps_cnv[,1],cps_cnv[,2]), drop = T)
cnv_high_expr$CFB.x = factor(cnv_high_expr$CFB.x, levels = c('Amplication', 'Diploid', 'Deletion'), ordered = T)
ggviolin(cnv_high_expr, 'CFB.x', 'CFB.y', color = 'black', fill = 'CFB.x',
         palette = c("red", "orange", "blue"),
         add = c('boxplot'), add.params = list(fill='white')) +
  theme(axis.title.x = element_blank(),
        legend.position = 'none') +
  labs(y=paste0('CFB mRNA')) +
  stat_compare_means(comparisons = cps_cnv, label = 'p.signif')
#*gene=CLIC2
cps_cnv = t(combn(unique(cnv_high_expr$CLIC2.x), 2))
cps_cnv = split(cps_cnv,list(cps_cnv[,1],cps_cnv[,2]), drop = T)
cnv_high_expr$CLIC2.x = factor(cnv_high_expr$CLIC2.x, levels = c('Amplication', 'Diploid', 'Deletion'), ordered = T)
ggviolin(cnv_high_expr, 'CLIC2.x', 'CLIC2.y', color = 'black', fill = 'CLIC2.x',
         palette = c("red", "orange", "blue"),
         add = c('boxplot'), add.params = list(fill='white')) +
  theme(axis.title.x = element_blank(),
        legend.position = 'none') +
  labs(y=paste0('CLIC2 mRNA')) +
  stat_compare_means(comparisons = cps_cnv, label = 'p.signif')
#*gene=UBE2L6
cps_cnv = t(combn(unique(cnv_high_expr$UBE2L6.x), 2))
cps_cnv = split(cps_cnv,list(cps_cnv[,1],cps_cnv[,2]), drop = T)
cnv_high_expr$UBE2L6.x = factor(cnv_high_expr$UBE2L6.x, levels = c('Amplication', 'Diploid', 'Deletion'), ordered = T)
ggviolin(cnv_high_expr, 'UBE2L6.x', 'UBE2L6.y', color = 'black', fill = 'UBE2L6.x',
         palette = c("red", "orange", "blue"),
         add = c('boxplot'), add.params = list(fill='white')) +
  theme(axis.title.x = element_blank(),
        legend.position = 'none') +
  labs(y=paste0('UBE2L6 mRNA')) +
  stat_compare_means(comparisons = cps_cnv, label = 'p.signif')

##3.高低TIPRGPI组间hallmark通路富集差异（GSVA）、GSEA验证GSVA一致结果通路、K-M曲线评估一致通路的OS预后价值
#*GSVA
gsvaSet_hallmark = getGmt('D:/GSEA/msigdb_v7.5_GMTs/h.all.v7.5.symbols.gmt')
gsvaEs_hallmark = gsva(expr=as.matrix(tcga_rna_tumor_mean), gset.idx.list=gsvaSet_hallmark, parallel.sz=16, 
                       method = "ssgsea",
                       kcdf = "Poisson",
                       abs.ranking = TRUE,
                       verbose = FALSE)
gsvaEs_hallmark2 = as.data.frame(t(gsvaEs_hallmark[,train_lasso$sample])) %>% mutate('sample'=rownames(.)) %>%
  merge(train_lasso[,c('sample', 'cluster', 'OS', 'OS.time')], ., by='sample')
rownames(gsvaEs_hallmark2) = gsvaEs_hallmark2$sample
ttest_gsva = data.frame()
for (pathway in rownames(gsvaEs_hallmark)) {
  aa_result = t.test(gsvaEs_hallmark2 %>% filter(cluster=='High TIPRGPI') %>% pull(pathway),
                     gsvaEs_hallmark2 %>% filter(cluster=='Low TIPRGPI') %>% pull(pathway))
  ttest_gsva = rbind(ttest_gsva,
                     data.frame('statistics'=aa_result$statistic,
                                'Pvalue'=aa_result$p.value))
}
rownames(ttest_gsva) = rownames(gsvaEs_hallmark)
ttest_gsva = ttest_gsva %>% arrange(Pvalue) %>% tibble::rownames_to_column('pathway')
ttest_gsva_signif = ttest_gsva %>% filter(Pvalue<0.05)
ggplot(ttest_gsva_signif, aes(x=reorder(pathway, statistics), y=statistics, fill=statistics>0)) +
  geom_bar(stat = 'identity') +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = 'transparent'),
        legend.position = 'none') +
  scale_fill_manual(values = c('cyan4', 'red')) +
  coord_flip() +
  ylim(-20, 10) +
  labs(x=NULL, y='t value of GSVA score,\n High TIPRGPI vs Low TIPRGPI')
#*GSEA
gseaSet_hallmark = read.gmt('D:/GSEA/msigdb_v7.5_GMTs/h.all.v7.5.symbols.gmt')
compare_2limma_gsea = function(a){
  rna_mean_i = tcga_rna_tumor_mean[,train_lasso$sample]
  group_i = ifelse(train_lasso %>% tibble::column_to_rownames('sample') %>% .[colnames(rna_mean_i),] %>% .$cluster==a, 'high', 'low')
  design = model.matrix(~0+factor(group_i))
  colnames(design) = c('high', 'low')
  rownames(design) = colnames(rna_mean_i)
  #CS-CRS需要随分组更改
  contrast.matrix <- makeContrasts(contrasts = paste0(c('high','low'),collapse = '-'), levels=design)
  fit = lmFit(rna_mean_i, design)
  fit1 = contrasts.fit(fit, contrast.matrix)
  fit2 = eBayes(fit1)
  dif = topTable(fit2, coef=1, adjust="BH", number=ncol(rna_mean_i))
  dif_new = dif[(dif$adj.P.Val < 0.05 & abs(dif$logFC) >= log2(1.5)),]
  #write.table(dif_new, file=paste0('output/4.1_',a,'_',b,'_limma_dif.txt'), sep='\t', quote=F)
  return(list(dif, dif_new))
}
high_low = compare_2limma_gsea('High TIPRGPI')[1][[1]] #457个
high_low_dif = compare_2limma_gsea('High TIPRGPI')[2][[1]] #414个,0个up,13个down
high_low_dif = high_low_dif %>% mutate('gene'=rownames(.))
genelist_gsea = 2**high_low_dif$logFC %>% sort(., decreasing = T) %>% "names<-"(., high_low_dif$gene)
gsea_ttest <- GSEA(genelist_gsea, TERM2GENE = gseaSet_hallmark)
gsea_ttest2 = gsea_ttest@result
intersect(ttest_gsva_signif$pathway, gsea_ttest2$ID)
enrichplot::gseaplot2(gsea_ttest, geneSetID = 1:2, color="red",pvalue_table = F,base_size=10,ES_geom="line", rel_heights = c(4,1,1))
#*生存曲线
gsvaEs_hallmark3 = gsvaEs_hallmark2 %>% select(sample, OS, OS.time, gsea_ttest2$ID)
ann_colors = c(ann_colors, list(cluster=c('High Score'='red', 'Low Score'='blue')))
for (pathway in gsea_ttest2$ID) {
  gsvaEs_hallmark3$cluster = ifelse(gsvaEs_hallmark3[,pathway]>median(gsvaEs_hallmark3[,pathway]), 'High Score', 'Low Score')
  pdf(paste0("output/5.2_survival_", pathway, ".pdf"), width=6, height=5, onefile = F)
  plot_surv(gsvaEs_hallmark3, c('High Score', 'Low Score'), 5, '')
  dev.off()
  tiff(paste0("output/5.2_survival_", pathway, ".tiff"), res=300, width=6, height=5, compression="lzw", units="in")
  plot_surv(gsvaEs_hallmark3, c('High Score', 'Low Score'), 5, '')
  dev.off()
}

##Section6：TIPRGPI关联免疫过程
##1.高低TIPRGPI风险分组间肿瘤微环境细胞浸润差异（区分适应性免疫、固有免疫细胞、基质细胞）
#*SSGSEA
write.table(tcga_rna_tumor_mean, 'output/6.1_tcga_rna_tumor_mean.txt', col.names = T, row.names = T, sep = '\t', quote = F)
write.table(es, 'output/6.1_es.txt', col.names = T, row.names = T, sep = '\t', quote = F)
es_result2 = merge(es %>% mutate(., 'sample'=rownames(.)),
                   train_lasso[,c('sample', 'cluster')], by='sample')
es_result3 = melt(es_result2, id.vars = c('sample', 'cluster'), variable.name = 'cell', value.name = 'ES')
cell_class = list('Stromal'=c('Endothelial cells', 'Fibroblasts'),
                  'Innate'=c('Macrophages M0', 'Macrophages M1', 'Macrophages M2', 'Eosinophils', 'Dendritic cells activated', 'Dendritic cells resting', 'Mast cells activated', 'Mast cells resting', 'Monocytes', 'Neutrophils', 'NK cells activated', 'NK cells resting', 'Plasma cells'),
                  'Adaptive'=c('B cells memory', 'B cells naive', 'T cells CD4 memory activated', 'T cells CD4 memory resting', 'T cells CD4 naive', 'T cells CD8', 'T cells follicular helper', 'T cells gamma delta', 'T cells regulatory (Tregs)'))
cell_class = data.frame('cell'=unlist2(cell_class),
                        'type'=names(unlist2(cell_class)))
rownames(cell_class) = cell_class$cell
es_result3 = merge(es_result3, cell_class, by='cell')
es_result3$cell = factor(es_result3$cell, levels = cell_class$cell, ordered = T)
ggplot(es_result3, aes(x=cell,y=ES, fill=cluster))+
  geom_violin() +
  geom_boxplot(width=0.1, position = position_dodge(0.9))+
  theme_bw()+
  theme(legend.position="top",
        axis.text.x=element_text(colour=Recode(cell_class$type, "'Stromal'='blue';'Innate'='green4';'Adaptive'='red'"),family="Times",size=12, angle = 45, vjust=1, hjust = 1),
        axis.text.y=element_text(family="Times",size=14,face="plain"),
        axis.title.y=element_text(family="Times",size = 14,face="plain"),
        axis.title.x=element_text(family="Times",size = 14,face="plain"),
        plot.title = element_text(family="Times",size=15,face="bold",hjust = 0.5,),
        axis.text.x.bottom = element_text(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  xlab("")+ylab("Cell Infiltration") +
  scale_fill_manual(values = c("red", "blue")) +
  stat_compare_means(aes(group = cluster), label='p.signif')

##2.计算TIPRGPI与TME细胞的相关性分析
es_score = merge(es %>% mutate('sample'=rownames(.)),
                 train_lasso[,c('sample', 'TIPRGPI', 'OS', 'OS.time')], by='sample')
es_score_cor = rcorr(es_score %>% select(-c('sample')) %>% as.matrix())
es_score_p = data.frame('cell'=colnames(es),
                        'p'=es_score_cor$P['TIPRGPI', colnames(es)])
es_score_p = es_score_p %>% filter(p<0.05)
es_score_r = data.frame('cell'=colnames(es),
                        'r'=es_score_cor$r['TIPRGPI', colnames(es)])
es_score_cor2 = merge(es_score_p, es_score_r, by='cell')
es_score_cor2 = merge(es_score_cor2, cell_class, by='cell')
es_score_cor2$cell = factor(es_score_cor2$cell, levels = cell_class$cell, ordered = T)
ggplot(es_score_cor2 %>% rename('correlation'='r'), aes(x=cell, y=correlation)) +
  geom_bar(stat = 'identity', width = 0.1) +
  geom_point(aes(size=abs(correlation), color=p)) +
  scale_color_continuous(low='red', high='blue') +
  theme_bw() +
  coord_flip() +
  theme(axis.text.y=element_text(colour=Recode(cell_class$type[-1], "'Stromal'='blue';'Innate'='green4';'Adaptive'='red'"),family="Times", size = 11),
        axis.text.x=element_text(size = 11),
        axis.title.x = element_text(size = 12)) +
  labs(x=NULL)

##3.TIPRGPI与TME最相关的细胞筛选
#*生存相关的细胞类型
cell_cox_result = data.frame()
for (i in colnames(es)) {
  print(i)
  es_score$'x' = es_score[,i]
  model = coxph(Surv(OS.time, OS) ~ x, data=es_score[!is.na(es_score[,i]),])
  print(c(colnames(es_score)[i], cox.zph(model)$table))
  coef = cbind(summary(model)$coefficients, summary(model)$conf.int)
  #rownames(coef)[1] = colnames(pheno15)[i]
  cell_cox_result = rbind(cell_cox_result, coef) #254基因*9列cox属性
}
rownames(cell_cox_result) = colnames(es)
write.table(cell_cox_result, file='output/6.3_cell_cox_result.txt', row.names=T, col.names=T, quote=F, sep='\t')
#*ssGSEA差异分析、相关分析、生存分析基因统计venn图
cells_dif = setdiff(colnames(es), 'Fibroblasts')
cells_relative = es_score_p$cell
cells_surv = rownames(cell_cox_result)[cell_cox_result$`Pr(>|z|)`<0.05]
ggvenn(list('Differential analysis'=cells_dif,'Correlation analysis'=cells_relative,'Survival analysis'=cells_surv), show_percentage = F, set_name_size = 4.3, fill_color = c('green4', 'blue', 'red'))
cells_stat = data.frame('type'=c('Differential analysis','Correlation analysis','Survival analysis'),
                        'num'=c(length(cells_dif), length(cells_relative), length(cells_surv)))
cells_stat$type = factor(cells_stat$type, levels = cells_stat$type, ordered = T)
ggplot(cells_stat, aes(x=type, y=num)) +
  geom_bar(stat = 'identity', aes(fill=type)) +
  theme(panel.background = element_blank(),
        axis.line = element_line(),
        legend.position = 'none',
        plot.title = element_text(hjust = 0.5, face = 'italic')) +
  labs(title = 'Size of each list', x=NULL) +
  scale_fill_manual(values = c('green4', 'blue', 'red'))
cells_top = Reduce(intersect, list(cells_dif, cells_relative, cells_surv))
  
##4.TIPRGPI组成基因与TME最相关的细胞相关性分析（区分适应性免疫、固有免疫细胞、基质细胞，第二张图可放附图）
es_hubs_cor = rcorr(t(tcga_rna_tumor_mean[hubs_lasso,]), as.matrix(es[,cells_top]))es_hubs_cor2 = es_hubs_cor
es_hubs_cor2$r = es_hubs_cor2$r[hubs_lasso, cells_top]
es_hubs_cor2$P = es_hubs_cor2$P[hubs_lasso, cells_top]
es_hubs_cor2$P = ifelse(es_hubs_cor2$P<=0.0001, '****', ifelse(es_hubs_cor2$P<=0.001, '***', ifelse(es_hubs_cor2$P<0.01, '**', ifelse(es_hubs_cor2$P<0.05, '*', ''))))
ph_es_hubs_cor = pheatmap::pheatmap(as.data.frame(es_hubs_cor2$r), color = colorRampPalette(colors = c("blue","white","red"))(100), display_numbers = es_hubs_cor2$P, cluster_rows = F, cluster_cols = F, angle_col = '45', breaks = unique(c(seq(-3,3, length=100))), border_color = 'white')
ph_es_hubs_cor$gtable$grobs[[2]]$gp = gpar(col=Recode(cell_class[cells_top,]$type, "'Stromal'='blue';'Innate'='green4';'Adaptive'='red'"))
ph_es_hubs_cor

##5.TIPRGPI风险高低组间不同的免疫基因集表达差异
pheatmap(tcga_rna_tumor_mean[immunecell_genes, annotation_col_tcga$sample], cluster_row = F, cluster_col = F, scale = 'row', show_colnames = F, annotation_col=annotation_col_tcga[c('TIPRGPI')], color = colorRampPalette(colors = c("royalblue2","white","red2"))(100), show_rownames=T, annotation_colors = ann_colors[3], border_color = NA, breaks = unique(c(seq(-2, 2, length=100))), gaps_col = length(tcga_tumor_ids_high))

pheatmap(tcga_rna_tumor_mean[immunesuppression_genes, annotation_col_tcga$sample], cluster_row = F, cluster_col = F, scale = 'row', show_colnames = F, annotation_col=annotation_col_tcga[c('TIPRGPI')], color = colorRampPalette(colors = c("royalblue2","white","red2"))(100), show_rownames=T, annotation_colors = ann_colors[3], border_color = NA, breaks = unique(c(seq(-2, 2, length=100))), gaps_col = length(tcga_tumor_ids_high))

pheatmap(tcga_rna_tumor_mean[setdiff(innateimmunity_genes, 'IFNA1'), annotation_col_tcga$sample], cluster_row = F, cluster_col = F, scale = 'row', show_colnames = F, annotation_col=annotation_col_tcga[c('TIPRGPI')], color = colorRampPalette(colors = c("royalblue2","white","red2"))(100), show_rownames=T, annotation_colors = ann_colors[3], border_color = NA, breaks = unique(c(seq(-2, 2, length=100))), gaps_col = length(tcga_tumor_ids_high))

pheatmap(tcga_rna_tumor_mean[adaptiveimmunity_genes, annotation_col_tcga$sample], cluster_row = F, cluster_col = F, scale = 'row', show_colnames = F, annotation_col=annotation_col_tcga[c('TIPRGPI')], color = colorRampPalette(colors = c("royalblue2","white","red2"))(100), show_rownames=T, annotation_colors = ann_colors[3], border_color = NA, breaks = unique(c(seq(-2, 2, length=100))), gaps_col = length(tcga_tumor_ids_high))

pheatmap(tcga_rna_tumor_mean[setdiff(app_genes, c("HLA-DQA2", "LILRA3")), annotation_col_tcga$sample], cluster_row = F, cluster_col = F, scale = 'row', show_colnames = F, annotation_col=annotation_col_tcga[c('TIPRGPI')], color = colorRampPalette(colors = c("royalblue2","white","red2"))(100), show_rownames=T, annotation_colors = ann_colors[3], border_color = NA, breaks = unique(c(seq(-2, 2, length=100))), gaps_col = length(tcga_tumor_ids_high))

pheatmap(tcga_rna_tumor_mean[cytotox_genes, annotation_col_tcga$sample], cluster_row = F, cluster_col = F, scale = 'row', show_colnames = F, annotation_col=annotation_col_tcga[c('TIPRGPI')], color = colorRampPalette(colors = c("royalblue2","white","red2"))(100), show_rownames=T, annotation_colors = ann_colors[3], border_color = NA, breaks = unique(c(seq(-2, 2, length=100))), gaps_col = length(tcga_tumor_ids_high))

pheatmap(tcga_rna_tumor_mean[setdiff(inflammation_genes, 'IL17A'), annotation_col_tcga$sample], cluster_row = F, cluster_col = F, scale = 'row', show_colnames = F, annotation_col=annotation_col_tcga[c('TIPRGPI')], color = colorRampPalette(colors = c("royalblue2","white","red2"))(100), show_rownames=T, annotation_colors = ann_colors[3], border_color = NA, breaks = unique(c(seq(-2, 2, length=100))), gaps_col = length(tcga_tumor_ids_high))

immune_signatures = data.frame('Gene.Signature'=c('Immune cells recruitment', 'Immune suppression genes', 'Innate immunity', 'Adaptive immunity', 'Antigen presentation and processing', 'Cytotoxicity/killing of cancer cells', 'Inflammation'),
                               'Platform'='Gene expression profiling',
                               'Genes'=c(paste0(immunecell_genes, collapse = ','), paste0(immunesuppression_genes, collapse = ','),
                                         paste0(setdiff(innateimmunity_genes, 'IFNA1'), collapse = ','), paste0(adaptiveimmunity_genes, collapse = ','),
                                         paste0(setdiff(app_genes, c("HLA-DQA2", "LILRA3")), collapse = ','), paste0(cytotox_genes, collapse = ','),
                                         paste0(setdiff(inflammation_genes, 'IL17A'), collapse = ',')))
write.table(immune_signatures, 'output/6.5_immune_signatures.txt', sep = '\t', quote = F, row.names = F, col.names = T)
rownames(immune_signatures) = immune_signatures$Gene.Signature
write.gmt <- function(df,file_path){
  sink(file_path)
  sapply(immune_signatures$Gene.Signature, function(i){
    cat(paste(c(immune_signatures[i,'Gene.Signature'], immune_signatures[i,'Platform'], strsplit2(immune_signatures[i,'Genes'], split = ',')),collapse='\t'))
    cat('\n')
  })
  sink()
}
write.gmt(immune_signatures, 'output/6.5_immune_signatures.gmt')
immune_signatures_gmt = getGmt('output/6.5_immune_signatures.gmt')
gsvaEs_immune_signatures = gsva(expr=as.matrix(tcga_rna_tumor_mean), gset.idx.list=immune_signatures_gmt, parallel.sz=32)
gsvaEs_immune_signatures2 = as.data.frame(t(gsvaEs_immune_signatures)) %>% mutate(., 'sample'=rownames(.))
gsvaEs_immune_signatures2 = merge(gsvaEs_immune_signatures2, train_lasso[,c('sample', 'cluster')], by='sample')
wilcox_results_gsva_immune = data.frame()
for (i in immune_signatures$Gene.Signature) {
  print(i)
  wilcox_result = wilcox.test(as.numeric(gsvaEs_immune_signatures2[gsvaEs_immune_signatures2$cluster=='High TIPRGPI', i]), 
                              as.numeric(gsvaEs_immune_signatures2[gsvaEs_immune_signatures2$cluster=='Low TIPRGPI', i]))
  wilcox_results_gsva_immune = rbind(wilcox_results_gsva_immune,
                                 data.frame('statistic'=wilcox_result$statistic,
                                            'pvalue'=wilcox_result$p.value,
                                            'item'=i))
}
rownames(wilcox_results_gsva_immune) = wilcox_results_gsva_immune$item
wilcox_results_gsva_immune = wilcox_results_gsva_immune[order(wilcox_results_gsva_immune$pvalue),]
write.table(wilcox_results_gsva_immune, 'output/6.5_wilcox_results_gsva_immune.txt', sep = '\t', quote = F, row.names = F, col.names = T)

##Section7：TIPRGPI可作为免疫治疗敏感预测因子
##1.高低TIPRGPI风险分组间免疫检查点的表达差异；干扰素-γ通路marker表达差异; 调控因子表达差异；IPS评分差异；TIDE与TIPRGPI相关性
#*icgs
tcga_rna_tumor_mean_icgs = tcga_rna_tumor_mean[icgs,]
tcga_rna_tumor_mean_icgs = tcga_rna_tumor_mean_icgs[complete.cases(tcga_rna_tumor_mean_icgs),]
tcga_rna_tumor_mean_icgs = merge(as.data.frame(t(tcga_rna_tumor_mean_icgs)) %>% tibble::rownames_to_column('sample'),
                                 train_lasso[,c('sample', 'cluster')], by='sample')
tcga_rna_tumor_mean_icgs2 = melt(tcga_rna_tumor_mean_icgs, id.vars = c('sample', 'cluster'), variable.name = 'gene', value.name = 'expr')
pdf(file='output/7.1_icg.pdf', width=20, height=5, onefile = F)
tiff("output/7.1_icg.tiff", res=300, width=20, height=5, compression="lzw", units="in")
ggplot(tcga_rna_tumor_mean_icgs2, aes(x=gene, y=expr, fill=cluster)) +
  geom_boxplot(outlier.color="white") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(x=NULL, y='expression level') +
  guides(fill=guide_legend(title = NULL)) +
  scale_fill_manual(values = ann_colors[[2]]) +
  stat_compare_means(aes(group = cluster), label='p.signif')
dev.off()
#*IFN
tcga_rna_tumor_mean_ifn = tcga_rna_tumor_mean[ifn_genes,]
tcga_rna_tumor_mean_ifn = tcga_rna_tumor_mean_ifn[complete.cases(tcga_rna_tumor_mean_ifn),]
tcga_rna_tumor_mean_ifn = merge(as.data.frame(t(tcga_rna_tumor_mean_ifn)) %>% tibble::rownames_to_column('sample'),
                                train_lasso[,c('sample', 'cluster')], by='sample')
tcga_rna_tumor_mean_ifn2 = melt(tcga_rna_tumor_mean_ifn, id.vars = c('sample', 'cluster'), variable.name = 'gene', value.name = 'expr')
ggplot(tcga_rna_tumor_mean_ifn2, aes(x=gene, y=expr, fill=cluster)) +
  geom_boxplot(outlier.color="white") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(x=NULL, y='expression level') +
  guides(fill=guide_legend(title = NULL)) +
  scale_fill_manual(values = ann_colors[[2]]) +
  stat_compare_means(aes(group = cluster), label='p.signif')
#*
tcga_rna_tumor_mean_ = tcga_rna_tumor_mean[,]
tcga_rna_tumor_mean_ = tcga_rna_tumor_mean_[complete.cases(tcga_rna_tumor_mean_),]
tcga_rna_tumor_mean_ = merge(as.data.frame(t(tcga_rna_tumor_mean_)) %>% tibble::rownames_to_column('sample'),
                                train_lasso[,c('sample', 'cluster')], by='sample')
tcga_rna_tumor_mean_2 = melt(tcga_rna_tumor_mean_, id.vars = c('sample', 'cluster'), variable.name = 'gene', value.name = 'expr')
ggplot(tcga_rna_tumor_mean_2, aes(x=gene, y=expr, fill=cluster)) +
  geom_boxplot(outlier.color="white") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(x=NULL, y='expression level') +
  guides(fill=guide_legend(title = NULL)) +
  scale_fill_manual(values = ann_colors[[2]]) +
  stat_compare_means(aes(group = cluster), label='p.signif')
##IPS评分差异（下载：https://tcia.at/patients）
ips_raw = read.delim('output/7.1_TCIA-ClinicalData.tsv', sep = '\t')
ips = merge(ips_raw[,c('barcode', 'ips_ctla4_neg_pd1_neg', 'ips_ctla4_neg_pd1_pos', 'ips_ctla4_pos_pd1_neg', 'ips_ctla4_pos_pd1_pos')],
            train_lasso[,c('sample', 'cluster')] %>% mutate('barcode'=substr(sample, 1, 12)),
            by='barcode')
ips = ips %>% rename('IPS'='ips_ctla4_neg_pd1_neg', 'IPS-CTLA4-and PD1/PDL1/PDL2 blocker'='ips_ctla4_pos_pd1_pos',
                     'IPS-CTLA4 blocker'='ips_ctla4_pos_pd1_neg', 'IPS-PD1/PDL1/PDL2 blocker'='ips_ctla4_neg_pd1_pos')
ggplot(ips, aes(cluster, IPS)) + 
  ggdist::stat_halfeye(aes(color=cluster,fill=cluster),adjust = .5, width = .3, .width = 0, justification = -.3, point_colour = NA) + 
  geom_boxplot(aes(color=cluster),width = .1, outlier.shape = NA) +
  gghalves::geom_half_point(aes(color=cluster),side = "l", range_scale = .4, alpha = .5) +
  ggsci::scale_color_nejm()+
  ggsci::scale_fill_nejm() +
  labs(x=NULL) +
  theme_bw() +
  theme(legend.position = 'none') +
  scale_fill_manual(values = ann_colors[[2]]) +
  scale_color_manual(values = ann_colors[[2]]) +
  stat_compare_means(comparisons = cps2, label='p.signif', label.y = 11)
ggplot(ips, aes(cluster, `IPS-CTLA4-and PD1/PDL1/PDL2 blocker`)) + 
  ggdist::stat_halfeye(aes(color=cluster,fill=cluster),adjust = .5, width = .3, .width = 0, justification = -.3, point_colour = NA) + 
  geom_boxplot(aes(color=cluster),width = .1, outlier.shape = NA) +
  gghalves::geom_half_point(aes(color=cluster),side = "l", range_scale = .4, alpha = .5) +
  ggsci::scale_color_nejm()+
  ggsci::scale_fill_nejm() +
  labs(x=NULL) +
  theme_bw() +
  theme(legend.position = 'none') +
  scale_fill_manual(values = ann_colors[[2]]) +
  scale_color_manual(values = ann_colors[[2]]) +
  stat_compare_means(comparisons = cps2, label='p.signif', label.y = 11)
ggplot(ips, aes(cluster, `IPS-CTLA4 blocker`)) + 
  ggdist::stat_halfeye(aes(color=cluster,fill=cluster),adjust = .5, width = .3, .width = 0, justification = -.3, point_colour = NA) + 
  geom_boxplot(aes(color=cluster),width = .1, outlier.shape = NA) +
  gghalves::geom_half_point(aes(color=cluster),side = "l", range_scale = .4, alpha = .5) +
  ggsci::scale_color_nejm()+
  ggsci::scale_fill_nejm() +
  labs(x=NULL) +
  theme_bw() +
  theme(legend.position = 'none') +
  scale_fill_manual(values = ann_colors[[2]]) +
  scale_color_manual(values = ann_colors[[2]]) +
  stat_compare_means(comparisons = cps2, label='p.signif', label.y = 11)
ggplot(ips, aes(cluster, `IPS-PD1/PDL1/PDL2 blocker`)) + 
  ggdist::stat_halfeye(aes(color=cluster,fill=cluster),adjust = .5, width = .3, .width = 0, justification = -.3, point_colour = NA) + 
  geom_boxplot(aes(color=cluster),width = .1, outlier.shape = NA) +
  gghalves::geom_half_point(aes(color=cluster),side = "l", range_scale = .4, alpha = .5) +
  ggsci::scale_color_nejm()+
  ggsci::scale_fill_nejm() +
  labs(x=NULL) +
  theme_bw() +
  theme(legend.position = 'none') +
  scale_fill_manual(values = ann_colors[[2]]) +
  scale_color_manual(values = ann_colors[[2]]) +
  stat_compare_means(comparisons = cps2, label='p.signif', label.y = 11)
#TIDE
tide2 = merge(train_lasso[,c('sample', 'TIPRGPI', 'cluster')], tide[,c('Patient', 'TIDE')], by.y='Patient', by.x='sample')
tide_scatter = ggscatter(tide2, x = "TIPRGPI", y = 'TIDE',
          color = "cluster", size = 3,# 点的颜色与大小
          add = "reg.line",  # 添加回归线
          add.params = list(color = "blue", fill = "gray"), # 回归线的颜色设置为红色，区间颜色设置为灰色
          conf.int = TRUE, # 添加回归线的置信区间
          cor.coef = TRUE, # 添加相关系数
          cor.coeff.args = list(method = "pearson", label.sep = "\t"),#选择Pearson相关,
          margin.params = list(fill = "blue3"),
          ggtheme = theme_bw(),
          xlab = 'TIPRGPI', ylab = 'TIDE',
          repel = F) +
  rremove('legend') +
  scale_color_manual(values = ann_colors[[2]])
tide_box = ggplot(tide2, aes(x = cluster, y = TIDE)) +
  geom_boxplot(aes(fill=cluster), outlier.color = 'black') +
  scale_fill_manual(values = ann_colors[[2]]) +
  theme(panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()) +
  #rremove('legend') +
  stat_compare_means(comparisons = cps2, label='p.signif') +
  guides(fill=guide_legend(title = NULL))
cowplot::plot_grid(tide_scatter, tide_box, ncol=2, align = 'hv', rel_widths = c(1,1))
  
##Section8：TIPRGPI相关核心靶点识别和候选分子预测
##1.高低TIPRGP组筛选差异基因构建PPI互作网络根据节点度筛选核心因子
high_low_dif2 = high_low_dif %>% filter(adj.P.Val<0.01 & abs(logFC)>1.5)
write.table(rownames(high_low_dif2), file='output/8.1_high_low_dif2_genes.txt', row.names=F, col.names=F, quote=F, sep='\t')
ppi_out = read.table('output/8.1_string_interactions_short.tsv', sep = '\t', header = T)












































