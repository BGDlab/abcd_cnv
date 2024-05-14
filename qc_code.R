###########################
# sample and CNV level qc #
###########################
load('/shaz/abcd_cnv/preprocess/cnv_rawdata_merge.RData')
load('/shaz/abcd_cnv/preprocess/qc_rawdata.RData')

# remove CNVs on ChrX
cnv_rawdata_merge=cnv_rawdata_merge[-which(cnv_rawdata_merge$Chr==23),]

# remove samples from batch461. as there is something wrong with genotyping
batch=read.table('/shaz/abcd_cnv/preprocess/ABCD_release3.0_.batch_info.txt',header=T)
batch461=batch$abcd.id_redcap [which(batch$Axiom_Plate==461)]
cnv_rawdata_merge=cnv_rawdata_merge[-which(match(cnv_rawdata_merge$sub_id,batch461)>0),]

# exclude individuals with lrr_sd>0.35 or baf_sd>0.1 or |wf|>0.05
tem2=which(qc_rawdata$lrr_sd>0.35 | qc_rawdata$baf_sd>0.1 | abs(qc_rawdata$wf)>0.05)
tem3=qc_rawdata$sub_id[tem2]

tem4=which(match(cnv_rawdata_merge$sub_id,tem3)>0)

cnv_rawdata_merge=cnv_rawdata_merge[-tem4,]

# exclude individuals with the number of CNVs>=500
tem5=unique(cnv_rawdata_merge$SampleID)

tem6=c()
for (i in 1:length(tem5)) {
  tem6[i]=length(which(cnv_rawdata_merge$FID==tem5[i]))
}

tem7=tem5[which(tem6>=500)] 

tem8=which(match(cnv_rawdata_merge$FID,tem7)>0)
cnv_rawdata_merge=cnv_rawdata_merge[-tem8,]

# we only include CNVs with 70% overlapping between two algorithms, only include autosomal CNVs
cnv=cnv_rawdata_merge[ which(cnv_rawdata_merge$X.Two.Algs>70), ]

# exclude CNVs<50kb
cnv=cnv[-which(cnv$Size<50000),]

# exclude CNVs with confidence score<30 and numsnp<20
cnv=cnv[-which(cnv$Max.Score<30 | cnv$Max.SNP<20),]

# exclude CNVs 50% overlapping with segmental duplication
load('/shaz/abcd_cnv/preprocess/segmental_duplication_hg19.RData')

library(IRanges)
cnv_bind=c()
for (i in 1:22) {
  tem10=seg_dup[which(seg_dup$chrom==i),2]
  tem11=seg_dup[which(seg_dup$chrom==i),3]
  tem12=IRanges(start=tem10,end=tem11)
  
  cnv_index=which(cnv$Chr==i)
  tem13=cnv[which(cnv$Chr==i),4]
  tem14=cnv[which(cnv$Chr==i),5]
  tem15=IRanges(start=tem13,end=tem14)
  
  tem16=as.matrix(findOverlaps(tem15, tem12))
  
  tem17=c()
  for (f in 1:nrow(tem16)) {
    index1=tem16[f,1] 
    index2=tem16[f,2] 
    
    start1=tem13[index1] 
    end1=tem14[index1] 
    
    start2=tem10[index2] 
    end2=tem11[index2] 
    
    overlap=min(end1,end2)-max(start1,start2)
    cnv_length=end1-start1
    ref_length=end2-start2
    tem17=c(tem17,min(overlap/cnv_length,overlap/ref_length))
  }
  
  tem18=matrix(tem16[which(tem17 >= .5),],ncol = 2)
  
  bad.rows.chr=unique(tem18[,1])
  
  cnv_bad=cnv_index[bad.rows.chr]
  
  cnv_bind=c(cnv_bind,cnv_bad)
}

cnv=cnv[-cnv_bind,]

# exclude CNVs 50% overlapping with centrometic regions
load('/shaz/abcd_cnv/preprocess/centromeric_regions_hg19.RData')

cnv_bind=c()
for (i in 1:22) {
  tem10=cen_reg[which(cen_reg$chr==i),2]
  tem11=cen_reg[which(cen_reg$chr==i),3]
  tem12=IRanges(start=tem10,end=tem11)
  
  cnv_index=which(cnv$Chr==i)
  tem13=cnv[which(cnv$Chr==i),4]
  tem14=cnv[which(cnv$Chr==i),5]
  tem15=IRanges(start=tem13,end=tem14)
  
  tem16=as.matrix(findOverlaps(tem15, tem12))
  
  if (nrow(tem16)!=0) {
  
  tem17=c()
  for (f in 1:nrow(tem16)) {
    index1=tem16[f,1] 
    index2=tem16[f,2] 
    
    start1=tem13[index1] 
    end1=tem14[index1] 
    
    start2=tem10[index2] 
    end2=tem11[index2] 
    
    overlap=min(end1,end2)-max(start1,start2)
    cnv_length=end1-start1
    ref_length=end2-start2
    tem17=c(tem17,min(overlap/cnv_length,overlap/ref_length))
  }
  
  tem18=matrix(tem16[which(tem17 >= .5),],ncol = 2)
  
  bad.rows.chr=unique(tem18[,1])
  
  cnv_bad=cnv_index[bad.rows.chr]
  
  cnv_bind=c(cnv_bind,cnv_bad)
  }
}

cnv=cnv[-cnv_bind,]

# exclude CNVs 50% overlapping with telomeric regions
load('/shaz/abcd_cnv/preprocess/telomeric_regions_hg19.RData')

cnv_bind=c()
for (i in 1:22) {
  tem10=tel_reg[which(tel_reg$chr==i),2]
  tem11=tel_reg[which(tel_reg$chr==i),3]
  tem12=IRanges(start=tem10,end=tem11)
  
  cnv_index=which(cnv$Chr==i)
  tem13=cnv[which(cnv$Chr==i),4]
  tem14=cnv[which(cnv$Chr==i),5]
  tem15=IRanges(start=tem13,end=tem14)
  
  tem16=as.matrix(findOverlaps(tem15, tem12))
  
  if (nrow(tem16)!=0) {
    
    tem17=c()
    for (f in 1:nrow(tem16)) {
      index1=tem16[f,1] 
      index2=tem16[f,2] 
      
      start1=tem13[index1] 
      end1=tem14[index1] 
      
      start2=tem10[index2] 
      end2=tem11[index2] 
      
      overlap=min(end1,end2)-max(start1,start2)
      cnv_length=end1-start1
      ref_length=end2-start2
      tem17=c(tem17,min(overlap/cnv_length,overlap/ref_length))
    }
    
    tem18=matrix(tem16[which(tem17 >= .5),],ncol = 2)
    
    bad.rows.chr=unique(tem18[,1])
    
    cnv_bad=cnv_index[bad.rows.chr]
    
    cnv_bind=c(cnv_bind,cnv_bad)
  }
}

cnv=cnv[-cnv_bind,]

# exclude CNVs 50% overlapping with MHC
cnv_bind=c()
tem10=28477797
tem11=33448354
tem12=IRanges(start=tem10,end=tem11)

cnv_index=which(cnv$Chr==6)
tem13=cnv[which(cnv$Chr==6),4]
tem14=cnv[which(cnv$Chr==6),5]
tem15=IRanges(start=tem13,end=tem14)

tem16=as.matrix(findOverlaps(tem15, tem12))
tem17=c()
for (f in 1:nrow(tem16)) {
  index1=tem16[f,1] 
  index2=tem16[f,2] 
  
  start1=tem13[index1] 
  end1=tem14[index1] 
  
  start2=tem10[index2] 
  end2=tem11[index2] 
  
  overlap=min(end1,end2)-max(start1,start2)
  cnv_length=end1-start1
  ref_length=end2-start2
  tem17=c(tem17,min(overlap/cnv_length,overlap/ref_length))
}

tem18=matrix(tem16[which(tem17 >= .5),],ncol = 2)

bad.rows.chr=unique(tem18[,1])

cnv_bad=cnv_index[bad.rows.chr]

cnv_bind=c(cnv_bind,cnv_bad)

cnv=cnv[-cnv_bind,]

# remove subjects that only have genetic data, but do not have phenotypic data
load('/shaz/abcd_cnv/preprocess/cnv_qced.RData')

pheno=read.table('/mnt/isilon/bgdlab_resnas03/Data/ABCD/release4/abcd_cbcls01.txt',header=T)
pheno=pheno[which(pheno$eventname=='baseline_year_1_arm_1'),]

test1=intersect(unique(cnv$sub_id),pheno$src_subject_id)
test2=setdiff(cnv$sub_id,test1)
test3=which(match(cnv$sub_id,test2)>0)

cnv_clean=cnv[-test3,]

# remove CNVs with the type of '0,1,3' or '1,3'
cnv_clean=cnv_clean[ -which(cnv_clean$Type=='0, 1, 3' | cnv_clean$Type=='1, 3'), ]

# for families, include only one member
fam=read.table('/mnt/isilon/bgdlab_resnas03/Data/ABCD/release4/acspsw03.txt',header=T)
fam=fam[which(fam$eventname=='baseline_year_1_arm_1'),]

cnv_clean$rel_family_id=fam$rel_family_id[ match(cnv_clean$sub_id,fam$src_subject_id) ]

# rel_relationship: 
# 0:single 
# 1:sibling 
# 2:twin 
# 3:triple
cnv_clean$rel_relationship=fam$rel_relationship[ match(cnv_clean$sub_id,fam$src_subject_id) ]

cnv_del=c()
for (i in 1:nrow(cnv_clean)) {
  if (cnv_clean$rel_relationship[i]>0) {
    tem1=which(cnv_clean$rel_family_id==cnv_clean$rel_family_id[i])
    tem2=cnv_clean$sub_id[tem1]
    
    if (length(unique(tem2))>1) {
      tem3=tem1[which(tem2!=unique(tem2)[1])]
      cnv_del=c(cnv_del,tem3)
    }
  }
}
cnv_clean=cnv_clean[-unique(cnv_del),]









