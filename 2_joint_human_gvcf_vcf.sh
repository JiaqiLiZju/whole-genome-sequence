#!/bin/bash

##define parameters
sample_ID=$1
samples=$(echo $sample | tr "," "\n")
in_dir=$out_dir
out_name=$2

##mkdir
if [ ! -d "$out_dir/population"  ]; then
    mkdir -p $out_dir/population
fi

###joint genotyping
sample_gvcfs=""
for sample in $samples; do
	sample_gvcfs=${sample_gvcfs}"-V $out_dir/$sample.gvcf.gz \\"\n
done
time $gatk CombineGVCFs \
	-R $reference ${sample_gvcfs} \
	-O $out_dir/population/$out_name.gvcf.gz \
	&& echo "$out_name.gvcf done"
time $gatk GenotypeGVCFs \
	-R $reference \
	-V $out_dir/population/$out_name.gvcf.gz \
	-O $out_dir/population/$out_name.vcf.gz \

###VQSR
##SNP model
time $gatk variantRecalibrator \
	-R $reference \
	-V $out_dir/population/$out_name.vcf.gz \
	-resource:hapmap,known=false,training=true,truth=true,prior=15.0 $GATK_bundle/hapmap_3.3.hg38.vcf \
	-resource:omini,known=false,training=true,truth=false,prior=12.0 $GATK_bundle/1000G_omni2.5.hg38.vcf \
	-resource:1000G,known=false,training=true,truth=false,prior=10.0 $GATK_bundle/1000G_phase1.snps.high_confidence.hg38.vcf \
	-resource:dbsnp,known=true,training=false,truth=false,prior=6.0 $GATK_bundle/dbsnp_146.hg38.vcf \
	-an DP -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
	-mode SNP \
	-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 95.0 -tranche 90.0 \
	-rscriptFile $out_dir/population/$out_name.snp.plot.R \
	--tranches-file $out_dir/population/$out_name.snp.tranches \
	-O $out_dir/population/$out_name.snp.recal &&\
time $gatk ApplyVQSR \
	-R $reference \
	-V $out_dir/population/$out_name.vcf.gz \
	--ts_filter_level 99.0 \
	--tranches-file $out_dir/population/$out_name.snp.tranches \
	-recalFile $out_dir/population/$out_name.snp.recal \
	-mode SNP \
	-O $out_dir/population/$out_name.snp.VQSR.vcf.gz && echo "snp.VQSR done"

#Indel mode
time $gatk variantRecalibrator \
	-R $reference \
	-input $out_dir/population/$out_name.snp.VQSR.vcf.gz \
	-resource:mills,known=true,training=true,truth=true,prior=12.0 $GATK_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf \
	-an DP -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
	-mode INDEL \
	--max-gaussian 6 \
	-rscriptFile $out_dir/population/$out_name.indel.plot.R \
	--tranches-file $out_dir/population/$out_name.indel.tranches \
	-O $out_dir/population/$out_name.indel.recal && \
time $gatk ApplyVQSR \
	-R $reference \
	-input $out_dir/population/$out_name.snp.VQSR.vcf.gz \
	--ts_filter_level 99.0 \
	--tranches-file $out_dir/population/$out_name.indel.tranches \
	-recalFile $out_dir/population/$out_name.indel.recal \
	-mode INDEL \
	-O $out_dir/population/$out_name.indel.VQSR.vcf.gz && echo "indel VQSR done"
