#BSUB -q normal
#BSUB -J tumor1_gatk
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -R "span[ptile=36]"
#BSUB -n 96

cd /share/home/guoguoji/RAWDATA/YUCX/20180711tumor
 
/share/home/guoguoji/tools/STAR-2.5.2a/bin/Linux_x86_64_static/STAR \
	--runThreadN 20 --genomeDir /share/home/guoguoji/Desktop/xyy/zebrafish_reference.10/genomeDir \
	--twopassMode Basic --readFilesIn out_H_R2.fastq --outSAMtype BAM SortedByCoordinate \
	--outFileNamePrefix ./ 

java -jar /share/home/guoguoji/tools/picard-tools-1.119/AddOrReplaceReadGroups.jar \
	I=Aligned.sortedByCoord.out.bam O=_RG_added_sorted.bam \
	SO=coordinate RGID=tumor1 RGLB=rna RGPL=illumina RGPU=hiseq RGSM=tumor1

java -jar /share/home/guoguoji/tools/picard-tools-1.119/MarkDuplicates.jar \
	I=_RG_added_sorted.bam O=_dedup.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=_dedup.metrics

java -jar /share/home/guoguoji/Desktop/callsnptools/GenomeAnalysisTK.jar -T SplitNCigarReads -R /share/home/guoguoji/Desktop/xyy/zebrafish_reference.10/Danio_rerio.GRCz10.dna.toplevel.fa -I _dedup.bam -o _dedup_split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS

java -jar /share/home/guoguoji/Desktop/callsnptools/GenomeAnalysisTK.jar -T RealignerTargetCreator -R /share/home/guoguoji/Desktop/xyy/zebrafish_reference.10/Danio_rerio.GRCz10.dna.toplevel.fa -I _dedup_split.bam -o _realign_interval.list 

java -jar /share/home/guoguoji/Desktop/callsnptools/GenomeAnalysisTK.jar -T IndelRealigner -R /share/home/guoguoji/Desktop/xyy/zebrafish_reference.10/Danio_rerio.GRCz10.dna.toplevel.fa -I _dedup_split.bam -o _realign.bam -targetIntervals _realign_interval.list

java -jar /share/home/guoguoji/Desktop/callsnptools/GenomeAnalysisTK.jar -T BaseRecalibrator -R /share/home/guoguoji/Desktop/xyy/zebrafish_reference.10/Danio_rerio.GRCz10.dna.toplevel.fa -I _realign.bam --knownSites /share/home/guoguoji/Desktop/callsnptools/zebrafish_reference/danio_rerio.vcf.gz -o recal_data.table

java -jar /share/home/guoguoji/Desktop/callsnptools/GenomeAnalysisTK.jar -T PrintReads -R /share/home/guoguoji/Desktop/xyy/zebrafish_reference.10/Danio_rerio.GRCz10.dna.toplevel.fa -I _realign.bam -BQSR recal_data.table -o _BQSR.bam

java -jar /share/home/guoguoji/Desktop/callsnptools/GenomeAnalysisTK.jar -T MuTect2 -I:tumor /share/home/guoguoji/RAWDATA/JUNQING/sandai/WT/_recal_data.table/_BQSR.bam --dbsnp /share/home/guoguoji/Desktop/callsnptools/human_reference/dbsnp_138.hg38.vcf.gz --output_mode EMIT_VARIANTS_ONLY -o _mutect2.vcf.gz -L 5 -R /share/home/guoguoji/tools/STAR_Reference_Human/Homo_sapiens.GRCh38.fa
