#/bin/bash

# star 1-pass index
##make STAR reference
STAR --runMode genomeGenerate --runThreadN 12 \
	--genomeDir reference/star_index_1pass/ --genomeFastaFiles reference/tea_genome.fa \
	--sjdbGTFfile tea.gtf --sjdbOverhang 99 --genomeChrBinNbits 15

##alignment
STAR --runThreadN 12 --twopassMode Basic \
	--genomeDir ../reference/ \
	--readFilesIn ../Rawdata/Tea_C413-T01_good_1.fq.gz ../Rawdata/Tea_C413-T01_good_2.fq.gz --readFilesCommand zcat \
	--outSAMtype BAM SortedByCoordinate --outFileNamePrefix ./

# star 2-pass index
STAR --runThreadN 8 --runMode genomeGenerate \
	--genomeDir reference/star_index_2pass/ \
	--genomeFastaFiles reference/tea_genome.fa \
	--sjdbFileChrStartEnd T1/SJ.out.tab T2/SJ.out.tab T3/SJ.out.tab T4/SJ.out.tab 
	
# star 2-pass align
STAR --runThreadN 8 --genomeDir ../reference/star_index_2pass/ \
	--readFilesIn ../Rawdata/Tea_C413-T03_good_1.fq.gz ../Rawdata/Tea_C413-T03_good_2.fq.gz --readFilesCommand zcat \
	--outSAMtype BAM SortedByCoordinate --outFileNamePrefix ./