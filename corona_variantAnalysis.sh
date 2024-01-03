#store reference genome index in ref_genome variable
ref_genome=$1

##create directories named as bam, mark, raw_vcf filter_vcf and annotated_vcf
mkdir bam markduplicate raw_vcf filter_vcf annotated_vcf

##
#Run a for loop (all the filtered fastq file should be present in QC folder)
for f in QC/*R1*.fastq.gz
do
#Run below command to map filter reads on reference genome (It is assumed that index file of refence genome is alredy available)
#if index file is not available use this command to generate index file (bwa index  NC_045512.2.fasta)

bwa mem -t 10 -Y -R '@RG\tID:$(basename $f _R1.clean.fastq.gz)\tLB:$(basename $f _R1.clean.fastq.gz)\tPL:ILLUMINA\tPM:Novaseq\tSM:$(basename $f _R1.clean.fastq.gz)' $ref_genome $f QC/$(basename $f _R1.clean.fastq.gz)_R2.clean.fastq.gz | samtools view -@ 5 -1 - > bam/$(basename $f _R1.clean.fastq.gz).aligned.bam 2>bam/$(basename $f _R1.clean.fastq.gz).aligned.log

##sorting alignment file ###

samtools sort -O bam -@ 10 -o bam/$(basename $f _R1.clean.fastq.gz).aligned.sort.bam bam/$(basename $f _R1.clean.fastq.gz).aligned.bam

##mark duplicates###

gatk MarkDuplicates --INPUT bam/$(basename $f _R1.clean.fastq.gz).aligned.sort.bam --OUTPUT markduplicate/$(basename $f _R1.clean.fastq.gz).aligned.dedup.bam --METRICS_FILE markduplicate/$(basename $f _R1.clean.fastq.gz).matrix.txt --VALIDATION_STRINGENCY SILENT --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 --ASSUME_SORT_ORDER "queryname"

##collect alignment and insert size metrices

java -jar picard.jar CollectAlignmentSummaryMetrics R=$ref_genome I=markduplicate/$(basename $f _R1.clean.fastq.gz).aligned.dedup.bam O=markduplicate/$(basename $f _R1.clean.fastq.gz).alignment_metrices.txt

java -jar picard.jar CollectInsertSizeMetrics INPUT=markduplicate/$(basename $f _R1.clean.fastq.gz).aligned.dedup.bam OUTPUT=markduplicate/$(basename $f _R1.clean.fastq.gz).alignment_metrices.txt HISTOGRAM_FILE=markduplicate/$(basename $f _R1.clean.fastq.gz).insert_size_histogram.pdf

samtools depth -a markduplicate/$(basename $f _R1.clean.fastq.gz).aligned.dedup.bam > markduplicate/$(basename $f _R1.clean.fastq.gz).depth_out.txt

### index dedup bam file ####
samtools index markduplicate/$(basename $f _R1.clean.fastq.gz).aligned.dedup.bam

###variant calling ##
gatk HaplotypeCaller -R $ref_genome -I markduplicate/$(basename $f _R1.clean.fastq.gz).aligned.dedup.bam -O raw_vcf/$(basename $f _R1.clean.fastq.gz).raw_variants.vcf

##Filter variants (see DQ_filter tag in vcf) ##
gatk VariantFiltration -R $ref_genome -V raw_vcf/$(basename $f _R1.clean.fastq.gz).raw_variants.vcf -O filter_vcf/$(basename $f _R1.clean.fastq.gz).filtered_variants.vcf -filter-name "QD_filter" -filter "QD < 10.0"

###annotation filtered variants 
java -Xmx8g -jar /gpfs/data/user/rajesh/software/tools/snpEff/snpEff.jar corona_virus filter_vcf/$(basename $f _R1.clean.fastq.gz).filtered_variants.vcf -htmlStats annotated_vcf/$(basename $f _R1.clean.fastq.gz).html -csvStats annotated_vcf/$(basename $f _R1.clean.fastq.gz).csv > annotated_vcf/$(basename $f _R1.clean.fastq.gz).filtered_variants.anno.vcf

done
