### Mapping
# 比对到参考基因组
bwa mem -M -k 48 -t 7 -R "@RG\tID:Hiseq\tPL:Illumina\tSM:SOMAD_0001" /zs32/data-analysis/reflib/genomes/human_UCSC_hg19/chrAll.fasta /zs32/home/shliu/sampling/SOMAD_0001.fq.gz | samtools view -b -S - | samtools sort -m 500000000 - align/SOMAD_0001.sorted


# 去除PCR duplication
/zs32/home/yjiang/bin/java -jar /opt/tools/seq-analysis/picard-tools-1.119/MarkDuplicates.jar INPUT=align/SOMAD_0001.sorted.bam OUTPUT=rmdup2/SOMAD_0001.rmdup.sorted.bam METRICS_FILE=rmdup2/SOMAD_0001.metrics REMOVE_DUPLICATES=true PROGRAM_RECORD_ID=null

samtools index rmdup2/SOMAD_0001.rmdup.sorted.bam

# indelrealigner
/usr/bin/java -Xms14g -Xmx18g -jar /zs32/home/yjiang/softwares/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar  -T RealignerTargetCreator -R /zs32/data-analysis/reflib/genomes/human_UCSC_hg19/chrAll.fasta -o gatk2/SOMAD_0001.bam.list -I rmdup2/SOMAD_0001.rmdup.sorted.bam -nt 6

/usr/bin/java -Xms14g -Xmx18g -Djava.io.tmpdir=gatk2/SOMAD_0001.tmpdir -jar /zs32/home/yjiang/softwares/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar  -I rmdup2/SOMAD_0001.rmdup.sorted.bam -R /zs32/data-analysis/reflib/genomes/human_UCSC_hg19/chrAll.fasta -T IndelRealigner -targetIntervals gatk2/SOMAD_0001.bam.list -o gatk2/SOMAD_0001.realigned.bam --maxReadsForRealignment 30000 --maxReadsInMemory 1000000

samtools index gatk2/SOMAD_0001.realigned.bam

### BQSR
/usr/bin/java -Xms14g -Xmx18g -jar /zs32/home/yjiang/softwares/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar  -l INFO -R /zs32/data-analysis/reflib/genomes/human_UCSC_hg19/chrAll.fasta -T BaseRecalibrator --knownSites:dbsnp,VCF /zs32/home/yjiang/database/dbsnp147.GRCH37/All_20160601.vcf.gz -I gatk2/SOMAD_0001.realigned.bam -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov ContextCovariate -o gatk2/SOMAD_0001.recal_data.csv --default_platform illumina -nct 6

/usr/bin/java -Xms14g -Xmx18g -jar /zs32/home/yjiang/softwares/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar  -l INFO -R /zs32/data-analysis/reflib/genomes/human_UCSC_hg19/chrAll.fasta -I gatk2/SOMAD_0001.realigned.bam -T PrintReads -o gatk2/SOMAD_0001.realigned.recal.bam -BQSR gatk2/SOMAD_0001.recal_data.csv

samtools index gatk2/SOMAD_0001.realigned.recal.bam

# 查看BAM文件
samtools view -h gatk/8521602001262.realigned.recal.bam|less -S

### Genotyping
# call genotypes (HaplotypeCaller)
/usr/bin/java -Xmx20g -jar /zs32/home/yjiang/softwares/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar -l INFO -R /zs32/data-analysis/reflib/genomes/human_UCSC_hg19/chrAll.fasta -T HaplotypeCaller -nct 6 --dbsnp /zs32/home/yjiang/database/dbsnp147.GRCH37/All_20160601.vcf.gz -o vcf/$sampleName.interval.raw.vcf -L /zs32/data-analysis/liucy_group/shareData/TGEN/TGEN/bams/CDS.chr1.bed  -I SOMAD_0001.realigned.recal.bam -I SOMAD_0002.realigned.recal.bam -I SOMAD_0003.realigned.recal.bam

### VQSR(变异位点质量值重新校正:VQSLOD一个位点是真实的概率比上这个位点可能是假阳性的概率的log odds ratio（对数差异比），因此，可以定性的认为，这个值越大就越好)
# snp recalibration
/usr/bin/java -Xmx20g -jar /zs32/home/yjiang/softwares/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar -T VariantRecalibrator -R /zs32/data-analysis/reflib/genomes/human_UCSC_hg19/chrAll.fasta -input vcf2/interval.raw.vcf -resource:dbsnp,known=true,training=true,truth=true,prior=6.0 /zs32/home/yjiang/database/dbsnp147.GRCH37/All_20160601.vcf.gz -an QD  -an MQRankSum -an ReadPosRankSum -an FS -an MQ -mode SNP -recalFile vcf2/snp.recal -tranchesFile vcf2/snp.tranches -rscriptFile vcf2/snp.plots.R -nt 6 --maxGaussians 6

/usr/bin/java -Xmx20g -jar /zs32/home/yjiang/softwares/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar -T ApplyRecalibration -R /zs32/data-analysis/reflib/genomes/human_UCSC_hg19/chrAll.fasta -input vcf2/interval.raw.vcf --ts_filter_level 99.0 -tranchesFile vcf2/snp.tranches -recalFile vcf2/snp.recal -mode SNP -o vcf2/snp.recal.vcf -nt 6

# indel recalibration
/usr/bin/java -Xmx20g -jar /zs32/home/yjiang/softwares/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar -T VariantRecalibrator -R /zs32/data-analysis/reflib/genomes/human_UCSC_hg19/chrAll.fasta -input vcf2/snp.recal.vcf -resource:dbsnp,known=true,training=true,truth=true,prior=6.0 /zs32/home/yjiang/database/dbsnp147.GRCH37/All_20160601.vcf.gz -an QD  -an MQRankSum -an ReadPosRankSum -an FS -an MQ -mode INDEL -recalFile vcf2/both.recal -tranchesFile vcf2/both.tranches -rscriptFile vcf2/both.plots.R -nt 6 --maxGaussians 6

/usr/bin/java -Xmx20g -jar /zs32/home/yjiang/softwares/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar -T ApplyRecalibration -R /zs32/data-analysis/reflib/genomes/human_UCSC_hg19/chrAll.fasta -input vcf2/snp.recal.vcf -ts_filter_level 99.0 -tranchesFile vcf2/both.tranches -recalFile vcf2/both.recal -mode INDEL -o vcf2/both.recal.vcf -nt 6


#####################################
### IMPUTATION
#####################################

## BEAGLE
/zs32/home/yjiang/bin/java -jar /zs32/home/yjiang/softwares/beagle/beagle.27Jul16.86a.jar map=/zs32/home/yjiang/softwares/beagle/map/plink.chr1.GRCh37.map gl=interval.raw.vcf out=combine136.both.filter.gl(vcf文件的chr要删除，看map文件就知道)

# ## SHAPEIT2
# ls /zs32/data-analysis/liucy_group/jiangyi/WGS/gatk/*.bam|awk 'BEGIN{FS="[/.]"}{print $8,$0,"chr1"}' > combine136.chr1.bamlist
# awk '$0~/^#/||$5!~/,/' ../vcf/combine136.chr1.both.filter.vcf > combine136.chr1.biallelic.vcf
# ~/softwares/extractPIRs.v1.r68.x86_64/extractPIRs --bam combine136.chr1.bamlist --vcf combine136.chr1.biallelic.vcf --out combine136.chr1.PIRsList
# shapeit -check --input-vcf combine136.chr1.biallelic.vcf -R ~/softwares/impute_v2.3.2_x86_64_dynamic/Reference/1000GP_Phase3/1000GP_Phase3_chr1.hap.gz ~/softwares/impute_v2.3.2_x86_64_dynamic/Reference/1000GP_Phase3/1000GP_Phase3_chr1.legend.gz ~/softwares/impute_v2.3.2_x86_64_dynamic/Reference/1000GP_Phase3/1000GP_Phase3.sample --output-log combine136.chr1.AlignmentChecks
# shapeit -assemble --input-vcf combine136.chr1.biallelic.vcf --input-pir combine136.chr1.PIRsList -R ~/softwares/impute_v2.3.2_x86_64_dynamic/Reference/1000GP_Phase3/1000GP_Phase3_chr1.hap.gz ~/softwares/impute_v2.3.2_x86_64_dynamic/Reference/1000GP_Phase3/1000GP_Phase3_chr1.legend.gz ~/softwares/impute_v2.3.2_x86_64_dynamic/Reference/1000GP_Phase3/1000GP_Phase3.sample -O combine136.chr1.Haplotype --exclude-snp combine136.chr1.AlignmentChecks.snp.strand.exclude --no-mcmc
# python3.5 statistics.phase.py -i combine136.chr1.AlignmentChecks.log -s combine136 -c chr1 -d null --head > statistics.phase.xls
# 
# ## IMPUTE2
# for((i=5000000;i<=254250620;i=$i+5000000)); do let ii=$i-4999999; nohup impute2 -filt_rules_l 'ALL<0.001' 'TYPE==LOWCOV' -use_prephased_g -h ~/softwares/impute_v2.3.2_x86_64_dynamic/Reference/1000GP_Phase3/1000GP_Phase3_chr1.hap.gz -l ~/softwares/impute_v2.3.2_x86_64_dynamic/Reference/1000GP_Phase3/1000GP_Phase3_chr1.legend.gz -m ~/softwares/impute_v2.3.2_x86_64_dynamic/Reference/1000GP_Phase3/genetic_map_chr1_combined_b37.txt -known_haps_g combine136.chr1.Haplotype.haps -int $ii $i -o combine136.chr1.$ii-$i.impute2 -phase & sleep 30; done

## annotation
/opt/tools/seq-analysis/annovar201602/table_annovar.pl /zs32/home/shliu/gatk2 /zs32/data-analysis/reflib/annovar_humandb -outfile combine136.$chr.both.filter.gl -buildver hg19 -protocol refGene,cytoBand,avsnp147,1000g2015aug_all,gnomad_genome -operation g,r,f,f,f --remove --otherinfo --verbose --vcfinput -nastring .
