#####################################
### Data preprocess
#####################################

###Remove adapter and prune reads according to base quality（Trimmomatic）
java -jar /opt/tools/seq-analysis/Trimmomatic-0.36/trimmomatic-0.36.jar SE(PE) -threads 20 -phred33 *.fastq  ILLUMINACLIP:/opt/tools/seq-analysis/Trimmomatic-0.36/adapters/TruSeq2-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36(50) #单端一般36 双端一般50


###QC
/opt/tools/FastQC/fastqc -t 8 -o /zs32/home/shliu/FastQC *.fastq


#####################################
### Call genotypes
#####################################

### Mapping
# Map to reference genome
bwa mem -M -k 48 -t 7 -R "@RG\tID:Hiseq\tPL:Illumina\tSM:SOMAD_0004" /zs32/data-analysis/reflib/genomes/human_UCSC_hg19/chrAll.fasta /zs32/home/shliu/chr_all/SOMAD_0004.fq.gz | samtools view -b -S - | samtools sort -m 500000000 - align/SOMAD_0004.sorted

#Statistical mapping information
samtools flagstat *.bam>*.stat 

# Remove PCR duplication
/zs32/home/yjiang/bin/java -jar /opt/tools/seq-analysis/picard-tools-1.119/MarkDuplicates.jar INPUT=align/SOMAD_0004.sorted.bam OUTPUT=rmdup/SOMAD_0004.rmdup.sorted.bam METRICS_FILE=rmdup/SOMAD_0004.metrics REMOVE_DUPLICATES=true PROGRAM_RECORD_ID=null

samtools index rmdup/SOMAD_0004.rmdup.sorted.bam

# indelrealigner
/usr/bin/java -Xms14g -Xmx18g -jar /zs32/home/yjiang/softwares/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar  -T RealignerTargetCreator -R /zs32/data-analysis/reflib/genomes/human_UCSC_hg19/chrAll.fasta -o gatk/SOMAD_0004.bam.list -I rmdup/SOMAD_0004.rmdup.sorted.bam -nt 6

/usr/bin/java -Xms14g -Xmx18g -Djava.io.tmpdir=gatk/SOMAD_0004.tmpdir -jar /zs32/home/yjiang/softwares/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar  -I rmdup/SOMAD_0004.rmdup.sorted.bam -R /zs32/data-analysis/reflib/genomes/human_UCSC_hg19/chrAll.fasta -T IndelRealigner -targetIntervals gatk/SOMAD_0004.bam.list -o gatk/SOMAD_0004.realigned.bam --maxReadsForRealignment 30000 --maxReadsInMemory 1000000

samtools index gatk/SOMAD_0004.realigned.bam

### BQSR
/usr/bin/java -Xms14g -Xmx18g -jar /zs32/home/yjiang/softwares/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar  -l INFO -R /zs32/data-analysis/reflib/genomes/human_UCSC_hg19/chrAll.fasta -T BaseRecalibrator --knownSites:dbsnp,VCF /zs32/home/yjiang/database/dbsnp147.GRCH37/All_20160601.vcf.gz -I gatk/SOMAD_0004.realigned.bam -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov ContextCovariate -o gatk/SOMAD_0004.recal_data.csv --default_platform illumina -nct 6

/usr/bin/java -Xms14g -Xmx18g -jar /zs32/home/yjiang/softwares/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar  -l INFO -R /zs32/data-analysis/reflib/genomes/human_UCSC_hg19/chrAll.fasta -I gatk/SOMAD_0004.realigned.bam -T PrintReads -o gatk/SOMAD_0004.realigned.recal.bam -BQSR gatk/SOMAD_0004.recal_data.csv

samtools index gatk/SOMAD_0003.realigned.recal.bam

### Genotyping
# call genotypes (HaplotypeCaller; per chr a time)
/usr/bin/java -Xmx20g -jar /zs32/home/yjiang/softwares/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar -l INFO -R /zs32/data-analysis/reflib/genomes/human_UCSC_hg19/chrAll.fasta -T HaplotypeCaller -nct 6 --dbsnp /zs32/home/yjiang/database/dbsnp147.GRCH37/All_20160601.vcf.gz -o vcf/high.interval.raw.vcf -L chr1  -I SOMAD_0001.realigned.recal.bam -I SOMAD_0003.realigned.recal.bam -I SOMAD_0004.realigned.recal.bam

### VQSR
# snp recalibration
/usr/bin/java -Xmx20g -jar /zs32/home/yjiang/softwares/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar -T VariantRecalibrator -R /zs32/data-analysis/reflib/genomes/human_UCSC_hg19/chrAll.fasta -input vcf/high.interval.raw.vcf -resource:dbsnp,known=true,training=true,truth=true,prior=6.0 /zs32/home/yjiang/database/dbsnp147.GRCH37/All_20160601.vcf.gz -an QD  -an MQRankSum -an ReadPosRankSum -an FS -an MQ -mode SNP -recalFile vcf/snp.recal -tranchesFile vcf/snp.tranches -rscriptFile vcf/snp.plots.R -nt 6 --maxGaussians 6

/usr/bin/java -Xmx20g -jar /zs32/home/yjiang/softwares/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar -T ApplyRecalibration -R /zs32/data-analysis/reflib/genomes/human_UCSC_hg19/chrAll.fasta -input vcf/high.interval.raw.vcf --ts_filter_level 99.0 -tranchesFile vcf/snp.tranches -recalFile vcf/snp.recal -mode SNP -o vcf/snp.recal.vcf -nt 6

# indel recalibration
/usr/bin/java -Xmx20g -jar /zs32/home/yjiang/softwares/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar -T VariantRecalibrator -R /zs32/data-analysis/reflib/genomes/human_UCSC_hg19/chrAll.fasta -input vcf/snp.recal.vcf -resource:dbsnp,known=true,training=true,truth=true,prior=6.0 /zs32/home/yjiang/database/dbsnp147.GRCH37/All_20160601.vcf.gz -an QD  -an MQRankSum -an ReadPosRankSum -an FS -an MQ -mode INDEL -recalFile vcf/both.recal -tranchesFile vcf/both.tranches -rscriptFile vcf/both.plots.R -nt 6 --maxGaussians 6

/usr/bin/java -Xmx20g -jar /zs32/home/yjiang/softwares/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar -T ApplyRecalibration -R /zs32/data-analysis/reflib/genomes/human_UCSC_hg19/chrAll.fasta -input vcf/snp.recal.vcf -ts_filter_level 99.0 -tranchesFile vcf/both.tranches -recalFile vcf/both.recal -mode INDEL -o vcf/both.recal.vcf -nt 6



#####################################
### IMPUTATION
#####################################

## BEAGLE(none-ref-based or ref-based)

for i in {1..22} X; do time awk '$1~/^#/||$7=="PASS"' chr$i.both.recal.vcf  > chr$i.both.filter.vcf & done

cat $chr.both.filter.vcf |sed -e 's/^chr//' > $chr.both.filter.vcf

#none-ref-based imputation

/opt/tools/java_v8/jre1.8.0_121/bin/java -jar beagle.27Jan18.7e1.jar map=/zs32/home/yjiang/softwares/beagle/map/plink.chr1.GRCh37.map gl=chr1.vcf out=chr1.both.filter.gl

# BEAGLE ref-based imputation
java -jar /zs32/home/yjiang/softwares/beagle/beagle.27Jul16.86a.jar map=/zs32/home/yjiang/softwares/beagle/map/plink.chr1.GRCh37.map ref=/zs32/home/yjiang/softwares/beagle/reference/chr1.1kg.phase3.v5a.bref gt=ATAC.chr1.raw.gl .vcf.gz out=ATAC.chr1.raw.bref

#Removes snps which AR2<0.8
zcat chr1.both.filter.gl.vcf.gz | awk '{if($1~/^#/){print}else{flag=0;split($8,a,";");for(i in a){if(a[i]~/^AR2=/&&substr(a[i],5)>=0.8){flag=1}}if(flag==1){print}}}'|gzip >chr1.beagle.result.vcf.gz



# ## SHAPEIT2
# ls /zs32/data-analysis/liucy_group/jiangyi/WGS/gatk/*.bam|awk 'BEGIN{FS="[/.]"}{print $8,$0,"chr2"}' > combine136.chr2.bamlist
# awk '$0~/^#/||$5!~/,/' ../vcf/combine136.chr2.both.filter.vcf > combine136.chr2.biallelic.vcf
# /zs32/home/yjiang/softwares/extractPIRs.v1.r68.x86_64/extractPIRs --bam combine136.chr2.bamlist --vcf combine136.chr2.biallelic.vcf --out combine136.chr2.PIRsList
# shapeit -check --input-vcf combine136.chr2.biallelic.vcf -R /zs32/home/yjiang/softwares/softwares/impute_v2.3.2_x86_64_dynamic/Reference/1000GP_Phase3/1000GP_Phase3_chr2.hap.gz /zs32/home/yjiang/softwares/softwares/impute_v2.3.2_x86_64_dynamic/Reference/1000GP_Phase3/1000GP_Phase3_chr2.legend.gz /zs32/home/yjiang/softwares/softwares/impute_v2.3.2_x86_64_dynamic/Reference/1000GP_Phase3/1000GP_Phase3.sample --output-log combine136.chr2.AlignmentChecks
# shapeit -assemble --input-vcf combine136.chr2.biallelic.vcf --input-pir combine136.chr2.PIRsList -R /zs32/home/yjiang/softwares/softwares/impute_v2.3.2_x86_64_dynamic/Reference/1000GP_Phase3/1000GP_Phase3_chr2.hap.gz /zs32/home/yjiang/softwares/softwares/impute_v2.3.2_x86_64_dynamic/Reference/1000GP_Phase3/1000GP_Phase3_chr2.legend.gz /zs32/home/yjiang/softwares/softwares/impute_v2.3.2_x86_64_dynamic/Reference/1000GP_Phase3/1000GP_Phase3.sample -O combine136.chr2.Haplotype --exclude-snp combine136.chr2.AlignmentChecks.snp.strand.exclude --no-mcmc
# python3.5 statistics.phase.py -i combine136.chr2.AlignmentChecks.log -s combine136 -c chr2 -d null --head > statistics.phase.xls
# 
# ## IMPUTE2
# for((i=5000000;i<=254250620;i=$i+5000000)); do let ii=$i-4999999; nohup impute2 -filt_rules_l 'ALL<more0.3' 'TYPE==highCOV' -use_prephased_g -h /zs32/home/yjiang/softwares/softwares/impute_v2.3.2_x86_64_dynamic/Reference/1000GP_Phase3/1000GP_Phase3_chr2.hap.gz -l /zs32/home/yjiang/softwares/softwares/impute_v2.3.2_x86_64_dynamic/Reference/1000GP_Phase3/1000GP_Phase3_chr2.legend.gz -m /zs32/home/yjiang/softwares/softwares/impute_v2.3.2_x86_64_dynamic/Reference/1000GP_Phase3/genetic_map_chr2_combined_b37.txt -known_haps_g combine136.chr2.Haplotype.haps -int $ii $i -o combine136.chr2.$ii-$i.impute2 -phase & sleep 30; done

## annotation
/opt/tools/seq-analysis/annovar201602/table_annovar.pl chrall.beagle.result.vcf /zs32/data-analysis/reflib/annovar_humandb -outfile /zs32/home/shliu/annovar-highdeep.vcf -buildver hg19 -protocol refGene,cytoBand,avsnp147,exac03,esp6500siv2_all,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_amr,1000g2015aug_eas,1000g2015aug_eur,1000g2015aug_sas,cg69,hrcr1,kaviar_20150923,clinvar_20161128,genomicSuperDups,gerp++elem,gerp++gt2,caddgt20,cadd13gt20,revel,mcap,ljb26_all -operation g,r,f,f,f,f,f,f,f,f,f,f,f,f,f,r,r,f,f,f,f,f,f --remove --otherinfo --verbose --vcfinput -nastring.

#####################################
### eQTL
#####################################

# FastQTL
# input: genotypes (VCF.gz) & gene expression data(bed.gz)
# Linux

# 校正(lm)：
# 已知因素：性别、年龄.....
# 未知因素：SVA (top ## PCs),PEER

mod = model.matrix( ~ sample$BrainWeight + sample$RIN + sample$YearAutopsy + sample$Ethnicity + sample$SequencingPlatform + sample$TissueState + sample$Sex + sample$Hemisphere)
n.sv = num.sv(Lmbeta, mod)   # n.sv=4
svobj = sva(Lmbeta, mod, n.sv=n.sv)
modSv = cbind(mod, svobj$sv)


###FASTQTL
#Keep snps which maf > 0.05 or hwe>0.000001
for i in {1..22} X; do time vcftools --gzvcf ../BEAGLE/combine136.chr$i.both.filter.AFgt0.05.gl.vcf.gz --plink --out combine136.chr$i.both.filter.AFgt0.05.gl & sleep 1; done
for i in {1..22} X; do time plink --file combine136.chr$i.both.filter.AFgt0.05.gl --out combine136.chr$i.both.filter.AFgt0.05.gl --hardy --noweb & sleep 1; done
for i in {1..22} X; do time plink --file combine136.chr$i.both.filter.AFgt0.05.gl --out combine136.chr$i.both.filter.AFgt0.05.gl --frq --noweb & sleep 1; done

for i in {1..22} X
do
cat chr$i.frq | awk '$5>0.05{print $2}' >chr$i.frqID.txt
cat chr$i.hwe | awk '$3=="ALL"&&$9>0.000001{print $2}' >chr$i.hweID.txt
cat chr$i.hweID.txt chr$i.frqID.txt >chr$i.all.txt
cat chr$i.all.txt | awk 'a[$0]++' >chr$i.final.txt
cat ../BEAGLE/combine151.chr$i.rs.vcf | awk 'NR==FNR{a[$1]=0}NR>FNR&&FNR<=11{print $0}NR>FNR&&FNR>=12&&(($3 in a))' chr$i.final.txt - >../../eQTL/combine151.chr$i.vcf
done

#bgzip and tabix
cd eQTL
for i in {1..22} X 
do
bgzip combine151.chr$i.vcf && tabix -p vcf combine151.chr$i.vcf.gz
done

#FastQTL(normial and permutation)
for i in {1..22} 
do
/zs32/home/yjiang/softwares/FastQTL/bin/fastQTL.static --vcf combine151.chr$i.vcf.gz --bed phenotype.bed.gz --region "$i":10000-250000000 --out "$i"_nominals.default.txt.gz --cov cov.txt
/zs32/home/yjiang/softwares/FastQTL/bin/fastQTL.static --vcf combine151.chr$i.vcf.gz --bed phenotype.bed.gz --region "$i":10000-250000000 --permute 1000 --out "$i"_permutations.default.txt.gz --cov cov.txt
done
