

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
