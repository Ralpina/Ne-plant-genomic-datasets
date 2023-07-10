module load vcftools/0.1.16
# 27 genomic scaffolds:
vcftools --gzvcf ventoux.recode.vcf.gz --chr Fsylvatica_scaffold_4 --chr Fsylvatica_scaffold_6 --chr Fsylvatica_scaffold_11 --chr Fsylvatica_scaffold_7 --chr Fsylvatica_scaffold_5 --chr Fsylvatica_scaffold_8 --chr Fsylvatica_scaffold_13 --chr Fsylvatica_scaffold_10 --chr Fsylvatica_scaffold_9 --chr Fsylvatica_scaffold_18 --chr Fsylvatica_scaffold_14 --chr Fsylvatica_scaffold_12 --chr Fsylvatica_scaffold_15 --chr Fsylvatica_scaffold_16 --chr Fsylvatica_scaffold_17 --chr Fsylvatica_scaffold_19 --chr Fsylvatica_scaffold_21 --chr Fsylvatica_scaffold_24 --chr Fsylvatica_scaffold_20 --chr Fsylvatica_scaffold_23 --chr Fsylvatica_scaffold_25 --chr Fsylvatica_scaffold_27 --chr Fsylvatica_scaffold_22 --chr Fsylvatica_scaffold_26 --chr Fsylvatica_scaffold_28 --chr Fsylvatica_scaffold_30 --chr Fsylvatica_scaffold_31 --recode --out 27scaffolds
# output file: 27scaffolds.recode.vcf

# 12 genomic scaffolds:
vcftools --gzvcf ventoux.recode.vcf.gz --chr Fsylvatica_scaffold_1 --chr Fsylvatica_scaffold_2 --chr Fsylvatica_scaffold_3 --chr Fsylvatica_scaffold_4 --chr Fsylvatica_scaffold_5 --chr Fsylvatica_scaffold_6 --chr Fsylvatica_scaffold_7 --chr Fsylvatica_scaffold_8 --chr Fsylvatica_scaffold_9 --chr Fsylvatica_scaffold_10 --chr Fsylvatica_scaffold_11 --chr Fsylvatica_scaffold_12 --recode --out 12scaffolds
# output file: 12scaffolds.recode.vcf

