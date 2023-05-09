# Ne estimation in GONE

The scripts presented here describe the analyses carried out on four plant genomic datasets to estimate effective population size (Ne) using the software GONE.
The aims of these analyses are described in Gargiulo et al. 2023. Please cite this article when using the scripts described here.

* Contact: r.gargiulo@kew.org | [@RobertaGargiu10](https://twitter.com/RobertaGargiu10)

## Datasets used:
- *Symphonia globulifera* in [Schmitt et al. 2021](https://onlinelibrary.wiley.com/doi/10.1111/mec.16116)
- *Mercurialis annua* in [González-Martínez et al. 2017](https://www.sciencedirect.com/science/article/pii/S0960982217308655)
- *Fagus sylvatica* (Ivan Scotti, pers. comm.)
- *Prunus armeniaca* in [Groppi et al. 2021](https://www.nature.com/articles/s41467-021-24283-6)

## Software required
vcftools: https://vcftools.github.io/index.html
bcftools
GONE: [Santiago et al. 2020]

## Directory conventions
- data
- results
- indlist
- snpslist


## *Symphonia globulifera*

## *Mercurialis annua*

## *Fagus sylvatica*

## *Prunus armeniaca*
#### Dataset: full dataset from Groppi et al., 2021:
data/apricot_collection_2019_marouch_v3.1.vcf.gz

We first generate subsets of individuals, based on the population structure (Q-values) found in Groppi et al. (2021) using fastStructure.
We are interested in two main clusters (red or Southern gene pool, and yellow or Northern gene pool), with Q-values per individual ranging from >70% to >99%.
The lists of individuals will be in the directory "indlist", where the number corresponds to Q-value for each genetic cluster.

- indlist/yellow99
- indlist/yellow95
- indlist/yellow90
- indlist/yellow80
- indlist/yellow70
- indlist/red99
- indlist/red95
- indlist/red90
- indlist/red80
- indlist/red70
- indlist/all (=Red70+Yellow70)

We first remove all the indels, as we only need SNPs; simultaneously, we also extract only the individuals we need, to create smaller files!

```
cd data
vcftools --gzvcf apricot_collection_2019_marouch_v3.1.vcf.gz --remove-indels --keep ./indlist/all --recode --recode-INFO-all --stdout | gzip -c > armeniaca.SNPs.vcf.gz 
```

Our new dataset is data/armeniaca.SNPs.vcf.gz

To count the number of unique chromosomes or contigs:
```
bcftools view -H armeniaca.SNPs.vcf.gz | cut -f 1 | uniq | awk '{print $0"\t"$0}' | wc -l
8
```
8 chromosomes!

let's see how many SNPs per chromosome:
```
bcftools view -H armeniaca.SNPs.vcf.gz | cut -f 1 | sort | uniq -c > sites_per_contigs
cat sites_per_contigs
6203547 chr1
4110478 chr2
3388721 chr3
3935530 chr4
2329589 chr5
3886518 chr6
3071933 chr7
2941746 chr8
```
They are too many for GONE, which accepts 10 million SNPs (or 1 million per chromosome) 


















