# Ne estimation in GONE

The scripts presented here describe the analyses carried out on four plant genomic datasets to estimate effective population size (Ne) using the software GONE [Santiago et al. 2020](https://doi.org/10.1093/molbev/msaa169)
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
plink v1.9:

## Directory conventions
- data
- results
- indlist
- snpslist


## *Symphonia globulifera*

## *Mercurialis annua*

## *Fagus sylvatica*

## *Prunus armeniaca*
#### Dataset: the full dataset from [Groppi et al. 2021](https://www.nature.com/articles/s41467-021-24283-6)
data/apricot_collection_2019_marouch_v3.1.vcf.gz

We first generate subsets of individuals, based on the population structure (Q-values) found in Groppi et al. (2021) using fastStructure.
We are interested in two main clusters (red or Southern gene pool, and yellow or Northern gene pool), with Q-values per individual ranging from >70% to >99%.
The lists of individuals will be in the directory "indlist", where the number corresponds to Q-value for each genetic cluster.

```sh
# indlist/yellow99
# indlist/yellow95
# indlist/yellow90
# indlist/yellow80
# indlist/yellow70
# indlist/red99
# indlist/red95
# indlist/red90
# indlist/red80
# indlist/red70
# indlist/all (=Red70+Yellow70)
```

We first remove all the indels, as we only need SNPs; simultaneously, we also extract only the individuals we need, to create smaller files!

```sh
cd data
vcftools --gzvcf apricot_collection_2019_marouch_v3.1.vcf.gz --remove-indels --keep ./indlist/all --recode --recode-INFO-all --stdout | gzip -c > armeniaca.SNPs.vcf.gz 
```

Our new dataset is data/armeniaca.SNPs.vcf.gz

To count the number of unique chromosomes or contigs:
```sh
bcftools view -H armeniaca.SNPs.vcf.gz | cut -f 1 | uniq | awk '{print $0"\t"$0}' | wc -l
```
8 chromosomes!

let's see how many SNPs per chromosome:
```sh
bcftools view -H armeniaca.SNPs.vcf.gz | cut -f 1 | sort | uniq -c > sites_per_contigs
cat sites_per_contigs
# 6203547 chr1
# 4110478 chr2
# 3388721 chr3
# 3935530 chr4
# 2329589 chr5
# 3886518 chr6
# 3071933 chr7
# 2941746 chr8
```

We will need to subsample the number of sites eventually, as GONE accepts a maximum of 10 million SNPs (or 1 million per chromosome)


#### 1. Influence of number of SNPs on Ne estimation. 

We select only the individuals with Q-value > 99% (77 individuals) from the Northern gene pool (yellow gene pool), to avoid the influence of population structure on the results obtained.

```sh
module load vcftools/0.1.16
vcftools --gzvcf armeniaca.SNPs.vcf.gz --keep ../indlist/yellow99 --plink --out yellow99.SNPs
```
From this, we obtain a map and a ped file.
We use the map file to extract the list of chromosomes and their variant sites.
```sh
cut -f 2 yellow99.SNPs.map > ./snpslist/yellow99.snps           
wc -l yellow99.snps
# 23974765 yellow99.snps
```
We want to subsample SNPs to understand how Ne estimates change depending on the number of SNPs. 
We also want to create confidence intervals, so we replicate each sampling procedure several times (50), to obtain different estimates across subsampled datasets.
We can use the ```shuf``` function to sample without replacement:
```sh
cd ../
for i in {1..50}; do
 shuf -n 10000000 ./snpslist/yellow99.snps | sort > ./snpslist/subset10M$i.snps       
done
```
where ```-n``` tells the shuffle command to sample the specified maximum number of lines.
```sort``` will reorder the list 50 SNPs subsets were created this way.

We repeat the above also by subsampling 50 times 10000000, 7000000, 3500000, 500000, 300000, 150000, 80000 and 40000 SNPs (every time, we change the file names accordingly).
If needed, we can check how many SNPs per chromosome have been subsampled, for example in the first subset:
```sh
grep -c "chr1" ./snpslist/subset10M1.snps 
```

Now we need to create the input files for GONE, and we can use the initial ped and map files generated, specified in the command below by "yellow99.SNPs" without extensions:
```sh
cd ../
module load plink/1.90
for i in {1..50}; do
  plink --file ./data/yellow99.SNPs --recode --extract ./snpslist/subset10M$i.snps --out ./data/yellow99.subset10M$i
done
```
This produces subsampled ped and map files, that we need to convert to the format required by GONE.
First, we need to convert the delimiter from space to tab:
```sh
for i in {1..50}; do
 tr ' ' \\t < ./data/yellow99.subset10M$i.ped > ./data/temp$i
done
```
then we rename the temp files:
```sh
for i in {1..50}; do
 mv ./data/temp$i ./data/yellow99.subset10M$i.ped
done
```
We then prepare a six-columns, tab-delimited file called ```indlist/Y99```, that looks like this:
```sh
1	IND1	0	0	1	-9
1	IND2	0	0	1	-9
1	IND3	0	0	1	-9
1	IND4	0	0	1	-9
1	IND5	0	0	1	-9
...
1	IND77	0	0	1	-9
```
in order to replace the first six columns of the ped files obtained above (GONE likes individual names in the format IND followed by a number).
We cut away the first seven columns of our ped files and replace them with those in the six-columns file:
```sh
for i in {1..50}; do
 cut -f7- ./data/yellow99.subset10M$i.ped > ./data/yellow99.subset10M${i}Geno
 paste ./indlist/Y99 ./data/yellow99.subset10M${i}Geno > ./data/yellow99.subset10M$i.ped
done
```
then we can remove the temporary files:
```sh
for i in {1..50}; do
 rm ./data/yellow99.subset10M${i}Geno
done
```
We then transfer the obtained ped and map files to the directory where we will run GONE ```results/yellow99```
```sh
mv ./data/yellow99.subset10M* ./results/yellow99
```
And finally we run GONE:
```sh
cd ./results/yellow99
for i in {1..50}; do
 bash script_GONE.sh yellow99.subset10M$i
done
```
Once the analyses in GONE are complete, we want to calculate the average number of SNPs per chromosome used by GONE:
```sh
for i in {1..50}; do
 grep "NSNP_calculations=" ./results/yellow99/OUTPUT_yellow99.subset10M$i | cut -d "=" -f 2 
done | awk 'BEGIN{s=0;}{s+=$1;}END{print s/NR;}'  
```
For the dataset with 10 million SNPs, for example, the average number of SNPs per chromosome is ```20777.5```.
We will report these in the plot (see below).
Then, we need to extract the Ne estimates and use them to calculate median values and confidence intervals (for Ne estimates relating to the last generation). 
For example:
```sh
for i in {1..50}; do
 awk 'NR==3' ./results/yellow99/Output_Ne_yellow99.subset10M$i | cut -f 2
done >> ./results/yellow99/Ne10M.txt
```

We calculate median values and 95% confidence intervals using R
```sh
module load R/4.1.0
R
Ne10M <- read.delim("./results/yellow99/Ne10M.txt")
median(Ne10M[[1]])
quantile(Ne10M[[1]], 0.025, type = 1)
quantile(Ne10M[[1]], 0.975, type = 1)
```

We can now prepare the input file for plotting the results.
Plotting the results in R:
``` 
SNPsubsets <- read.delim("SNPsubsets.txt")
r <-ggplot(SNPsubsets, aes(x=SNPs, y=Ne)) + geom_point(colour="#00C08D", size = 10) + geom_line(colour="#00C08D", linewidth = 6) + ylim(0,10000)

rp <- r+geom_ribbon(aes(ymin=LowCI, ymax=HighCI), linetype=2, alpha=0.3, bg = "#00C08D") +
  theme_classic() + theme(axis.text.x = element_text(face="bold", size=24), axis.text.y = element_text(face="bold", size=24)) +
  theme(axis.title.x = element_text(size=24, face="bold"), axis.title.y = element_text(size=24, face="bold")) +
  scale_y_continuous(labels=comma) + scale_x_continuous(labels=comma) + labs(x = "average number of SNPs per chromosome (used by GONE)", y = "Ne (Geometric mean)") 
```


  
  
















