# *N*<sub>e</sub> estimation in GONE

The analyses presented here describe effective population size (*N*<sub>e</sub>) estimation carried out for four plant genomic datasets,  using the software GONE [Santiago et al. 2020](https://doi.org/10.1093/molbev/msaa169).
The aims of these analyses are described in Gargiulo et al. 2023. Please cite the article and this repository when using the code described here. Most scripts are included in this README; a few additional scripts are available in the scripts directory of this repository.

* Contact: r.gargiulo@kew.org | [@RobertaGargiu10](https://twitter.com/RobertaGargiu10)

## Table of Contents
1. [Datasets used](https://github.com/Ralpina/Ne-plant-genomic-datasets#datasets-used)  
2. [Software programmes and tools required](https://github.com/Ralpina/Ne-plant-genomic-datasets#software-programmes-and-tools-required)  
3. [Directory conventions](https://github.com/Ralpina/Ne-plant-genomic-datasets#directory-conventions)  
4. [*Prunus armeniaca*](https://github.com/Ralpina/Ne-plant-genomic-datasets#prunus-armeniaca)
   - [Preparing the dataset](https://github.com/Ralpina/Ne-plant-genomic-datasets#preparing-the-dataset)
   - [1. Influence of number of SNPs on *N*<sub>e</sub> estimation](https://github.com/Ralpina/Ne-plant-genomic-datasets#1-influence-of-number-of-snps-on-ne-estimation)
   - [2. Influence of sample size on *N*<sub>e</sub> estimation](https://github.com/Ralpina/Ne-plant-genomic-datasets#2-influence-of-sample-size-resampling-only-individuals-with-the-smallest-proportion-of-admixture-on-ne-estimation)
   - [3. Influence of population structure on *N*<sub>e</sub> estimation](https://github.com/Ralpina/Ne-plant-genomic-datasets#3-influence-of-population-structure-on-ne-estimation)
   - [4. Effect of using genomic scaffolds instead of chromosomes on *N*<sub>e</sub> estimation](https://github.com/Ralpina/Ne-plant-genomic-datasets#4-effect-of-using-genomic-scaffolds-instead-of-chromosomes-on-ne-estimation)
   - [5. Influence of missing data on *N*<sub>e</sub> estimation](https://github.com/Ralpina/Ne-plant-genomic-datasets#5-influence-of-missing-data-on-ne-estimation)
5. [*Symphonia globulifera*](https://github.com/Ralpina/Ne-plant-genomic-datasets#symphonia-globulifera)
   - [Preparing the datasets and running GONE](https://github.com/Ralpina/Ne-plant-genomic-datasets#preparing-the-datasets-and-running-gone)
6. [*Mercurialis annua*](https://github.com/Ralpina/Ne-plant-genomic-datasets#mercurialis-annua)
   - [Preparing the datasets and running GONE](https://github.com/Ralpina/Ne-plant-genomic-datasets#preparing-the-datasets-and-running-gone-1)
7. [*Fagus sylvatica*](https://github.com/Ralpina/Ne-plant-genomic-datasets#fagus-sylvatica)
   - [Preparing the dataset](https://github.com/Ralpina/Ne-plant-genomic-datasets#preparing-the-dataset-1)
   - [Influence of missing data on *N*<sub>e</sub> estimation](https://github.com/Ralpina/Ne-plant-genomic-datasets#influence-of-missing-data-on-ne-estimation)
8. [References](https://github.com/Ralpina/Ne-plant-genomic-datasets#references)
  
## Datasets used:
- *Prunus armeniaca* in [Groppi et al. 2021](https://www.nature.com/articles/s41467-021-24283-6)
- *Symphonia globulifera* in [Schmitt et al. 2021](https://onlinelibrary.wiley.com/doi/10.1111/mec.16116)
- *Mercurialis annua* in [González-Martínez et al. 2017](https://www.sciencedirect.com/science/article/pii/S0960982217308655)
- *Fagus sylvatica* in [Lesur-Kupin & Scotti 2023](https://doi.org/10.57745/FJRYI1)

## Software programmes and tools required
[vcftools](https://vcftools.github.io/index.html)  
[bcftools](http://samtools.github.io/bcftools/)  
[GONE](https://github.com/esrud/GONE)  
[plink v1.9](https://www.cog-genomics.org/plink2/)  
[R](https://www.r-project.org/)  

## Directory conventions
All scripts include relative paths assuming that directories are organised as follows:
- data
- results
- indlist
- snpslist  

## *Prunus armeniaca*
#### Preparing the dataset
The full dataset from [Groppi et al. 2021](https://www.nature.com/articles/s41467-021-24283-6) is:  
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


#### 1. Influence of number of SNPs on *N*<sub>e</sub> estimation. 

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
We want to subsample SNPs to understand how *N*<sub>e</sub> estimates change depending on the number of SNPs. 
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
Then, we need to extract the *N*<sub>e</sub> estimates and use them to calculate median values and confidence intervals (for *N*<sub>e</sub> estimates relating to the last generation). 
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
setwd("./results/yellow99")
Ne10M <- read.delim("Ne10M.txt")
median(Ne10M[[1]])
quantile(Ne10M[[1]], 0.025, type = 1)
quantile(Ne10M[[1]], 0.975, type = 1)
```

We can now prepare the input file for plotting the results. The file ```SNPsubsets.txt``` looks like this:  
```
Ne	LowCI	HighCI	SNPs	totalSNPs
2828.0	2374.4	3128.9	20777.5	10000000
# ...etc.
```
Plotting the results in R:
``` 
library (ggplot2)
SNPsubsets <- read.delim("SNPsubsets.txt")
r <-ggplot(SNPsubsets, aes(x=SNPs, y=Ne)) + geom_point(colour="#00C08D", size = 10) + geom_line(colour="#00C08D", linewidth = 6) 

rp <- r+geom_ribbon(aes(ymin=LowCI, ymax=HighCI), linetype=2, alpha=0.3, bg = "#00C08D") +
  theme_classic() + theme(axis.text.x = element_text(face="bold", size=24, colour = "black"), axis.text.y = element_text(face="bold", size=24, colour = "black")) +
  theme(axis.title.x = element_text(size=24, face="bold", colour = "black"), axis.title.y = element_text(size=24, face="bold", colour = "black")) +
  scale_y_continuous(labels=comma, limits=c(0, 7000)) + scale_x_continuous(labels=comma) + labs(x = "Average number of SNPs per chromosome (used by GONE)", y = "Ne") 
 ```
#### 2. Influence of sample size (resampling only individuals with the smallest proportion of admixture) on *N*<sub>e</sub> estimation

We only select the individuals with Q-values > 99% (77 individuals) from the Northern gene pool (yellow gene pool), to avoid the influence of population structure on the results obtained.

We generate 50 permuted subsets of 15, 30, 45, 60, 75 individuals 
The ```shuf``` function will permute individuals and sample them without replacement.

```sh
for i in {1..50}; do
 shuf -n 15 ./indlist/yellow99 > ./indlist/indset15.$i
done

for i in {1..50}; do
 shuf -n 30 ./indlist/yellow99 > ./indlist/indset30.$i
done

for i in {1..50}; do
 shuf -n 45 ./indlist/yellow99 > ./indlist/indset45.$i
done

for i in {1..50}; do
 shuf -n 60 ./indlist/yellow99 > ./indlist/indset60.$i
done

for i in {1..50}; do
 shuf -n 75 ./indlist/yellow99 > ./indlist/indset75.$i
done
```

We also need to subsample SNPs because they are too many for GONE.
To get a list of SNPs to sample randomly, we can do the following (similarly to what we've done above):
```sh
module load vcftools/0.1.16
vcftools --gzvcf ./data/armeniaca.SNPs.vcf.gz --keep ./indlist/yellow99 --plink --out ./data/yellow99.SNPs
```
(it doesn't actually matter which individuals we include, because we only need the list of SNPs and that won't change among individual subsets).
From the map file generated above from vcftools, we extract the list of SNPs from which we sample 3.5 million SNPs:
```sh
cut -f 2 ./data/yellow99.SNPs.map > ./snpslist/yellow99.SNPs          
shuf -n 3500000 ./snpslist/yellow99.SNPs | sort  > ./snpslist/yellow3.5subset.snps       
```
Now we can generate the set of ped and map files for GONE, with resampled individuals:
```sh
module load vcftools/0.1.16
module load plink/1.90

for i in {1..50}; do 
 vcftools --gzvcf ./data/armeniaca.SNPs.vcf.gz --plink --keep ./indlist/indset15.$i --chrom-map ./data/chrom-map --out ./data/yellow99.15inds.$i
done

echo "subsampling individuals done"

for i in {1..50}; do 
 plink --file ./data/yellow99.15inds.$i --recode --extract ./snpslist/yellow3.5subset.snps --out ./data/15inds.$i
done

echo "subsampling SNPs done"

for i in {1..50}; do 
 rm ./data/yellow99.15inds.$i.ped 
done

for i in {1..50}; do 
 rm ./data/yellow99.15inds.$i.map 
done

echo "raw ped and map files created"

# to change spaces to tab delimiters:
for i in {1..50}; do
	tr ' ' \\t < ./data/15inds.$i.ped > ./data/test$i
done

for i in {1..50}; do
	mv ./data/test$i ./data/15inds.$i.ped
done

# to replace the first six columns:
for i in {1..50}; do
	cut -f7- ./data/15inds.$i.ped > ./data/15inds.${i}Geno
done

for i in {1..50}; do
	paste ./indlist/15 ./data/15inds.${i}Geno > ./data/15inds.$i.ped
done

# where ./indlist/15 includes the first six columns in the format required by GONE, for 15 individuals, as in: 
# 1	IND1	0	0	1	-9

echo "ped files for GONE obtained"

# We then transfer the obtained ped and map files to the directory where we will run GONE ```results/yellow99```

mv ./data/15inds.* ./results/yellow99
```

And finally we run GONE:
```sh
cd ./results/yellow99
for i in {1..50}; do
 bash script_GONE.sh 15inds.$i
done
```

Then repeat the above for each individual subsampling (30, 45, 60 and 75 individuals).

Preparing results to plot:
```sh
for i in {1..50}; do
 awk 'NR==3' ./results/yellow99/Output_Ne_15inds.$i | cut -f 2
done  >> ./results/yellow99/Ne15inds.txt
awk 'BEGIN{ FS = OFS = "\t" } { print $0, "15" }' ./results/yellow99/Ne15inds.txt > tmp && mv tmp ./results/yellow99/Ne15inds.txt

for i in {1..50}; do
 awk 'NR==3' ./results/yellow99/Output_Ne_30inds.$i | cut -f 2
done >> ./results/yellow99/Ne30inds.txt
awk 'BEGIN{ FS = OFS = "\t" } { print $0, "30" }' ./results/yellow99/Ne30inds.txt > tmp && mv tmp ./results/yellow99/Ne30inds.txt

for i in {1..50}; do
 awk 'NR==3' ./results/yellow99/Output_Ne_45inds.$i | cut -f 2
done  >> ./results/yellow99/Ne45inds.txt
awk 'BEGIN{ FS = OFS = "\t" } { print $0, "45" }' ./results/yellow99/Ne45inds.txt > tmp && mv tmp ./results/yellow99/Ne45inds.txt

for i in {1..50}; do
 awk 'NR==3' ./results/yellow99/Output_Ne_60inds.$i | cut -f 2
done  >> ./results/yellow99/Ne60inds.txt
awk 'BEGIN{ FS = OFS = "\t" } { print $0, "60" }' ./results/yellow99/Ne60inds.txt > tmp && mv tmp ./results/yellow99/Ne60inds.txt

for i in {1..50}; do
 awk 'NR==3' ./results/yellow99/Output_Ne_75inds.$i | cut -f 2 
done >> ./results/yellow99/Ne75inds.txt
awk 'BEGIN{ FS = OFS = "\t" } { print $0, "75" }' ./results/yellow99/Ne75inds.txt > tmp && mv tmp ./results/yellow99/Ne75inds.txt

cat ./results/yellow99/Ne15inds.txt ./results/yellow99/Ne30inds.txt ./results/yellow99/Ne45inds.txt ./results/yellow99/Ne60inds.txt ./results/yellow99/Ne75inds.txt > ./results/yellow99/NeIndividuals.txt
```
The file ```./results/yellow99/NeIndividuals.txt``` includes *N*<sub>e</sub> estimates in the first column and the subset of individuals considered in the second column. We can add the header or not.  
The file ```NeVsInds.txt``` looks like this:  
```
Ne	inds
10054.7	15
847.069	15
6232.01	15
# ...
23444.2	30
3729.25	30
15170.3	30
# ... end of file:
3415.5	75
3605.97	75
3513.67	75
```
Plotting the results in R:
```
library (ggplot2)
library (scales)
indsubsets <- read.delim("NeVsInds.txt")
p1 <-ggplot(indsubsets, aes(x=inds, y=Ne), color = "#00C08D", fill = "#00C08D") + 
   stat_summary(fun.data = median_hilow, mapping = aes(group=inds), 
    geom = "pointrange", color = "#00C08D", size = 3, linewidth = 1, fun.args = list(conf.int = 0.95)) +
	theme_test() + 
	theme(axis.title.x = element_text(size=24, face="bold"), axis.title.y = element_text(size=24, face="bold"), 
	axis.text.x = element_text(size=24, colour = "black"),axis.text.y = element_text(size=24, colour = "black")) + 
	scale_x_continuous(limits=c(10,75), breaks=seq(0,75,15)) + scale_y_continuous(labels=comma) +
	labs(x ="Number of individuals sampled", y = "Ne (Geometric mean)")

# to display zoom on point *N*<sub>e</sub> 
p2 <- ggplot(indsubsets, aes(x=inds, y=Ne), color = "#00C08D", fill = "#00C08D") + 
   stat_summary(fun = median, mapping = aes(group=inds), 
    geom = "point", color = "#00C08D", size = 10) +
	theme_test() + 
	theme(axis.title.x = element_text(size=16, face="bold"), axis.title.y = element_text(size=16, face="bold"), 
	axis.text.x = element_text(size=16, colour = "black"),axis.text.y = element_text(size=16, colour = "black")) + 
	scale_x_continuous(limits=c(10,75), breaks=seq(0,75,15)) + scale_y_continuous(labels=comma) +
	labs(x ="Number of individuals sampled", y = "Ne (Geometric mean)") + ylim(0, 5000)

# to create an insert within the main plot.	
p1 + annotation_custom(ggplotGrob(p2), xmin = 40, xmax = 75, 
                       ymin = 30000, ymax = 75000)		
```
#### 3. Influence of population structure on *N*<sub>e</sub> estimation
We first generate a list of SNPs that we need to subsample the dataset):
```sh
zgrep -v "#" ./data/armeniaca.SNPs.vcf.gz | cut -f 1,2 > ./snpslist/SNPlist
shuf ./snpslist/SNPlist | head -n 3500000 | sort > ./snpslist/3.5M.SNPs

module load vcftools/0.1.16
vcftools --gzvcf ./data/armeniaca.SNPs.vcf.gz --positions ./snpslist/3.5M.SNPs --recode --out ./data/armeniaca.3.5SNPs
```
this generates the file ```armeniaca.3.5SNPs.recode.vcf``` that we will need later.  
Then we resample individuals for each Q-value group (for both the Southern and the Northern gene pool), repeating the resampling procedure 50 times:
```sh
# we resample 77 individuals in the Northern gene pool (the minimum number of individual in a Q-value group):
for i in {1..50}; do
 shuf -n 77 ./indlist/yellow99 > ./indlist/yellow99.$i
done

for i in {1..50}; do
 shuf -n 77 ./indlist/yellow95 > ./indlist/yellow95.$i
done

for i in {1..50}; do
 shuf -n 77 ./indlist/yellow90 > ./indlist/yellow90.$i
done

for i in {1..50}; do
 shuf -n 77 ./indlist/yellow80 > ./indlist/yellow80.$i
done

for i in {1..50}; do
 shuf -n 77 ./indlist/yellow70 > ./indlist/yellow70.$i
done

# we resample 21 individuals in the Northern gene pool (the minimum number of individual in a Q-value group):
for i in {1..50}; do
 shuf -n 21 ./indlist/red99 > ./indlist/red99.$i
done

for i in {1..50}; do
 shuf -n 21 ./indlist/red95 > ./indlist/red95.$i
done

for i in {1..50}; do
 shuf -n 21 ./indlist/red90 > ./indlist/red90.$i
done

for i in {1..50}; do
 shuf -n 21 ./indlist/red80 > ./indlist/red80.$i
done

for i in {1..50}; do
 shuf -n 21 ./indlist/red70 > ./indlist/red70.$i
done
```
We also resample the entire dataset (255 individuals), to see the influence of mixing two different gene pools.
```sh
for i in {1..50}; do
 shuf -n 77 ./indlist/all > ./indlist/all.$i
done
```
As the scripts above have only generated lists of individuals, now we generate the actual SNPs dataset, in ped and map formats.
```sh
module load vcftools/0.1.16

for i in {1..50}; do 
 vcftools --vcf ./data/armeniaca.3.5SNPs.recode.vcf --plink --keep ./indlist/yellow99.$i --chrom-map ./data/chrom-map --out ./data/yellow99.$i
done

for i in {1..50}; do 
 vcftools --vcf ./data/armeniaca.3.5SNPs.recode.vcf --plink --keep ./indlist/yellow95.$i --chrom-map ./data/chrom-map --out ./data/yellow95.$i
done

for i in {1..50}; do 
 vcftools --vcf ./data/armeniaca.3.5SNPs.recode.vcf --plink --keep ./indlist/yellow90.$i --chrom-map ./data/chrom-map --out ./data/yellow90.$i
done

for i in {1..50}; do 
 vcftools --vcf ./data/armeniaca.3.5SNPs.recode.vcf --plink --keep ./indlist/yellow80.$i --chrom-map ./data/chrom-map --out ./data/yellow80.$i
done

for i in {1..50}; do 
 vcftools --vcf ./data/armeniaca.3.5SNPs.recode.vcf --plink --keep ./indlist/yellow70.$i --chrom-map ./data/chrom-map --out ./data/yellow70.$i
done

for i in {1..50}; do 
 vcftools --vcf ./data/armeniaca.3.5SNPs.recode.vcf --plink --keep ./indlist/red99.$i --chrom-map ./data/chrom-map --out ./data/red99.$i
done

for i in {1..50}; do 
 vcftools --vcf ./data/armeniaca.3.5SNPs.recode.vcf --plink --keep ./indlist/red95.$i --chrom-map ./data/chrom-map --out ./data/red95.$i
done

for i in {1..50}; do 
 vcftools --vcf ./data/armeniaca.3.5SNPs.recode.vcf --plink --keep ./indlist/red90.$i --chrom-map ./data/chrom-map --out ./data/red90.$i
done

for i in {1..50}; do 
 vcftools --vcf ./data/armeniaca.3.5SNPs.recode.vcf --plink --keep ./indlist/red80.$i --chrom-map ./data/chrom-map --out ./data/red80.$i
done

for i in {1..50}; do 
 vcftools --vcf ./data/armeniaca.3.5SNPs.recode.vcf --plink --keep ./indlist/red70.$i --chrom-map ./data/chrom-map --out ./data/red70.$i
done


for i in {1..50}; do 
 vcftools --vcf ./data/armeniaca.3.5SNPs.recode.vcf --plink --keep ./indlist/all.$i --chrom-map ./data/chrom-map --out ./data/all.$i
done
```
We then modify ped and map files to the format required by GONE, as done previously, for example, for the Q-value 70% dataset:
```sh
for i in {1..50}; do
	cut -f7- ./data/yellow70.$i.ped > ./data/yellow70.${i}Geno
done

for i in {1..50}; do
	paste ./indlist/77 ./data/yellow70.${i}Geno > ./data/yellow70.$i.ped
done

rm ./data/*Geno
```
We then move the ped and map files to the respective folders where we will run GONE,
for example in ```./results/yellow70/```:
```sh
for i in {1..50}; do
 bash script_GONE.sh yellow70.$i   
done
```
Plotting the results; let's first extract last-generation *N*<sub>e</sub> estimates for the Northern gene pool datasets:
```sh
for i in {1..50}; do
 awk 'NR==3' ./results/yellow99/Output_Ne_yellow99.$i 
done | cut -f 2 > ./results/yellow99/Ne99.txt

for i in {1..50}; do
 awk 'NR==3' ./results/yellow95/Output_Ne_yellow95.$i 
done | cut -f 2 > ./results/yellow95/Ne95.txt

for i in {1..50}; do
 awk 'NR==3' ./results/yellow90/Output_Ne_yellow90.$i 
done | cut -f 2 > ./results/yellow90/Ne90.txt

for i in {1..50}; do
 awk 'NR==3' ./results/yellow80/Output_Ne_yellow80.$i 
done | cut -f 2 > ./results/yellow80/Ne80.txt

for i in {1..50}; do
 awk 'NR==3' ./results/yellow70/Output_Ne_yellow70.$i 
done | cut -f 2 > ./results/yellow70/Ne70.txt

# to add two columns with labels:
sed $'s/$/\tyellow99\tyellow/' ./results/yellow99/Ne99.txt > ./results/yellow99/Ne99yellow.txt
sed $'s/$/\tyellow95\tyellow/' ./results/yellow95/Ne95.txt > ./results/yellow95/Ne95yellow.txt
sed $'s/$/\tyellow90\tyellow/' ./results/yellow90/Ne90.txt > ./results/yellow90/Ne90yellow.txt
sed $'s/$/\tyellow80\tyellow/' ./results/yellow80/Ne80.txt > ./results/yellow80/Ne80yellow.txt
sed $'s/$/\tyellow70\tyellow/' ./results/yellow70/Ne70.txt > ./results/yellow70/Ne70yellow.txt

cat ./results/yellow99/Ne99yellow.txt ./results/yellow95/Ne95yellow.txt ./results/yellow90/Ne90yellow.txt ./results/yellow80/Ne80yellow.txt ./results/yellow70/Ne70yellow.txt > ./results/NeYellow.txt
```
Then for the Southern gene pool dataset:
```sh
for i in {1..50}; do
 awk 'NR==3' ./results/red99/Output_Ne_red99.$i 
done | cut -f 2 > ./results/red99/Ne99.txt

for i in {1..50}; do
 awk 'NR==3' ./results/red95/Output_Ne_red95.$i 
done | cut -f 2 > ./results/red95/Ne95.txt

for i in {1..50}; do
 awk 'NR==3' ./results/red90/Output_Ne_red90.$i 
done | cut -f 2 > ./results/red90/Ne90.txt

for i in {1..50}; do
 awk 'NR==3' ./results/red80/Output_Ne_red80.$i 
done | cut -f 2 > ./results/red80/Ne80.txt

for i in {1..50}; do
 awk 'NR==3' ./results/red70/Output_Ne_red70.$i 
done | cut -f 2 > ./results/red70/Ne70.txt

# to add two columns with labels:
sed $'s/$/\tred99\tred/' ./results/red99/Ne99.txt > ./results/red99/Ne99red.txt
sed $'s/$/\tred95\tred/' ./results/red95/Ne95.txt > ./results/red95/Ne95red.txt
sed $'s/$/\tred90\tred/' ./results/red90/Ne90.txt > ./results/red90/Ne90red.txt
sed $'s/$/\tred80\tred/' ./results/red80/Ne80.txt > ./results/red80/Ne80red.txt
sed $'s/$/\tred70\tred/' ./results/red70/Ne70.txt > ./results/red70/Ne70red.txt

cat ./results/red99/Ne99red.txt ./results/red95/Ne95red.txt ./results/red90/Ne90red.txt ./results/red80/Ne80red.txt ./results/red70/Ne70red.txt > ./results/NeRed.txt 
```
We concatenate results for both gene pools:
```sh
cat ./results/NeYellow.txt ./results/NeRed.txt > ./results/NeVsStructure.txt
```
We also want to get the *N*<sub>e</sub> estimates for the last 25 generations:
```sh
for i in {1..50}; do
 awk 'NR>=3 && NR<=27' ./results/yellow99/Output_Ne_yellow99.$i 
done | cut -f 1,2 > ./results/yellow99/Ne99_25gen.txt

for i in {1..50}; do
  awk 'NR>=3 && NR<=27' ./results/yellow95/Output_Ne_yellow95.$i 
done | cut -f 1,2 > ./results/yellow95/Ne95_25gen.txt

for i in {1..50}; do
  awk 'NR>=3 && NR<=27' ./results/yellow90/Output_Ne_yellow90.$i 
done | cut -f 1,2 > ./results/yellow90/Ne90_25gen.txt

for i in {1..50}; do
  awk 'NR>=3 && NR<=27' ./results/yellow80/Output_Ne_yellow80.$i 
done | cut -f 1,2 > ./results/yellow80/Ne80_25gen.txt

for i in {1..50}; do
  awk 'NR>=3 && NR<=27' ./results/yellow70/Output_Ne_yellow70.$i 
done | cut -f 1,2 > ./results/yellow70/Ne70_25gen.txt

# to add two columns with labels:

sed $'s/$/\tyellow99\tyellow/' ./results/yellow99/Ne99_25gen.txt > ./results/yellow99/Ne99yellow_25gen.txt
sed $'s/$/\tyellow95\tyellow/' ./results/yellow95/Ne95_25gen.txt > ./results/yellow95/Ne95yellow_25gen.txt
sed $'s/$/\tyellow90\tyellow/' ./results/yellow90/Ne90_25gen.txt > ./results/yellow90/Ne90yellow_25gen.txt
sed $'s/$/\tyellow80\tyellow/' ./results/yellow80/Ne80_25gen.txt > ./results/yellow80/Ne80yellow_25gen.txt
sed $'s/$/\tyellow70\tyellow/' ./results/yellow70/Ne70_25gen.txt > ./results/yellow70/Ne70yellow_25gen.txt

cat ./results/yellow99/Ne99yellow_25gen.txt ./results/yellow95/Ne95yellow_25gen.txt ./results/yellow90/Ne90yellow_25gen.txt ./results/yellow80/Ne80yellow_25gen.txt ./results/yellow70/Ne70yellow_25gen.txt > ./results/NeYellow_25gen.txt
```
We do the same for the Southern gene pool, then we concatenate the results in:
```sh
cat ./results/NeYellow_25gen.txt ./results/NeRed_25gen.txt > ./results/NeVsStructure_25gen.txt
```
For the total dataset, we keep results in separate files, but we can also combine them in the same plot if needed (not shown):
```sh
# *N*<sub>e</sub> during the last generation across 50 datasets:
for i in {1..50}; do
 awk 'NR==3' ./results/all/Output_Ne_all.$i 
done | cut -f 2 > ./results/all/allNe.txt
# last 25 generations across 50 datasets:
for i in {1..50}; do
 awk 'NR>=3 && NR<=27' ./results/all/Output_Ne_all.$i 
done | cut -f 1,2 > ./results/all/all_Ne25gen.txt
```
The file ```NeVsStructure.txt``` looks like this:  
```
Ne	group	pop
3399.45	yellow99	yellow
3673.48	yellow99	yellow
3454.42	yellow99	yellow
# ...
6856.61	yellow70	yellow
5603.9	yellow70	yellow
3774.62	yellow70	yellow
178.688	red99	red
170.765	red99	red
196.732	red99	red
# ... end of file:
4118.63	all	red
1233.03	all	red
1784.73	all	red
```
The file ```NeVsStructure_25gen.txt``` looks like this:  
```
gen	Ne	group	pop
1	3399.45	p99	yellow
2	3399.45	p99	yellow
3	3399.45	p99	yellow
# ...
24	4471.06	p80	yellow
25	4408.21	p80	yellow
1	6023.96	p80	yellow
2	6023.96	p80	yellow
3	6023.96	p80	yellow
# ... end of file:
24	5070.48	allR	red
25	5087.75	allR	red
```
Plotting results in R:
```
library (ggplot2)
library (scales)
library(RColorBrewer)
library(viridis)

structure <- read.delim("NeVsStructure.txt")

t <-ggplot(structure, aes(x=group, y=Ne, colour = pop)) + 
   stat_summary(fun.data = median_hilow, mapping = aes(group=group), geom = "pointrange", 
   size = 3, linewidth = 2, fun.args = list(conf.int = 0.95)) +
   theme_test() + scale_y_continuous(labels=comma) +
   theme(axis.title.x = element_text(size=24, face="bold"), axis.title.y = element_text(size=24, face="bold"), 
   axis.text.x = element_text(size=24, colour = "black"),axis.text.y = element_text(size=24, colour = "black"), legend.title = element_text(size=14), legend.text = element_text(size=14)) + 
    labs(x ="Proportion of genetic membership", y = "Ne (Geometric mean)")

color_scale <- scale_color_manual(values = c("#FB8072","#FFFFB3"))
t <- t + color_scale
w <- t + facet_wrap(. ~ pop, scales="free_x") + theme(strip.background = element_rect(colour="black", fill="white", 
                                       linewidth=0.5, linetype="solid"), strip.text.x = element_text(size=24, face="bold"))
```	
In the last 25 generations:
```
structure25 <- read.delim("NeVsStructure_25gen.txt")
s <- ggplot(structure25, aes(x=gen, y=Ne, colour = group)) + 
  geom_line(stat="summary", fun.data=median_hilow, aes(x=gen, y=Ne, group=group), 
            size = 2, linetype="solid") +   
  scale_color_viridis_d(option="magma") + 
  theme_test() + scale_y_continuous(labels=comma) +
  theme(axis.title.x = element_text(size=24, face="bold"), 
        axis.title.y = element_text(size=24, face="bold"), 
        axis.text.x = element_text(size=24, colour = "black"),
        axis.text.y = element_text(size=24, colour = "black"), legend.title = element_text(size=14), legend.text = element_text(size=14)) + 
  labs(x ="Generation", y = "Ne (Geometric mean)") + ylim(0,6000)

v <- s + facet_wrap(. ~ pop, scales="free_x") + 
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        linewidth=0.5, linetype="solid"), 
        strip.text.x = element_text(size=24, face="bold"))

library(patchwork)
w / v
```

#### 4. Effect of using genomic scaffolds instead of chromosomes on *N*<sub>e</sub> estimation
We work on 45 random individuals from the Northern gene pool dataset, using a subset of 3.5 M SNPs:
```sh
module load vcftools/0.1.16
module load plink/1.90
vcftools --gzvcf ./data/armeniaca.SNPs.vcf.gz --plink --keep ./indlist/indset45.1 --chrom-map ./data/chrom-map --out ./data/yellow99.45inds.1
# As above, we get ped and map files
plink --file ./data/yellow99.45inds.1 --recode --extract ./snpslist/yellow3.5subset.snps --out ./data/45inds.8chr
# Then we need to edit the ped file for GONE
tr ' ' \\t < ./data/45inds.8chr.ped > ./data/test8chr
mv ./data/test8chr ./data/45inds.8chr.ped
cut -f7- ./data/45inds.8chr.ped > ./data/45inds.8chr.Geno
paste ./indlist/45 ./data/45inds.8chr.Geno > ./data/45inds.8chr.ped
```
In this analysis, we want to try to assign SNPs to a number of chromosomes greater than the true one (8 chromosomes), simulating what happens when we have information about genomic scaffolds or short linkage blocks but we don't have SNPs mapped to full chromosomes. We want to assign SNPs to either:  
-16 chromosomes  
-32 chromosomes  
-64 chromosomes  
-128 chromosomes  
Approximately, we want to get for each of these "assumed" chromosomes:  
3500000/8 = 437,500 SNPs  
3500000/16 = 218,750 SNPs  
3500000/32 = 109,375 SNPs  
3500000/64 = 54,688 SNPs  
3500000/128 = 27,344 SNPs  
Of course, these are averaged values, but in reality some chromosomes will have more and some fewer SNPs.

We then modify the file ```45inds.8chr.map``` to redistribute the SNPs to a progressively higher number of chromosomes.  
The file includes 3500000 rows and four columns (tab separated). The first row looks like this:  
"1	chr1:85	0	85"  
The last row looks like this:    
"8	chr8:20012701	0	20012701"  
We create a script that replaces the number in the first column and the number next to the string "chr" (before the ":") starting from row 218750, by progressively increasing the number from 1 to 16, for each and every 218750 rows. The script is also available in the scripts folder of this repository as "replace_chromosomes.sh" (every attempt to improve it is very welcome!).
```sh
# setting the input and output file paths
input_file="45inds.8chr.map"
output_file="45inds.16chr.map"
# setting the starting row and increment
start_row=218750
increment=1
# loop through each line in the input file
line_number=0
while read line; do
    line_number=$((line_number + 1))

    # checking if this line should have the numbers replaced
    if [ $line_number -ge $start_row ]; then
        # calculating the new values for the first two columns
        new_first_col=$(( (line_number - start_row) / start_row + increment + 1 ))
        new_second_col="chr$new_first_col:"
        
        # replacing the first column and the number after "chr"
        modified_line=$(echo "$line" | sed "s/^[0-9]\+\s/""$new_first_col""\t/")
        modified_line=$(echo "$modified_line" | sed "s/chr[0-9]\+:/""$new_second_col""/")

        # appending the modified line to the output file
        echo "$modified_line" >> "$output_file"
    else
        # appending the original line to the output file
        echo "$line" >> "$output_file"
    fi
done < "$input_file"
```

Then I do the same for the other redistributions (32 to 128 chromosomes), every time changing the script according to the values required in this parts of the script:
```
input_file="45inds.8chr.map"
output_file="45inds.16chr.map"
```
and
```
#setting the starting row and increment
start_row=218750
increment=1
```
At the end, we will have the following chromosome maps:  
```
45inds.8chr.map
45inds.16chr.map
45inds.32chr.map
45inds.64chr.map
45inds.128chr.map
```
Warning: 45inds.16chr.map and 45inds.32chr.map have one chromosome more in the last line (17 and 33 respectively); I just replace the values manually.  
We do not change the ped file, but we rename it every time, as GONE needs ped and map file to have the same prefix.  
```sh
cp 45inds.8chr.ped 45inds.16chr.ped
cp 45inds.8chr.ped 45inds.32chr.ped
cp 45inds.8chr.ped 45inds.64chr.ped
cp 45inds.8chr.ped 45inds.128chr.ped
```
We then run GONE in the appropriate folder (we cannot run more than one analysis in the same folder, as GONE overwrites output files), e.g.:
```sh
bash script_GONE.sh 45inds.16chr
```
To prepare the results for plotting:
```sh
awk 'NR>=3 && NR<=27' ./results/yellow99/Output_Ne_45inds.8chr > ./results/yellow99/Ne8chr.txt
awk 'NR>=3 && NR<=27' ./results/yellow99/Output_Ne_45inds.16chr > ./results/yellow99/Ne16chr.txt
awk 'NR>=3 && NR<=27' ./results/yellow99/Output_Ne_45inds.32chr > ./results/yellow99/Ne32chr.txt
awk 'NR>=3 && NR<=27' ./results/yellow99/Output_Ne_45inds.64chr > ./results/yellow99/Ne64chr.txt
awk 'NR>=3 && NR<=27' ./results/yellow99/Output_Ne_45inds.128chr > ./results/yellow99/Ne128chr.txt

cat ./results/yellow99/Ne8chr.txt ./results/yellow99/Ne16chr.txt ./results/yellow99/Ne32chr.txt ./results/yellow99/Ne64chr.txt ./results/yellow99/Ne128chr.txt > ./results/yellow99/NeByChrom.txt
```
We can modify the file "NeByChrom.txt" manually to include information about chromosomes. The file will look like this:
```
gen	Ne	chrom	facet
1	4566.61	1-8chr	facet3
2	4566.61	1-8chr	facet3
3	4566.61	1-8chr	facet3
# ...
1	13229.1	3-32chr	facet3
2	13229.1	3-32chr	facet3
3	13229.1	3-32chr	facet3
# ...
23	761005	5-128chr	facet1
24	769090	5-128chr	facet1
25	784894	5-128chr	facet1
```
Facets facilitate displaying very different *N*<sub>e</sub> values over the y axis, see below.  
In R:
```
library(ggplot2)
library(viridis)
chromsubs <- read.delim("NeByChrom.txt")

c <-ggplot(chromsubs, aes(x=gen, y=Ne, group=chrom)) + 
  geom_line(aes(color=chrom), linewidth = 1) +
  geom_point(aes(fill=chrom), color="grey20", size = 5, alpha = 0.8, shape=21)+
    scale_fill_viridis(discrete=TRUE) + theme_bw() + theme_classic() +
    scale_color_viridis(discrete=TRUE) + 
	theme(axis.title.x = element_text(size=24, face="bold"), axis.title.y = element_text(size=24, face="bold"), 
	axis.text.x = element_text(size=24, colour = "black"),axis.text.y = element_text(size=24, colour = "black"),
	legend.title = element_text(size=20, colour = "black"),legend.text = element_text(size=20, colour = "black")) + 
	labs(x ="Generation", y = "Ne (Geometric mean)") + scale_y_continuous(labels=comma) 
	
c + facet_grid(facet ~ ., scales="free_y") + theme(strip.text.y = element_blank())
```

#### 5. Influence of missing data on *N*<sub>e</sub> estimation
Starting from the Northern gene pool dataset with 77 individuals and 3.5 million SNPs:
```sh
module load vcftools/0.1.16
vcftools --vcf ./data/armeniaca.3.5SNPs.recode.vcf --keep ./indlist/yellow99 --recode --out ./data/yellow3.5
```
we'll get the file ```yellow3.5.recode.vcf```. Now we want to randomly input missing data using a custom python v3.10 script (script name "missing.py")

```py
import random

input_file = "yellow3.5.recode.vcf"
output_file = "yellow3.5.recodeMissing20.vcf"
missing_rate = 0.2

with open(input_file, 'r') as input_vcf, open(output_file, 'w') as output_vcf:
    for line in input_vcf:
        if line.startswith("#"):  # Write header lines
            output_vcf.write(line)
        else:
            fields = line.split("\t")
            genotypes = fields[9:]  # Get genotype fields for all individuals
            for i in range(len(genotypes)):
                if random.random() < missing_rate:
                    genotypes[i] = "./."  # Introduce missing genotype
            fields[9:] = genotypes
            output_vcf.write("\t".join(fields))
```

We want to test what happens when we change the missing data proportion from 20% to 80%, and therefore we edit the script to generate the datasets accordingly (i.e. missing_rate, ranging from 0.2 to 0.8 and output_file name).  
To check that everything is correct, we can calculate the proportion of missing data per individual for each file generated:
```sh
module load vcftools/0.1.16
vcftools --vcf yellow3.5.recodeMissing20.vcf --missing-indv --out miss20
vcftools --vcf yellow3.5.recodeMissing40.vcf --missing-indv --out miss40
vcftools --vcf yellow3.5.recodeMissing60.vcf --missing-indv --out miss60
vcftools --vcf yellow3.5.recodeMissing80.vcf --missing-indv --out miss80
vcftools --vcf yellow3.5.recode.vcf --missing-indv --out  miss3.5recode77inds
```
For each of the dataset we create ped and map files:
```sh
vcftools --vcf ./data/yellow3.5.recodeMissing20.vcf --plink --chrom-map ./data/chrom-map --out ./data/yellow3.5Missing20
vcftools --vcf ./data/yellow3.5.recodeMissing40.vcf --plink --chrom-map ./data/chrom-map --out ./data/yellow3.5Missing40
vcftools --vcf ./data/yellow3.5.recodeMissing60.vcf --plink --chrom-map ./data/chrom-map --out ./data/yellow3.5Missing60
vcftools --vcf ./data/yellow3.5.recodeMissing80.vcf --plink --chrom-map ./data/chrom-map --out ./data/yellow3.5Missing80
```
We then modify them and run GONE, for example:
```sh
cut -f7- ./data/yellow3.5Missing20.ped > ./data/yellow3.5Missing20Geno
paste ./indlist/77 ./data/yellow3.5Missing20Geno > ./data/yellow3.5Missing20.ped
bash script_GONE.sh yellow3.5Missing20
```
Notice that *N*<sub>e</sub> estimation in GONE failed for the datasets yellow3.5Missing60 and yellow3.5Missing80.  
Let's extract *N*<sub>e</sub> estimates in the last 25 generations and prepare the results to plot:
```sh
awk 'NR>=3 && NR<=27' ./results/yellow99/Output_Ne_yellow99.subset3.5M1 | cut -f 1,2 >> ./results/yellow99/NeMissing.txt
awk 'NR>=3 && NR<=27' ./results/yellow99/Output_Ne_yellow3.5Missing20 | cut -f 1,2 >> ./results/yellow99/NeMissing.txt
awk 'NR>=3 && NR<=27' ./results/yellow99/Output_Ne_yellow3.5Missing40 | cut -f 1,2 >> ./results/yellow99/NeMissing.txt
```
The file NeMissing.txt will look like:
```
gen	Ne	missing
1	3744.29	3-8% 
2	3744.29	3-8%
3	3744.29	3-8%
# ...
1	1.46E+06	2-20%
2	1.46E+06	2-20%
3	1.46E+06	2-20%
# ...
1	4.14E+06	1-40%
2	4.14E+06	1-40%
3	4.14E+06	1-40%
```
Where 8% in the first rows corresponds to the original dataset (true percentage of missing data = 8%).  
In R:
```
library(ggplot2)
library(viridis)

miss <- read.delim("NeVsMiss.txt")

m <-ggplot(miss, aes(x=gen, y=Ne, group=missing)) + 
  geom_line(aes(color=missing), linewidth = 1) +
  geom_point(aes(fill=missing, shape = missing), color="grey20", size = 5, alpha = 0.8) +
    scale_fill_viridis(discrete=TRUE) + theme_bw() + theme_classic() +
    scale_color_viridis(discrete=TRUE) + 
	scale_shape_manual(values=c(21, 22, 23)) +
	theme(axis.title.x = element_text(size=24, face="bold"), axis.title.y = element_text(size=24, face="bold"), 
	axis.text.x = element_text(size=24, colour = "black"),axis.text.y = element_text(size=24, colour = "black"),
	legend.title = element_text(size=20, colour = "black"),legend.text = element_text(size=20, colour = "black")) + 
	labs(x ="Generation", y = "Ne") + scale_y_continuous(labels=comma) 
	
m + facet_grid(missing ~ ., scales="free_y") + theme(strip.text.y = element_blank())
```

## *Symphonia globulifera*
#### Preparing the datasets and running GONE
We use the dataset from [Schmitt et al. 2021](https://onlinelibrary.wiley.com/doi/10.1111/mec.16116), also available [here](https://doi.org/10.5281/zenodo.4727831). We need to group individuals based on their Q-values (according to the analyses in Schmitt et al. 2021). 
```sh
# obtaining Q-values from original data:
paste symcapture.all.biallelic.snp.filtered.nonmissing.paracou.fam symcapture.all.biallelic.snp.filtered.nonmissing.paracou.3.Q > Individuals_Q
# extracting individuals based on a Q-value >= 95%:
awk '$7 >= 0.95' Individuals_Q > Species1
awk '$8 >= 0.95' Individuals_Q > Species2
awk '$9 >= 0.95' Individuals_Q > Species3
# we can also count them:
awk '$7 >= 0.95' Individuals_Q | wc -l
# 228
awk '$8 >= 0.95' Individuals_Q | wc -l
# 107
awk '$9 >= 0.95' Individuals_Q | wc -l
# 30
mkdir Ind_lists
mv Species1 ./Ind_lists
mv Species2 ./Ind_lists
mv Species3 ./Ind_lists
mv Individuals_Q ./Ind_lists
# getting only the list of individuals
cut -f 2 ./Ind_lists/Species1 > Ind_lists/Inds1
cut -f 2 ./Ind_lists/Species2 > Ind_lists/Inds2
cut -f 2 ./Ind_lists/Species3 > Ind_lists/Inds3
```
We now want to select the longest contigs (we do not have SNPs mapping or chromosome information for this species):
```sh
module load bcftools/1.13
bcftools view -H symcapture.all.biallelic.snp.filtered.nonmissing.paracou.vcf | cut -f 1 | sort | uniq -c > ContigList
```
Manual step: I sorted the contigs obtained in the file by number of SNPs and extracted the 125 contigs with the largest number of SNPs. Then I created the vcf file with the script "extract_contigs_Symphonia.sh" (see scripts folder). The dataset generated is ```symcapture.all.biallelic.snp.filtered.nonmissing.125Contigs.paracou```. We now need to obtain ped and map files for the three gene pools (Species 1, Species 2 and Species 3). 
```sh 
# To get the map file, we need to create a chromosome map in which we assume that contigs = chromosomes.
module load bcftools/1.13
bcftools view -H symcapture.all.biallelic.snp.filtered.nonmissing.125Contigs.paracou.recode.vcf | cut -f 1 | uniq | awk '{print $0"\t"$0}' > chrom-map
# replacing contig names with numbers, as required by GONE:
nl chrom-map > cc
cut -f 1 cc > c1
cut -f 2 cc > c2
paste c2 c1 > chrom-map
# now preparing the files in plink format for GONE
module load vcftools/0.1.16
vcftools --vcf symcapture.all.biallelic.snp.filtered.nonmissing.125Contigs.paracou.recode.vcf --keep ./Ind_lists/Inds1 --plink --chrom-map chrom-map --out Species1
vcftools --vcf symcapture.all.biallelic.snp.filtered.nonmissing.125Contigs.paracou.recode.vcf --keep ./Ind_lists/Inds2 --plink --chrom-map chrom-map --out Species2
vcftools --vcf symcapture.all.biallelic.snp.filtered.nonmissing.125Contigs.paracou.recode.vcf --keep ./Ind_lists/Inds3 --plink --chrom-map chrom-map --out Species3
# then we modify the ped file according to GONE requirements:
cut -f7- Species3.ped > Sp3
paste Species3 Sp3 > Species3.ped
cut -f7- Species2.ped > Sp2
paste Species2 Sp2 > Species2.ped
cut -f7- Species1.ped > Sp1
paste Species1 Sp1 > Species1.ped
# where Species1 Species2 and Species3 are tab delimited files with 6 columns and a number of row equal to the number of individuals in the respective dataset for the three species. Their first row looks like this:
# 1 IND1 0 0 0 -9 A T (...etc.)
```
Once the datasets are ready, we can move them to the appropriate directory and run GONE (remember that GONE overwrites output files, so one analysis can be run at a given time in one directory).
```sh
bash script_GONE.sh Species1
```
As explained in Gargiulo et al. 2023, the analyses run on each gene pool from *Symphonia* returned the error "Too few SNPs".

## *Mercurialis annua*
#### Preparing the datasets and running GONE
We use the dataset from [González-Martínez et al. (2017)](https://www.sciencedirect.com/science/article/pii/S0960982217308655). For *Mercurialis*, we do not have information about the physical location of SNPs on chromosomes . As previously done for *Symphonia*, we extract information about the number of SNPs in each contig and sort the list of contigs manually (in excel), to extract the identity of the contigs with the largest number of SNPs. We also group individuals based on the gene pools found in González-Martínez et al. (2017) and divide them in different datasets: "Core", "Atlantic" and "Mediterranean". For the final analyses, we extract the 48 contigs that have a number of SNPs > 500 (to avoid the issue of having too few SNPs, observed above with *Symphonia*) and generate a new vcf file (see script "extract_contigs_Mercurialis.sh" in the scripts folder). The newly obtained dataset is ```annua.48Contigs.recode.vcf```. We then obtain ped and map files as follows: 
```sh
module load vcftools/0.1.16
vcftools --vcf annua.48Contigs.recode.vcf --keep ./Core --plink --chrom-map chrom-map48 --out CoreReduced
vcftools --vcf annua.48Contigs.recode.vcf --keep ./Atlantic --plink --chrom-map chrom-map48 --out AtlanticReduced
vcftools --vcf annua.48Contigs.recode.vcf --keep ./Mediter --plink --chrom-map chrom-map48 --out MediterReduced
# we create the tab-delimited files "C", "A" and "M", which contains the correct format required for GONE inputs (as explained for Symphonia and Prunus above) 
cut -f7- CoreReduced.ped > CoreReducedGeno
paste C CoreReducedGeno > CoreReduced.ped

cut -f7- AtlanticReduced.ped > AtlanticReducedGeno
paste A AtlanticReducedGeno > AtlanticReduced.ped

cut -f7- MediterReduced.ped > MediterReducedGeno
paste M MediterReducedGeno > MediterReduced.ped

# we then run GONE for the three datasets, starting with:
bash script_GONE.sh CoreReduced
```

## *Fagus sylvatica*
#### Preparing the dataset
Starting from [this dataset](https://doi.org/10.57745/FJRYI1), which includes 167 individual genotypes and 406 genomic scaffolds (SNPs not mapped to chromosomes), we want to evaluate how *N*<sub>e</sub> estimation in GONE changes depending on the percentage of missing data. We first extract the first 27 genomic scaffolds with the largest number of SNPs (see script "extract_scaffolds_Fagus.sh"). However, we exclude scaffolds with more than one million SNPs, to prevent GONE from crashing (the maximum number of SNPs accepted by GONE is 10 millions). We then create a chromosome map, to replace scaffold names by numbers, as required by GONE:
```sh
nl chrom-map > cc
cut -f 1 cc > c1
cut -f 2 cc > c2
paste c2 c1 > chrom-map
rm cc c1 c2
```
#### Influence of missing data on *N*<sub>e</sub> estimation 
We can use vcftools to obtain missing data proportions per individual:
```sh
module load vcftools/0.1.16
vcftools --vcf 27scaffolds.recode.vcf --missing-indv --out miss
```
Most individuals have a missing data proportion > 95%. We exclude these and keep 67 individuals with missing data proportion < 95% for the following analyses.
We then generate 150 datasets by permuting 35 individuals (out of the 67 individuals retained):
```sh
module load vcftools/0.1.16
for i in {1..150}; do
 shuf -n 35 ./indlist/all > ./indlist/indset$i
done
# generating datasets from the lists of individuals above:
for i in {1..150}; do
 vcftools --gzvcf ./data/27scaffolds.recode.vcf.gz --keep ./indlist/indset$i --recode --stdout | gzip > ./data/Fagus$i.recode.vcf.gz
done

# computing missing data proportion in each dataset:
for i in {1..150}; do
 vcftools --gzvcf ./data/Fagus$i.recode.vcf.gz --missing-indv --out ./data/Fagus$i
done

# extracting missing data proportions per individual and printing the results in missingdata.txt
for i in {1..150}; do
 cut -f 1,5 ./data/Fagus$i.imiss | grep -v "INDV" >> ./data/missingdata.txt
done

# generating input files for GONE:
for i in {1..150}; do
 vcftools --gzvcf ./data/27scaffolds.recode.vcf.gz --plink --keep ./indlist/indset$i --chrom-map ./data/chrom-map --out ./data/Fagus$i
done

# Then we need to replace the first columns of the ped file to get the format accepted by GONE (see manual):
# The file ./indlist/F35 looks like:
# 1	IND1	0	0	1	-9
# 1	IND2	0	0	1	-9
# 1	IND3	0	0	1	-9
# 1	IND4	0	0	1	-9
# 1	IND5	0	0	1	-9
# 1	IND6	0	0	1	-9
# 1	IND7	0	0	1	-9
# ...
# 1	IND35	0	0	1	-9

# editing the ped files:
for i in {1..150}; do
 cut -f7- ./data/Fagus$i.ped > ./data/Fagus${i}Geno
done

for i in {1..150}; do
 paste ./indlist/F35 ./data/Fagus${i}Geno > ./data/Fagus$i.ped
 rm ./data/Fagus${i}Geno
done

# moving all files .ped and .map to the directory where running GONE ../results/27contigs/
for i in {1..50}; do
 cp ./data/Fagus$i.ped ./results/27contigs/
done

for i in {1..50}; do
 cp ./data/Fagus$i.map ./results/27contigs/
done

for i in {1..50}; do
 rm ./data/Fagus$i.ped 
 rm ./data/Fagus$i.map 
done

# make sure GONE files are there too (script_GONE.sh; PROGRAMMES directory and INPUT_PARAMETERS_FILE)
# now running GONE by creating the following bash script: "slurmGONE":

for i in {1..150}; do
 bash script_GONE.sh Fagus$i
done
```
Notice that the output file "outfileHWD" will be overwritten, so don't run GONE using a loop if you need that file.
The other output files will have their specific code associated with the input file, so they will be fine.
We can also remove .ped and .map no longer needed (check results have been obtained first!)
```sh
rm  ./results/27contigs/*.map
rm  ./results/27contigs/*.ped
```
Extract *N*<sub>e</sub> estimates in the last generation from all the GONE output files:
```sh
for i in {1..150}; do
 awk 'NR==3' ./results/27contigs/Output_Ne_Fagus$i 
done | cut -f 2 > ./results/27contigs/Ne.txt

# Then let's combine missing data estimates and *N*<sub>e</sub> values with the following script:

# reading the input file
input_file="Ne.txt"
input_data=$(cat "$input_file")

# creating an array with each value repeated 35 times
repeated_data=()
while read -r line; do
    for i in $(seq 1 35); do
        repeated_data+=("$line")
    done
done <<< "$input_data"

# writing the repeated data to the output file
output_file="Ne_repeated.txt"
printf '%s\n' "${repeated_data[@]}" > "$output_file"

# reading the contents of the missingdata.txt and Ne_repeated.txt files into arrays
readarray -t missingdata < missingdata.txt
readarray -t nerepeated < Ne_repeated.txt

# looping over the missingdata array and add the corresponding element from nerepeated to each line
for i in "${!missingdata[@]}"; do
    echo -e "${missingdata[$i]}\t${nerepeated[$i]}"
done > missingvsNe.txt
```
The files missingdata.txt and Ne.txt will be used to build the R script. The file "missingvsNe.txt" will then look like this:
```
ind	miss	Ne
AAABOSDC	0.204339	578.102
AAACOSDC	0.251473	578.102
AAAGOSDC	0.924035	578.102
# ...
```
and the file "missing.txt":
```
INDV	F_MISS	order
AAFLOSDC	0.191767	1
AAFNOSDC	0.198035	2
AABYOSDC	0.199144	3
```
Plotting results in R:
```
library(ggplot2)
library(patchwork)

missing <- read.delim("missingvsNe.txt")
missing$Ne <- as.factor(missing$Ne)

f <- ggplot(missing, aes(x = miss, y = Ne))
g <- f + stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), 
    geom = "linerange", color = "purple", linewidth = 0.3)  + 
	stat_summary(fun.data = median_hilow, mapping = aes(group=Ne), 
    geom = "point", color = "purple", size = 2) +	
	theme_classic() + 
	theme(axis.title.x = element_text(size=24, face="bold"), axis.title.y = element_text(size=24, face="bold"), 
	axis.text.x = element_text(size=24, colour = "black"),axis.text.y = element_text(size=24, colour = "black")) 
# Change tick mark labels
h <- g + scale_y_discrete(breaks = levels(missing$Ne)[c(T, rep(F, 20))],
        labels=c("200","600", "1400", "12000", "260000", "4937000", "4991000", "4994000")) + 
		labs(x ="Proportion of missing data per individual", y = "Ne (Geometric mean)")
# barplot:
missing <- read.delim("missing.txt")
M <- ggplot(missing, aes(x=order, y=F_MISS))+
  geom_bar(stat="identity", fill="purple")+
  theme_classic() + 
	theme(axis.title.x = element_text(size=24, face="bold"), axis.title.y = element_text(size=24, face="bold"), 
	axis.text.x = element_text(size=24, colour = "black"),axis.text.y = element_text(size=24, colour = "black")) +
	labs(x ="Individual", y = "Proportion of missing data")

h / M
```

#### *N*<sub>e</sub> estimation with fewer genomic scaffolds
We can also check what happens when we use fewer genomic scaffolds, for example the 12 scaffolds with the largest number of SNPs. 
First, see script "extract_scaffolds_Fagus.sh" to generate a vcf file with a subset of scaffolds.
```sh
module load vcftools/0.1.16
# generate vcf files with a reduced number of individuals (and also thinned dataset to run in NeEstimator):
vcftools --vcf ./data/12scaffolds.recode.vcf --keep ./indlist/35inds --recode --stdout | gzip > ./data/Fagus12.recode.vcf.gz
vcftools --vcf ./data/12scaffolds.recode.vcf --keep ./indlist/35inds --thin 10000 --recode --out ./data/Fagus12.thinned35
vcftools --vcf ./data/12scaffolds.recode.vcf --keep ./indlist/35inds --thin 5000 --recode --out ./data/Fagus12.thinned2 
# creating chromosome map for GONE analyses
bcftools view -H Fagus12.recode.vcf.gz | cut -f 1 | uniq | awk '{print $0"\t"$0}' > chrom-map12
# we need to replace scaffold name with number, as required by GONE:
nl chrom-map12 > cc
cut -f 1 cc > c1
cut -f 2 cc > c2
paste c2 c1 > chrom-map12
rm cc c1 c2
# generating ped and map files
vcftools --gzvcf ./data/Fagus12.recode.vcf.gz --plink --chrom-map ./data/chrom-map12 --out ./data/Fagus12
cut -f7- ./data/Fagus12.ped > ./data/Fagus12Geno
paste ./indlist/F35 ./data/Fagus12Geno > ./data/Fagus12.ped
rm ./data/Fagus12Geno
# moving all files .ped and .map to the directory where running GONE ../results/12contigs/
cp ./data/Fagus12.ped ./results/12contigs/Fagus12.ped
cp ./data/Fagus12.map ./results/12contigs/Fagus12.map
# running GONE
bash script_GONE.sh Fagus12
```

## References
Gargiulo, R., Decroocq, V., González-Martínez, S.C., Lesur-Kupin, I., Paz-Vinas, I., Schmitt, S., Scotti, I., & Heuertz, M. (2023). Estimation of contemporary effective population size in plant populations: limitations of genomic datasets. In preparation for *Evolutionary Applications*.  

Gargiulo, R., Heuertz, M., & Decroocq, V. (2023). *Prunus armeniaca* SNPs dataset [Data set]. Zenodo. https://doi.org/10.5281/zenodo.8124822  

González-Martínez, S.C., Ridout, K., Pannell, J.R. (2017). Range expansion compromises adaptive evolution in an outcrossing plant. *Current Biology* 27, 2544–2551. https://doi.org/10.1016/j.cub.2017.07.007  

Groppi, A., Liu, S., Cornille, A.,  Decroocq, S ... &  Decroocq, V. (2021). Population genomics of apricots unravels domestication history and adaptive events. *Nature Communications* 12,3956. https://doi.org/10.1038/s41467-021-24283-6  

Lesur-Kupin, I. & Scotti, I. (2023). *Fagus sylvatica* supporting information. https://doi.org/10.57745/FJRYI1  

Schmitt, S., Tysklind, N., Hérault, B., & Heuertz, M. (2021). Topography drives microgeographic adaptations of closely related species in two tropical tree species complexes. *Molecular Ecology* 30, 5080-5093. https://doi.org/10.1111/mec.16116  

Schmitt, S., Tysklind, N., Hérault, B., & Heuertz, M. (2021). Topography drives microgeographic adaptations of closely-related species in two tropical tree species complexes [Data set]. In *Molecular Ecology* (0.1.0). Zenodo. https://doi.org/10.5281/zenodo.4727831    






Shield: [![CC BY 4.0][cc-by-shield]][cc-by]

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg
