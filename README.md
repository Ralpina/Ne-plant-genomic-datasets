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
setwd("./results/yellow99")
Ne10M <- read.delim("Ne10M.txt")
median(Ne10M[[1]])
quantile(Ne10M[[1]], 0.025, type = 1)
quantile(Ne10M[[1]], 0.975, type = 1)
```

We can now prepare the input file for plotting the results.
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
#### 2. Influence of sample size (resampling only individuals with the smallest proportion of admixture) 

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
The file ```./results/yellow99/NeIndividuals.txt``` includes Ne estimates in the first column and the subset of individuals considered in the second column. We can add the header or not.  

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

# to display zoom on point Ne 
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


#### 3. Influence of population structure
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
Plotting the results; let's first extract last-generation Ne estimates for the Northern gene pool datasets:
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
We also want to get the Ne estimates for the last 25 generations:
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
# Ne during the last generation across 50 datasets:
for i in {1..50}; do
 awk 'NR==3' ./results/all/Output_Ne_all.$i 
done | cut -f 2 > ./results/all/allNe.txt
# last 25 generations across 50 datasets:
for i in {1..50}; do
 awk 'NR>=3 && NR<=27' ./results/all/Output_Ne_all.$i 
done | cut -f 1,2 > ./results/all/all_Ne25gen.txt
```

  
















