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
