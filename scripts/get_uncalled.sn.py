# get uncalled from vcf

import gzip

uncalled_dict = {}
with gzip.open(snakemake.input[0]) as vcf:
    for line in vcf:
        if line.startswith("#CHROM"):
            sample_line = line.split()[9:]
            for sample in sample_line:
                uncalled_dict[sample] = []
        elif not line.startswith("#"):
            variant_line = line.split()
            chrom = variant_line[0]
            pos = int(variant_line[1])
            genotypes = variant_line[9:]
            for i, gt in enumerate(genotypes):
                sample = sample_line[i]
                if "./." in gt.split(":")[0]:
                    uncalled_dict[sample].append("pos")

with open(snakemake.output[1], 'w') as uncalled:
    for sample in sample_line:
        uncalled.write(f'{sample}/t{uncalled_dict[sample].join("/t")}/n')
