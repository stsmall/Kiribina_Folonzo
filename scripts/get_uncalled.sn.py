import gzip
import sys

uncalled_dict = {}
with gzip.open(sys.argv[1], 'rb') as vcf:
    for line in vcf:
        line = line.decode()
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

with open(sys.argv[2], 'wt') as uncalled:
    for sample in sample_line:
        tab_line = "/t".join(uncalled_dict[sample])
        uncalled.write(f'{sample}/t{tab_line}/n')
