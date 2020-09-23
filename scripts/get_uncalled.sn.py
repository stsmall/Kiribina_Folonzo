import gzip
import sys
from collections import defaultdict
uncalled_dict = defaultdict(list)

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
            pos = variant_line[1]
            genotypes = variant_line[9:]
            for i, gt in enumerate(genotypes):
                sample = sample_line[i]
                if "./." in gt.split(":")[0] or ".|." in gt.split(":")[0]:
                    uncalled_dict[sample].append(pos)

with gzip.open(f"{sys.argv[2]}.gz", 'wt') as uncalled:
    for sample in sample_line:
        tab_line = "\t".join(uncalled_dict[sample])
        uncalled.write(f'{sample}\t{chrom}\t{tab_line}\n')
