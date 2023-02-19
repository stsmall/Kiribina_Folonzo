
import glob
import pickle as pk
import numpy as np
from collections import defaultdict


def get_freq(pop, snp_dict):
    """ """
    times = None
    for chrom in ["X", "2", "3"]:
        freq_files = glob.glob(f"{pop}.{chrom}*.freq")
        for freq_file in freq_files:
            with open(freq_file, 'r') as ff:
                line = next(ff)
                header = line.split()[3:-2]
            if times:
                assert times == header
            else:
                times = header
            for line in ff:
                x = line.split()
                snp = x[0]
                freqs = np.array(x[3:-2], dtype=np.float)
                assert len(times) == len(freqs)
                snp_dict[chrom][snp] = freqs

        lin_files = glob.glob(f"{pop}.{chrom}*.lin")
        for lin_file in lin_files:
            with open(lin_file, 'r') as lf:
                line = next(lf)
                header = line.split()[3:-1]
                assert times == header
                for line in lf:
                    x = line.split()
                    snp = x[0]
                    lins = np.array(x[3:-1], dtype=np.float)
                    assert len(times) == len(lins)
                    freqs = snp_dict[chrom][snp]
                    assert len(lins) == len(freqs)
                    snp_dict[chrom][snp] = freqs/lins
    return snp_dict, times


snp_dict_f, times_F = get_freq("F", defaultdict(dict))
with open('F.freq.sel.txt', 'wb') as handle:
    pk.dump(snp_dict_f, handle, protocol=pk.HIGHEST_PROTOCOL)

snp_dict_k, times_K = get_freq("K", defaultdict(dict))
with open('K.freq.sel.txt', 'wb') as handle:
    pk.dump(snp_dict_k, handle, protocol=pk.HIGHEST_PROTOCOL)

assert times_F == times_K
times = list(map(float, times_F))
with open('KF.times.txt', 'wb') as handle:
    pk.dump(times, handle, protocol=pk.HIGHEST_PROTOCOL)
