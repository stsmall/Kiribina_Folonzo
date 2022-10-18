
import pickle as pk
import sys

with open(sys.argv[1], 'rb') as handle:  # pkb file of freqs
    snp_dt = pk.load(handle)

with open(sys.argv[2], 'rb') as handle:  # pkb file of times
    times = pk.load(handle)

snp = None
sel = 0

freq_file = open("freq_file.txt", 'w')

with open(f"{sys.argv[3]}.newfrq.txt", 'w') as f:
    with open(sys.argv[3], 'r') as infile:  # infile of pop specific data
        line = next(infile)
        f.write(line)
        header = line.split()
        try:
            sel_ix = header.index("sel_time")
            frq_ix = header.index("freq")
            age_ix = header.index("age")
            chrom_ix = header.index("chromosome")
            snp_ix = header.index("snp_pos")
        except ValueError as ve:
            print(ve)
            sys.exit("check headers")
        for line in infile:
            # add zero time
            if snp and float(line.split()[sel_ix]) > sel:
                x[sel_ix] = "0"
                x[sel_ix + 1] = "1"
                rel_frq_0 = snp_dt[chrom][snp][-1]
                try:
                    assert rel_frq_0 == frq_0
                except AssertionError:
                    freq_file.write(f"{snp}\t{rel_frq_0}\t{frq_0}\n")
                x[frq_ix] = str(rel_frq_0)
                x[age_ix] = str(age)
                f.write("{}\n".format("\t".join(x)))
            # split
            x = line.split()
            # gather
            chrom = x[chrom_ix]
            snp = x[snp_ix]
            sel = float(x[sel_ix])
            frq_0 = float(x[frq_ix])
            age = float(x[age_ix])
            # replace
            time_ix = times.index(sel)
            new_frq = snp_dt[chrom][snp][time_ix]
            x[frq_ix] = str(new_frq)
            x[age_ix] = str(round((age - sel), -1))
            # write
            f.write("{}\n".format("\t".join(x)))
