#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
from Bio import Phylo as phy
import csv

trees = phy.parse(sys.argv[1], 'newick')
sp2 = sys.argv[2]
sp3 = sys.argv[3]
sp4 = sys.argv[4]
DistanceABAB = []
DistanceBCBC = []
DistanceACAB = []
DistanceACBC = []

for tree in trees:  # for each genealogy
    if tree.distance(sp3,sp4) < tree.distance(sp2,sp3) and tree.distance(sp3,sp4) < tree.distance(sp2,sp4): #if the topology is AB
        DistanceABABi = tree.distance(sp3,sp4)  # AB given AB
        DistanceACABi = tree.distance(sp2,sp4)  # AC given AB
        DistanceABAB.append(DistanceABABi)
        DistanceACAB.append(DistanceACABi)
    if tree.distance(sp2,sp3) < tree.distance(sp3,sp4) and tree.distance(sp2,sp3) < tree.distance(sp2,sp4): #if the topology is BC
        DistanceBCBCi = tree.distance(sp2,sp3)  # BC given BC
        DistanceACBCi = tree.distance(sp2,sp4)  # AC given BC
        DistanceBCBC.append(DistanceBCBCi)
        DistanceACBC.append(DistanceACBCi)

D1 = (sum(DistanceABAB)/len(DistanceABAB))-(sum(DistanceBCBC)/len(DistanceBCBC))  # Estimate D1
D2 = (sum(DistanceACAB)/len(DistanceACAB))-(sum(DistanceACBC)/len(DistanceACBC))  # Estimate D2

# Write to csv
with open(sys.argv[5], 'a') as statsfile:
    wr = csv.writer(statsfile)
    wr.writerow([D1, D2])
