__author__ = 'Nina'

import sys

marker = sys.argv[1]

file = open(marker,"r")
species = ""
speciesHash = {}
for line in file:
    if line.startswith(">"):
        species = line.split("\t")[0][1:]
    elif line.startswith("#"):
        chrom = line.rstrip("\n")[2:]
    elif not line == "\n":
        mark = line.rstrip("\n")[:-2].split(" ")
       # print mark
        if species in speciesHash:
            speciesHash[species] += mark
        else:
            speciesHash[species] = mark


refList = speciesHash["hg18"]
for elem in refList:
    for species in speciesHash:
        negative = "-"+elem
        if not (elem in speciesHash[species] or negative in speciesHash[species]):
            print species
            print elem
