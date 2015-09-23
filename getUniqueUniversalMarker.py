__author__ = 'Nina'
import sys


markerFile = sys.argv[1]

def removeMarker(markerFile):
    markerDict = {}
    speciesDict = {}
    speciesOrder = []
    species = ""
    file = open(markerFile)
    for line in file:
        if line.startswith(">"):
            if not species == "":
                speciesDict[species] = chromosomes
            species = line
            speciesOrder.append(species)
            chromosomes = {}
        elif line.startswith("#"):
            chrom = line
        elif not line == "\n":
            order = line[:-2].split(" ")
            chromosomes[chrom] = order
            for elem in order:
                if "-" in elem:
                    elem = elem[1:]
                if elem in markerDict:
                    markerDict[elem].append(species)
                else:
                    markerDict[elem] = [species]
    speciesDict[species] = chromosomes
    file.close()

    newspeciesDict = {}
    for speci in speciesDict:
        chroms = speciesDict[speci]
        newchroms = {}
        for chro in chroms:
            markorder = chroms[chro]
            newmarkorder = []
            for mark in markorder:
                minus = False
                if "-" in mark:
                    mark = mark[1:]
                    minus = True
                if len(markerDict[mark]) == 8:
                    if minus:
                        newmarkorder.append("-"+mark)
                    else:
                        newmarkorder.append(mark)
            newchroms[chro] = newmarkorder
        newspeciesDict[speci] = newchroms

    return newspeciesDict,speciesOrder



def outputMarker(speciesdict,speciesOrder):
    out = open("universal.txt","w")
    for speci in speciesOrder:
        out.write(speci)
        chroms = speciesdict[speci]
        for chro in chroms:
            out.write(chro)
            out.write(" ".join(chroms[chro])+" $\n")
    out.close()



newDicts,order = removeMarker(markerFile)
outputMarker(newDicts,order)