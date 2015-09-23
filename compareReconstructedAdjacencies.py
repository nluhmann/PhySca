__author__ = 'Nina'
import sys

file1 = sys.argv[1]
file2 = sys.argv[2]


#read list of reconstructed adjacencies
def readList(file):
    adjacencyHash = {}
    f = open(file,"r")
    for line in f:
        if line.startswith(">"):
            species = line[1:].split(" ")[0]
            adjacencyHash[species] = []
        elif not line == " " and not line == "\n":
            adj = line.rstrip("\n")
            adjacencyHash[species].append(adj)
    f.close()
    return adjacencyHash

adjHash1 = readList(file1)
adjHash2 = readList(file2)
print ">"+str(file1)+" "+str(file2)
#compare hashes!
for spec in adjHash1:
    adjList1 = adjHash1[spec]
    adjList2 = adjHash2[spec]
    elems1 = []

    for elem in adjList1:
        if not elem in adjList2:
            elems1.append(elem)
    if not elems1 == []:
        print spec
        print "adjacencies only in "+str(file1)
        print "\n".join(elems1)
    elems2 = []

    for elem in adjList2:
        if not elem in adjList1:
            elems2.append(elem)
    if not elems2 == []:
        print spec
        print "adjacencies only in "+str(file2)
        print "\n".join(elems2)