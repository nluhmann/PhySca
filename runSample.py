from ete2 import Tree
import getAdjacencies
import globalAdjacencyGraph
import SR
import scaffolding
import argparse
import time

from calculate_SCJ import calculate_SCJ
import multiprocessing

def runSample(lock, ccs, tree, extantAdjacencies, adjacencyProbs, alpha, i,
              extantAdjacencies_species_adj, reconstructedMarkerCount, allSampleReconstructionStatistic,dict_SCJ):
    t1 = time.time()

    jointLabels, first = SR.enumJointLabelings(ccs)
    validLabels, validAtNode = SR.validLabels(jointLabels, first)

    print "###################### " + str(i) + " ######################"
    lock.acquire()
    topDown = SR.sampleLabelings(tree, ccs, validAtNode, extantAdjacencies, adjacencyProbs, alpha)
    lock.release()
    reconstructedAdj = SR.reconstructedAdjacencies(topDown)
    SR.outputReconstructedAdjacencies(reconstructedAdj, "reconstructed_adjacencies_" + str(i))
    lock.acquire()
    for node in reconstructedAdj:
        print node
        print "Number of reconstructed adjacencies: " + str(len(reconstructedAdj[node]))
        # count for each adjaency on each internal node, how often this adjacencies over all samples occurs there
        for adjacency in reconstructedAdj[node]:
            #if node in allSampleReconstructionStatistic:
            #    # print allSampleReconstructionStatistic[node]
            #    if adjacency in allSampleReconstructionStatistic[node]:
            #        allSampleReconstructionStatistic[node][adjacency] += 1
            #    else:
            #        dict_adj = {adjacency: 1}
            #        allSampleReconstructionStatistic[node].update(dict_adj)
            #        # print allSampleReconstructionStatistic[node][adjacency]
            #else:
            #    allSampleReconstructionStatistic.update({node: {adjacency: 1}})
            if (node,adjacency) in allSampleReconstructionStatistic:
                allSampleReconstructionStatistic[(node,adjacency)] += 1
            else:
                allSampleReconstructionStatistic.update({(node,adjacency):1})
    lock.release()
    scaffolds = scaffolding.scaffoldAdjacencies(reconstructedAdj)
    undoubled = scaffolding.undoubleScaffolds(scaffolds)
    scaffolding.outputUndoubledScaffolds(undoubled, "undoubled_scaffolds_" + str(i))
    scaffolding.outputScaffolds(scaffolds, "doubled_scaffolds_" + str(i))
    scaffolding.sanityCheckScaffolding(undoubled)

    for node in undoubled:
        print node
        markerCounter = 0
        for scaffold in undoubled[node]:
            first = scaffold[0]
            last = scaffold[-1]
            if not first == last:
                markerCounter = markerCounter + len(scaffold)
            else:
                markerCounter = markerCounter + len(scaffold) - 1
        print node + " number of reconstructed undoubled marker in scaffolds: " + str(markerCounter)
        # number of reconstructed markerIds
        #reconstructedMarkerCount = len(reconstructedMarker)
        # singleton scaffolds number / number of not reconstructed marker
        notReconstructedMarkerCount = reconstructedMarkerCount - markerCounter
        # Vergleich/Differenz markerCounter-rekonstruierte Markeranzahl=singleton scaffolds Anzahl
        # number of all scaffolds
        allScaffoldCount = markerCounter + notReconstructedMarkerCount
        # Summe markerCounter +singleton scaffolds Anzahl= Gesamt Scaffold anzahl
        print node + " number of singleton scaffolds (not reconstructed marker): " + str(
            notReconstructedMarkerCount)
        print node + " number of scaffolds: " + str(allScaffoldCount)
    print time.time() - t1, "seconds process time"
    lock.acquire()
    scj = calculate_SCJ(tree, reconstructedAdj, extantAdjacencies_species_adj)
    dict_SCJ.update({'Sample_' + str(i): scj})
    lock.release()
    #return reconstructedAdj, scj