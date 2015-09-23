from unittest import TestCase
import getAdjacencies
import globalAdjacencyGraph
from ete2 import Tree
import networkx as nx

__author__ = 'Nina'


class TestGraph(TestCase):
  chrom = {}
  chrom["one"] = ["3","4","5","6","7","8"]
  chrom["two"] = ["1","2"]
  species = {}
  species["A"] = chrom
  chrom = {}
  chrom["one"] = ["-4","5","-3","6","7","8"]
  chrom["two"] = ["1","2"]
  species["B"] = chrom
  chrom = {}
  chrom["one"] = ["3","4","5","6","7","8"]
  chrom["two"] = ["1","2"]
  species["E"] = chrom
  chrom = {}
  chrom["one"] = ["4","5","-3","6","7","8"]
  chrom["two"] = ["1","2"]
  species["D"] = chrom

  adj = getAdjacencies.findAdjacencies(species)
  tree = Tree("(A:1,(B:1,(E:1,D:1)Int_1:0.5[&&NHX:ancient=1])Int_2:0.5[&&NHX:ancient=0])Root;", format=1)
  paths = getAdjacencies.findTreePaths(tree)
  internal,adjacenciesAncestral = getAdjacencies.assignAncestralAdjacencies(paths,adj,tree)
  graphs = globalAdjacencyGraph.createGraph(adj,adjacenciesAncestral)


  def test_createGraph(self):
    TestCase.assertEqual(self,len(self.graphs),5)
    TestCase.assertEqual(self,nx.number_of_nodes(self.graphs[0]),5)
    TestCase.assertEqual(self,nx.number_of_edges(self.graphs[0]),4)


  def test_analyseConnectedComponents(self):
    globalAdjacencyGraph.analyseConnectedComponents(self.graphs)
    #passt



