__author__ = 'Nina'

import unittest
import getAdjacencies
from ete2 import Tree


class MyTestCase(unittest.TestCase):
    chrom = {}
    chrom["one"] = ["3","4","5"]
    chrom["two"] = ["1","2"]
    species = {}
    species["A"] = chrom
    chrom = {}
    chrom["one"] = ["-4","5","-3"]
    chrom["two"] = ["1","2"]
    species["B"] = chrom
    tree = Tree("(A:1,(B:1,(E:1,D:1)Int_1:0.5[&&NHX:ancient=1])Int_2:0.5[&&NHX:ancient=0])Root;", format=1)
    paths = getAdjacencies.findTreePaths(tree)
    adj = getAdjacencies.findAdjacencies(species)
    leaveA = tree.get_leaves_by_name("A")[0]
    leaveB = tree.get_leaves_by_name("B")[0]
    leaveD = tree.get_leaves_by_name("D")[0]
    leaveE = tree.get_leaves_by_name("E")[0]
    int1 = leaveE.get_common_ancestor(leaveD)
    int2 = leaveE.get_common_ancestor(leaveB)
    int3 = leaveE.get_common_ancestor(leaveA)


    def test_doubleMarker(self):
        marker = "3"
        doubled = getAdjacencies.doubleMarker(marker)
        unittest.TestCase.assertEqual(self,doubled,(5,6))
        marker = "-3"
        doubled = getAdjacencies.doubleMarker(marker)
        unittest.TestCase.assertEqual(self,doubled,(6,5))


    def test_findAdjacencies(self):
        unittest.TestCase.assertIn(self,(6,7),self.adj)
        unittest.TestCase.assertIn(self,(8,9),self.adj)
        unittest.TestCase.assertIn(self,(2,3),self.adj)
        unittest.TestCase.assertIn(self,(7,9),self.adj)
        unittest.TestCase.assertIn(self,(10,6),self.adj)

    def test_findTreePath(self):
        unittest.TestCase.assertIn(self,(self.leaveA,self.leaveB),self.paths)
        pathAB = [self.int3, self.int2]
        pAB = self.paths[(self.leaveA,self.leaveB)]
        unittest.TestCase.assertEqual(self,pathAB,pAB)
        pathDE = [self.int1]
        pDE = self.paths[(self.leaveE,self.leaveD)]
        unittest.TestCase.assertEqual(self,pathDE,pDE)

    def test_assignAncestralAdjacencies(self):
        internal,adjacenciesAncestral = getAdjacencies.assignAncestralAdjacencies(self.paths,self.adj,self.tree)
        internal1 = internal[self.int1]
        unittest.TestCase.assertEqual(self,internal1,set([]))
        internal2 = internal[self.int2]
        internal3 = internal[self.int3]
        unittest.TestCase.assertEqual(self,internal2,internal3)
        shared = set([(2,3)])
        unittest.TestCase.assertEqual(self,internal2,shared)
        shared = set([(7,9)])
        unittest.TestCase.assertNotEqual(self,internal2,shared)
        anc = {(2,3):set([self.int3,self.int2])}
        unittest.TestCase.assertEqual(self,anc,adjacenciesAncestral)

    def test_deCloneProbabilities(self):
        extant = {(1,2):[("bosTau3","chr18"),("hg18","chr1"), ("mm9","chr2")]}
        getAdjacencies.deCloneProbabilities(extant,"","")


if __name__ == '__main__':
    unittest.main()
