# PhySca

PhySca samples solutions to the SCJ small parsimony problem for weighted gene adjacencies.


### Requirements

PhySca is composed of a set of Python scripts. It requires the following to be available:

* python (2.7)
* python packages: networkx and ete2
* software DeClone (https://github.com/yannponty/DeClone)


### Usage
```
Main.py [-h] [-tree <treefile>] [-alpha <a>] (-m <markerfile> | -a <adjacencies>)
               (-x <x> | -w <weights>) [-gx <weights>] [-kT <kT>] [-s <Z>]
               
```


  *-h*, --help   show this help message and exit
  
  *-tree \<treefile>*            tree file in newick format, if DeClone is used, the tree has to be specified in Extended Newick 
        (NHX) format
  
  *-alpha \<a>*          alpha parameter in objective function, [0,1]
  
  *-m , --marker \<markerfile>*       marker order of extant genomes
                        
  *-a , --adjacencies \<adjacencies>*
                        adjacencies in extant genomes
                        
  *-x <x>*                 Assign potential adjacencies by weight threshold, [0,1]
                        
  *-w , --weight \<weightfile>*
                        weights for adjacencies at specific internal nodes,
                        adjacencies at other nodes are assigned a weight of 0
                        
  *-gx \<weightfile>*                weights for adjacencies at specific internal nodes,
                        adjacencies at other nodes have weights computed by
                        Declone, therefore parameter -x must be given!
                        
  *-kT \<kT>*                deClone constant
  
  *-s , --sampling <Z>*
                        sample Z solutions for given set of parameters

#### Input files
##### Marker orders for extant genomes
To account for orientation, marker are signed.

Example:
```
>genome 1   #marker     #chromosomes
# chr1
1 2 3 4 5 6 7 8 $
# chr2
10 -11 -14 15 16 $
```

##### Adjacencies in extant genomes

Markers were doubled to account for orientation: marker X consists of two extremities, with respective IDs 2\*X (for the head of the marker) and 2\*X-1 (for the tail of the marker). Then we have the following line for each adjacency A=(extremity1, extremity2), followed by a list of extant genomes that contain this adjacency.

Example:
```
A0:extremity1 extremity2  #genome1:start-stop - #genome2:start-stop  #genome3:start-stop  
A1: ...
```

##### Given adjacency weights for specific internal nodes of the tree

Example:
```
>extrem1   extrem2    gapID   weight in [0,1]     genome1:start-stop       genome2:start-stop   
```


#### Output files

* *ancestral_assigned_adjacencies_with_weight* --- contains all potential adjacencies at each internal node after filtering, together with its weight value

* *conflicts* ---
  contains all pairwise conflicting adjacencies in the set of potential adjacencies at each internal node of the tree

* *reconstructed_adjacencies* --- contains all reconstructed adjacencies at each internal node of the tree

* *doubled_scaffolds, undoubled_scaffolds* --- contains the set of scaffolds/CARs at each internal node according to the marker order file format described above, either as sequences of marker extremities or signed marker




