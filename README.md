# PhySca

PhySca samples solutions to the Single-Cut-and-Join (SCJ) small parsimony problem for weighted gene adjacencies.


### Requirements

PhySca is composed of a set of Python scripts. It requires the following to be available:

* python (2.7)
* python packages: networkx and ete2
* software DeClone (https://github.com/yannponty/DeClone)


### Quick Start

If you want to quickly make a whole run through of the program including preprocessing you can use the script 'parallel_main.sh'.

### parallel_main.sh


```
parallel_main.sh [-h](-nhx <nhx_tree> | -nf <newick_tree>) 
		 (-a <adjacencies>/ -m <markers>) 
		 [-out <output>] 
		 [-pN <number of processes>] 
		 [-alpha <alpha>] [-s <Z>] [-x <x>] [-kT <kT>]
		 [-phySca <phySca>]
		

```

  *-h*,--help		 show this help message and exit

  *-nhx \<nhx_tree>*		tree file in newick (NHX) format

  *-nf \<newick_tree>*		path to the file with NEWICK-tree

  *-a \<adjacencies>*	adjacencies in extant genomes

  *-m \<markers>*			marker order of extant genomes

  *-output \<output>*		path to directory for preprocessing output, default: ./testlauf

  *-pN, \<processnumber>*		number of processes used for sampling. Max: [number of cpu], default=1

  *-alpha \<alpha>*		alpha parameter in objective function, [0,1]
	
  *-x \<x>*                  assign potential adjacencies by weight threshold, [0,1]
	
  *-kT \<kT>*                deClone constant

  *-s \<Z>*		sample Z solutions for given set of parameters

  *-phySca \<phySca>*	path to directory with main program (PhySca), default: ./

### Usage

### 1) Preprocessing with weightingWithDeClone


```
weightingWithDeClone.py [-h] (-nhx <nhx_tree> | -nf <newick_tree>) [-i]
                               [-sm <minimum>] (-a <adjacencies> | -m <markers>)
                               [-kT <kT>] [-jP] [-out <output>]

```

  *-h*, --help            show this help message and exit

  *-nhx , --nhx_Tree \<nhx_tree>*		tree file in newick (NHX) format

  *-nf , --Newick \<newick_tree>*		path to the file with NEWICK-tree

  *-i , --ignore_weights*	if set, for either ignore edge length/weights, when parsing Newick Tree into nhx Tree

  *-sm , --set_minimum \<minimum>*		minimal value for every edge length, when parsing Newick Tree into nhx Tree

  *-a , --adjacencies \<adjacencies>*	adjacencies in extant genomes

  *-m, --markers \<markers>*			marker order of extant genomes

  *-kT \<kT>*                deClone constant

  *-jP , --just_Parse*     if set, just parse the Newick file into nhx file and not run DeClone after it.

  *-out , --output \<output>*		output directory, current directory as default

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

* *nhx_tree* --- contains given tree in nhx format

* *single_leaf_adjacencies* --- contains all adjacencies, which only ocur at one leaf/ in one extant genome

* *weighted_extant_adjacencies* --- contains all adjacencies in extant genomes

* *weighted_internal_adjacencies* --- contains all adjacencies at internal nodes


```
Structure of weighted_extant_adjacencies:
>node	adjacency
```

```
Structure of weighted_internal_adjacencies:
>node	adjacency	weight
```

### 2) Main

```
Main.py [-h] [-tree <treefile>] [-alpha <alpha>] [-extant <extant>]
               [-internal <internal>] [-x <x>] [-s <Z>] [-pN <processnumber>] [-out <output>] [-skip <int>]

```

  *-h*, --help			show this help message and exit

  *-tree \<treefile>*		tree file in newick or nhx format

  *-alpha \<alpha>*		alpha parameter in objective function, [0,1]

  *-extant \<extant>*		file with precomputed weighted adjacencies for
                        external nodes

  *-internal \<internal>*	    file with precomputed weighted adjacencies for
                        internal nodes

  *-x \<x>*                  (optional) assign potential adjacencies by weight threshold,
                        [0,1]

  *-s , --sampling \<Z>*		(optional) sample Z solutions for given set of parameters

  *-pN, --processNumber \<processnumber>*		(optional) number of processes used for sampling. Max: [number of cpu], default=1

  *-out , --output \<output>*		(optional) output directory, current directory as default

  *-skip , --skip \<int>*       (optional) skip all connected components that consist of more edges than <int>. Note that these adjacencies will not be reconstructed.

#### Input files (precomputed by weightingWithDeClone)

* *weighted_extant_adjacencies* --- contains all adjacencies in extant genomes

* *weighted_internal_adjacencies* --- contains all adjacencies at internal nodes
 

```
Structure of weighted_extant_adjacencies:
>node	adjacency
```

```
Structure of weighted_internal_adjacencies:
>node	adjacency	weight
```



#### Output files

* *conflicts* ---
  contains all pairwise conflicting adjacencies in the set of potential adjacencies at each internal node of the tree

* *reconstructed_adjacencies* --- contains all reconstructed adjacencies at each internal node of the tree

* *doubled_scaffolds, undoubled_scaffolds* --- contains the set of scaffolds/CARs at each internal node according to the marker order file format described above, either as sequences of marker extremities or signed marker

* *SCJ_distances* --- contains the Single Cut or Join distance for each sampled solution

* *statistic_allSampled_Reconstructed_Adjacencies* --- contains for each internal node how often each adjacency occured at this node over all samples 


### Reference
Nina Luhmann, Manuel Lafond, Annelyse Thevenin, Aida Ouangraoua, Roland Wittler, Cedric Chauve, "The SCJ Small Parsimony Problem for Weighted Gene Adjacencies", IEEE/ACM Transactions on Computational Biology and Bioinformatics, [doi:10.1109/TCBB.2017.2661761 ](http://ieeexplore.ieee.org/document/7837680/)
