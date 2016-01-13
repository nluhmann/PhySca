# PhySca

PhySca samples solutions to the SCJ small parsimony problem for weighted gene adjacencies.


## Requirements

PhySca is composed of a set of Python scripts. It requires the following to be available:

* python (2.7)
* python packages: networkx and ete2
* software DeClone (https://github.com/yannponty/DeClone)


## Usage

Main.py [-h] [-tree \<treefile>] \[-alpha \<a>] (-m \<markerfile> | -a \<adjacencies>)
               (-x <x> | -w \<weights>) [-gx \<weights>] [-kT \<kT>] [-s \<Z>]
               



  -h, --help              show this help message and exit
  
  -tree \<treefile>            tree file in newick format, if DeClone is used, the tree has to be specified in Extended Newick 
        (NHX) format
  
  -alpha \<a>          alpha parameter in objective function, [0,1]
  
  -m , --marker \<markerfile>  
                        marker order of extant genomes
                        
  -a , --adjacencies \<adjacencies>
                        adjacencies in extant genomes
                        
  -x <x>                  Assign potential adjacencies by weight threshold, [0,1]
                        
  -w , --weight \<weightfile>
                        weights for adjacencies at specific internal nodes,
                        adjacencies at other nodes have a weight of 0
                        
  -gx \<weightfile>                weights for adjacencies at specific internal nodes,
                        adjacencies at other nodes have weights computed by
                        Declone, therefore parameter -x must be given!
                        
  -kT \<kT>                deClone constant
  
  -s , --sampling <Z>
                        sample Z solutions for given set of parameters

### Specific input file formats


## Output
