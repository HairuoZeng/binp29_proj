# binp29_proj
## Introduction 
This project is used to provide some simple tools to assemble the mitochondrial genome based on Needleman-Wrench algorithm with greedy method or topological sort by directed graph.

### Summary 
There are four python scripts for different function.

reads_generation.py can simulate the process of reads generation in genome sequencing.

greedy_method.py can overlap the best match reads each round based on Needleman-Wrench algorithm

greedy_method_notNW.py can overlap the best match reads based on number of full match bases between sequences

Directed_graph.py provides two algorithm (Dijkstra_algorithm or Breadth-First Search) for genome assembly based on directed graph

### Usage
reads_generation.py  fasta_file(only including one full length sequence)  simulated_sequencing_depth  simulated_length_of_reads
*output file will be named automaticallly based on the input file.

greedy_method(_notNW).py forwards_direction_reads(generated from reads_generation.py)  reverse_direction_reads(generated from reads_generation.py) output_file_name

Directed_graph.py  forwards_direction_reads(generated from reads_generation.py)  reverse_direction_reads(generated from reads_generation.py)  threshold_of_match_bases  output_file_name
*User can will be asked to choose one algorithm for assembly.BFS is set as default choice.


## Requirements 
build in Python 3.7 
