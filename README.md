# GetSPAdesContigGraph

This tool produces a fastg graph from either the contigs or scaffolds from a [SPAdes](http://bioinf.spbau.ru/spades) assembly.  It uses the sequences from the contigs/scaffolds fasta file and adds connections that it infers from the assembly_graph.fastg file and the paths file.

The sequences in a contigs/scaffolds file made by SPAdes contains longer sequences than the assembly_graph.fastg file, due to repeat resolution and scaffolding.  However, the assembly_graph.fastg file contains connections between sequences, which can be useful for analysing assemblies in [Bandage](http://rrwick.github.io/Bandage/).  This tool aims to combine the best of each to create a fastg file with both longer sequences and connections.

GetSPAdesContigGraph requires an assembly from SPAdes 3.6.2 or later.


## Usage

```GetSPAdesContigGraph.py [-h] [-o OUTPUT] graph contigs paths```

Example usage:
`GetSPAdesContigGraph.py assembly_graph.fastg contigs.fasta contigs.paths -o contigs.fastg`

If an output file isn't specified, the graph will be outputted to stdout.

GetSPAdesContigGraph has no dependencies other than Python.



## License

GNU General Public License, version 3
