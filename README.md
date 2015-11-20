# GetSPAdesContigGraph

This tool produces a fastg graph from either the contigs or scaffolds from a [SPAdes](http://bioinf.spbau.ru/spades) assembly.  It uses the sequences from the contigs/scaffolds fasta file and adds connections that it infers from the assembly_graph.fastg file and the paths file.

The resulting graph has advantages over both the SPAdes assembly_graph.fastg file and the contigs/scaffolds fasta file:
* Sequences have been processed by SPAdes' repeat resolution, like the contigs/scaffolds fasta files but unlike the assembly_graph.fastg file.
* Connections are present between sequences, like the assembly_graph.fastg file but unlike the contigs/scaffolds files.

GetSPAdesContigGraph requires an assembly from SPAdes 3.6.2 or later.


## Usage

```GetSPAdesContigGraph.py [-h] [-o OUTPUT] graph contigs paths```

Example usage:
`GetSPAdesContigGraph.py assembly_graph.fastg contigs.fasta contigs.paths -o contigs.fastg`

If an output file isn't specified, the graph will be outputted to stdout.

GetSPAdesContigGraph has no dependencies other than Python.



## License

GNU General Public License, version 3
