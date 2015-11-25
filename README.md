# GetSPAdesContigGraph

This tool produces a fastg graph from either the contigs or scaffolds from a [SPAdes](http://bioinf.spbau.ru/spades) assembly.  It uses the sequences from the contigs/scaffolds fasta file and adds connections that it infers from the assembly_graph.fastg file and the paths file.

The sequences in a contigs/scaffolds file made by SPAdes contains longer sequences than the assembly_graph.fastg file, due to repeat resolution and scaffolding.  However, the assembly_graph.fastg file contains connections between sequences, which can be useful for analysing assemblies in [Bandage](http://rrwick.github.io/Bandage/).  This tool combines the two to create a fastg file with both longer sequences and connections.


# How it works

GetSPAdesContigGraph looks at the paths file to determine which graph segment each contigs starts and ends with.  It then uses the graph file to see if there are any connections between these starting/ending segments.  These connections are added to the contigs to create a contig graph.

GetSPAdesContigGraph requires an assembly from SPAdes 3.6.2 or later, as this is the version of SPAdes which introduced the paths files.


# Limitations

It is important to note that the contig graph which results from this program will NOT necessary have all of the connections present in the full assembly graph.

To illustrate why, consider the following example (described by SPAdes author Anton Korobeynikov):
* An assembly graph has three segments: A, B and C.
* There are two connections in the graph: A -> B and C -> B.
* When producing contigs from the graph, SPAdes makes two: AB and C.

If you then used GetSPAdesContigGraph on this example, there would be no connections.  The A -> B connection is already included in a contig and the C -> B connection has nowhere to attach.

Therefore, consider graphs made by GetSPAdesContigGraph as containing SOME connections but not all.  You must view the SPAdes' assembly_graph.fastg file to see all assembly connections.


## Usage

```GetSPAdesContigGraph.py [-h] [-o OUTPUT] graph contigs paths```

Example usage:

`GetSPAdesContigGraph.py assembly_graph.fastg contigs.fasta contigs.paths -o contigs.fastg`

If an output file isn't specified, the graph will be outputted to stdout.

GetSPAdesContigGraph has no dependencies other than Python.


## License

GNU General Public License, version 3
