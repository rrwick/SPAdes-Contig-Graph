# GetSPAdesContigGraph

This tool produces a fastg graph from either the contigs or scaffolds from a [SPAdes](http://bioinf.spbau.ru/spades) assembly.  It uses the sequences from the contigs/scaffolds fasta file and adds connections that it infers from the assembly_graph.fastg file and the paths file.

The sequences in a contigs/scaffolds file made by SPAdes contains longer sequences than the assembly_graph.fastg file, due to repeat resolution and scaffolding.  However, the assembly_graph.fastg file contains connections between sequences, which can be useful for analysing assemblies in [Bandage](http://rrwick.github.io/Bandage/).  This tool combines the two to create a fastg file with both longer sequences and connections.

GetSPAdesContigGraph requires an assembly from SPAdes 3.6.2 or later, as this is the version of SPAdes which introduced the paths files.



# Tools modes

### Length priority

When this script is run in length priority mode (with the `-l` flag), it will produce a graph with the exact same sequences as the input contig file, but with some connections added.  It adds whatever connections are possible between the ends of contigs.

In this mode, the program looks at the paths file to determine which graph segment each contig starts and ends with.  It then uses the graph file to see if there are any connections between these starting/ending segments.  These connections are added to the contigs to create a contig graph.

### Connection priority

When run in connection priority mode (with the `-c` flag), this script will attempt to maintain the connections in the assembly graph, and may split contigs to accomplish this.  Compared to length priority mode, this will produce results more and shorter contigs/scaffolds, but there will be fewer erroneous dead ends in the graph.

In this mode, the program first produces a contig graph like in length priority mode.  It then looks at the assembly graph to see which graph connections are not represented in the contig graph.  It will then split contigs as necessary to allow for these connections to be added to the contig graph.



## Usage

```GetSPAdesContigGraph.py [-h] [-p PATHS_OUT] (-c | -l) graph contigs paths output```

When run in connection priority mode, this script requires BLAST to be installed (specifically it uses `makeblastdb` and `blastn`).

Example usage:

Produce a length-priority contig graph from SPAdes contigs:
`GetSPAdesContigGraph.py -l assembly_graph.fastg contigs.fasta contigs.paths -o contigs.fastg`

Produce a connection-priority contig graph from SPAdes scaffolds and also save the paths file for this graph:
`GetSPAdesContigGraph.py -c -p scaffolds.paths assembly_graph.fastg scaffolds.fasta scaffolds.paths scaffolds.fastg



# Limitations

It is important to note that the contig graphs made by this program may NOT have all of the connection information present in the assembly graph.  This is especially true if the program is run in length priority mode, where many assembly graph connections cannot be added to the contig graph.

To illustrate why, consider the following example (described by SPAdes author Anton Korobeynikov):
* An assembly graph has three segments: A, B and C.
* There are two connections in the graph: A -> B and C -> B.
* SPAdes makes two contigs from the graph: AB and C.

If you then used GetSPAdesContigGraph in length priority mode on this example, there would be no connections.  The A -> B connection is already included in a contig and the C -> B connection has nowhere to attach.  If, however, you used the program in connection priority mode, the AB contig would be split into two contigs (A and B) to allow for C to attach.  In this simple case, the result is no better than the original assembly graph.

Therefore, when using this tool, please keep in mind that the contig graphs it generates may not contain the full complexity of the assembly graph.  In complex or ambiguous cases, it is recommended to use the SPAdes assembly graph.



## License

GNU General Public License, version 3
