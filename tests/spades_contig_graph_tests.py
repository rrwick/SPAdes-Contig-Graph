#!/usr/bin/env python
"""
SPAdes Contig Graph tests

This script runs automated tests on the SPAdes Contig Graph tool.

Author: Ryan Wick
email: rrwick@gmail.com
"""

from __future__ import print_function
import os
import sys
import inspect
import subprocess

sys.path.insert(1, os.path.join(sys.path[0], '..'))

from spades_contig_graph import *

def main():
    print('Running SPAdes Contig Graph tests\n')
    length_priority_test()
    connection_priority_test()

def length_priority_test():
    """
    Run a test on the full execution of spades_contig_graph.py, in length
    priority mode.
    """
    print('Length priority test............... ', end='')
    sys.stdout.flush()

    command = ['../spades_contig_graph.py', 'test_assembly_graph.fastg', 'test_contigs.fasta', 'test_contigs.paths', 'out.fastg', '-l']
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = process.communicate()

    segment_sequences, segment_depths = load_graph_sequences_2('out.fastg')
    links = load_graph_links_2('out.fastg')
    link_tuples = get_link_tuples(links)

    one_test('Segments', len(segment_sequences), 4)
    one_test('Links', len(links), 4)
    one_test('Link tuples', link_tuples, [('2-', '2+')])
    one_test('Node 1 sequence', segment_sequences['1+'], 'CGTTATTCGCGCCCACTCTCCCATTTATCCGCGCAAGCGGATGCGATGCGATTGCCCGCTAAGATATTCTTACCATTCTCGACATGCTGAGCTGAGACGGCGTCGATGCATAGCGGACTTTCGGTCAGTCGCAATTCCTCACGAGACTGGTCCTGTTGACTACAATGGGCCCAACTCAATCACAGCTCGAGCGCCTTGAATAACATACTCATCTCTATACATTCTCGACATGCTGAGCTGAGACGGCGTCGATGCATAGCGGACTTTCGGTCAGTCGCAATTCCTCACGAGACTGGTCCTGTTACAGAGCTGGCGTACGCGTTGAACACTTCACAGATGATAGGGATTCGGGTAAAGAGCGTGTCATTGGGGGCTT')
    one_test('Node 2 sequence', segment_sequences['2+'], 'ATGGCAAGGTACTTCCGGTCTTAATGAATGGCCGGGAAAGGTACGCACGCGGTATGGGGGGGTGAAGGGGCGAATAGACAGGCTNNNNNNNTAAAAATGACAGTGGTTGGTGCTCTAAACTTCATTTGGTTAACTCGTGTATCAGCGCGATAGGCTGTTAGAGGTTTAATATTGTTTAATCCAATTCCCTCATTTAGGACCCTACCAAGTCAACATTGGTATATGAATGCGACCTCGAAGAGGCCGCCAAAGAACAAAGGCTTACTGTGCGCAGAGGAACGCCCATTTAGCGGCTGGCGTTTTGAATCCTTTTAATATTGTTTAATCCAATTCCCTCATTTAGGACCCTACCAAGTCAACATTGGTATATGAATGCGACCTCGAAGAGGCCGCCACCGTTTTAGGGGGGGAAGGTTGAAGATCTCCTCTTCTCATGACTGAACTCGCGAGGGCCGTATTGGGGGCTTCATACATAGAGCAAGGGCGTCGAACGGTCGTGAAAGTCTTAGTACCGCACGTACCAACTTACTGAGGATATTG')

    os.remove('out.fastg')
    print('pass')

def connection_priority_test():
    """
    Run a test on the full execution of spades_contig_graph.py, in connection
    priority mode.
    """
    print('Connection priority test........... ', end='')
    sys.stdout.flush()

    command = ['../spades_contig_graph.py', 'test_assembly_graph.fastg', 'test_contigs.fasta', 'test_contigs.paths', 'out.fastg', '-c']
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = process.communicate()

    segment_sequences, segment_depths = load_graph_sequences_2('out.fastg')
    links = load_graph_links_2('out.fastg')
    link_tuples = get_link_tuples(links)

    one_test('Segments', len(segment_sequences), 6)
    one_test('Links', len(links), 6)
    one_test('Link tuples', link_tuples, [('1+', '3+'), ('1-', '1+'), ('2+', '3+'), ('3-', '1-'), ('3-', '2-')])
    one_test('Node 1 sequence', segment_sequences['1+'], 'ATGGCAAGGTACTTCCGGTCTTAATGAATGGCCGGGAAAGGTACGCACGCGGTATGGGGGGGTGAAGGGGCGAATAGACAGGCTNNNNNNNTAAAAATGACAGTGGTTGGTGCTCTAAACTTCATTTGGTTAACTCGTGTATCAGCGCGATAGGCTGTTAGAGGTTTAATATTGTTTAATCCAATTCCCTCATTTAGGACCCTACCAAGTCAACATTGGTATATGAATGCGACCTCGAAGAGGCCGCCAAAGAACAAAGGCTTACTGTGCGCAGAGGAACGCCCATTTAGCGGCTGGCGTTTTGAATCCTTTTAATATTGTTTAATCCAATTCCCTCATTTAGGACCCTACCAAGTCAACATTGGTATATGAATGCGACCTCGAAGAGGCCGCCACCGTTTTAGGGGGGGAAGGTTGAAGATCTCCTCTTCTCATGACTGAACTCGCGAGGGCCGTATTGGGGGCTT')
    one_test('Node 2 sequence', segment_sequences['2+'], 'CGTTATTCGCGCCCACTCTCCCATTTATCCGCGCAAGCGGATGCGATGCGATTGCCCGCTAAGATATTCTTACCATTCTCGACATGCTGAGCTGAGACGGCGTCGATGCATAGCGGACTTTCGGTCAGTCGCAATTCCTCACGAGACTGGTCCTGTTGACTACAATGGGCCCAACTCAATCACAGCTCGAGCGCCTTGAATAACATACTCATCTCTATACATTCTCGACATGCTGAGCTGAGACGGCGTCGATGCATAGCGGACTTTCGGTCAGTCGCAATTCCTCACGAGACTGGTCCTGTTACAGAGCTGGCGTACGCGTTGAACACTTCACAGATGATAGGGATTCGGGTAAAGAGCGTGTCATTGGGGGCTT')
    one_test('Node 3 sequence', segment_sequences['3+'], 'ATTGGGGGCTTCATACATAGAGCAAGGGCGTCGAACGGTCGTGAAAGTCTTAGTACCGCACGTACCAACTTACTGAGGATATTG')

    os.remove('out.fastg')
    print('pass')

def one_test(test_name, actual, expected):
    """
    Carry out a single test, comparing actual value to expected value.  If they
    are different, print an error message and quit the program.
    """
    if expected != actual:
        print('fail')
        print('Failed test:', test_name)
        print('  Expected:', expected)
        print('  Actual:  ', actual)
        quit()

def get_link_tuples(link_dict):
    """
    Use the link dictionary to produce a sorted list of link tuples in the from
    of (start_segment, end_segment).
    """
    link_tuples_set = set()
    for start, ends in link_dict.iteritems():
        for end in ends:
            link_tuples_set.add((start, end))
    return sorted(list(link_tuples_set))

if __name__ == '__main__':
    main()
