#!/usr/bin/env python

print "\nReplace tip names (accession or names) in newick tree with a new name." 
print "Contributor: Laura Wegener Parfrey"
print "November 2016"
print "\nrun: python replace_in_tree.py -h for help.\n"

# author: Laura Wegener Parfrey
# date: November 2016
# goal: Replace tip names (accession or names) in newick tree with a new name 
# run: python replace_in_tree.py -h for help.

from sys import argv
import re
import argparse

parser = argparse.ArgumentParser(
    description='Replace tip name in tree')
parser.add_argument(
    '-t',
    '--input_tree',
    help='Path to input newick tree',
    required=True)
parser.add_argument(
    '-o', 
    '--output_tree',
    help='Path to output newick tree with renamed tips',
    required=True)
parser.add_argument(
    '-i',
    '--tip_metadata',
    help='Path to input tab delimited file current tip names and updated names.',
    required=True)
parser.add_argument(
    '--initial_column',
    help='number of the column with initial tip names. Default = 1 (first column)',
    required=False)
parser.add_argument(
    '--rename_column',
    help='Number of the colume with new names. Default = 2.',
    required=False)



args = parser.parse_args()

input_tree = args.input_tree
output_tree = args.output_tree
tip_metadata = args.tip_metadata
initial_column = args.initial_column
rename_column = args.rename_column

if initial_column is not None:
	initial = int(initial_column) - 1
else: 
	initial = 0
	
if rename_column is not None:
	rename = int(rename_column) - 1
else:
	rename = 1
	
# tree 
tree = open(input_tree, "U").read()

# open tab delimited file to with current (column 1) and replacement (column 2) tree tip names. 
for line in open(tip_metadata, "U"):
	entry = line.strip().split('\t')
	tree = tree.replace(entry[initial],entry[rename].strip()) 

out_tree = open(output_tree, "w")

out_tree.write(tree)
out_tree.close()

