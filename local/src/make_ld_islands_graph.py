#!/usr/bin/env python
from __future__ import with_statement

from sys import stdin, stderr
from optparse import OptionParser
import networkx as nx

def main():
	parser = OptionParser(usage="make_ld_islands_graph -c cutoff -r chr < STDIN")
	
	parser.add_option('-c', '--cutoff', type=float, dest='cutoff', default=0.8, help='R2 cutoff used to define islands [default: %default]')
	parser.add_option('-r', '--chr', type=str, dest='chr', help='chr')

	options, args = parser.parse_args()
	
	if len(args) != 0:
		exit('Unexpected argument number.')
	if options.chr == None:
		exit("-r is compulsory")

	# undirected graph (with self loops but who cares?).
	ldGraph = nx.Graph()
	for line in stdin:
		tokens = line.split('\t')
		node1 = tokens[2] + "_" + tokens[0]
		node2 = tokens[3] + "_" + tokens[1]
		r2 = float(tokens[4])
		if r2 >= options.cutoff:
			ldGraph.add_edge(node1, node2)
		else:
			ldGraph.add_node(node1)
			ldGraph.add_node(node2)
	
	islands = nx.connected_components(ldGraph)
	for cc in islands:
		# cc is a set
		#print "\t".join(cc)
		start = 3079843747 # yay
		end = -1
		snps = []
		for snp in cc:
			explode = snp.split("_")
			snps.append(explode[0])
			coord = int(explode[1])
			if coord < start:
				start = coord
			if coord > end: #not elif for islands made of a single snp
				end = coord
		start = int(start) - 1 
		print(options.chr + "\t" + str(start) + "\t" + str(end) + "\t" + "_".join(snps)) # convert coords from 1 to 0 based

	
if __name__ == '__main__':
	main()

