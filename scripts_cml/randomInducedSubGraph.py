# __author__ = 'Simon Dirmeier'
# __email__  = 'simon.dirmeier@bsse.ethz.ch'
# __date__   = 04/04/16

from __future__ import print_function, absolute_import
import sys
import argparse


class Graph:
    def __init__(self, n):
        self.n = int(n)
        self.genes = set()
        self.edges = set()

    def add_node(self, nodes):
        for n in nodes:
            if len(self.genes) < self.n:
                    self.genes.add(n)

    def add_edge(self, nodes):
        try:
            if nodes[0] in self.genes and nodes[1] in self.genes:
                nodes.sort()
                self.edges.add(",".join(nodes))
        except IndexError:
            print(nodes)

    def __iter__(self):
        for e in self.edges:
            yield e


def parse_options(args):
    parser = argparse.ArgumentParser(description='Get an induced subgraph from a csv.')
    parser.add_argument('-f', type=str, help='input file', required=True, metavar='input-file')
    parser.add_argument('-s', type=str, help='size of graph', required=True, metavar='size')
    opts = parser.parse_args(args)
    return opts.f, opts.s


def main(args):
    fl, si = parse_options(args[1:])
    g = Graph(si)
    with open(fl, "r") as infile:
        for ln in infile:
            toks = ln.rstrip().split("\t")
            if len(toks) > 1:
                g.add_node(toks[0:2])
                g.add_edge(toks[0:2])
    for e in g:
        print(e)

if __name__ == "__main__":
    main(sys.argv)
