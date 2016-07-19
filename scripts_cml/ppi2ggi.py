# __author__ = 'Simon Dirmeier'
# __email__  = 'simon.dirmeier@bsse.ethz.ch'
# __date__   = 06/04/16

from __future__ import print_function, absolute_import
import sys
import argparse
import re


def parse_options(args):
    parser = argparse.ArgumentParser(description='Get a mapping from PPI to GGI')
    parser.add_argument('-p', type=str, help='ppi file (e.g. homo_sapiens_ppi.v10.tsv)', required=True,
                        metavar='ppi-file')
    parser.add_argument('-m', type=str, help='mapping file that maps proteins to genes (e.g. entrez2string.v10.tsv)', required=True,
                        metavar='map-file')
    opts = parser.parse_args(args)
    return opts.p, opts.m


def multi_map(mf, sp):
    hash = {}
    with open(mf, "r") as infile:
        for line in infile:
            entrez, ensemble = line.strip().split(sp)
            ensemble = re.sub(r"^\d+.", "", ensemble)
            hash[ensemble] = entrez
    return hash


def get_induced_subgraph(ppif, map, sp):
    print("Entrez\tEntrez\tScore")
    with open(ppif, "r") as infile:
        for line in infile:
            ens1, ens2, score = line.strip().split(sp)
            ens1 = re.sub(r"^\d+.", "", ens1)
            ens2 = re.sub(r"^\d+.", "", ens2)
            if ens1 in map and ens2 in map:
                print(map[ens1], "\t", map[ens2], "\t", score)



def main(args):
    ppif, mf = parse_options(args[1:])
    map = multi_map(mf, "\t")
    get_induced_subgraph(ppif, map, "\t")

if __name__ == "__main__":
    main(sys.argv)
