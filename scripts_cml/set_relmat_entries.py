# __author__ = 'Simon Dirmeier'
# __email__  = 'simon.dirmeier@bsse.ethz.ch'
# __date__   = 01/06/16

import sys
import argparse
import numpy as np


def parse_options(args):
    parser = argparse.ArgumentParser(description='Set siRNA-entrez pairs to a specific values')
    parser.add_argument('-f', type=str, help='target relation matrix', required=True,
                        metavar='matrix')
    parser.add_argument('-r', type=str, help='rownames of target relation matrix', required=True,
                        metavar='rownames')
    parser.add_argument('-c', type=str, help='colnames of relation matrix', required=True,
                        metavar='colnames')
    parser.add_argument('-m', type=str, help='siRNA-entrez mappings', required=True,
                        metavar='map-file')
    parser.add_argument('-o', type=str, help='outfile', required=True,
                        metavar='map-file')
    opts = parser.parse_args(args)
    return opts.f, opts.r, opts.c, opts.m, opts.o


def read_sirna_ent_list(mf):
    m = []
    with open(mf, "r") as f:
        f.readline().strip().split("\t")
        for l in f.readlines():
            tok = l.strip().split("\t")
            m.append(tok)
    return m


def read_list(lf):
    with open(lf, "r") as f:
        l = [line.strip() for line in f]
    return l


def read_matrix(mf):
    x = np.loadtxt(mf, skiprows=1, delimiter="\t")
    return x


def idx(arr, val):
    try:
        idx = arr.index(val)
    except ValueError:
        print("Could not find", val)
        idx = -1
    return idx


def set_relmat_els(X, rownames, colnames, sirna_entrez_list, of):
    for el in sirna_entrez_list:
        ent = el[0]
        sirnas = el[1].split(",")
        ent_idx = idx(colnames, ent)
        for sirna in sirnas:
            sirna_idx = idx(rownames, sirna)
            if ent_idx != -1 and sirna_idx != -1:
                X[sirna_idx, ent_idx] = 0.75
    np.savetxt(of, X, delimiter="\t")


def run(mat, rs, cls, ma, of):
    se_l = read_sirna_ent_list(ma)
    rownames = read_list(rs)
    colnames = read_list(cls)
    X = read_matrix(mat)
    set_relmat_els(X, rownames, colnames, se_l, of)
#

def main(args):
    mat, rs, cls, ma, of = parse_options(args[1:])
    run(mat, rs, cls, ma, of)


if __name__ == "__main__":
    main(sys.argv)
