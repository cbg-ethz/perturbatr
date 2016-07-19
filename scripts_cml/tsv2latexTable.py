# __author__ = 'Simon Dirmeier'
# __email__  = 'simon.dirmeier@bsse.ethz.ch'
# __date__   = 22/03/16


from __future__ import print_function, absolute_import
import argparse
import os

__DATA_KEY__ = "d"
__HEADER_KEY__ = "h"


def parse_options():
    parser = argparse.ArgumentParser(description='Parse an R table file to a latex table as pdf.')
    parser.add_argument('-f', type=str, help='input file as table with header', required=True, metavar='input-file')
    parser.add_argument('-o', type=str, help='output folder', required=True, metavar='output-folder')
    parser.add_argument('-s', type=str, help='split at group in index', required=True, metavar='rowidx-split')
    opts = parser.parse_args()
    return opts.f, opts.o, opts.s


def read_table(infile, spl):
    header = []
    spl = int(spl)
    data = dict()
    with open(infile, "r") as inhandle:
        header = inhandle.readline().strip().split("\t")
        for line in inhandle.readlines():
            toks = line.strip().split("\t")
            vir = toks[spl]
            if vir not in data:
                data[vir] = []
            data[vir].append(" & ".join(toks))
    m = {__DATA_KEY__: data, __HEADER_KEY__: header}
    return m


def write_latex(m, outfolder):
    if not os.path.isdir(outfolder):
        os.makedirs(outfolder)
    tex = outfolder + "/table.tex"
    with open(tex, "w") as tmp:
        tmp.write("\\documentclass[9pt]{article}\n")
        tmp.write("\\usepackage{booktabs}\n")
        tmp.write("\\newcommand{\\ra}[1]{\\renewcommand{\\arraystretch}{#1}}\n")
        tmp.write("\\begin{document}\n")
        for spl in m[__DATA_KEY__].keys():
            write_split(tmp, m[__HEADER_KEY__], m[__DATA_KEY__][spl])
        tmp.write("\\end{document}\n")


def write_split(tmp, header, data):
    tmp.write("\\begin{table} \centering\n")
    tmp.write("\\ra{1.3}\n")
    tmp.write("\\begin{tabular}{@{}")
    for i in header:
        tmp.write("l")
    tmp.write("@{}}\n")
    tmp.write("\\toprule\n")
    tmp.write(" & ".join(header))
    tmp.write(" \\\\ \n \\midrule \small \n")
    tmp.write(" \\\ \n".join(data))
    tmp.write(" \\\ \n \\bottomrule\n")
    tmp.write("\\end{tabular}\n")
    tmp.write("\\end{table}\n")


def write_pdf(infile, outfolder, spl):
    m = read_table(infile, spl)
    write_latex(m, outfolder)


def main():
    infile, outfolder, spl = parse_options()
    write_pdf(infile, outfolder, spl)


if __name__ == "__main__":
    main()
