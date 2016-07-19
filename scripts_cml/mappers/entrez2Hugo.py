# __author__ = 'Simon Dirmeier'
# __email__  = 'simon.dirmeier@bsse.ethz.ch'
# __date__   = 03/03/16

from __future__ import print_function, absolute_import
from os import listdir
from xlrd import open_workbook
import re
import sys


def parse(fls):
    for f in fls:
        parse_file(f)


def parse_file(layout_file):
    patr = re.compile("^E\w+\d{2}\.xls$")
    files = (layout_file + "/" + f for f in listdir(layout_file) if patr.match(f))
    for f in files:
        read_xls(f, 0)


def read_xls(file, sheet):
    book = open_workbook(file, 'r')
    sheet = book.sheet_by_index(sheet)
    for i in range(sheet.nrows):
        row = sheet.row(i)
        if not str(row[0].value).startswith("EXT"):
            continue
        l = parse_row(row)
        print("\t".join(l))


def parse_row(row):
    toks = [f.value for f in row]
    ret = [str(cast(toks[3])), toks[2].strip()]
    return ret


def cast(val):
    try:
        return int(val)
    except ValueError:
        return val


def main(fls):
    parse(fls)


if __name__ == "__main__":
    main(sys.argv[1:])
