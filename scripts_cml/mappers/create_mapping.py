# __author__ = 'Simon Dirmeier'
# __email__  = 'simon.dirmeier@bsse.ethz.ch'
# __date__   = 03/03/16


__PATH__    = "/Users/simondi/PHD/data/data/sysvirdrug/maps/"
__GENES__   = __PATH__ + "hugo_gene_list.tsv"
__NCBI__    = __PATH__ + "hugo2entrez_ncbi.tsv"
__BIOMART__ = __PATH__ + "hugo2entrez_biomart.tsv"
__SCREEN__  = __PATH__ + "hugo2entrez_virusdata.tsv"
__OUT_COMPL__ =  __PATH__ + "hugo2entrez_complete_final.tsv"
__OUT_SCREEN__ =  __PATH__ + "hugo2entrez_screen_final.tsv"


class Gene:
    def __init__(self, hugo=None, entrez=None):
        self.hugo_ = hugo
        self.entrez_ = entrez

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.hugo_ == other.hugo_
        return False

    def __hash__(self):
        return hash(str(self.hugo_))

    def __str__(self):
        return self.hugo_

    def __repr__(self):
        return self.__str__()


def _read_genes(fl):
    gene_set = set()
    with open(fl) as f:
        for line in f:
            gene_set.add(line.strip())
    return gene_set


def _read_mapping(fl):
    with open(fl) as f:
        gene_map = [line.rstrip().split("\t") for line in f if not line.startswith("Entrez")]
    return gene_map


def _init_gene_set(hugo2entrez_maps):
    h_e_mapping = {}
    for map in hugo2entrez_maps:
        for k, v in map:
            if k != "NA" and k != ""  != "NA" and v != "":
                if v in h_e_mapping and k != h_e_mapping[v]:
                    h_e_mapping[v].add(k)
                else:
                    h_e_mapping[v] = set()
                    h_e_mapping[v].add(k)
    return h_e_mapping


def _init_genes(h_e_mapping):
    genes = []
    for g in h_e_mapping:
        els = list(h_e_mapping[g])
        els.insert(0, g)
        genes.append(els)
    genes = sorted(genes, key=lambda x: len(x), reverse=True)
    return genes


def _write_all(genes):
    with open(__OUT_COMPL__, "w") as f:
            for g in genes:
                f.write("\t".join(g))
                f.write("\n")


def _add(mapping, hugo, entrez):
    if hugo not in mapping:
        mapping[hugo] = entrez
    else:
        entry = mapping[hugo]
        if entry[1] == "ncbi" and entrez[1] == "data":
            mapping[hugo] = entrez

def _init_data_map(h_e_mapping):
    screen_list = _read_mapping(__SCREEN__)
    mapping = {}
    for s in screen_list:
        hugo = s[1]
        entrez = s[0]
        if entrez == "NA" or entrez == "":
            if hugo in h_e_mapping:
                entrez =  [list(h_e_mapping[hugo])[0], "ncbi"]
                _add(mapping, hugo, entrez)

        else:
            entrez = [entrez, "data"]
            _add(mapping, hugo, entrez)
    return mapping


def _write_screen_genes(dat_mapping):
    with open(__OUT_SCREEN__, "w") as f:
        for g in dat_mapping.keys():
            el = dat_mapping[g]
            el.insert(0, g)
            f.write("\t".join(el))
            f.write("\n")


def map_all():
    hugo2entrez_maps = [_read_mapping(x) for x in [__NCBI__, __BIOMART__, __SCREEN__]]
    h_e_mapping = _init_gene_set(hugo2entrez_maps)
    dat_mapping = _init_data_map(h_e_mapping)
    _write_all(_init_genes(h_e_mapping))
    _write_screen_genes(dat_mapping)



def main():
    map_all()


if __name__ == "__main__":
    main()
