package net.digital_alexandria;

import net.digital_alexandria.diffusion.NetworkPropagator;
import net.digital_alexandria.diffusion.NetworkPropagatorFactory;
import net.digital_alexandria.param.FlagType;
import net.digital_alexandria.param.Param;
import net.digital_alexandria.param.ParamList;
import net.digital_alexandria.param.ParamsParser;
import net.digital_alexandria.structs.Hit;
import net.digital_alexandria.structs.enums.DataType;
import net.digital_alexandria.structs.enums.HeaderType;
import net.digital_alexandria.structs.graph.Edge;
import net.digital_alexandria.structs.graph.Graph;
import net.digital_alexandria.structs.graph.GraphFactory;
import net.digital_alexandria.structs.graph.Node;
import net.digital_alexandria.structs.maps.MapFactory;
import net.digital_alexandria.structs.maps.MultiMap;
import net.digital_alexandria.util.FileReader;
import net.digital_alexandria.util.annotations.ParamAnnotation;
import net.digital_alexandria.util.reflections.FieldParser;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.*;

/**
 * @author Simon Dirmeier {@literal simon.dirmeier@gmx.de}
 */
final class Controller
{
    private static final boolean DEBUG = false;
    private final static Logger _LOGGER = LoggerFactory.getLogger(Controller.class);

    @ParamAnnotation(value = "-t")
    private String _DATA_TABLE_FILE;
    @ParamAnnotation(value = "-g")
    private String _GGI_FILE;
    @ParamAnnotation(value = "-l")
    private String _GENE_SCORES;
    @ParamAnnotation(value = "-o")
    private String _OUT_FILE;
    // the protein protein interaction graph using gene IDs
    private Graph<Integer, Integer> _ggiGraph;
    // hugo symbol to entrez mapping from data
    private MultiMap<String, Integer> _entrez2HugoDualMap;
    // the markov random walk object


    Controller(String[] args)
    {
        // parses the command line arguments
        _parse(args);
        // a mapping from hugoSymbol <-> entrezID (from screening data)
        _initMappings();
        // create a GGI graph using the entrez IDs from the data
        _initGGI();
    }

    private void _parse(String... args)
    {
        ParamsParser parser = _setArgs();
        parser.parse(args);
        FieldParser.parse(this, Controller.class, parser);
    }

    private void _initMappings()
    {
        this._entrez2HugoDualMap = MapFactory.newMultiMap(DataType.DataTable,
                                                          this._DATA_TABLE_FILE);
    }

    private void _initGGI()
    {
        // The set of entrez keys in the study
        Set<Integer> availabelEntrezKeys =
            new HashSet<>(this._entrez2HugoDualMap.secondKeys());
        // create ggi using data genes only
        // read ggi graph to memory
        this._ggiGraph = GraphFactory.newGraph(
            this._GGI_FILE, availabelEntrezKeys,
            _entrez2HugoDualMap);
        // set and  gene effects
        this._ggiGraph.setGeneScores(
            this._entrez2HugoDualMap,
            FileReader.readTSV(this._GENE_SCORES,
                               HeaderType.HeaderOneLine,
                               new String[]{"GeneID", "Effect"},
                               "\t"));

        if (DEBUG)
        {
            System.out.println("Edges");
            for (Node<?> n : _ggiGraph.nodes())
            {
                for (Edge e : n.edges())
                {
                    System.out.println(n + " " + e.other(n));
                }
            }
        }
    }

    final void run()
    {
        NetworkPropagator<Integer> mrw;
        mrw = NetworkPropagatorFactory.newMarkovRandomWalk(this._ggiGraph);
        List<Hit<Integer>> m =  mrw.diffuse();
        _print(m);
    }

    private void _print(List<Hit<Integer>> hits)
    {
        Collections.sort(hits);
        for (Hit<Integer> h : hits)
        {
            Set<String> hugos = this._entrez2HugoDualMap.getFrom(h.id());
            System.out.println(h.id() + "\t" + h.score() + "\t" + hugos);
        }
    }

    private ParamsParser _setArgs()
    {
        ParamsParser parser = ParamsParser
            .newInstance(
                ParamList.newInstance(
                    new Param.Builder()
                        .flag("-g")
                        .desc("the ggi file")
                        .flagType(FlagType.NEEDED_ARG)
                        .build(),
                    new Param.Builder()
                        .flag("-l")
                        .desc("the file of genes scores")
                        .flagType(FlagType.NEEDED_ARG)
                        .build(),
                    new Param.Builder()
                        .flag("-t")
                        .desc("rnai data table")
                        .flagType(FlagType.NEEDED_ARG)
                        .build(),
                    new Param.Builder()
                        .flag("-o")
                        .desc("outfile")
                        .flagType(FlagType.NEEDED_ARG)
                        .build()),
                "java -jar network_diffusion.jar");
        return parser;
    }

}
