package net.digital_alexandria.mvc.controller;

import edu.uci.ics.jung.graph.UndirectedGraph;
import edu.uci.ics.jung.graph.UndirectedSparseGraph;
import net.digital_alexandria.diffusion.diffusors.VizMarkovRandomWalk;
import net.digital_alexandria.param.FlagType;
import net.digital_alexandria.param.Param;
import net.digital_alexandria.param.ParamList;
import net.digital_alexandria.param.ParamsParser;
import net.digital_alexandria.structs.graph.GraphFactory;
import net.digital_alexandria.structs.table.DataTable;
import net.digital_alexandria.structs.graph.Edge;
import net.digital_alexandria.structs.graph.Graph;
import net.digital_alexandria.structs.graph.Node;
import net.digital_alexandria.util.annotations.ParamAnnotation;

import java.lang.reflect.Field;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * @author Simon Dirmeier {@literal simon.dirmeier@gmx.de}
 */
public class Controller
{
    @ParamAnnotation(value = "-f")
    private String _FILE_NAME;
    private final VizMarkovRandomWalk<Integer> _MRW;
    private final DataTable _TABLE;
    private final Graph<Integer, Integer> _GRAPH;
    private UndirectedGraph _vizGraph;


    public Controller(String[] args)
    {
        _parse(args);
        this._TABLE = net.digital_alexandria.util.FileReader.readDataTable(_FILE_NAME);
        this._GRAPH = GraphFactory.newGraph(_TABLE);
        this._MRW = new VizMarkovRandomWalk(_GRAPH);
        this._vizGraph = _initGraph(this._GRAPH);
    }

    public Controller(Graph g)
    {
        this._GRAPH = g;
        this._vizGraph = _initGraph(this._GRAPH);
        this._MRW = null;
        this._TABLE = null;
    }

    private UndirectedGraph _initGraph(Graph graph)
    {
        UndirectedGraph<Node, Edge> g = new UndirectedSparseGraph<>();
        _GRAPH.nodes().forEach(g::addVertex);
        for (Edge<Integer> e : _GRAPH.edges())
            g.addEdge(e, e.from(), e.to());
        return g;
    }

    public Node<Integer> walk()
    {
        this._MRW.walk();
        if (this._MRW.steps() > 3)
            this._GRAPH.setFrequentNodes();
        return _MRW.state();

    }

    private void _parse(String... args)
    {

        ParamsParser parser = ParamsParser
            .newInstance(
                ParamList.newInstance(
                    new Param.Builder()
                        .flag("-f")
                        .desc("data-matrix")
                        .flagType(FlagType.NEEDED_ARG)
                        .build()),
                "java -jar network_diffusion.jar");
        parser.parse(args);
        for (Field f : Controller.class.getDeclaredFields())
        {
            if (f.isAnnotationPresent(ParamAnnotation.class))
            {
                ParamAnnotation ann = f.getAnnotation(ParamAnnotation.class);
                try
                {
                    f.set(this, parser.getArgument(ann.value()));
                }
                catch (IllegalAccessException e)
                {
                    Logger.getLogger(Controller.class.getSimpleName())
                          .log(Level.SEVERE,
                               "Reflections error parsing argument: " + ann.value());
                }
            }
        }
    }

    public edu.uci.ics.jung.graph.Graph vizGraph()
    {
        return _vizGraph;
    }

    public Graph graph()
    {
        return _GRAPH;
    }

    public void stats()
    {
        this._MRW.stats();
    }
}
