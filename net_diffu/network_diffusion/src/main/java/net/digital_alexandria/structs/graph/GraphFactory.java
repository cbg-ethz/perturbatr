package net.digital_alexandria.structs.graph;

import net.digital_alexandria.diffusion.kernels.GaussianKernel;
import net.digital_alexandria.diffusion.kernels.Kernel;
import net.digital_alexandria.structs.maps.MultiMap;
import net.digital_alexandria.structs.table.DataTable;
import net.digital_alexandria.structs.table.DoubleRow;
import net.digital_alexandria.util.math.Sets;
import org.slf4j.LoggerFactory;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Set;

/**
 * @author Simon Dirmeier {@literal simon.dirmeier@gmx.de}
 */
public class GraphFactory
{

    private GraphFactory() { }

    private final static org.slf4j.Logger _LOGGER = LoggerFactory.getLogger(GraphFactory.class);

    public static Graph<Integer, Integer> newGraph(String tsvFile, Set<Integer> relNodes,
                                                   MultiMap<String, Integer> _altIDs)
    {
        return newGraphFromTSV(tsvFile, relNodes, _altIDs);
    }

    private static Graph<Integer, Integer> newGraphFromTSV(String tsvFile, Set<Integer> relNodes,
                                                           MultiMap<String, Integer> _altIDs)
    {

        Graph<Integer, Integer> g = _initGraph(tsvFile, relNodes, _altIDs);
        return g;
    }

    private static Graph<Integer, Integer> _initGraph(String tsvFile, Set<Integer> relNodes,
                                                      MultiMap<String, Integer> _altIDs)
    {
        Graph<Integer, Integer> g = new Graph<>();
        BufferedReader bR;
        try
        {
            bR = new BufferedReader(new FileReader(new File(tsvFile)));
            String line;
            int run = 0;
            while ((line = bR.readLine()) != null)
            {
                if (!line.matches("^\\d+.*\\d+.*\\d+$")) continue;
                String[] toks = line.split("\t");
                int node1 = Integer.parseInt(toks[0].trim());
                int node2 = Integer.parseInt(toks[1].trim());
                if (!relNodes.contains(node1) ||
                    !relNodes.contains(node2))
                {
                    continue;
                }
                double weight = Double.parseDouble(toks[2]);
                _addToGraph(g, node1, node2, run++, weight);
                g.node(node1).addAlternativeIDs(_altIDs.getFrom(node1));
                g.node(node2).addAlternativeIDs(_altIDs.getFrom(node2));
            }
            bR.close();
        }
        catch (IOException e)
        {
            _LOGGER.error("IO-error");
        }
        Set<Integer> missingEls = Sets.difference(relNodes, g.nodeIDs());
        for (int mis : missingEls)
        {
           // _LOGGER.warn("Could not find: " + mis + " for GGI!");
        }
        return g;
    }

    private static <T extends Comparable, U extends Comparable>
    void _addToGraph(Graph<T, U> g, T toNodeLabel, T fromNodeLabel, U edgeLabel, double edgeWeight)
    {

        Node<T> n1 = g.getOrPut(toNodeLabel);
        Node<T> n2 = g.getOrPut(fromNodeLabel);
        Edge<U> e = new Edge<>(edgeLabel, n1, n2, edgeWeight);
        g.add(toNodeLabel, n1);
        g.add(fromNodeLabel, n2);
        g.add(edgeLabel, e);
    }

    public static Graph<Integer, Integer> newGraph(DataTable table)
    {
        Graph<Integer, Integer> g = new Graph<>();
        Kernel kernel = GaussianKernel.getInstance();
        int ncnt = 1;
        int ecnt = 1;
        for (DoubleRow row : table)
        {
            Node<Integer> n = new Node<>(ncnt, row);
            g.put(ncnt, n);
            ncnt++;
        }
        for (Node<Integer> n : g.nodes())
        {
            if (n.degree() >= 3) continue;
            for (Node m : g.nodes())
            {
                if (m.degree() >= 5) continue;
                if (n == m) continue;
                Edge<Integer> e = new Edge<>(ecnt,
                                             n, m, kernel.map(n.values(), m.values()));
                if (!g.contains(e))
                {
                    g.put(ecnt, e);
                    n.addOutgoing(e);
                    m.addOutgoing(e);
                    ecnt++;
                }
            }
        }
        g.normalizeEdgeWeights();
        return g;
    }
}
