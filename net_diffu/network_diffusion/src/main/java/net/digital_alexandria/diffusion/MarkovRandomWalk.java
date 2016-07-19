package net.digital_alexandria.diffusion;

import net.digital_alexandria.structs.Hit;
import net.digital_alexandria.structs.graph.Edge;
import net.digital_alexandria.structs.graph.Graph;
import net.digital_alexandria.structs.graph.Node;
import net.digital_alexandria.util.math.Matrix;
import net.digital_alexandria.util.math.Vector;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author Simon Dirmeier {@literal simon.dirmeier@gmx.de}
 */
final class MarkovRandomWalk<T extends Comparable, U extends Comparable> implements NetworkPropagator<T>
{

    private final static Logger _LOGGER = LoggerFactory.getLogger(Graph.class);
    private static final boolean DEBUG = false;

    private static final double _RESTART_PROBABILITY = 0.5;
    private static final double _STEP_PROBABILITY = 1 - _RESTART_PROBABILITY;
    private static final double _THRESHOLD = 1e-5;

    private final Graph<T, U> _G;
    private final Map<T, Integer> _NODE_INDEXES;
    private final List<Node<T>> _NODES;
    private final double _P0[];
    private double _P_INF[];
    private final double _W[][];

    MarkovRandomWalk(Graph<T, U> g)
    {
        this._G = g;
        this._NODES = new ArrayList<>(this._G.nodes());
        this._NODE_INDEXES = _setNodeMap(this._NODES);
        this._P0 = _setStartDistribution(this._NODES);
        this._W = _setTransitionMatrix(this._NODES);

        if (DEBUG)
            _dprint();

    }

    private void _dprint()
    {
        System.out.println("Node mapping");
        for (Node<T> n : new ArrayList<>(this._G.nodes()))
        {
            System.out.println(this._NODE_INDEXES.get(n.id()) + " " + n + " " + n.alternativeIDs());
        }
        System.out.println("Starting distribution");
        for (int i = 0; i < this._P0.length; i++)
        {
            System.out.print(_P0[i] + " ");
        }
        System.out.println("\nTransition matrix");
        for (int i = 0; i < _W.length; i++)
        {
            for (int j = 0; j < _W[i].length; j++)
            {
                System.out.print(_W[i][j] + " ");
            }
            System.out.println();
        }
    }

    private Map<T, Integer> _setNodeMap(List<Node<T>> nodes)
    {
        Map<T, Integer> map = new HashMap<>();
        int sz = nodes.size();
        for (int i = 0; i < sz; i++)
        {
            map.put(nodes.get(i).id(), i);
        }
        return map;
    }

    private double[] _setStartDistribution(List<Node<T>> nodes)
    {
        int sz = nodes.size();
        double[] po = new double[sz];
        for (int i = 0; i < sz; i++)
        {
            po[i] = nodes.get(i).weight();
        }
        po = Vector.normalize(po);
        return po;
    }

    private double[][] _setTransitionMatrix(List<Node<T>> nodes)
    {
        int sz = nodes.size();
        double[][] w = new double[sz][sz];
        for (int i = 0; i < sz; i++)
        {
            Node<T> n = nodes.get(i);
            for (Edge e : n.edges())
            {
                Node<?> other = e.other(n);
                int j = this._NODE_INDEXES.get(other.id());
                w[i][j] = e.weight();
            }
        }
        w = Matrix.columnNormalize(w);
        return w;
    }

    @Override
    public List<Hit<T>> diffuse()
    {
        this._P_INF = new double[this._P0.length];
        System.arraycopy(this._P0, 0, this._P_INF, 0, this._P_INF.length);
        double[] pold = new double[this._P0.length];
        do
        {
            System.arraycopy(this._P_INF, 0, pold, 0, this._P_INF.length);
            _P_INF = walk(this._W, pold, this._RESTART_PROBABILITY, this._P0);
            int k = 2;
        }
        while (!Vector.converges(this._P_INF, pold, this._THRESHOLD));
        if (DEBUG)
        {
            System.out.println("Stationary:");
            for (int i = 0; i < this._P_INF.length; i++)
            {
                System.out.println(_P_INF[i]);
            }
        }

        return map(this._P_INF);
    }

    private double[] walk(final double[][] W, final double[] pt, final double r, final double[] p0)
    {
        // r * W %*% p_t
        double[] step = Vector.multiply(_STEP_PROBABILITY, Matrix.multiply(W, pt));
        double[] restart = Vector.multiply(_RESTART_PROBABILITY, p0);
        double[] ptnew = Vector.add(step, restart);
        return ptnew;
    }

    private List<Hit<T>> map(double[] v)
    {
        List<Hit<T>> l = new ArrayList<>();
        for (int i = 0; i < v.length; i++)
            l.add(new Hit<>(this._NODES.get(i).id(), v[i]));
        return l;
    }

}
