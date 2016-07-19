package net.digital_alexandria.diffusion.diffusors;

import net.digital_alexandria.structs.graph.Graph;
import net.digital_alexandria.structs.graph.Node;


import java.util.Random;


/**
 * @author Simon Dirmeier {@literal simon.dirmeier@gmx.de}
 */
public class VizMarkovRandomWalk<T extends Comparable>
{
    private static final double _THRESHOLD = 0.00000000000001;
    private static final int _SEED = 23;
    private static final Random _RND = new Random(_SEED);
    private Node<T> _currNode;
    private final Node<T> _START_NODE;
    private int[] _counts;
    private final double[] _P0;
    private double[] _Pt;
    private double[] _Pold;
    private int nsteps;
    private final static double _RESTART_PROB = 0.2;

    public VizMarkovRandomWalk(Graph graph)
    {
        this._counts = new int[graph.size()];
        this._currNode = graph.random();
        this._START_NODE = _currNode;
        this._START_NODE.setStart();
        this._P0 = new double[this._counts.length];
        this._P0[0] = 1.0;
    }

    public VizMarkovRandomWalk walk()
    {
        nsteps++;
        _Pold = _stationaryDistribution();
        if (_RND.nextDouble() < _RESTART_PROB)
            _currNode = _START_NODE;
        else
            _currNode = _currNode.random();
        _currNode.incWalkBy();
        //_counts[_currNode.id() - 1]++;
        _Pt = _stationaryDistribution();
        return this;
    }

    public Node<T> state()
    {
        return _currNode;
    }

    public void stats()
    {
        _stationaryDistribution();
    }

    private double[] _stationaryDistribution()
    {
        double[] dist = new double[_counts.length];
        int sum = 0;
        for (int i = 0; i < _counts.length; i++)
            sum += _counts[i];
        for (int i = 0; i < dist.length; i++)
        {
            dist[i] = _counts[i] * 1.0 / sum;
        }
        return dist;
    }

    public boolean converges()
    {
        double sum = 0.0;
        for (int i = 0; i < this._Pt.length; i++)
        {
            sum += Math.pow(_Pt[i] - _Pold[i], 2);
        }
        return sum < _THRESHOLD;
    }

    public void printStationaryDistribution()
    {
        for (int i = 0; i < _Pt.length; i++)
        {
            System.out.print(_Pt[i]+ " ");
        }
        System.out.println();
    }

    public int steps()
    {
        return nsteps;
    }
}
