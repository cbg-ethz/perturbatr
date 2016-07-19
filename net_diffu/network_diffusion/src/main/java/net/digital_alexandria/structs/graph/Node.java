package net.digital_alexandria.structs.graph;

import net.digital_alexandria.structs.table.DoubleRow;
import net.digital_alexandria.util.StdRnd;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

/**
 * @author Simon Dirmeier {@literal simon.dirmeier@gmx.de}
 */
public class Node<T extends Comparable>
{
    private final T _ID;
    private Set<String> _altIDs;
    private double[] _vals;
    private ArrayList<Edge> _edges;
    private boolean _isStart;
    private int cntWalkedBy;
    private boolean _frequent;
    private double _effect;

    Node(T id)
    {
        this._ID = id;
        this._altIDs = new HashSet<>();
        this._edges = new ArrayList<>();
        this._effect = 0.0;
    }

    Node(T id, DoubleRow row)
    {
        this(id);
        this._vals = row.values();
    }


    public ArrayList<Edge> edges()
    {
        return this._edges;
    }


    public void incWalkBy()
    {
        this.cntWalkedBy++;
    }

    @Override
    public boolean equals(Object o)
    {
        return o instanceof Node && this._ID.equals(((Node) o)._ID);
    }

    public String toString()
    {
        return String.valueOf(this._ID);
    }

    @Override
    public int hashCode()
    {
        return this._ID.hashCode() % 7;
    }

    void addOutgoing(Edge e)
    {
        this._edges.add(e);
    }

    void normalizeEdgeWeights()
    {
        double s = 0.0;
        for (Edge e : this._edges)
            s += e.weight();
        for (Edge e : this._edges)
            e.weight(e.weight() / s);
        int k = 2;
    }

    public T id()
    {
        return _ID;
    }

    public int degree()
    {
        return this._edges.size();
    }

    public double[] values()
    {
        return this._vals;
    }

    public Node random()
    {
        double rnd = StdRnd.runif();
        double stp = 0.0;
        for (Edge e : _edges)
        {
            stp += e.weight();
            if (rnd <= stp) return e.other(this);
        }
        return _edges.get(_edges.size() - 1).other(this);
    }

    public void setStart()
    {
        this._isStart = true;
    }

    public boolean isStart()
    {
        return _isStart;
    }

    public int walkedBy()
    {
        return this.cntWalkedBy;
    }

    public void setFrequent(boolean b)
    {
        this._frequent = b;
    }

    public boolean isFrequent()
    {
        return _frequent;
    }

    public void weight(double effect)
    {
        this._effect = effect;
    }

    public double weight()
    {
        return _effect;
    }

    public void addAlternativeIDs(Set<String> ids)
    {
        this._altIDs.addAll(ids);
    }

    public Set<String> alternativeIDs()
    {
        return _altIDs;
    }
}
