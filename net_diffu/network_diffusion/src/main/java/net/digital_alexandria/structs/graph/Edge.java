package net.digital_alexandria.structs.graph;

import java.util.Comparator;

/**
 * @author Simon Dirmeier {@literal simon.dirmeier@gmx.de}
 */
public class Edge<T extends Comparable>  implements Comparable<Edge>
{
    private final T _ID;
    private final Node _FROM;
    private final Node _TO;
    private double _weight;

    Edge(T id, Node from, Node to, double weight)
    {
        this._ID = id;
        this._FROM = from;
        this._TO = to;
        this._weight = weight;
    }

    @Override
    public boolean equals(Object o)
    {
        if (o instanceof Edge)
        {
            Edge e = (Edge) o;
            return (e._FROM.id() == _FROM.id() && e._TO.id() == _TO.id() )||
                   (e._FROM.id() == _TO.id() && e._TO.id() == _FROM.id());
        }
        return false;
    }

    @Override
    public int hashCode()
    {
        return (_FROM.hashCode() + _TO.hashCode()) % 13;
    }

    public double weight()
    {
        return _weight;
    }

    public void weight(double w)
    {
        this._weight = w;
    }

    @Override
    public int compareTo(Edge e)
    {
        return this._weight - e._weight < 0 ? 1 : -1;
    }

    public T id()
    {
        return _ID;
    }

    public Node to()
    {
        return _TO;
    }

    public Node from() { return _FROM; }

    public Node other(Node node)
    {
        return node == _FROM ? _TO : _FROM;
    }
}
