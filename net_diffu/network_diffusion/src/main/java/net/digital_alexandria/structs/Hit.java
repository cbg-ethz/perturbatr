package net.digital_alexandria.structs;

/**
 * @author Simon Dirmeier {@literal simon.dirmeier@gmx.de}
 */
public class Hit<T extends Comparable> implements Comparable<Hit<T>>
{
    private T _id;
    private double _score;


    public Hit(T id, double score)
    {
        this._id = id;
        this._score = score;
    }

    @Override
    public int compareTo(Hit<T> other)
    {
        if (this._score < other._score) return 1;
        else if (this._score > other._score) return -1;
        return 0;
    }

    public double score()
    {
        return this._score;
    }

    public T id()
    {
        return _id;
    }
}
