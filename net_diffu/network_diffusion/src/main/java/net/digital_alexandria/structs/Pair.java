package net.digital_alexandria.structs;

/**
 * @author Simon Dirmeier {@literal simon.dirmeier@gmx.de}
 */
public final class Pair<T, U>
{
    public final T t;
    public final U u;

    public Pair(T t, U u)
    {
        this.t = t;
        this.u = u;
    }
}
