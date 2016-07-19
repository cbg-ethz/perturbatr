package net.digital_alexandria.structs.maps;

import java.util.*;

/**
 * @author Simon Dirmeier {@literal simon.dirmeier@gmx.de}
 */
public class MultiMap<T, U> implements Iterable<Map.Entry<U, Set<T>>>
{
    private Map<T, Set<U>> _toMap;
    private HashMap<U, Set<T>> _fromMap;

    MultiMap()
    {
        this._fromMap = new HashMap<>();
        this._toMap = new HashMap<>();
    }

    void put(T t, U u)
    {
        if (!_toMap.containsKey(t))
        {
            _toMap.put(t, new HashSet<>());
        }
        if (!_fromMap.containsKey(u))
        {
            _fromMap.put(u, new HashSet<>());
        }

        this._toMap.get(t).add(u);
        this._fromMap.get(u).add(t);
    }

    @Override
    public Iterator<Map.Entry<U, Set<T>>> iterator()
    {
        return  _fromMap.entrySet().iterator();
    }

    public boolean containsFrom(U key)
    {
        return _fromMap.containsKey(key);
    }

    public Set<T> getFrom(U key)
    {
        return _fromMap.get(key);
    }

    public Set<U> secondKeys()
    {
        return _fromMap.keySet();
    }

    public Set<U> firstValues(T key)
    {
        return _toMap.get(key);
    }

    public Set<T> firstKeys()
    {
        return _toMap.keySet();
    }

    public Set<U> getTo(T k)
    {
        return this._toMap.get(k);
    }
}
