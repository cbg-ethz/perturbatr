package net.digital_alexandria.util.math;

import java.util.HashSet;
import java.util.Set;

/**
 * @author Simon Dirmeier {@literal simon.dirmeier@gmx.de}
 */
public class Sets
{

    private Sets()
    {

    }

    public static <T> Set<T> difference(Set<T> a, Set<T> b)
    {
        Set<T> missing = new HashSet<>();
        for (T el : a)
        {
            if (!b.contains(el)) missing.add(el);
        }
        return missing;
    }
}
