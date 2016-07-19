package net.digital_alexandria.util;

import java.util.Random;

/**
 * @author Simon Dirmeier {@literal simon.dirmeier@gmx.de}
 */
public class StdRnd
{

    private static Random rnd;

    static {
        rnd = new Random(23);
    }

    public static double runif()
    {
        return rnd.nextDouble();
    }
}
