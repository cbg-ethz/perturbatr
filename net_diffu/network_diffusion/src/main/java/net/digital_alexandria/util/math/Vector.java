package net.digital_alexandria.util.math;

/**
 * @author Simon Dirmeier {@literal simon.dirmeier@gmx.de}
 */
public final class Vector
{
    private Vector() {}

    public static double[] normalize(double[] v)
    {
        double sum = 0.0;
        for (int i = 0; i < v.length; i++)
        {
            sum += v[i];
        }
        for (int i = 0; i < v.length; i++)
        {
            v[i] /= sum;
        }
        return v;
    }

    public static boolean converges(final double[] a, final double[] b,  final double threshold)
    {
        double euc = 0.0;
        for (int i = 0; i < a.length; i++)
        {
            double dif = a[i] - b[i];
            euc += dif * dif;
        }
        return Math.sqrt(euc) < threshold;
    }


    public static double[] multiply(final double r, final double[] v)
    {
        double[] w = new double[v.length];
        for (int i = 0; i < v.length; i++)
        {
            w[i] = v[i] * r;
        }
        return w;
    }

    public static double[] add(final double[] a, final double[] b)
    {
        double[] n = new double[a.length];
        for (int i = 0; i < n.length; i++)
        {
            n[i] = a[i] + b[i];
        }
        return n;
    }
}
