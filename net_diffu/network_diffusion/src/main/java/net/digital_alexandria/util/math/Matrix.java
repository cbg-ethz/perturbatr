package net.digital_alexandria.util.math;

/**
 * @author Simon Dirmeier {@literal simon.dirmeier@gmx.de}
 */
public final class Matrix
{

    private Matrix () {}

    public static double[][] columnNormalize(double[][] m)
    {
        int nrow = m.length;
        int ncol = m[0].length;
        for (int j = 0; j < ncol; j++)
        {
            double colSum = 0.0;
            for (int i = 0; i < nrow; i++)
            {
                    colSum += m[i][j];
            }
            for (int i = 0; i < nrow; i++)
            {
                m[i][j] /= colSum;
            }
        }
        return m;
    }

    public static double[] multiply(final double[][] M, final double[] v)
    {
        double w[] = new double[v.length];
        for (int i = 0; i < M.length; i++)
        {
            double dotp = 0.0;
            for (int j = 0; j < v.length; j++)
            {
                dotp += M[i][j] * v[j];
            }
            w[i] = dotp;
        } 
        return w;
    }
}
