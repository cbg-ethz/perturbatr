package net.digital_alexandria.diffusion.kernels;

/**
 * @author Simon Dirmeier {@literal simon.dirmeier@gmx.de}
 */
public class GaussianKernel implements Kernel
{

    private static GaussianKernel kernel;

    private GaussianKernel()
    {
    }

    public static GaussianKernel getInstance()
    {
        if (kernel == null)
            kernel = new GaussianKernel();
        return kernel;
    }

    @Override
    public double map(double[] x1, double[] x2)
    {

        double ret = 0.0;
        for (int i = 0; i < x1.length; i++)
        {
            double sub = x1[i] - x2[i];
            ret += sub*sub;
        }
        ret /= 2.0;
        ret = Math.exp(-ret);
        return ret;
    }
}
