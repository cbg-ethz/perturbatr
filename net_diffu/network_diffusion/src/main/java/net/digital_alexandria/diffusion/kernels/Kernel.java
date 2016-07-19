package net.digital_alexandria.diffusion.kernels;

/**
 * @author Simon Dirmeier {@literal simon.dirmeier@gmx.de}
 */
public interface Kernel
{
    public double map(double[] x1, double[] x2);
}
