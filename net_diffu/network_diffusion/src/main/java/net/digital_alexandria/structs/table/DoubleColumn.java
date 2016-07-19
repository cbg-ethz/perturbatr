package net.digital_alexandria.structs.table;

import java.util.Arrays;

/**
 * @author Simon Dirmeier {@literal simon.dirmeier@gmx.de}
 */
public class DoubleColumn
{

    private static final int _INIT_NROW = 100;
    private static final double _REALLOC_SZ = 1.6;

    private double[] _vals;
    private int _nrow;

    public DoubleColumn()
    {
        this._vals = new double[_INIT_NROW];
    }

    public void add(String tok)
    {
        if (this._nrow == this._vals.length - 1) _resize();
        _add(tok);
    }

    private void _add(String tok)
    {
        this._vals[_nrow] = Double.parseDouble(tok);
        this._nrow++;
    }

    public double get(int idx)
    {
        return this._vals[idx];
    }

    private void _resize()
    {
        this._vals = Arrays.copyOf(this._vals, (int) (this._vals.length * _REALLOC_SZ));
    }
}
