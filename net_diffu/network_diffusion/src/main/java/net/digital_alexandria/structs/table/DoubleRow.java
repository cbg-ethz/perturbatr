package net.digital_alexandria.structs.table;

/**
 * @author Simon Dirmeier {@literal simon.dirmeier@gmx.de}
 */
public class DoubleRow
{

    private String _id;
    private double[] _vals;

    public DoubleRow(int le)
    {
        this._vals = new double[le];
    }

    public void set(String[] rowNames, DoubleColumn[] readouts, int idx)
    {
        this._id = rowNames[idx];
        for (int i = 0; i < readouts.length; i++)
        {
            this._vals[i] = readouts[i].get(idx);
        }
    }

    public String toString()
    {
        StringBuilder sb = new StringBuilder(this._id).append("\t");
        for (int i = 0; i < this._vals.length; i++)
        {
            sb.append(this._vals[i]);
            if (i < this._vals.length - 1)
                sb.append(" ");
        }
        return sb.toString();
    }

    public double[] values()
    {
        return _vals;
    }

    public String id()
    {
        return _id;
    }
}
