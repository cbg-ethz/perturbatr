package net.digital_alexandria.structs.table;

import java.util.Arrays;
import java.util.Iterator;

/**
 * @author Simon Dirmeier {@literal simon.dirmeier@gmx.de}
 */
public class DataTable implements Iterable<DoubleRow>
{

    private static final int _INIT_NROW = 100;
    private static final double _REALLOC_SZ = 1.6;

    private DoubleColumn[] _readouts;
    private String[] _rowNames;
    private int _nrow;

    public static DataTable newInstance()
    {
        return new DataTable();
    }

    private DataTable()
    {
    }

    public void rbind(String[] toks)
    {
        if (_readouts == null) _init(toks.length - 1);
        if (_nrow == this._rowNames.length - 1) _resize();
        _rbind(toks);
    }

    private void _init(int cols)
    {
        this._rowNames = new String[_INIT_NROW];
        this._readouts = new DoubleColumn[cols];
        for (int i = 0; i < this._readouts.length; i++)
            this._readouts[i] = new DoubleColumn();
    }

    private void _rbind(String[] toks)
    {
        String geneName = toks[0];
        for (int i = 0; i < toks.length - 1; i++)
            this._readouts[i].add(toks[i+1]);
        this._rowNames[this._nrow] = geneName;
        this._nrow++;
    }

    private void _resize()
    {
        this._rowNames = Arrays.copyOf(this._rowNames,
                                       (int) (this._rowNames.length * _REALLOC_SZ));
    }

    @Override
    public Iterator<DoubleRow> iterator()
    {
        return new Iterator<DoubleRow>() {

            private final int _SZ = _nrow;
            private DoubleRow _row = new DoubleRow(_readouts.length);
            private int _run = 0;

            @Override
            public boolean hasNext()
            {
                return _run < _nrow;
            }

            @Override
            public DoubleRow next()
            {
                _row.set(_rowNames, _readouts, _run);
                _run++;
                return _row;
            }
        };
    }
}
