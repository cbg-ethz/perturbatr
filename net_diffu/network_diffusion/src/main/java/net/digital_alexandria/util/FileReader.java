package net.digital_alexandria.util;

import net.digital_alexandria.structs.enums.HeaderType;
import net.digital_alexandria.structs.table.DataTable;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.HashMap;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


/**
 * @author Simon Dirmeier {@literal simon.dirmeier@gmx.de}
 */
public class FileReader
{

    private static Logger logger = LoggerFactory.getLogger(FileReader.class);

    private FileReader()
    {

    }

    public static DataTable readDataTable(String inFile)
    {
        DataTable tab = DataTable.newInstance();
        BufferedReader bR;
        try
        {
            bR = new BufferedReader(new java.io.FileReader(new File(inFile)));
            String line;
            while ((line = bR.readLine()) != null)
                    tab.rbind(line.split("\t"));
            bR.close();
        }
        catch (IOException e)
        {
            logger.error("I cant read this file!");
        }

        return tab;
    }

    public static HashMap<String, Double> readTSV(String file, HeaderType headType, String[] cols, String spl)
    {
        switch (headType)
        {
            case HeaderOneLine:
                return _readOneLineHeaderTSV(file, cols, spl);
            default:
                throw new IllegalArgumentException("Ay, where is the rub!");

        }
    }

    private static HashMap<String, Double> _readOneLineHeaderTSV(String file, String[] cols, String spl)
    {
        HashMap<String, Double> map = new HashMap<>();
        BufferedReader bR;
        try
        {
            bR = new BufferedReader(new java.io.FileReader(new File(file)));
            String line = bR.readLine();
            String headerToks[] = line.split(spl);
            int idx[] = new int[cols.length];
            for (int i = 0; i < headerToks.length; i++)
            {
                for (int j = 0; j < cols.length; j++)
                {
                    if (headerToks[i].equals(cols[j]))
                    {
                        idx[j] = i;
                        j = cols.length;
                    }
                }
            }
            while ((line = bR.readLine()) != null)
            {
                String[] toks = line.split("\t");
                map.put(toks[idx[0]], Double.parseDouble(toks[idx[1]]));
            }
            bR.close();
        }
        catch (IOException e)
        {
            logger.error("I cant read this file!");
        }
        return map;
    }
}
