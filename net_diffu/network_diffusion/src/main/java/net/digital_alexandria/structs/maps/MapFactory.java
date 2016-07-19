package net.digital_alexandria.structs.maps;

import net.digital_alexandria.structs.enums.DataType;
import org.slf4j.LoggerFactory;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;


/**
 * @author Simon Dirmeier {@literal simon.dirmeier@gmx.de}
 */
public class MapFactory
{

    private final static org.slf4j.Logger _LOGGER = LoggerFactory.getLogger(MapFactory.class);

    private MapFactory()
    {
    }

    public static MultiMap<String, Integer> newMultiMap(DataType type, String file)
    {
        switch (type)
        {
            case DataTable:
                return newMappingFromDataTable(file);
            case FLAT:
                return newMappingFromFlatFile(file);
            default:
                return null;
        }
    }

    private static MultiMap<String, Integer> newMappingFromFlatFile(String flatFile)
    {
        MultiMap<String, Integer> map = new MultiMap<>();
        BufferedReader bR;
        try
        {
            bR = new BufferedReader(new FileReader(new File(flatFile)));
            String line;
            while ((line = bR.readLine()) != null)
            {
                if (line.startsWith("#")) continue;
                String[] toks = line.split("\t");
                String ens = toks[1].replace("9606.", "");
                int entrez = Integer.parseInt(toks[0]);
                map.put(ens, entrez);
            }
            bR.close();
        }
        catch (IOException e)
        {
            Logger.getLogger(MapFactory.class.getSimpleName())
                  .log(Level.SEVERE, "IO-error");
        }
        return map;
    }

    private static MultiMap<String, Integer> newMappingFromDataTable(String dataTableFile)
    {
        MultiMap<String, Integer> map = new MultiMap<>();
        BufferedReader bR;
        try
        {
            bR = new BufferedReader(new FileReader(new File(dataTableFile)));
            String line;
            while ((line = bR.readLine()) != null)
            {
                if (line.startsWith("Virus")) continue;
                String[] toks = line.split("\t");
                if (toks[5].equals("NA") || toks[15].equals("NA")) continue;
                String hugo = toks[5];
                int entrez = Integer.parseInt(toks[15]);
                map.put(hugo, entrez);
            }
            bR.close();
        }
        catch (IOException e)
        {
            _LOGGER.error("IO-error");
        }
        return map;
    }
}
