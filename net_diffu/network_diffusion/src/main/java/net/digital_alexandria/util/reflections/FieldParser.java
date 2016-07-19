package net.digital_alexandria.util.reflections;

import net.digital_alexandria.util.annotations.ParamAnnotation;

import java.lang.reflect.Field;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * @author Simon Dirmeier {@literal simon.dirmeier@gmx.de}
 */
public class FieldParser
{

    private FieldParser() {}

    /**
     * Parses command line arguments and sets the field values of object
     *
     * @param o
     * @param c
     * @param p
     */
    public static void parse(Object o,
                             Class c,
                             net.digital_alexandria.param.ParamsParser p)
    {
        for (Field f : c.getDeclaredFields())
        {
            f.setAccessible(true);
            if (f.isAnnotationPresent(ParamAnnotation.class))
            {

                ParamAnnotation ann = f.getAnnotation(ParamAnnotation.class);
                try
                {
                    f.set(o, p.getArgument(ann.value()));
                }
                catch (IllegalAccessException e)
                {
                    Logger.getLogger(c.getSimpleName())
                          .log(Level.SEVERE,
                               "Reflections error parsing argument: " + ann.value());
                }
            }
            f.setAccessible(false);
        }
    }
}
