package net.digital_alexandria.diffusion;

import net.digital_alexandria.structs.Hit;

import java.util.List;

/**
 * @author Simon Dirmeier {@literal simon.dirmeier@gmx.de}
 */
public interface NetworkPropagator<T extends Comparable>
{
    List<Hit<T>> diffuse();
}
