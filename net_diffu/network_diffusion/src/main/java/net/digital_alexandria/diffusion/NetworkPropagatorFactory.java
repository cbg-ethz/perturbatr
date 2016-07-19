package net.digital_alexandria.diffusion;

import net.digital_alexandria.structs.graph.Graph;

/**
 * @author Simon Dirmeier {@literal simon.dirmeier@gmx.de}
 */
public class NetworkPropagatorFactory
{
    private NetworkPropagatorFactory(){ }

    public static MarkovRandomWalk<Integer, Integer> newMarkovRandomWalk(Graph<Integer, Integer> g)
    {
        return new MarkovRandomWalk<>(g);
    }
}
