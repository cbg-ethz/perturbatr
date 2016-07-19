package net.digital_alexandria.structs.graph;


import net.digital_alexandria.structs.maps.MultiMap;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.*;
import java.util.stream.Collectors;

/**
 * @author Simon Dirmeier {@literal simon.dirmeier@gmx.de}
 */
public class Graph<T extends Comparable, U extends Comparable> implements Iterable<Node<T>>
{

    private final static Logger _LOGGER = LoggerFactory.getLogger(Graph.class);

    private Map<T, Node<T>> _nodes;
    private Map<U, Edge<U>> _edges;

    Graph()
    {
        this._nodes = new HashMap<>();
        this._edges = new HashMap<>();
    }

    private void _normalizeEdgeWeights()
    {
        for (Node<T> n : this._nodes.values())
        {
            n.normalizeEdgeWeights();
        }
    }

    public int size()
    {
        return _nodes.size();
    }

    @Override
    public Iterator<Node<T>> iterator()
    {
        return this._nodes.values().iterator();
    }

    public Set<Edge<U>> edges()
    {
        return new HashSet<>(_edges.values());
    }

    public Set<Node<T>> nodes()
    {
        return new HashSet<>(_nodes.values());
    }

    public Node random()
    {
        return new ArrayList<>(_nodes.values())
            .get(new Random().nextInt(_nodes.values().size()));
    }

    void add(T key, Node<T> value)
    {
        if (!this._nodes.containsKey(key))
            this._nodes.put(key, value);
    }

    void add(U key, Edge<U> value)
    {
        value.from().addOutgoing(value);
        value.to().addOutgoing(value);
        this._edges.put(key, value);
    }

    public void normalizeEdgeWeights()
    {
        _normalizeEdgeWeights();
    }

    void put(T k, Node<T> n)
    {
        this._nodes.put(k, n);
    }

    void put(U k, Edge<U> n)
    {
        this._edges.put(k, n);
    }

    boolean contains(Edge<U> u)
    {
        return this._edges.containsValue(u);
    }

    public void setFrequentNodes()
    {
        ArrayList<Node<T>> nodes = new ArrayList<>(_nodes.values());
        Collections.sort(nodes, new Comparator<Node<T>>()
        {
            @Override
            public int compare(Node<T> o1, Node<T> o2)
            {
                if (o1.walkedBy() < o2.walkedBy())
                {
                    return 1;
                }
                else if (o1.walkedBy() > o2.walkedBy())
                {
                    return -1;
                }
                return 0;
            }
        });
        nodes.get(0).setFrequent(true);
        nodes.get(1).setFrequent(true);
        nodes.get(2).setFrequent(true);
        for (int i = 3; i < nodes.size(); i++)
        {
            nodes.get(i).setFrequent(false);
        }
    }

    Set<T> nodeIDs()
    {
        Set<T> ids = _nodes.values().stream().map(Node::id).collect(Collectors.toSet());
        return ids;
    }

    public void setGeneScores(MultiMap<String, Integer> entrezHugoMap,
                              HashMap<String, Double> scoreMap)
    {
        _setWeights(entrezHugoMap, scoreMap);
    }

    private void _setWeights(MultiMap<String, Integer> entrezHugoMap, HashMap<String, Double> scoreMap)
    {
        // sets scores for every node in the ggi
        for (Map.Entry<String, Double> e : scoreMap.entrySet())
        {
            String hugo = e.getKey();
            Set<Integer> entrezIDs = entrezHugoMap.getTo(hugo);
            if (entrezIDs.size() != 1)
                _LOGGER.warn("Multiple|No entrez IDs found for hugo symbol: " + hugo);
            int ent = new ArrayList<>(entrezIDs).get(0);
            double score = e.getValue();
            if (this._nodes.containsKey(ent))
            {
                this._nodes.get(ent).weight(score);

            }
            else
            {
                //  if (entrezIDs.size() != 0)
                // _LOGGER.warn("Could not get entrez node '" + ent + "' " +
                // sou   "to set gene effects!");
                //  else
                //  _LOGGER.warn("Could not get hugo node '" + hugo + "' " +
                //               "to set gene effects and entrez ID!");
            }
        }
    }

    private void _normalizeNodeWeights()
    {
        double sumOfEffects = 0.0;
        for (Node<T> n : this._nodes.values())
        {
            sumOfEffects += n.weight();
        }
        for (Node<T> n : this._nodes.values())
        {
            n.weight(n.weight() / sumOfEffects);
        }
    }

    Node<T> node(T k)
    {
        return this._nodes.get(k);
    }

    Node<T> getOrPut(T t)
    {
        if (this._nodes.containsKey(t))
            return this._nodes.get(t);
        else return new Node<>(t);
    }

    public void normalizeNodeWeights()
    {
        _normalizeNodeWeights();
    }
}
