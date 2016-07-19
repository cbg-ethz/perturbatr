package net.digital_alexandria.mvc.view;

import edu.uci.ics.jung.algorithms.layout.CircleLayout;
import edu.uci.ics.jung.visualization.BasicVisualizationServer;
import edu.uci.ics.jung.visualization.decorators.ToStringLabeller;
import net.digital_alexandria.mvc.controller.Controller;
import net.digital_alexandria.structs.graph.Edge;
import net.digital_alexandria.structs.graph.Graph;
import net.digital_alexandria.structs.graph.Node;
import org.apache.commons.collections15.Transformer;

import javax.swing.*;
import java.awt.*;

/**
 * @author Simon Dirmeier {@literal simon.dirmeier@gmx.de}
 */
public class GraphForm extends javax.swing.JFrame
{
    private Controller cntr;
    private static final int maxSteps = 1000;
    private static final int DIM = 600;

    private BasicVisualizationServer<Node, Edge> _vs;

    private GraphForm(String[] args)
    {
        this.cntr = new Controller(args);
        _display();
    }

    private GraphForm(Graph g)
    {
        this.cntr = new Controller(g);
        _display();
    }

    private void _display()
    {

        setDefaultCloseOperation(javax.swing.WindowConstants.EXIT_ON_CLOSE);
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        setTitle("Markov Random Walk");
        _set_graph_layout();
        pack();
        setVisible(true);
    }

    private void _set_graph_layout()
    {

        _vs = new BasicVisualizationServer<>
            (new CircleLayout<>(cntr.vizGraph()), new Dimension(DIM, DIM));
        _vs.getRenderContext().setVertexLabelTransformer(new ToStringLabeller<>());
        _vs.getRenderContext().setVertexFillPaintTransformer(integer -> Color.BLACK);
        if (cntr.graph().size() <= 10)
            _vs.getRenderContext().setEdgeLabelTransformer(e -> String.format("%.2f", e.weight()));
        getContentPane().add(_vs);
    }

    private void _recolor(Node<Integer> id)
    {
        try
        {
            Transformer<Node, Paint> vertexColor = i -> {
                if (i.id() == id.id()) return Color.GREEN;
                else if (i.isStart()) return Color.BLUE;
                else if (i.isFrequent()) return Color.CYAN;
                return Color.BLACK;
            };
            _vs.getRenderContext().setVertexFillPaintTransformer(vertexColor);
            try
            {
                invalidate();
                validate();
                repaint();
            }
            catch (Exception ee)
            {
                System.err.println("maaansss");
            }
        }
        catch (Exception e)
        {
            System.err.println("maaan");
        }

    }

    public static void run(String[] args)
    {
        try
        {
            for (javax.swing.UIManager.LookAndFeelInfo info : javax.swing.UIManager.getInstalledLookAndFeels())
            {
                if ("Nimbus".equals(info.getName()))
                {
                    javax.swing.UIManager.setLookAndFeel(info.getClassName());
                    break;
                }
            }
        }
        catch (ClassNotFoundException |
            InstantiationException |
            javax.swing.UnsupportedLookAndFeelException |
            IllegalAccessException ex)
        {
            java.util.logging.Logger.getLogger(GraphForm.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        }
        GraphForm g = new GraphForm(args);
        g._startWalk();
    }

    private void _startWalk()
    {

        int cnt = 0;
        while (cnt++ < maxSteps)
        {
            try
            {
                if (maxSteps <= 20)
                    Thread.sleep(3000);
                else
                    Thread.sleep(400);
            }
            catch (InterruptedException ignored)
            {
            }
            System.out.println("Doing step: " + cnt);
            _recolor(this.cntr.walk());
        }
    }

    public static void show(Graph g)
    {
        new GraphForm(g);
    }
}
