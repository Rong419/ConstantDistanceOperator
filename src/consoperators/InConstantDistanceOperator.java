package consoperators;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.operators.TreeOperator;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;

import java.text.DecimalFormat;
import java.util.List;

@Description("For internal nodes: propose a new node time")
public class InConstantDistanceOperator extends TreeOperator {
    public final Input<Double> twindowSizeInput =
            new Input<>("twindowSize", "the size of the window when proposing new node time", Input.Validate.REQUIRED);
    public final Input<BranchRateModel.Base> branchRateModelInput = new Input<>("branchRateModel",
            "A model describing the rates on the branches of the beast.tree.");
    final public Input<RealParameter> rateInput = new Input<>("rates", "the rates associated with nodes in the tree for sampling of individual rates among branches.", Input.Validate.REQUIRED);

    private double twindowSize;
    private RealParameter rates;
    private double hastingsRatio;
    protected BranchRateModel.Base branchRateModel;


    @Override
    public void initAndValidate() {
        twindowSize = twindowSizeInput.get();
        branchRateModel = branchRateModelInput.get();
        rates = rateInput.get();
    }

    @Override
    public double proposal() {
        final Tree tree = treeInput.get(this);
        //the chosen node to work on
        Node node;

        //the original node time
        double t_x;
        double t_j;
        double t_k;

        //the proposed node time
        double t_x_;

        double r_x;
        double r_i;

        double d_i;
        double d_x;

        //Step 1: randomly select an internal node, denoted by node x
       int nodeCount = tree.getNodeCount();//return the number of nodes in the tree
        do {
            final int nodeNr = nodeCount / 2 + 1 + Randomizer.nextInt(nodeCount / 2);
            node = tree.getNode(nodeNr);
        } while (node.isRoot() || node.isLeaf());

        t_x = node.getHeight();//get the time of this node
        // son
        Node son = node.getChild(0);//get the left child of this node, i.e. son
        t_j = son.getHeight();//node time of son
        r_i = branchRateModel.getRateForBranch(son);
        d_i = r_i * (t_x - t_j);
        // daughter
        Node daughter = node.getChild(1);//get the right child of this node, i.e. daughter
        t_k = daughter.getHeight();//node time of daughter
        r_x = branchRateModel.getRateForBranch(daughter);
        d_x = r_x * (t_x - t_k);

        double a = Randomizer.uniform(-twindowSize, twindowSize);

        //to propose a new node time for this node
        t_x_ = t_x + a;

        double upper = node.getParent().getHeight();
        double lower = Math.max(t_j, t_k);

        if (t_x_<= lower || t_x_ >= upper) {
            return Double.NEGATIVE_INFINITY;
        }
        node.setHeight(t_x_);

        // there are three rates in total
       double r_node = branchRateModel.getRateForBranch(node);//branch rate for this node

       // propose the new rates
       double r_node_ = r_node * (upper - t_x) / (upper - t_x_);
       double r_i_ = d_i / (t_x_ - t_j);
       double r_x_ = d_x / (t_x_ - t_k);

       // set the proposed new rates
       rates.setValue(node.getNr(), r_node_);
       rates.setValue(son.getNr(), r_i_);
       rates.setValue(daughter.getNr(), r_x_);

       return hastingsRatio = 0.0;
    }


    @Override
    public double getCoercableParameterValue() {
        return twindowSize;
    }


    @Override
    public void setCoercableParameterValue(double value) {
        twindowSize = value;
    }

    /**
     * called after every invocation of this operator to see whether
     * a parameter can be optimised for better acceptance hence faster
     * mixing
     *
     * @param logAlpha difference in posterior between previous state & proposed state + hasting ratio
     */

    @Override
    public void optimize(double logAlpha) {
        // must be overridden by operator implementation to have an effect
        double delta = calcDelta(logAlpha);

        delta += Math.log(twindowSize);
        twindowSize = Math.exp(delta);
    }

    @Override
    public final String getPerformanceSuggestion() {
        double prob = m_nNrAccepted / (m_nNrAccepted + m_nNrRejected + 0.0);
        double targetProb = getTargetAcceptanceProbability();

        double ratio = prob / targetProb;
        if (ratio > 2.0) ratio = 2.0;
        if (ratio < 0.5) ratio = 0.5;

        // new scale factor
        double newWindowSize = twindowSize * ratio;

        DecimalFormat formatter = new DecimalFormat("#.###");
        if (prob < 0.10) {
            return "Try setting window size to about " + formatter.format(newWindowSize);
        } else if (prob > 0.40) {
            return "Try setting window size to about " + formatter.format(newWindowSize);
        } else return "";
    }


}

