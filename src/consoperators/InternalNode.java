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
import java.util.ArrayList;
import java.util.List;

@Description("For internal nodes: propose a new node time")
public class InternalNode extends TreeOperator {
    public final Input<Double> twindowSizeInput =
            new Input<>("twindowSize", "the size of the window when proposing new node time", Input.Validate.REQUIRED);
    final public Input<RealParameter> rateInput = new Input<>("rates", "the rates associated with nodes in the tree for sampling of individual rates among branches.", Input.Validate.REQUIRED);
    public final Input<BranchRateModel.Base> branchRateModelInput = new Input<>("branchRateModel",
            "A model describing the rates on the branches of the beast.tree.");

    private double twindowSize;
    private RealParameter rates;
    Tree tree;
    protected BranchRateModel.Base branchRateModel;

    @Override
    public void initAndValidate() {
        twindowSize = twindowSizeInput.get();
        rates = rateInput.get();
        tree = treeInput.get(this);
        branchRateModel = branchRateModelInput.get();
    }

    @Override
    public double proposal() {

        //the chosen node to work on
        Node node;
        //the original node times
        double t_x;
        double t_j;
        double t_k;
        //the original distances
        double d_i;
        double d_k;
        //the original rates
        double r_k;
        double r_i;

        //the proposed node time
        double t_x_;

        //Step 1: randomly select an internal node, denoted by node x

        List<Node> In_node = tree.getInternalNodes();

        int Count = In_node.size();

        //is the nodeCount equal to rateCount???
        int nodeCount = Randomizer.nextInt(Count);

        node = In_node.get(nodeCount);

        //node
        t_x = node.getHeight();//get the time of this node
        double r_node = rates.getValue(nodeCount);//branch rate for this node
        // son
        Node son = node.getChild(0);//get the left child of this node, i.e. son
        t_j = son.getHeight();//node time of son
        int sonCount = son.getNr();
        r_i = rates.getValue(sonCount);
        d_i = r_i * (t_x - t_j);
        // daughter
        Node daughter = node.getChild(1);//get the right child of this node, i.e. daughter
        t_k = daughter.getHeight();//node time of daughter
        int dauCount = daughter.getNr();
        r_k = rates.getValue(dauCount);
        d_k = r_k * (t_x - t_k);

        double a = Randomizer.uniform(-twindowSize, twindowSize);

        //Step2: to propose a new node time for this node
        t_x_ = t_x + a;

        double upper = node.getParent().getHeight();
        double lower = Math.max(t_j, t_k);

        if (t_x_<= lower || t_x_ >= upper) {
            return Double.NEGATIVE_INFINITY;
        }
        node.setHeight(t_x_);



       //Step3: propose the new rates
       //there are three rates in total
       //r_node, r_i, r_x
       double r_node_ = r_node * (upper - t_x) / (upper - t_x_);
       double r_i_ = d_i / (t_x_ - t_j);
       double r_k_ = d_k / (t_x_ - t_k);
       //if (r_i_ < 0.1 || r_k_  <0.1 || r_node_ <0.1) {
           //System.out.println("t_x_=" + t_x_);
           //System.out.println("r_node_ =" + r_node_ + ",r_i_=" + r_i_ + ",r_k_=" + r_k_);
           //System.out.println("d_node=" + r_node * (upper - t_x) + ",d_i=" + d_i + ",d_k=" + d_k);
           //return Double.NEGATIVE_INFINITY;
       //}
       // set the proposed new rates
       rates.setValue(node.getNr(), r_node_);
       rates.setValue(son.getNr(), r_i_);
       rates.setValue(daughter.getNr(), r_k_);

       return 0.0;
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

