package consoperators;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.Operator;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.tree.Tree;
import beast.base.inference.distribution.LogNormalDistributionModel;
import beast.base.inference.distribution.ParametricDistribution;
import beast.base.inference.distribution.Uniform;
import beast.base.util.Randomizer;
import beast.base.evolution.tree.Node;
import org.apache.commons.math.MathException;
import beast.base.evolution.operator.TreeOperator;

import java.text.DecimalFormat;
import java.util.List;

@Description("Small pulley version2: Propose a new genetic distance and a new root time")
public class SmallPulleyV2 extends TreeOperator {
    //public final Input<Tree> treeInput = new Input<>("tree", "the rooted time tree for the operator to work on");
    public final Input<Double> twindowSizeInput =
            new Input<>("twindowSize", "the size of the window for proposing new node time", Input.Validate.REQUIRED);
    public final Input<Double> dwindowSizeInput =
            new Input<>("dwindowSize", "the size of the window for proposing new genetic distance");
    public final Input<BranchRateModel.Base> branchRateModelInput = new Input<>("branchRateModel",
            "A model describing the rates on the branches of the beast.tree.");
    final public Input<RealParameter> rateInput = new Input<>("rates", "the rates associated with nodes in the tree for sampling of individual rates among branches.", Input.Validate.REQUIRED);

    private double twindowSize;
    private double dwindowSize;
    private Tree tree;
    private RealParameter rates;
    private double hastingsRatio;
    protected BranchRateModel.Base branchRateModel;

    @Override
    public void initAndValidate() {
        twindowSize = twindowSizeInput.get();
        dwindowSize = dwindowSizeInput.get();
        tree = treeInput.get();
        branchRateModel = branchRateModelInput.get();
        rates = rateInput.get();
    }

    @Override
    public double proposal() {
        final Tree tree = treeInput.get();

        //original rates
        double r_i;
        double r_j;
        double r_k;
        double r_x;
        //proposed rates
        double r_i_;
        double r_j_;
        double r_k_;
        double r_x_;
        //the chosen node to work on
        Node node;
        //nodes below this chosen node
        //Node node1 = new Node(); //the higher child node
        //Node node2 = new Node(); //the lower child node
        double a; //t factor in the operator, i.e. to move the node time
        //the original node time
        double t_x;
        double t_j;
        double t_k;
        double t_i;
        double t1 = 1.0;
        double t2 = 0.0;
        //the proposed node time
        double t_x_;

        // Step 1: get the root of the tree
        node = tree.getRoot();

        // this node
        t_x = node.getHeight();//get the time of this node
        //int nodeN01 = node.getNr();//node number of this node
        // son
        Node son = node.getChild(0);//get the left child of this node, i.e. son
        t_j = son.getHeight();//node time of son
        int nodeN02 = son.getNr();//node number of son
        // daughter
        Node daughter = node.getChild(1);//get the right child of this node, i.e. daughter
        t_k = daughter.getHeight();//node time of daughter
        int nodeN03 = daughter.getNr();//node number of daughter


        double lower = Math.max(t_j, t_k);

        // Step 2: generate a random number from a Uniform distribution
        //         to propose a new node time for this node
        a = Randomizer.uniform(-twindowSize, twindowSize);
        t_x_ = t_x + a;

        if (t_x_ <= lower) {
            return Double.NEGATIVE_INFINITY;
        }

        node.setHeight(t_x_);

        // get the rates on branches linked to root
        r_j = branchRateModel.getRateForBranch(son);//branch rate for son
        r_k = branchRateModel.getRateForBranch(daughter);//branch rate for daughter

        double d = r_k * (t_x - t_k);//d3, which should be logged
        double D = r_j * (t_x - t_j) + d;//d3+d4
        double b = Randomizer.uniform(-dwindowSize, dwindowSize);
        double d_ = d + b;
        if (d_ <= 0.0 || d_ >= D) {
            return Double.NEGATIVE_INFINITY;
        }

        //Step 3: make changes on the rates
        r_j_ = (D - d_) / (t_x_ - t_j);
        r_k_ = d_ / (t_x_ - t_k);

        //Step 4: set the proposed new rates
        rates.setValue(nodeN02, r_j_);
        rates.setValue(nodeN03, r_k_);

        return hastingsRatio = 0.0; //in log space
    }
    /**
     * automatic parameter tuning *
     */
    @Override
    public double getCoercableParameterValue() {
        return dwindowSize;
    }

    @Override
    public void setCoercableParameterValue(double value) {
        dwindowSize = value;
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

        delta += Math.log(dwindowSize);
        dwindowSize = Math.exp(delta);
    }

    @Override
    public final String getPerformanceSuggestion() {
        double prob = m_nNrAccepted / (m_nNrAccepted + m_nNrRejected + 0.0);
        double targetProb = getTargetAcceptanceProbability();

        double ratio = prob / targetProb;
        if (ratio > 2.0) ratio = 2.0;
        if (ratio < 0.5) ratio = 0.5;

        // new scale factor
        double newWindowSize = dwindowSize * ratio;

        DecimalFormat formatter = new DecimalFormat("#.###");
        if (prob < 0.10) {
            return "Try setting window size to about " + formatter.format(newWindowSize);
        } else if (prob > 0.40) {
            return "Try setting window size to about " + formatter.format(newWindowSize);
        } else return "";
    }
}





