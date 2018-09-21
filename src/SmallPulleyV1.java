package myoperator;

import beast.core.Description;
import beast.core.Input;
import beast.core.Operator;
import beast.core.parameter.RealParameter;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.tree.Tree;
import beast.math.distributions.LogNormalDistributionModel;
import beast.math.distributions.ParametricDistribution;
import beast.math.distributions.Uniform;
import beast.util.Randomizer;
import beast.evolution.tree.Node;
import org.apache.commons.math.MathException;
import beast.evolution.operators.TreeOperator;

import java.util.List;

@Description("Small pulley version1: Propose a new genetic distance")
public class SmallPulleyV1 extends TreeOperator {
    //public final Input<Tree> treeInput = new Input<>("tree", "the rooted time tree for the operator to work on");
    public final Input<Double> dwindowSizeInput =
            new Input<>("dwindowSize", "the size of the window in Big Pulley");
    public final Input<BranchRateModel.Base> branchRateModelInput = new Input<>("branchRateModel",
            "A model describing the rates on the branches of the beast.tree.");
    final public Input<RealParameter> rateInput = new Input<>("rates", "the rates associated with nodes in the tree for sampling of individual rates among branches.", Input.Validate.REQUIRED);


    private double dwindowSize;
    private Tree tree;
    private RealParameter rates;
    private double hastingsRatio;
    protected BranchRateModel.Base branchRateModel;


    @Override
    public void initAndValidate() {
        dwindowSize = dwindowSizeInput.get();
        tree = treeInput.get();
        branchRateModel = branchRateModelInput.get();
        rates = rateInput.get();

    }

    @Override
    public double proposal() {
        final Tree tree = treeInput.get(this);

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

        // Step1: get the root of the tree
        node = tree.getRoot();


        // original root time
        t_x = node.getHeight();//get the time of this node

        // son
        Node son = node.getChild(0);//get the left child of this node, i.e. son
        t_j = son.getHeight();//node time of son
        int nodeN02 = son.getNr();//node number of son
        // daughter
        Node daughter = node.getChild(1);//get the right child of this node, i.e. daughter
        t_k = daughter.getHeight();//node time of daughter
        int nodeN03 = daughter.getNr();//node number of daughter

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
        r_j_ = (D - d_) / (t_x - t_j);
        r_k_ = d_ / (t_x - t_k);

        //Step 4: set the proposed new rates
        rates.setValue(nodeN02, r_j_);
        rates.setValue(nodeN03, r_k_);

        return hastingsRatio = 0.0;
    }
}

   /*
                 * Strategy one: Simple Distance
                 */
//node.setHeight(t_x_);
//Step 3: make changes on the rates
//r_j_ = r_j * (t_x - t_j) / (t_x_ - t_j);
//r_k_ = r_k * (t_x - t_k) / (t_x_ - t_k);
                /*
                 * Strategy two: Small Pulley
                 */



