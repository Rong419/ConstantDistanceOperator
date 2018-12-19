package consoperators;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.operators.TreeOperator;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import beast.evolution.tree.Node;

@Description("For internal nodes: propose two node times and rates that are associated with the node times")
public class TwoNodeOperator extends TreeOperator {
    final public Input<RealParameter> rateInput = new Input<>("rates", "the rates associated with nodes in the tree for sampling of individual rates among branches.", Input.Validate.REQUIRED);
    public final Input<BranchRateModel.Base> branchRateModelInput = new Input<>("branchRateModel",
            "A model describing the rates on the branches of the beast.tree.");

    RealParameter rates;
    BranchRateModel.Base branchRateModel;
    JacobianMatrixDeterminant JD = new JacobianMatrixDeterminant();

    @Override
    public void initAndValidate() {
        rates = rateInput.get();
        branchRateModel = branchRateModelInput.get();
    }

    @Override
    public double proposal() {
        //the internal node to operate on
        Node node;
        //the tree to work with
        Tree tree = treeInput.get(this);
        //the proposed node times
        double t_x_; double t_p_;
        //the random number to propose new node times
        double u1;double u2;
        //the proposed rates
        double r1_;double r2_;double r3_;double r4_;
        //determinant
        double Det;

        /*
        Step1: randomly select one internal node in the tree
         */
        int nodeCount = tree.getNodeCount();
        do {
            final int nodeNr = nodeCount / 2 + 1 + Randomizer.nextInt(nodeCount / 2);
            node = tree.getNode(nodeNr);
        } while (node.isRoot() || node.isLeaf());
        //node time of the selected node
        double t_x = node.getHeight();
        //get the details of the child nodes of the selected node
        Node Ch1 = node.getChild(0);
        Node Ch2 = node.getChild(1);
        double t1 = Ch1.getHeight();
        double t2 = Ch2.getHeight();
        double t_O = Math.max(t1,t2);

        /*
        Step2: get the parent of the selected node
         */
        Node P = node.getParent();
        //node time of parental node
        double t_p = P.getHeight();

        //get the original rates
        double r1 = branchRateModel.getRateForBranch(Ch1);
        double r2 = branchRateModel.getRateForBranch(Ch2);
        double r3 = branchRateModel.getRateForBranch(P.getChild(0));
        double r4 = branchRateModel.getRateForBranch(P.getChild(1));


        /*
        Step3: Propose node times and rates
         */
        if (P.isRoot()) {
          //Case1: P is the root of the tree
            //propose node times
            u1 = Randomizer.uniform(t_O, Double.POSITIVE_INFINITY);
            u2 = Randomizer.uniform(t_O, Double.POSITIVE_INFINITY);
            t_x_ = Math.min(u1,u2);
            t_p_ = Math.max(u1,u2);
            node.setHeight(t_x);
            P.setHeight(t_p_);

            //propose new rates
            r1_ = r1 * (t_x - t1) / (t_x_ - t1);
            r2_ = r2 * (t_x - t2) / (t_x_ - t2);
            r3_ = r3 * (t_p - t_x) / (t_p_ - t_x_);
            r4_ = r4 * (t_p - P.getChild(1).getHeight()) / (t_p_ - P.getChild(1).getHeight());
            rates.setValue(Ch1.getNr(),r1_);
            rates.setValue(Ch2.getNr(),r2_);
            rates.setValue(node.getNr(),r3_);
            rates.setValue(P.getChild(1).getNr(),r4_);

            //calculate HastingsRatio
            double [][] J1 = new double[4][4];
            Det = JD.Determinant(J1,3);
        }
        else {
          //Case2: P is an internal node
            //get the grandparent of the selected node
            Node GP = P.getParent();
            double t_gp = GP.getHeight();

            double r5 = branchRateModel.getRateForBranch(P);
            //propose node times
            u1 = Randomizer.uniform(t_O,t_gp);
            u2 = Randomizer.uniform(t_O,t_gp);
            t_x_ = Math.min(u1,u2);
            t_p_ = Math.max(u1,u2);
            node.setHeight(t_x);
            P.setHeight(t_p_);

            //propose new rates
            r1_ = r1 * (t_x - t1) / (t_x_ - t1);
            r2_ = r2 * (t_x - t2) / (t_x_ - t2);
            r3_ = r3 * (t_p - t_x) / (t_p_ - t_x_);
            r4_ = r4 * (t_p - P.getChild(1).getHeight()) / (t_p_ - P.getChild(1).getHeight());
            double r5_ = r5 * (t_gp - t_p) / (t_gp - t_p);
            rates.setValue(Ch1.getNr(),r1_);
            rates.setValue(Ch2.getNr(),r2_);
            rates.setValue(node.getNr(),r3_);
            rates.setValue(P.getChild(1).getNr(),r4_);
            rates.setValue(GP.getNr(),r5_);

            //calculate HastingsRatio
            double [][] J2 = new double[5][5];
            Det = JD.Determinant(J2,4);
        }

        return Math.log(Det);
    }
}
