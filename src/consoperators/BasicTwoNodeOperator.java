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
public class BasicTwoNodeOperator extends TreeOperator {
    final public Input<RealParameter> rateInput = new Input<>("rates", "the rates associated with nodes in the tree for sampling of individual rates among branches.", Input.Validate.REQUIRED);
    public final Input<BranchRateModel.Base> branchRateModelInput = new Input<>("branchRateModel",
            "A model describing the rates on the branches of the beast.tree.");

    RealParameter rates;
    BranchRateModel.Base branchRateModel;


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
        //tree.getInternalNodes()
        do {
            final int nodeNr = nodeCount / 2 + 1 + Randomizer.nextInt(nodeCount / 2);
            node = tree.getNode(nodeNr);

        } while (node.isRoot() || node.isLeaf());//rule out the leaf and root



        //get the details of the child nodes of the selected node
        //get parent node
        Node C1 = node.getChild(0);//always work on this node!!!
        if (C1.isLeaf()){
            return Double.NEGATIVE_INFINITY;
        }
        Node C2 = node.getChild(1);
        Node GC1 = C1.getChild(0);
        Node GC2 = C1.getChild(1);

        //node time of the nodes
        double t_p = node.getHeight();
        double t_x = C1.getHeight();
        double t3 = C2.getHeight();
        double t1 = GC1.getHeight();
        double t2 = GC2.getHeight();

        // /grand parent node of C1
        Node P = node.getParent();
        double t_GP = P.getHeight();

        //Lower bound in the Uniform distribution for the node times
        double L = Math.max(t1,t2);
        double H;
        /*
        Step2: get the original rates
         */
        double r1 = branchRateModel.getRateForBranch(GC1);
        double r2 = branchRateModel.getRateForBranch(GC2);
        double r3 = branchRateModel.getRateForBranch(C1);
        double r4 = branchRateModel.getRateForBranch(C2);

        /*
        Step3: Propose node times and rates
        */

            double r5 = branchRateModel.getRateForBranch(node);
            //propose node times
            H = t_GP;
            u1 = Randomizer.uniform(L,H);
            u2 = Randomizer.uniform(L,H);
            if (u1 == u2){
                return Double.NEGATIVE_INFINITY;
            }
            t_x_ = Math.min(u1,u2);
            t_p_ = Math.max(u1,u2);
            C1.setHeight(t_x_);
            node.setHeight(t_p_);

            //propose new rates
            r1_ = r1 * (t_x - t1) / (t_x_ - t1);
            r2_ = r2 * (t_x - t2) / (t_x_ - t2);
            r3_ = r3 * (t_p - t_x) / (t_p_ - t_x_);
            r4_ = r4 * (t_p - t3) / (t_p_ - t3);
            double r5_ = r5 * (t_GP - t_p) / (t_GP - t_p_);
            rates.setValue(GC1.getNr(),r1_);
            rates.setValue(GC2.getNr(),r2_);
            rates.setValue(C1.getNr(),r3_);
            rates.setValue(C2.getNr(),r4_);
            rates.setValue(node.getNr(),r5_);

            //calculate determinant of Jacobian

            double nu = (t_x - t1)* (t_x - t2) * (t_p -t3) * (t_p - t_x) * (t_GP - t_p);
            double de = (t_x_ - t1)* (t_x_ - t2) * (t_p_ -t3) * (t_p_ - t_x) * (t_GP - t_p_);
            Det = nu / de;


        //return Hastings Ratio in log space
        return Math.log(Det);
    }
}

        /*
        Step3: Propose node times and rates
         */
        /*
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
        }
        */