package consoperators;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.operator.TreeOperator;
import beast.base.evolution.tree.Tree;
import beast.base.util.Randomizer;
import beast.base.evolution.tree.Node;

@Description("For internal nodes: propose two node times and rates that are associated with the node times")
public class SpecialTwoNodeOperator extends TreeOperator {
    final public Input<RealParameter> rateInput = new Input<>("rates", "the rates associated with nodes in the tree for sampling of individual rates among branches.", Input.Validate.REQUIRED);
    public final Input<BranchRateModel.Base> branchRateModelInput = new Input<>("branchRateModel",
            "A model describing the rates on the branches of the beast.tree.");

    RealParameter rates;
    BranchRateModel.Base branchRateModel;
    Node node_P; Node node_C;

    @Override
    public void initAndValidate() {
        rates = rateInput.get();
        branchRateModel = branchRateModelInput.get();
    }

    @Override
    public double proposal() {
        //the tree to work with
        Tree tree = treeInput.get();
        //the internal node to operate on
        Node node;
        //hastings ratio to return
        double hastingsratio;

        /*
        Step1: randomly select one node in the tree
         */
        int nodeCount = tree.getNodeCount();
        int nodeNr = nodeCount / 2 + 1 + Randomizer.nextInt(nodeCount / 2);
        node = tree.getNode(nodeNr);

        /*
        Case1: A leaf node
              get the grandparent node
              P: The left child node of grandparent
              C: The left child of parent
              propose two node times, i.e. P and C
              propose five rates
         */
        if (node.isLeaf()){
           node_P = node.getParent().getParent();
           node_C = node_P.getChild(0);

            hastingsratio = FiveRates(node_C,node_P);
        }
        /*
        Case2: Root of the tree
               get the left child node (can not be a leaf node)
               propose two node time, i.e. root and left child
               propose four rates
         */
        else if (node.isRoot()) {
             node_P = node;
             node_C = node.getChild(0);

             if (node_C.isLeaf()) {
                 return Double.NEGATIVE_INFINITY;
             }

             hastingsratio = FourRates(node_C,node_P);
        }
        /*
        Case3: An internal node
               get the left child node (can't be a leaf node)
               propose two node times, i.e. the internal node and its left child
               propose five rates
         */
        else {
            node_P = node;
            node_C = node.getChild(0);

            if (node_C.isLeaf()) {
                return Double.NEGATIVE_INFINITY;
            }

            hastingsratio = FiveRates(node_C, node_P);

        }
        return Math.log(hastingsratio);
    }

    /*
    N1: the child node
    N2: the parent node
    N1 is the left child node of N2
    t_N1: the new node time for N1
    t_N2: the new node time for N2
    this method proposes rates for s1, s2, s3 and N1
               N2
               /\
           N1 /  \
             /\   \
            /  \   \
           s1  s2  s3
     */
    public double BasicRateChange (Node N1, Node N2, double t_N1, double t_N2) {
        //child nodes of N1
        Node s1 = N1.getChild(0);
        Node s2 = N1.getChild(1);
        //right child node of N2
        Node s3 = N2.getChild(1);

        //original node times
        double TN1 = N1.getHeight();
        double TN2 = N2.getHeight();
        double t_s1 = s1.getHeight();
        double t_s2 = s2.getHeight();
        double t_s3 = s3.getHeight();
        //original rates
        double r1 = branchRateModel.getRateForBranch(s1);
        double r2 = branchRateModel.getRateForBranch(s2);
        double r3 = branchRateModel.getRateForBranch(s3);
        double r4 = branchRateModel.getRateForBranch(N1);

        //propose new rates
        double r1_ = r1 * (TN1 - t_s1) / (t_N1 - t_s1);
        double r2_ = r2 * (TN1 - t_s2) / (t_N1 - t_s2);
        double r3_ = r3 * (TN2 - t_s3) / (t_N2 - t_s3);
        double r4_ = r4 * (TN2 - TN1) / (t_N2 - t_N1);
        //update the new rates
        rates.setValue(s1.getNr(),r1_);
        rates.setValue(s2.getNr(),r2_);
        rates.setValue(s3.getNr(),r3_);
        rates.setValue(N1.getNr(),r4_);

        //the determinant of the Jacobian matrix
        double nu = (TN1 - t_s1)*(TN1 - t_s2)*(TN2 - t_s3)*(TN2 - TN1);
        double de = (t_N1 - t_s1)*(t_N1 - t_s2)*(t_N2 - t_s3)*(t_N2 - t_N1);

        return nu/de;
    }

    /*
    The method of changing four rates when the root is involved
    N1: child node
    N2: parent node
    N1 is the left child of N2

    The method proposes four rates for s1,s2,s3 and N1
     */
    public double FourRates (Node N1, Node N2) {
        //child nodes of N1
        Node s1 = N1.getChild(0);
        Node s2 = N1.getChild(1);

        //original node times
        double t_s1 = s1.getHeight();
        double t_s2 = s2.getHeight();

        double H1 = 10.0; //upper bound: for the root time, the upper bound is positive infinity theoretically
        double L1 = Math.max(t_s1,t_s2);//lower bound
        //propose new node times
        double T [] = TwoTime(L1, H1);
        double t_N1 = T[0]; double t_N2 = T[1];

        //propose four rates
        double ratio = BasicRateChange(N1,N2,t_N1,t_N2);

        //set new node times to N1 and N2
        N1.setHeight(t_N1);
        N2.setHeight(t_N2);

        return ratio;
    }

    /*
    N1: the child node
    N2: the parent node
    N1 is the left child node of N2
    The method proposes five rates for s1,s2,s3, N1 and N2
     */
    public double FiveRates (Node N1, Node N2) {
        //Here, getting the node s1,s2,s4 and their node times is used to
        //propose new node times for N1 and N2
        //because it needs boundaries for Uniform distribution

        //child nodes of N1
        Node s1 = N1.getChild(0);
        Node s2 = N1.getChild(1);
        //parent node of N2, i.e. grandparent of N1
        Node s4 = N2.getParent();

        //original node times
        double t_s1 = s1.getHeight();
        double t_s2 = s2.getHeight();
        double TN2 = N2.getHeight();
        double H2 = s4.getHeight(); //upper bound for new times

        double L2 = Math.max(t_s1,t_s2);//lower bound for new times
        //propose new node times for N1 and N2
        //denoted by t_N1 and t_N2
        double T [] = TwoTime(L2, H2);
        double t_N1 = T[0]; double t_N2 = T[1];

        //change the four rates
        double ratio = BasicRateChange(N1,N2,t_N1,t_N2);

        //set the new node times for node N1 and N2
        //must after BasicRateChange
        //because original rates will be required in the method
        N1.setHeight(t_N1);
        N2.setHeight(t_N2);

        //propose additional rate above N2
        double r5 = branchRateModel.getRateForBranch(N2);
        double r5_ = r5 * (H2 - TN2) / (H2 - t_N2);
        rates.setValue(N2.getNr(),r5_);

        return ratio * (H2 - TN2) / (H2 - t_N2);
    }
    /*
     This method generates two new node times from Uniform distribution
     (1) get two random numbers from the Uniform distribution
     (2) save the smaller one as t_C, i.e. the first element in Time
     (3) save the larger one as t_P, i.e. the second element in Time
     */
    public double[] TwoTime (double lower, double upper) {
        double u1 = Randomizer.uniform(lower,upper);
        double u2 = Randomizer.uniform(lower,upper);
        double t_C = Math.min(u1, u2);
        double t_P = Math.max(u1, u2);
        double Time []= {0.0,0.0};
        Time[0] =t_C; Time[1] =t_P;
        return Time;
    }


}
