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
public class WideTwoNodeOperator extends TreeOperator {
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
        //the tree to work with
        Tree tree = treeInput.get();
        //the internal node to operate on
        Node node;
        //hastings ratio to return
        double hastingsratio;

        /*
        Step1: randomly select one internal node in the tree
         */
        int nodeCount = tree.getNodeCount();
        //tree.getInternalNodes()
        do {
            final int nodeNr = nodeCount / 2 + 1 + Randomizer.nextInt(nodeCount / 2);
            node = tree.getNode(nodeNr);

        } while (node.isRoot() || node.isLeaf());//rule out the leaf and root

          /*
        Step2: get the parent node and child nodes of this node
         */
          Node P = node.getParent();
          Node C1 = node.getChild(0);
          Node C2 = node.getChild(1);

          /*
        Step3: get the sibling of this node
         */
          Node S;
          if (node == P.getChild(0)) {
              S = P.getChild(1);
          } else {
              S = P.getChild(0);
          }


         //Step 4: Propose new node times and rates
         /*
        Case1: P is the root of tree
               change 4 rates above C1,C2,S and node
                P
               /\
         node /  \
             /\   \
            /  \   \
           C1  C2  S
         */
         if (P.isRoot()) {
              hastingsratio = RateChange(C1,C2,S,node,P);
         }

         /*
        Case2: P is not the root
               firstly, change the 4 rates above C1, C2, S and node
               secondly, change the rate above P, i.e. r5 ---> r5'

                 GP
                 /\
              P /  \
               /\   \
         node /  \   \
             /\   \   \
            /  \   \   \
           C1  C2  S
         */
         else {
             double t_GP = P.getParent().getHeight();
             double r5 = branchRateModel.getRateForBranch(P);
             //original node time of P
             double t5 = P.getHeight();
             //propose four rates on branches above C1, C2, S and node
             double r = RateChange(C1,C2,S,node,P);
             //get the proposed node time of P
             double t5_ = P.getHeight();
             //To Propose a new rate for branch above P
             double r5_ = r5 * (t_GP - t5) / (t_GP - t5_);
             rates.setValue(P.getNr(),r5_);

             hastingsratio = r * (t_GP - t5) / (t_GP - t5_);

         }

        return Math.log(hastingsratio);
    }


    /*
    This method changes four rates on the branches above N1, N2, S, C
    In the meantime, two node times are proposed for C and P
                P
               /\
            C /  \
             /\   \
            /  \   \
           N1  N2  S
     */
    public double RateChange (Node N1, Node N2, Node S, Node C, Node P) {
        //original node times
        double t1 = N1.getHeight();
        double t2 = N2.getHeight();
        double t3 = S.getHeight();
        double tc = C.getHeight();
        double tp = P.getHeight();

        //original rates
        double r1 = branchRateModel.getRateForBranch(N1);
        double r2 = branchRateModel.getRateForBranch(N2);
        double r3 = branchRateModel.getRateForBranch(S);
        double r4 = branchRateModel.getRateForBranch(C);

        //propose new node times
        double T [] = NodeTime(C, P);
        double tc_ = T[0]; double tp_ = T[1];

        //propose new rates
        double r1_ = r1 * (tc - t1) / (tc_ - t1);
        double r2_ = r2 * (tc - t2) / (tc_ - t2);
        double r3_ = r3 * (tp - t3) / (tp_ - t3);
        double r4_ = r4 * (tp - tc) / (tp_- tc_);

        //update the new rates
        rates.setValue(N1.getNr(),r1_);
        rates.setValue(N2.getNr(),r2_);
        rates.setValue(S.getNr(),r3_);
        rates.setValue(C.getNr(),r4_);

        //the determinant of the Jacobian matrix
        double nu = (tc - t1)*(tc - t2)*(tp - t3)*(tp - tc);
        double de = (tc_ - t1)*(tc_ - t2)*(tp_ - t3)*(tp_ - tc_);

        return nu/de;
    }

    /*
    This method proposes two node times
    And also give the new values to the corresponding nodes

            A____________ t_A
            /\
        D  /__\__________ t_D
          /\   \
         /  \   \

    returns proposed times by Time = [t_D'  t_A']
     */
    public double[] NodeTime (Node D, Node A) {
        double upper;
        if (A.isRoot()) {
            upper = 10.0;
        } else {
            upper = A.getParent().getHeight();
        }
        double lower = Math.max(D.getChild(0).getHeight(),D.getChild(1).getHeight());

        double u1 = Randomizer.uniform(lower,upper);
        double u2 = Randomizer.uniform(lower,upper);
        double t_D = Math.min(u1, u2);
        double t_A = Math.max(u1, u2);

        //set the new node times
        D.setHeight(t_D);
        A.setHeight(t_A);

        //return the new node times for proposing new rates
        double Time []= {0.0,0.0};
        Time[0] =t_D; Time[1] =t_A;
        return Time;
    }


}
