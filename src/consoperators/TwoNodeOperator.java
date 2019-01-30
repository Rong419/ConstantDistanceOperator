package consoperators;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.operators.TreeOperator;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import beast.evolution.tree.Node;

import java.text.DecimalFormat;

@Description("For internal nodes: propose two node times and rates that are associated with the node times")
public class TwoNodeOperator extends TreeOperator {
    final public Input<RealParameter> rateInput = new Input<>("rates", "the rates associated with nodes in the tree for sampling of individual rates among branches.", Input.Validate.REQUIRED);
    public final Input<Double> twindowSizeInput =
            new Input<>("twindowSize", "the size of the window when proposing new node time", Input.Validate.REQUIRED);
    public final Input<BranchRateModel.Base> branchRateModelInput = new Input<>("branchRateModel",
            "A model describing the rates on the branches of the beast.tree.");

    RealParameter rates;
    double twindowSize;
    BranchRateModel.Base branchRateModel;

    @Override
    public void initAndValidate() {
        rates = rateInput.get();
        twindowSize = twindowSizeInput.get();
        branchRateModel = branchRateModelInput.get();

    }

    @Override
    public double proposal() {
        //the tree to work with
        Tree tree = treeInput.get(this);
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

        //set the new node times
        if (P.isRoot()) {
            if (tc_ < Math.max(t1, t2)) {
                return Double.NEGATIVE_INFINITY;
            } else {
                C.setHeight(tc_);
                P.setHeight(tp_);
            }
        } else {
            if (tc_ < Math.max(t1, t2) || tp_ > P.getParent().getHeight()) {
                return Double.NEGATIVE_INFINITY;
            } else {
                C.setHeight(tc_);
                P.setHeight(tp_);
            }
        }

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

        //generate two random variables from Uniform distribution

        double u1 = Randomizer.uniform(-twindowSize,twindowSize);
        double u2 = Randomizer.uniform(-twindowSize,twindowSize);

        double a = D.getHeight() + u1;
        double b = A.getHeight() + u2;

        double t_D = Math.min(a,b);
        double t_A = Math.max(a,b);

        //set the new node times
        //D.setHeight(t_D);
        //A.setHeight(t_A);

        //return the new node times for proposing new rates
        double Time []= {0.0,0.0};
        Time[0] =t_D; Time[1] =t_A;
        return Time;
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
