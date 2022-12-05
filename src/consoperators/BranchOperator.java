package consoperators;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.operator.TreeOperator;
import beast.base.evolution.tree.Tree;
import beast.base.util.Randomizer;
import beast.base.evolution.tree.Node;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.*;

@Description("For internal nodes: move the branch without changing the rate on it")
public class BranchOperator extends TreeOperator {
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
        Tree tree = treeInput.get();

        //hastings ratio to return
        double hastingsratio;

        /*
        Step1: select the internal node in the tree
         */
        Node node = GetInternalNode(tree);

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
               change 3 rates above C1,C2 and S
                P
               /\
         node /  \
             /\   \
            /  \   \
           C1  C2  S
         */
         if (P.isRoot()) {
              hastingsratio = RateChangeV2(C1,C2,S,node,P);
         }

         /*
        Case2: P is not the root
               firstly, change the 3 rates above C1, C2 and S
               secondly, change the rate above P, i.e. r4 ---> r4'

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
             double r4 = branchRateModel.getRateForBranch(P);
             //original node time of P
             double t4 = P.getHeight();
             //propose three rates on branches above C1, C2, S
             double r = RateChangeV2(C1,C2,S,node,P);
             //get the proposed node time of P
             double t4_ = P.getHeight();
             //To Propose a new rate for branch above P
             double r4_ = r4 * (t_GP - t4) / (t_GP - t4_);
             rates.setValue(P.getNr(),r4_);

             hastingsratio = r * (t_GP - t4) / (t_GP - t4_);
         }

        return Math.log(hastingsratio);
        //return 0.0;
    }

    /*
    This method is used to assign each internal node a weight
    that is the minimum space for the branch to move
     */
    public Node GetInternalNode (Tree tree){
        Node Nx;
        int nodeCount = tree.getNodeCount();
        //tree.getInternalNodes()
        do {
            final int nodeNr = nodeCount / 2 + 1 + Randomizer.nextInt(nodeCount / 2);
            Nx = tree.getNode(nodeNr);

        } while (Nx.isRoot() || Nx.isLeaf());//rule out the leaf and root

        /*
        List<Double> W = new ArrayList<>(); //the list of weight for each internal node
        List<Integer> N = new ArrayList<>(); //the list of node number in the tree
        for (Node I : tree.getInternalNodes()) {
            int Nr = I.getNr();
            double t = I.getHeight();
            double T = I.getParent().getHeight();

            //node time of the grand parent node
            double tG;
            if (I.getParent().isRoot()) {
                tG = Double.POSITIVE_INFINITY;
            } else {
                tG = I.getParent().getParent().getHeight();
            }

            //node time of the older child
            double tc = Math.max(I.getChild(0).getHeight(),I.getChild(1).getHeight());

            //the minimum space for the internal node and its parent node to move
            double space = Math.min(tG - T, t - tc);
            W.add(space);
            N.add(Nr);

        }
        //to select the internal node
        double wx = Collections.max(W);
        int nx = N.get(W.indexOf(wx));
        Nx= tree.getNode(nx);
        */
        return Nx;
    }

    /*
    This method changes three rates on the branches above N1, N2, S
    In the meantime, two node times are proposed for C and P
                P
               /\
            C /  \
             /\   \
            /  \   \
           N1  N2  S
     */
    public double RateChangeV2 (Node N1, Node N2, Node S, Node C, Node P) {
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

        //propose new node times
        double T [] = NodeTime(C, P);
        double tc_ = T[0]; double tp_ = T[1];

        //set the new node times
        if (P.isRoot()) {
            if (tc_ <= Math.max(t1, t2)) {
                return Double.NEGATIVE_INFINITY;
            } else {
                C.setHeight(tc_);
                P.setHeight(tp_);
            }
        } else {
            if (tc_ <= Math.max(t1, t2) || tp_ >= P.getParent().getHeight()) {
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

        //update the new rates
        rates.setValue(N1.getNr(),r1_);
        rates.setValue(N2.getNr(),r2_);
        rates.setValue(S.getNr(),r3_);

        //the determinant of the Jacobian matrix
        double nu = (tc - t1)*(tc - t2)*(tp - t3);
        double de = (tc_ - t1)*(tc_ - t2)*(tp_ - t3);

        return nu/de;
    }

    /*
    This method proposes two node times
    Use the same random number to propose the node times
    Which makes the difference between there two node times constant
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
        double u = Randomizer.uniform(-twindowSize,twindowSize);

        double t_D = D.getHeight() + u;
        double t_A = A.getHeight() + u;

        //do{} while()
        //to deal with invalid proposals
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
