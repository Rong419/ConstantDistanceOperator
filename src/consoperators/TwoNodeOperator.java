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
    Tree tree;

    @Override
    public void initAndValidate() {
        rates = rateInput.get();
        twindowSize = twindowSizeInput.get();
        branchRateModel = branchRateModelInput.get();
        tree = treeInput.get(this);
    }

    @Override
    public double proposal() {
        //the internal node to operate on
        Node node;
        //hastings ratio to return
        double hastingsratio;

        /*
        Step1: randomly select one internal node in the tree
         */
        int nodeCount = tree.getNodeCount();
        do {
            int nodeNr = nodeCount / 2 + 1 + Randomizer.nextInt(nodeCount / 2);
            node = tree.getNode(nodeNr);

        } while (node.isRoot() || node.isLeaf());//rule out the leaf and root

          /*
        Step2: get the parent node and child nodes of this node
         */
          Node P = node.getParent();
          if (P.isRoot()){
              return Double.NEGATIVE_INFINITY;
          }
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
          double ts = S.getHeight();

         /*
         Step 4: Propose new node times
         */
        double tx = node.getHeight();
        double tp = P.getHeight();
        double t1 = C1.getHeight();
        double t2 = C2.getHeight();

        double u1 = Randomizer.uniform(-twindowSize,twindowSize);
        double u2 = Randomizer.uniform(-twindowSize,twindowSize);
/*
        double a = tx + u1;
        double b = tp + u2;
        double tx_ = Math.min(a,b);
        double tp_ = Math.max(a,b);
*/

        double tx_ = tx + u1;
        double tp_ = tp + u2;
        if (tx_ >= tp_){
            return Double.NEGATIVE_INFINITY;
        }

        if (P.isRoot()) {
            if (tx_ < Math.max(t1, t2)) {
                return Double.NEGATIVE_INFINITY;
            } else {
                node.setHeight(tx_);
                P.setHeight(tp_);
            }
        } else {
            if (tx_ <= Math.max(t1, t2) || tp_ >= P.getParent().getHeight()) {
                return Double.NEGATIVE_INFINITY;
            } else {
                node.setHeight(tx_);
                P.setHeight(tp_);
            }
        }

         /*
         Step 5: Propose new rates
         */
        //original rates
        double r1 = branchRateModel.getRateForBranch(C1);
        double r2 = branchRateModel.getRateForBranch(C2);
        double r3 = branchRateModel.getRateForBranch(S);
        double r4 = branchRateModel.getRateForBranch(node);
        //propose new rates
        double r1_ = r1 * (tx - t1) / (tx_ - t1);
        double r2_ = r2 * (tx - t2) / (tx_ - t2);
        double r3_ = r3 * (tp - ts) / (tp_ - ts);
        double r4_ = r4 * (tp - tx) / (tp_- tx_);
        //update the new rates
        rates.setValue(C1.getNr(),r1_);
        rates.setValue(C2.getNr(),r2_);
        rates.setValue(S.getNr(),r3_);
        rates.setValue(node.getNr(),r4_);

        //the determinant of the Jacobian matrix
        double nu = (tx - t1)*(tx - t2)*(tp - ts)*(tp - tx);
        double de = (tx_ - t1)*(tx_ - t2)*(tp_ - ts)*(tp_ - tx_);
        double r = nu / de;

        if (!P.isRoot()) {
            double t_GP = P.getParent().getHeight();
            double r5 = branchRateModel.getRateForBranch(P);
            //To Propose a new rate for branch above P
            double r5_ = r5 * (t_GP - tp) / (t_GP - tp_);
            rates.setValue(P.getNr(),r5_);
            hastingsratio = r * (t_GP - tp) / (t_GP - tp_);
         } else {
            hastingsratio = r;
        }
        return Math.log(2*hastingsratio);
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
