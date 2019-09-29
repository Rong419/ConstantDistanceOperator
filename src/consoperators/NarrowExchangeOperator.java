package consoperators;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.operators.TreeOperator;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import java.text.DecimalFormat;


@Description("Implements the NARROW variety in Exchange operator in beas2" +
             "but maintain the genetic distances by proposing branch rates")
public class NarrowExchangeOperator extends TreeOperator {
    public final  Input<RealParameter> rateInput = new Input<>("rates", "the rates associated with nodes in the tree for sampling of individual rates among branches.", Input.Validate.REQUIRED);
    public final Input<Double> rwindowSizeInput =
            new Input<>("rwindowSize", "the size of the window when proposing branch rate for node 'parent'", Input.Validate.REQUIRED);

    private RealParameter rates;
    private double rwindowSize;

    @Override
    public void initAndValidate() {
        rates = rateInput.get();
        rwindowSize = rwindowSizeInput.get();
    }

    @Override
    public double proposal() {
        final Tree tree = treeInput.get(this);

        // Step 1: randomly select a valid grandParent node
        // get the number of internal nodes
        final int internalNodes = tree.getInternalNodeCount();
        if (internalNodes <= 1) {
            return Double.NEGATIVE_INFINITY;
        }

        // child nodes of a valid grandParent node cannot be both leaf nodes
        Node grandParent = tree.getNode(internalNodes + 1 + Randomizer.nextInt(internalNodes));
        while (grandParent.getChild(0).isLeaf() && grandParent.getChild(1).isLeaf()) {
            grandParent = tree.getNode(internalNodes + 1 + Randomizer.nextInt(internalNodes));
        }

        // Step2: access to the child nodes of grandParent
        // the older child node is denoted by parent
        // the younger child nodes is denoted by uncle
        Node parent = grandParent.getChild(0);
        Node uncle = grandParent.getChild(1);
        if (parent.getHeight() < uncle.getHeight()) {
            parent = grandParent.getChild(1);
            uncle = grandParent.getChild(0);
        }

        // randomly select a child node of parent, denoted by child
        final Node child = (Randomizer.nextBoolean() ? parent.getChild(0) : parent.getChild(1));

        // get the node times
        double tp = parent.getHeight();
        double tc = child.getHeight();
        double tu = uncle.getHeight();
        double tgp = grandParent.getHeight();

        // avoid the parent is a node with tip date
        if( parent.isLeaf() ) {
            return Double.NEGATIVE_INFINITY;
        }

        // the probability of making this proposal
        // = the probability of picking the random grandParent
        // = 1 / (the number of all valid grandParent nodes in the current tree)
        int validGP = 0;

        for(int i = internalNodes + 1; i < 1 + 2*internalNodes; ++i) {
            validGP += isg(tree.getNode(i));
        }

        final int c2 = sisg(parent) + sisg(uncle);


        // rates in original state
        double ru = rates.getValue(uncle.getNr());
        double rc = rates.getValue(child.getNr());
        double rp = rates.getValue(parent.getNr());

        // Step3: exchange child and uncle, child.e. child <--> uncle
        exchangeNodes(child, uncle, parent, grandParent);

        // Step4: propose new rates
        double upper = rc * (tp - tc) / (tgp - tp);
        double rp_ = Randomizer.uniform(-rwindowSize,rwindowSize) + rp;

        if (rp_ <= 0 || rp_ >= upper) {
            return Double.NEGATIVE_INFINITY;
        }

        double ru_ = (ru * (tgp - tu) + rp * (tgp - tp)) / (tp - tu);

        double rc_ = (rc * (tp - tc) - rp_ * (tgp - tp)) / (tgp - tc);


        // set the new rates
        rates.setValue(child.getNr(),rc_);
        rates.setValue(uncle.getNr(),ru_);
        rates.setValue(parent.getNr(),rp_);



        // Step4: return the hastings ratio in log space
        // the probability of going back to the original state
        // = 1 / (the number of all valid grandParent nodes in the proposed tree)
        final int validGPafter = validGP - c2 + sisg(parent) + sisg(uncle);

        // the ratio of exchange child and uncle
        double hastingsRatio = (float)validGP/validGPafter;

        // the determinant of jacobian matrix
        double Jacobian = ((tp - tc) * (tgp - tu)) / ((tgp - tc) * (tp - tu));

        // the Green's corrected hastings ratio
        return Math.log(hastingsRatio * Jacobian);
    }

    private int isg(final Node n) {
        return (n.getLeft().isLeaf() && n.getRight().isLeaf()) ? 0 : 1;
    }

    private int sisg(final Node n) {
        return n.isLeaf() ? 0 : isg(n);
    }

    private void exchangeNodes(Node i, Node j,
                                 Node p, Node jP) {
        // precondition p -> i & jP -> j
        replace(p, i, j);
        replace(jP, j, i);
        // postcondition p -> j & p -> i
    }

    /*
     Tuning the parameter: rwindowsize
    */
    @Override
    public double getCoercableParameterValue() {
        return rwindowSize;
    }


    @Override
    public void setCoercableParameterValue(double value) {
        rwindowSize = value;
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

        delta += Math.log(rwindowSize);
        rwindowSize = Math.exp(delta);
    }

    @Override
    public final String getPerformanceSuggestion() {
        double prob = m_nNrAccepted / (m_nNrAccepted + m_nNrRejected + 0.0);
        double targetProb = getTargetAcceptanceProbability();

        double ratio = prob / targetProb;
        if (ratio > 2.0) ratio = 2.0;
        if (ratio < 0.5) ratio = 0.5;

        // new scale factor
        double newWindowSize = rwindowSize * ratio;

        DecimalFormat formatter = new DecimalFormat("#.###");
        if (prob < 0.10) {
            return "Try setting window size to about " + formatter.format(newWindowSize);
        } else if (prob > 0.40) {
            return "Try setting window size to about " + formatter.format(newWindowSize);
        } else return "";
    }

}
