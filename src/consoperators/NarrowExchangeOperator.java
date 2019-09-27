package consoperators;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.operators.TreeOperator;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import com.sun.org.apache.xerces.internal.util.SynchronizedSymbolTable;

import java.text.DecimalFormat;
import java.util.Arrays;

@Description("Implements the NARROW variety in Exchange operator in beas2" +
             "but maintain the genetic distances by proposing branch rates")
public class NarrowExchangeOperator extends TreeOperator {
    final public Input<RealParameter> rateInput = new Input<>("rates", "the rates associated with nodes in the tree for sampling of individual rates among branches.", Input.Validate.REQUIRED);

    private RealParameter rates;

    @Override
    public void initAndValidate() {
        rates = rateInput.get();
    }

    @Override
    public double proposal() {
        final Tree tree = treeInput.get(this);


        // get the number of internal nodes
        final int internalNodes = tree.getInternalNodeCount();
        if (internalNodes <= 1) {
            return Double.NEGATIVE_INFINITY;
        }

        // randomly select a valid grandParent node
        // whose child nodes cannot be both leaf nodes
        Node grandParent = tree.getNode(internalNodes + 1 + Randomizer.nextInt(internalNodes));
        while (grandParent.getChild(0).isLeaf() && grandParent.getChild(1).isLeaf()) {
            grandParent = tree.getNode(internalNodes + 1 + Randomizer.nextInt(internalNodes));
        }

        // access to the child nodes of grandParent
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
        Double [] r = rates.getValues();
        double ru = rates.getValue(uncle.getNr());
        double rc = rates.getValue(child.getNr());
        double rp = rates.getValue(parent.getNr());

        // exchange child and uncle, child.e. child <--> uncle
        exchangeNodes(child, uncle, parent, grandParent);

        // propose new rates
        double upper = rc * (tp - tc) / (tgp - tp);
        double rp_ = Randomizer.uniform(0,upper);

        if (rp_ <= 0 || rp_ >= upper) {
            return Double.NEGATIVE_INFINITY;
        }

        double ru_ = (ru * (tgp - tu) + rp * (tgp - tp)) / (tp - tu);

        double rc_ = (rc * (tp - tc) - rp_ * (tgp - tp)) / (tgp - tc);


        // set the new rates
        rates.setValue(child.getNr(),rc_);
        rates.setValue(uncle.getNr(),ru_);
        rates.setValue(parent.getNr(),rp_);

        // the probability of going back to the original state
        // = 1 / (the number of all valid grandParent nodes in the proposed tree)
        final int validGPafter = validGP - c2 + sisg(parent) + sisg(uncle);

        // return the hastings ratio in log space
        double R1 = (float)validGP/validGPafter;
        double R2 = (rc * rp * (tp - tc)) / (rp_ * ru_ * (tp - tu));

        System.out.println("r = "+ Arrays.toString(r));
        System.out.println("rc="+rc+",ru="+ru+",rp="+rp);
        System.out.println("upper="+upper+",rp_="+rp_);
        System.out.println("rc'="+rc_+",ru'="+ru_+",rp'="+rp_);
        System.out.println("r'=" + Arrays.toString(rates.getValues()));
        System.out.println("R2="+R2);

        //System.out.println("R1="+R1+",R2="+R2);
        double J = ((tp - tc) * (tgp - tu)) / ((tgp - tc) * (tp - tu));
        return Math.log(R1*R2*J);
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
}
