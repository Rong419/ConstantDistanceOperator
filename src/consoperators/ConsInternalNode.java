package consoperators;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.operator.TreeOperator;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.util.Randomizer;
import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.BetaDistribution;
import org.apache.commons.math.distribution.BetaDistributionImpl;

import java.text.DecimalFormat;

@Description("For internal nodes: propose a new node time")
public class ConsInternalNode extends TreeOperator {
    public final Input<Double> twindowSizeInput =
            new Input<>("twindowSize", "the size of the window when proposing new node time");
    final public Input<RealParameter> rateInput =
            new Input<>("rates", "the rates associated with nodes in the tree for sampling of individual rates among branches.", Input.Validate.REQUIRED);
    final public Input<Boolean> useRandomInput =
            new Input<>("random", "use random walk proposal for node time.");
    final public Input<Boolean> useUniformInput =
            new Input<>("uniform", "use uniform proposal for node time.");
    final public Input<Boolean> useBactrianInput =
            new Input<>("bactrian", "use bactrian proposal for node time.");
    final public Input<Boolean> useBetaInput =
            new Input<>("beta", "use beta proposal for node time.");

    private double twindowSize;
    private RealParameter rates;

    enum Style {
        RANDOM_WALK, UNIFORM, BACTRIAN, BETA
    }

    private Boolean RANDOM_WALK;
    private Boolean UNIFORM;
    private Boolean BACTRIAN;
    private Boolean BETA;
    Style style = Style.RANDOM_WALK;

    @Override
    public void initAndValidate() {
        rates = rateInput.get();
        UNIFORM = useUniformInput.get();
        RANDOM_WALK = useRandomInput.get();
        BACTRIAN = useBactrianInput.get();
        BETA = useBetaInput.get();
        twindowSize = twindowSizeInput.get();

        if(RANDOM_WALK == null) {
            if (UNIFORM != null) {
                style = Style.UNIFORM;
            } else if (BACTRIAN != null) {
                style = Style.BACTRIAN;
            } else {
                style = Style.BETA;
            }
        } else{
                style = Style.RANDOM_WALK;
            }
        }


    private BetaDistribution betaDistr = new BetaDistributionImpl(1,1);
    private BetaDistribution new_betaDistr = new BetaDistributionImpl(1,1);

    @Override
    public double proposal() {
        final Tree tree = treeInput.get();
        int nodeCount = tree.getNodeCount(); //return the number of nodes in the tree
        int branchCount = nodeCount - 1; //the number of branches of the tree


        //the chosen node to work on
        Node node;

        //the original node times
        double t_x;
        double t_j;
        double t_k;
        //the original distances
        double d_j;
        double d_k;
        //the original rates
        double r_k;
        double r_j;

        double hastingsRatio = 0.0;

        //the proposed node time
        double t_x_ = 1.0;

        //Step 1: randomly select an internal node, denoted by node x.
        // Avoid fake nodes used to connect direct ancestors into tree.
        do {
            final int nodeNr = nodeCount / 2 + 1 + Randomizer.nextInt(nodeCount / 2);
            node = tree.getNode(nodeNr);
        } while (node.isRoot() || node.isLeaf() || node.isFake());

        // the number of this node
        int nodeNr = node.getNr();
        // if this node has max number, then use the free index stored in root node to get rate.
        if (nodeNr == branchCount) {
            nodeNr = node.getTree().getRoot().getNr();
        }

        //rate and time for this node
        t_x = node.getHeight();
        double r_node = rates.getValues()[nodeNr];

        //Step 2: Access to the child nodes of this node
        // son
        Node son = node.getChild(0);//get the left child of this node, i.e. son
        t_j = son.getHeight();//node time of son

        int sonNr = son.getNr();// node number of son
        if (sonNr == branchCount) {
            sonNr = son.getTree().getRoot().getNr();
        }

        r_j = rates.getValues()[sonNr]; // rate of branch above son
        d_j = r_j * (t_x - t_j); // distance of branch above son


        // daughter
        Node daughter = node.getChild(1);//get the right child of this node, i.e. daughter
        t_k = daughter.getHeight();//node time of daughter

        int dauNr = daughter.getNr(); // node time of daughter
        if (dauNr == branchCount) {
            dauNr = daughter.getTree().getRoot().getNr();
        }

        r_k = rates.getValues()[dauNr];// rate of branch above daughter
        d_k = r_k * (t_x - t_k);// distance of branch above daughter


        //Step3: to propose a new node time
        double upper = node.getParent().getHeight();
        double lower = Math.max(t_j, t_k);
        switch (style) {
            case RANDOM_WALK: {
                t_x_ = t_x + Randomizer.uniform(-twindowSize, twindowSize);
            }
            break;

            case UNIFORM: {
                    t_x_ = Randomizer.uniform(lower, upper);
            }
            break;

            case BACTRIAN: {
                Node oldChild;
                if (t_j >= t_k) {
                    oldChild = son;
                } else {
                    oldChild = daughter;
                }

                double branchLength = node.getLength() + oldChild.getLength();
                double mu1 = t_x + (branchLength * twindowSize);
                double mu2 = t_x - (branchLength * twindowSize);
                double sigma = branchLength * 0.5 * twindowSize;

                // index of going up or going down
                double bactrianIndex = Randomizer.uniform(0, 1);

                if (bactrianIndex <= 0.5 ) {
                    // half of probability to go up
                    t_x_ = sigma * Randomizer.nextGaussian() + mu1;
                } else {
                    // half of probability to go down
                    t_x_ = sigma * Randomizer.nextGaussian() + mu2;
                }

                //t_x_ = 0.5 * (sigma * Randomizer.nextGaussian() + mu1) + 0.5 * (sigma * Randomizer.nextGaussian() + mu2);
            }
            break;

            case BETA: {
                double oldChildHeight = Math.max(t_j, t_k);
                double m = (t_x - oldChildHeight) / (upper - oldChildHeight);
                double aBeta = twindowSize * m + 1.0;
                double bBeta = twindowSize * (1.0 - m) + 1.0;
                betaDistr.setAlpha(aBeta);
                betaDistr.setBeta(bBeta);
                double new_m = 0;
                try {
                    new_m = betaDistr.inverseCumulativeProbability(Randomizer.nextDouble());
                } catch (MathException e) {
                    e.printStackTrace();
                }
                t_x_ = (upper - oldChildHeight) * new_m + oldChildHeight;

                // hastings ratio
                double forward = betaDistr.density(new_m);
                double new_aBeta = twindowSize * new_m + 1.0;
                double new_bBeta = twindowSize * (1.0 - new_m) + 1.0;
                new_betaDistr.setAlpha(new_aBeta);
                new_betaDistr.setBeta(new_bBeta);
                double backward = new_betaDistr.density(m);
                hastingsRatio = Math.log(backward / forward);
            }
        }

        //reject the proposal if meets the boundary
        if (t_x_ == lower || t_x_ == upper) {
            return Double.NEGATIVE_INFINITY;
        }
        /*
        do {
            if (t_x_ <= lower) {
                t_x_ = 2 * lower - t_x_;
            }
            if (t_x_ >= upper) {
            t_x_ = 2 * upper - t_x_;
            }
        } while (t_x_<= lower || t_x_ >= upper);
        */
        double e; double n; double r;
        if (t_x_ > upper) {
            e = t_x_ - upper;
            n = Math.floor(e / (upper - lower));
            r = e - n * (upper - lower);
            if (n % 2 == 0) {
                t_x_ = upper - r;
            } else {
                t_x_ = lower + r;
            }
        } else if (t_x_ < lower) {
            e = lower - t_x_;
            n = Math.floor(e / (upper - lower));
            r = e - n * (upper - lower);
            if (n % 2 == 0) {
                t_x_ = lower + r;
            } else {
                t_x_ = upper - r;
            }
        }

        node.setHeight(t_x_);


       //Step4: propose the new rates
       //there are three rates in total
       //r_node, r_i, r_x
       double r_node_ = r_node * (upper - t_x) / (upper - t_x_);
       double r_j_ = d_j / (t_x_ - t_j);
       double r_k_ = d_k / (t_x_ - t_k);

       // set the proposed new rates
       rates.setValue(nodeNr, r_node_);
       rates.setValue(sonNr, r_j_);
       rates.setValue(dauNr, r_k_);


       //Step4: calculate the Hastings ratio
        /*
         *input:t_x,r_node,r_j,r_k
         *
         *f:the function
         *f(t_x,r_node,r_j,r_k)
         *
         *output:t_x_,r_node_,r_j_,r_k_
         *t_x_ = t_x + a
         *r_j_ = r_j * (t_x - t_j) / (t_x_ - t_j)
         *r_k_ = r_k * (t_x - t_k) / (t_x_ - t_k)
         *r_node_ = r_node * (upper - t_x) / (upper - t_x_)
         */
       double nu =(upper - t_x) * (t_x - t_j) * (t_x - t_k) ;
       double de = (upper - t_x_) * (t_x_ - t_j) * (t_x_ - t_k);
       hastingsRatio = hastingsRatio + Math.log(nu / de);
       return hastingsRatio;
}



    /*
    Tuning the parameter: twindowsize represents the range of Uniform distribution
     */
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

