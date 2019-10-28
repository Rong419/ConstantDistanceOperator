package consoperators;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.branchratemodel.UCRelaxedClockModel;
import beast.evolution.operators.TreeOperator;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.math.distributions.LogNormalDistributionModel;
import beast.math.distributions.ParametricDistribution;
import beast.util.Randomizer;
import org.apache.commons.math.MathException;
import org.apache.commons.math3.util.FastMath;

import java.text.DecimalFormat;
import java.util.List;

@Description("For internal nodes: propose a new node time")
public class InConstantDistanceOperator extends TreeOperator {
    public final Input<Double> twindowSizeInput =
            new Input<>("twindowSize", "the size of the window when proposing new node time", Input.Validate.REQUIRED);
    final public Input<RealParameter> rateInput = new Input<>("rates", "the rates associated with nodes in the tree for sampling of individual rates among branches.", Input.Validate.REQUIRED);
    final public Input<RealParameter> quantileInput = new Input<>("quantiles", "the quantiles of each branch rate.", Input.Validate.XOR,rateInput);
    final public Input<ParametricDistribution> rateDistInput = new Input<>("distr", "the distribution governing the rates among branches.");

    private double twindowSize;
    private RealParameter rates;
    private RealParameter quantiles;
    private LogNormalDistributionModel LN ;
    protected ParametricDistribution distribution;

    @Override
    public void initAndValidate() {
        twindowSize = twindowSizeInput.get();

        distribution = rateDistInput.get();

        if (quantileInput.get() == null) {
            rates = rateInput.get();
        } else {
            quantiles = quantileInput.get();
        }
    }

    @Override
    public double proposal() {
        final Tree tree = treeInput.get(this);
        int nodeCount = tree.getNodeCount(); //return the number of nodes in the tree
        int branchCount = nodeCount - 1; //the number of branches of the tree

        if(distribution instanceof LogNormalDistributionModel) {
            LN = (LogNormalDistributionModel) distribution;
        }


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

        //the proposed node time
        double t_x_;

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
       //double r_node = rates.getValues()[nodeNr];
        double r_node = getRateForNode(nodeNr);

       //Step 2: Access to the child nodes of this node
       // son
       Node son = node.getChild(0);//get the left child of this node, i.e. son
       t_j = son.getHeight();//node time of son

       int sonNr = son.getNr();// node number of son
       if (sonNr == branchCount) {
           sonNr = son.getTree().getRoot().getNr();
        }

       //r_j = rates.getValues()[sonNr]; // rate of branch above son
       r_j = getRateForNode(sonNr);
       d_j = r_j * (t_x - t_j); // distance of branch above son


       // daughter
       Node daughter = node.getChild(1);//get the right child of this node, i.e. daughter
       t_k = daughter.getHeight();//node time of daughter

       int dauNr = daughter.getNr(); // node time of daughter
       if (dauNr == branchCount) {
            dauNr = daughter.getTree().getRoot().getNr();
       }

       //r_k = rates.getValues()[dauNr];// rate of branch above daughter
       r_k = getRateForNode(dauNr);
       d_k = r_k * (t_x - t_k);// distance of branch above daughter

       double q_node = 1.0; double q_j = 1.0; double q_k = 1.0;
       if (quantileInput.get()!=null) {
          q_node = quantiles.getValues()[nodeNr];
          q_j = quantiles.getValues()[sonNr];
          q_k = quantiles.getValues()[dauNr];
       }

       //Step3: to propose a new node time for this node
       double a = Randomizer.uniform(-twindowSize, twindowSize);
       t_x_ = t_x + a;

       //reject the proposal if exceeds the boundary
       double upper = node.getParent().getHeight();
       double lower = Math.max(t_j, t_k);

       if (t_x_<= lower || t_x_ >= upper) {
            return Double.NEGATIVE_INFINITY;
        }
        node.setHeight(t_x + a);


       //Step4: propose the new rates
       //there are three rates in total
       //r_node, r_i, r_x
       double r_node_ = r_node * (upper - t_x) / (upper - t_x_);
       double r_j_ = d_j / (t_x_ - t_j);
       double r_k_ = d_k / (t_x_ - t_k);

       // set the proposed new rates
       //rates.setValue(nodeNr, r_node_);
       //rates.setValue(sonNr, r_j_);
       //rates.setValue(dauNr, r_k_);
        assignRateForNode(r_node_,nodeNr);
        assignRateForNode(r_j_,sonNr);
        assignRateForNode(r_k_,dauNr);


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
       double JD = Math.log(nu /de);
       if (quantileInput.get() != null) {
           JD = JD + calculateHastingsRatio(r_node_,q_node) + calculateHastingsRatio(r_j_,q_j) + calculateHastingsRatio(r_k_, q_k);
       }
       return JD;
}

    private double getRateForNode(int nodeNr) {
        if (quantileInput.get() == null) {
            return rates.getValues()[nodeNr];
        }
        else {
            return getRealRate(quantiles.getValues()[nodeNr]);
        }
    }

    private void assignRateForNode(double rate, int nodeNr) {
        if (quantileInput.get() == null) {
            rates.setValue(nodeNr,rate);
        } else {
            quantiles.setValue(nodeNr,getRateQuantiles(rate));
        }
    }

    private double getRateQuantiles (double r) {
        try {
            return LN.cumulativeProbability(r);
        } catch (MathException e) {
            throw new RuntimeException("Failed to compute cumulative probability because rate =" + r);
        }
    }

    private double getRealRate (double q) {
        try {
            return LN.inverseCumulativeProbability(q);
        } catch (MathException e) {
            throw new RuntimeException("Failed to compute inverse cumulative probability because quantile = " + q);
        }
    }

    private double calculateHastingsRatio(double r, double q) {
        double stdev = LN.SParameterInput.get().getValue();
        double variance = FastMath.exp(stdev * stdev) -1; // sigma square of lognormal
        double miu = - 0.5 * FastMath.log(1 + variance); // miu of lognormal

        double a = erfInv(2 * q - 1);
        double b = FastMath.log(r);
        double c = 2 * stdev * stdev;
        double x = b - miu;
        double x_sq = x * x / c;
        double d = Math.sqrt(c);
        return -b - (x_sq/c) + miu + (d * a) + (a * a);
    }

    
    public static double calculateHastingsRatio(double r, double q, double stdev) {
        double variance = FastMath.exp(stdev * stdev) -1; // sigma square of lognormal
        double miu = - 0.5 * FastMath.log(1 + variance); // miu of lognormal

        double a = erfInv(2 * q - 1);
        double b = FastMath.log(r);
        double c = 2 * stdev * stdev;
        double x = b - miu;
        double x_sq = x * x / c;
        double d = Math.sqrt(c);
        return -b - (x_sq/c) + miu + (d * a) + (a * a);
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

    private static double erfInv(final double x) {
        double w = - FastMath.log((1.0 - x) * (1.0 + x));
        double p;

        if (w < 6.25) {
            w -= 3.125;
            p = -3.6444120640178196996e-21;
            p = -1.685059138182016589e-19 + p * w;
            p = 1.2858480715256400167e-18 + p * w;
            p = 1.115787767802518096e-17 + p * w;
            p = -1.333171662854620906e-16 + p * w;
            p = 2.0972767875968561637e-17 + p * w;
            p = 6.6376381343583238325e-15 + p * w;
            p = -4.0545662729752068639e-14 + p * w;
            p = -8.1519341976054721522e-14 + p * w;
            p = 2.6335093153082322977e-12 + p * w;
            p = -1.2975133253453532498e-11 + p * w;
            p = -5.4154120542946279317e-11 + p * w;
            p = 1.051212273321532285e-09 + p * w;
            p = -4.1126339803469836976e-09 + p * w;
            p = -2.9070369957882005086e-08 + p * w;
            p = 4.2347877827932403518e-07 + p * w;
            p = -1.3654692000834678645e-06 + p * w;
            p = -1.3882523362786468719e-05 + p * w;
            p = 0.0001867342080340571352 + p * w;
            p = -0.00074070253416626697512 + p * w;
            p = -0.0060336708714301490533 + p * w;
            p = 0.24015818242558961693 + p * w;
            p = 1.6536545626831027356 + p * w;
        } else if (w < 16.0) {
            w = FastMath.sqrt(w) - 3.25;
            p = 2.2137376921775787049e-09;
            p = 9.0756561938885390979e-08 + p * w;
            p = -2.7517406297064545428e-07 + p * w;
            p = 1.8239629214389227755e-08 + p * w;
            p = 1.5027403968909827627e-06 + p * w;
            p = -4.013867526981545969e-06 + p * w;
            p = 2.9234449089955446044e-06 + p * w;
            p = 1.2475304481671778723e-05 + p * w;
            p = -4.7318229009055733981e-05 + p * w;
            p = 6.8284851459573175448e-05 + p * w;
            p = 2.4031110387097893999e-05 + p * w;
            p = -0.0003550375203628474796 + p * w;
            p = 0.00095328937973738049703 + p * w;
            p = -0.0016882755560235047313 + p * w;
            p = 0.0024914420961078508066 + p * w;
            p = -0.0037512085075692412107 + p * w;
            p = 0.005370914553590063617 + p * w;
            p = 1.0052589676941592334 + p * w;
            p = 3.0838856104922207635 + p * w;
        } else if (!Double.isInfinite(w)) {
            w = FastMath.sqrt(w) - 5.0;
            p = -2.7109920616438573243e-11;
            p = -2.5556418169965252055e-10 + p * w;
            p = 1.5076572693500548083e-09 + p * w;
            p = -3.7894654401267369937e-09 + p * w;
            p = 7.6157012080783393804e-09 + p * w;
            p = -1.4960026627149240478e-08 + p * w;
            p = 2.9147953450901080826e-08 + p * w;
            p = -6.7711997758452339498e-08 + p * w;
            p = 2.2900482228026654717e-07 + p * w;
            p = -9.9298272942317002539e-07 + p * w;
            p = 4.5260625972231537039e-06 + p * w;
            p = -1.9681778105531670567e-05 + p * w;
            p = 7.5995277030017761139e-05 + p * w;
            p = -0.00021503011930044477347 + p * w;
            p = -0.00013871931833623122026 + p * w;
            p = 1.0103004648645343977 + p * w;
            p = 4.8499064014085844221 + p * w;
        } else {
            p = Double.POSITIVE_INFINITY;
        }
        return p * x;
    }

}

