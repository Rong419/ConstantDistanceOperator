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
import org.apache.commons.math3.special.Erf;
import org.apache.commons.math3.util.FastMath;
import java.text.DecimalFormat;
import java.util.List;

@Description("For internal nodes: propose a new node time")
public class InternalOperator extends TreeOperator {
    public final Input<Double> twindowSizeInput =
            new Input<>("twindowSize", "the size of the window when proposing new node time", Input.Validate.REQUIRED);
    final public Input<RealParameter> rateInput = new Input<>("rates", "the rates associated with nodes in the tree for sampling of individual rates among branches.",Input.Validate.REQUIRED);
    final public Input<RealParameter> quantileInput = new Input<>("rateQuantiles", "the rate quantiles associated with nodes in the tree for sampling of individual rates among branches.",Input.Validate.XOR, rateInput);
    final public Input<ParametricDistribution> rateDistInput = new Input<>("distr", "the distribution governing the rates among branches.", Input.Validate.REQUIRED);
    final public Input<RealParameter> stdevInput = new Input<>("stdev", "the distribution governing the rates among branches.");
    public final Input<Double> scaleFactorInput = new Input<>("scaleFactor", "scaling factor: larger means more bold proposals", 1.0);



    enum Mode {
        quantiles,
        rates
    }
    private Mode mode = Mode.rates;
    private double twindowSize;
    private RealParameter rates;
    private RealParameter quantiles;
    private RealParameter ucldStdev;
    double new_stdev;
    double stdev;
    private double m_fScaleFactor;

    protected ParametricDistribution distribution;

    @Override
    public void initAndValidate() {
        twindowSize = twindowSizeInput.get();
        distribution = rateDistInput.get();
        if (quantileInput.get() != null) {
            mode = Mode.quantiles;
            quantiles = quantileInput.get();
        } else {
            rates = rateInput.get();
        }
        stdev = stdevInput.get().getValue();
        m_fScaleFactor = scaleFactorInput.get();
    }

    @Override
    public double proposal() {
        final Tree tree = treeInput.get(this);
        ucldStdev=stdevInput.get();
        //double m_fScaleFactor = 0.5;
        double b =  (m_fScaleFactor + (Randomizer.nextDouble() * ((1.0 / m_fScaleFactor) - m_fScaleFactor)));
        new_stdev = stdev * b;
        ucldStdev.setValue(new_stdev);

        //the original rates
        double r_node;
        double r_k;
        double r_j;

        //Step 1: randomly select an internal node, denoted by node x
        //the chosen node to work on
        Node node;
        int nodeCount = tree.getNodeCount();//return the number of nodes in the tree
        do {
            final int nodeNr = nodeCount / 2 + 1 + Randomizer.nextInt(nodeCount / 2);
            node = tree.getNode(nodeNr);
        } while (node.isRoot() || node.isLeaf());


       //Step 2: Access to the child nodes of this node
       //get the left child of this node, i.e. son
       Node son = node.getChild(0);
       //get the right child of this node, i.e. daughter
       Node daughter = node.getChild(1);


       //node times for these nodes
       double t_x = node.getHeight();
       double t_j = son.getHeight();
       double t_k = daughter.getHeight();

       //original rates for these nodes
       if (mode == Mode.quantiles) {
           // use quantiles
           double q_node = quantiles.getValues()[node.getNr()];
           double q_j = quantiles.getValues()[son.getNr()];
           double q_k = quantiles.getValues()[daughter.getNr()];
           //System.out.println("q1="+q_node+",q2="+q_j+",q3="+q_k);
           r_node = getRealRate(q_node);
           r_j = getRealRate(q_j);
           r_k = getRealRate(q_k);
           //System.out.println("r1="+r_node+",r2="+r_j+",r3="+r_k);
       } else {
           // use real rates
           r_node = rates.getValues()[node.getNr()];
           r_j = rates.getValues()[son.getNr()];
           r_k = rates.getValues()[daughter.getNr()];
       }

       if (r_node == 0.0 || r_j == 0.0 || r_k == 0.0) {
           return Double.NEGATIVE_INFINITY;
       }
       //System.out.println("r1="+r_node+",r2="+r_j+",r3="+r_k);

       //Step3: to propose a new node time for this node
       double a = Randomizer.uniform(-twindowSize, twindowSize);
       double t_x_ = t_x + a;

       //reject the proposal if exceeds the boundary
       double upper = node.getParent().getHeight();
       double lower = Math.max(t_j, t_k);

       if (t_x_<= lower || t_x_ >= upper) {
            return Double.NEGATIVE_INFINITY;
        }
        node.setHeight(t_x + a);


       //Step4: propose the new rates
       double r_node_ = r_node * (upper - t_x) / (upper - t_x_);
       double r_j_ = r_j * (t_x - t_j) / (t_x_ - t_j);
       double r_k_ = r_k * (t_x - t_k) / (t_x_ - t_k);

        //System.out.println("r1'="+r_node_+",r2'="+r_j_+",r3'="+r_k_);
       // set the proposed new rates
        if (mode == Mode.rates) {
            rates.setValue(node.getNr(), r_node_);
            rates.setValue(son.getNr(), r_j_);
            rates.setValue(daughter.getNr(), r_k_);
        }
        else {
            //System.out.println("q1="+getRateQuantiles(r_node_)+",q2="+getRateQuantiles(r_j_)+",q3="+getRateQuantiles(r_k_));
            quantiles.setValue(node.getNr(),getRateQuantiles(r_node_));
            quantiles.setValue(son.getNr(),getRateQuantiles(r_j_));
            quantiles.setValue(daughter.getNr(),getRateQuantiles(r_k_));
        }

       //Step4: calculate the Hastings ratio
       double nu =(upper - t_x) * (t_x - t_j) * (t_x - t_k) ;
       double de = (upper - t_x_) * (t_x_ - t_j) * (t_x_ - t_k);
       double JD =  b * getDicdf(r_node) * getDicdf(r_j) * getDicdf(r_k) * nu / de;

       return Math.log(JD);
}

    private double getRateQuantiles (double r) {
        try {
            return distribution.cumulativeProbability(r);
        } catch (MathException e) {
            throw new RuntimeException("Failed to compute inverse cumulative probability because rate =" + r);
        }
    }

    private double getRealRate (double q) {
        try {
            return distribution.inverseCumulativeProbability(q);
        } catch (MathException e) {
            throw new RuntimeException("Failed to compute inverse cumulative probability because quantile = " + q);
        }
    }

    private double getDicdf(double r) {

        double miu = Math.log(1/(Math.sqrt(1+stdev)));
        double miu_ = Math.log(1/(Math.sqrt(1+new_stdev)));
        double jita = Math.sqrt(Math.log(1 + stdev * stdev));
        double a = 1/(jita * r);
        double b = 1 - Math.pow((Math.log(r) - miu) / (Math.sqrt(2) * jita),2);
        double c = erfInv(2 * getRateQuantiles(r) -1);
        double d = miu_ + Math.sqrt(2)*Math.E*c + c*c + b;

        return a * Math.exp(d);
    }

    /*
    Tuning the parameter: twindowsize represents the range of Uniform distribution
     */
    @Override
    /*public double getCoercableParameterValue() { return twindowSize; }*/
    public double getCoercableParameterValue() {
        return m_fScaleFactor;
    }

    @Override
    /*public void setCoercableParameterValue(double value) {
        twindowSize = value;
    }
    */
    public void setCoercableParameterValue(final double value) {
        double upper = 1.0 - 1e-8;
        double lower = 1e-8;
        m_fScaleFactor = Math.max(Math.min(value, upper), lower);
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

        /*double delta = calcDelta(logAlpha);
        delta += Math.log(twindowSize);
        twindowSize = Math.exp(delta);
        */

        double delta = calcDelta(logAlpha);
        delta += Math.log(1.0 / m_fScaleFactor - 1.0);
        setCoercableParameterValue(1.0 / (Math.exp(delta) + 1.0));
    }

    @Override
    public final String getPerformanceSuggestion() {
        double prob = m_nNrAccepted / (m_nNrAccepted + m_nNrRejected + 0.0);
        double targetProb = getTargetAcceptanceProbability();

        double ratio = prob / targetProb;
        if (ratio > 2.0) ratio = 2.0;
        if (ratio < 0.5) ratio = 0.5;

        //double newWindowSize = twindowSize * ratio;
        final double sf = Math.pow(m_fScaleFactor, ratio);

        DecimalFormat formatter = new DecimalFormat("#.###");
        if (prob < 0.10) {
            //return "Try setting window size to about " + formatter.format(newWindowSize);
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else if (prob > 0.40) {
            //return "Try setting window size to about " + formatter.format(newWindowSize);
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else return "";
    }

   private static double erfInv(final double x) {
    // beware that the logarithm argument must be
    // commputed as (1.0 - x) * (1.0 + x),
    // it must NOT be simplified as 1.0 - x * x as this
    // would induce rounding errors near the boundaries +/-1
             double w = - FastMath.log((1.0 - x) * (1.0 + x));
             double p;

             if (w < 6.25) {
                        w -= 3.125;
                        p =  -3.6444120640178196996e-21;
                        p =   -1.685059138182016589e-19 + p * w;
                        p =   1.2858480715256400167e-18 + p * w;
                        p =    1.115787767802518096e-17 + p * w;
                        p =   -1.333171662854620906e-16 + p * w;
                        p =   2.0972767875968561637e-17 + p * w;
                        p =   6.6376381343583238325e-15 + p * w;
                        p =  -4.0545662729752068639e-14 + p * w;
                        p =  -8.1519341976054721522e-14 + p * w;
                        p =   2.6335093153082322977e-12 + p * w;
                        p =  -1.2975133253453532498e-11 + p * w;
                        p =  -5.4154120542946279317e-11 + p * w;
                        p =    1.051212273321532285e-09 + p * w;
                        p =  -4.1126339803469836976e-09 + p * w;
                        p =  -2.9070369957882005086e-08 + p * w;
                        p =   4.2347877827932403518e-07 + p * w;
                        p =  -1.3654692000834678645e-06 + p * w;
                        p =  -1.3882523362786468719e-05 + p * w;
                        p =    0.0001867342080340571352 + p * w;
                        p =  -0.00074070253416626697512 + p * w;
                        p =   -0.0060336708714301490533 + p * w;
                        p =      0.24015818242558961693 + p * w;
                        p =       1.6536545626831027356 + p * w;
                    } else if (w < 16.0) {
                        w = FastMath.sqrt(w) - 3.25;
                        p =   2.2137376921775787049e-09;
                        p =   9.0756561938885390979e-08 + p * w;
                        p =  -2.7517406297064545428e-07 + p * w;
                        p =   1.8239629214389227755e-08 + p * w;
                        p =   1.5027403968909827627e-06 + p * w;
                        p =   -4.013867526981545969e-06 + p * w;
                        p =   2.9234449089955446044e-06 + p * w;
                      p =   1.2475304481671778723e-05 + p * w;
                        p =  -4.7318229009055733981e-05 + p * w;
                        p =   6.8284851459573175448e-05 + p * w;
                        p =   2.4031110387097893999e-05 + p * w;
                        p =   -0.0003550375203628474796 + p * w;
                        p =   0.00095328937973738049703 + p * w;
                        p =   -0.0016882755560235047313 + p * w;
                        p =    0.0024914420961078508066 + p * w;
                        p =   -0.0037512085075692412107 + p * w;
                        p =     0.005370914553590063617 + p * w;
                        p =       1.0052589676941592334 + p * w;
                        p =       3.0838856104922207635 + p * w;
                    } else if (!Double.isInfinite(w)) {
                        w = FastMath.sqrt(w) - 5.0;
                        p =  -2.7109920616438573243e-11;
                        p =  -2.5556418169965252055e-10 + p * w;
                        p =   1.5076572693500548083e-09 + p * w;
                        p =  -3.7894654401267369937e-09 + p * w;
                        p =   7.6157012080783393804e-09 + p * w;
                        p =  -1.4960026627149240478e-08 + p * w;
                        p =   2.9147953450901080826e-08 + p * w;
                        p =  -6.7711997758452339498e-08 + p * w;
                        p =   2.2900482228026654717e-07 + p * w;
                        p =  -9.9298272942317002539e-07 + p * w;
                        p =   4.5260625972231537039e-06 + p * w;
                        p =  -1.9681778105531670567e-05 + p * w;
                        p =   7.5995277030017761139e-05 + p * w;
                        p =  -0.00021503011930044477347 + p * w;
                       p =  -0.00013871931833623122026 + p * w;
                    p =       1.0103004648645343977 + p * w;
                       p =       4.8499064014085844221 + p * w;
                  } else {
                   p = Double.POSITIVE_INFINITY;
                 }

               return p * x;

          }


}

