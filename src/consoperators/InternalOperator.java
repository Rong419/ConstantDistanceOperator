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


    enum Mode {
        quantiles,
        rates
    }
    private Mode mode = Mode.rates;
    private double twindowSize;
    private RealParameter rates;
    private RealParameter quantiles;
    private RealParameter ucldStdev;

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
        ucldStdev = stdevInput.get();
    }

    @Override
    public double proposal() {
        final Tree tree = treeInput.get(this);
        double stdev = ucldStdev.getValue();
        double m_fScaleFactor = 0.5;
        double b =  (m_fScaleFactor + (Randomizer.nextDouble() * ((1.0 / m_fScaleFactor) - m_fScaleFactor)));
        double new_stdev = stdev * b;
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
           System.out.println("q1="+q_node+",q2="+q_j+",q3="+q_k);
           r_node = getRealRate(q_node);
           r_j = getRealRate(q_j);
           r_k = getRealRate(q_k);
           System.out.println("r1="+r_node+",r2="+r_j+",r3="+r_k);
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
       double icdfq = getDicdf(getRateQuantiles(r_node_)) * getDicdf(getRateQuantiles(r_j_)) * getDicdf(getRateQuantiles(r_k_));
       double JD =  b * nu /de;

       return Math.log(JD);
}

    private double getRateQuantiles (double r) {
        try {
            LogNormalDistributionModel mylognormal=(LogNormalDistributionModel)distribution;
            double S=mylognormal.SParameterInput.get().getValue();
            //System.out.println("S="+S);
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

    private double getDicdf(double q) {
        LogNormalDistributionModel mylognormal=(LogNormalDistributionModel)distribution;
        double S=mylognormal.SParameterInput.get().getValue();
        double miu = Math.log(1/(Math.sqrt(1+S)));
        return miu + 4 * Math.sqrt(2) * Math.E *((2*q-1)* (Erf.erf(2*q-1)) + (Math.exp(-(2*q-1)*(2*q-1)))/(Math.sqrt(Math.PI)));
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

