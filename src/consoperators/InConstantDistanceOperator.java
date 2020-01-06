package consoperators;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.branchratemodel.UCRelaxedClockModel;
import beast.evolution.operators.TreeOperator;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.math.distributions.LogNormalDistributionModel;
import beast.math.distributions.ParametricDistribution;
import beast.math.distributions.PiecewiseLinearDistribution;
import beast.util.Randomizer;
import org.apache.commons.math.MathException;
import java.text.DecimalFormat;


@Description("For internal nodes: propose a new node time")
public class InConstantDistanceOperator extends TreeOperator {
    public final Input<Double> twindowSizeInput =
            new Input<>("twindowSize", "the size of the window when proposing new node time", Input.Validate.REQUIRED);
    final public Input<RealParameter> rateInput = new Input<>("rates", "the rates associated with nodes in the tree for sampling of individual rates among branches.", Input.Validate.REQUIRED);
    final public Input<RealParameter> quantileInput = new Input<>("quantiles", "the quantiles of each branch rate.", Input.Validate.XOR,rateInput);
    final public Input<UCRelaxedClockModel> clockModelInput = new Input<>("clockModel", "relaxed clock model used to deal with quantiles", Input.Validate.REQUIRED);

    private double twindowSize;
    private RealParameter rates;
    private RealParameter quantiles;
    private enum rateMode {
        quantiles,
        rates
    }
    private rateMode mode = rateMode.rates;

    @Override
    public void initAndValidate() {
        twindowSize = twindowSizeInput.get();
        if (rateInput.get() == null) {
            quantiles = quantileInput.get();
            mode = rateMode.quantiles;
        } else {
            rates = rateInput.get();
            mode = rateMode.rates;
        }
    }

    @Override
    public double proposal() {
        final Tree tree = treeInput.get(this);
        ParametricDistribution rateDistribution = clockModelInput.get().rateDistInput.get();
        int nodeCount = tree.getNodeCount(); //return the number of nodes in the tree
        int branchCount = nodeCount - 1; //the number of branches of the tree

        // the chosen node to work on
        Node node;

        // the original node times
        double t_x, t_j, t_k;

        // the original rates
        double r_x, r_k, r_j;

        // the original quantiles
        double q_x = 0.5; double q_k = 0.5; double q_j = 0.5;

        // the proposed quantiles
        double q_x_ = 0.5; double q_k_ = 0.5; double q_j_ = 0.5;

        double hastingsRatio = 0.0;

        // Step 1: randomly select an internal node, denoted by node x.
        // avoid fake nodes used to connect direct ancestors into tree.
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

       // time for this node
       t_x = node.getHeight();


       // Step 2: Access to the child nodes of this node
       // son
       Node son = node.getChild(0); // get the left child of this node, i.e. son
       t_j = son.getHeight(); // node time of son
       int sonNr = son.getNr(); // node number of son
       if (sonNr == branchCount) {
           sonNr = son.getTree().getRoot().getNr();
        }

       // daughter
       Node daughter = node.getChild(1); //get the right child of this node, i.e. daughter
       t_k = daughter.getHeight(); // node time of daughter
       int dauNr = daughter.getNr(); // node number of daughter
       if (dauNr == branchCount) {
            dauNr = daughter.getTree().getRoot().getNr();
       }

       switch (mode) {
           case rates: {
               r_x = rates.getValues()[nodeNr]; // rate of branch above this node
               r_j = rates.getValues()[sonNr]; // rate of branch above son
               r_k = rates.getValues()[dauNr]; // rate of branch above daughter
               break;
           }

           case quantiles: {
               q_x = quantiles.getValues()[nodeNr];
               q_j = quantiles.getValues()[sonNr];
               q_k = quantiles.getValues()[dauNr];
               try {
                   r_x = rateDistribution.inverseCumulativeProbability(q_x);
                   r_j = rateDistribution.inverseCumulativeProbability(q_j);
                   r_k = rateDistribution.inverseCumulativeProbability(q_k);
               } catch (MathException e) {
                   e.printStackTrace();
                   return Double.NEGATIVE_INFINITY;
               }
               break;
           }

           default: {
               return Double.NEGATIVE_INFINITY;
           }
       }


       // Step3: to propose a new node time for this node
       double a = Randomizer.uniform(-twindowSize, twindowSize);
       double t_x_ = t_x + a;

       // deal with the boundary cases
       double upper = node.getParent().getHeight();
       double lower = Math.max(t_j, t_k);

        if (t_x_ == lower || t_x_ == upper) {
            return Double.NEGATIVE_INFINITY;
        }
        // fold the proposed node time
        double err; double n; double r;
        if (t_x_ > upper) {
            err = t_x_ - upper;
            n = Math.floor(err / (upper - lower));
            r = err - n * (upper - lower);
            if (n % 2 == 0) {
                t_x_ = upper - r;
            } else {
                t_x_ = lower + r;
            }
        } else if (t_x_ < lower) {
            err = lower - t_x_;
            n = Math.floor(err / (upper - lower));
            r = err - n * (upper - lower);
            if (n % 2 == 0) {
                t_x_ = lower + r;
            } else {
                t_x_ = upper - r;
            }
        }

        // set the proposed node time
        node.setHeight(t_x_);


       // Step4: propose the new rates
       // there are three rates in total
       // r_x_, r_j_, r_k_
       double r_x_ = r_x * (upper - t_x) / (upper - t_x_);
       double r_j_ = r_j * (t_x - t_j) / (t_x_ - t_j);
       double r_k_ = r_k * (t_x - t_k) / (t_x_ - t_k);


       // set the proposed new rates or quantiles
        switch (mode) {
            case rates: {
                // set rates directly
                rates.setValue(nodeNr, r_x_);
                rates.setValue(sonNr, r_j_);
                rates.setValue(dauNr, r_k_);
                break;
            }

            case quantiles: {
                try {
                    // new quantiles of proposed rates
                    q_x_ = rateDistribution.cumulativeProbability(r_x_);
                    q_j_ = rateDistribution.cumulativeProbability(r_j_);
                    q_k_ = rateDistribution.cumulativeProbability(r_k_);

                    // set quantiles
                    quantiles.setValue(nodeNr, q_x_);
                    quantiles.setValue(sonNr, q_j_);
                    quantiles.setValue(dauNr, q_k_);

                } catch (MathException e) {
                    e.printStackTrace();
                    return Double.NEGATIVE_INFINITY;
                }
                break;
            }

            default: {

            }

        }


        // Step4: calculate the Hastings ratio
        switch (mode) {
            case rates: {
                double nu =(upper - t_x) * (t_x - t_j) * (t_x - t_k) ;
                double de = (upper - t_x_) * (t_x_ - t_j) * (t_x_ - t_k);
                hastingsRatio = Math.log(nu / de);
                break;
            }

            case quantiles: {
                if (rateDistribution instanceof LogNormalDistributionModel) {
                    hastingsRatio = ConsOperatorUtils.getHRForLN(r_x_, q_x, rateDistribution)
                                  + ConsOperatorUtils.getHRForLN(r_j_, q_j, rateDistribution)
                                  + ConsOperatorUtils.getHRForLN(r_k_, q_k, rateDistribution);
                }

                else if (rateDistribution instanceof PiecewiseLinearDistribution) {
                    hastingsRatio = ConsOperatorUtils.getHRForPieceWise(r_x_, q_x, q_x_, rateDistribution)
                                  + ConsOperatorUtils.getHRForPieceWise(r_j_, q_j, q_j_, rateDistribution)
                                  + ConsOperatorUtils.getHRForPieceWise(r_k_, q_k, q_k_, rateDistribution);
                }

                else {
                    hastingsRatio = ConsOperatorUtils.getHRUseNumericApproximation(r_x_, q_x, rateDistribution)
                                  + ConsOperatorUtils.getHRUseNumericApproximation(r_j_, q_j, rateDistribution)
                                  + ConsOperatorUtils.getHRUseNumericApproximation(r_k_, q_k, rateDistribution);
                }

                break;
            }

            default: {

            }

        }
        return hastingsRatio;
}

    // Tuning the parameter: twindowsize represents the range of Uniform distribution
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

