package consoperators;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.Operator;
import beast.base.inference.parameter.RealParameter;
import beast.base.core.Log;
import beast.base.inference.distribution.Exponential;
import beast.base.inference.distribution.LogNormalDistributionModel;
import beast.base.inference.distribution.ParametricDistribution;
import beast.base.util.Randomizer;
import org.apache.commons.math.MathException;


@Description("IndependentOperator: randomly propose a new rate to replace the current rate")
public class IndependentOperator extends Operator {
    final public Input<RealParameter> rateInput = new Input<>("rates",
            "the rates associated with nodes", Input.Validate.REQUIRED);
    final public Input<ParametricDistribution> rateDistInput = new Input<>("distr", "the distribution g" +
            "overning the rates among branches." , Input.Validate.REQUIRED);

    ParametricDistribution distribution; //the distribution of the rates
    RealParameter rates;
    double M;double S;boolean MeanInRealSpace;
    double lamda;
    double r;int index;

    @Override
    public void initAndValidate() {
        distribution = rateDistInput.get();
        rates = rateInput.get();
    }

    @Override
    public double proposal() {

       double P_old;double P_new;

        //randomly select a rate to work on
       index =Randomizer.nextInt(rates.getDimension());

       double old_r = rates.getValue(index);

        //for rates follow a LogNormal distribution
        if(distribution instanceof LogNormalDistributionModel) {
            LogNormalDistributionModel L = (LogNormalDistributionModel) distribution;
            M = L.MParameterInput.get().getArrayValue();
            S = L.SParameterInput.get().getArrayValue();
            MeanInRealSpace = L.hasMeanInRealSpaceInput.get();
            r  = Randomizer.nextLogNormal(M,S,MeanInRealSpace);
        }

        //for rates follow an Exponential distribution
        else if(distribution instanceof LogNormalDistributionModel) {
            Exponential E = (Exponential) distribution;
            lamda = E.lambdaInput.get().getArrayValue();
            r = Randomizer.nextExponential(lamda);
        }
        else {
            Log.info.println("IndependentOperator: the required distribution is unavailable.");
        }


        if (r <= 0.0){
            return Double.NEGATIVE_INFINITY;
        }

        //set the proposed value for the selected rate
        rates.setValue(index, r);


        //A Gibbs sampler: the proposal is always accepted
        // ????
        //return Double.MAX_VALUE;

        try {
            P_new = distribution.cumulativeProbability(r);
            P_old = distribution.cumulativeProbability(old_r);
        } catch (MathException e) {
            throw new RuntimeException("Failed to compute cumulative probability!");
        }

        return Math.log(P_new / P_old);
    }

}
