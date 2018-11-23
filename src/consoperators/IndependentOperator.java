package consoperators;

import beast.core.Description;
import beast.core.Input;
import beast.core.Operator;
import beast.core.parameter.RealParameter;
import beast.evolution.operators.TreeOperator;
import beast.evolution.tree.Tree;
import beast.math.distributions.LogNormalDistributionModel;
import beast.math.distributions.ParametricDistribution;
import beast.util.Randomizer;
import beast.evolution.tree.Node;

@Description("IndependentOperator: propose a new rate to replace the current rate")
public class IndependentOperator extends Operator {
    final public Input<RealParameter> rateInput = new Input<>("rates",
            "the rates associated with nodes", Input.Validate.REQUIRED);
    final public Input<ParametricDistribution> rateDistInput = new Input<>("distr", "the distribution g" +
            "overning the rates among branches." , Input.Validate.REQUIRED);

    ParametricDistribution distribution; //the distribution of the rates
    RealParameter rates;
    //private Tree tree;
    double M;double S;boolean MeanInRealSpace;

    @Override
    public void initAndValidate() {
        distribution = rateDistInput.get();
        rates = rateInput.get();
        //tree = treeInput.get();
    }

    @Override
    public double proposal() {
        LogNormalDistributionModel L=(LogNormalDistributionModel)distribution;
        M = L.MParameterInput.get().getValue();
        S = L.SParameterInput.get().getValue();
        MeanInRealSpace = L.hasMeanInRealSpaceInput.get();

        int index =Randomizer.nextInt(rates.getDimension());

        double r = Randomizer.nextLogNormal(M,S,MeanInRealSpace);
        if (r <= 0.0){
            return Double.NEGATIVE_INFINITY;
        }

        rates.setValue(index, r);

    return  0.0;
    }

}
