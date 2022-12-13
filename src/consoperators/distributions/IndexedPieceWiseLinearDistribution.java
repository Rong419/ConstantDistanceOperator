package consoperators.distributions;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.math.distribution.ContinuousDistribution;
import beast.base.inference.Distribution;
import beast.base.inference.distribution.LogNormalDistributionModel;
import beast.base.inference.distribution.ParametricDistribution;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;

@Description("Selects a parametricdistribution from among its inputs based on index variable and "
		+ "approximates distribution by piecewise linear approximation.")
public class IndexedPieceWiseLinearDistribution extends PiecewiseLinearDistribution {
	final public Input<IntegerParameter> indexInput = new Input<>("index", "index of the distribution to be selected", Validate.REQUIRED);
	final public Input<List<ParametricDistribution>> distrsInput = new Input<>("paramdistr", "one or more parametric distributions to be chosen from", new ArrayList<>(), Validate.REQUIRED);
		
	IntegerParameter index;
	List<ParametricDistribution> distrs;
	

    public IndexedPieceWiseLinearDistribution() {
    	distrInput.setRule(Validate.FORBIDDEN);
    }
    
    /** construct log-normal distribution with mean = 1 **/ 
    public IndexedPieceWiseLinearDistribution(RealParameter S) {
    	LogNormalDistributionModel lnormal = new LogNormalDistributionModel();
    	lnormal.initByName("S", S, "M", "1.0", "meanInRealSpace", true);
    	initByName("distr", lnormal);
    }
    
    @Override
	public void initAndValidate() {
    	int dim = numberOfDiscreteRates.get();
    	rates = new double[dim];
    	storedRates = new double[dim];

    	for (ParametricDistribution distribution : distrsInput.get()) {
    		org.apache.commons.math.distribution.Distribution d = distribution.getDistribution();
	    	if (d instanceof ContinuousDistribution) {
	    		underlyingDistr = (ContinuousDistribution) d;
	    	} else {
	    		throw new IllegalArgumentException("Expected parametric distribution that is continuous");
	    	}
    	}
        refresh();

        distrs = distrsInput.get();
		
		index = indexInput.get();
		index.setLower(0);
		index.setUpper(distrs.size()-1);
		underlyingDistr = distrsInput.get().get(index.getValue());
		
		cutOffEnd = cutOffEndInput.get();

    	limit_0 = limitInput.get();
    	if (limit_0 <= 0 || limit_0 >= 1) {
    		throw new IllegalArgumentException("limit should be between 0 and 1");
    	}
		limitLow = limit_0 / dim;
    	limitUp = 1.0 - limit_0 / dim;

		if (offsetInput.get() != 0.0) {
			throw new IllegalArgumentException("Offset should be set on inputs of this IndexedDistribution, not on the IndexedDistribution itself");
		}
    }

    @Override
    public ContinuousDistribution getUnderlyingDistr() {
    	return distrs.get(index.getValue());
	}
    
    @Override
    public boolean requiresRecalculation() {
    	
		underlyingDistr = distrs.get(index.getValue());
    	if (index.somethingIsDirty() || ((ParametricDistribution) underlyingDistr).isDirtyCalculation()) {
    		Arrays.fill(rates, 0.0);
    		return true;
    	}
    	
    	return super.requiresRecalculation();
    }
	
}
