package consoperators;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math.MathException;
import beast.base.inference.Distribution;

import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.Operator;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.core.Log;
import beast.base.evolution.branchratemodel.UCRelaxedClockModel;
import beast.base.inference.distribution.ParametricDistribution;
import beast.base.util.Randomizer;
import consoperators.distributions.IndexedPieceWiseLinearDistribution;
import consoperators.distributions.OneParameterMeanOneDistribution;

public class IndexChangeConstantRate extends Operator {

	
	final public Input<IntegerParameter> indexInput = new Input<>("index", "the clock model indicator.", Validate.REQUIRED);
	final public Input<RealParameter> quantileInput = new Input<>("quantiles", "the quantiles of each branch rate.", Validate.REQUIRED);
    final public Input<IndexedPieceWiseLinearDistribution> distrInput = new Input<>("distr", "the clock model approximation function.", Validate.REQUIRED);
    
    
	final boolean DEBUG = false;
    UCRelaxedClockModel clockModel;
    IndexedPieceWiseLinearDistribution distr; 
    RealParameter quantiles;
    IntegerParameter index;
    
    List<OneParameterMeanOneDistribution> clockModelDistributions = new ArrayList<>();
    int numEligibleDistributions = 0;
    double[] rates_original;
	
	@Override
	public void initAndValidate() {
		distr = distrInput.get();
		
		String thisName = this.getID() == "" ? "" + this.getClass() : this.getID();
		
		
		// Find all the eligible clock models. Strict clock is not eligible.
		List<ParametricDistribution> distributions = distr.distrsInput.get();
		for (ParametricDistribution distribution : distributions) {
			try {
				OneParameterMeanOneDistribution d = (OneParameterMeanOneDistribution) distribution;
				
		    	if (d.modeInput.get() == OneParameterMeanOneDistribution.Mode.strict) {
		    		Log.warning(thisName + " will not propose new rates when transiting to/from the strict clock model");
		    	}
		    	
		    	clockModelDistributions.add(d);
		    	numEligibleDistributions ++;
				
				
			}catch(Exception e) {
				e.printStackTrace();
				throw new IllegalArgumentException(thisName + " can not operate on " + distribution.getClass() + " because it is not a OneParameterMeanOneDistribution clock distribution");
	    	}
	    		
	    	
	    	
	    	
		}
		
		
		if (numEligibleDistributions < 2) {
			Log.warning(thisName + " will not recalculate rates because there are less than 2 eligible clock models");
		}
		
		index = indexInput.get();
		quantiles = quantileInput.get();
		
		rates_original = new double[quantiles.getDimension()];
		
		
	}

	@Override
	public double proposal() {
		
		if (numEligibleDistributions < 2) return Double.NEGATIVE_INFINITY;
		
		
		//if (DEBUG) System.out.println("Applying operator " + this.getClass());
		
		
		// Hastings ratio
		double logHR = 0;
		
		
		// The original clock model
		OneParameterMeanOneDistribution originalDist = (OneParameterMeanOneDistribution) distr.getUnderlyingDistr();
		int originalIndex = index.getValue();
		
		
		// Get the rates under the original clock model
		if (originalDist.modeInput.get() != OneParameterMeanOneDistribution.Mode.strict) {
			for (int i = 0; i < quantiles.getDimension(); i ++) {
				
				double q = quantiles.getArrayValue(i);
				
				try {
					rates_original[i] = distr.inverseCumulativeProbability(q);
					logHR += Math.log(distr.getDerivativeAtQuantile(q));
					
					if (distr.getDerivativeAtQuantile(q) < 0) {
						double x = distr.getDerivativeAtQuantile(q);
						if (DEBUG) System.out.println( originalDist.modeInput.get() + " rate = " + rates_original[i] + " distr.getDerivativeAtQuantile(q) = " + distr.getDerivativeAtQuantile(q) + " q =  " + q + " logHR = " + logHR);
						
					}
				
				} catch (MathException e) {
					Log.warning(e.getMessage());
					return Double.NEGATIVE_INFINITY;
				}
			}
		}
		
		
		// New random index. Sample the number of models to move away from the current one (and wrap around using modulus function)
		// Will not re-sample the current model
		int proposedModelIndex = (index.getValue() + 1 + Randomizer.nextInt(clockModelDistributions.size() - 1)) % clockModelDistributions.size();
		
		
		// Get the new clock model distribution
		index.setValue(proposedModelIndex);
		
		OneParameterMeanOneDistribution proposedDist = (OneParameterMeanOneDistribution) distr.getUnderlyingDistr();
		distr.requiresRecalculation();
		
		
		// If either the original or the proposed model is a strict clock, then finish the proposal now
		if (originalDist.modeInput.get() == OneParameterMeanOneDistribution.Mode.strict) return 0;
		if (proposedDist.modeInput.get() == OneParameterMeanOneDistribution.Mode.strict) return 0;
		

		
		// Propose new quantiles such that the rates remain constant
		try {

			double rmin = distr.getRangeMin();
			double rmax = distr.getRangeMax();
			//if (DEBUG) System.out.println("rmin = " + rmin + ", rmax = " + rmax);
			
			// First check that all of the rates are within range of the new distribution
			for (int i = 0; i < quantiles.getDimension(); i ++) {
				
				double rate = rates_original[i];
				if (rate <= rmin || rate >= rmax) {
					//if (DEBUG) System.out.println("Rates out of range: " + rmin + " !< " + rate + " !< " + rmax);
					index.setValue(originalIndex);
					distr.requiresRecalculation();
					return Double.NEGATIVE_INFINITY;
				}
			}
			
			if (DEBUG) System.out.println("Transiting from " + originalDist.modeInput.get() + " to " + proposedDist.modeInput.get());
			
			
			for (int i = 0; i < quantiles.getDimension(); i ++) {
				
				// The rate will stay constant.
				double rate = rates_original[i];
				
				// Propose a new quantile such that the rate remains constant
				double q = distr.cumulativeProbability(rate);
				
				if (q <= 0 || q >= 1) {
					index.setValue(originalIndex);
					distr.requiresRecalculation();
					return Double.NEGATIVE_INFINITY;
				}
				
				logHR += Math.log(distr.getDerivativeAtQuantileInverse(rate, q));
				//if (DEBUG) System.out.println("q_old = " + quantiles.getArrayValue(i) + ", rate = " + rate + ", q_new = " + q);
				quantiles.setValue(i, q);
				
				
				//double temp = distr.getDerivativeAtQuantile(q);
				//double temp2 = 1 / distr.getDerivativeAtQuantileInverse(rate, q);
				//if (DEBUG) System.out.println("temp = " + temp + ", temp2 = " + temp2 + " logHR = " + logHR);
				
				
			}
			
			
			distr.requiresRecalculation();
			
		} catch (MathException e) {
			Log.warning(e.getMessage());
			index.setValue(originalIndex);
			distr.requiresRecalculation();
			return Double.NEGATIVE_INFINITY;
		}

		if (DEBUG) System.out.println("Hastings ratio: " + logHR);
		return logHR;
		

	}

}










