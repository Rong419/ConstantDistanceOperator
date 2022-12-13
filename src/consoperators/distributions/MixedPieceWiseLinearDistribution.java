package consoperators.distributions;


import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.ContinuousDistribution;
import beast.base.inference.Distribution;
import beast.base.inference.distribution.LogNormalDistributionModel;
import beast.base.inference.distribution.ParametricDistribution;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.RealParameter;


@Description("Mixture of parametricdistribution weighted among its inputs by a weight vector and "
		+ "approximates distribution by piecewise linear approximation.")
public class MixedPieceWiseLinearDistribution extends ParametricDistribution {
	final public Input<RealParameter> weightInput = new Input<>("weight", "weights of the distribution to be selected, should have same dimension as number of parameters", new RealParameter("1.0"));
	final public Input<List<ParametricDistribution>> distrsInput = new Input<>("paramdistr", "one or more parametric distributions to be chosen from", new ArrayList<>(), Validate.REQUIRED);
    final public Input<Integer> numberOfDiscreteRates = new Input<>("bins", "number of bins used to approximate the distribution piecewise linearly. (default 100)", 100);
		
	RealParameter weight;
	List<ParametricDistribution> distrs;
	
	ContinuousDistribution mixtureDistr;
    protected double[] rates; // the output rates of mixture distribution
    protected double[] storedRates; // stored rates

    protected double[][] drates; // rates of component distributions
    protected double[][] storedDRates; // stored drates

    public MixedPieceWiseLinearDistribution() {
    }
    
    /** construct log-normal distribution with mean = 1 **/ 
    public MixedPieceWiseLinearDistribution(RealParameter S) {
    	LogNormalDistributionModel lnormal = new LogNormalDistributionModel();
    	lnormal.initByName("S", S, "M", "1.0", "meanInRealSpace", true);
    	initByName("distr", lnormal);
    }
    
    @Override
	public void initAndValidate() {
        distrs = distrsInput.get();
        int dim = numberOfDiscreteRates.get();
        
    	rates = new double[dim];
    	drates = new double[distrs.size()][dim + 2]; // + 2 to deal with boundary conditions
    	storedRates = new double[dim];
    	storedDRates = new double[distrs.size()][dim + 2];

    	for (ParametricDistribution distribution : distrsInput.get()) {
    		org.apache.commons.math.distribution.Distribution d = distribution.getDistribution();
	    	if (!(d instanceof ContinuousDistribution)) {
	    		throw new IllegalArgumentException("Expected parametric distribution that is continuous");
	    	}
    	}

		
		weight = weightInput.get();
		weight.setLower(0.0);
		weight.setUpper(1.0);
		
		if (offsetInput.get() != 0.0) {
			throw new IllegalArgumentException("Offset should be set on inputs of this MixtureDistribution, not on the MixtureDistribution itself");
		}
		
		if (weight.getDimension() != distrs.size()) {
			throw new IllegalArgumentException("There should be one weight per distribution, but there are " + 
					weight.getDimension() + " weights and " + distrs.size() + " distributions. Change the dimension of the " +
			        "weight parameter to fix this.");
		}
		
		mixtureDistr = new MixtureLinearPiecewiseImpl();
		refresh(true);
    }

    public class MixtureLinearPiecewiseImpl implements ContinuousDistribution {

        public MixtureLinearPiecewiseImpl() {
        }

        @Override
        public double cumulativeProbability(double x) throws MathException {
        	// Return exact cdf using piecewise linear approximation
            int i = Arrays.binarySearch(rates, x);
            if (i < 0) {
            	i = -i-1;
            }
            if (i < rates.length-1) {
            	double cdf = (i + (x - rates[i]) / (rates[i+1] - rates[i])) / (rates.length-1);
            	cdf = Math.max(0, cdf);
            	return cdf;
            }
            return 1.0;
        }
        
		@Override
        public double cumulativeProbability(double x0, double x1) throws MathException {
            return cumulativeProbability(x1) - cumulativeProbability(x0);
        }

        @Override
        public double inverseCumulativeProbability(double q) throws MathException {
            double v = q * (rates.length - 1);
            int i = (int) v;
            
            // make sure cached rates are calculated
            if (rates[i] == 0.0) {
            	refresh(false);
            }
            
            // return piecewise linear approximation
            double r = rates[i];
            if (i < rates.length - 1) {
            	r += (rates[i+1] - rates[i]) * (v - i);
            }
            return r;
        }

        @Override
        public double density(double x) {
        	throw new IllegalArgumentException("not implemented yet");
        }

        @Override
        public double logDensity(double x) {
        	return Math.log(density(x));
        }

    } // class MixtureLinearPiecewiseImpl

    
    void refresh(boolean initial) {
    	int k = 0;
		int n1 = rates.length-1;
		double n = rates.length - 1;
		int m = distrs.size();
    	double [] weights = weight.getDoubleValues();

    	// sanity check
    	double sum = 0;
    	for (double d : weights) {
    		sum += d;
    	}
    	if (Math.abs(sum - 1.0) > 1e-10) {
    		throw new IllegalArgumentException("Weights should sum to 1 (not " + sum +")");
    	}
    	
    	try {
        	// collect rates from underlying distributions
	    	for (ParametricDistribution distr : distrs) {
	    		if (initial || distr.isDirtyCalculation()) {
	    			// only update at the start and when parametric distributions change
	    			double [] rates = drates[k];
					rates[0] = Double.NEGATIVE_INFINITY;
					rates[1] = distr.inverseCumulativeProbability(0.1/n);
	    			for (int i = 1; i < n1; i++) {
	    				rates[i+1] = distr.inverseCumulativeProbability(i/n);
	    			}
					rates[rates.length - 2] = distr.inverseCumulativeProbability(1.0 - 0.1/n);
					rates[rates.length - 1] = Double.POSITIVE_INFINITY;
	    		}
	    		k++;
	    	}
	    	
	    	// create rate mixture stored in allrates and allweights
	    	double [] allrates = new double[m * rates.length];
	    	double [] allweights = new double[m * rates.length]; // weight of interval
	    	//double [] r = new double[m]; // current rate
	    	int [] j2 = new int[m]; // index of current rate
	    	
	    	// find minimum rate
	    	allrates[0] = drates[0][1];
    		for (int a = 1; a < m; a++) {
    			allrates[0] = Math.min(allrates[0], drates[a][1]);
    		}
	    	for (int i = 1; i < allrates.length; i++) {
	    		double min = drates[0][j2[0]+1];
	    		int minIndex = 0;
	    		for (int a = 1; a < m; a++) {
	    			if (min > drates[a][j2[a]+1]) {
	    				min = drates[a][j2[a]+1];
	    				minIndex = a;
	    			}
	    		}
	    		allrates[i] = min;
	    		for (int a = 0; a < m; a++) {
	    			if (drates[a][j2[a]+1] != drates[a][j2[a]]) {
	    				allweights[i] += weights[a] * (allrates[i] - allrates[i-1]) / (drates[a][j2[a]+1] - drates[a][j2[a]]);
	    			} else {
	    				allweights[i] += weights[a];
	    			}
	    		}
	    		j2[minIndex]++;
	    	}
	    	
	    	// represent rate mixture by piecewise linear approximation
	    	rates[0] = allrates[0];
	    	double target = 0;
	    	double cumProb = 0;
	    	int j3 = 0;
	    	for (int i = 1; i < rates.length; i++) {
	    		target = i / n;
	    		while (cumProb < target) {
	    			cumProb += allweights[j3]/n;
	    			j3++;
		    		if (j3 == allrates.length) {
		    			break;
		    		}
	    		}
	    		double p = cumProb - allweights[j3-1];
	    		if (j3 == allrates.length) {
	    			j3--;
	    		}
	    		rates[i] = allrates[j3-1] + (allrates[j3] - allrates[j3-1]) * (target - p)/(cumProb - p);
	    		if (Double.isNaN(rates[i])) {
	    			System.out.println("Panick: NaN rate " + rates[i]);
	    		}
	    		if (rates[i] < 0) {
	    			System.out.println("Panick: negative rate " + rates[i]);
	    		}
	    	}
	    	rates[rates.length - 1] = allrates[allrates.length - 1];
	    	
	    	
        } catch (MathException e) {
            throw new RuntimeException("Failed to compute inverse cumulative probability!");
        }
    }
    
    @Override
    protected double getMeanWithoutOffset() {
    	throw new IllegalArgumentException("not implemented yet");
    }
    
    @Override
    protected void store() {
        System.arraycopy(rates, 0, storedRates, 0, rates.length);
        for (int i = 0; i < drates.length; i++) {
            System.arraycopy(drates[i], 0, storedDRates[i], 0, drates[i].length);
        }
        super.store();
    }
    
    @Override
    protected void restore() {
        double[] tmp = rates;
        rates = storedRates;
        storedRates = tmp;
        
        double[][] tmp2 = drates;
        drates = storedDRates;
        storedDRates = tmp2;

        super.restore();
    }
    
    @Override
    protected boolean requiresRecalculation() {
    	
    	boolean hasDirtyDistr = false;
    	for (ParametricDistribution distr : distrs) {
    		if (distr.isDirtyCalculation()) {
    			hasDirtyDistr = true;
    			break;
    		}
    	}
    	if (weight.somethingIsDirty() || hasDirtyDistr) {
    		Arrays.fill(rates, 0.0);
    		return true;
    	}
    	
    	return super.requiresRecalculation();
    }

	@Override
	public ContinuousDistribution getDistribution() {
		return mixtureDistr;
	}
	
}
