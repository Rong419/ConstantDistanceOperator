package beast.math.distributions;

import java.util.Arrays;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.ContinuousDistribution;
import org.apache.commons.math.distribution.Distribution;
import org.apache.commons.math.distribution.NormalDistributionImpl;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;



@Description("Approximates parametric distribution by piecewise linear approximation.")
public class PiecewiseLinearDistribution extends ParametricDistribution {
    final public Input<ParametricDistribution> distrInput = new Input<>("distr", "Underlying parametric distribution that is approximated.", Validate.REQUIRED);
    final public Input<Integer> numberOfDiscreteRates = new Input<>("bins", "number of bins used to approximate the distribution piecewise linearly. (default 100)", 100);
    
    final public Input<Double> limitInput = new Input<>("limit", "fraction of bins at end of distribution to be cut off -- should be between 0 and 1", 0.1);
    final public Input<Boolean> cutOffEndInput = new Input<>("cutOffEnd", "approximate values below for quantiles close to 0 and 1 with rate at limit/bins and (1-limit/bins) respectively."
    		+ "If false, for quantiles below limit/bins and above (1-limit/bins) the inverseCumulativeProbability methods of the underlying parametric distribution is called"
    		+ "(which is slower, but should not happen very often).", true);

    LinearPiecewiseImpl dist = new LinearPiecewiseImpl();
    ContinuousDistribution underlyingDistr = new NormalDistributionImpl();
    
    // ParametricDistribution distribution;
    
    protected double[] rates; //the output rates
    protected double[] storedRates; //
    
    private double limit_0 = 0.1; // limit_0 / (numberOfDiscreteRates-1) is the minimum quantile
    private double limitLow, limitUp;
    private boolean cutOffEnd;


    public PiecewiseLinearDistribution() {}
    
    /** construct log-normal distribution with mean = 1 **/ 
    public PiecewiseLinearDistribution(RealParameter S) {
    	LogNormalDistributionModel lnormal = new LogNormalDistributionModel();
    	lnormal.initByName("S", S, "M", "1.0", "meanInRealSpace", true);
    	initByName("distr", lnormal);
    }
    
    @Override
	public void initAndValidate() {
    	int dim = numberOfDiscreteRates.get();
    	rates = new double[dim];
    	storedRates = new double[dim];
    	ParametricDistribution distribution = distrInput.get();
    	Distribution d = distribution.getDistribution();
    	if (d instanceof ContinuousDistribution) {
    		underlyingDistr = (ContinuousDistribution) d;
    	} else {
    		throw new IllegalArgumentException("Expected parametric distribution that is continuous");
    	}
    	limit_0 = limitInput.get();
    	if (limit_0 <= 0 || limit_0 >= 1) {
    		throw new IllegalArgumentException("limit should be between 0 and 1");
    	}
    	limitLow = limit_0 / dim;
    	limitUp = 1.0 - limit_0 / dim;
    	cutOffEnd = cutOffEndInput.get();
        refresh();
    }

    /**
     * make sure internal state is up to date *
     */
    void refresh() {
    }

    @Override
    public Distribution getDistribution() {
        refresh();
        return dist;
    }

    public class LinearPiecewiseImpl implements ContinuousDistribution {

        public LinearPiecewiseImpl() {
        }


        @Override
        public double cumulativeProbability(double x) throws MathException {
        	// Return exact cdf using piecewise linear approximation
            int i = getIntervalFor(x);

            if (!cutOffEnd) {
            	// TODO: needs testing
            	if (i == 0) {
            		if (x < rates[0]) {
                		return underlyingDistr.cumulativeProbability(x);
            		}
                	double cdf = (limitLow + (1-limit_0) * (x - rates[i]) / (rates[i+1] - rates[i])) / (rates.length-1);
            		return cdf;
            	} else if (i == rates.length-1) {
            		return underlyingDistr.cumulativeProbability(x);
            	} else if (i == rates.length-2) {
                	double cdf = (i + (1-limit_0) * (x - rates[i]) / (rates[i+1] - rates[i])) / (rates.length-1);
            		return cdf;
            	}
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
        	if (!cutOffEnd && (q < limitLow  || q > limitUp)) {
        		return underlyingDistr.inverseCumulativeProbability(q);
        	}
            double v = q * (rates.length - 1);
            int i = (int) v;
            
            // make sure cached rates are calculated
            if (rates[i] == 0.0) {
    	        try {
    	        	if (i > 0) {
    	        		rates[i] = underlyingDistr.inverseCumulativeProbability(((double)i) / (rates.length-1));
    	        	} else {
    	        		rates[i] = underlyingDistr.inverseCumulativeProbability(limit_0 / (rates.length-1));
    	        	}
    	        } catch (MathException e) {
    	            throw new RuntimeException("Failed to compute inverse cumulative probability!");
    	        }
            }
            if (i < rates.length - 1 && rates[i + 1] == 0.0) {
    	        try {
    	        	if (i < rates.length - 2) {
    	        		rates[i + 1] = underlyingDistr.inverseCumulativeProbability(((double)(i + 1)) / (rates.length-1));
    	        	} else {
    	        		rates[i + 1] = underlyingDistr.inverseCumulativeProbability((rates.length - 1 - limit_0) / (rates.length-1));
    	        	}
    	        } catch (MathException e) {
    	            throw new RuntimeException("Failed to compute inverse cumulative probability!");
    	        }
            }
            
            if (cutOffEnd) {
            	if (i == 0) {
            		v = (v-limit_0) / (1-limit_0);
            	} else if (i == rates.length-1) {
            		v = i + (v - i) / (1-limit_0);
            	}
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


        protected int getIntervalFor(double x) throws MathException {
            double q = underlyingDistr.cumulativeProbability(x);
            int i = getIntervalFor(x, q);
			return i;
		}

        protected int getIntervalFor(double r, double qNew) {
	        double v = qNew * (rates.length - 1);
	        int i = (int) v;
	        if (rates[i] == 0.0) {
		        try {
		        	if (i > 0) {
		        		rates[i] = underlyingDistr.inverseCumulativeProbability(((double)i) / (rates.length - 1));
		        	} else {
		        		rates[i] = underlyingDistr.inverseCumulativeProbability(limit_0 / (rates.length - 1));
		        	}
		        } catch (MathException e) {
		            throw new RuntimeException("Failed to compute inverse cumulative probability!");
		        }
	        }
	        if (i < rates.length - 1 && rates[i + 1] == 0.0) {
		        try {
		        	if (i < rates.length - 2) {
		        		rates[i + 1] = underlyingDistr.inverseCumulativeProbability(((double)(i + 1)) / (rates.length - 1));
		        	} else {
		        		rates[i + 1] = underlyingDistr.inverseCumulativeProbability((rates.length - 1 - limit_0) / (rates.length - 1));
		        	}
		        } catch (MathException e) {
		            throw new RuntimeException("Failed to compute inverse cumulative probability!");
		        }
	        }
	        
	        // test boundary: r should be between rates[i] and rates[i+1]
	        // but due to numerical errors in the piecewise linear approximation could
	        // fall just outside, so make sure the condition is met, and rates are calculated
	        while (i > 0 && r < rates[i]) {
	        	i--;
	            if (i > 0 && rates[i] == 0.0) {
	    	        try {
	    	        	if (i > 0) {
	    	        		rates[i] = underlyingDistr.inverseCumulativeProbability(((double)i) / (rates.length - 1));
	    	        	} else {
	    	        		rates[i] = underlyingDistr.inverseCumulativeProbability(limit_0 / (rates.length - 1));
	    	        	}
	    	        } catch (MathException e) {
	    	            throw new RuntimeException("Failed to compute inverse cumulative probability!");
	    	        }
	            }
	        }
	        while (i < rates.length - 1 && r > rates[i+1]) {
	        	i++;
	            if (i < rates.length - 1 && rates[i + 1] == 0.0) {
	    	        try {
	    	        	if (i < rates.length - 2) {
	    	        		rates[i + 1] = underlyingDistr.inverseCumulativeProbability(((double)(i + 1)) / (rates.length - 1));
	    	        	} else {
	    	        		rates[i + 1] = underlyingDistr.inverseCumulativeProbability((rates.length - limit_0) / (rates.length - 1));
	    	        	}
	    	        } catch (MathException e) {
	    	            throw new RuntimeException("Failed to compute inverse cumulative probability!");
	    	        }
	            }
	        }
	        return i;
		}
    } // class LinearPiecewiseImpl

    @Override
    protected double getMeanWithoutOffset() {
    	throw new IllegalArgumentException("not implemented yet");
    }
    
    @Override
    protected void store() {
        System.arraycopy(rates, 0, storedRates, 0, rates.length);
        super.store();
    }
    
    @Override
    protected void restore() {
        double[] tmp = rates;
        rates = storedRates;
        storedRates = tmp;
        
        underlyingDistr = getUnderlyingDistr();
        
    	super.restore();
    }

    protected ContinuousDistribution getUnderlyingDistr() {
    	return (ContinuousDistribution) distrInput.get().getDistribution();
	}

	@Override
    protected boolean requiresRecalculation() {
    	
    	if (distrInput.get() != null && distrInput.get().isDirtyCalculation()) {
    		Arrays.fill(rates, 0.0);
    		underlyingDistr = getUnderlyingDistr();
    		return true;
    	}
    	
    	return super.requiresRecalculation();
    }
    
	/**
	 * Assumes quantile parameterisation of clockModel, so clockModel must be specified
	 * @param q quantile
	 * @return derivative of rate distribution at quantile q
	 */
	final static private double EPSILON = 1e-8;
	public double getDerivativeAtQuantile(double q) {
    	if (!cutOffEnd && (q < limitLow  || q > limitUp)) {
        	// TODO: needs testing
			try {
	    		double r = underlyingDistr.inverseCumulativeProbability(q);
	    		double rPlusH = underlyingDistr.inverseCumulativeProbability(q + EPSILON);
	    		double dR = (rPlusH - r) / EPSILON;
	    		return dR;
			} catch (MathException e) {
	            throw new RuntimeException("Failed to compute inverse cumulative probability!" + e.getMessage());
			}
    	}

		// use cached rates
        double v = q * (rates.length - 1);
        int i = (int) v;
        if (rates[i] == 0.0) {
	        try {
	        	if (i > 0) {
	        		rates[i] = underlyingDistr.inverseCumulativeProbability(((double)i) /(rates.length - 1));
	        	} else {
	        		rates[i] = underlyingDistr.inverseCumulativeProbability(limit_0 / (rates.length - 1));
	        	}
	        } catch (MathException e) {
	            throw new RuntimeException("Failed to compute inverse cumulative probability!");
	        }
        }
        if (i < rates.length - 1 && rates[i + 1] == 0.0) {
	        try {
	        	if (i < rates.length - 2) {
	        		rates[i + 1] = underlyingDistr.inverseCumulativeProbability(((double)(i + 1)) / (rates.length - 1));
	        	} else {
	        		rates[i + 1] = underlyingDistr.inverseCumulativeProbability((rates.length - 1 - limit_0) / (rates.length - 1));
	        	}
	        } catch (MathException e) {
	            throw new RuntimeException("Failed to compute inverse cumulative probability!");
	        }
        }
        if (!cutOffEnd) {
        	if (i == 0 || i == rates.length - 2) {
            	double r = (rates[i+1] - rates[i]) / ((1.0 - limit_0)/(rates.length - 1));
            	return r;        		
            }
        }
        if (i < rates.length - 1) {
        	double r = (rates[i+1] - rates[i]) / (1.0/(rates.length - 1));
        	return r;
        }
        return 0;
    }
    
	/**
	 * Assumes quantile parameterisation of clockModel, so clockModel must be specified
	 * @param r rate
	 * @param qNew: quantile associated with rate r
	 * @return derivative of quantile distribution at rate r
	 */
    public double getDerivativeAtQuantileInverse(double r, double qNew) {
    	// TODO: verify the following is correct
    	// return 1.0 / getDerivativeAtQuantile(qNew);
    	
    	// TODO: implement !cutOffEnd version
    	
    	int i = dist.getIntervalFor(r, qNew);
        if (i < rates.length - 1) {
            double derivative = (1.0/(rates.length - 1))/(rates[i+1] - rates[i]);
            return derivative;        	
        }
        return 0;
    }
    
    
    
    /**
     * @return The minimum rate supported by the function 
     * @throws MathException 
     */
    public double getRangeMin() throws MathException {
    	if (!cutOffEnd) {
    		return underlyingDistr.inverseCumulativeProbability(0.0);
    	}
    	if (rates[0] == 0.0) {
    		rates[0] = underlyingDistr.inverseCumulativeProbability(limit_0 / (rates.length - 1));
    	}
    	return rates[0];
    	//return getDerivativeAtQuantile(limit_0 / (rates.length - 1));
    }
    
    
    /**
     * @return The maximum rate supported by the function 
     * @throws MathException 
     */
    public double getRangeMax() throws MathException {
    	if (!cutOffEnd) {
    		return underlyingDistr.inverseCumulativeProbability(1.0);
    	}
    	if (rates[rates.length-1] == 0.0) {
    		rates[rates.length-1] = underlyingDistr.inverseCumulativeProbability((rates.length - limit_0 - 1) / (rates.length - 1));
    	}
    	return rates[rates.length-1];
    	//return getDerivativeAtQuantile((rates.length - 1 - limit_0) / (rates.length - 1));
    }
    
    
    
    
 }
