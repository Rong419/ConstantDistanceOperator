package consoperators;


import org.apache.commons.math.MathException;
import org.apache.commons.math3.util.FastMath;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.branchratemodel.UCRelaxedClockModel;
import beast.base.evolution.operator.TreeOperator;
import beast.base.inference.distribution.LogNormalDistributionModel;
import beast.base.inference.distribution.ParametricDistribution;
import consoperators.distributions.PiecewiseLinearDistribution;

@Description("Operator on tree and rates for relaxed clock")
abstract public class RateAndTreeOperator extends TreeOperator {
    final public Input<UCRelaxedClockModel> clockModelInput = new Input<>("clockModel", "relaxed clock model used to deal with quantiles", Validate.REQUIRED);

    private RealParameter rates;
    protected RealParameter quantiles;
    protected ParametricDistribution distribution;
    protected UCRelaxedClockModel clockModel;

	@Override
	public void initAndValidate() {
        clockModel = clockModelInput.get();
        rates = clockModel.rateInput.get();
        quantiles = clockModel.quantileInput.get();
        if (clockModel != null) {
        	distribution = clockModel.rateDistInput.get();
        }
	}

	protected void setRateValue(int nodeNr, double rate) {
		if (rates != null) {
			rates.setValue(nodeNr, rate);
		} else {			
			try {
				double q = distribution.cumulativeProbability(rate);
				quantiles.setValue(nodeNr, q);
			} catch (MathException e) {
				e.printStackTrace();
			}
		}
		
	}

	protected double getRateValue(int nodeNr) {
		if (rates != null) {
			return rates.getArrayValue(nodeNr);
		} else {
			double q = quantiles.getArrayValue(nodeNr);
			double r = 1;
			try {
				 r = distribution.inverseCumulativeProbability(q);
			} catch (MathException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			return r;
		}
	}

    // step size in numeric approximation of derivatives
	final static private double EPSILON = 1e-8;
    final private static double SQRT2 = Math.sqrt(2.0);

	/**
	 * calculate contribution to HR due to quantile => rate => quantile transformations
	 * Special cases:
	 * 1. rates cached in clock model, if used for piecewise linear approximation for distribution
	 * 2. for log normal, special case is worked out efficiently
	 * 3. for other distributions, numerical approximation is used
	 * 
	 * @param rNew proposed rate
	 * @param qOld original quantile
	 * @param qNew proposed quantile
	 * @return log of HR
	 */
    public double calculateHastingsRatio(double rNew, double qOld, double qNew) {
    	double [] rates = clockModel.getRates();

    	if (rates == null) {
    		// special case for log normal distribution
    		if (distribution instanceof LogNormalDistributionModel) {
    	        double stdev = ((LogNormalDistributionModel) distribution).SParameterInput.get().getArrayValue();
    	        double miu = - 0.5 * stdev * stdev;

    	        double b = FastMath.log(rNew);
    	        double c = 2 * stdev * stdev;
    	        double x = b - miu;
    	        double x_sq = x * x / c;
    	        double rateHR = b + x_sq;

    	        double a = erfInv(2 * qOld - 1);
    	        double quantileHR = miu + SQRT2 * stdev * a + a * a;
    	        return quantileHR - rateHR;
    		}
        	// use numeric approximation of derivative
    		double logHR = 0;
	        try {
	        	double r0 = distribution.inverseCumulativeProbability(qOld);
	        	double r0h = distribution.inverseCumulativeProbability(qOld + EPSILON);
	        	logHR += FastMath.log((r0h - r0) / EPSILON);
	        	
	        	double q0 = distribution.cumulativeProbability(rNew);
        		double q0h = distribution.cumulativeProbability(rNew + EPSILON);
        		logHR += FastMath.log((q0h - q0) / EPSILON);
	        } catch (MathException e) {
	            throw new RuntimeException("Failed to compute inverse cumulative probability!");
	        }
	        return -logHR;
    	}

    	if (distribution instanceof PiecewiseLinearDistribution) {
			PiecewiseLinearDistribution pld = (PiecewiseLinearDistribution) distribution;
			double logHR = Math.log(pld.getDerivativeAtQuantile(qOld));
			logHR +=  Math.log(pld.getDerivativeAtQuantileInverse(rNew, qNew));
			return logHR;
		}
    	
    	throw new IllegalArgumentException("rates should be specified or PiecewiseLinearDistribution should be used");
    }
	
	
	/**
	 * Assumes quantile parameterisation of clockModel, so clockModel must be specified
	 * @param q quantile
	 * @return derivative of rate distribution at quantile q
	 */
	private double getDerivativeAtQuantile(double q, double [] rates) {
        // use cached rates
        double v = q * (rates.length - 1);
        int i = (int) v;
        if (i == rates.length) {
        	int h = 3;
        	h++;
        }
        if (rates[i] == 0.0) {
	        try {
	        	if (i > 0) {
	        		rates[i] = distribution.inverseCumulativeProbability(((double)i) /(rates.length - 1));
	        	} else {
	        		rates[i] = distribution.inverseCumulativeProbability(0.1 / (rates.length - 1));
	        	}
	        } catch (MathException e) {
	            throw new RuntimeException("Failed to compute inverse cumulative probability!");
	        }
        }
        if (i < rates.length - 1 && rates[i + 1] == 0.0) {
	        try {
	        	if (i < rates.length - 2) {
	        		rates[i + 1] = distribution.inverseCumulativeProbability(((double)(i + 1)) / (rates.length - 1));
	        	} else {
	        		rates[i + 1] = distribution.inverseCumulativeProbability((rates.length - 1 - 0.1) / (rates.length - 1));
	        	}
	        } catch (MathException e) {
	            throw new RuntimeException("Failed to compute inverse cumulative probability!");
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
	 * @return derivative of quantile distribution at rate r
	 */
    private double getDerivativeAtQuantileInverse(double r, double qNew, double [] rates) {
    	// find quantile for new rate
    	// TODO: pass this as argument to this method, since it will be required for setting quantiles anyway?
    	
        double v = qNew * (rates.length - 1);
        int i = (int) v;
        if (rates[i] == 0.0) {
	        try {
	        	if (i > 0) {
	        		rates[i] = distribution.inverseCumulativeProbability(((double)i) / (rates.length - 1));
	        	} else {
	        		rates[i] = distribution.inverseCumulativeProbability(0.1 / (rates.length - 1));
	        	}
	        } catch (MathException e) {
	            throw new RuntimeException("Failed to compute inverse cumulative probability!");
	        }
        }
        if (i < rates.length - 1 && rates[i + 1] == 0.0) {
	        try {
	        	if (i < rates.length - 2) {
	        		rates[i + 1] = distribution.inverseCumulativeProbability(((double)(i + 1)) / (rates.length - 1));
	        	} else {
	        		rates[i + 1] = distribution.inverseCumulativeProbability((rates.length - 1 - 0.1) / (rates.length - 1));
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
    	        		rates[i] = distribution.inverseCumulativeProbability(((double)i) / (rates.length - 1));
    	        	} else {
    	        		rates[i] = distribution.inverseCumulativeProbability(0.1 / (rates.length - 1));
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
    	        		rates[i + 1] = distribution.inverseCumulativeProbability(((double)(i + 1)) / (rates.length - 1));
    	        	} else {
    	        		rates[i + 1] = distribution.inverseCumulativeProbability((rates.length - 0.1) / (rates.length - 1));
    	        	}
    	        } catch (MathException e) {
    	            throw new RuntimeException("Failed to compute inverse cumulative probability!");
    	        }
            }
        }

        if (i < rates.length - 1) {
            double derivative = (1.0/(rates.length - 1))/(rates[i+1] - rates[i]);
            return derivative;        	
        }
        return 0;
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
