package beast.math.distributions;



import java.util.Arrays;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.*;
import org.apache.commons.math.special.Gamma;
import org.apache.commons.math3.util.FastMath;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.math.distributions.ParametricDistribution;

@Description("Distribution with mean one gouverned by a single variance parameter")
public class OneParameterMeanOneDistribution extends ParametricDistribution {
    final public Input<RealParameter> shapeInput = new Input<>("sigma", "variance parameter that governs the distribution");
    
    public enum Mode {weibull, gamma, invgamma, lognormal, exp, strict};
    final public Input<Mode> modeInput = new Input<>("mode", "which distribution to use. One of " + Arrays.toString(Mode.values()), 
    		Mode.weibull, Mode.values());

    ContinuousDistribution dist;
    RealParameter shape;
    
    double [] weibullCache;
    
	@Override
	public void initAndValidate() {
		shape = shapeInput.get();
		switch (modeInput.get()) {
		case weibull:
			dist = new WeibullDistributionImpl(1.0, 1.0);
			break;
		case gamma:
			dist = new GammaDistributionImpl(1.0, 1.0);
			break;
		case exp:
			dist = new ExponentialDistributionImpl(1.0);
			break;
		case invgamma:
			dist = new InverseGammaDistributionImpl(2.0, 1.0);
			break;
		case lognormal:
			dist = new LogNormalImpl(1.0, 1.0);
			break;
		case strict:
			dist = new DeltaDistrImpl();
			break;
		}
		
		// stores S in cache for values of k from 0.0 to 10.0
		weibullCache = new double[1000];
		for (int i = 0; i < 1000; i++) {
			double k = (i+2) / 100.0;
			double x = org.apache.commons.math3.special.Gamma.gamma(1.0 + 1.0 / k);
			weibullCache[weibullCache.length - i - 1] = org.apache.commons.math3.special.Gamma.gamma(1.0 + 2.0 / k) / (x*x) - 1;
		}
		
		refresh();
	}
	
	/** sets other parameters of the distribution in order to maintain a mean of 1 **/
	@SuppressWarnings("deprecation")
	protected void refresh() {
		double S = this.shape.getValue();
		
		switch (modeInput.get()) {
		case weibull:
			double k = 0;
			if (S > weibullCache[weibullCache.length - 1]) {
				k = 10;
			} else {
				int i = Arrays.binarySearch(weibullCache, S);
				if (i > 0) {
					k = 10.0 - i/100.0;
				} else {
					i = -i-1;
					k = 10.0 - (i/100.0 + 0.01 * (S - weibullCache[i]) / (weibullCache[i+1] - weibullCache[i]));
				}
			}
			double scale = 1.0 / org.apache.commons.math3.special.Gamma.gamma(1.0 + 1.0/k);
			((WeibullDistributionImpl)dist).setShape(k);
			((WeibullDistributionImpl)dist).setScale(scale);
			break;
		case gamma:
			double alphaG = 1.0/S;
			double betaG = S;
			((GammaDistributionImpl)dist).setAlpha(alphaG);
			((GammaDistributionImpl)dist).setBeta(betaG);
			break;
		case exp:
			break;
		case invgamma:
			double alphaIG = 2.0 + 1.0 / S;
			double betaIG = alphaIG - 1.0;
			((InverseGammaDistributionImpl)dist).setAlphaBeta(alphaIG, betaIG);
			break;
		case lognormal:
			double s = Math.sqrt(Math.log(S + 1.0));
			double m =  - (0.5 * s * s);
			((LogNormalImpl)dist).setMeanAndStdDev(m, s);
			break;
		case strict:
			break;
		}
	}

	public void setShape(double value) {
		shape.setValue(value);
		refresh();
	}
	
	public double cumulativeProbability(double x) throws MathException {
		return dist.cumulativeProbability(x);
	}
	
	public double inverseCumulativeProbability(double p) throws MathException {
		return dist.inverseCumulativeProbability(p);
	}
	
		
    public class DeltaDistrImpl implements ContinuousDistribution {
		@Override
		public double cumulativeProbability(double x) throws MathException {
			if (x < 1)
				return 0;
			return 1;
		}

		@Override
		public double cumulativeProbability(double x0, double x1) throws MathException {
			return cumulativeProbability(x1) - cumulativeProbability(x0);
		}

		@Override
		public double inverseCumulativeProbability(double p) throws MathException {
			return 1;
		}

		@Override
		public double density(double x) {
			if (x == 1) {
				return 1;
			}		
			return 0;
		}

		@Override
		public double logDensity(double x) {
			return Math.log(density(x));
		}
    	
    }
    
    public class LogNormalImpl implements ContinuousDistribution {
        double m_fMean;
        double m_fStdDev;
        NormalDistributionImpl m_normal = new NormalDistributionImpl(0, 1);

        public LogNormalImpl(double mean, double stdDev) {
            setMeanAndStdDev(mean, stdDev);
        }

        @SuppressWarnings("deprecation")
		void setMeanAndStdDev(double mean, double stdDev) {
            m_fMean = mean;
            m_fStdDev = stdDev;
            m_normal.setMean(mean);
            m_normal.setStandardDeviation(stdDev);
        }

        @Override
        public double cumulativeProbability(double x) throws MathException {
            return m_normal.cumulativeProbability(Math.log(x));
        }

        @Override
        public double cumulativeProbability(double x0, double x1) throws MathException {
            return cumulativeProbability(x1) - cumulativeProbability(x0);
        }

        @Override
        public double inverseCumulativeProbability(double p) throws MathException {
        	return Math.exp(m_fMean + Math.sqrt(2*m_fStdDev*m_fStdDev) * erfInv(2*p-1)); 
            // return Math.exp(m_normal.inverseCumulativeProbability(p));
        }

        private double erfInv(final double x) {
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


        @Override
        public double density(double x) {
            if( x <= 0 ) {
                return 0;
            }
            return m_normal.density(Math.log(x)) / x;
        }

        @Override
        public double logDensity(double x) {
            if( x <= 0 ) {
                return  Double.NEGATIVE_INFINITY;
            }
            return m_normal.logDensity(Math.log(x)) - Math.log(x);
        }
    } // class LogNormalImpl

    public class InverseGammaDistributionImpl extends AbstractContinuousDistribution {
		private static final long serialVersionUID = 1L;
		double alpha;
        double beta;
        // log of the constant beta^alpha/Gamma(alpha)
        double C;

        InverseGammaDistributionImpl(double alpha, double beta) {
            setAlphaBeta(alpha, beta);
        }

        public void setAlphaBeta(double alpha, double beta) {
            this.alpha = alpha;
            this.beta = beta;
            C = alpha * Math.log(this.beta) - org.apache.commons.math.special.Gamma.logGamma(alpha);
        }

        @Override
        public double cumulativeProbability(double x) throws MathException {
        	return Gamma.regularizedGammaP(alpha, beta/x) / Gamma.digamma(alpha); 
        }

        @Override
        public double inverseCumulativeProbability(double p) throws MathException {
            if (p == 0) {
                return 0d;
            }
            if (p == 1) {
                return Double.POSITIVE_INFINITY;
            }
            // this is what R thinks
            final org.apache.commons.math.distribution.GammaDistribution g = new GammaDistributionImpl(alpha, beta);
            return 1/g.inverseCumulativeProbability(1-p);
            // return super.inverseCumulativeProbability(p);
        }


        
        @Override
        public double density(double x) {
            double logP = logDensity(x);
            return Math.exp(logP);
        }

        @Override
        public double logDensity(double x) {
            double logP = -(alpha + 1.0) * Math.log(x) - (beta / x) + C;
            return logP;
        }

		@Override
		protected double getInitialDomain(double p) {
			return 1.0;
		}

		@Override
		protected double getDomainLowerBound(double p) {
			return 0;
		}

		@Override
		protected double getDomainUpperBound(double p) {
			return 1000;
		}
    } // class InverseGammaImpl

	@Override
	public Distribution getDistribution() {		
		return dist;
	}

	@Override
	protected double getMeanWithoutOffset() {
		// OneParameterMeanOneDistribution should have mean = 1...
		return 1;
	}
	
	@Override
	protected void restore() {
		requiresRecalculation(); // refreshes for any distribution not exp or strict
	}

	
	@Override
	protected boolean requiresRecalculation() {
		switch (modeInput.get()) {
		case exp:
		case strict:
			return false;
		default:
			refresh();
			return true;
		}
	}
}

