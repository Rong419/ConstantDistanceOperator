package beast.math.distributions;



import java.util.Arrays;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.*;
import org.apache.commons.math.special.Gamma;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.math.distributions.ParametricDistribution;

@Description("Distribution with mean one gouverned by a single variance parameter")
public class OneParameterMeanOneDistribution extends ParametricDistribution {
    final public Input<RealParameter> shapeInput = new Input<>("sigma", "variance parameter that gouverns the distribution");
    
    enum Mode {weibull, gamma, invgamma, lognormal};
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
		case invgamma:
			dist = new InverseGammaDistributionImpl(2.0, 1.0);
			break;
		case lognormal:
			dist = new LogNormalImpl(1.0, 1.0);
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
            return Math.exp(m_normal.inverseCumulativeProbability(p));
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
		switch (modeInput.get()) {
		case weibull:
			double shape = ((WeibullDistributionImpl)dist).getShape();
			double scale = ((WeibullDistributionImpl)dist).getScale();
			return scale * org.apache.commons.math3.special.Gamma.gamma(1 + 1.0 / shape);
		case gamma:
			double alpha = ((GammaDistributionImpl)dist).getAlpha();
			double beta = ((GammaDistributionImpl)dist).getBeta();
			return alpha * beta;
		case invgamma:
			double a = ((InverseGammaDistributionImpl)dist).alpha;
			double b = ((InverseGammaDistributionImpl)dist).beta;
			return b / (a-1);
		case lognormal:
			double m = ((LogNormalImpl)dist).m_fMean;
			double s = ((LogNormalImpl)dist).m_fStdDev;
    		return Math.exp(m + s * s/2.0);
		}
		return 0;
	}
	
}

