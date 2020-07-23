package beast.math.distributions;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.ContinuousDistribution;
import org.apache.commons.math.distribution.Distribution;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;

@Description("Caches parametric distribution for inverscumulative methods.")
public class CachedDistribution extends ParametricDistribution {
    final public Input<ParametricDistribution> distrInput = new Input<>("distr", "Underlying parametric distribution that is approximated.", Validate.REQUIRED);
    
    ContinuousDistribution dist;
    
    // ParametricDistribution distribution;
    
    Map<Double, Double> cache = new HashMap<>(); 
    Map<Double, Double> storedcache = new HashMap<>(); 
    
    @Override
	public void initAndValidate() {
    	dist = new CachedImpl();
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

    public class CachedImpl implements ContinuousDistribution {
    	ContinuousDistribution underlyingDistr;
    	
        public CachedImpl() {
        	underlyingDistr = getUnderlyingDistr();
        }

        @Override
        public double cumulativeProbability(double x) throws MathException {
        	// System.err.println(x);
        	if (x <= 0) {
        		return 0;
        	}
        	if (Double.isInfinite(x)) {
        		return 1;
        	}
        	return underlyingDistr.cumulativeProbability(x);
        }
        
		@Override
        public double cumulativeProbability(double x0, double x1) throws MathException {
            return cumulativeProbability(x1) - cumulativeProbability(x0);
        }

        @Override
        public double inverseCumulativeProbability(double q) throws MathException {
        	if (q >= 1) {
        		return Double.POSITIVE_INFINITY;
        	}
        	if (Double.isNaN(q)) {
        		return 0;
        	}
        	final Double qD = q;
        	if (cache.containsKey(qD)) {
        		return cache.get(qD);
        	}
        	double x = underlyingDistr.inverseCumulativeProbability(q);
        	cache.put(qD, x);
        	return x;
//        	return underlyingDistr.inverseCumulativeProbability(q);
        }

        @Override
        public double density(double x) {
        	return underlyingDistr.density(x);
        }

        @Override
        public double logDensity(double x) {
        	return Math.log(density(x));
        }
    } // class CachedImpl

    @Override
    protected double getMeanWithoutOffset() {
    	throw new IllegalArgumentException("not implemented yet");
    }
    
    @Override
    protected void store() {
    	storedcache.clear();
    	for (Double d : cache.values()) {
    		storedcache.put(d, cache.get(d));
    	}
        super.store();
    }
    
    @Override
    protected void restore() {
    	super.restore();
    }

    protected ContinuousDistribution getUnderlyingDistr() {
    	return (ContinuousDistribution) distrInput.get().getDistribution();
	}

	@Override
    protected boolean requiresRecalculation() {
    	
    	if (distrInput.get() != null && distrInput.get().isDirtyCalculation()) {
    		cache.clear();
    		return true;
    	}
    	
    	return super.requiresRecalculation();
    }
    
 }
