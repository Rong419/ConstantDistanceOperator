package test.beast.math.distributions;

import org.apache.commons.math.MathException;
import org.junit.Test;

import beast.base.core.Log;
import beast.base.inference.distribution.LogNormalDistributionModel;
import consoperators.distributions.OneParameterMeanOneDistribution;
import junit.framework.TestCase;

public class OneParameterMeanOneDistributionTest extends TestCase {
	
	@Test
	public void testOneParameterLogNormal() throws MathException {
		// implementation of inverseCumulativeProbability optimised for speed
		OneParameterMeanOneDistribution p = new OneParameterMeanOneDistribution();
		p.initByName("mode", "lognormal", "sigma", "1.0");
		
		// the standard, but slow implementation
		LogNormalDistributionModel base = new LogNormalDistributionModel();
		double s = Math.sqrt(Math.log(1.0 + 1.0));
		double m =  - (0.5 * s * s);
		base.initByName("M", m + "", "S", s + "");
		
		double q = 0.01;
		while (q < 1.0) {
			double actual = p.inverseCumulativeProbability(q);
			double expected = base.inverseCumulativeProbability(q);
			assertEquals(expected, actual, 1e-6);
			q += 0.01;
		}
	}

}
