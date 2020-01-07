package test.beast.math.distributions;

import org.apache.commons.math.MathException;
import org.junit.Test;

import beast.core.parameter.RealParameter;
import beast.math.distributions.MixedPieceWiseLinearDistribution;
import beast.math.distributions.Normal;
import junit.framework.TestCase;

public class MixedPieceWiseLinearDistributionTest extends TestCase {

		@Test
		public void testMixtureOfSameDistributions() throws MathException {
			Normal normal1 = new Normal();
			normal1.initByName("mean", "1.0", "sigma", "1.0");
			Normal normal2 = new Normal();
			normal2.initByName("mean", "1.0", "sigma", "1.0");
			
			RealParameter weights = new RealParameter("0.5 0.5");
			
			MixedPieceWiseLinearDistribution distr = new MixedPieceWiseLinearDistribution();
			distr.initByName("weight", weights, "paramdistr", normal1, "paramdistr", normal2);
			
			double p1 = normal1.cumulativeProbability(1.0);
			double p2 = distr.cumulativeProbability(1.0);
			assertEquals(p1, p2, 1e-2);
		}
	
		@Test
		public void testMixtureOfTwoNormalDistributions() throws MathException {
			Normal normal1 = new Normal();
			normal1.initByName("mean", "1.0", "sigma", "1.0");
			Normal normal2 = new Normal();
			normal2.initByName("mean", "3.0", "sigma", "1.0");
			
			RealParameter weights = new RealParameter("0.5 0.5");
			
			MixedPieceWiseLinearDistribution distr = new MixedPieceWiseLinearDistribution();
			distr.initByName("weight", weights, "paramdistr", normal1, "paramdistr", normal2);
			
			double p2 = distr.cumulativeProbability(2.0);
			assertEquals(0.5, p2, 1e-2);
		}
}
