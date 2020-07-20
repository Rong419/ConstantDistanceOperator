package test.beast.math.distributions;

import org.apache.commons.math.MathException;
import org.junit.Test;

import beast.math.distributions.LogNormalDistributionModel;
import beast.math.distributions.ParametricDistribution;
import beast.math.distributions.PiecewiseLinearDistribution;
import junit.framework.TestCase;

public class PieceWiseLinearDistributionTest extends TestCase {
	
	@Test
	public void testDerivative() throws MathException {
		testDerivative(false);
		//testDerivative(true);
	}
	
	public void testDerivative(boolean cutoff) throws MathException {

		ParametricDistribution distr = new LogNormalDistributionModel();
		distr.initByName("M", "1.0", "S", "1.0", "meanInRealSpace", true);
		
		PiecewiseLinearDistribution pwld = new PiecewiseLinearDistribution();
		pwld.initByName("distr", distr, "cutOffEnd", cutoff);

		double x = 0.00001;
		while (x < 0.01) {
			System.out.println(x + " " + distr.inverseCumulativeProbability(x) + " " + pwld.inverseCumulativeProbability(x));
			x += 0.0005;
		}
		while (x < 0.99) {
			System.out.println(x + " " + distr.inverseCumulativeProbability(x) + " " + pwld.inverseCumulativeProbability(x));
			x += 0.01;
		}
		while (x < 1) {
			System.out.println(x + " " + distr.inverseCumulativeProbability(x) + " " + pwld.inverseCumulativeProbability(x));
			x += 0.0005;
		}
		
		
		
		testDerivativeAt(0.5, distr, pwld, 1e-3);		

		testDerivativeAt(0.25, distr, pwld, 1e-2);
		testDerivativeAt(0.75, distr, pwld, 7e-2);
		testDerivativeAt(0.1, distr, pwld, 1e-2);
		testDerivativeAt(0.9, distr, pwld, 1e-0);
		testDerivativeAt(0.01, distr, pwld, 1.3);
		testDerivativeAt(0.0011, distr, pwld, 4.5);

		// cases outside cut-off
		testDerivativeAt(0.0001, distr, pwld, 1e-2);
		testDerivativeAt(0.9999, distr, pwld, 4e-1);

		testDerivativeAt(0.989, distr, pwld, 550);
		testDerivativeAt(0.9991, distr, pwld, 550);
	}

	private void testDerivativeAt(double q, ParametricDistribution distr, PiecewiseLinearDistribution pwld, double accuracy) throws MathException {
		double h = 1e-8;
		double x = distr.inverseCumulativeProbability(q);
		double xh = distr.inverseCumulativeProbability(q+h);
		double derivativeE = (xh-x)/h;
		
		double derivativeA = pwld.getDerivativeAtQuantile(q);
		
		assertEquals(derivativeE, derivativeA, accuracy);
	}

	@Test
	public void testInverseDerivative() throws MathException {
		testInverseDerivative(false);		
		// testInverseDerivative(true);
	}
	
	public void testInverseDerivative(boolean cutoff) throws MathException {
		ParametricDistribution distr = new LogNormalDistributionModel();
		distr.initByName("M", "1.0", "S", "1.0", "meanInRealSpace", true);
		
		PiecewiseLinearDistribution pwld = new PiecewiseLinearDistribution();
		pwld.initByName("distr", distr, "cutOffEnd", cutoff);

		
		testInverseDerivativeAt(0.5, distr, pwld, 1e-3);
		testInverseDerivativeAt(0.01, distr, pwld, 1.3);
		testInverseDerivativeAt(0.25, distr, pwld, 1e-2);
		testInverseDerivativeAt(0.1, distr, pwld, 1e-2);
		testInverseDerivativeAt(0.75, distr, pwld, 7e-2);
		testInverseDerivativeAt(0.9, distr, pwld, 1e-0);

		// cases outside cut-off
		testInverseDerivativeAt(0.0001, distr, pwld, 1e-5);
		testInverseDerivativeAt(0.9999, distr, pwld, 1e-5);

		testInverseDerivativeAt(0.9991, distr, pwld, 1e-3);
		testInverseDerivativeAt(0.0011, distr, pwld, 4.5);
		testInverseDerivativeAt(0.99, distr, pwld, 550);
	}
	
	private void testInverseDerivativeAt(double q, ParametricDistribution distr, PiecewiseLinearDistribution pwld, double accuracy) throws MathException {
		double h = 1e-8;
		double x = distr.inverseCumulativeProbability(q);
		double xh = distr.inverseCumulativeProbability(q+h);
		double derivativeE = 1.0/((xh-x)/h);
		
		double derivativeA = pwld.getDerivativeAtQuantileInverse(x, q);
		
		assertEquals(derivativeE, derivativeA, accuracy);
	}

}
