package test;

import java.text.DecimalFormat;

import consoperators.ConsOperatorUtils;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.interpolation.*;
import org.apache.commons.math3.util.FastMath;
import org.junit.Test;

import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.branchratemodel.UCRelaxedClockModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.distribution.LogNormalDistributionModel;
import beast.base.util.Randomizer;
import consoperators.InConstantDistanceOperator;
import consoperators.RateAndTreeOperator;
import consoperators.distributions.PiecewiseLinearDistribution;
import junit.framework.TestCase;
import test.beast.BEASTTestCase;

public class RateAndTreeOperatorTest extends TestCase {
	
	// check that HR = 0 when rate == invCummulative(quantile)
	@Test
	public void testHR() throws Exception {
//		int n = 10000000;
//		long start = System.currentTimeMillis();
//		double sum = 0;
//		for (int i = 0; i < n; i++) {
//			sum += Math.exp(Randomizer.nextDouble());
//		}
//		long end = System.currentTimeMillis();
//		System.out.println(sum + " " + (start-end) + " ms");
//		sum = 0;
//		for (int i = 0; i < n; i++) {
//			sum += FastMath.exp(Randomizer.nextDouble());
//		}
//		long end2 = System.currentTimeMillis();
//		System.out.println(sum + " " + (end-end2) + " ms");
		
		
		
		double stdev0 = 0.5;
        Alignment data = BEASTTestCase.getAlignment();
        Tree tree = BEASTTestCase.getTree(data);
        RealParameter quantiles = new RealParameter("0.0");

        RealParameter stdev = new RealParameter(stdev0 + "");
        LogNormalDistributionModel distr = new LogNormalDistributionModel();
        distr.initByName("M", "1.0", "meanInRealSpace", true, "S", stdev);
    
        UCRelaxedClockModel clockModel = new UCRelaxedClockModel();
        clockModel.initByName("tree", tree, "distr", distr, "clock.rate", "1.0", "rateQuantiles", quantiles);//, "numberOfDiscreteRates", 100);

        //double q = 0.5/100;
        double q = 0.01;
        Node node = tree.getNode(0);
        while (q < 1) {
        	quantiles.setValue(0, q);
        	double rate = clockModel.getRateForBranch(node);
        
        	double HR = ConsOperatorUtils.calculateHastingsRatio(rate, q, stdev0);
        	assertEquals(0, HR, 1e-9);
        	q = q * 2;
        }
    }	
	
	
	public void testPWLA() throws Exception {
		double stdev = 0.25;
        PiecewiseLinearDistribution distr = new PiecewiseLinearDistribution();
        LogNormalDistributionModel lognormal = new LogNormalDistributionModel();
        lognormal.initByName("M", "1.0", "meanInRealSpace", true, "S", stdev + "");
        distr.initByName("distr", lognormal);
        
        for (int i = 0; i < 1000; i++) {
        	double q = Randomizer.nextDouble();
        	double q2 = distr.cumulativeProbability(distr.inverseCumulativeProbability(q));
        	assertEquals(q, q2, 1e-9);
        }
	}

	
	final static int N = 32;
	// CubicSplineFast cubicApprox = new CubicSplineFast(N);
	// QuadraticInterpolation quadraticApprox = new QuadraticInterpolation(N);
	

	/** compare log normal distribution with piecewise linear approximation in UCRelaxedClock **/
	@Test
	public void testTolerance() throws Exception {
		double stdev = 0.01;
		while (stdev < 1.5) {
			double maxRelErr = testTolerance(stdev);
			// allow up to 5% error
			assertTrue(maxRelErr < 0.05);
			System.err.println("Stdev:" + stdev + " error: " + maxRelErr);
			stdev *= 2;
		}
	}
	
	private double testTolerance(double stdev0) throws Exception {
        Alignment data = BEASTTestCase.getAlignment();
        Tree tree = BEASTTestCase.getTree(data);
        RealParameter quantiles = new RealParameter("0.0");

        RealParameter stdev = new RealParameter(stdev0 + "");
        PiecewiseLinearDistribution distr = new PiecewiseLinearDistribution();
        LogNormalDistributionModel lognormal = new LogNormalDistributionModel();
        lognormal.initByName("M", "1.0", "meanInRealSpace", true, "S", stdev);
        distr.initByName("distr", lognormal);
    
        UCRelaxedClockModel clockModel = new UCRelaxedClockModel();
        clockModel.initByName("tree", tree, "distr", distr, "clock.rate", "1.0", "rateQuantiles", quantiles, "numberOfDiscreteRates", 100);

        
        double [] x = new double[N];
        double [] y  = new double[N];
        x[0] = 0.1 / N / N;
        for (int i = 1; i < N/3; i++) {
        	x[i] = 0.1 * ((double)i) / (N/3 - 1);
        }
        for (int i = 0; i < N/3; i++) {
        	x[N/3 + i] = 0.1 + 0.7 * ((double)i+1) / (N/3);
        }
        int k = N/3; k = 2 * k;
        for (int i = k; i < N; i++) {
        	x[i] = 0.8 + 0.2 * ((double)i-k) / (N - k);
        }
        x[N-1] = (N - 0.1 / N) / N;
        for (int i = 0; i < N; i++) {
        	y[i] = distr.inverseCumulativeProbability(x[i]);
        }        
        //quadraticApprox.resetData(x, y);
        //cubicApprox.calcDeriv();
        
        UnivariateInterpolator fi = new LinearInterpolator();
        UnivariateFunction f = fi.interpolate(x, y);
        
        
        
        double q = 1.0/100;
        double maxErr = 0, maxRelErr = 0;
        double qAtMaxErr = 0, qAtMaxRelErr = 0;

        double maxErr2 = 0, maxRelErr2 = 0;
        double qAtMaxErr2 = 0, qAtMaxRelErr2 = 0;
        Node node = tree.getNode(0);
        DecimalFormat formatter = new DecimalFormat("#.#####");
        System.out.println("q    r      rApprox  err    relErr   ");
        while (q < 1 - 1.0/100) {
        	quantiles.setValue(0, q);
        	double approxRate = clockModel.getRateForBranch(node);
        	double exactRate = distr.inverseCumulativeProbability(q);
        	double absErr = Math.abs(approxRate - exactRate);
        	double relErr = absErr/exactRate;
        	if (absErr > maxErr) {
        		maxErr = absErr;
        		qAtMaxErr = q;
        	}
        	if (relErr > maxRelErr) {
        		maxRelErr = relErr;
        		qAtMaxRelErr = q;
        	} 

        	// double approxRate2 = quadraticApprox.interpolate(q);
        	double approxRate2 = f.value(q);
        	double absErr2 = Math.abs(approxRate2 - exactRate);
        	double relErr2 = absErr2/exactRate;
        	if (absErr2 > maxErr2) {
        		maxErr2 = absErr2;
        		qAtMaxErr2 = q;
        	}
        	if (relErr2 > maxRelErr2) {
        		maxRelErr2 = relErr2;
        		qAtMaxRelErr2 = q;
        	} 
            System.out.println(formatter.format(q) +  
//            		" " + formatter.format(exactRate) + " " + formatter.format(approxRate) +
//            		" " + formatter.format(absErr) + " " + formatter.format(relErr));
    		" " + formatter.format(exactRate) + " " + formatter.format(approxRate2) +
    		" " + formatter.format(absErr2) + " " + formatter.format(relErr2));
        	q += 0.00109;
        }
        System.out.println("\nmaxErr: " + maxErr + " @q=" + qAtMaxErr);
        System.out.println("maxRelErr: " + maxRelErr + " @q=" + qAtMaxRelErr);
        System.out.println("\nmaxErr2: " + maxErr2 + " @q=" + qAtMaxErr2);
        System.out.println("maxRelErr2: " + maxRelErr2 + " @q=" + qAtMaxRelErr2);
        return maxRelErr;
	}


	@Test
	public void testDerivativeTolerance() throws Exception {
		Randomizer.setSeed(127);
		double stdev = 0.01;
		while (stdev < 1.0) {
			double maxRelErr = testDerivativeTolerance(stdev);
			System.err.println("Stdev:" + stdev + " error: " + maxRelErr);
			// allow up to 5% error
			// assertTrue(maxRelErr < 0.05);
			stdev *= 2;
		}
		System.err.println("testDerivativeTolerance done");	
	}

	private double testDerivativeTolerance(double stdev0) throws Exception {
        Alignment data = BEASTTestCase.getAlignment();
        Tree tree = BEASTTestCase.getTree(data);
        RealParameter quantiles = new RealParameter("0.0");

        RealParameter stdev = new RealParameter(stdev0 + "");
        PiecewiseLinearDistribution distr = new PiecewiseLinearDistribution();
        LogNormalDistributionModel lognormal = new LogNormalDistributionModel();
        lognormal.initByName("M", "1.0", "meanInRealSpace", true, "S", stdev);
        distr.initByName("distr", lognormal);
    
        UCRelaxedClockModel clockModel = new UCRelaxedClockModel();
        clockModel.initByName("tree", tree, "distr", distr, "clock.rate", "1.0", "rateQuantiles", quantiles, "numberOfDiscreteRates", 100);

        RateAndTreeOperator operator1 = new RateAndTreeOperator() {
        	@Override
        	public double proposal() {
        		return 0;
        	}
        };
        operator1.initByName("tree", tree, "clockModel", clockModel, "weight", 1.0);

        
        UCRelaxedClockModel clockModel2 = new UCRelaxedClockModel();
        clockModel2.initByName("tree", tree, "distr", distr, "clock.rate", "1.0", "rateQuantiles", quantiles);

        RateAndTreeOperator operator2 = new RateAndTreeOperator() {
        	@Override
        	public double proposal() {
        		return 0;
        	}
        };
        operator2.initByName("tree", tree, "clockModel", clockModel2, "weight", 1.0);

        double q = 1.0/100;
        double maxErr = 0, maxRelErr = 0;
        double qAtMaxErr = 0, qAtMaxRelErr = 0;
        Node node = tree.getNode(0);
//        DecimalFormat formatter = new DecimalFormat("#.#####");
//        System.out.println("q    r      rApprox  err    relErr   ");
        while (q < 1 - 1.0/100) {
        	quantiles.setValue(0, q);
        	double qNew = Randomizer.nextDouble() * 0.9 + 0.05;
        	double r = distr.inverseCumulativeProbability(qNew);
        	double approxLogHR1 = operator1.calculateHastingsRatio(r, q, qNew);
        	double exactLogHR = operator2.calculateHastingsRatio(r, q, qNew);
        	
        	double absErr = Math.abs(approxLogHR1 - exactLogHR);
        	double relErr = absErr/Math.abs(exactLogHR);
        	if (absErr > maxErr) {
        		maxErr = absErr;
        		qAtMaxErr = q;
        	}
        	if (relErr > maxRelErr) {
        		maxRelErr = relErr;
        		qAtMaxRelErr = q;
        	} 
//            System.out.println(formatter.format(q) +  
//            		" " + formatter.format(exactRate) + " " + formatter.format(approxRate) +
//            		" " + formatter.format(absErr) + " " + formatter.format(relErr));
        	q += 0.0109;
        }
        System.out.println("\nmaxErr: " + maxErr + " @q=" + qAtMaxErr);
        System.out.println("maxRelErr: " + maxRelErr + " @q=" + qAtMaxRelErr);
        return maxRelErr;
	}


	@Test
	public void testDerivativeAt0() throws Exception {
		double stdev = 0.01;
		while (stdev < 1.5) {
			double maxRelErr = testDerivativeAt0(stdev);
			System.err.println("Stdev:" + stdev + " error: " + maxRelErr);
			// allow up to 5% error
		    assertTrue(maxRelErr < 1e-5);
			stdev *= 2;
		}
		System.err.println("testDerivativeAt0 done");	

	}
	
	private double testDerivativeAt0(double stdev0) throws Exception {
        Alignment data = BEASTTestCase.getAlignment();
        Tree tree = BEASTTestCase.getTree(data);
        RealParameter quantiles = new RealParameter("0.0");

        RealParameter stdev = new RealParameter(stdev0 + "");
        PiecewiseLinearDistribution distr = new PiecewiseLinearDistribution();
        LogNormalDistributionModel lognormal = new LogNormalDistributionModel();
        lognormal.initByName("M", "1.0", "meanInRealSpace", true, "S", stdev);
        distr.initByName("distr", lognormal);
    
        UCRelaxedClockModel clockModel = new UCRelaxedClockModel();
        clockModel.initByName("tree", tree, "distr", distr, "clock.rate", "1.0", "rateQuantiles", quantiles, "numberOfDiscreteRates", 100);

        RateAndTreeOperator operator1 = new RateAndTreeOperator() {
        	@Override
        	public double proposal() {
        		return 0;
        	}
        };
        operator1.initByName("tree", tree, "clockModel", clockModel, "weight", 1.0);

        
        UCRelaxedClockModel clockModel2 = new UCRelaxedClockModel();
        clockModel2.initByName("tree", tree, "distr", distr, "clock.rate", "1.0", "rateQuantiles", quantiles);

        RateAndTreeOperator operator2 = new RateAndTreeOperator() {
        	@Override
        	public double proposal() {
        		return 0;
        	}
        };
        operator2.initByName("tree", tree, "clockModel", clockModel2, "weight", 1.0);

        double q = 1.0/100;
        double maxErr = 0, maxRelErr = 0;
        double qAtMaxErr = 0, qAtMaxRelErr = 0;
        Node node = tree.getNode(0);
//        DecimalFormat formatter = new DecimalFormat("#.#####");
//        System.out.println("q    r      rApprox  err    relErr   ");
        while (q < 1 - 1.0/100) {
        	quantiles.setValue(0, q);
        	double qNew = q;
        	double r = clockModel.getRateForBranch(node); //distr.inverseCumulativeProbability(qNew);
        	double approxLogHR1 = operator1.calculateHastingsRatio(r, q, qNew);
        	double exactLogHR = operator2.calculateHastingsRatio(r, q, qNew);
        	
        	double absErr = Math.abs(approxLogHR1);// - exactLogHR);
        	double relErr = absErr/Math.abs(exactLogHR);
        	if (absErr > maxErr) {
        		maxErr = absErr;
        		qAtMaxErr = q;
        	}
        	if (relErr > maxRelErr) {
        		maxRelErr = relErr;
        		qAtMaxRelErr = q;
        	} 
//            System.out.println(formatter.format(q) +  
//            		" " + formatter.format(exactRate) + " " + formatter.format(approxRate) +
//            		" " + formatter.format(absErr) + " " + formatter.format(relErr));
        	q += 0.0109;
        }
        System.out.println("\nmaxErr: " + maxErr + " @q=" + qAtMaxErr);
        System.out.println("maxRelErr: " + maxRelErr + " @q=" + qAtMaxRelErr);
        return maxRelErr;
	}

}
