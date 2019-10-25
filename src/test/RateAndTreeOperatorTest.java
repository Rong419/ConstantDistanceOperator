package test;

import java.text.DecimalFormat;

import org.junit.Test;

import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.branchratemodel.UCRelaxedClockModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.math.distributions.LogNormalDistributionModel;
import junit.framework.TestCase;
import test.beast.BEASTTestCase;

public class RateAndTreeOperatorTest extends TestCase {

	@Test
	public void testTolerance() throws Exception {
		double stdev = 0.01;
		while (stdev < 1.5) {
			double maxRelErr = testTolerance(stdev);
			// allow up to 5% error
			assertTrue(maxRelErr < 0.05);
			System.err.println("Stdev:" + stdev + " error" + maxRelErr);
			stdev *= 2;
		}
	
	}
	
	private double testTolerance(double stdev0) throws Exception {
        Alignment data = BEASTTestCase.getAlignment();
        Tree tree = BEASTTestCase.getTree(data);
        RealParameter quantiles = new RealParameter("0.0");

        RealParameter stdev = new RealParameter(stdev0 + "");
        LogNormalDistributionModel distr = new LogNormalDistributionModel();
        distr.initByName("M", "1.0", "meanInRealSpace", true, "S", stdev);
    
        UCRelaxedClockModel clockModel = new UCRelaxedClockModel();
        clockModel.initByName("tree", tree, "distr", distr, "clock.rate", "1.0", "rateQuantiles", quantiles, "numberOfDiscreteRates", 100);

        double q = 1.0/100;
        double maxErr = 0, maxRelErr = 0;
        double qAtMaxErr = 0, qAtMaxRelErr = 0;
        Node node = tree.getNode(0);
//        DecimalFormat formatter = new DecimalFormat("#.#####");
//        System.out.println("q    r      rApprox  err    relErr   ");
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
//            System.out.println(formatter.format(q) +  
//            		" " + formatter.format(exactRate) + " " + formatter.format(approxRate) +
//            		" " + formatter.format(absErr) + " " + formatter.format(relErr));
        	q += 0.00109;
        }
        System.out.println("\nmaxErr: " + maxErr + " @q=" + qAtMaxErr);
        System.out.println("maxRelErr: " + maxRelErr + " @q=" + qAtMaxRelErr);
        return maxRelErr;
	}


}
