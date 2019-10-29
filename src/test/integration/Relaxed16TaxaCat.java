package test.integration;

// Requirements: the file analysis.xml must be in the directory
// from where this script is run.

import beast.evolution.substitutionmodel.*;
import beast.evolution.tree.Tree;
import beast.math.distributions.LogNormalDistributionModel;
import beast.evolution.sitemodel.*;
import beast.evolution.alignment.*;
import beast.evolution.branchratemodel.UCRelaxedClockModel;
import beast.util.*;
import beast.app.seqgen.MergeDataWith;
import beast.app.seqgen.SequenceSimulator;
import beast.core.*;
import beast.core.parameter.*;

import java.io.File;
import java.io.IOException;
import java.util.List;

import beagle.BeagleFlag;

public class Relaxed16TaxaCat {
	// This script performs a simulation study for the
	// bModelTest model by sampling from the prior, simulate sequences and
	// running ananlyses to recover the model.


	static int N = 100;


	static void process(String dir) throws IllegalArgumentException, IllegalAccessException, IOException, XMLParserException {

		if (!(new File(dir).exists())) {
			new File(dir).mkdirs();
		}
		System.err.print("Processing " + dir);

		
		LogAnalyser trace = new LogAnalyser(wdir + "/cat16.log", 0, false, false);
		Double [] kappa = trace.getTrace("kappa");
		Double [] freqParameter1 = trace.getTrace("freqParameter.1");
		Double [] freqParameter2 = trace.getTrace("freqParameter.2");
		Double [] freqParameter3 = trace.getTrace("freqParameter.3");
		Double [] freqParameter4 = trace.getTrace("freqParameter.4");
		Double [] s = trace.getTrace("ucldStdev");
		
		Double [][] rateCategoriesX = new Double[30][];
		for (int i = 0; i < 30; i++) {
			rateCategoriesX[i] = trace.getTrace("rateCategories." + (i+1));
		}
		
	for (int i = 1; i < N+1; i++) {

		// set up model to draw samples from
		Sequence A=new Sequence(); A.initByName("taxon","sequence01","value","?");
		Sequence B=new Sequence(); B.initByName("taxon","sequence02","value","?");
		Sequence C=new Sequence(); C.initByName("taxon","sequence03","value","?");
		Sequence D=new Sequence(); D.initByName("taxon","sequence04","value","?");
		Sequence E=new Sequence(); E.initByName("taxon","sequence05","value","?");
		Sequence F=new Sequence(); F.initByName("taxon","sequence06","value","?");
		Sequence G=new Sequence(); G.initByName("taxon","sequence07","value","?");
		Sequence I=new Sequence(); I.initByName("taxon","sequence08","value","?");
		Sequence H=new Sequence(); H.initByName("taxon","sequence09","value","?");
		Sequence J=new Sequence(); J.initByName("taxon","sequence10","value","?");
		Sequence K=new Sequence(); K.initByName("taxon","sequence11","value","?");
		Sequence L=new Sequence(); L.initByName("taxon","sequence12","value","?");
		Sequence M=new Sequence(); M.initByName("taxon","sequence13","value","?");
		Sequence N=new Sequence(); N.initByName("taxon","sequence14","value","?");
		Sequence O=new Sequence(); O.initByName("taxon","sequence15","value","?");
		Sequence P=new Sequence(); P.initByName("taxon","sequence16","value","?");
		Alignment data = new Alignment();
		data.initByName("sequence", A, "sequence", B, "sequence", C, "sequence", D, "sequence", E,
				"sequence", F, "sequence", G, "sequence", H, "sequence", I, "sequence", J,
				"sequence", K, "sequence", L, "sequence", M, "sequence", N, "sequence", O,
				"sequence", P
				);
//		tree = new beast.util.TreeParser(newick="(((A:0.1,B:0.1),C:0.15):0.05,(D:0.1,E:0.1):0.1)", taxa=data, IsLabelledNewick=true);
		Tree tree = trees.get(i);


		
		RealParameter freqs=new RealParameter(freqParameter1[i] + " " + freqParameter2[i] + " " + freqParameter3[i] + " " + freqParameter4[i]);
		Frequencies f = new Frequencies();
		f.initByName("frequencies",freqs);
		
		HKY hky = new beast.evolution.substitutionmodel.HKY();
		hky.initByName("frequencies", f, 
			"kappa", kappa[i]+""
		);
		LogNormalDistributionModel distr = new beast.math.distributions.LogNormalDistributionModel();
		distr.initByName("S", s[i] + "", "meanInRealSpace", true, "M", "1.0");
		
		
		String rates = "";
		for (int j = 0; j < 30; j++) {
			rates += (int)((double)rateCategoriesX[j][i]) + " ";
		}
		IntegerParameter  rateCategories = new beast.core.parameter.IntegerParameter(rates);
		
		UCRelaxedClockModel clockmodel = new beast.evolution.branchratemodel.UCRelaxedClockModel();
		clockmodel.initByName("distr", distr, "rateCategories", rateCategories, "tree", tree);



	    SiteModel sitemodel = new SiteModel();
		sitemodel.initByName("gammaCategoryCount", 1, "substModel", hky, "shape", "1.0", "proportionInvariant", "0.0");
		
		MergeDataWith mergewith = new beast.app.seqgen.MergeDataWith();
		mergewith.initByName("template", wdir + "/analysis16relaxedCat.xml", "output", dir + "/analysis-out" + i + ".xml");

		SequenceSimulator sim = new beast.app.seqgen.SequenceSimulator();
		sim.initByName("data", data, "tree", tree, "sequencelength", 1000, "outputFileName", 
				"gammaShapeSequence.xml", "siteModel", sitemodel, "branchRateModel", clockmodel, 
				"merge", mergewith);
		// produce gammaShapeSequence.xml and merge with analysis.xml to get analysis-out.xml
		sim.run();
		System.err.print('.');
	  }
	System.err.println();
	}


	static List<Tree> trees;
	static String wdir = "/Users/remco/workspace/ConstantDistanceOperator/validation/quantiles";

			public static void main(String[] args) throws IOException, IllegalArgumentException, IllegalAccessException, XMLParserException {
		
		Logger.FILE_MODE = beast.core.Logger.LogFileMode.overwrite;

		// set up flags for BEAGLE -- YMMV
		long beagleFlags = BeagleFlag.VECTOR_SSE.getMask() | BeagleFlag.PROCESSOR_CPU.getMask();
		System.setProperty("beagle.preferred.flags", Long.toString(beagleFlags));

		
		NexusParser parser = new beast.util.NexusParser();
		File fin = new File(wdir + "/cat16.trees");
		parser.parseFile(fin);
		trees = parser.trees;

		process(wdir + "/cat16/");
		
		System.err.println("Done");

	}

}
