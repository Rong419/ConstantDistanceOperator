package test;

import consoperators.NarrowExchangeOperator;
import org.junit.Test;

import beast.base.inference.State;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.util.Randomizer;
import beast.base.evolution.tree.TreeParser;
import junit.framework.TestCase;


public class NarrowExchangeOperatorTest extends TestCase {

    @Test
    public void test4Taxa() throws Exception {

        int runs = 1;
        Randomizer.setSeed(666);
        // test that going from source tree to target tree
        // is as likely as going the other way around
        // taking the HR in account.
        Sequence A = new Sequence("A", "A");
        Sequence B = new Sequence("B", "A");
        Sequence C = new Sequence("C", "A");
        Sequence D = new Sequence("D", "A");

        Alignment data = new Alignment();
        data.initByName("sequence", A, "sequence", B, "sequence", C, "sequence", D,
                "dataType", "nucleotide"
        );
        String sourceTree = "((A:2.0,B:2.0):1.0,(C:1.0,D:1.0):2.0):0.0"; // ((A,B),(C,D))
        testNarrowExchange(sourceTree,runs, data);
    }

    void testNarrowExchange(String sourceTree, int runs, Alignment data) throws Exception {

        // first test going from source to target
        double match = 0;
        for (int i = 0; i < runs; i++) {
            TreeParser tree = new TreeParser();
            tree.initByName("taxa", data, "newick", sourceTree, "IsLabelledNewick", true);
            State state = new State();
            state.initByName("stateNode", tree);
            state.initialise();
            NarrowExchangeOperator operator = new NarrowExchangeOperator();
            operator.initByName("isNarrow", true, "tree", tree, "weight", 1.0);

            double logHR = operator.proposal();
            String treeString = tree.getRoot().toNewick();
        }
    }

}