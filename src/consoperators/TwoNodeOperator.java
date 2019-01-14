package consoperators;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.operators.TreeOperator;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import beast.evolution.tree.Node;

@Description("For internal nodes: propose two node times and rates that are associated with the node times")
public class TwoNodeOperator extends TreeOperator {
    final public Input<RealParameter> rateInput = new Input<>("rates", "the rates associated with nodes in the tree for sampling of individual rates among branches.", Input.Validate.REQUIRED);
    public final Input<BranchRateModel.Base> branchRateModelInput = new Input<>("branchRateModel",
            "A model describing the rates on the branches of the beast.tree.");

    RealParameter rates;
    BranchRateModel.Base branchRateModel;


    @Override
    public void initAndValidate() {
        rates = rateInput.get();
        branchRateModel = branchRateModelInput.get();
    }

    @Override
    public double proposal() {
        //the tree to work with
        Tree tree = treeInput.get(this);
        //the internal node to operate on
        Node node;

        /*
        Step1: randomly select one node in the tree
         */
        int nodeCount = tree.getNodeCount();
        int nodeNr = nodeCount / 2 + 1 + Randomizer.nextInt(nodeCount / 2);
        node = tree.getNode(nodeNr);

        if (node.isLeaf()){

        }

        else if (node.isRoot()) {

        }
        else {

        }

        return 0;
    }

    public void FourRates () {

    }

    public void FiveRates () {

    }


}
