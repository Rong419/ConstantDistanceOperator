package consoperators;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.operators.TreeOperator;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import java.text.DecimalFormat;
import java.util.List;

@Description("For internal nodes: propose a new node time")
public class InConstantDistanceOperator extends TreeOperator {
    public final Input<Double> twindowSizeInput =
            new Input<>("twindowSize", "the size of the window when proposing new node time", Input.Validate.REQUIRED);
    public final Input<BranchRateModel.Base> branchRateModelInput = new Input<>("branchRateModel",
            "A model describing the rates on the branches of the beast.tree.");
    final public Input<RealParameter> rateInput = new Input<>("rates", "the rates associated with nodes in the tree for sampling of individual rates among branches.", Input.Validate.REQUIRED);

    private double twindowSize;
    private RealParameter rates;

    protected BranchRateModel.Base branchRateModel;
    JacobianMatrixDeterminant JD = new JacobianMatrixDeterminant();

    @Override
    public void initAndValidate() {
        twindowSize = twindowSizeInput.get();
        branchRateModel = branchRateModelInput.get();
        rates = rateInput.get();
    }

    @Override
    public double proposal() {
        final Tree tree = treeInput.get(this);
        //the chosen node to work on
        Node node;
        //the original node times
        double t_x;
        double t_j;
        double t_k;
        //the original distances
        double d_j;
        double d_k;
        //the original rates
        double r_k;
        double r_j;

        //the proposed node time
        double t_x_;

        //Step 1: randomly select an internal node, denoted by node x
       int nodeCount = tree.getNodeCount();//return the number of nodes in the tree

       //Node [] node = tree.getNodesAsArray();
       // List nodes = tree.getInternalNodes();
       //int nodeCount = nodes.size();
       //int nodeNr = Randomizer.nextInt(nodeCount);
       //node = tree.getNode(nodeNr);
       //if (node.isRoot() || node.isLeaf()){
           //return Double.NEGATIVE_INFINITY;
       //}

        do {
            final int nodeNr = nodeCount / 2 + 1 + Randomizer.nextInt(nodeCount / 2);
            node = tree.getNode(nodeNr);
        } while (node.isRoot() || node.isLeaf());

        //node
        t_x = node.getHeight();//get the time of this node
        double r_node = branchRateModel.getRateForBranch(node);//branch rate for this node
        // son
        Node son = node.getChild(0);//get the left child of this node, i.e. son
        t_j = son.getHeight();//node time of son
        r_j = branchRateModel.getRateForBranch(son);
        d_j = r_j * (t_x - t_j);
        // daughter
        Node daughter = node.getChild(1);//get the right child of this node, i.e. daughter
        t_k = daughter.getHeight();//node time of daughter
        r_k = branchRateModel.getRateForBranch(daughter);
        d_k = r_k * (t_x - t_k);

        double a = Randomizer.uniform(-twindowSize, twindowSize);

        //Step2: to propose a new node time for this node
        t_x_ = t_x + a;

        double upper = node.getParent().getHeight();
        double lower = Math.max(t_j, t_k);

        if (t_x_<= lower || t_x_ >= upper) {
            return Double.NEGATIVE_INFINITY;
        }
        node.setHeight(t_x_);



       //Step3: propose the new rates
       //there are three rates in total
       //r_node, r_i, r_x
       double r_node_ = r_node * (upper - t_x) / (upper - t_x_);
       double r_j_ = d_j / (t_x_ - t_j);
       double r_k_ = d_k / (t_x_ - t_k);
       //if (r_i_ < 0.1 || r_k_  <0.1 || r_node_ <0.1) {
           //System.out.println("t_x_=" + t_x_);
           //System.out.println("r_node_ =" + r_node_ + ",r_i_=" + r_i_ + ",r_k_=" + r_k_);
           //System.out.println("d_node=" + r_node * (upper - t_x) + ",d_i=" + d_i + ",d_k=" + d_k);
           //return Double.NEGATIVE_INFINITY;
       //}
       // set the proposed new rates
       rates.setValue(node.getNr(), r_node_);
       rates.setValue(son.getNr(), r_j_);
       rates.setValue(daughter.getNr(), r_k_);


       //Step4: calculate the Hastings ratio
       double [][] J = new double[4][4];
       J[0][0] = 1.0;
       J[1][0] = r_j / (t_x_ - t_j);
       J[2][0] = r_k / (t_x_ - t_k);
       J[3][0] = (-r_node * upper) / (upper - t_x_);
       J[1][1] = (t_x - t_j) / (t_x_ - t_j);
       J[2][2] = (t_x - t_k) / (t_x_ - t_k);
       J[3][3] = (upper - t_x) / (upper - t_x_);

       double Det = JD.Determinant(J,3);

       return Math.log(Det);
    }
/*
    public static double Determinant(double[][] Matrix, int N)
    {
        int T0;
        int T1;
        int T2;
        double Num;
        int Cha;
        double[][] B;
        if (N > 0) {
            Cha = 0;
            B = new double[N][N];
            Num = 0;
            if (N == 1) {
                return Matrix[0][0] * Matrix[1][1] - Matrix[0][1] * Matrix[1][0];
            }
            for (T0 = 0; T0 <= N; T0++)
            {
                for (T1 = 1; T1 <= N; T1++)
                {
                    for (T2 = 0; T2 <= N - 1; T2++)
                    {
                        if (T2 == T0) {
                            Cha = 1;
                        }
                        B[T1 - 1][T2] = Matrix[T1][T2 + Cha];
                    }
                    Cha = 0;
                }
                Num = Num + Matrix[0][T0] * Determinant(B, N - 1) * Math.pow((-1), T0);
            }
            return Num;
        } else if (N == 0) {
            return Matrix[0][0];
        }
        return 0;
    }
    */

    @Override
    public double getCoercableParameterValue() {
        return twindowSize;
    }


    @Override
    public void setCoercableParameterValue(double value) {
        twindowSize = value;
    }

    /**
     * called after every invocation of this operator to see whether
     * a parameter can be optimised for better acceptance hence faster
     * mixing
     *
     * @param logAlpha difference in posterior between previous state & proposed state + hasting ratio
     */

    @Override
    public void optimize(double logAlpha) {
        // must be overridden by operator implementation to have an effect
        double delta = calcDelta(logAlpha);

        delta += Math.log(twindowSize);
        twindowSize = Math.exp(delta);
    }

    @Override
    public final String getPerformanceSuggestion() {
        double prob = m_nNrAccepted / (m_nNrAccepted + m_nNrRejected + 0.0);
        double targetProb = getTargetAcceptanceProbability();

        double ratio = prob / targetProb;
        if (ratio > 2.0) ratio = 2.0;
        if (ratio < 0.5) ratio = 0.5;

        // new scale factor
        double newWindowSize = twindowSize * ratio;

        DecimalFormat formatter = new DecimalFormat("#.###");
        if (prob < 0.10) {
            return "Try setting window size to about " + formatter.format(newWindowSize);
        } else if (prob > 0.40) {
            return "Try setting window size to about " + formatter.format(newWindowSize);
        } else return "";
    }


}

