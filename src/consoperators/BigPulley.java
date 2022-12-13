package consoperators;

import beast.base.core.Citation;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.operator.TreeOperator;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.util.Randomizer;

import java.text.DecimalFormat;
import java.util.List;

@Description("Big pulley: Propose a new genetic distance and a new root time.")
@Citation(value =
        "Zhang, R., Drummond, A. (2020) Improving the performance of Bayesian phylogenetic inference\n" +
                "  under relaxed clock models. BMC Evol Biol 20, 54", DOI = "https://doi.org/10.1186/s12862-020-01609-4",
        year = 2020, firstAuthorSurname = "Zhang")
public class BigPulley extends TreeOperator {
    public final Input<Double> twindowSizeInput =
            new Input<>("twindowSize", "the size of the window for proposing new node time", Input.Validate.REQUIRED);
    public final Input<Double> dwindowSizeInput =
            new Input<>("dwindowSize", "the size of the window for proposing new genetic distance");
    final public Input<RealParameter> rateInput = new Input<>("rates", "the rates associated with nodes in the tree for sampling of individual rates among branches.", Input.Validate.REQUIRED);

    private double twindowSize;
    private double dwindowSize;
    private Tree tree;
    private int branchCount;
    private RealParameter rates;
    private double hastingsRatio;
    private JacobianMatrixDeterminant JD = new JacobianMatrixDeterminant();

    @Override
    public void initAndValidate() {
        twindowSize = twindowSizeInput.get();
        dwindowSize = dwindowSizeInput.get();
        rates = rateInput.get();
        tree = treeInput.get();
    }

    @Override
    public double proposal() {
        final Tree tree = treeInput.get();
        branchCount = tree.getNodeCount() - 1; //the number of branches of the tree

        //the chosen node to work on
        Node node;
        Node nodeOlder = new Node(); //the older child node
        Node nodeYounger = new Node(); //the younger child node
        Node son1 ;
        Node daughter1;
        Node son2;
        Node daughter2;
        //the original node time
        double t_x;
        double t_j;
        double t_k;
        double tOlder = 0.0;
        double tYounger = 0.0;
        double newtOlder = 1.0;
        //the proposed node time
        double t_x_;
        double t_j_;
        double t_k_;

        double r_x;
        double r_i;
        double r_j;
        double r_k;

        double d_i;
        double d_x;
        double d_j;
        double d_k;
        double dOld = 1.0;
        double dOld_ = 1.0;
        double dYoung = 1.0;
        List C1;List C2;List E;
        double ParaA = 1.0; double ParaB =1.0; double ParaC;
        double det = 1.0;

        // Step 1: randomly select a node, denoted by node x
        node = tree.getRoot();
        t_x = node.getHeight();//get the time of this node

        // information of son
        Node son = node.getChild(0);//get the left child of this node, i.e. son
        t_j = son.getHeight();//node time of son
        r_i = rates.getValues()[getNodeNr(son)];
        d_i = r_i * (t_x - t_j);

        // information of daughter
        Node daughter = node.getChild(1);//get the right child of this node, i.e. daughter
        t_k = daughter.getHeight();//node time of daughter
        r_x = rates.getValues()[getNodeNr(daughter)];
        d_x = r_x * (t_x - t_k);

        // random numbers used in the proposal
        double a = Randomizer.uniform(-twindowSize, twindowSize);//for the root time
        double a1 = Randomizer.uniform(-twindowSize, twindowSize);//for the son's node time
        double a2 = Randomizer.uniform(-twindowSize, twindowSize);//for the daughter's node time
        double u1 = Randomizer.uniform(0, 1);//which node to pick, i.e. the son or the daughter of the root
        double u2 = Randomizer.uniform(0, 1);//which tree topology to go, i.e. which node to exchange
        double d0 = Randomizer.uniform(-dwindowSize, dwindowSize);//for distance of branch above son
        double d1 = Randomizer.uniform(-dwindowSize, dwindowSize);//for distance of branch above daughter

        // Step 2: to propose a new values
        t_x_ = t_x + a;  // new root time
        t_j_ = t_j + a1; // new node time of  son
        t_k_ = t_k + a2; // new node time of daughter
        double d_i_ = d_i + d0; // new distance of branch above son
        double d_x_ = d_x + d1; // new distance of branch above daughter


        // to prepare the constant numbers for jacobian
        // which is dependent on how the tree topology is proposed
        if (t_j > t_k) {
            ParaA = 0.5;
            ParaB = 1.0;
        }
        if (t_j < t_k) {
            ParaA = 1.0;
            ParaB = 0.5;
        }
        if (t_j == t_k) {
            ParaA = 1.0;
            ParaB = 1.0;
        }

        // get grandchildren of the root
        // to decide which tree shape it is
        C1 = son.getChildren();
        C2 = daughter.getChildren();

        // symmetric tree shapes
        // i.e. both son and daughter have child nodes
        if ((C1.size() != 0) && (C2.size() != 0)) {
            // access to child nodes of son and daughter
            son1 = son.getChild(0);
            daughter1 = son.getChild(1);
            son2 = daughter.getChild(0);
            daughter2 = daughter.getChild(1);

            double t3 = son1.getHeight();
            double t4 = daughter1.getHeight();
            double t5 = son2.getHeight();
            double t6 = daughter2.getHeight();

            r_j = rates.getValues()[getNodeNr(son1)];
            r_k = rates.getValues()[getNodeNr(daughter1)];
            double r_m = rates.getValues()[getNodeNr(son2)];
            double r_n = rates.getValues()[getNodeNr(daughter2)];

            d_j = r_j * (t_j - t3);
            d_k = r_k * (t_j - t4);
            double d_m = r_m * (t_k - t5);
            double d_n = r_n * (t_k - t6);

            // with 0.5 probability to pick the son and change the daughter
            if (u1 > 0.5) {
                if (t_j_ <= t_k) {
                    return Double.NEGATIVE_INFINITY;
                }
                if (t_x_ <= t_j_) {
                    return Double.NEGATIVE_INFINITY;
                }
                node.setHeight(t_x + a);
                son.setHeight(t_j_);

                //Case1: daughter <-> son1
                if (u2 > 0.5) {
                    if (d_j - d_i_ <= 0) {
                        return Double.NEGATIVE_INFINITY;
                    }
                    TreeChange(son1, daughter1, daughter, son, node, d_j, d_k, d_x, d_i, d_i_, t3, t4, t_k, t_j_, t_x_);
                    det = Det(r_j,r_k,r_x,r_i,t_x,t_x_,t_j,t_j_,t3,t4,t_k);
                    //det = Det2 (t_x_,t_j_,t3,t4,t_k);
                    E = son1.getChildren();
                    if (E.size() != 0) {
                        ParaA = 0.25;
                    }
                }

                //Case2: daughter <-> daughter1
                else if (u2 < 0.5) {

                    if (d_k - d_i_ <= 0) {
                        return Double.NEGATIVE_INFINITY;
                    }
                    TreeChange(daughter1, son1, daughter, son, node, d_k, d_j, d_x, d_i, d_i_, t4, t3, t_k, t_j_, t_x_);
                    det = Det(r_k,r_j,r_x,r_i,t_x,t_x_,t_j,t_j_,t4,t3,t_k);
                    //det = Det2 (t_x_,t_j_,t4,t3,t_k);
                    E = daughter1.getChildren();
                    if (E.size() != 0) {
                        ParaA = 0.25;
                    }
                } else {
                    return Double.NEGATIVE_INFINITY;
                }
                //hastingsRatio = Math.log(ParaA / 0.25);
                hastingsRatio = Math.log(4 * det * ParaA);
            }

            // with 0.5  probability to pick the daughter and change the son
            else if (u1 < 0.5) {
                if (t_k_ <= t_j) {
                    return Double.NEGATIVE_INFINITY;
                }
                if (t_x_ <= t_k_) {
                    return Double.NEGATIVE_INFINITY;
                }
                node.setHeight(t_x + a);
                daughter.setHeight(t_k_);

                //Case3: son <-> son2
                if (u2 > 0.5) {
                    if (d_m - d_x_ <= 0) {
                        return Double.NEGATIVE_INFINITY;
                    }
                    TreeChange(son2, daughter2, son, daughter, node, d_m, d_n, d_i, d_x, d_x_, t5, t6, t_j, t_k_, t_x_);
                    det = Det(r_m,r_n,r_i,r_x,t_x,t_x_,t_k,t_k_,t5,t6,t_j);
                    //det = Det2 (t_x_,t_k_,t5,t6,t_j);
                    E = son2.getChildren();
                    if (E.size() != 0) {
                        ParaB = 0.25;
                    }
                }

                //Case4: son <-> daughter2
                else if (u2 < 0.5) {
                    if (d_n - d_x_ <= 0) {
                        return Double.NEGATIVE_INFINITY;
                    }
                    TreeChange(daughter2, son2, son, daughter, node, d_n, d_m, d_i, d_x, d_x_, t6, t5, t_j, t_k_, t_x_);
                    det = Det(r_n,r_m,r_i,r_x,t_x,t_x_,t_k,t_k_,t6,t5,t_j);
                    //det = Det2 (t_x_,t_k_,t6,t5,t_j);
                    E = daughter2.getChildren();
                    if (E.size() != 0) {
                        ParaB = 0.25;
                    }
                } else {
                    return Double.NEGATIVE_INFINITY;
                }
                //hastingsRatio = Math.log(ParaB / 0.25);
                hastingsRatio = Math.log(4 * det * ParaB);
            } else {
                return Double.NEGATIVE_INFINITY;
            }
        }

        // asymmetric tree shape
        // i.e. either son or daughter is a tip node
        if (((C1.size() == 0) && (C2.size() != 0)) || ((C1.size() != 0) && (C2.size() == 0))) {
            // to identify which node is a tip (younger one)
            if (C2.size() == 0) {
                tOlder = t_j;
                tYounger = t_k;
                newtOlder = t_j_;
                dOld = d_i;
                dOld_ = d_i_;
                dYoung = d_x;
                nodeOlder = son;
                nodeYounger = daughter;
            }
            if (C1.size() == 0) {
                tOlder = t_k;
                tYounger = t_j;
                newtOlder = t_k_;
                nodeOlder = daughter;
                nodeYounger = son;
                dOld = d_x;
                dOld_ = d_x_;
                dYoung = d_i;
            }

            if ((t_x_ <= newtOlder) || (newtOlder <= tYounger)) {
                return Double.NEGATIVE_INFINITY;
            }
            node.setHeight(t_x + a);
            nodeOlder.setHeight(newtOlder);

            // access to child nodes of the node that is not a tip
            Node Child1 = nodeOlder.getChild(0);
            Node Child2 = nodeOlder.getChild(1);
            double t_Ch1 = Child1.getHeight();
            double t_Ch2 = Child2.getHeight();

            if ((newtOlder <= t_Ch1) || (newtOlder <= t_Ch2)) {
                return Double.NEGATIVE_INFINITY;
            }

            double r_Ch1 = rates.getValues()[getNodeNr(Child1)];
            double r_Ch2 = rates.getValues()[getNodeNr(Child2)];
            double d_Ch1 = r_Ch1 * (tOlder - t_Ch1);
            double d_Ch2 = r_Ch2 * (tOlder - t_Ch2);

            double r_Y = rates.getValues()[getNodeNr(nodeYounger)];
            double r_O = rates.getValues()[getNodeNr(nodeOlder)];

            if ((newtOlder > Math.max(t_Ch1, t_Ch2)) || (t_Ch1 == t_Ch2)) {

                //Case5: Child1 <-> nodeYounger
                if (u2 > 0.5) {
                    if (d_Ch1 - dOld_ <= 0) {
                        return Double.NEGATIVE_INFINITY;
                    }
                    TreeChange(Child1, Child2, nodeYounger, nodeOlder, node, d_Ch1, d_Ch2, dYoung, dOld, dOld_, t_Ch1, t_Ch2, tYounger, newtOlder, t_x_);
                    det = Det(r_Ch1,r_Ch2,r_Y,r_O,t_x,t_x_,tOlder,newtOlder,t_Ch1,t_Ch2,tYounger);
                    //det = Det2 (t_x_,newtOlder,t_Ch1,t_Ch2,tYounger);
                    List C = Child1.getChildren();
                    if (C.size() == 0) {
                        ParaC = 0.5;
                    } else {
                        ParaC = 0.25;
                    }
                }

                //Case6: Child2 <-> nodeYounger
                else if (u2 < 0.5) {
                    if (d_Ch2 - dOld_ <= 0) {
                        return Double.NEGATIVE_INFINITY;
                    }
                    TreeChange(Child2, Child1, nodeYounger, nodeOlder, node, d_Ch2, d_Ch1, dYoung, dOld, dOld_, t_Ch2, t_Ch1, tYounger, newtOlder, t_x_);
                    det = Det(r_Ch2,r_Ch1,r_Y,r_O,t_x,t_x_,tOlder,newtOlder,t_Ch2,t_Ch1,tYounger);
                    //det = Det2 (t_x_,newtOlder,t_Ch2,t_Ch1,tYounger);
                    List D = Child2.getChildren();
                    if (D.size() == 0) {
                        ParaC = 0.5;
                    } else {
                        ParaC = 0.25;
                    }
                } else {
                    return Double.NEGATIVE_INFINITY;
                }
                //hastingsRatio = Math.log(ParaC / 0.5);
                hastingsRatio = Math.log(2 * det * ParaC);
            }

            else if ((newtOlder > Math.min(t_Ch1, t_Ch2)) && (newtOlder <= Math.max(t_Ch1, t_Ch2))) {

                // /Case7: Child1 <-> nodeYounger
                if (t_Ch1 > t_Ch2) {
                    if (d_Ch1 - dOld_ <= 0) {
                        return Double.NEGATIVE_INFINITY;
                    }
                    TreeChange(Child1, Child2, nodeYounger, nodeOlder, node, d_Ch1, d_Ch2, dYoung, dOld, dOld_, t_Ch1, t_Ch2, tYounger, newtOlder, t_x_);
                    det = Det(r_Ch1,r_Ch2,r_Y,r_O,t_x,t_x_,tOlder,newtOlder,t_Ch1,t_Ch2,tYounger);
                    //det = Det2 (t_x_,newtOlder,t_Ch1,t_Ch2,tYounger);
                }
                if (t_Ch1 < t_Ch2) {
                    //Case8: Child2 <-> nodeYonger
                    if (d_Ch2 - dOld_ <= 0) {
                        return Double.NEGATIVE_INFINITY;
                    }
                    TreeChange(Child2, Child1, nodeYounger, nodeOlder, node, d_Ch2, d_Ch1, dYoung, dOld, dOld_, t_Ch2, t_Ch1, tYounger, newtOlder, t_x_);
                    det = Det(r_Ch2,r_Ch1,r_Y,r_O,t_x,t_x_,tOlder,newtOlder,t_Ch2,t_Ch1,tYounger);
                    //det = Det2 (t_x_,newtOlder,t_Ch2,t_Ch1,tYounger);
                }
                //hastingsRatio = Math.log(0.25);
                hastingsRatio = Math.log(0.25 * det);
            } else {
                return Double.NEGATIVE_INFINITY;
            }
        }

        return hastingsRatio;
    }

    private void TreeChange(Node ChangeA, Node SiblingA, Node ChangeC, Node ParentA, Node Root,
                              double dChangeA, double dSiblingA, double dChangC, double dParentA, double newdParentA,
                              double ChangeATime, double SiblingATime, double ChangeCTime,double newParentATime, double newRootTime) {
        //propose new rates
        double r_3_ = (dChangC + dParentA) / (newParentATime - ChangeCTime);
        double r_2_ = dSiblingA / (newParentATime - SiblingATime);
        double r_1_ = (dChangeA - newdParentA) / (newRootTime - ChangeATime);
        double r_4_ = newdParentA / (newRootTime - newParentATime);

        //exchange node A and C
        exchangeNodes(ChangeA, ChangeC, ParentA, Root);
        tree.setRoot(Root);

        //assign new rates to the nodes
        rates.setValue(getNodeNr(ChangeC), r_3_);
        rates.setValue(getNodeNr(SiblingA), r_2_);
        rates.setValue(getNodeNr(ChangeA), r_1_);
        rates.setValue(getNodeNr(ParentA), r_4_);


    }

    //calculate the determinant of the Jacobian matrix
    private double Det (double r1, double r2, double r3, double r4, double T, double T_, double t,double t_, double t1, double t2, double t3) {
        double [][] J = new double[6][6];
        //f(T,t,r1,r2,r3,r4) ---> T_,t_,r1_,r2_,r3_,r4_
        J[0][0] = 1;
        J[1][1] = 1;
        J[2][0] = -r4 / (T_ - t1);
        J[2][1] = (r1 + r4) / (T_ - t1);
        J[2][2] = (t- t1) / (T_ - t1);
        J[2][5] = (t - T) / (T_ - t1);
        J[3][1] = r2 / (t_ - t2);
        J[3][3] = (t - t2) / (t_ - t2);
        J[4][0] = (r3 + r4) / (t_ - t3);
        J[4][1] = -r4 / (t_ - t3);
        J[4][4] = (T - t3) / (t_ - t3);
        J[4][5] = (T - t) / (t_ - t3);
        J[5][0] = r4 / (T_ - t_);
        J[5][1] = -r4 / (T_ - t_);
        J[5][5] = (T - t) / (T_ - t_);

        return JD.Determinant(J,5);
    }

    private void exchangeNodes(Node i, Node j, Node p, Node jP) {
        //precondition p -> i & jP -> j
        replace(p, i, j);

        //postcondition p -> j & p -> i
        replace(jP, j, i);

    }

   private int getNodeNr (Node node) {
       int nodenumber = node.getNr();
       if (nodenumber == branchCount) {
           nodenumber =node.getTree().getRoot().getNr();
       }
       return nodenumber;
   }



    /*
     * automatic parameter tuning *
     */

    @Override
    public double getCoercableParameterValue() {
        return dwindowSize;
    }

    @Override
    public void setCoercableParameterValue(double value) {
        dwindowSize = value;
    }

    /*
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

        delta += Math.log(dwindowSize);
        dwindowSize = Math.exp(delta);
    }

    @Override
    public final String getPerformanceSuggestion() {
        double prob = m_nNrAccepted / (m_nNrAccepted + m_nNrRejected + 0.0);
        double targetProb = getTargetAcceptanceProbability();

        double ratio = prob / targetProb;
        if (ratio > 2.0) ratio = 2.0;
        if (ratio < 0.5) ratio = 0.5;

        // new scale factor
        double newWindowSize = dwindowSize * ratio;

        DecimalFormat formatter = new DecimalFormat("#.###");
        if (prob < 0.10) {
            return "Try setting window size to about " + formatter.format(newWindowSize);
        } else if (prob > 0.40) {
            return "Try setting window size to about " + formatter.format(newWindowSize);
        } else return "";
    }

}

