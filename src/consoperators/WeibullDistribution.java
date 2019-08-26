package consoperators;

import beast.math.distributions.ParametricDistribution;
import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.ContinuousDistribution;
import org.apache.commons.math.distribution.Distribution;
import org.apache.commons.math.distribution.NormalDistributionImpl;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import org.apache.commons.math.distribution.WeibullDistributionImpl;


public class WeibullDistribution extends ParametricDistribution {
    final public Input<RealParameter> kParameterInput = new Input<>("K", "Shape parameter.");
    final public Input<RealParameter> lambdaParameterInput = new Input<>("Lambda", "Scale parameter.");


    org.apache.commons.math.distribution.WeibullDistribution m_dist = new WeibullDistributionImpl(1.0, 1.0);

    @Override
    public void initAndValidate() {
        refresh();
    }

    /**
     * make sure internal state is up to date *
     */
    @SuppressWarnings("deprecation")
    void refresh() {
        double shape;
        double scale;
        if (kParameterInput.get() == null) {
            shape = 1;
        } else {
            shape = kParameterInput.get().getValue();
        }

        if (lambdaParameterInput.get() == null) {
            scale = 1;
        } else {
            scale = lambdaParameterInput.get().getValue();
        }
        m_dist.setScale(scale);
        m_dist.setShape(shape);
    }

    @Override
    public ContinuousDistribution getDistribution() {
        refresh();
        return m_dist;
    }

    @Override
    public double cumulativeProbability(double x) throws MathException {
        return m_dist.cumulativeProbability(x);
    }

    @Override
    public double inverseCumulativeProbability(double p) throws MathException {
        return m_dist.inverseCumulativeProbability(p);
    }

    @Override
    protected double getMeanWithoutOffset() {
        refresh();
        return m_dist.getScale() * m_dist.getShape();
    }


}