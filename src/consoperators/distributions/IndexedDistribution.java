package consoperators.distributions;

import java.util.*;

import beast.base.inference.Distribution;
import beast.base.inference.distribution.ParametricDistribution;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.IntegerParameter;

@Description("Selects a distribution from among its inputs based on index variable")
public class IndexedDistribution extends ParametricDistribution {
	final public Input<IntegerParameter> indexInput = new Input<>("index", "index of the distribution to be selected", Validate.REQUIRED);
	final public Input<List<ParametricDistribution>> distrInput = new Input<>("distr", "one or more distributions to be chosen from", new ArrayList<>(), Validate.REQUIRED);

	IntegerParameter index;
	List<ParametricDistribution> distrs;
	
	@Override
	public void initAndValidate() {
		distrs = distrInput.get();
		
		index = indexInput.get();
		index.setLower(0);
		index.setUpper(distrs.size()-1);
		
		if (offsetInput.get() != 0.0) {
			throw new IllegalArgumentException("Offset should be set on inputs of this IndexedDistribution, not on the IndexedDistribution itself");
		}
	}

	@Override
	public org.apache.commons.math.distribution.Distribution getDistribution() {
		return distrs.get(index.getValue()).getDistribution();
	}

	@Override
	public double getOffset() {
		return distrs.get(index.getValue()).getOffset();
	}
}
