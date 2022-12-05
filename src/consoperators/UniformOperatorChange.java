package consoperators;

import beast.base.core.Input;
import beast.base.inference.Operator;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.Parameter;
import beast.base.util.Randomizer;


// A uniform operator which acts on integers and always samples a new value
public class UniformOperatorChange extends Operator {

	final public Input<IntegerParameter> parameterInput = new Input<>("parameter", "an integer parameter to sample individual values for", Validate.REQUIRED, Parameter.class);

    IntegerParameter parameter;
    double lower, upper;

    @Override
    public void initAndValidate() {
        parameter = parameterInput.get();
    }

    @Override
    public double proposal() {
    	parameter = parameterInput.get();
        int lowerIndex = parameter.getLower();
        int upperIndex = parameter.getUpper();
        int range = upperIndex - lowerIndex + 1;
    	
    	// (index.getValue() + 1 + Randomizer.nextInt(clockModelDistributions.size() - 1)) % clockModelDistributions.size();
		if (range <= 1) return Double.NEGATIVE_INFINITY;
    	
        int index = Randomizer.nextInt(parameter.getDimension());
        int initalValue = parameter.getValue(index);
        
        
       
        
        int delta = initalValue + Randomizer.nextInt(range-1) + 1 - lowerIndex; // Increase in index
        delta = delta % range; // If the index goes beyond the range, then wrap it back under
        int newValue = delta + lowerIndex;
        parameter.setValue(index, newValue);
        
        //System.out.println("Changing from " + initalValue + " to " + newValue + " lower = " + lowerIndex + " upper = " + upperIndex);
        return 0.0;
        
    }

}
