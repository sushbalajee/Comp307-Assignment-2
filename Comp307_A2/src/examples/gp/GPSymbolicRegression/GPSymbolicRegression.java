package examples.gp.GPSymbolicRegression;

import org.jgap.*;
import org.jgap.gp.*;
import org.jgap.gp.terminal.*;
import org.jgap.gp.impl.*;
import org.jgap.gp.function.*;

public class GPSymbolicRegression extends GPProblem{

	  public static Variable vx;

	  protected static Float[] x = {-2.0f, -1.75f, -1.50f, -1.25f, -1.00f, -0.75f, -0.50f, -0.25f, 0.00f, 0.25f, 0.50f, 0.75f, 1.00f, 1.25f, 1.50f, 1.75f, 2.00f, 2.25f, 2.50f, 2.75f};    
	  protected static Float[] y = {37.00000f, 24.16016f, 15.06250f, 8.91016f, 5.00000f, 2.72266f, 1.56250f, 1.09766f, 1.00000f, 1.03516f, 1.06250f, 1.03516f, 1.00000f, 1.09766f, 1.56250f, 2.72266f, 5.00000f, 8.91016f, 15.06250f, 24.16016f};

	  private static int genCount = 800;
	  
	  public GPSymbolicRegression(GPConfiguration a_conf)
	      throws InvalidConfigurationException {
	    super(a_conf);
	  }

	  /**
	   * This method is used for setting up the commands and terminals that can be
	   * used to solve the problem.
	   * Please notice, that the variables types, argTypes and nodeSets correspond
	   * to each other: they have the same number of elements and the element at
	   * the i'th index of each variable corresponds to the i'th index of the other
	   * variables!
	   *
	   * @return GPGenotype
	   * @throws InvalidConfigurationException
	   */
	  public GPGenotype create()
	      throws InvalidConfigurationException {
	    GPConfiguration conf = getGPConfiguration();
	    // At first, we define the return type of the GP program.
	    // ------------------------------------------------------
	    @SuppressWarnings("rawtypes")
		Class[] types = {
	        // Return type of result-producing chromosome
	        CommandGene.FloatClass};
	    // Then, we define the arguments of the GP parts. Normally, only for ADF's
	    // there is a specification here, otherwise it is empty as in first case.
	    // -----------------------------------------------------------------------
	    @SuppressWarnings("rawtypes")
		Class[][] argTypes = {
	        // Arguments of result-producing chromosome: none
	        {},
	    };
	    // Next, we define the set of available GP commands and terminals to use.
	    // ----------------------------------------------------------------------
	    CommandGene[][] nodeSets = {{
	        // We use a variable that can be set in the fitness function.
	        // ----------------------------------------------------------
	        vx = Variable.create(conf, "X", CommandGene.FloatClass),
	        new Multiply(conf, CommandGene.FloatClass),
	        new Add(conf, CommandGene.FloatClass),
	        new Divide(conf, CommandGene.FloatClass),
	        new Subtract(conf, CommandGene.FloatClass),
	        new Terminal(conf, CommandGene.FloatClass, 1.0d, 10.0d, true),
	    }};
	    // ------------------------------------------------------------------------
	    return GPGenotype.randomInitialGenotype(conf, types, argTypes, nodeSets,
	        20, true);
	  }
	
	  /**
	   * Fitness function for evaluating the produced fomulas, represented as GP
	   * programs. The fitness is computed by calculating the result (Y) of the
	   * function/formula for integer inputs 0 to 20 (X). The sum of the differences
	   * between expected Y and actual Y is the fitness, the lower the better (as
	   * it is a defect rate here).
	   */
	  @SuppressWarnings("serial")
	public static class FormulaFitnessFunction
	      extends GPFitnessFunction {
	    protected double evaluate(final IGPProgram a_subject) {
	      return computeRawFitness(a_subject);
	    }

	    public double computeRawFitness(final IGPProgram ind) {
	      double error = 0.0f;
	      Object[] noargs = new Object[0];
	      // Evaluate function for input numbers 0 to 20.
	      // --------------------------------------------
	      for (int i = 0; i < 20; i++) {
	        // Provide the variable X with the input number.
	        // See method create(), declaration of "nodeSets" for where X is
	        // defined.
	        // -------------------------------------------------------------
	        vx.set(x[i]);
	        try {
	          // Execute the GP program representing the function to be evolved.
	          // As in method create().
	          // ----------------------------------------------------------------
	          double result = ind.execute_float(0, noargs);
	          // Sum up the error between actual and expected result to get a defect
	          // rate.
	          // -------------------------------------------------------------------
	          error += Math.abs(result - y[i]);
	          // If the error is too high, stop evlauation and return worst error
	          // possible.
	          // ----------------------------------------------------------------
	          if (Double.isInfinite(error)) {
	            return Double.MAX_VALUE;
	          }
	        } catch (ArithmeticException ex) {
	          // This should not happen, some illegal operation was executed.
	          // ------------------------------------------------------------
	          System.out.println("x = " + x[i].floatValue());
	          System.out.println(ind);
	          throw ex;
	        }
	      }
	      // In case the error is small enough, consider it perfect.
	      // -------------------------------------------------------
	      if (error < 0.001) {
	        error = 0.0d;
	      }
	      return error;
	    }
	  }
	  

	  public static void main(String[] args)
	      throws Exception {
	    System.out.println("Formula to discover: X^4 + 2X^3 + X^2 - X");
	    GPConfiguration config = new GPConfiguration();
	    config.setGPFitnessEvaluator(new DeltaGPFitnessEvaluator());
	    config.setFitnessFunction(new FormulaFitnessFunction());
	    config.setMaxInitDepth(4);
	    config.setPopulationSize(1000);
	    config.setMaxCrossoverDepth(8);
	    config.setStrictProgramCreation(false);
	    config.setCrossoverProb(75.0f);
        config.setMutationProb(25.0f);
        config.setReproductionProb(0.2f);
        
	    GPProblem problem = new GPSymbolicRegression(config);
	    GPGenotype gp = problem.create();
	    gp.setVerboseOutput(true);
	    
	    double bestFit = -1.0d;
	    // Termination criteria determined when number of generations reaches max or when fitness reaches 0
	    for (int gen = 1; gen <= genCount; gen++) {
	    	System.out.println("Generation Count: " + gen);
            gp.evolve();
            gp.calcFitness();
            double fitness = gp.getAllTimeBest().getFitnessValue();
            bestFit = fitness;
            if (bestFit == 0.0d) {
                break;
            }
        }
	  
	    System.out.println("---------------------------------------------------BEST SOLUTION--------------------------------------------------");
	    gp.outputSolution(gp.getAllTimeBest());
	    System.out.println("------------------------------------------------------------------------------------------------------------------");
	  }

	
}
