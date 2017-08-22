package examples.gp.GPClassification;

import java.io.FileNotFoundException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Scanner;

import org.jgap.InvalidConfigurationException;
import org.jgap.gp.CommandGene;
import org.jgap.gp.GPFitnessFunction;
import org.jgap.gp.GPProblem;
import org.jgap.gp.IGPProgram;
import org.jgap.gp.function.Add;
import org.jgap.gp.function.Add3;
import org.jgap.gp.function.Divide;
import org.jgap.gp.function.Multiply;
import org.jgap.gp.function.Subtract;
import org.jgap.gp.impl.DefaultGPFitnessEvaluator;
import org.jgap.gp.impl.GPConfiguration;
import org.jgap.gp.impl.GPGenotype;
import org.jgap.gp.terminal.Terminal;
import org.jgap.gp.terminal.Variable;

public class GPClassification extends GPProblem{

	private static final String trainingFile = "training.txt";
    private static final String testFile = "test.txt";
    
    private static Data training;
    
    private static int generationCount = 800;
    
    static List<Instance> patients = new ArrayList<Instance>();
    
    private static Variable CT;
    private static Variable USz;
    private static Variable UShp;
    private static Variable MA;
    private static Variable SESz;
    private static Variable BN;
    private static Variable BC;
    private static Variable NN;
    private static Variable M;
    
    public GPClassification(Data train, GPConfiguration conf) throws InvalidConfigurationException {
		super(conf);
	}

//--------------------------------------------------------------------------------------------------------------------\\ 
//--------------------------------------------------------------------------------------------------------------------\\  
    
    public static class Data{
    	public Data(String file) throws FileNotFoundException{
    		readFile(file);
    	}
    	
	    private static void readFile(String file) throws FileNotFoundException {
	    	
	    	Scanner scan = new Scanner(new InputStreamReader(ClassLoader.getSystemResourceAsStream(file)));
	    	
			HashMap<String, Integer> patientDetails = new HashMap<String, Integer>();
			
			while(scan.hasNextLine()){
				String line = scan.nextLine();
				String[] instanceAttributes = line.split(","); 
				int patientAttributes[] = new int[instanceAttributes.length]; 
				for(int i = 0; i < instanceAttributes.length; i++){
					if(instanceAttributes[i].equals("?")){
						instanceAttributes[i] = "-1";
					}
					 patientAttributes[i] = Integer.parseInt(instanceAttributes[i]);
				}
				patientDetails.put("ID", patientAttributes[0]);
				patientDetails.put("clumpThickness", patientAttributes[1]);
				patientDetails.put("uniformityOfCellSize", patientAttributes[2]);
				patientDetails.put("uniformityOfCellShape", patientAttributes[3]);
				patientDetails.put("marginalAdhesion", patientAttributes[4]);
				patientDetails.put("singleEpithelialCellSize", patientAttributes[5]);
				patientDetails.put("bareNuclei", patientAttributes[6]);
				patientDetails.put("blandChromatin", patientAttributes[7]);
				patientDetails.put("normalNucleoli", patientAttributes[8]);
				patientDetails.put("mitoses", patientAttributes[9]);
				patientDetails.put("cancerClassification", patientAttributes[10]);
				
				patients.add(new Instance(patientDetails));
				
				patientDetails.clear();
			}
			scan.close();
		}
	    public static List<Instance> getInstances() {
			return patients;
		} 
    }

//--------------------------------------------------------------------------------------------------------------------\\ 
//--------------------------------------------------------------------------------------------------------------------\\    
    
	@Override
	public GPGenotype create() 
			throws InvalidConfigurationException {
		GPConfiguration conf = getGPConfiguration();
		@SuppressWarnings("rawtypes")
		Class[] types = {
			CommandGene.IntegerClass,CommandGene.IntegerClass
		}; 
		@SuppressWarnings("rawtypes")
		Class[][] argTypes = {
			{}, {
			CommandGene.IntegerClass, CommandGene.IntegerClass, CommandGene.IntegerClass}
		};
		CommandGene[][] nodeSets = {{
			
			CT = Variable.create(conf, "CT", CommandGene.IntegerClass),
            USz = Variable.create(conf, "USz", CommandGene.IntegerClass),
            UShp = Variable.create(conf, "UShp", CommandGene.IntegerClass),
            MA = Variable.create(conf, "MA", CommandGene.IntegerClass),
            SESz = Variable.create(conf, "SESz", CommandGene.IntegerClass),
            BN = Variable.create(conf, "BN", CommandGene.IntegerClass),
            BC = Variable.create(conf, "BC", CommandGene.IntegerClass),
            NN = Variable.create(conf, "NN", CommandGene.IntegerClass),
            M = Variable.create(conf, "M", CommandGene.IntegerClass),
            
            new Multiply(conf, CommandGene.IntegerClass),
			new Add(conf, CommandGene.IntegerClass),
			new Divide(conf, CommandGene.IntegerClass),
			new Subtract(conf, CommandGene.IntegerClass),
			new Terminal(conf, CommandGene.IntegerClass, -10.0d, 10.0d, true),
		},{
            new Add3(conf, CommandGene.IntegerClass),
    }};
		return GPGenotype.randomInitialGenotype(conf, types, argTypes, nodeSets,
			20, true);
	}

	 @SuppressWarnings("serial")
	public static class FormulaFitnessFunction 
		extends GPFitnessFunction {
		 protected double evaluate(final IGPProgram a_subject) {
			 return computeRawFitness(a_subject); 
		 }

		 public static double computeRawFitness(final IGPProgram prog) {
			 int correct = 0;
			 for(Instance i : Data.getInstances()){
				 
				 CT.set(i.getClumpThickness());
	             USz.set(i.getUniformityOfCellSize());
	             UShp.set(i.getUniformityOfCellShape());
	             MA.set(i.getMarginalAdhesion());
	             SESz.set(i.getSingleEpithelialCellSize());
	             BN.set(i.getBareNuclei());
	             BC.set(i.getBlandChromatin());
	             NN.set(i.getNormalNucleoli());
	             M.set(i.getMitoses());
	             
				 int result = prog.execute_int(0, new Object[0]);
				 int classification = 0;
				 
				 if(result < 0){
					 classification = 4;}
				 else{
					 classification = 2;
				 }
				 if(classification == i.getCancerClassification()){
					 correct++;
				 }
			}
			double accuracy = ((double)correct/(double)Data.getInstances().size())*100;
			return accuracy; 
		 }
	 }
	 
	 public static void main(String[] args) throws Exception {
		 
			GPConfiguration config = new GPConfiguration();
			 config.setGPFitnessEvaluator(new DefaultGPFitnessEvaluator());	
			 config.setMaxInitDepth(4);
			 config.setPopulationSize(1000);
			 config.setMaxCrossoverDepth(8);
			 config.setFitnessFunction(new GPClassification.FormulaFitnessFunction());
			 config.setStrictProgramCreation(true);
			 config.setCrossoverProb(0.9f);
		     config.setMutationProb(0.1f);
		     config.setReproductionProb(0.1f);
		 
		 training = new Data(trainingFile);
		 GPProblem problem = new GPClassification(training, config);
		 GPGenotype gp = problem.create();
		 gp.setVerboseOutput(true);
		 for (int gen = 0; gen <= generationCount; gen++) {
			System.out.println("Generation Count: " + gen);
			gp.evolve();
			gp.calcFitness();
			if(gp.getAllTimeBest().getFitnessValue() >= 100)
				break;
		 }
		 System.out.println("---------------------------------------------------BEST SOLUTION--------------------------------------------------");
		 gp.outputSolution(gp.getAllTimeBest());
		 System.out.println("Correct Training Classifications: " + testResult(gp.getAllTimeBest(), trainingFile));
		 System.out.println("Correct Test Classifications: " + testResult(gp.getAllTimeBest(), testFile));
		 System.out.println("------------------------------------------------------------------------------------------------------------------");
	 }

	private static double testResult(IGPProgram allTimeBest, String file) throws FileNotFoundException{
		Data data = new Data(file);
		int correct = 0;
		for(Instance i : Data.getInstances()){
			CT.set(i.getClumpThickness());
	        USz.set(i.getUniformityOfCellSize());
	        UShp.set(i.getUniformityOfCellShape());
	        MA.set(i.getMarginalAdhesion());
	        SESz.set(i.getSingleEpithelialCellSize());
	        BN.set(i.getBareNuclei());
	        BC.set(i.getBlandChromatin());
	        NN.set(i.getNormalNucleoli());
	        M.set(i.getMitoses());
			int result = allTimeBest.execute_int(0, new Object[0]);
			int classification = 0;
			if(result < 0){
				classification = 4;
			}else{
				 classification = 2;
			}
			 if(classification == i.getCancerClassification()){
				 correct++;
			 }
		}
		double accuracy = ((double)correct/(double)Data.getInstances().size())*100;
		return accuracy;
	}
}

//--------------------------------------------------------------------------------------------------------------------\\ 
//--------------------------------------------------------------------------------------------------------------------\\   

    class Instance {
	private int ID;
	private int clumpThickness;
	private int uniformityOfCellSize;
	private int uniformityOfCellShape;
	private int marginalAdhesion;
	private int singleEpithelialCellSize;
	private int bareNuclei;
	private int blandChromatin;
	private int normalNucleoli;
	private int mitoses;
	private int cancerClassification;
	
	public Instance(HashMap<String, Integer> patientDetails) {
        this.ID = patientDetails.get("ID");
        this.clumpThickness = patientDetails.get("clumpThickness");
        this.uniformityOfCellSize = patientDetails.get("uniformityOfCellSize");
        this.uniformityOfCellShape = patientDetails.get("uniformityOfCellShape");
        this.marginalAdhesion = patientDetails.get("marginalAdhesion");
        this.singleEpithelialCellSize = patientDetails.get("singleEpithelialCellSize");
        this.bareNuclei = patientDetails.get("bareNuclei");
        this.blandChromatin = patientDetails.get("blandChromatin");
        this.normalNucleoli = patientDetails.get("normalNucleoli");
        this.mitoses = patientDetails.get("mitoses");
        this.cancerClassification = patientDetails.get("cancerClassification");
    }

	public int getID() {
		return ID;
	}
	public int getClumpThickness() {
		return clumpThickness;
	}
	public int getUniformityOfCellSize() {
		return uniformityOfCellSize;
	}
	public int getUniformityOfCellShape() {
		return uniformityOfCellShape;
	}
	public int getMarginalAdhesion() {
		return marginalAdhesion;
	}
	public int getSingleEpithelialCellSize() {
		return singleEpithelialCellSize;
	}
	public int getBareNuclei() {
		return bareNuclei;
	}
	public int getBlandChromatin() {
		return blandChromatin;
	}
	public int getNormalNucleoli() {
		return normalNucleoli;
	}
	public int getMitoses() {
		return mitoses;
	}
	public int getCancerClassification() {
		return cancerClassification;
	}
	
//--------------------------------------------------------------------------------------------------------------------\\ 
//--------------------------------------------------------------------------------------------------------------------\\   
	
}
