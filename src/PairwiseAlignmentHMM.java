/**
 * Created by user on 3/8/2017.
 */
//import com.sun.jmx.snmp.internal.SnmpSubSystem;

import java.lang.String;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.io.File;
import java.io.FileInputStream;

public class PairwiseAlignmentHMM {

	private  String File;
    private static int SeqLength = 0;
    private static int SeqXLength = 0;
	private static double[][] a_prob = new double[][]{    //a_prob[from][to]
		{0.838, 0.08, 0.08, 0.002},                     //from B to M, X, Y, E
		{0.838, 0.08, 0.08, 0.002},                     //from M to M, X, Y, E
		{0.648, 0.35, 0,    0.002},                     //from X to M, X, Y, E
		{0.648, 0,    0.35, 0.002}                      //from Y to M, X, Y, E
	};
	private static Double[] terminationScoresList = new Double[1000];

	private static double[] q_a = new double[]{             // q from qp.txt
			5.99e-02,
			5.60e-02,
			4.82e-02,
			5.22e-02,
			3.31e-02,
			5.02e-02,
			5.15e-02,
			5.97e-02,
			3.62e-02,
			5.96e-02,
			6.08e-02,
			5.33e-02,
			5.15e-02,
			5.02e-02,
			5.07e-02,
			5.47e-02,
			5.46e-02,
			1.93e-02,
			4.01e-02,
			5.83e-02};

	private static double[][] p_ab = new double[][]{        //p from qp.txt
		{1.36e-02, 2.26e-03, 1.37e-03, 1.49e-03, 1.89e-03, 2.02e-03, 2.07e-03, 3.40e-03, 1.03e-03, 2.40e-03, 2.45e-03, 2.15e-03, 2.07e-03, 1.43e-03, 2.04e-03, 4.41e-03, 3.11e-03, 3.89e-04, 1.14e-03, 3.32e-03},
		{2.26e-03, 1.69e-02, 2.57e-03, 1.39e-03, 6.23e-04, 3.78e-03, 2.74e-03, 1.59e-03, 1.93e-03, 1.12e-03, 1.62e-03, 5.68e-03, 1.94e-03, 9.46e-04, 1.35e-03, 2.06e-03, 2.06e-03, 3.64e-04, 1.07e-03, 1.10e-03},
		{1.37e-03, 2.57e-03, 1.77e-02, 3.39e-03, 5.37e-04, 2.30e-03, 2.36e-03, 2.74e-03, 2.35e-03, 9.67e-04, 9.86e-04, 2.45e-03, 1.18e-03, 8.14e-04, 1.16e-03, 3.55e-03, 2.51e-03, 2.21e-04, 9.21e-04, 9.46e-04},
		{1.49e-03, 1.39e-03, 3.39e-03, 2.07e-02, 5.81e-04, 2.49e-03, 5.11e-03, 2.10e-03, 1.27e-03, 1.05e-03, 7.55e-04, 1.87e-03, 9.04e-04, 8.81e-04, 1.78e-03, 2.71e-03, 1.92e-03, 2.40e-04, 7.05e-04, 1.02e-03},
		{1.89e-03, 6.23e-04, 5.37e-04, 5.81e-04, 2.36e-02, 5.59e-04, 4.05e-04, 6.65e-04, 4.03e-04, 1.33e-03, 1.35e-03, 5.93e-04, 1.15e-03, 7.90e-04, 5.64e-04, 1.22e-03, 1.22e-03, 3.04e-04, 6.32e-04, 1.30e-03},
		{2.02e-03, 3.78e-03, 2.30e-03, 2.49e-03, 5.59e-04, 1.36e-02, 4.92e-03, 1.43e-03, 1.73e-03, 1.01e-03, 1.45e-03, 3.60e-03, 2.46e-03, 8.48e-04, 1.71e-03, 2.61e-03, 1.85e-03, 4.61e-04, 1.36e-03, 1.39e-03},
		{2.07e-03, 2.74e-03, 2.36e-03, 5.11e-03, 4.05e-04, 4.92e-03, 1.43e-02, 1.46e-03, 1.77e-03, 1.03e-03, 1.05e-03, 3.69e-03, 1.26e-03, 8.69e-04, 1.76e-03, 2.68e-03, 1.89e-03, 3.34e-04, 9.83e-04, 1.43e-03},
		{3.40e-03, 1.59e-03, 2.74e-03, 2.10e-03, 6.65e-04, 1.43e-03, 1.46e-03, 2.71e-02, 1.03e-03, 8.47e-04, 8.64e-04, 1.51e-03, 1.03e-03, 1.01e-03, 1.44e-03, 3.11e-03, 1.55e-03, 5.49e-04, 8.07e-04, 1.17e-03},
		{1.03e-03, 1.93e-03, 2.35e-03, 1.27e-03, 4.03e-04, 1.73e-03, 1.77e-03, 1.03e-03, 1.99e-02, 7.26e-04, 7.40e-04, 1.30e-03, 8.87e-04, 1.22e-03, 8.73e-04, 1.33e-03, 9.41e-04, 3.32e-04, 2.77e-03, 7.10e-04},
		{2.40e-03, 1.12e-03, 9.67e-04, 1.05e-03, 1.33e-03, 1.01e-03, 1.03e-03, 8.47e-04, 7.26e-04, 1.35e-02, 6.90e-03, 1.07e-03, 4.13e-03, 2.85e-03, 1.02e-03, 1.55e-03, 2.19e-03, 3.87e-04, 1.61e-03, 9.36e-03},
		{2.45e-03, 1.62e-03, 9.86e-04, 7.55e-04, 1.35e-03, 1.45e-03, 1.05e-03, 8.64e-04, 7.40e-04, 6.90e-03, 1.41e-02, 1.54e-03, 5.96e-03, 2.90e-03, 1.04e-03, 1.58e-03, 2.24e-03, 5.58e-04, 1.64e-03, 4.77e-03},
		{2.15e-03, 5.68e-03, 2.45e-03, 1.87e-03, 5.93e-04, 3.60e-03, 3.69e-03, 1.51e-03, 1.30e-03, 1.07e-03, 1.54e-03, 1.53e-02, 1.85e-03, 9.00e-04, 1.82e-03, 2.77e-03, 1.96e-03, 3.46e-04, 1.02e-03, 1.48e-03},
		{2.07e-03, 1.94e-03, 1.18e-03, 9.04e-04, 1.15e-03, 2.46e-03, 1.26e-03, 1.03e-03, 8.87e-04, 4.13e-03, 5.96e-03, 1.85e-03, 1.43e-02, 2.46e-03, 1.24e-03, 1.89e-03, 1.89e-03, 6.69e-04, 1.39e-03, 4.04e-03},
		{1.43e-03, 9.46e-04, 8.14e-04, 8.81e-04, 7.90e-04, 8.48e-04, 8.69e-04, 1.01e-03, 1.22e-03, 2.85e-03, 2.90e-03, 9.00e-04, 2.46e-03, 1.92e-02, 6.05e-04, 1.31e-03, 1.31e-03, 1.30e-03, 5.42e-03, 1.97e-03},
		{2.04e-03, 1.35e-03, 1.16e-03, 1.78e-03, 5.64e-04, 1.71e-03, 1.76e-03, 1.44e-03, 8.73e-04, 1.02e-03, 1.04e-03, 1.82e-03, 1.24e-03, 6.05e-04, 2.77e-02, 1.86e-03, 1.86e-03, 2.33e-04, 6.85e-04, 1.41e-03},
		{4.41e-03, 2.06e-03, 3.55e-03, 2.71e-03, 1.22e-03, 2.61e-03, 2.68e-03, 3.11e-03, 1.33e-03, 1.55e-03, 1.58e-03, 2.77e-03, 1.89e-03, 1.31e-03, 1.86e-03, 1.14e-02, 4.02e-03, 3.55e-04, 1.04e-03, 1.52e-03},
		{3.11e-03, 2.06e-03, 2.51e-03, 1.92e-03, 1.22e-03, 1.85e-03, 1.89e-03, 1.55e-03, 9.41e-04, 2.19e-03, 2.24e-03, 1.96e-03, 1.89e-03, 1.31e-03, 1.86e-03, 4.02e-03, 1.61e-02, 5.02e-04, 1.04e-03, 3.03e-03},
		{3.89e-04, 3.64e-04, 2.21e-04, 2.40e-04, 3.04e-04, 4.61e-04, 3.34e-04, 5.49e-04, 3.32e-04, 3.87e-04, 5.58e-04, 3.46e-04, 6.69e-04, 1.30e-03, 2.33e-04, 3.55e-04, 5.02e-04, 1.60e-02, 1.47e-03, 3.79e-04},
		{1.14e-03, 1.07e-03, 9.21e-04, 7.05e-04, 6.32e-04, 1.36e-03, 9.83e-04, 8.07e-04, 2.77e-03, 1.61e-03, 1.64e-03, 1.02e-03, 1.39e-03, 5.42e-03, 6.85e-04, 1.04e-03, 1.04e-03, 1.47e-03, 1.74e-02, 1.58e-03},
		{3.32e-03, 1.10e-03, 9.46e-04, 1.02e-03, 1.30e-03, 1.39e-03, 1.43e-03, 1.17e-03, 7.10e-04, 9.36e-03, 4.77e-03, 1.48e-03, 4.04e-03, 1.97e-03, 1.41e-03, 1.52e-03, 3.03e-03, 3.79e-04, 1.58e-03, 1.30e-02}
	};

	//constructor

	public PairwiseAlignmentHMM(){

		this.File = this.getCurDirecoty() + "/src/2017-01-16_uniprot.fasta";


	}

	private static HashMap<String, Integer> hashMap = new HashMap<String, Integer>();

	//=================================================================================================================

    public static ArrayList globalViterbi(String[] Seq, String[] Seqx, int seqNum){
        //=====initialization=============
    	int xAxisLenth = Seqx.length+1;
    	int yAxisLenth = Seq.length+1;
    			
        String [][] Vm = new String[yAxisLenth][xAxisLenth];         //Vm[seq1][seq2]    [score : 0:M/ 1:X/ 2:Y]
        String [][] Vx = new String[yAxisLenth][xAxisLenth];         //Vx[seq1][seq2]
        String [][] Vy = new String[yAxisLenth][xAxisLenth];         //Vy[seq1][seq2]
        Vm[0][0] = "1: ";
        Vx[0][0] = "0: ";
        Vy[0][0] = "0: ";

        //Initialize the first column
		int product =1;
		for(int i=0; i<Seq.length;i++ ){
			if(Seq[i]!= null) {
				int ref1000 = hashMap.get(Seq[i]);
				double ini = Math.log(0.08 * (Math.pow(0.35, i - 1)) * (product * q_a[ref1000]));
				// round ini before toString()
				Vm[i+1][0] = String.valueOf(0).concat(":  ");
				Vx[i+1][0] = String.valueOf(ini).concat(":  ");
				Vy[i+1][0] = String.valueOf(0).concat(":  ");
			}else{
				break;
			}
		}
		//Initialize the first row
		product =1;
		for(int i=0; i<Seqx.length;i++ ){
			if(Seqx[i]!= null) {
				int refX = hashMap.get(Seqx[i]);
				double ini = Math.log(0.08 * (Math.pow(0.35, i - 1)) * (product * q_a[refX]));
				Vm[0][i+1] = String.valueOf(0).concat(":  ");
				Vx[0][i+1] = String.valueOf(0).concat(":  ");
				Vy[0][i+1] = String.valueOf(ini).concat(":  ");
			}else{
				break;
			}
		}

        //=====Recurrence=================
        for(int i =1; i<yAxisLenth; i++) {
		    for(int j=1; j<xAxisLenth; j++) {
		    	System.out.println("P" + i + " "+ j);
                if (Seq[i-1] != null && Seqx[j-1] != null) {
                    int ref1000 = hashMap.get(Seq[i-1]);
                    int refX = hashMap.get(Seqx[j-1]);
                    //=======Vm=====================
                    System.out.println("Vm");
                    double valvm;
                    double valvx;
                    double valvy;

                    valvm = Double.parseDouble(Vm[i - 1][j - 1].split(":")[0]);        //get the score Vm
                    valvx = Double.parseDouble(Vx[i - 1][j - 1].split(":")[0]);        //get the score Vx
                    valvy = Double.parseDouble(Vy[i - 1][j - 1].split(":")[0]);        //get the score Vy

                    double tempM = Math.log(Math.abs(Math.max(valvm * a_prob[1][0], Math.max(valvx * a_prob[2][0], valvy * a_prob[3][0]))));
                    if (tempM == Math.log(Math.abs(valvm * a_prob[1][0]))) {
                        double finalM = Math.log(p_ab[ref1000][refX]) + tempM;
                        Vm[i][j] = String.valueOf(finalM).concat(": 0");                    //store score and pointer to Vm
                    } else if (tempM == Math.log(Math.abs(valvx * a_prob[2][0]))) {
                        double finalX = Math.log(p_ab[ref1000][refX]) + tempM;
                        Vm[i][j] = String.valueOf(finalX).concat(": 1");                    //store score and pointer to Vx
                    } else if (tempM == Math.log(Math.abs(valvy * a_prob[3][0]))) {
                        double finalY = Math.log(p_ab[ref1000][refX]) + tempM;
                        Vm[i][j] = String.valueOf(finalY).concat(": 2");                    //store score and pointer to Vy
                    }
//                    System.out.println("valm" + valvm);            //-4.442079991122024
//                    System.out.println("valx" + valvx);
//                    System.out.println("valy" + valvy);
                    System.out.println("Vm: " + Vm[i][j]);

                    //=======Vx=====================
                    System.out.println("Vx");
                    valvm = Double.parseDouble(Vm[i - 1][j].split(":")[0]);        //get the score Vm
                    valvx = Double.parseDouble(Vx[i - 1][j].split(":")[0]);        //get the score Vx

                    double tempX = Math.log(Math.abs(Math.max(valvm * a_prob[1][1], valvx * a_prob[2][1])));
                    if (tempX == Math.log(Math.abs(valvm * a_prob[1][1]))) {
                        double finalM = Math.log(q_a[ref1000]) + tempX;
                        Vx[i][j] = String.valueOf(finalM).concat(": 0");                    //store score and pointer to Vm
                    } else if (tempX == Math.log(Math.abs(valvx * a_prob[2][1]))) {
                        double finalX = Math.log(q_a[ref1000]) + tempX;
                        Vx[i][j] = String.valueOf(finalX).concat(": 1");                    //store score and pointer to Vx
                    }
                    System.out.println("Vx: " + Vx[i][j]);
                    //=======Vy=====================
                    System.out.println("Vy");
                    valvm = Double.parseDouble(Vm[i][j - 1].split(":")[0]);        //get the score Vm
                    valvy = Double.parseDouble(Vy[i][j - 1].split(":")[0]);        //get the score Vy

                    double tempY = Math.log(Math.abs(Math.max(valvm * a_prob[1][2], valvy * a_prob[3][2])));
                    if (tempY == Math.log(Math.abs(valvm * a_prob[1][2]))) {
                        double finalM = Math.log(q_a[refX]) + tempY;
                        Vy[i][j] = String.valueOf(finalM).concat(": 0");                    //store score and pointer to Vm
                    } else if (tempY == Math.log(Math.abs(valvy * a_prob[3][2]))) {
                        double finalY = Math.log(q_a[refX]) + tempY;
                        Vy[i][j] = String.valueOf(finalY).concat(": 2");                    //store score and pointer to Vy
                    }
                    System.out.println("Vy: " + Vy[i][j]);

                }else{
                    break;
                }
            }
        }

        //=====Termination================
        double maxFromVm = Double.parseDouble(Vm[Seq.length][Seqx.length].split(":")[0]);
        double maxFromVx = Double.parseDouble(Vx[Seq.length][Seqx.length].split(":")[0]);
        double maxFromVy = Double.parseDouble(Vy[Seq.length][Seqx.length].split(":")[0]);
        
        double termination = Math.max(maxFromVm, Math.max(maxFromVx,maxFromVy))* 0.002;
        terminationScoresList[seqNum]=termination;
        System.out.println("termination");

        ArrayList res = new ArrayList(3);
        res.add(Vm);
        res.add(Vx);
        res.add(Vy);
        
        return res;

    }

	public static String getCurDirecoty (){

		String workingDir = System.getProperty("user.dir");
		System.out.println(workingDir);
		return workingDir;

	}
	
	
	public static ArrayList<Integer> indexesOfTopElements(Double[]terminationScoresList, int nummax) {
        Double[] copy =  Arrays.copyOf(terminationScoresList,1000);
        Arrays.sort(terminationScoresList);
        Double[] honey = Arrays.copyOfRange(copy,copy.length - nummax, copy.length);
        ArrayList<Integer> result = new ArrayList<Integer>(nummax);
        int resultPos = 0;
        for(int i = 0; i < terminationScoresList.length; i++) {
            double onTrial = terminationScoresList[i];
            int index = Arrays.binarySearch(honey,onTrial);
            if(index < 0) continue;
            resultPos++;
            result.set(resultPos, i);
        }
        return result;
    }
	
	public static String[] doTraceback(String[] Seq, String seq2, int seqNum){
		String[] seq2CharList = seq2.split("");
		StringBuilder templateAlignment = new StringBuilder();
		StringBuilder sequenceAlignment = new StringBuilder();
		ArrayList matricesObject =  globalViterbi(Seq,seq2CharList, seqNum);
		String [][] Vm = (String[][]) matricesObject.get(0);         //Vm[seq1][seq2]    [score : 0:M/ 1:X/ 2:Y]
	    String [][] Vx = (String[][]) matricesObject.get(1);         //Vx[seq1][seq2]
	    String [][] Vy = (String[][]) matricesObject.get(2);
		
		// TODO: do the traceback here jumping bwtween the three matrices
	    double maxFromVm = Double.parseDouble(Vm[Seq.length][seq2CharList.length].split(":")[0]);
	    String FromWhichVm = Vm[Seq.length][seq2CharList.length].split(":")[1].trim();
        double maxFromVx = Double.parseDouble(Vx[Seq.length][seq2CharList.length].split(":")[0]);
        String FromWhichVx = Vx[Seq.length][seq2CharList.length].split(":")[1].trim();
        double maxFromVy = Double.parseDouble(Vy[Seq.length][seq2CharList.length].split(":")[0]);
        String FromWhichVy = Vy[Seq.length][seq2CharList.length].split(":")[1].trim();
       
        double termination = Math.max(maxFromVm, Math.max(maxFromVx,maxFromVy));
        String priorMatrix = "";
        int i = seq2CharList.length-1;
        int j = Seq.length-1;
        //Initialize the alignment
        if (termination == maxFromVm){
        	templateAlignment.append(Seq[j]);
        	sequenceAlignment.append(seq2CharList[i]);
        	i -=1;
        	j -=1;
        	priorMatrix = FromWhichVm;
        	
        }else if (termination == maxFromVx){ //Delection xi,-
        	templateAlignment.append(Seq[j]);
        	sequenceAlignment.append("-");
        	j -=1;
        	priorMatrix = FromWhichVx;
        }else if (termination == maxFromVy){
        	templateAlignment.append("-");
        	sequenceAlignment.append(seq2CharList[i]);
        	i-=1;
        	priorMatrix = FromWhichVy;
        }
        
        
        while (i!=0 || j !=0){
        	if (i != 0) {
        		if (j != 0) {
        			if (priorMatrix.equals("0")){
                    	templateAlignment.append(Seq[j]);
                    	sequenceAlignment.append(seq2CharList[i]);
                    	i -=1;
                    	j -=1;
                    	priorMatrix = Vm[j][i].split(":")[1].trim();
                    	
                    }else if (priorMatrix.equals("1")){ //Delection xi,-
                    	templateAlignment.append(Seq[j]);
                    	sequenceAlignment.append("-");
                    	j -=1;
                    	priorMatrix = Vx[j][i].split(":")[1].trim();
                    }else if (priorMatrix.equals("2")){
                    	templateAlignment.append("-");
                    	sequenceAlignment.append(seq2CharList[i]);
                    	i-=1;
                    	priorMatrix = Vy[j][i].split(":")[1].trim();
                    }
        		} else { // i != 0 and j == 0
        			templateAlignment.append("-");
                	sequenceAlignment.append(seq2CharList[i]);
                	i-=1;
        		}
        	} else { // i == 0 and j != 0
				templateAlignment.append(Seq[j]);
				sequenceAlignment.append("-");
				j -=1;
        	}
        	
        	
        }
        
        String[] answer = new String[2];
        answer[0] = sequenceAlignment.toString();
        answer[1] = templateAlignment.toString();
		return answer;
		
	}

	//=================================================================================================================

	public static void main(String[]args){
		
		PairwiseAlignmentHMM pwa = new PairwiseAlignmentHMM();
		
		BufferedReader br = null;
		FileReader fr = null;
		String[] Seq = new String[1001];
		String[] Header = new String[1001];

		hashMap.put("A", 0);
		hashMap.put("R", 1);
		hashMap.put("N", 2);
		hashMap.put("D", 3);
		hashMap.put("C", 4);
		hashMap.put("Q", 5);
		hashMap.put("E", 6);
		hashMap.put("G", 7);
		hashMap.put("H", 8);
		hashMap.put("I", 9);
		hashMap.put("L", 10);
		hashMap.put("K", 11);
		hashMap.put("M", 12);
		hashMap.put("F", 13);
		hashMap.put("P", 14);
		hashMap.put("S", 15);
		hashMap.put("T", 16);
		hashMap.put("W", 17);
		hashMap.put("Y", 18);
		hashMap.put("V", 19);

		//Read File and add Seq to array
		try{
			fr = new FileReader(pwa.File);        //change to File
			br = new BufferedReader(fr);
			

			String cl;
			boolean found = false;
			br = new BufferedReader(new FileReader(pwa.File));         //change to File
			int j=-1;
			while((cl = br.readLine()) != null){
				if(cl.startsWith(">")) {
					found = false;
					j+=1;
				}
				if(found){
					if(Seq[j]!=null){
						String temp = Seq[j];
						Seq[j]=temp.concat(cl);
					}else{
						Seq[j] = cl;
					}
					//                    System.out.println(Seq[j]);

				}
				if(j>1000) break;
				if(cl.startsWith(">"))
				{   if(found) break;
				Header[j] = cl;
				found=true;
				//                    System.out.println(Header[j]);

				}

			}
			br.close();
		}catch(IOException e){
			e.printStackTrace();
		}finally {
			try{
				if(br != null)
					br.close();
				if(fr != null)
					fr.close();
			}catch(IOException ex){
				ex.printStackTrace();
			}

//			String Seq1000 = Seq[1000];
  /*          String Seq1000 = Seq[59];
			System.out.println(Seq1000);
			String[] Sequence1000 = new String[Seq1000.length()];
			int i =0;
			for (String retval: Seq1000.split("")) {
				Sequence1000[i] = retval;
				i++;
				SeqLength++;
			}
			//call functions
			for(int j=0; j<1000; j++){
				String SeqIndX = Seq[j];
				String[] SequenceX = new String[SeqIndX.length()];
				int k=0;
				for(String ret: SeqIndX.split("")){
					SequenceX[k] = ret;
					k++;
					SeqXLength++;
				}
				System.out.println("sequence100 "+ Sequence1000.length);
				System.out.println("SequenceX "+ SequenceX.length);

				//            System.out.println(k);
				pwa.globalViterbi(Sequence1000, SequenceX, j);


			}

			ArrayList<Integer> maxThreeAlignmentIndeces = new ArrayList<Integer>();
			
			maxThreeAlignmentIndeces = pwa.indexesOfTopElements(pwa.terminationScoresList,3);
			
			for (int j = 0; i < maxThreeAlignmentIndeces.size();i++){
				int index = maxThreeAlignmentIndeces.get(j);
				String name = Header[index];
				double value = pwa.terminationScoresList[index];
				System.out.println("Index=" + index + " Name=" + name +  " ln Pr="+ value);
				
				String alignment = pwa.doTraceback(Sequence1000, Seq[index], j);
			}
*/

  String[] sampleX = new String[]{"T", "A", "P","P", "A","C"};
  String[] sampleY = new String[]{"T","A","A","C"};;
  String sampleYString = "TAAC";
  pwa.globalViterbi(sampleX,sampleY,1);
  String[] alignment = pwa.doTraceback(sampleX, sampleYString, 1);
		}
	}

}
