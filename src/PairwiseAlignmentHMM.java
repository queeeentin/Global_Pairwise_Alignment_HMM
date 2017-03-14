/**
 * Created by user on 3/8/2017.
 */
import java.lang.String;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;

public class PairwiseAlignmentHMM {

	private  String File;
	private  double[][] a_prob = new double[][]{    //a_prob[from][to]
		{0.838, 0.08, 0.08, 0.002},                     //from B to M, X, Y, E
		{0.838, 0.08, 0.08, 0.002},                     //from M to M, X, Y, E
		{0.648, 0.35, 0,    0.002},                     //from X to M, X, Y, E
		{0.648, 0,    0.35, 0.002}                      //from Y to M, X, Y, E
	};

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
		this.File = this.getCurDirecoty() + "\\2017-01-16 uniprot.fasta";
	}

	private static HashMap<String, Integer> hashMap = new HashMap<String, Integer>();

	//=================================================================================================================

	public void globalViterbi(String[] Seq, String[] Seqx){
		//=====initialization=============
		String [][] Vm = new String[Seq.length][Seqx.length];         //Vm[seq1][seq2]    [score : 0:M/ 1:X/ 2:Y]
		String [][] Vx = new String[Seq.length][Seqx.length];         //Vx[seq1][seq2]
		String [][] Vy = new String[Seq.length][Seqx.length];         //Vy[seq1][seq2]
		Vm[0][0] = "1";

		//=====Recurrence=================

		//=====Termination================
	}

	public static String getCurDirecoty (){

		String workingDir = System.getProperty("user.dir");
		System.out.println(workingDir);
		return workingDir;

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

			String Seq1000 = Seq[900];
			String[] Sequence1000 = new String[900];
			int i =0;
			for (String retval: Seq1000.split("")) {
				Sequence1000[i] = retval;
				i++;
			}
			//call functions
			for(int j=0; j<1000; j++){
				String SeqIndX = Seq[j];
				String[] SequenceX = new String[900];
				int k=0;
				for(String ret: SeqIndX.split("")){
					SequenceX[k] = ret;
					k++;
				}
				//            System.out.println(k);
				//new Assg4Q3().globalViterbi(Sequence1000, SequenceX);

			}
		}
	}

}
