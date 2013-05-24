import java.util.*;
import java.io.*;
import java.lang.*;
public class AlignThreeSequence {

  public static int [][] tranpose(int [][] M) {
    int nrows = M.length;
    int ncols = M[0].length;
    
    int [][] N = new int [ncols][nrows];
    for (int i = 0; i < nrows; ++i){
      for (int j = 0; j < ncols; ++j) {
        N[j][i] = M[i][j];
      }
    }
    return N;
  }

  public static int [][] ClockWiseRotate90Degree(int [][]M){
    int nrows = M.length;
    int ncols = M[0].length;
    int [][] N = new int [ncols][nrows];
    for (int i = 0; i < nrows; ++i){
      for (int j = 0; j < ncols; ++j) {
        N[j][N[0].length - 1 - i] = M[i][j];
      }
    }
    return N;
  }  
  
  public static String join(int [] L, String delimiter){
	  StringBuilder result = new StringBuilder();
	  for(int i = 0; i < L.length - 1; ++i){
		  result.append(L[i]);
		  result.append(delimiter);
	  }
	  result.append(L[L.length - 1]);
	  return result.toString(); 
  }	

  public static String print2DMatrix(int [][] M){
	  int nrows = M.length;
	  int ncols = M[0].length;
	  String result = "";
	  for (int i = 0; i < M.length; ++i){
		  result += join(M[i], "\t") + "\n";
	  } 
	  return result;
  }

  public static String print3DMatrix(int [][][] M){
	  int nrows = M.length;
	  int ncols = M[0].length;
	  String result = "";
	  for (int i = 0; i < M.length; ++i){
		  result += "#" + i + ":\n";
		  for (int j = 0; j < M[0].length; ++j){
			  result += join(M[i][j], "\t") + "\n";
		  }
		  result += "\n";
	  } 
	  return result;
  }
  
  public static void PrintHash(HashMap m){
    Set set = m.entrySet(); 
    Iterator i = set.iterator();
    while (i.hasNext()){
      Map.Entry me = (Map.Entry) i.next();
      System.out.println(me.getKey() + " => " + me.getValue());
    }
  } 

  public static HashMap ScoreMatrixTable() {
    int mismatch = -4;
    int match = 5;
    int indel = -8;
    HashMap ScoreTwoChar = new HashMap();
    HashMap ScoreThreeChar = new HashMap();

    String [] alphabet = {"A", "C", "G", "T", "-"};
    for (String a: alphabet){
      for (String b: alphabet){
        if (a == b)
          if (a == "-")
            ScoreTwoChar.put(a+b, 0);
          else
            ScoreTwoChar.put(a+b, match);
        else
          if (a == "-" || b == "-")
           ScoreTwoChar.put(a+b, indel); 
          else
           ScoreTwoChar.put(a+b, mismatch); 
     }
    }

    for (String a: alphabet){
      for (String b: alphabet){
        for (String c: alphabet){
          if (a == b && b == c && c == "-")
            continue;
          else {
           int sum = (Integer)ScoreTwoChar.get(a+b) 
             + (Integer)ScoreTwoChar.get(a+c) 
             + (Integer)ScoreTwoChar.get(b+c);
            ScoreThreeChar.put(a+b+c, sum);
          }
            
        }
      }
    }

//    PrintHash(ScoreThreeChar);
    return ScoreThreeChar;
  }  
  
  public static void testScoreTwoSeq(String A, String B){
	  NamedObject2D result = scoreTwoSeqAB(A, B);
	  System.out.println("Score Matrix\n" + print2DMatrix(result.Score));
  }
  
  public static void testAlignThreeSeq(String A, String B, String C){
	  
	  NamedObject3D result = alignThreeSequence(A, B, C);
	  System.out.format("%s\n%s\n%s\n", result.A_alignment, result.B_alignment, result.C_alignment);
  }
  
  public static void test(){
    int [][] test = new int [3][4];
    int k = 0;
//    System.out.println("dimension of Array = (" + test.length + ", " + test[0].length + ")");
    for(int i = 0; i < test.length; ++i){
      for(int j = 0; j < test[0].length; ++j){
        test[i][j] = k++;
      }
    }
//    System.out.println("original:\n" + printMatrix(test));
//    System.out.println("Transposed:\n" + printMatrix(tranpose(test)));
//    System.out.println("ClockWiseRotate90Degree:\n" + printMatrix(ClockWiseRotate90Degree(test)));
    PrintHash(ScoreMatrixTable());
    testScoreTwoSeq("ATC", "ATC");
    testScoreTwoSeq("ATC", "ATG");
    testAlignThreeSeq("ATC", "ATC", "ATC");
  }
  
  public static class NamedObject2D {
    // encapsulate class
	  public final int [][] Score;
	  public final int [][][] Path;
	  public NamedObject2D(int [][] s, int [][][] p){
		  this.Score = s;
		  this.Path = p;
	  }
  }
  
  public static class NamedObject3D {
    // encapsulate class
    public final String A_alignment;
    public final String B_alignment;
    public final String C_alignment;
    public final int  Score;
    
    public NamedObject3D (String A, String B, String C, int  S) {
      this.A_alignment= A;
      this.B_alignment = B;
      this.C_alignment = C;
      this.Score = S;
    }
  }
  
 public static int [][] scoreThreeSeq(String A, String B, String C){
	 HashMap score = ScoreMatrixTable();
	 int new_indel = (Integer) score.get("A" + "-" + "-");
	 int new_match = (Integer) score.get("A" + "A" + "-");
	 int new_mismatch = (Integer) score.get("A" + "C" + "-");
	 NamedObject2D result = scoreTwoSeqBC(B, C);
	 int [][] prev = result.Score;
	 int [][] current = new int [B.length() + 1] [C.length() + 1];
	 
	 for (int a = 1; a < A.length() + 1; ++a){
		 current[0][0] = prev[0][0] + new_indel;
		 // fill dummy row when C == 0
		 for(int b = 1; b < B.length() + 1; ++b){
			 int one = current[b-1][0] + new_indel;
			 int diagonal = prev[b-1][0] + 
					 (Integer) score.get(String.valueOf(B.charAt(b-1)) 
					 + String.valueOf(A.charAt(a-1)) + "-");
			 int two = prev[b][0] + new_indel;
			 ArrayList<Integer> tmp = new ArrayList<Integer>();
			 tmp.add(one); tmp.add(two); tmp.add(diagonal);
			 current[b][0] = Collections.max(tmp);
		 }
			 // fill dummy row when B == 0
		 for(int c = 1; c < C.length() + 1; ++c){
			 int one = current[0][c-1] + new_indel;
			 int diagonal = prev[0][c-1] + 
					 (Integer) score.get(String.valueOf(C.charAt(c-1)) 
					 + String.valueOf(A.charAt(a-1)) + "-");
			 int two = prev[0][c] + new_indel;
			 ArrayList<Integer> tmp = new ArrayList<Integer>();
			 tmp.add(one); tmp.add(two); tmp.add(diagonal);
			 current[0][c] = Collections.max(tmp);
		 }	 
		 for(int b = 1; b < B.length() + 1; ++b){	 
			 for(int c = 1; c < C.length() + 1; ++c){
				 int seven = prev[b-1][c-1] +  
						 (Integer) score.get(String.valueOf(A.charAt(a-1)) + 
								 String.valueOf(B.charAt(b-1)) + 
								 String.valueOf(C.charAt(c-1)));
				 int six   = prev[b-1][c] + 
						 (Integer) score.get(String.valueOf(A.charAt(a-1)) + 
						 String.valueOf(B.charAt(b-1)) + "-");
				 int five  = prev[b][c-1] + 
						 (Integer) score.get(String.valueOf(A.charAt(a-1)) + "-" 
						 + String.valueOf(C.charAt(c-1)));
				 int three = prev[b][c] + 
						 (Integer) score.get(A.charAt(a-1) + "-" + "-");
				 int four  = current[b-1][c-1] + 
						 (Integer) score.get("-" + String.valueOf(B.charAt(b-1)) 
						 + String.valueOf(C.charAt(c-1)));
				 int two   = current[b-1][c] + 
						 (Integer) score.get("-" +  String.valueOf(B.charAt(b-1)) + "-");
				 int one   = current[b][c-1] + 
						 (Integer) score.get("-" +  "-" + String.valueOf(C.charAt(c-1)));
				 ArrayList<Integer> tmp = new ArrayList<Integer>();
                 tmp.add(seven); tmp.add(six); tmp.add(five); tmp.add(four); tmp.add(three);
                 tmp.add(two); tmp.add(one);
                 current[b][c] = Collections.max(tmp);
			 }
		 }
		 prev = current;
	 }
	 return current;
 }
 
 public static NamedObject2D scoreTwoSeqBC(String B, String C){
	 int [][] M = new int [B.length() + 1][C.length() + 1];
	 int [][][] P = new int [B.length() + 1][C.length() + 1][3];
	 HashMap score = ScoreMatrixTable();
	 int new_indel = (Integer) score.get("A" + "-" + "-");
	 int new_match = (Integer) score.get("A" + "A" + "-");
	 int new_mismatch = (Integer) score.get("A" + "C" + "-");
	 for (int i = 1; i < C.length() + 1; ++i){
 		 M[0][i] = M[0][i-1] + new_indel;
 		 P[0][i][0] = 0;
 		 P[0][i][1] = 0;
 		 P[0][i][2] = i-1;
 	 }
	 // initialize the dummy column
 	 for (int i = 1; i < B.length() + 1; ++i){ 
 		 M[i][0] = M[i-1][0] + new_indel;
 		 P[0][i][0] = 0;
 		 P[0][i][1] = i-1;
 		 P[0][i][2] = 0;
 	 }
 	 
 	 for (int i = 1; i < B.length() + 1; ++i){
 		 for (int j = 1; j < C.length() + 1; ++j){
 			 if (B.charAt(i-1) == C.charAt(j-1)){
 				 M[i][j] =  Math.max(Math.max(M[i-1][j-1] + new_match, 
 						 M[i-1][j] + new_indel), 
 						 M[i][j-1] + new_indel);

 			 } 
 			 else{
 				 M[i][j] =  Math.max(Math.max(M[i-1][j-1] + new_mismatch, 
 						 M[i-1][j] + new_indel), 
 						 M[i][j-1] + new_indel);				
 			 }
 			 int first, second;
 			 if (M[i][j] == M[i-1][j-1] + new_match || M[i][j] == M[i-1][j-1] + new_mismatch){
 				 first = i-1; second = j-1;
 			 }
 			 else if (M[i][j] == M[i-1][j] + new_indel){
 				 first = i-1; second = j;
 			 }
 			 else{ 
 				 first = i; second = j-1;
 			 }
				 P[i][j][0] = 0;
				 P[i][j][1] = first;
				 P[i][j][2] = second;
 		 }
 	 }
 	return new NamedObject2D(M, P);
 	 
 }
 	public static NamedObject2D scoreTwoSeqAB(String A, String B){
	 int [][] M = new int [A.length() + 1][B.length() + 1];
	 // P is used to store the 3D coordinate for previous position, where the alignment
	 // comes from
	 int [][][] P = new int [A.length() + 1][B.length() + 1][3];
	 HashMap score = ScoreMatrixTable();
	 int new_indel = (Integer) score.get("A" + "-" + "-");
	 int new_match = (Integer) score.get("A" + "A" + "-");
	 int new_mismatch = (Integer) score.get("A" + "C" + "-");
	 
	// initialize the dummy row 
 	 for (int i = 1; i < B.length() + 1; ++i){
 		 M[0][i] = M[0][i-1] + new_indel;
 		 P[0][i][0] = 0;
 		 P[0][i][1] = i-1;
 		 P[0][i][2] = 0;
 	 }
	 // initialize the dummy column
 	 for (int i = 1; i < A.length() + 1; ++i){ 
 		 M[i][0] = M[i-1][0] + new_indel;
 		 P[0][i][0] = i-1;
 		 P[0][i][1] = 0;
 		 P[0][i][2] = 0;
 	 }
 	 
 	 for (int i = 1; i < A.length() + 1; ++i){
 		 for (int j = 1; j < B.length() + 1; ++j){
 			 if (A.charAt(i-1) == B.charAt(j-1)){
 				 M[i][j] =  Math.max(Math.max(M[i-1][j-1] + new_match, 
 						 M[i-1][j] + new_indel), 
 						 M[i][j-1] + new_indel);

 			 } 
 			 else{
 				 M[i][j] =  Math.max(Math.max(M[i-1][j-1] + new_mismatch, 
 						 M[i-1][j] + new_indel), 
 						 M[i][j-1] + new_indel);				
 			 }
 			 int first, second;
 			 if (M[i][j] == M[i-1][j-1] + new_match || M[i][j] == M[i-1][j-1] + new_mismatch){
 				 first = i-1; second = j-1;
 			 }
 			 else if (M[i][j] == M[i-1][j] + new_indel){
 				 first = i-1; second = j;
 			 }
 			 else{ 
 				 first = i; second = j-1;
 			 }
 			 P[i][j][0] = first;
 			 P[i][j][1] = second;
 			 P[i][j][2] = 0; 		
 		 }
 	 }
 	 return new NamedObject2D(M, P);
  }
 	
	public static NamedObject2D scoreTwoSeqAC(String A, String C){
		 int [][] M = new int [A.length() + 1][C.length() + 1];
		 // P is used to store the 3D coordinate for previous position, where the alignment
		 // comes from
		 int [][][] P = new int [A.length() + 1][C.length() + 1][3];
		 HashMap score = ScoreMatrixTable();
		 int new_indel = (Integer) score.get("A" + "-" + "-");
		 int new_match = (Integer) score.get("A" + "A" + "-");
		 int new_mismatch = (Integer) score.get("A" + "C" + "-");
		 
		// initialize the dummy row 
	 	 for (int i = 1; i < C.length() + 1; ++i){
	 		 M[0][i] = M[0][i-1] + new_indel;
	 		 P[0][i][0] = 0;
	 		 P[0][i][1] = 0;
	 		 P[0][i][2] = i-1;
	 	 }
		 // initialize the dummy column
	 	 for (int i = 1; i < A.length() + 1; ++i){ 
	 		 M[i][0] = M[i-1][0] + new_indel;
	 		 P[0][i][0] = i-1;
	 		 P[0][i][1] = 0;
	 		 P[0][i][2] = 0;
	 	 }
	 	 
	 	 for (int i = 1; i < A.length() + 1; ++i){
	 		 for (int j = 1; j < C.length() + 1; ++j){
	 			 if (A.charAt(i-1) == C.charAt(j-1)){
	 				 M[i][j] =  Math.max(Math.max(M[i-1][j-1] + new_match, 
	 						 M[i-1][j] + new_indel), 
	 						 M[i][j-1] + new_indel);

	 			 } 
	 			 else{
	 				 M[i][j] =  Math.max(Math.max(M[i-1][j-1] + new_mismatch, 
	 						 M[i-1][j] + new_indel), 
	 						 M[i][j-1] + new_indel);				
	 			 }
	 			 int first, second;
	 			 if (M[i][j] == M[i-1][j-1] + new_match || M[i][j] == M[i-1][j-1] + new_mismatch){
	 				 first = i-1; second = j-1;
	 			 }
	 			 else if (M[i][j] == M[i-1][j] + new_indel){
	 				 first = i-1; second = j;
	 			 }
	 			 else{ 
	 				 first = i; second = j-1;
	 			 }
	 			 P[i][j][0] = first;
	 			 P[i][j][1] = 0;
	 			 P[i][j][2] = second; 		
	 		 }
	 	 }
	 	 return new NamedObject2D(M, P);
	  }
  public static NamedObject3D alignThreeSequence (String A, String B, String C) {
	 HashMap score = ScoreMatrixTable(); 
	 int [][][] M = new int [A.length() + 1] [B.length() + 1][C.length() + 1];
	 int [][][][] P = new int [A.length() + 1] [B.length() + 1][C.length() + 1][3];
	 
	 // convert 2D tuple to 3D tuple dummy surface A
	 NamedObject2D result = scoreTwoSeqBC(B, C);
	 M[0] = result.Score;
	 P[0] = result.Path;

	 
	 // convert 2D tuple to 3D tuple dummy surface B
	 NamedObject2D AC = scoreTwoSeqAC(A, C);
	 for (int i = 0; i < A.length() + 1; ++i){
		 for (int j = 0; j < C.length() + 1; ++j){
			 M[i][0][j] = AC.Score[i][j];
			 P[i][0][j] = AC.Path[i][j];
		 }
	 }
	 // convert 2D tuple to 3D tuple dummy surface C
	 NamedObject2D AB = scoreTwoSeqAB(A, B);
	 for (int i = 0; i < A.length() + 1; ++i){
		 for (int j = 0; j < B.length() + 1; ++j){
			 M[i][j][0] = AB.Score[i][j];
			 P[i][j][0] = AB.Path[i][j];
		 }
	 }
	 
	 for (int a = 1; a < A.length() + 1; ++a){
		 for (int b = 1; b < B.length() + 1; ++b){
			 for (int c = 1; c < C.length() + 1; ++c){
                 int seven = M[a-1][b-1][c-1] +  (Integer) score.get(String.valueOf(A.charAt(a-1)) + 
                		 String.valueOf(B.charAt(b-1)) + String.valueOf(C.charAt(c-1)));
                 int six   = M[a-1][b-1][c] + (Integer) score.get(String.valueOf(A.charAt(a-1)) + 
                		 String.valueOf(B.charAt(b-1)) + "-");
                 int five  = M[a-1][b][c-1] + (Integer) score.get(String.valueOf(A.charAt(a-1)) + "-" 
                		 + String.valueOf(C.charAt(c-1)));
                 int three = M[a-1][b][c] + (Integer) score.get(A.charAt(a-1) + "-" + "-");
                 int four  = M[a][b-1][c-1] + (Integer) score.get("-" + String.valueOf(B.charAt(b-1)) 
                		 + String.valueOf(C.charAt(c-1)));
                 int two   = M[a][b-1][c] + (Integer) score.get("-" +  String.valueOf(B.charAt(b-1)) + "-");
                 int one   = M[a][b][c-1] + (Integer) score.get("-" +  "-" + String.valueOf(C.charAt(c-1)));			 

                 ArrayList<Integer> tmp = new ArrayList<Integer>();
                 tmp.add(seven); tmp.add(six); tmp.add(five); tmp.add(four); tmp.add(three);
                 tmp.add(two); tmp.add(one);
			     M[a][b][c] = Collections.max(tmp);
//			     System.out.println("Here is A-----: ");
			     int first, second, third;
                 if ( M[a][b][c] == seven){
                	first = a-1; second = b-1; third = c-1;
                 }
                 else if ( M[a][b][c] == six){
                  	first = a-1; second = b-1; third = c;
                 }
                 else if (M[a][b][c] == four){
                   	first = a; second = b-1; third = c-1;
                 }
                 else if (M[a][b][c] == five){
                	first = a-1; second = b; third = c-1;
                 }
                 else if (M[a][b][c]  == one ){
                  	first = a; second = b; third = c-1;
                 }
                 else if (M[a][b][c] == two ){
                   	first = a; second = b-1; third = c;
                 }
                 else { // three
                	first = a-1; second = b; third = c;
                 }
                 P[a][b][c][0] = first;
                 P[a][b][c][0] = second;
                 P[a][b][c][0] = third;
			 }
		 }
	 }

	 // backtracing
    StringBuilder A_aln = new StringBuilder();
    StringBuilder B_aln = new StringBuilder();
    StringBuilder C_aln = new StringBuilder();
    int a = A.length(); int b = B.length(); int c = C.length();
    
    int max_score = M[a][b][c];
    // continue until meet the start point
    while (!(P[a][b][c][0] == 0 && P[a][b][c][1] == 0 && P[a][b][c][2] == 0)){
    	
      int new_a = P[a][b][c][0];
      int new_b = P[a][b][c][1];
      int new_c = P[a][b][c][2];
      System.out.println("Here is B------------: ");
      A_aln.insert(0, a - new_a > 0 ? A.charAt(a-1) : '-');
      B_aln.insert(0, b - new_b > 0 ? B.charAt(b-1) : '-');
      C_aln.insert(0, c - new_c > 0 ? C.charAt(c-1) : '-');
      a = new_a; b = new_b; c = new_c;
    }
    System.out.println("Here is A: " + A_aln.toString());
    System.out.println("Here is B: " + B_aln.toString());
    System.out.println("Here is C: " + C_aln.toString());
	return new NamedObject3D (A_aln.toString(), B_aln.toString(), C_aln.toString(),max_score); 
  }
  
  public static String reverseStr(String A){
	  return new StringBuffer(A).reverse().toString();
  }
  
  public static int[] partionBC(int[][] upper, int [][]down){
	  // find the 2D coordinates that maximizes the sum of two matrices
	  // store 2D coordinates into result
	  int [][] rotated180 = ClockWiseRotate90Degree(ClockWiseRotate90Degree(down));
	  int [] result = new int [2];
	  int max = -1000000;
	  for (int i = 0; i < upper.length; ++i){
		  for (int j = 0; j < upper[0].length; ++j){
			 if (upper[i][j] + rotated180[i][j] > max){
				 max = upper[i][j] + rotated180[i][j];
				 result[0] = i;
				 result[1] = j;
			 }
		  }
	  }
	  return result;
  }
  
  public static NamedObject3D recursive_call(String A, String B, String C){
	  NamedObject3D result;
	  if (A.length() <= 1 || B.length() <= 1 || C.length() <= 1){
		 result = alignThreeSequence(A, B, C);
	  }
	  else{
		  int xmid = A.length()/2;
		  int [][] upper = scoreThreeSeq(A.substring(0, A.length()/2), B, C);
		  int [][] down = scoreThreeSeq(reverseStr(A.substring(A.length()/2)),
				  reverseStr(B), reverseStr(C));
		  int [] coordinates = partionBC(upper, down);
		  int b = coordinates[0]; int c = coordinates[1];
		  NamedObject3D r_upper = recursive_call(A.substring(0, A.length()/2),
				  B.substring(0, b), C.substring(0, c));
  		  NamedObject3D r_down = recursive_call(A.substring(A.length()/2),
				  B.substring(b), C.substring(c));
  		  
  		  result = new NamedObject3D(r_upper.A_alignment + r_down.A_alignment,
  				  r_upper.B_alignment + r_down.B_alignment, 
  				  r_upper.C_alignment + r_down.C_alignment,
  				  r_upper.Score + r_down.Score);
	  }
	  return result;
  }
  public static void main(String[] args) {
    test();
  }
}

