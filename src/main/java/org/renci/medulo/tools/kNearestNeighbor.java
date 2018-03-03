package org.renci.medulo.tools;


import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.SortedSet;
import java.util.TreeSet;
import org.apache.commons.math3.linear.*;
import org.apache.commons.math3.stat.StatUtils;


	/**
	 * Created by Leon Aburime on 7/15/2014.
	 * Built using IntelliJ
	 */
	public class kNearestNeighbor {

	    //Data from http://www.saedsayad.com/k_nearest_neighbors.htm
	    //{Age, Loan Amount}
	    public double[][] input = {
	            {25, 40000},
	            {35, 60000},
	            {45, 80000},
	            {20, 20000},
	            {35, 120000},
	            {52, 18000},
	            {23, 95000},
	            {40, 62000},
	            {60, 100000},
	            {48, 220000},
	            {33, 150000}
	    };

	    //If Defaulted
	    public char[] output = {
	            'N',
	            'N',
	            'N',
	            'N',
	            'N',
	            'N',
	            'Y',
	            'Y',
	            'Y',
	            'Y',
	            'Y'
	    };

	    //Value we are trying to predict whether default will happen
	    public double[] predict = {48, 142000};
	    public double[] predict_standardize = {};
	    public int k = 3;//Num of nearest Neighbors

	    public static void main(String[] args) {
	        kNearestNeighbor knn = new kNearestNeighbor();

	        //Lets turn our imput intoa matrix
	        RealMatrix m = MatrixUtils.createRealMatrix(knn.input);
	        knn.predict_standardize = standardizeArray( m, knn.predict );//Standardize prediction values

	        //Standardize my matrix values
	        m = standardizeMatrix(m);

	        List<Character> answers = euclideanDistance(m, knn.predict_standardize, knn.output);

	        char prediction = getMostOccurringElement(answers, knn.k);
	        System.out.println("Predicted value of " + Arrays.toString( knn.predict ) + " is " + prediction);

	    }

	    //Slice the list from 0 to k and return the character that occurs the most
	    public static Character getMostOccurringElement( List <Character>answers, int k){
	        int max = 0;
	        Character c = ' ';
	        List <Character>sublist = answers.subList(0,k);

	        for (int i=0; i<sublist.size(); i++){
	            int occurrences = Collections.frequency(answers.subList(0,k), sublist.get(i));

	            //If counted occurrences are greater than what I have update to current values
	            if(occurrences > max) {
	                c = (Character) sublist.get(i);
	                max = occurrences;
	            }
	        }

	        System.out.println("Sorted Array by distance : " + Arrays.toString( sublist.toArray()));
	        System.out.println(c + " occurs most frequently with value " + max +".");

	        return c;
	    }

	    //Perform Euclidean distance formula to find out the distance
	    //between our prediction value and each row in the matrix
	    public static List <Character>euclideanDistance( RealMatrix m, double[] y,
	                                            char[] output ){

	        Map  <Double,Character>map = new HashMap <Double,Character>();
	        //Lets turn out 'y' value or label into vector for easier math operations
	        RealVector Y = MatrixUtils.createRealVector( y);

	        for (int i=0; i<m.getRowDimension(); i++){
	            RealVector vec = m.getRowVector(i);
	            RealVector sub = vec.subtract( Y );

	            //Take square root of sum of square values that were subtracted a line above
	            double distance = Math.sqrt(StatUtils.sumSq(sub.toArray()));
	            //Use the distance to each data point(or row) as key with the 'default' option as value
	            map.put( distance  , output[i] );
	        }

	        //Now lets sort the map's keys into a set
	        SortedSet<Double> keys = new TreeSet<Double>(map.keySet());
	        List<Character> neighbors = new ArrayList<Character>();

	        //For each key add the values in that order into the list
	        for (Double key : keys) {
	            neighbors.add((Character) map.get(key));
	        }


	        return neighbors;
	    }

	    //Standardize our matrix i.e. ..
	    //(each_value - min_val_in_column)/(max_val_in_column - min_val_in_column
	    public static RealMatrix standardizeMatrix(RealMatrix m){
	        for (int i=0; i<m.getColumnDimension(); i++){
	            //Get each column as Vector
	            RealVector vec = m.getColumnVector(i);

	            //Lets get the max and the min of the column
	            double max = vec.getMaxValue();
	            double min = vec.getMinValue();

	            vec.mapSubtractToSelf(min);
	            vec.mapDivideToSelf(max - min);

	            //Now lets reset the value back into the matrix
	            m.setColumnVector(i, vec);
	        }
	        return m;
	    }

	    //Standardize our vector by the matrix i.e. ..
	    //(each_array_value - min_val_in_column)/(max_val_in_column - min_val_in_column
	    public static double[] standardizeArray(RealMatrix m, double[] arr) {
	        //Create new array to store the our standardized values in
	        double[] new_arr = new double[arr.length];

	        //Iterate each of the columns of the matrix to find max and min.
	        //Lets use those for calculations on each individual entry of the vector
	        for (int i = 0; i < m.getColumnDimension(); i++) {
	            //Get each column as Vector
	            RealVector vec = m.getColumnVector(i);

	            //Lets get the max and the min of the column
	            double max = vec.getMaxValue();
	            double min = vec.getMinValue();

	            //Store each value in the vector
	            new_arr[i] = (arr[i] - min) / (max - min);
	        }

	        return new_arr;
	    }


	}


