package aco;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Random;

import aco.TSP_ACO;
import aco.Ants.Ant;

import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.CellStyle;
import org.apache.poi.ss.usermodel.Font;
import org.apache.poi.ss.usermodel.IndexedColors;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.xssf.usermodel.XSSFSheet;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;

/**
 * ACO algorithms for the TSP
 * 
 * This code is based on the ACOTSP project of Thomas Stuetzle.
 * It was initially ported from C to Java by Adrian Wilke.
 * 
 * Project website: http://adibaba.github.io/ACOTSPJava/
 * Source code: https://github.com/adibaba/ACOTSPJava/
 */
public class Utilities {
	
	/*
     * ################################################
     * ########## ACO algorithms for the TSP ##########
     * ################################################
     * 
     * Version: 1.0
     * File: utilities.c
     * Author: Thomas Stuetzle
     * Purpose: some additional useful procedures
     * Check: README and gpl.txt
     * Copyright (C) 2002 Thomas Stuetzle
     */

    /***************************************************************************
     * Program's name: acotsp
     * 
     * Ant Colony Optimization algorithms (AS, ACS, EAS, RAS, MMAS, BWAS) for the
     * symmetric TSP
     * 
     * Copyright (C) 2004 Thomas Stuetzle
     * 
     * This program is free software; you can redistribute it and/or modify
     * it under the terms of the GNU General Public License as published by
     * the Free Software Foundation; either version 2 of the License, or
     * (at your option) any later version.
     * 
     * This program is distributed in the hope that it will be useful,
     * but WITHOUT ANY WARRANTY; without even the implied warranty of
     * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
     * GNU General Public License for more details.
     * 
     * You should have received a copy of the GNU General Public License
     * along with this program; if not, write to the Free Software
     * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
     * 
     * email: stuetzle no@spam informatik.tu-darmstadt.de
     * mail address: Universitaet Darmstadt
     * Fachbereich Informatik
     * Hochschulstr. 10
     * D-64283 Darmstadt
     * Germany
     ***************************************************************************/

    private static Random random;

    static int seed;
    static int rowIndex = 1;
    
    //private static String filePath = "output/Experiments Fuzzy cMeans.xlsx";
                                                                                                      
    private static String filePath = "../../../Jurnale/jurnal Swarm Intelligence (Springer)/Experimente/Rulari ACO MinMax_" + TSP_ACO.instanceName + ".xlsx";
    private static String filePath1 = "../../../Jurnale/jurnal Swarm Intelligence (Springer)/Experimente/Fronturi Pareto algoritmi (pe parcurs)/ACO MinMax_voiajor prob/ParetoFront_" + TSP_ACO.instanceName + " (m=" + MTsp.m + ")_amplitude_2200 iter (ACO MinMax_voiajor prob)" + ".xlsx";
    private static String filePath2 = "../../../Jurnale/jurnal Swarm Intelligence (Springer)/Experimente/fisiere intermediare/Rulari ACO MinMax_voiajor probabilist_" + TSP_ACO.instanceName + " (m=" + MTsp.m + ")_total, subtour costs and amplitude.txt";
    
    //for setting the content to be inserted in the Excel file
    private static String tours[];
    private static int nrCities[];
    private static double subtoursCost[];
    private static double totalCost;
    private static double longestSubtour;
    private static double amplitude;
    
    //for verbose output: save at each 5 iteration the best (minimum) total cost of all m subtours 
    private static ArrayList<Double> iterTotalCost;
    
    //for verbose output: save at each 5 iteration the cost of the longest tour and the number of the corresponding iteration
    private static ArrayList<Double> iterLongestCost;
    private static ArrayList<Integer> iterNumber;

    //auxiliary routine for sorting an integer array
    static void swap2(double v[], int v2[], int i, int j) {
    	double tmp1;
		int tmp2;
	
		tmp1 = v[i];
		v[i] = v[j];
		v[j] = tmp1;
		
		tmp2 = v2[i];
		v2[i] = v2[j];
		v2[j] = tmp2;
    }

    //recursive routine (quicksort) for sorting one array; second array does the same sequence of swaps
    static void sort2(double v[], int v2[], int left, int right) {
		    int k, last;
		
			if (left >= right)
			    return;
			swap2(v, v2, left, (left + right) / 2);
			last = left;
			for (k = left + 1; k <= right; k++)
			    if (v[k] < v[left])
			    	swap2(v, v2, ++last, k);
			swap2(v, v2, left, last);
			sort2(v, v2, left, last);
			sort2(v, v2, last + 1, right);
    }
    
    //auxiliary routine for sorting an integer array and a double array
    static void swap2_(Double v[], Integer v2[], int i, int j) {
    	double tmp1;
		int tmp2;
	
		tmp1 = v[i];
		v[i] = v[j];
		v[j] = tmp1;
		
		tmp2 = v2[i];
		v2[i] = v2[j];
		v2[j] = tmp2;
    }

    //recursive routine (quicksort) for sorting one array; second array does the same sequence of swaps
    static void sort2_(Double v[], Integer v2[], int left, int right) {
		    int k, last;
		
			if (left >= right)
			    return;
			swap2_(v, v2, left, (left + right) / 2);
			last = left;
			for (k = left + 1; k <= right; k++)
			    if (v[k] < v[left])
			    	swap2_(v, v2, ++last, k);
			swap2_(v, v2, left, last);
			sort2_(v, v2, left, last);
			sort2_(v, v2, last + 1, right);
    }

    //generate a random number that is uniformly distributed in [0,1]
    static double random01() {
		if (random == null) {
		    random = new Random();
		}
	
		return random.nextDouble();
    }
    
    static void writeExcel(int n, int m, int result) {
    	//the file already exists; we should add a new row as the last one in the Excel file
	    if (new File(filePath).canRead()) {
	    	//System.out.println("File already exists..");
	    	try {
		    	FileInputStream file = new FileInputStream(new File(filePath));
		    	
		    	//Create Workbook instance holding reference to .xlsx file
	            XSSFWorkbook workbook1 = new XSSFWorkbook(file);
	 
	            //Get first/desired sheet from the workbook
	            XSSFSheet sheet1 = workbook1.getSheetAt(2);
		    	int countRows = sheet1.getLastRowNum() + 1;
		    	Row newRow = sheet1.createRow(countRows++);
		    	
		    	int cellnum = 0;
		    	Cell cell = newRow.createCell(cellnum++);
		    	cell.setCellValue(n);
		    	cell = newRow.createCell(cellnum++);
		    	cell.setCellValue(m);
		    	cell = newRow.createCell(cellnum++);
		    	cell.setCellValue(result);
			    
				//Write the workbook in file system
			    FileOutputStream out = new FileOutputStream(new File(filePath));
			    workbook1.write(out);
			    out.close();
			    
			    //System.out.println("Written successfully on disk.");
	    	}
			catch (Exception e) {
			    e.printStackTrace();
			}

	    }
	    else {
	    	//Blank workbook
			XSSFWorkbook workbook2 = new XSSFWorkbook(); 
			
			//Create a blank sheet
			XSSFSheet sheet2 = workbook2.createSheet("Results - 51 cities");
		 
			//Iterate over data and write to sheet
			int rownum = 0, cellnum = 0;
			Row row = sheet2.createRow(rownum++);
			Cell cell = row.createCell(cellnum++);
			cell.setCellValue(n);
	    	cell = row.createCell(cellnum++);
	    	cell.setCellValue(m);
	    	cell = row.createCell(cellnum++);
	    	cell.setCellValue(result);
			
			try {
				//Write the workbook in file system
			    FileOutputStream out = new FileOutputStream(new File(filePath));
			    workbook2.write(out);
			    out.close();
			    
			    //System.out.println("Written successfully on disk.");
			} 
			catch (Exception e) {
			    e.printStackTrace();
			}
			    
	    }
    }
    
    static void writeResultsExcel(int trialNumber, boolean saveIterCosts) {
    	Row r, r1;
    	Cell c;
    	int index1 = 0, index2 = 0, index3 = 0, index4 = 0, index5 = 0;
    	//int index6 = 0;
    	
    	//the file already exists; we should add a new row as the last one in the Excel file
	    if (new File(filePath).canRead()) {
	    	//System.out.println("File already exists..");
	    	try {
		    	FileInputStream file = new FileInputStream(new File(filePath));
		    	
		    	//Create Workbook instance holding reference to .xlsx file
	            XSSFWorkbook workbook1 = new XSSFWorkbook(file);
	            
	            int startIndex = 0, rowIndex = 0;
	            switch (MTsp.m) {
	            	case 2: 
	            		startIndex = 0;
	            		rowIndex = 4;
	            		break;
	            	case 3: 
	            		startIndex = 2;
	            		rowIndex = 5;
	            		break;
	            	case 5: 
	            		startIndex = 4;
	            		rowIndex = 7;
	            		break;
	            	case 7: 
	            		startIndex = 6;
	            		rowIndex = 9;
	            		break;
	            	default:
	            		System.out.println("Unknown value for m");
	            		break;         		
	            }
	            
	            //Get desired sheet from the workbook
	            XSSFSheet sheet1 = workbook1.getSheetAt(startIndex);  //for tours
	            /*XSSFSheet sheet2 = workbook1.getSheetAt(startIndex + 1);  //for number of assigned cities
	            XSSFSheet sheet3 = workbook1.getSheetAt(startIndex + 2);  //for cost of individual subtours
	            XSSFSheet sheet4 = workbook1.getSheetAt(startIndex + 3);  //for total cost of subtours
	            XSSFSheet sheet5 = workbook1.getSheetAt(startIndex + 4);  //for verbose output of total cost at each 5 iteration
	            */
	            XSSFSheet sheet2 = workbook1.getSheetAt(startIndex + 1);  //for verbose output of longest cost at each 5 iteration
	            
	            //define a cell style for bold font
	            CellStyle style = workbook1.createCellStyle();
	            Font font = workbook1.createFont();
	            font.setBoldweight(Font.BOLDWEIGHT_BOLD);
	            style.setFont(font);
	            
	            //define style with bold font and blue color for font
	            CellStyle styleBoldBlue = workbook1.createCellStyle();
	            font = workbook1.createFont();
	            font.setBoldweight(Font.BOLDWEIGHT_BOLD);
	            font.setColor(IndexedColors.BLUE.index);
	            styleBoldBlue.setFont(font);
	            
	            index1 = 18;
	            if (!saveIterCosts) {
	            	//write only once the name of the algorithm that was run
	            	if (trialNumber == 0) {
			            r = sheet1.getRow(index1); 
			            if (r == null) {
			               // First cell in the row, create
			               //System.out.println("Empty row, create new one");
			               r = sheet1.createRow(index1);
			            }
		
			            c = r.getCell(0); 
			            if (c == null) {
			                // New cell
			            	//System.out.println("Empty cell, create new one");
			                c = r.createCell(0);
			            }
			            c.setCellValue("Obtained solutions (values) after running vers. III (ACO MinMax global, voiajor ales probabilist) - vers. veche de alg. ACO");
			            c.setCellStyle(styleBoldBlue);
	            	}
		            
		            //write only once the table header
		            index1 = index1 + 3;
		            r = sheet1.getRow(index1); 
		            if (r == null) {
		               // First cell in the row, create
		               //System.out.println("Empty row, create new one");
		               r = sheet1.createRow(index1);
		            }
	
		            c = r.getCell(0); 
		            if (c == null) {
		                // New cell
		            	//System.out.println("Empty cell, create new one");
		                c = r.createCell(0);
		            }
		            c.setCellValue("Run #");
		            c.setCellStyle(style);
		            
		            c = r.getCell(1); 
		            if (c == null) {
		                // New cell
		            	//System.out.println("Empty cell, create new one");
		                c = r.createCell(1);
		            }
		            c.setCellValue("MinMax (cost of longest subtour)");
		            c.setCellStyle(style);
		            
		            c = r.getCell(2); 
		            if (c == null) {
		                // New cell
		            	//System.out.println("Empty cell, create new one");
		                c = r.createCell(2);
		            }
		            c.setCellValue("Total Cost");
		            c.setCellStyle(style);
		            
		            c = r.getCell(3); 
		            if (c == null) {
		                // New cell
		            	//System.out.println("Empty cell, create new one");
		                c = r.createCell(3);
		            }
		            c.setCellValue("Amplitude");
		            c.setCellStyle(style);
		            
		            //write number of run
		            index1 = 22 + trialNumber;
		            r = sheet1.getRow(index1); 
		            if (r == null) {
		               // First cell in the row, create
		               //System.out.println("Empty row, create new one");
		               r = sheet1.createRow(index1);
		            }
	
		            c = r.getCell(0); 
		            if (c == null) {
		                // New cell
		            	//System.out.println("Empty cell, create new one");
		                c = r.createCell(0);
		            }
		            c.setCellValue(trialNumber + 1);

		            //write MinMax (cost of longest subtour)
		            double longestSubtour = getLongestSubtour();
		            c = r.getCell(1); 
		            if (c == null) {
		                // New cell
		            	//System.out.println("Empty cell, create new one");
		                c = r.createCell(1);
		            }
		            c.setCellValue(longestSubtour); 
		            
		            //write total cost
		            double totalCost = getTotalCost();
		            c = r.getCell(2); 
		            if (c == null) {
		                // New cell
		            	//System.out.println("Empty cell, create new one");
		                c = r.createCell(2);
		            }
		            c.setCellValue(totalCost); 
		            
		            //write amplitude
		            double amplitude = getAmplitude();
		            c = r.getCell(3); 
		            if (c == null) {
		                // New cell
		            	//System.out.println("Empty cell, create new one");
		                c = r.createCell(3);
		            }
		            c.setCellValue(amplitude); 
            	}
	            
	            index5 = 1;
	            if (saveIterCosts) {
	            	//write only once the name of the algorithm that was run
	            	if (trialNumber == 0) {
			            r = sheet2.getRow(index5); 
			            if (r == null) {
			               // First cell in the row, create
			               //System.out.println("Empty row, create new one");
			               r = sheet2.createRow(index5);
			            }
		
			            c = r.getCell(0); 
			            if (c == null) {
			                // New cell
			            	//System.out.println("Empty cell, create new one");
			                c = r.createCell(0);
			            }
			            c.setCellValue("Longest cost of subtour at each 5 iteration after running  vers. III (ACO MinMax global, voiajor ales probabilist) - vers. veche de alg. ACO");
			            c.setCellStyle(styleBoldBlue);
			          
			            int tempIndex = index5 + 3;
			            r = sheet2.getRow(tempIndex);
			            if (r == null) {
				               // First cell in the row, create
				               //System.out.println("Empty row, create new one");
				               r = sheet2.createRow(tempIndex);
				        }
			            ArrayList<Integer> iterNumber = getIterNumber();
			            
			            c = r.getCell(0); 
			            if (c == null) {
			                // New cell
			            	//System.out.println("Empty cell, create new one");
			                c = r.createCell(0);
			            }
			            c.setCellValue("Nr Iter");
			            c.setCellStyle(style); 
			            
			            int indexTemp = 0;
			            for (int j = 0; j < iterNumber.size(); j++) {
			            	indexTemp = tempIndex + 1 + j;
				            r1 = sheet2.getRow(indexTemp); 
				            if (r1 == null) {
				               // First cell in the row, create
				               //System.out.println("Empty row, create new one");
				               r1 = sheet2.createRow(indexTemp);
				            }
			
				            c = r1.getCell(0); 
				            if (c == null) {
				                // New cell
				            	//System.out.println("Empty cell, create new one");
				                c = r1.createCell(0);
				            }
				            c.setCellValue(iterNumber.get(j));     
		            	}            
	            	}
	            	
	            	index5 = index5 + 3;
		            r = sheet2.getRow(index5); 
		            if (r == null) {
		               // First cell in the row, create
		               //System.out.println("Empty row, create new one");
		               r = sheet2.createRow(index5);
		            }
		            
		            //for each trial run save at each 5 iteration the best longest cost of a subtour so far
		            ArrayList<Double> iterLongestCost = getIterLongestCost();
		            int index;
		            
		            //for each run write the table header cell
		            c = r.getCell(trialNumber + 1); 
		            if (c == null) {
		                // New cell
		            	//System.out.println("Empty cell, create new one");
		                c = r.createCell(trialNumber + 1);
		            }
		            c.setCellValue("Run " + (trialNumber + 1));
		            c.setCellStyle(style); 
		            	
		            for (int j = 0; j < iterLongestCost.size(); j++) {
	            		index = index5 + 1 + j;
			            r1 = sheet2.getRow(index); 
			            if (r1 == null) {
			               // First cell in the row, create
			               //System.out.println("Empty row, create new one");
			               r1 = sheet2.createRow(index);
			            }
		
			            c = r1.getCell(trialNumber + 1); 
			            if (c == null) {
			                // New cell
			            	//System.out.println("Empty cell, create new one");
			                c = r1.createCell(trialNumber + 1);
			            }
			            c.setCellValue(iterLongestCost.get(j));     
	            	}
	            }

				//Write the workbook in file system
			    FileOutputStream out = new FileOutputStream(new File(filePath));
			    workbook1.write(out);
			    out.close();
			    
			    int nrOfRun = trialNumber + 1;
			    System.out.println("\nRun #" + nrOfRun + " written successfully on disk.\n");
	    	}
			catch (Exception e) {
			    e.printStackTrace();
			}

	    }
	    else {
	    	//Blank workbook
	    	System.out.println("File " + filePath + " doesn't exists.."); 
			    
	    }
    }
    
    static void writeParetoSet(ArrayList<Ant> bestSoFarPareto, int trial) {
    	Row r;
    	Cell c;
    	int lineNumber = 0;
    	
    	//the file already exists; we should add a new row as the last one in the Excel file
	    if (new File(filePath1).canRead()) {
	    	//System.out.println("File already exists..");
	    	try {
		    	FileInputStream file = new FileInputStream(new File(filePath1));
		    	
		    	//Create Workbook instance holding reference to .xlsx file
	            XSSFWorkbook workbook1 = new XSSFWorkbook(file);
	 
	            //Get first/desired sheet from the workbook
	            XSSFSheet sheet1 = workbook1.getSheetAt(trial);
	            
	            //write table header cells
	            r = sheet1.getRow(lineNumber); 
	            if (r == null) {
	               // First cell in the row, create
	               r = sheet1.createRow(lineNumber);
	            }
	            c = r.getCell(0); 
	            if (c == null) {
	                // New cell
	                c = r.createCell(0);
	            }
	            c.setCellValue("Point #");   
	            c = r.getCell(1); 
	            if (c == null) {
	                // New cell
	                c = r.createCell(1);
	            }
	            c.setCellValue("Total tours length");            
	            c = r.getCell(2); 
	            if (c == null) {
	                // New cell
	                c = r.createCell(2);
	            }
	            c.setCellValue("Amplitude of tours");	            
	            c = r.getCell(3); 
	            if (c == null) {
	                // New cell
	                c = r.createCell(3);
	            }
	            c.setCellValue("List with cost of subtours");	
	            
	            lineNumber++;
	            for (int i = 0; i < bestSoFarPareto.size(); i++) {
	            	r = sheet1.getRow(i + lineNumber); 
		            if (r == null) {
		               // First cell in the row, create
		               r = sheet1.createRow(i + lineNumber);
		            }
		            //write point id
		            c = r.getCell(0); 
		            if (c == null) {
		                // New cell
		                c = r.createCell(0, Cell.CELL_TYPE_NUMERIC);
		            }
		            c.setCellValue(i + 1);
		            //write total cost and amplitude
	            	for (int indexObj = 0; indexObj < 2; indexObj++) {
	            		c = r.getCell(indexObj + 1); 
			            if (c == null) {
			                // New cell
			                c = r.createCell(indexObj + 1, Cell.CELL_TYPE_NUMERIC);
			            }
			            c.setCellValue(bestSoFarPareto.get(i).costObjectives[indexObj]);
	            	}
	            	//write cost of each individual subtour
		            for (int j = 0; j < bestSoFarPareto.get(i).tour_lengths.length; j++) {
			            c = r.getCell(j + 3); 
			            if (c == null) {
			                // New cell
			                c = r.createCell(j + 3);
			            }
			            c.setCellValue(bestSoFarPareto.get(i).tour_lengths[j]);
		            }
	            }
    
				//Write the workbook in file system
			    FileOutputStream out = new FileOutputStream(new File(filePath1));
			    workbook1.write(out);
			    out.close();
			    
			    //System.out.println("\nWritten Pareto front points successfully on disk.\n");
			    int nrOfRun = trial + 1;
			    System.out.println("\nRun #" + nrOfRun + " written Pareto front points successfully on disk.\n");
	    	}
			catch (Exception e) {
			    e.printStackTrace();
			}

	    }
	    else {
	    	System.out.println(" File " + filePath1 + " doesn't exists" );
	    }

    }

    //save in a .txt output file the best solution resulted after a run to be later used when
    //computing the Pareto front
    static void writeParetoSolutions(ArrayList<Ant> bestSoFarPareto) {
    	File f = new File(filePath2);
    	double[] objValues = new double[2];
    	
    	try {
    		BufferedWriter bw = new BufferedWriter(new FileWriter(f, true));
    			
    		for (int i = 0; i < bestSoFarPareto.size(); i++) {
    			if (rowIndex > 0 && i != 0) {
    				bw.newLine();
    			}
    			bw.write(rowIndex + "\t");
    			//get values total cost and amplitude
    			for (int indexObj = 0; indexObj < 2; indexObj++) {
    				objValues[indexObj] = bestSoFarPareto.get(i).costObjectives[indexObj];
    			}
    			bw.write(objValues[0] + "\t");
    			//write cost of each individual subtour
    			for (int j = 0; j < bestSoFarPareto.get(i).tour_lengths.length; j++) {
    				bw.write(bestSoFarPareto.get(i).tour_lengths[j] + "\t");
    			}
    			bw.write(objValues[1] + "\t");
    			
    			rowIndex++;
    		}
    		bw.newLine();
    		bw.close();
        }
    	catch (IOException e) {
    		System.out.println("error writing file");
        }
    }

	public static String[] getTours() {
		return tours;
	}

	public static void setTours(String[] tours) {
		Utilities.tours = tours;
	}

	public static int[] getNrCities() {
		return nrCities;
	}

	public static void setNrCities(int[] nrCities) {
		Utilities.nrCities = nrCities;
	}

	public static double[] getSubtoursCost() {
		return subtoursCost;
	}

	public static void setSubtoursCost(double[] subtoursCost) {
		Utilities.subtoursCost = subtoursCost;
	}

	public static double getTotalCost() {
		return totalCost;
	}

	public static void setTotalCost(double totalCost) {
		Utilities.totalCost = totalCost;
	}

	public static ArrayList<Double> getIterTotalCost() {
		return iterTotalCost;
	}

	public static void setIterTotalCost(ArrayList<Double> iterTotalCost) {
		Utilities.iterTotalCost = iterTotalCost;
	}

	public static ArrayList<Double> getIterLongestCost() {
		return iterLongestCost;
	}
	
	public static ArrayList<Integer> getIterNumber() {
		return iterNumber;
	}

	public static void setIterLongestCost(ArrayList<Double> iterLongestCost) {
		Utilities.iterLongestCost = iterLongestCost;
	}
	
	public static void setIterNumber(ArrayList<Integer> iterNumber) {
		Utilities.iterNumber = iterNumber;
	}

	public static double getLongestSubtour() {
		return longestSubtour;
	}

	public static void setLongestSubtour(double longestSubtour) {
		Utilities.longestSubtour = longestSubtour;
	}

	public static double getAmplitude() {
		return amplitude;
	}

	public static void setAmplitude(double amplitude) {
		Utilities.amplitude = amplitude;
	}
    
	
    
    
}
