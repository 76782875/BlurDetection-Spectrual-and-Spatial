package image_scaling_UI;

import java.awt.Color;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Set;
import java.util.TreeSet;

import javax.imageio.ImageIO;

import Jama.Matrix;
import helper.Complex;
import helper.FFT;

import static jeigen.Shortcuts.*;

import jeigen.DenseMatrix;

public class imProcess {
	BufferedImage origin_img;
	//	int[] RGBs;
	//	double[] gray_matrix;
	//	float[][] rgb_matrix;

	public imProcess() {
		////nothing
	}

	/**
	 * return and "return" some matrices of image features
	 * @param gray_matrix_int 
	 */
	public void bufferedImage2matrix(BufferedImage origin_img, double[] gray_matrix, float[][] rgb_matrix, int[] RGBs, int[] gray_matrix_int) {
		int width = origin_img.getWidth();
		int height = origin_img.getHeight();
		//		this.gray_matrix = new double[width*height];
		//		this.rgb_matrix = new float[width*height][3];
		//		this.RGBs = new int[width*height];
		try {
			int count = 0;
			//			alpha = (origin_img.getRGB(0,0)>>24)&0xff0000;
			for(int i=0; i< width; i++) {
				for(int j=0; j < height; j++) {
					RGBs[count] = origin_img.getRGB(i,j);
					rgb_matrix[count] = getRgb(origin_img.getRGB(i,j));
					gray_matrix[count] = grayScale(origin_img.getRGB(i,j));
					gray_matrix_int[count++] = grayScale(origin_img.getRGB(i,j));
				}
			}
		}
		catch(Exception e) {
			System.out.println("error from matrix2bufferedImage: " + e);
		}
	}

	/**
	 * execute the FFT to image matrix and return the realPart
	 */
	public static double[] real(double[] gray_matrix) {
		Complex[] origin_re_im = new Complex[gray_matrix.length];
		for(int i = 0; i < origin_re_im.length; i++) {
			origin_re_im[i] = new Complex(gray_matrix[i], 0);
		}
		System.out.println("finish constructing complex numbers");

		Complex[] DTF_re_im = FFT.fft(origin_re_im);
		System.out.println("finish fft");

		double[] real = new double[origin_re_im.length];
		for(int i = 0; i < real.length; i++) {
			real[i] = DTF_re_im[i].re();
		}
		return real;
	}

	/**
	 * execute the logorithm to realPart and return log_output
	 */
	public static double[] log(double[] realPart) {
		shiftToPositive(realPart);
		double[] result_log = new double[realPart.length];
		//		imProcess.printdouble(realPart, 256, 256);

		for(int i = 0; i < realPart.length; i ++) {
			result_log[i] = Math.log(realPart[i]);
			//			if(realPart[i] == 0) {
			//				System.err.println(i + ": " + realPart[i]);
			//			}
			//			System.err.println("65534 " + realPart[65534]);
			//			System.out.println(realPart[i] + "->" + result_log[i]);
			//note: if the value is 0, log(0) is infinity, but result_log[i] keeps 0.0;
		}
		//		imProcess.printdouble(result_log, 256, 256);
		return result_log;
	}

	/**
	 * add the min in the array to all the numbers, return the result array
	 */
	private static void shiftToPositive(double[] input) {
		double[] shifted = input.clone();
		Arrays.sort(shifted);
		double min = shifted[0];
		for(int i = 0; i < shifted.length; i++) {
			input[i] -= min;
			//			input[i] += 1;//make sure no log(0)
			if(input[i] == 0) {
				input[i] += 1;
				//				System.out.println("i: " + i + " = " + min);
			}
		}
		//		//debugging
		//		boolean bad = false;
		//		for(int i = 0; i < shifted.length; i++) {
		//			if(input[i] <= 0) {
		//				System.out.println("bad one: " + i +":" + input[i]);
		//				bad = true;
		//			}
		//		}
		//		System.out.println("bad: " + bad);
	}


	/**
	 * perform the average filter and "return" areraged_log
	 * @return 
	 */
	public static double[] A_f(int filter_size, int width, int height, double[] log_output) {
		int kernel_pixels = filter_size*filter_size;
		int[][][]mask = correlationKernel(filter_size);
		int index;
		double[] smoothed = new double[log_output.length];
		int mask_ind_w, mask_ind_h, cur_ind_for_masking;
		double sum = 0;
		for(int i = 0; i < width; i++) {
			for(int j = 0; j < height; j++) {
				index = i*width + j;
				for(int m = 0; m < filter_size; m++) {
					for(int n = 0; n < filter_size; n++) {
						mask_ind_w = i + mask[m][n][0];
						mask_ind_h = j + mask[m][n][1];
						//System.out.println("cur maks_w index: " + mask_ind_w + " cur mask h index: " + mask_ind_h);
						if(!(mask_ind_w < 0 || mask_ind_w > width - 1 || mask_ind_h < 0 || mask_ind_h > height - 1)) {
							cur_ind_for_masking = mask_ind_w*width + mask_ind_h;
							sum += log_output[cur_ind_for_masking];
						}
					}
				}
				//calculate the new_gray
				smoothed[index] = sum*1.0/kernel_pixels;
				sum = 0;
			}
		}
		return smoothed;
	}

	/**
	 * return spectrum_residual = log(realOutput - averaged_log)
	 */
	public static double[] S_f(double[] realOutput, double[] averaged_log) {
		double[] spectrum_residual = new double[averaged_log.length];
		for(int i = 0; i < averaged_log.length; i++) {
			spectrum_residual[i] = realOutput[i] - averaged_log[i];
		}
		return spectrum_residual;
	}

	/**
	 * ifft operation
	 */
	public static double[] B_I(double[] S_f) {
		Complex[] origin_re_im = new Complex[S_f.length];
		for(int i = 0; i < origin_re_im.length; i++) {
			origin_re_im[i] = new Complex(S_f[i], 0);
		}
		System.out.println("finish constructing complex numbers");

		Complex[] DTF_re_im = FFT.ifft(origin_re_im);
		//debugging
		//		for(int i = 0; i < DTF_re_im.length; i++) {
		//			if(i%10 == 0) {
		//				System.out.println(origin_re_im[i] + ": " + DTF_re_im[i].im());
		//			}
		//			if(i%256 == 0) {
		//				System.out.println();
		//			}
		//		}
		System.out.println();
		System.out.println("finish ifft");

		double[] B_I = new double[origin_re_im.length];
		for(int i = 0; i < B_I.length; i++) {
			//debugging
			//			if(i%10 == 0) {
			//				System.out.println(DTF_re_im[i] + ": " + DTF_re_im[i].re());
			//			}
			//			if(i%256 == 0) {
			//				System.out.println();
			//			}
			B_I[i] = DTF_re_im[i].re();
		}
		return B_I;
	}

	/**
	 * logistic regression function, modify the input array literally 
	 */
	public static void B_I_bar(double[] input) {
		for(int i = 0; i < input.length; i++) {
			input[i] = 1/(1+Math.exp(input[i]));
		}
	}

	public static String showCounterImage(int[] RGBs, int[] klabels, int width, int height, SLIC slic, int color) {
		slic.DrawContoursAroundSegments(RGBs, klabels, width, height, color);
		BufferedImage image_new = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);

		String image_name = "contour.jpg";
		int index;
		for(int i = 0; i < width; i++) {
			for(int j = 0; j < height; j++) {
				index = i*width + j;
				image_new.setRGB(i, j, 0xff000000 | RGBs[index]);
			}
		}
		imProcess.saveImage(image_new, image_name);
		return image_name;
	}

	/**
	 * @param cs_mean
	 * @param slic
	 * @param CIElabs
	 * @param RGBs
	 * @param klabels
	 * @param num_of_labels
	 */
	public static void rbg2CIElab( SLIC slic, double[][] CIElabs, int[] RGBs, int[] klabels) {
		double[] lvec = new double[RGBs.length], avec = new double[RGBs.length], bvec = new double[RGBs.length];
		slic.DoRGBtoLABConversion(RGBs, lvec, avec, bvec);
		CIElabs[0] = lvec;
		CIElabs[1] = avec;
		CIElabs[2] = bvec;
		//debug
		//		for(int i = 0; i < lvec.length; i++) {
		//			if(i < 70000 && i > 40000)
		//				System.out.println(i + ": " + lvec[i] + "  "+ avec[i] + " " + bvec[i]);
		//		}
	}
	public static void CIElabregionAverage(float[][] cs_mean, SLIC slic, double[][] CIElabs, int[] klabels, int num_of_labels) {
		//1. mean color of each region
		//get the mean of CIElab
		int[] np = new int[num_of_labels];
		int last_label = 0;//klabel starts with 0.
		//sum up for each labeled area
		for(int i = 0; i < CIElabs[0].length; i++) {
			if(klabels[i] == last_label) {
				//same region
				np[last_label] ++;
				cs_mean[last_label][0] += CIElabs[0][i];
				cs_mean[last_label][1] += CIElabs[1][i];
				cs_mean[last_label][2] += CIElabs[2][i];
				//				if(last_label == 0 || last_label == 34) {
				//					System.out.println(last_label + ": " + CIElabs[0][i]);
				//				}

			}else {
				//				if(last_label == 0 || last_label == 34) {
				//					System.out.println(np[last_label] + " " + last_label + ": " + cs_mean[last_label][0] + " " + cs_mean[last_label][1] + " " + cs_mean[last_label][2]);
				//				}
				//new region
				last_label = klabels[i];
				i--;//for not missing this label and resume the sumup in the next round
			}
		}
		//averaging
		for(int i = 0; i < np.length; i++) {
			cs_mean[i][0] /= np[i];
			cs_mean[i][1] /= np[i];
			cs_mean[i][2] /= np[i];
			//			if(i == 89) {
			//				System.out.println(i + " np: " + np[i]);
			//				System.out.println(cs_mean[i][0] + " " + cs_mean[i][1] + " " + cs_mean[i][2]);
			//			}

			//			if(i == 0)
			//				System.out.println(cs_mean[i][0] + " " + cs_mean[i][1] + " " + cs_mean[i][2]);
		}

		//debugging
		//		boolean b = false;
		//		//verify
		//		for(int j = 0; j < num_of_labels; j++) {
		//			if(cs_mean[j][0] > 255 || cs_mean[j][0] < 0) {
		//				System.out.println(j + "l :" + cs_mean[j][0]);
		//				b = true;
		//			}
		//			if(cs_mean[j][1] > 128 || cs_mean[j][1] < -128) {
		//				System.out.println(j + "a :" + cs_mean[j][1]);
		//				b = true;
		//			}
		//			if(cs_mean[j][2] > 128 || cs_mean[j][2] < -128) {
		//				System.out.println(j + "b :" + cs_mean[j][2]);
		//				b = true;
		//			}
		//		}
		//		System.out.println("b: " + b);
	}


	/**
	 * 
	 * @param num_of_labels
	 * @param cs_neighbor
	 * @param width
	 * @param height
	 * @param klabels
	 */
	public static void getNeighbors(int num_of_labels, HashMap<Integer, Set<Integer>> cs_neighbor, int width, int height, int[] klabels) {
		int index;
		//initialize the cs_neighbor
		for(int i = 0; i < num_of_labels; i++) {
			cs_neighbor.put(i, new TreeSet<Integer>());
		}
		//look to the right and below
		for(int i = 0; i < width-1; i++) {
			for(int j = 0; j < height-1; j++) {
				index = i*width + j;
				cs_neighbor.get(klabels[index]).add(klabels[index + 1]);
				cs_neighbor.get(klabels[index+1]).add(klabels[index]);


				cs_neighbor.get(klabels[index]).add(klabels[index + width]);
				cs_neighbor.get(klabels[index + width]).add(klabels[index]);
			}
		}
		//right edge of image
		index = - 1;
		for(int j = 0; j < height-1; j++) {
			index += width;
			cs_neighbor.get(klabels[index]).add(klabels[index + width]);
		}
		//bottom edge of image
		index = (height-1)*width - 1;
		for(int i = 0; i < width-1; i++) {
			index += 1;
			cs_neighbor.get(klabels[index]).add(klabels[index + 1]);
		}
	}

	public static void colorSimilarity(int num_of_labels, HashMap<Integer, Set<Integer>> cs_neighbor, float[][] cs_mean, double sigma, float[][] cs_ij) {
		for(int i = 0; i < num_of_labels; i ++) {
			for(int j = 0; j < num_of_labels; j++) {
				if(i == j) {
					cs_ij[i][j] = 0;
				}else if(cs_neighbor.get(i).contains(j)) {
					cs_ij[i][j] = (float) Math.exp(-dis_CIElab(cs_mean[i], cs_mean[j])/sigma);
					//					System.out.println(i + " " + j + ": " + Arrays.toString(cs_mean[i]) + " " + Arrays.toString(cs_mean[j]) + " "+ -dis_CIElab(cs_mean[i], cs_mean[j]) + " " + Math.exp(-dis_CIElab(cs_mean[i], cs_mean[j])/sigma_sqr));
					//					System.out.print ("(" + i + "," + j + ": " + cs_ij[i][j] + " dis: " + -dis_CIElab(cs_mean[i], cs_mean[j]) + ") ");
					//					System.out.print (i + ":" + cs_mean[i][0] + " " +cs_mean[i][1] + " " +cs_mean[i][2] + " VS ");
					//					System.out.println(j + ":" + cs_mean[j][0] + " " +cs_mean[j][1] + " " +cs_mean[j][2]);

				}else {
					//not my neighbor
					cs_ij[i][j] = -1;
				}
			}
		}

	}

	public static String drawDebugCIE(double[][] CIE, int width, int height, String filename) {
		BufferedImage image_new = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
		//		CIE2RGB cie2rgb = new CIE2RGB();

		int index, L, a, b, color;
		double rgb[];
		for(int i = 0; i < width; i++) {
			for(int j = 0; j < height; j++) {
				index = i*width + j;
				L = (int) Math.round(CIE[0][index]);
				a = (int) Math.round(CIE[1][index]);
				b = (int) Math.round(CIE[2][index]);
				//				color = cie2rgb.rgb(L, a, b);

				rgb = CIE2RGB.cie2rgb2(L, a, b);
				//				System.out.println(Arrays.toString(rgb));
				color= 0xff000000 | // hardcode alpha
						((((int)rgb[0])<<16)&0xff0000) |
						((((int)rgb[1])<<8)&0xff00) |
						((int)rgb[2]) ;
				image_new.setRGB(i, j, color);
			}
		}

		String image_name = filename + ".png";
		saveImage(image_new, image_name);
		return image_name;
	}

	public static String drawMean(float[][] cs_mean, int width, int height, int[] klabels) {
		BufferedImage image_new = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
		//		CIE2RGB cie2rgb = new CIE2RGB();

		int index, L, a, b, color;
		double[] rgb;
		for(int i = 0; i < width; i++) {
			for(int j = 0; j < height; j++) {
				index = i*width + j;
				L = (int) Math.round(cs_mean[klabels[index]][0]);
				a = (int) Math.round(cs_mean[klabels[index]][1]);
				b = (int) Math.round(cs_mean[klabels[index]][2]);
				rgb = CIE2RGB.cie2rgb2(L, a, b);

				//				System.out.println(klabels[index] + ": " + L +  " " + a + " " + b + "->"+ Arrays.toString(rgb));
				color= 0xff000000 | // hardcode alpha
						((((int)rgb[0])<<16)&0xff0000) |
						((((int)rgb[1])<<8)&0xff00) |
						((int)rgb[2]) ;

				image_new.setRGB(i, j, color);
			}
			//			System.out.println();
		}

		String image_name = "cie mean.jpg";
		saveImage(image_new, image_name);
		return image_name;
	}

	private static double dis_CIElab(float[] cielab1, float[] cielab2) {
		return Math.sqrt(Math.pow(cielab1[0] - cielab2[0], 2) + Math.pow(cielab1[1] - cielab2[1], 2) + Math.pow(cielab1[2] - cielab2[2], 2));
	}
	private static float[] getRgb(int value) {
		return new float[] {(value&0xff0000)>> 16, (value&0xff00)>>8, value&0xff};
	}

	private static int grayScale(int value) {
		int red, green, blue;
		red = (value&0xff0000)>> 16;
			green = (value&0xff00)>>8;
						blue = value&0xff;
						value = (int)((red + green + blue)/3);//gray scale
						return value;
	}

	private static void saveImage(BufferedImage image_new, String file_name) {
		File outputfile = new File(file_name);
		try {
			ImageIO.write(image_new, "png", outputfile);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	//the leading color is always 0 and the tailing is 255;
	//	private int[] getColors(int amount) {
	//		int[] colors = new int[amount];
	//		int interval = 256/(amount-1);
	//		colors[0] = 0;
	//		for(int i = 1; i < amount - 1; i++) {
	//			colors[i] = i * interval;
	//		}
	//		colors[amount-1] = 255; 
	//		return colors;
	//	}
	//
	//	private int mapValue(float value, int[] range, int[] colors) {
	//		for(int i = range.length - 1; i >= 0; i--) {
	//			if(value > range[i]) {
	//				return colors[i+1];
	//			}else if(i == 0) {//falls to the last range
	//				return colors[0];
	//			}
	//		}
	//		return 0;
	//	}


	/**
	 * convolution, not correlation
	 * @param isXdirection
	 * @param grayValue
	 * @param gradients
	 * @param width
	 * @param height
	 * @return
	 */
	public static String getGradient_Xdirection(boolean isXdirection, int[] grayValue, int[] gradients, int width, int height) {
		int size = 3;
		String image_name;
		BufferedImage image_new = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
		int[][][]mask = convolutionKernal3();
//		int[][][]mask = correlationKernel(size);
		//make the filter 3x3
		int[][] filter = new int[size][];
		if(isXdirection == true) {
			filter[0] = new int[] {-1, 0, 1};
			filter[1] = new int[] {-2, 0, 2};
			filter[2] = new int[] {-1, 0, 1};
			image_name = "Gx.jpg";
		}else {
			//y direction
			filter[0] = new int[] {-1, -2, -1};
			filter[1] = new int[] {0, 0, 0};
			filter[2] = new int[] {1, 2, 1};
			image_name = "Gy.jpg";
		}

		int index, /*gray,*/ new_gray = 0;

		for(int i = 0; i < width; i++) {
			for(int j = 0; j < height; j++) {
				index = i*width + j;
				//				gray = grayValue[index];

				int mask_ind_w, mask_ind_h, cur_ind_for_masking = 0;
				for(int m = 0; m < size; m++) {
					for(int n = 0; n < size; n++) {
						mask_ind_w = i + mask[m][n][0];
						mask_ind_h = j + mask[m][n][1];
						cur_ind_for_masking = mask_ind_w*height + mask_ind_h;
						if(!(mask_ind_w < 0 || mask_ind_w > width - 1 || mask_ind_h < 0 || mask_ind_h > height - 1)) {
							//filtering
							new_gray += filter[m][n]*grayValue[cur_ind_for_masking];
						}
					}
				}

				gradients[index] = new_gray;
				//reset
				new_gray = 0;
			}
		}
		//normalization
		int[] gradients_normal = normalization(gradients, width*height);
		//		System.out.println(Arrays.toString(gradients));

		//set the new image
		for(int i = 0; i < gradients_normal.length; i++) {
			image_new.setRGB(i/width, i%width, gradients_normal[i]);
		}
		saveImage(image_new, image_name);
		return image_name;
	}

	public static String meanGradientOfRegion(int[] gX, int[] gY, int[] klabels, int[][] gradient_of_Region, int num_of_labels, int width, int height) {
		//get the mean of CIElab
		int last_label = 0;//klabel starts with 0.
		int[] np = new int[num_of_labels];
		for(int i = 0; i < gX.length; i++) {
			if(klabels[i] == last_label) {
				//same region
				//sum(x^2, y^2) and save int the spot]
				np[last_label] ++;
				//debugging
				if(last_label < 2) {
					//					System.out.println(last_label + ":" + np[last_label] + "-> "+ Math.pow(gX[i], 2) + " " + Math.pow(gY[i], 2));
				}
				gradient_of_Region[0][last_label] += gX[i];
				gradient_of_Region[1][last_label] += gY[i];
				//				System.out.println(sum + " = " + gX[i] + ", " + gY[i]);
			}else {
				//new region
				last_label = klabels[i];
				i--;//for not missing this label and resume the sumup in the next round
			}
		}
		for(int i = 0; i < np.length; i++) {
			//debugging
			//			if(i < 2) {
			//				System.out.print(i + ": " + np[i] + " before: " + gradient_of_Region[0][i] + " after: ");
			//			}
			gradient_of_Region[0][i] = (int)Math.round(gradient_of_Region[0][i] * 1.0 /np[i]);
			gradient_of_Region[1][i] = (int)Math.round(gradient_of_Region[1][i] * 1.0 /np[i]);

			//debugging
			//			if(i < 2) {
			//				System.out.println(gradient_of_Region[0][i] + " " + gradient_of_Region[1][i]);
			//			}
		}
		//		System.out.println("len: " + gradient_of_Region[0].length);
		int[] imageX = normalization(gradient_of_Region[0], gradient_of_Region[0].length);
		int[] imageY = normalization(gradient_of_Region[1], gradient_of_Region[1].length);
		Draws.drawGradientBasedOnAveragedGradient(imageX, klabels, width, height, true);
		return Draws.drawGradientBasedOnAveragedGradient(imageY, klabels, width, height, false);
		//		System.out.println(Arrays.toString(gradient_of_Region));
	}

	public static DenseMatrix covarientMatrix(int Gx, int Gy, int GxMean, int GyMean) {
		//form the gx gy matrix 
		DenseMatrix xy_mat = new DenseMatrix(new double[][]{{Gx},{Gy}});

		//form the mean marices
		DenseMatrix mean_mat = new DenseMatrix(new double[][]{{GxMean},{GyMean}});

		//operation
		DenseMatrix sub_mat = xy_mat.sub(mean_mat);//subtraction
		DenseMatrix sub_tran = sub_mat.t();//transpose
		DenseMatrix coMat = sub_mat.mmul(sub_tran);
		//debug
		//		System.out.println("xy_mat: "+ xy_mat);
		//		System.out.println("mean_mat: "+ mean_mat);
		//		System.out.println("sub     : "+ sub_mat);
		//		System.out.println("trans   : "+ sub_tran);
		//		System.out.println("mul: " + coMat);
		//debugging
		return coMat;
	}

	public static DenseMatrix[] regionCoMat(int[] klabels, int[][] gradient_of_Region, int[] Gx, int[] Gy, int num_of_labels) {
		DenseMatrix[] coMats = new DenseMatrix[num_of_labels];
		for(int i = 0; i < num_of_labels; i++) {
			coMats[i] = zeros(2,2);
		}
		int[] np = new int[num_of_labels];
		int last_label = 0;
		//sumup
		for(int i = 0; i < Gx.length; i++) {
			if(klabels[i] == last_label) {
				np[last_label] ++;
				//debugging
				//				if(i == 1) {
				//					System.out.println("before adding: " + coMats[last_label]);
				//					System.out.println("sub: " + covarientMatrix(Gx[i], Gy[i], gradient_of_Region[0][last_label], gradient_of_Region[1][last_label]));				
				//				}
				coMats[last_label] = coMats[last_label].add(covarientMatrix(Gx[i], Gy[i], gradient_of_Region[0][last_label], gradient_of_Region[1][last_label]));
				//debugging
				//					System.out.println("Gx[i]" + Gx[i]);
				//					System.out.println("Gy[i]" + Gy[i]);
				//
				//					System.out.println("meanX: " + gradient_of_Region[0][last_label]);
				//					System.out.println("meanY: " + gradient_of_Region[1][last_label]);
				//
				//					System.out.println("coMats[last_label]: " + coMats[last_label]);
				//debugging
				//				if(i == 1) {
				//					System.out.println(i + " label:" + last_label + " " + coMats[last_label]);
				//				}
				//				if(last_label == 1) {
				//					System.out.println(np[last_label] + " " + coMats[last_label]);
				//				}
			}else {
				last_label = klabels[i];
				i--;
			}
		}
		//average
		for(int i = 0; i < coMats.length; i++) {
			//debugging
			//			if(i == 1) {
			//				System.out.println("int avearging " + i + coMats[i]);
			//				System.out.println("counter: " + np[i]);
			//			}
			coMats[i] = coMats[i].div(np[i]*1.0);
		}
		return coMats;
	}

	public static int[] normalization(int[] px_new, int size) {
		int new_gray, gray;
		int[] image = new int[size];
		//sort the new pixels
		int[] px_new_sort = px_new.clone();
		Arrays.sort(px_new_sort);
		int min = px_new_sort[0];
		int max = px_new_sort[px_new_sort.length - 1];
		int diff = max - min;
		for(int i = 0; i < px_new.length; i++) {
			//normalized and set the new image
			new_gray = (int) ((px_new[i]-min) * 255.0/diff);
			////TODO: right? no bugs
			px_new[i] = new_gray;
			gray= 0xff000000 | // hardcode alpha
					((((int)new_gray)<<16)&0xff0000) |
					((((int)new_gray)<<8)&0xff00) |
					((int)new_gray) ;
			image[i] = gray;
		}

		return image;
	}

	public static DenseMatrix msqrt(DenseMatrix mat) {
		double A = mat.get(0, 0), B = mat.get(0, 1), C = mat.get(1, 0), D = mat.get(1, 1);
		double tao = A + D;
		double determinant = A*D - B*C;
		double s_positive = Math.sqrt(determinant);

		//		System.out.println("tao: " + tao);
		//		System.out.println("determinate: " + determinant);
		//		System.out.println("s_positive: " + s_positive);

		double t = Math.sqrt(s_positive * 2 + tao);
		DenseMatrix I = eye(2);
		DenseMatrix R = mat.add(I.mul(s_positive)).div(t);

		//		System.out.println("t: " + t);
		//		System.out.println("I.mul(s_positive) " + I.mul(s_positive));
		//		System.out.println("mat.add(I.mul(s_positive)) " + mat.add(I.mul(s_positive)));
		//		System.out.println("mat.add(I.mul(s_positive)).div(t)" + mat.add(I.mul(s_positive)).div(t));
		return R;
	}

	/**TODO: i and j are for debugging**/
	static boolean debug = false;
	public static DenseMatrix inverse(DenseMatrix mat) {
		double[][] arr = new double[][]{{mat.get(0, 0), mat.get(0, 1)}, {mat.get(1, 0), mat.get(1, 1)}};
		Matrix m_jama = new Matrix(arr);

		debug = false;
		if(m_jama.det() == 0) {
			m_jama.set(0, 0, m_jama.get(0, 0)+0.00001);
			m_jama.set(1, 0, m_jama.get(0, 0)+0.00001);

			//			System.out.print(/*i +*/ "\'nei is " /*+ j*/ + ":new sqrt: " + Arrays.deepToString(m_jama.getArray()));
			debug =true;
		}
		//		System.out.println("sqrt: " + Arrays.deepToString(m_jama.getArray()));

		Matrix m_inverse = m_jama.inverse();
		//		System.out.println("inverse " + Arrays.deepToString(m_inverse.getArray()));

		return new DenseMatrix( new double[][]{{m_inverse.get(0, 0), m_inverse.get(0, 1)}, {m_inverse.get(1, 0), m_inverse.get(1, 1)}});
	}

	public static DenseMatrix mulplyC1negHalfC2C1negHalf(DenseMatrix C1, DenseMatrix C2) {
		//		System.out.println("msqurt of C1:" + msqrt(C1));
		DenseMatrix C1_half = inverse(msqrt(C1));
		//		System.out.println("C1_half " + C1_half);
		//		System.out.println("C1_half.mul(C2): " + C1_half.mmul(C2));
		return C1_half.mmul(C2).mmul(C1_half);
	}

	public static double frobeniusForm(DenseMatrix mat, int m, int n) {
		double sum = 0;
		for(int i = 0; i < mat.rows; i++) {
			for(int j = 0; j < mat.cols; j++) {
				sum += Math.pow(mat.get(i, j), 2);
				//				if(m == 29 && (n == 1 || n == 28)) {
				//					System.out.println(i + ", " + j + ": " + mat.get(i, j));
				//				}
			}
		}
		return Math.sqrt(sum);
	}

	//	public static DenseMatrix trickMat(DenseMatrix log) {
	//		double[][] arr = new double[][]{{log.get(0, 0), log.get(0, 1)}, {log.get(1, 0), log.get(1, 1)}};
	//		Matrix m_jama = new Matrix(arr);
	//		if(m_jama.det() == 0) {
	//			m_jama.set(0, 0, m_jama.get(0, 0)+10);
	//		}
	////		DenseMatrix
	////		System.out.println("in trick: " + );
	//		return new DenseMatrix(m_jama.getArray());
	//	}

	public static double dAI(DenseMatrix C1, DenseMatrix C2, int i, int j) {
		DenseMatrix mul = mulplyC1negHalfC2C1negHalf(C1, C2);
		DenseMatrix log = mul.mlog();
		if(Double.isNaN(log.get(0, 0))) {
			log = new DenseMatrix("0 0;0 0");
		}
		//		if(i == 29 && (j == 1)) {
		//			System.out.println("mulllllllllllllllllllllllllllll " + mul);
		//			System.out.println("dAI log : " + log);
		//			System.out.println("dAI frobeniusForm: " + frobeniusForm(log, i, j));
		//		}
		//		if(debug == true) {
		//			System.out.println("debug is truuuuuuuuuuuuuuuuuueeeeeeeeeeeee");
		//			System.out.println(frobeniusForm(log, i, j));
		//		}
		return frobeniusForm(log, i, j);
	}

	public static void gradientSimilarity(int num_of_labels, HashMap<Integer, Set<Integer>> cs_neighbor, DenseMatrix[] coMatrices, double sigma, float[][] gs_ij) {
		for(int i = 0; i < num_of_labels; i ++) {
			for(int j = 0; j < num_of_labels; j++) {
				if(i == j) {
					gs_ij[i][j] = 0;
				}else if(cs_neighbor.get(i).contains(j)) {
					//					System.out.println("i: " + i + " = " +  coMatrices[i]);
					//					System.out.println("j: " + j + " = " +  coMatrices[j]);
					//					if(i == 29 && (j == 1 || j == 28)) {
					//						System.out.println("i: " + coMatrices[i]);
					//						System.out.println(j + " " + coMatrices[j]);
					//					}
					gs_ij[i][j] = (float) Math.exp(-dAI(coMatrices[i], coMatrices[j], i, j)/sigma);
				}else {
					//not my neighbor
					gs_ij[i][j] = -1;
				}
			}
		}
	}

	public static void similarites(float[][] cs_ij, float[][] gs_ij, float[][] s_ij, int num_of_labels, double beta) {
		double one_minus_beta = 1 - beta;
		for(int i = 0; i < num_of_labels; i ++) {
			for(int j = 0; j < num_of_labels; j++) {
				s_ij[i][j] = (float) (beta * cs_ij[i][j] + one_minus_beta * gs_ij[i][j]);
								if(i == 29 && (j == 1 || j == 28))
									System.out.println(i + ", " + j + "cs: " + cs_ij[i][j] + " gs: " + gs_ij[i][j] +" sim: " + s_ij[i][j]);
			}
			//			System.out.println();
		}
	}

	/**
	 * coh_i = 1/(1+max(s_ij)
	 * @param s_ij
	 * @param num_of_labels
	 * @param coh_i
	 */
	public static void coherence(float[][] s_ij, int num_of_labels, float[] coh_i) {
		float[] temp;
		int last_ind = num_of_labels - 1;
		for(int i = 0; i < num_of_labels; i ++) {
			temp = s_ij[i].clone();
			Arrays.sort(temp);
			coh_i[i] = 1/(1+temp[last_ind]);
			//debugging
			//			if(i == 29) {
			//				System.out.println(Arrays.toString(temp));
			//				System.out.println(temp);
			//				for(int j = 0; j < s_ij[i].length; j++) {
			//					if(Double.isNaN(s_ij[i][j]))
			//						System.out.println("in coherence: " + i + ": " + j + ": " + s_ij[i][j] + " coh:" + coh_i[i]);
			//				}
			//			}
		}
	}

	public static void initializeCoarseBlurMap(double[] B_I, double[] b_i_region, HashMap<Integer, Set<Integer>> cs_neighbor, int[] klabels) {
		int[] np = new int[b_i_region.length];
		int last_label = 0;//klabel starts with 0.
		//sum up for each labeled area
		for(int i = 0; i < B_I.length; i++) {
			if(klabels[i] == last_label) {
				//same region
				np[last_label] ++;
				b_i_region[last_label] += B_I[i];
				//				if(last_label == 0 || last_label == 34) {
				//					System.out.println(last_label + ": " + B_I[i] + " " + b_i_region[last_label]);
				//				}

			}else {
				//new region
				last_label = klabels[i];
				i--;//for not missing this label and resume the sumup in the next round
			}
		}
		//averaging
		for(int i = 0; i < np.length; i++) {
			b_i_region[i] /= np[i];
			//			if(last_label == 0 || last_label == 34) {
			//				System.out.println("average: " + last_label + ": " + b_i_region[last_label]);
			//			}
		}

		//debugging
		//		boolean b = false;
		//		//verify
		//		for(int j = 0; j < num_of_labels; j++) {
		//			if(cs_mean[j][0] > 255 || cs_mean[j][0] < 0) {
		//				System.out.println(j + "l :" + cs_mean[j][0]);
		//				b = true;
		//			}
		//			if(cs_mean[j][1] > 128 || cs_mean[j][1] < -128) {
		//				System.out.println(j + "a :" + cs_mean[j][1]);
		//				b = true;
		//			}
		//			if(cs_mean[j][2] > 128 || cs_mean[j][2] < -128) {
		//				System.out.println(j + "b :" + cs_mean[j][2]);
		//				b = true;
		//			}
		//		}
		//		System.out.println("b: " + b);
	}

	/**
	 * determine a pixel to be drew or not
	 * @param klabels
	 * @param determineRegion
	 * @param width
	 * @param height
	 * @return
	 */
	public static int[] getDetermine(int[] klabels, int[] determineRegion, int width, int height) {
		int[] d = new int[width*height];
		int index;
		for(int i = 0; i < width; i++) {
			for(int j = 0; j < height; j++) {
				index = i*width + j;
				d[index] = determineRegion[klabels[index]] == 1? 1 :0;
			}
		}
		return d;
	}

	/**
	 * determine a region to be drew or not
	 * @param rates
	 * @param factor
	 * @return
	 */
	public static int[] determinationSet(double[] rates, double factor) {
		//get the average and * factor as in the other paper
		//sum
		//		double sum = 0;
		//		for(int i = 0; i < rates.length; i++) {
		//			sum += rates[i];
		//			//			System.out.println("sum: " + sum);
		//		}
		//		//Average
		//		double ave = sum * 1.0/rates.length;
		double ave = imProcess.getAve(rates);
		System.out.println("ave: " + ave);
		double threshold = ave * factor;

		int[] d = new int[rates.length];
		for(int i = 0; i < rates.length; i++) {
			d[i] = rates[i] >= threshold? 1:0;
		}
		return d;
	}

	

	static double getAve(double rate[]) {
		double sum = 0;
		for(int i = 0; i < rate.length; i++) {
			sum += rate[i];
		}
		return sum*1.0/rate.length;
	}

	private static double[] getDis(double rate[], double ave) {
		double[] dis = new double[rate.length];
		for(int i = 0; i < rate.length; i++) {
			dis[i] = Math.pow((rate[i] - ave)*1.0,2)/ rate.length;
		}
		return dis;
	}

	/*******************************************/
	/*******************************************/

	private static double sum_s_ij_b_t(int i, float[][]s_ij, double[] B_I) {
		double sum = 0;
		for(int j = 0; j < s_ij.length; j++) {
			if(s_ij[i][j] != -1) {
				sum = s_ij[i][j] * B_I[j];
			}
		}
		return sum;
	}

	public static double[] iterationUpdate(float[] coh_i, double[] B_I, float[][]s_ij, int step) {
		double[] b_i = B_I.clone();
		for(int s = 0; s < step ; s++) {
			for(int i = 0; i < b_i.length; i ++) {
				b_i[i] = coh_i[i] * b_i[i] + (1 - coh_i[i]) * sum_s_ij_b_t(i, s_ij, b_i);
			}
		}
		return b_i;
	}

	public double[] smoothingFiltering(int width, int height, int kernel, double[] grayValue) {
		int kernel_pixels = kernel*kernel;
		int[][][]mask = correlationKernel(kernel);
		int index;
		double[] smoothed = new double[grayValue.length];
		for(int i = 0; i < width; i++) {
			for(int j = 0; j < height; j++) {
				index = i*width + j;
				int mask_ind_w, mask_ind_h, cur_ind_for_masking, sum = 0;
				for(int m = 0; m < kernel; m++) {
					for(int n = 0; n < kernel; n++) {
						mask_ind_w = i + mask[m][n][0];
						mask_ind_h = j + mask[m][n][1];
						//System.out.println("cur maks_w index: " + mask_ind_w + " cur mask h index: " + mask_ind_h);
						if(!(mask_ind_w < 0 || mask_ind_w > width - 1 || mask_ind_h < 0 || mask_ind_h > height - 1)) {
							cur_ind_for_masking = mask_ind_w*height + mask_ind_h;
							sum += grayValue[cur_ind_for_masking];
						}
					}
				}
				//calculate the new_gray
				smoothed[index] = (int)Math.round(sum*1.0/kernel_pixels);
				//				System.out.print(grayValue[index] + "(" + smoothed[index] + ")  ");
				sum = 0;
			}
			//			System.out.println();
		}
		return smoothed;
	}

	private static int[][][] convolutionKernal3(){
		int kernel = 3;
		int mask[][][] = correlationKernel(kernel);
		//flip relative to the center
		for(int i = 0; i < 2; i++) {
			int tmp = mask[2][2][i];
			mask[2][2][i] = mask[0][0][i];
			mask[0][0][i] = tmp;

			tmp = mask[0][2][i];
			mask[0][2][i] = mask[2][0][i];
			mask[2][0][i] = tmp;
			
			tmp = mask[1][2][i];
			mask[1][2][i] = mask[1][0][i];
			mask[1][0][i] = tmp;
			
			tmp = mask[0][1][i];
			mask[0][1][i] = mask[2][1][i];
			mask[2][1][i] = tmp;
		}
//		System.out.println();
//		for(int i = 0; i < kernel; i++) {
//			for(int j = 0; j < kernel; j++) {
//				System.out.print(Arrays.toString(mask[i][j]));
//				System.out.print(" ");
//			}
//			System.out.println();
//		}
		return mask;
	}


	private static int[][][] correlationKernel(int kernel){
		int mask[][][] = new int[kernel][kernel][2];
		int tmp_index;

		int half_kernel_pixels = kernel*kernel/2;
		for(int i = 0; i < kernel; i++) {
			for(int j = 0; j < kernel; j++) {
				tmp_index = i*kernel + j;
				int[] pair = new int[2];
				if(tmp_index != half_kernel_pixels) {
					pair[0] = i - kernel/2;
					pair[1] = j - kernel/2;
					mask[i][j] = pair;
				}else{
					pair[0] = 0;
					pair[1] = 0;
					mask[i][j] = pair;
				}
				//				System.out.print(Arrays.toString(pair));
				//				System.out.print(" ");
			}
			//			System.out.println();
		}
		return mask;
	}

	//	private static double[] createGrayValueMatric(int[] gray_matrix, int height, int width){
	//		double[] grayValue = new double[height*width];
	//		int index, gray, red, green, blue;
	//		for(int i = 0; i < width; i++) {
	//			for(int j = 0; j < height; j++) {
	//				index = i*width + j;
	//				gray = (int)gray_matrix[index];
	//				red = (gray&0xff0000)>> 16;
	//			green = (gray&0xff00)>>8;
	//			blue = gray&0xff;
	//			gray = (int)((red + green + blue)/3);//gray scale
	//			grayValue[index] = gray;
	//			//				System.out.print(red + "," + gray + " ");
	//			}
	//			//			System.out.println();
	//		}
	//		return grayValue;
	//	}


	public double[] normalization_image(double[] input) {
		double[] image = new double[input.length];
		//sort the new pixels
		double[] input_sort = input.clone();
		Arrays.sort(input_sort);
		double min = input_sort[0];
		double max = input_sort[input_sort.length - 1];
		double diff = max - min;
		for(int i = 0; i < input.length; i++) {
			//normalized and set the new image
			image[i] = (int) ((input[i]-min) * 255.0/diff);
			//			System.out.println("in normal: " + input[i] + " " + image[i]);
		}
		return image;
	}

	public static void printdouble(double[] to_print, int width, int height) {
		//		for(int i = 0; i < width; i += 11) {
		//			for(int j =0; j < height; j += 11) {
		//				System.out.print(to_print[i*width + j] + " ");
		//			}
		//			System.out.println();
		//		}
		int index;
		for(int i = 0; i < width; i ++) {
			for(int j =0; j < height; j ++) {
				index = i*width + j;
				if(Double.isNaN(to_print[index])) {
					System.err.print("NAN: " + index + " " + to_print[index] + " ");
				}else if(Double.isInfinite(to_print[index])) {
					System.err.print("Infinite: " + index + " " + to_print[index] + " ");
				}
			}
		}
	}

}
