package image_scaling_UI;

import java.awt.Color;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;

import javax.imageio.ImageIO;

public class Draws {
	private static double[] normalization_general(double[] input, double bottom, double up) {
		double[] image = new double[input.length];
		//sort the new pixels
		double[] input_sort = input.clone();
		Arrays.sort(input_sort);
		double min = input_sort[0];
		double max = input_sort[input_sort.length - 1];
		double diff = max - min;
		double range = up - bottom;
		if(range <= 0) {
			System.err.println("normalization range is <= 0");
		}
		System.out.println("min, max, diff - " + min + " " + max + " " + diff);
		for(int i = 0; i < input.length; i++) {
			//normalized and set the new image
			image[i] = (int) ((input[i]-min) * range/diff + bottom);
			//			System.out.println("in normal: " + input[i] + " " + image[i]);
		}
		return image;
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
	
	public static String drawDebugGray(double[] grayscale, int width, int height) {
		BufferedImage image_new = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
		int index, gray;
		for(int i = 0; i < width; i++) {
			for(int j = 0; j < height; j++) {
				index = i*width + j;
				gray = 0xff000000 | // hardcode alpha
						((((int)grayscale[index])<<16)&0xff0000) |
						((((int)grayscale[index])<<8)&0xff00) |
						((int)grayscale[index]) ;
				image_new.setRGB(i, j, gray);
			}
		}
		String image_name = "debug.png";
		saveImage(image_new, image_name);
		return image_name;
	}
	
	public static String drawGradientBasedOnAveragedGradient(int[] gradients, int[] klabels, int width, int height, boolean isX) {
		String image_name = isX==true ?"AveragedGradientX.png" : "AveragedGradientY.png";
		BufferedImage image_new = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
		//set the new image
		int index, gray;
		for(int i = 0; i < width; i++) {
			for(int j = 0; j < height; j++) {
				index = i*width + j;
				gray = gradients[klabels[index]];
				//				System.out.println(gray);
				image_new.setRGB(i, j, gray);
			}
		}
		saveImage(image_new, image_name);
		return image_name;
	}
	
	public static String drawImageWithRegressionBase(int[] determine, Color bgColor, Color fgColor, int width, int height) {
		BufferedImage image_new = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);

		int index;
		for(int i = 0; i < width; i++) {
			for(int j = 0; j < height; j++) {
				index = i*width + j;
				if(determine[index] == 1) {
					image_new.setRGB(i, j, bgColor.getRGB());
				}else {
					image_new.setRGB(i, j, fgColor.getRGB());
				}
			}
		}
		String image_name = "twoColorBase.jpg";
		saveImage(image_new, image_name);
		return image_name;
	}

	public static Color[] intervalColors(float angleFrom, float angleTo, int n) {
		float angleRange = angleTo - angleFrom;
		float stepAngle = angleRange / n;
		
		Color[] colors = new Color[n];
		for (int i = 0; i < n; i++) {
			float angle = angleFrom + i*stepAngle;
			colors[i] = Color.getHSBColor(angle, 1, 1);        
		}
		return colors;
	}
	
	public static String drawImageWithRegressionOnRegionHSB(double[]rate, int[] klabels, int width, int height, String image_name) {
		BufferedImage image_new = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
		//red = 0, magenta: 300, blue: 240, cyan = 180, green, 120, yellow: 60 || n=10
		int csteps = 200;
		Color[] colors = Draws.intervalColors(300, 240, csteps);//red = 0, magenta: 300, blue: 240, cyan = 180, green, 120, yellow: 60 || n=10
		int[] cindex = colorPickIndex(rate, imProcess.getAve(rate), 0, 100);
		//map the rate to color
		int index_k = 0;//pixel index
		int index_c = 0;//color index of the cluster
		for(int i = 0; i < width; i++) {
			for(int j = 0; j < height; j++) {
				index_c = cindex[klabels[index_k++]];//klabels[index_k++]: cluster index
				image_new.setRGB(i, j, colors[index_c].getRGB());
			}
		}
		saveImage(image_new, image_name);
		return image_name;
	}
	
	private static int[] colorPickIndex(double rate[], double ave, int normal_up, int normal_bottom) {
		int[] indexes = new int[rate.length];
		double[] normalized = normalization_general(rate, normal_up, normal_bottom);
		for(int i = 0; i < normalized.length; i++) {
			indexes[i] = (int) Math.floor(normalized[i]);
		}
		return indexes;
	}
	
	public static String drawImageOnRegionHSB2(double[]rate, int[] klabels, int width, int height) {
		BufferedImage image_new = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
		Color[] colors = DisplyHSB.intervalHSBs(rate);
		for(int i = 0; i < width; i++) {
			for(int j = 0; j < height; j++) {
				image_new.setRGB(i, j, colors[klabels[i*width + j]].getRGB());
			}
		}
		String image_name = "twoColorFinal.jpg";
		saveImage(image_new, image_name);
		return image_name;
	}

	/**
	 * draw an image based on the determination on a region
	 * @param determine
	 * @param bgColor: bg color
	 * @param fgColor: fg color
	 * @return
	 */
	public static String drawImageWithRegressionOnRegion(int[] determine, int[] klabels, Color bgColor, Color fgColor, int width, int height) {
		BufferedImage image_new = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);

		int index;
		for(int i = 0; i < width; i++) {
			for(int j = 0; j < height; j++) {
				index = i*width + j;
				if(determine[klabels[index]] == 1) {
					image_new.setRGB(i, j, fgColor.getRGB());
				}else {
					image_new.setRGB(i, j, bgColor.getRGB());
				}
			}
		}
		String image_name = "twoColorFinal.jpg";
		saveImage(image_new, image_name);
		return image_name;
	}
	
	private static double[] rateFromRegionToPixel(int[] klabel, double[] rate, int width, int height) {
		double[] pixel_rate = new double[width*height];
		for(int i = 0; i < pixel_rate.length; i++) {
			pixel_rate[i] = rate[klabel[i]];
		}
		return pixel_rate;
	}
	
	public static String drawImageWithRateForPixels(double[] rate, String image_name, int width, int height) {
		BufferedImage image_new = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
		int csteps = 110;
		Color[] colors = Draws.intervalColors(0, 120, csteps);//red = 0, magenta: 300, blue: 240, cyan = 180, green, 120, yellow: 60 || n=10
		int[] cindex = colorPickIndex(rate, imProcess.getAve(rate), 0, 100);
		System.out.println("size: " + cindex.length + " width: " + width + " height: " + height);
		for(int i = 0; i < width; i++) {
			for(int j = 0; j < height; j++) {
				image_new.setRGB(i, j, colors[cindex[i*width+j]].getRGB());
			}
		}
		saveImage(image_new, image_name);
		return image_name;
	}
	
	public static double[] floatssToDoubles(float[][] floats, int width, int height) {
		double[] doubles = new double[width*height];
		for(int i = 0; i < width; i++) {
			for(int j = 0; j < height; j++) {
				doubles[i*width+j] = floats[i][j];
			}
		}
		return doubles;
	}
	
	public static void main(String[] arg) {
		int[] klabel = {1, 1, 0};
		double[] rate = {0.1, 0.4};
		int w = 3, h = 1;
		System.out.println(Arrays.toString(rateFromRegionToPixel(klabel, rate, w, h)));
		
		double i = 1.245;
		System.out.println((float) i);
	}
}
