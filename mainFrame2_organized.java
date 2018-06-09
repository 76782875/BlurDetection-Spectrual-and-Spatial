package image_scaling_UI;
// Using layouts
import java.awt.BorderLayout;
import java.awt.Checkbox;
import java.awt.CheckboxGroup;
import java.awt.Color;
import java.awt.Container;
// Using AWT event classes and listener interfaces
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Set;

import javax.imageio.ImageIO;
// Using Swing components and containers
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSlider;
import javax.swing.JTextField;
import javax.swing.SwingUtilities;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

//DFT
import jeigen.DenseMatrix;


// A Swing GUI application inherits the top-level container javax.swing.JFrame
public class mainFrame2_organized extends JFrame {
	public static boolean RIGHT_TO_LEFT = false;
	//	String original = "/Users/xlei/Documents/off-jpl/digital_image_process/lena_gray.gif";
	//	String original = "src/image_scaling_UI/images/flower256.png";
	String original = "src/image_scaling_UI/images/balls.gif";
	//		String original = "src/image_scaling_UI/images/biking.gif";


	String processed = original, last_processed = original, algorithm;
	String colorSim = original, gradientSim = original, baseNonRegion = original, baseRegion = original;
	BufferedImage origin_bufferedImage = null;
	BufferedImage last_bufferedImage = null;
	BufferedImage processed_bufferedImage = null;
	BufferedImage operating_bufferedImage = null;
	BufferedImage colorSim_bufferedImage = null;
	BufferedImage gradientSim_bufferedImage = null;
	BufferedImage baseNonRegionImage = null;
	BufferedImage baseRegionImage = null;

	String method;

	CheckboxGroup cbg;
	Checkbox origin_cb;
	Checkbox process_cb;
	Checkbox selected_cb;
	FileChooser fc;

	JSlider stepSlider = new JSlider(0, 200, 5);
	JSlider kLabelSlider = new JSlider(0, 1500, 5);
	JTextField jTextStep = new JTextField("0", 2);
	JTextField jTextKsize = new JTextField("0", 3);

	JTextField jTextsigmaColorSim = new JTextField("0.1", 2);
	JTextField jTextsigmaGradientSim = new JTextField("0.1", 2);
	JTextField jTextBeta = new JTextField("0.4", 2);
	JButton OKbutton;

	JLabel label_origin;
	JLabel label_processed;
	JLabel label_colorSim;
	JLabel label_gradientSim;
	JLabel label_baseNonRegion;
	JLabel label_baseRegion;

	ImageIcon origin_ii;
	ImageIcon process_ii;
	ImageIcon colorSim_ii;
	ImageIcon gradientSim_ii;
	ImageIcon baseNonRegion_ii;
	ImageIcon baseRegion_ii;

	// Constructor to setup the GUI components and event handlers
	public mainFrame2_organized() throws IOException {
		// Retrieve the content-pane of the top-level container JFrame
		// All operations done on the content-pane
		Container cp = getContentPane();
		if (!(cp.getLayout() instanceof BorderLayout)) {
			cp.add(new JLabel("Container doesn't use BorderLayout!"));
			return;
		}
		if (RIGHT_TO_LEFT) {
			cp.setComponentOrientation(
					java.awt.ComponentOrientation.RIGHT_TO_LEFT);
		}
		JPanel pane_top = new JPanel();
		fc = new FileChooser();
		pane_top.add(fc, BorderLayout.LINE_START);

		JPanel checkBoxes = new JPanel();
		createListImageSource();
		checkBoxes.add(selected_cb);
		checkBoxes.add(origin_cb);
		checkBoxes.add(process_cb);
		pane_top.add(checkBoxes, BorderLayout.LINE_END);
		cp.add(pane_top, BorderLayout.PAGE_START);

		createOriginLabel();
		JScrollPane scrollPane = new JScrollPane(label_origin, 
				JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
				JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
		cp.add(scrollPane, BorderLayout.LINE_START);

		JPanel pane_processing = new JPanel();
		createProcessedLabel();
		scrollPane = new JScrollPane(label_processed, 
				JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
				JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
		pane_processing.add(scrollPane);

		createBaseNonRegionLabel();
		scrollPane = new JScrollPane(label_baseNonRegion, 
				JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
				JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
		pane_processing.add(scrollPane);

		createBaseRegionLabel();
		scrollPane = new JScrollPane(label_baseRegion, 
				JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
				JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
		pane_processing.add(scrollPane);

		createColorSimLabel();
		scrollPane = new JScrollPane(label_colorSim, 
				JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
				JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
		pane_processing.add(scrollPane);

		createGradientSimLabel();
		scrollPane = new JScrollPane(label_gradientSim, 
				JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
				JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
		pane_processing.add(scrollPane);


		cp.add(pane_processing);

		pane_top = new JPanel();
		createListAlgorithems();

		pane_top.add(new JLabel("step: 0~200"));
		pane_top.add(stepSlider);
		pane_top.add(jTextStep);

		pane_top.add(new JLabel("klabel: 0~1000"));
		pane_top.add(kLabelSlider);
		pane_top.add(jTextKsize);

		//add lisnter
		kLabelSlider.addChangeListener(createListner(jTextKsize));
		stepSlider.addChangeListener(createListner(jTextStep));

		//textfield
		pane_top.add(new JLabel("Color sigma"));
		pane_top.add(jTextsigmaColorSim);
		pane_top.add(new JLabel("Gradient sigma"));
		pane_top.add(jTextsigmaGradientSim);
		pane_top.add(new JLabel("Beta"));
		pane_top.add(jTextBeta);

		createOKbutton();
		pane_top.add(OKbutton);
		cp.add(pane_top, BorderLayout.PAGE_END);

		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);  // Exit program if close-window button clicked
		setTitle("Image Scaling"); // "super" Frame sets title
		setSize(1100, 1000);  // "super" Frame sets initial size
		setVisible(true);   // "super" Frame shows
	}

	private ChangeListener createListner(JTextField textField) {
		ChangeListener listener = new ChangeListener()
		{
			public void stateChanged(ChangeEvent event)
			{
				// update text field when the slider value changes
				JSlider source = (JSlider) event.getSource();
				textField.setText(""+source.getValue());
			}
		};
		return listener;
	}

	private void createOKbutton() {
		OKbutton = new JButton("OK");
		OKbutton.addActionListener( new ActionListener()
		{
			@Override
			public void actionPerformed(ActionEvent e)
			{
				System.out.println("\nOK button is pressed ****************");
				//get the selected item:
//				int step = stepSlider.getValue();
//				int klabel = kLabelSlider.getValue();
				int step = Integer.parseInt(jTextStep.getText());
				int klabel = Integer.parseInt(jTextKsize.getText());
				//jTextKsize.setText("klabel"+Integer.toString(klabel));
				//jTextStep.setText("step"+Integer.toString(step));
				System.out.println("You seleted: " + step + " steps");
				System.out.println("You seleted: " + klabel + " labels");

				double sigmaColor = Double.parseDouble(jTextsigmaColorSim.getText()); 
				double sigmaGradient = Double.parseDouble(jTextsigmaGradientSim.getText());
				double Beta = Double.parseDouble(jTextBeta.getText());


				//perform the process
				try {
					process(step, klabel, sigmaColor, sigmaGradient, Beta);
				} catch (IOException e1) {
					e1.printStackTrace();
				}
				//update components		
				//update origin icon
				origin_ii = new ImageIcon(last_processed);
				origin_ii.getImage().flush();
				label_origin.setIcon(origin_ii);
				System.out.println("last: " + last_processed);

				//update process icon
				process_ii = new ImageIcon(processed);
				process_ii.getImage().flush();
				label_processed.setIcon(process_ii);
				System.out.println("process: " + processed);

				//update 

				//update baseNonRegion icon
				baseNonRegion_ii = new ImageIcon(baseNonRegion);
				baseNonRegion_ii.getImage().flush();
				label_baseNonRegion.setIcon(baseNonRegion_ii);
				System.out.println("baseNonRegion: " + baseNonRegion);

				//update baseNonRegion icon
				baseRegion_ii = new ImageIcon(baseRegion);
				baseRegion_ii.getImage().flush();
				label_baseRegion.setIcon(baseRegion_ii);
				System.out.println("baseRegion: " + baseRegion);

				//update colorSim icon
				colorSim_ii = new ImageIcon(colorSim);
				colorSim_ii.getImage().flush();
				label_colorSim.setIcon(colorSim_ii);
				System.out.println("colorSim: " + colorSim);

				//update gradientSim icon
				gradientSim_ii = new ImageIcon(gradientSim);
				gradientSim_ii.getImage().flush();
				label_gradientSim.setIcon(gradientSim_ii);
				System.out.println("gradientSim: " + gradientSim);

				//swap the last_process and process for next round
				last_processed = processed;
			}
		});
	}

	private void createOriginLabel() throws IOException {
		last_bufferedImage = ImageIO.read(new File(original));
		origin_ii = new ImageIcon(last_bufferedImage);
		label_origin = new JLabel();
		label_origin.setIcon(origin_ii);
		origin_bufferedImage = last_bufferedImage;
	}

	private void createProcessedLabel() throws IOException {
		processed_bufferedImage = ImageIO.read(new File(processed));
		process_ii = new ImageIcon(processed_bufferedImage);
		label_processed = new JLabel();
		label_processed.setIcon(process_ii);
	}

	private void createColorSimLabel() throws IOException {
		colorSim_bufferedImage = ImageIO.read(new File(colorSim));
		colorSim_ii = new ImageIcon(colorSim_bufferedImage);
		label_colorSim = new JLabel();
		label_colorSim.setIcon(colorSim_ii);
	}

	private void createGradientSimLabel() throws IOException {
		gradientSim_bufferedImage = ImageIO.read(new File(gradientSim));
		gradientSim_ii = new ImageIcon(gradientSim_bufferedImage);
		label_gradientSim = new JLabel();
		label_gradientSim.setIcon(gradientSim_ii);
	}

	private void createBaseNonRegionLabel() throws IOException {
		baseNonRegionImage = ImageIO.read(new File(baseNonRegion));
		baseNonRegion_ii = new ImageIcon(baseNonRegionImage);
		label_baseNonRegion = new JLabel();
		label_baseNonRegion.setIcon(baseNonRegion_ii);
	}

	private void createBaseRegionLabel() throws IOException {
		baseRegionImage = ImageIO.read(new File(baseRegion));
		baseRegion_ii = new ImageIcon(baseRegionImage);
		label_baseRegion = new JLabel();
		label_baseRegion.setIcon(baseRegion_ii);
	}

	private void createListImageSource(){
		cbg = new CheckboxGroup();
		selected_cb = new Checkbox("new selected", cbg, false);
		selected_cb.addItemListener(new ItemListener() {

			@Override
			public void itemStateChanged(ItemEvent e) {
				last_processed = fc.getSelectedFile();				
			}
		});
		origin_cb = new Checkbox("original", cbg, true);
		origin_cb.addItemListener(new ItemListener() {

			@Override
			public void itemStateChanged(ItemEvent e) {
				last_processed = original;
			}
		});
		process_cb = new Checkbox("last processed",cbg, false);	
		process_cb.addItemListener(new ItemListener() {

			@Override
			public void itemStateChanged(ItemEvent e) {
				//				last_processed = last_processed;
			}
		});
	}

	private void createListAlgorithems(){
		String[] functions = new String[] {"coarse blurriness", "CIE color show", "CIE color average", "color similarity blurriness", "Gradient X", "Gradient Y", "Gradient Average", "gradient similarity blurriness", "final"};
	}

	private void setLastProcess() {
		String selected = cbg.getSelectedCheckbox().getLabel();
		switch (selected) {
		case "original":
			last_processed = original;
			break;
		case "new selected":
			System.out.println("selected new selected");
			if(fc.getSelectedFile() != null) {
				last_processed = fc.getSelectedFile();
			}else if(last_processed == null){
				last_processed = original;
			}
			System.out.println("then last_processed: " + last_processed);
			break;
		case "last processed":
			//			last_processed = last_processed;
			break;
		default:
			break;
		}
	}

	private void process(int step, int K_size, double sigmaColor, double sigmaGradient, double beta) throws IOException {
		setLastProcess();
		System.out.println("int process(): last_processed: " + last_processed);
		last_bufferedImage = ImageIO.read(new File(last_processed));
		imProcess ip = new imProcess();

		int width = last_bufferedImage.getWidth();
		int height = last_bufferedImage.getHeight();

		int[] gray_matrix_int = new int[width*height];
		double[] gray_matrix = new double[width*height];
		float[][] rgb_matrix = new float[width*height][3];
		int[] RGBs = new int[width*height];

		ip.bufferedImage2matrix(last_bufferedImage, gray_matrix, rgb_matrix, RGBs, gray_matrix_int);

		//		double[] inreal = gray_matrix;
		//		double[] inimag = new double[inreal.length];
		//		System.out.println(Arrays.toString(inreal));
		//		System.out.println(Ar rays.toString(inimag));

		/***** debug gray scale********/
		//draw the grayscale picture
		this.processed = Draws.drawDebugGray(gray_matrix, width, height);


		/**** 1.DFT *********/
		System.out.println("in DFT********************");
		double[] real_part = imProcess.real(gray_matrix);
		//		double[][] real_imag = ip.DFT(inreal, inimag);
		//debugging
		//		ip.printdouble(real_part, width, height);


		//		System.out.println(Arrays.toString(real_imag[0]));
		//		System.out.println(Arrays.toString(real_imag[1]));

		//			double[] inrealinv = {2.0, -2.0, 0.0, 4.0};
		//			double[] inimaginv = {0.0, -2.0, -2.0, 4.0};

		//		/*** 2. inverse****/
		//			real_imag = ip.inverseDFT(real_imag[0], real_imag[1]);
		//			System.out.println(Arrays.toString(real_imag[0]));
		//			System.out.println(Arrays.toString(real_imag[1]));

		/*** 2. log****/
		System.out.println("in log");
		double[] log_output = imProcess.log(real_part);
		imProcess.printdouble(log_output, width, height);
		System.out.println("finish log");

		/** 3. Average log with (3x3)**/
		System.out.println("Average filter(3x3)");
		int kernel = 3;
		double[] averaged_log = imProcess.A_f(kernel, width, height, log_output);
		//		ip.printdouble(averaged_log, width, height);
		System.out.println("finish average log");

		/** 4. S(f) = log(real(DFT(f))) - A(f)(whicaveraged_log)**/
		System.out.println("S_f");
		double[] S_f = imProcess.S_f(log_output, averaged_log);
		//				ip.printdouble(S_f, width, height);

		/** 5. B(I)** inverse S_f**/
		System.out.println("B(I)");
		double[] B_I = imProcess.B_I(S_f);
		//				ip.printdouble(B_I, width, height);

		System.out.println("finish constructing B_I");
		//		B_I = ip.normalization(B_I);

		/** 6. B(I)bar: simple logistic regression function**/
		System.out.println("B(I)bar");
		imProcess.B_I_bar(B_I);
		//				ip.printdouble(B_I, width, height);


		/** 7. SLIC **/
		SLIC slic = new SLIC();
		int numlabels;
		boolean LAB_space = false;
		int[] klabels = new int[width*height];
		numlabels = 1;
		//		K_size = 100;//required number of superpixels
		double m = 5;//weight given to spatial distance
		//int DoSuperpixelSegmentation_ForGivenK(
		//        int[] ubuff,
		//        int width, 
		//        int height,
		//        int[] klabels,
		//        int numlabels,
		//        int K,//required number of superpixels
		//        double m,//weight given to spatial distance
		//        boolean LAB_space) {
		slic.DoSuperpixelSegmentation_ForGivenK(RGBs, width, height, klabels, numlabels, K_size, m, LAB_space);

		/**** test base map ******/

		//		int[] determine = imProcess.determinationSet(B_I, 1);
		//		Color bg = Color.DARK_GRAY;
		//		Color fg = Color.MAGENTA;
		//		String image_name = Draws.drawImageWithRegressionOnRegion(determine, klabels, bg, fg, width, height);
		String image_name = Draws.drawImageWithRateForPixels(B_I, "BaseNonRegion.png", width, height);
		this.baseNonRegion = image_name;

		int contour_color = Color.WHITE.getRGB();
		int[] RGBs_for_contour = RGBs.clone();
		image_name = imProcess.showCounterImage(RGBs_for_contour, klabels, width, height, slic, contour_color);

		System.out.println(image_name);

		/**** 9. CIElab for each pixel *****/
		//0. get the CIElab for each pixel
		double[][] CIElabs = new double[3][];//region color mean
		imProcess.rbg2CIElab(slic, CIElabs, RGBs, klabels);
		System.out.println("done: getting CIElab");
		//1. mean color of each region
		//get the max label
		int[] temp_klabels = klabels.clone();
		Arrays.sort(temp_klabels);
		int max = temp_klabels[temp_klabels.length-1];
		//number of labels
		int num_of_labels = max + 1;

		//		//debuging cielab
		//		boolean b = false;
		//		for(int i = 0; i < 3; i++) {
		//			//verify
		//			for(int j = 0; j < CIElabs[i].length; j++) {
		//				if(CIElabs[0][j] > 255 || CIElabs[0][j] < 0) {
		//					System.out.println(j + " :" + CIElabs[0][j]);
		//					b = true;
		//				}
		//				if(CIElabs[1][j] > 128 || CIElabs[1][j] < -128) {
		//					System.out.println(j + " :" + CIElabs[1][j]);
		//					b = true;
		//				}
		//				if(CIElabs[2][j] > 128 || CIElabs[2][j] < -128) {
		//					System.out.println(j + " :" + CIElabs[2][j]);
		//					b = true;
		//				}
		//			}
		//		}
		//		System.out.println("b: " + b);
		//get the mean of CIElab
		float[][] cs_mean = new float[num_of_labels][3];//region color mean
		imProcess.CIElabregionAverage(cs_mean, slic, CIElabs, klabels, num_of_labels);
		//print klabels
		//				System.out.println("printing klables");
		//				int index;
		//				for(int i = 0; i < width; i += 11) {
		//					for(int j = 0; j < height; j += 11) {
		//						index = i*width + j;
		//						System.out.printf("%3d ", klabels[index]);
		//					}
		//					System.out.println();
		//				}

		//						for(int k = 0; k < klabels.length; k++) {
		//							if(k%width == 0) {
		//								System.out.println();
		//							}
		//							//miniature
		//							if(k%1 == 0) {
		//				//				String.format("%02d", klabels[k]);
		//								System.out.printf("%3d ", klabels[k]);
		//							}
		//						}

		//print cs_mean
		//		System.out.println("printing cs_mean: " + cs_mean.length);
		//		for(int k = 0; k < cs_mean.length; k++) {
		////			if(k%10 == 0) {
		////				System.out.println();
		////			}
		//			System.out.println("mean: " + cs_mean[k][0] + " " + cs_mean[k][1] + " " +  cs_mean[k][2]);
		//		}
		/*** 10. Neighbors ***/
		HashMap<Integer, Set<Integer>> cs_neighbor = new HashMap<Integer, Set<Integer>>();
		imProcess.getNeighbors(num_of_labels, cs_neighbor, width, height, klabels);

		//				System.out.println("printing neighbors");
		//				Iterator<Entry<Integer, Set<Integer>>> iter = cs_neighbor.entrySet().iterator();
		//				while (iter.hasNext()) {
		//					Entry<Integer, Set<Integer>> entry = iter.next();
		//					System.out.println(entry.getKey() + ": " + entry.getValue().toString());
		//				}

		image_name = imProcess.drawMean(cs_mean, width, height, klabels);
		image_name = imProcess.drawDebugCIE(CIElabs, width, height, "CIE");
		System.out.println(image_name);

		/*** 11. setup B_I for regions ****/
		double[] b_0_region = new double[num_of_labels];
		imProcess.initializeCoarseBlurMap(B_I, b_0_region, cs_neighbor, klabels);
		this.baseRegion = Draws.drawImageWithRegressionOnRegionHSB(B_I, klabels, width, height, "baseRegion.png");
		System.out.println("regionize the Blur map");


		/*** 12. cs_ij ***/
		float[][] cs_ij = new float[num_of_labels][num_of_labels];
		imProcess.colorSimilarity(num_of_labels, cs_neighbor, cs_mean, sigmaColor , cs_ij);
		//		//debug
		//		for(int i = 0; i < num_of_labels; i++) {
		//			for(int j = 0; j <num_of_labels; j++) {
		//				if(cs_ij[i][j] != -1 && cs_ij[i][j] != 0.0) {
		//					System.out.print ("(" + i + "," + j + ": " + cs_ij[i][j] + ") ");
		//					System.out.print (i + ":" + cs_mean[i][0] + " " +cs_mean[i][0] + " " +cs_mean[i][0] + " VS ");
		//					System.out.println(j + ":" + cs_mean[j][0] + " " +cs_mean[j][0] + " " +cs_mean[j][0]);
		//				}
		//			}
		//			//			System.out.println();
		//		}
		System.out.println("finish color similarity");


		/***** 13. gradient similarity ********/
		int[] Gx = new int[gray_matrix.length];
		image_name = imProcess.getGradient_Xdirection(true, gray_matrix_int, Gx, width, height);
		//						this.processed = image_name;

		int[] Gy = new int[gray_matrix.length];
		image_name = imProcess.getGradient_Xdirection(false, gray_matrix_int, Gy, width, height);
		//						this.processed = image_name;

		/**** 14. get mean of gradient in X direction and Y direction*******/
		int[][] mean_gradient = new int[2][klabels.length];
		image_name = imProcess.meanGradientOfRegion(Gx, Gy, klabels, mean_gradient, num_of_labels, width, height);
		//				this.processed = image_name;


		/***** 15. covariance matrix ******/
		DenseMatrix[] coMatrices= imProcess.regionCoMat(klabels, mean_gradient, Gx, Gy, num_of_labels);
		//debugging
		//				boolean check = false;
		//				for(int i = 0; i < coMatrices.length; i++) {
		//		//			if((coMatrices[i].get(0, 0) * coMatrices[i].get(1, 1) - coMatrices[i].get(0, 1) * coMatrices[i].get(1, 0))< 0) {
		//					if(i == 10 || i == 11) {
		//					System.out.println(i + ": " + coMatrices[i]);
		//					check = true;
		//					}
		//				}
		//				System.out.println("check: " + check);

		/**** 16. AI mulplyC1negHalfC2C1negHalf ***/
		float[][] gs_ij = new float[num_of_labels][num_of_labels];
		imProcess.gradientSimilarity(num_of_labels, cs_neighbor, coMatrices, sigmaColor, gs_ij);
		//		//debug
		//		for(int i = 0; i < num_of_labels; i++) {
		//			for(int j = 0; j <num_of_labels; j++) {
		////				if(gs_ij[i][j] != -1 && gs_ij[i][j] != 0.0) {
		//				if(cs_ij[i][j] != gs_ij[i][j])
		//					System.out.println("(" + i + "," + j + ": " + cs_ij[i][j] + " and " + gs_ij[i][j] + ") ");
		////			}
		//			}
		//						System.out.println();
		//		}
		System.out.println("finish gs_ij");

		/**** 17. total similarites ****/
		//draw colSim 
//		drawColor(num_of_labels, cs_ij, gs_ij, b_0_region, step, width, height, temp_klabels);

		//draw gradientSim 
//		drawGradient(num_of_labels, cs_ij, gs_ij, b_0_region, step, width, height, temp_klabels);


		float[][] s_ij = new float[num_of_labels][num_of_labels];
		imProcess.similarites(cs_ij, gs_ij, s_ij, num_of_labels, beta);

		/**** 18. coherence ****/
		float[] coh_i = new float[num_of_labels];
		imProcess.coherence(s_ij, num_of_labels, coh_i);
		//		//debug
		//		for(int i = 0; i < num_of_labels; i++) {
		//			System.out.println(coh_i[i]);
		//		}
		//		
		/*** 19. iteration ****/
		//		int step = 50;
		System.out.println("len of b_0_region: " + b_0_region.length);
		double[] blur_map = imProcess.iterationUpdate(coh_i, b_0_region, s_ij, step);

		//		System.out.println(Arrays.toString(blur_map));
		//		for(int i = 0; i < blur_map.length; i++) {
		//			System.out.print(blur_map[i] + " " );
		//			if(i%10 == 0) {
		//				System.out.println();
		//			}
		//		}	
		/**** 20. draw the final image!!!*****/
		//printing blurmap
		//		for(int i = 0; i < blur_map.length; i++) {
		//			if(i%50 == 0) {
		//				System.out.println();
		//			}
		////			if(Double.isNaN(blur_map[i]))
		//			System.out.println(i + ": blur_map:" + blur_map[i] + " " + "coh_i: " + coh_i[i] + " s_ij:" + s_ij[i] + " b_0_region:" + b_0_region[i]);
		//		}
		System.out.println("len of blur_map: " + blur_map.length);
		System.out.println("You seleted: " + step);

		//		this.processed = imProcess.drawImageWithRegressionOnRegion(determine,klabels, bg, fg, width, height);
		this.processed = Draws.drawImageWithRegressionOnRegionHSB(blur_map, klabels, width, height, "final_process");
	}

	private double[] simProcess(int num_of_labels, float[][] cs_ij, float[][] gs_ij, float beta, double[] b_0_region, int step, int width, int height, int[] klabels) {
		System.out.println("/////////// simProcess////////////////");
		float[][] s_ij = new float[num_of_labels][num_of_labels];
		imProcess.similarites(cs_ij, gs_ij, s_ij, num_of_labels, beta);

		/**** 18. coherence ****/
		float[] coh_i = new float[num_of_labels];
		imProcess.coherence(s_ij, num_of_labels, coh_i);

		/*** 19. iteration ****/
		System.out.println("len of b_0_region: " + b_0_region.length);
		double[] blur_map = imProcess.iterationUpdate(coh_i, b_0_region, s_ij, step);

		/**** 20. draw the final image!!!*****/
		System.out.println("len of blur_map: " + blur_map.length);
		System.out.println("You seleted: " + step);

		System.out.println("/////////// end of simProcess////////////////\n");
		return blur_map;
	}

	private void drawColor(int num_of_labels, float[][] cs_ij, float[][] gs_ij, double[] b_0_region, int step, int width, int height, int[] klabels) {
		double[] rate = simProcess(num_of_labels, cs_ij, gs_ij, 1, b_0_region, step, width, height, klabels);
		this.colorSim = Draws.drawImageWithRegressionOnRegionHSB(rate, klabels, width, height, "colSim.png");
	}

	private void drawGradient(int num_of_labels, float[][] cs_ij, float[][] gs_ij, double[] b_0_region, int step, int width, int height, int[] klabels) {
		double[] rate = simProcess(num_of_labels, cs_ij, gs_ij, 0, b_0_region, step, width, height, klabels);
		this.gradientSim = Draws.drawImageWithRegressionOnRegionHSB(rate, klabels, width, height, "gradientSim.png");
	}

	public static void main(String[] args) {
		// Run the GUI construction in the Event-Dispatching thread for thread-safety
		SwingUtilities.invokeLater(new Runnable() {
			@Override
			public void run() {
				try {
					new mainFrame2_organized();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		});
	}
}