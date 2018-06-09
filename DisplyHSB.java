package image_scaling_UI;

import java.awt.Color;

public class DisplyHSB {
	/**
	 * Remaps x from old-interval to new-interval.
	 * DoubleInterval just wraps double values a, b.
	 */
	public static double remap(double x, DoubleInterval oldDomain, DoubleInterval newDomain) {
	    double oldRange = oldDomain.size(); // oldDomain.b - oldDomain.a
	    double oldRangeValue = x - oldDomain.a; // lowerBound a is just an offset
	    double percentile = oldRangeValue / oldRange;

	    double newRange = newDomain.size(); // newDomain.b - newDomain.a
	    double newRangeValue = percentile * newRange;
	    double newVal = newRangeValue + newDomain.a;
	    return newVal;
	}

	/**
	 * Returns color from specified color-interval that matches normValue <0,1>.
	 * If normValue = 0, angleFrom = 0 (red), angleTo = 120 (green) then color = red. 
	 */
	public static Color intervalColor(double normValue, float angleFrom, float angleTo) {
	    double angleValue = remap(normValue, new DoubleInterval(0, 1), new DoubleInterval(angleFrom, angleTo));
	    return Color.getHSBColor((float) angleValue, 1, 1);        
	}

	/**
	 * Returns inversion of specified value in given domain.
	 * Example: if x = 0.3 and domain of x is <0,1> then inversion = 0.7
	 */
	public static double invert(double x, DoubleInterval domain) {
	    return (domain.b - x) + domain.a;
	}

	/**
	 * Returns color from specified color-interval that matches inverted normValue <0,1>.
	 * If normValue = 0 and angleFrom = 0 (red) and angleTo = 120 (green) then color = green. 
	 */
	public static Color invertedIntervalColor(float normValue, float angleFrom, float angleTo) {
	    double invNormValue = invert(normValue, new DoubleInterval(0, 1));
	    return intervalColor(invNormValue, angleFrom, angleTo);      
	}
	
	public static Color[] intervalHSBs(double[] rates) {
		Color[] c = new Color[rates.length];
		for(int i = 0; i < c.length; i++) {
			c[i] = intervalColor(rates[i], 0, 120);
		}
		return c;
	}
//	
//	public static void main(String[] args) {
//		
//
//	}
}
