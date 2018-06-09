package image_scaling_UI;

public class DoubleInterval {
		double a;
		double b;
		public DoubleInterval(double a, double b) {
			this.a = a;
			this.b = b;
		}
		public double size() {
			return this.b - this.a;
		}
}
