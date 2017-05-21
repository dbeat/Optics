package com.dbeat.rayCurvature;

import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.ClassicalRungeKuttaIntegrator;

// see Landau & Lifschitz Vol-2 The Classical Theory of Fields
class Ray implements FirstOrderDifferentialEquations {
	
	 private double c;
	 private double n;

	    public Ray(double n) {
	        this.c     = 3e8;
	        this.n     = n;
	    }

	public void computeDerivatives(double t, double[] y, double[] yDot)
			throws MaxCountExceededException, DimensionMismatchException {
		// y[0], y[1], and y[2] = kx, ky, kz and y[3], y[4], and y[5] = qx, qy, qz
		double normk = Math.sqrt(Math.pow(y[0],2)+Math.pow(y[1],2)+Math.pow(y[2],2));
		yDot[0] = 0.0;
		yDot[1] = 0.0;
		yDot[2] = 0.0;
		yDot[3] = c /(n * normk) * y[0];
		yDot[4] = c /(n * normk) * y[1];
		yDot[5] = c /(n * normk) * y[2];
		
	}

	public int getDimension() {
		// TODO Auto-generated method stub
		return 6;
	}
	
	public static void main(String[] args) {
    	FirstOrderIntegrator dp853 = new ClassicalRungeKuttaIntegrator(1.0e-13);
    	FirstOrderDifferentialEquations ode = new Ray( 1.0 );
    	double c = 3.0e8;
    	double n = 1.0;
    	double lambda0 = 660.0e-9;
    	double omega = 2*Math.PI*c/(n*lambda0);
    	double[] L0 = new double[] {1.0, 0.0, 0.0};
    	double normL0 = Math.sqrt(Math.pow(L0[0],2)+Math.pow(L0[1],2)+Math.pow(L0[2],2));
    	double[] k0 = new double[] {omega*n*L0[0]/(c*normL0),omega*n*L0[1]/(c*normL0),omega*n*L0[2]/(c*normL0)};
    	double[] q0 = new double[] {0.0, 0.0, 0.0};
    	double[] y = new double[] { k0[0], k0[1], k0[2], q0[0], q0[1], q0[2] }; // initial state
    	double t = 1.0e-9;
    	dp853.integrate(ode, 0.0, y, t, y); // now y contains final state at time t=1.0e-9
    	System.out.println("kx: "+y[0]+" ky: "+y[1]+" kz: "+y[2]); // should give 9.5200e6 rad/m, 0.0 rad/m, 0.0 rad/m
    	System.out.println("qx: "+y[3]+" qy: "+y[4]+" qz: "+y[5]); // should give 0.29979 m, 0.0 m, 0.0 m
    }
}
