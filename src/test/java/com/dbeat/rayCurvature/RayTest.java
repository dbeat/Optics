package com.dbeat.rayCurvature;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.ClassicalRungeKuttaIntegrator;

import com.dbeat.rayCurvature.Ray;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * Unit test for simple Ray.
 */
public class RayTest 
    extends TestCase
{
    /**
     * Create the test case
     *
     * @param testName name of the test case
     */
    public RayTest( String testName )
    {
        super( testName );
    }

    /**
     * @return the suite of tests being tested
     */
    public static Test suite()
    {
        return new TestSuite( RayTest.class );
    }

    /**
     * Rigourous Test :-)
     */
    public void testRay()
    {
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
    	dp853.integrate(ode, 0.0, y, t, y); // now y contains final state at time t=1.0 ns
        assertEquals(9.5200e6, y[0], 100.0); // x-component of the wavevector after 1 ns
        assertEquals(0.29979, y[3], 1.0e-3); // x-coordinate position of the ray after 1 ns
    }
}
