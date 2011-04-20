
using System;

namespace alglib
{
    public class testidwunit
    {
        /*************************************************************************
        Testing IDW interpolation
        *************************************************************************/
        public static bool testidw(bool silent)
        {
            bool result = new bool();
            double[,] xy = new double[0,0];
            int i = 0;
            int j = 0;
            double vx = 0;
            double vy = 0;
            double vz = 0;
            int d = 0;
            int dtask = 0;
            int nx = 0;
            int n = 0;
            int nq = 0;
            int nw = 0;
            int smalln = 0;
            int largen = 0;
            bool waserrors = new bool();
            bool idwerrors = new bool();

            idwerrors = false;
            smalln = 256;
            largen = 1024;
            nq = 10;
            nw = 18;
            
            //
            // Simple test:
            // * F = x^3 + sin(pi*y)*z^2 - (x+y)^2
            // * space is either R1=[-1,+1] (other dimensions are
            //   fixed at 0), R1^2 or R1^3.
            //* D = -1, 0, 1, 2
            //
            for(nx=1; nx<=2; nx++)
            {
                xy = new double[largen, nx+1];
                for(i=0; i<=largen-1; i++)
                {
                    for(j=0; j<=nx-1; j++)
                    {
                        xy[i,j] = 2*AP.Math.RandomReal()-1;
                    }
                    if( nx>=1 )
                    {
                        vx = xy[i,0];
                    }
                    else
                    {
                        vx = 0;
                    }
                    if( nx>=2 )
                    {
                        vy = xy[i,1];
                    }
                    else
                    {
                        vy = 0;
                    }
                    if( nx>=3 )
                    {
                        vz = xy[i,2];
                    }
                    else
                    {
                        vz = 0;
                    }
                    xy[i,nx] = vx*vx*vx+Math.Sin(Math.PI*vy)*AP.Math.Sqr(vz)-AP.Math.Sqr(vx+vy);
                }
                for(d=-1; d<=2; d++)
                {
                    testxy(ref xy, largen, nx, d, nq, nw, ref idwerrors);
                }
            }
            
            //
            // Another simple test:
            // * five points in 2D - (0,0), (0,1), (1,0), (-1,0) (0,-1)
            // * F is random
            // * D = -1, 0, 1, 2
            //
            nx = 2;
            xy = new double[5, nx+1];
            xy[0,0] = 0;
            xy[0,1] = 0;
            xy[0,2] = 2*AP.Math.RandomReal()-1;
            xy[1,0] = 1;
            xy[1,1] = 0;
            xy[1,2] = 2*AP.Math.RandomReal()-1;
            xy[2,0] = 0;
            xy[2,1] = 1;
            xy[2,2] = 2*AP.Math.RandomReal()-1;
            xy[3,0] = -1;
            xy[3,1] = 0;
            xy[3,2] = 2*AP.Math.RandomReal()-1;
            xy[4,0] = 0;
            xy[4,1] = -1;
            xy[4,2] = 2*AP.Math.RandomReal()-1;
            for(d=-1; d<=2; d++)
            {
                testxy(ref xy, 5, nx, d, nq, nw, ref idwerrors);
            }
            
            //
            // Degree test.
            //
            // F is either:
            // * constant (DTask=0)
            // * linear (DTask=1)
            // * quadratic (DTask=2)
            //
            // Nodal functions are either
            // * constant (D=0)
            // * linear (D=1)
            // * quadratic (D=2)
            //
            // When DTask<=D, we can interpolate without errors.
            // When DTask>D, we MUST have errors.
            //
            for(nx=1; nx<=3; nx++)
            {
                for(d=0; d<=2; d++)
                {
                    for(dtask=0; dtask<=2; dtask++)
                    {
                        testdegree(smalln, nx, d, dtask, ref idwerrors);
                    }
                }
            }
            
            //
            // Noisy test
            //
            testnoisy(ref idwerrors);
            
            //
            // report
            //
            waserrors = idwerrors;
            if( !silent )
            {
                System.Console.Write("TESTING INVERSE DISTANCE WEIGHTING");
                System.Console.WriteLine();
                System.Console.Write("* IDW:                                   ");
                if( !idwerrors )
                {
                    System.Console.Write("OK");
                    System.Console.WriteLine();
                }
                else
                {
                    System.Console.Write("FAILED");
                    System.Console.WriteLine();
                }
                if( waserrors )
                {
                    System.Console.Write("TEST FAILED");
                    System.Console.WriteLine();
                }
                else
                {
                    System.Console.Write("TEST PASSED");
                    System.Console.WriteLine();
                }
                System.Console.WriteLine();
                System.Console.WriteLine();
            }
            result = !waserrors;
            return result;
        }


        /*************************************************************************
        Unsets 2D array.
        *************************************************************************/
        private static void unset2d(ref AP.Complex[,] a)
        {
            a = new AP.Complex[0+1, 0+1];
            a[0,0] = 2*AP.Math.RandomReal()-1;
        }


        /*************************************************************************
        Unsets 1D array.
        *************************************************************************/
        private static void unset1d(ref double[] a)
        {
            a = new double[0+1];
            a[0] = 2*AP.Math.RandomReal()-1;
        }


        /*************************************************************************
        Testing IDW:
        * generate model using N/NX/D/NQ/NW
        * test basic properties
        *************************************************************************/
        private static void testxy(ref double[,] xy,
            int n,
            int nx,
            int d,
            int nq,
            int nw,
            ref bool idwerrors)
        {
            double threshold = 0;
            double lipschitzstep = 0;
            int i = 0;
            int j = 0;
            int i1 = 0;
            int i2 = 0;
            double v = 0;
            double v1 = 0;
            double v2 = 0;
            double t = 0;
            double l1 = 0;
            double l2 = 0;
            idwint.idwinterpolant z1 = new idwint.idwinterpolant();
            double[] x = new double[0];
            int i_ = 0;

            threshold = 1000*AP.Math.MachineEpsilon;
            lipschitzstep = 0.001;
            x = new double[nx];
            
            //
            // build
            //
            idwint.idwbuildmodifiedshepard(ref xy, n, nx, d, nq, nw, ref z1);
            
            //
            // first, test interpolation properties at nodes
            //
            for(i=0; i<=n-1; i++)
            {
                for(i_=0; i_<=nx-1;i_++)
                {
                    x[i_] = xy[i,i_];
                }
                idwerrors = idwerrors | (double)(idwint.idwcalc(ref z1, ref x))!=(double)(xy[i,nx]);
            }
            
            //
            // test Lipschitz continuity
            //
            i1 = AP.Math.RandomInteger(n);
            do
            {
                i2 = AP.Math.RandomInteger(n);
            }
            while( i2==i1 );
            l1 = 0;
            t = 0;
            while( (double)(t)<(double)(1) )
            {
                v = 1-t;
                for(i_=0; i_<=nx-1;i_++)
                {
                    x[i_] = v*xy[i1,i_];
                }
                v = t;
                for(i_=0; i_<=nx-1;i_++)
                {
                    x[i_] = x[i_] + v*xy[i2,i_];
                }
                v1 = idwint.idwcalc(ref z1, ref x);
                v = 1-(t+lipschitzstep);
                for(i_=0; i_<=nx-1;i_++)
                {
                    x[i_] = v*xy[i1,i_];
                }
                v = t+lipschitzstep;
                for(i_=0; i_<=nx-1;i_++)
                {
                    x[i_] = x[i_] + v*xy[i2,i_];
                }
                v2 = idwint.idwcalc(ref z1, ref x);
                l1 = Math.Max(l1, Math.Abs(v2-v1)/lipschitzstep);
                t = t+lipschitzstep;
            }
            l2 = 0;
            t = 0;
            while( (double)(t)<(double)(1) )
            {
                v = 1-t;
                for(i_=0; i_<=nx-1;i_++)
                {
                    x[i_] = v*xy[i1,i_];
                }
                v = t;
                for(i_=0; i_<=nx-1;i_++)
                {
                    x[i_] = x[i_] + v*xy[i2,i_];
                }
                v1 = idwint.idwcalc(ref z1, ref x);
                v = 1-(t+lipschitzstep/3);
                for(i_=0; i_<=nx-1;i_++)
                {
                    x[i_] = v*xy[i1,i_];
                }
                v = t+lipschitzstep/3;
                for(i_=0; i_<=nx-1;i_++)
                {
                    x[i_] = x[i_] + v*xy[i2,i_];
                }
                v2 = idwint.idwcalc(ref z1, ref x);
                l2 = Math.Max(l2, Math.Abs(v2-v1)/(lipschitzstep/3));
                t = t+lipschitzstep/3;
            }
            idwerrors = idwerrors | (double)(l2)>(double)(2.0*l1);
        }


        /*************************************************************************
        Testing IDW:
        * generate model using R-based model
        * test basic properties
        *************************************************************************/
        private static void testrxy(ref double[,] xy,
            int n,
            int nx,
            double r,
            ref bool idwerrors)
        {
            double threshold = 0;
            double lipschitzstep = 0;
            int i = 0;
            int j = 0;
            int i1 = 0;
            int i2 = 0;
            double v = 0;
            double v1 = 0;
            double v2 = 0;
            double t = 0;
            double l1 = 0;
            double l2 = 0;
            idwint.idwinterpolant z1 = new idwint.idwinterpolant();
            double[] x = new double[0];
            int i_ = 0;

            threshold = 1000*AP.Math.MachineEpsilon;
            lipschitzstep = 0.001;
            x = new double[nx];
            
            //
            // build
            //
            idwint.idwbuildmodifiedshepardr(ref xy, n, nx, r, ref z1);
            
            //
            // first, test interpolation properties at nodes
            //
            for(i=0; i<=n-1; i++)
            {
                for(i_=0; i_<=nx-1;i_++)
                {
                    x[i_] = xy[i,i_];
                }
                idwerrors = idwerrors | (double)(idwint.idwcalc(ref z1, ref x))!=(double)(xy[i,nx]);
            }
            
            //
            // test Lipschitz continuity
            //
            i1 = AP.Math.RandomInteger(n);
            do
            {
                i2 = AP.Math.RandomInteger(n);
            }
            while( i2==i1 );
            l1 = 0;
            t = 0;
            while( (double)(t)<(double)(1) )
            {
                v = 1-t;
                for(i_=0; i_<=nx-1;i_++)
                {
                    x[i_] = v*xy[i1,i_];
                }
                v = t;
                for(i_=0; i_<=nx-1;i_++)
                {
                    x[i_] = x[i_] + v*xy[i2,i_];
                }
                v1 = idwint.idwcalc(ref z1, ref x);
                v = 1-(t+lipschitzstep);
                for(i_=0; i_<=nx-1;i_++)
                {
                    x[i_] = v*xy[i1,i_];
                }
                v = t+lipschitzstep;
                for(i_=0; i_<=nx-1;i_++)
                {
                    x[i_] = x[i_] + v*xy[i2,i_];
                }
                v2 = idwint.idwcalc(ref z1, ref x);
                l1 = Math.Max(l1, Math.Abs(v2-v1)/lipschitzstep);
                t = t+lipschitzstep;
            }
            l2 = 0;
            t = 0;
            while( (double)(t)<(double)(1) )
            {
                v = 1-t;
                for(i_=0; i_<=nx-1;i_++)
                {
                    x[i_] = v*xy[i1,i_];
                }
                v = t;
                for(i_=0; i_<=nx-1;i_++)
                {
                    x[i_] = x[i_] + v*xy[i2,i_];
                }
                v1 = idwint.idwcalc(ref z1, ref x);
                v = 1-(t+lipschitzstep/3);
                for(i_=0; i_<=nx-1;i_++)
                {
                    x[i_] = v*xy[i1,i_];
                }
                v = t+lipschitzstep/3;
                for(i_=0; i_<=nx-1;i_++)
                {
                    x[i_] = x[i_] + v*xy[i2,i_];
                }
                v2 = idwint.idwcalc(ref z1, ref x);
                l2 = Math.Max(l2, Math.Abs(v2-v1)/(lipschitzstep/3));
                t = t+lipschitzstep/3;
            }
            idwerrors = idwerrors | (double)(l2)>(double)(2.0*l1);
        }


        /*************************************************************************
        Testing degree properties

        F is either:
        * constant (DTask=0)
        * linear (DTask=1)
        * quadratic (DTask=2)

        Nodal functions are either
        * constant (D=0)
        * linear (D=1)
        * quadratic (D=2)

        When DTask<=D, we can interpolate without errors.
        When DTask>D, we MUST have errors.
        *************************************************************************/
        private static void testdegree(int n,
            int nx,
            int d,
            int dtask,
            ref bool idwerrors)
        {
            double threshold = 0;
            int nq = 0;
            int nw = 0;
            int i = 0;
            int j = 0;
            double v = 0;
            double c0 = 0;
            double[] c1 = new double[0];
            double[,] c2 = new double[0,0];
            double[] x = new double[0];
            double[,] xy = new double[0,0];
            idwint.idwinterpolant z1 = new idwint.idwinterpolant();
            double v1 = 0;
            double v2 = 0;
            int i_ = 0;

            threshold = 1.0E6*AP.Math.MachineEpsilon;
            nq = 2*(nx*nx+nx+1);
            nw = 10;
            System.Diagnostics.Debug.Assert(nq<=n, "TestDegree: internal error");
            
            //
            // prepare model
            //
            c0 = 2*AP.Math.RandomReal()-1;
            c1 = new double[nx];
            for(i=0; i<=nx-1; i++)
            {
                c1[i] = 2*AP.Math.RandomReal()-1;
            }
            c2 = new double[nx, nx];
            for(i=0; i<=nx-1; i++)
            {
                for(j=i+1; j<=nx-1; j++)
                {
                    c2[i,j] = 2*AP.Math.RandomReal()-1;
                    c2[j,i] = c2[i,j];
                }
                do
                {
                    c2[i,i] = 2*AP.Math.RandomReal()-1;
                }
                while( (double)(Math.Abs(c2[i,i]))<=(double)(0.3) );
            }
            
            //
            // prepare points
            //
            xy = new double[n, nx+1];
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=nx-1; j++)
                {
                    xy[i,j] = 4*AP.Math.RandomReal()-2;
                }
                xy[i,nx] = c0;
                if( dtask>=1 )
                {
                    v = 0.0;
                    for(i_=0; i_<=nx-1;i_++)
                    {
                        v += c1[i_]*xy[i,i_];
                    }
                    xy[i,nx] = xy[i,nx]+v;
                }
                if( dtask==2 )
                {
                    for(j=0; j<=nx-1; j++)
                    {
                        v = 0.0;
                        for(i_=0; i_<=nx-1;i_++)
                        {
                            v += c2[j,i_]*xy[i,i_];
                        }
                        xy[i,nx] = xy[i,nx]+xy[i,j]*v;
                    }
                }
            }
            
            //
            // build interpolant, calculate value at random point
            //
            idwint.idwbuildmodifiedshepard(ref xy, n, nx, d, nq, nw, ref z1);
            x = new double[nx];
            for(i=0; i<=nx-1; i++)
            {
                x[i] = 4*AP.Math.RandomReal()-2;
            }
            v1 = idwint.idwcalc(ref z1, ref x);
            
            //
            // calculate model value at the same point
            //
            v2 = c0;
            if( dtask>=1 )
            {
                v = 0.0;
                for(i_=0; i_<=nx-1;i_++)
                {
                    v += c1[i_]*x[i_];
                }
                v2 = v2+v;
            }
            if( dtask==2 )
            {
                for(j=0; j<=nx-1; j++)
                {
                    v = 0.0;
                    for(i_=0; i_<=nx-1;i_++)
                    {
                        v += c2[j,i_]*x[i_];
                    }
                    v2 = v2+x[j]*v;
                }
            }
            
            //
            // Compare
            //
            if( dtask<=d )
            {
                idwerrors = idwerrors | (double)(Math.Abs(v2-v1))>(double)(threshold);
            }
            else
            {
                idwerrors = idwerrors | (double)(Math.Abs(v2-v1))<(double)(threshold);
            }
        }


        /*************************************************************************
        Noisy test:
         * F = x^2 + y^2 + z^2 + noise on [-1,+1]^3
         * space is either R1=[-1,+1] (other dimensions are
           fixed at 0), R1^2 or R1^3.
         * D = 1, 2
         * 4096 points is used for function generation,
           4096 points - for testing
         * RMS error of "noisy" model on test set must be
           lower than RMS error of interpolation model.
        *************************************************************************/
        private static void testnoisy(ref bool idwerrors)
        {
            double noiselevel = 0;
            int nq = 0;
            int nw = 0;
            int d = 0;
            int nx = 0;
            int ntrn = 0;
            int ntst = 0;
            int i = 0;
            int j = 0;
            double v = 0;
            double t = 0;
            double v1 = 0;
            double v2 = 0;
            double ve = 0;
            double[,] xy = new double[0,0];
            double[] x = new double[0];
            idwint.idwinterpolant z1 = new idwint.idwinterpolant();
            idwint.idwinterpolant z2 = new idwint.idwinterpolant();
            double rms1 = 0;
            double rms2 = 0;

            nq = 20;
            nw = 40;
            noiselevel = 0.2;
            ntrn = 256;
            ntst = 1024;
            for(d=1; d<=2; d++)
            {
                for(nx=1; nx<=2; nx++)
                {
                    
                    //
                    // prepare dataset
                    //
                    xy = new double[ntrn, nx+1];
                    for(i=0; i<=ntrn-1; i++)
                    {
                        v = noiselevel*(2*AP.Math.RandomReal()-1);
                        for(j=0; j<=nx-1; j++)
                        {
                            t = 2*AP.Math.RandomReal()-1;
                            v = v+AP.Math.Sqr(t);
                            xy[i,j] = t;
                        }
                        xy[i,nx] = v;
                    }
                    
                    //
                    // build interpolants
                    //
                    idwint.idwbuildmodifiedshepard(ref xy, ntrn, nx, d, nq, nw, ref z1);
                    idwint.idwbuildnoisy(ref xy, ntrn, nx, d, nq, nw, ref z2);
                    
                    //
                    // calculate RMS errors
                    //
                    x = new double[nx];
                    rms1 = 0;
                    rms2 = 0;
                    for(i=0; i<=ntst-1; i++)
                    {
                        ve = 0;
                        for(j=0; j<=nx-1; j++)
                        {
                            t = 2*AP.Math.RandomReal()-1;
                            x[j] = t;
                            ve = ve+AP.Math.Sqr(t);
                        }
                        v1 = idwint.idwcalc(ref z1, ref x);
                        v2 = idwint.idwcalc(ref z2, ref x);
                        rms1 = rms1+AP.Math.Sqr(v1-ve);
                        rms2 = rms2+AP.Math.Sqr(v2-ve);
                    }
                    idwerrors = idwerrors | (double)(rms2)>(double)(rms1);
                }
            }
        }


        /*************************************************************************
        Silent unit test
        *************************************************************************/
        public static bool testidwunit_test_silent()
        {
            bool result = new bool();

            result = testidw(true);
            return result;
        }


        /*************************************************************************
        Unit test
        *************************************************************************/
        public static bool testidwunit_test()
        {
            bool result = new bool();

            result = testidw(false);
            return result;
        }
    }
}
