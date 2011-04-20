
using System;

namespace alglib
{
    public class testkmeansunit
    {
        public static bool testkmeans(bool silent)
        {
            bool result = new bool();
            int nf = 0;
            int maxnf = 0;
            int nc = 0;
            int maxnc = 0;
            int passcount = 0;
            int pass = 0;
            bool waserrors = new bool();
            bool converrors = new bool();
            bool simpleerrors = new bool();
            bool complexerrors = new bool();
            bool othererrors = new bool();
            bool restartserrors = new bool();

            
            //
            // Primary settings
            //
            maxnf = 5;
            maxnc = 5;
            passcount = 10;
            waserrors = false;
            converrors = false;
            othererrors = false;
            simpleerrors = false;
            complexerrors = false;
            restartserrors = false;
            
            //
            //
            //
            for(nf=1; nf<=maxnf; nf++)
            {
                for(nc=1; nc<=maxnc; nc++)
                {
                    simpletest1(nf, nc, passcount, ref converrors, ref othererrors, ref simpleerrors);
                }
            }
            restartstest(ref converrors, ref restartserrors);
            
            //
            // Final report
            //
            waserrors = converrors | othererrors | simpleerrors | complexerrors | restartserrors;
            if( !silent )
            {
                System.Console.Write("K-MEANS TEST");
                System.Console.WriteLine();
                System.Console.Write("TOTAL RESULTS:                           ");
                if( !waserrors )
                {
                    System.Console.Write("OK");
                    System.Console.WriteLine();
                }
                else
                {
                    System.Console.Write("FAILED");
                    System.Console.WriteLine();
                }
                System.Console.Write("* CONVERGENCE:                           ");
                if( !converrors )
                {
                    System.Console.Write("OK");
                    System.Console.WriteLine();
                }
                else
                {
                    System.Console.Write("FAILED");
                    System.Console.WriteLine();
                }
                System.Console.Write("* SIMPLE TASKS:                          ");
                if( !simpleerrors )
                {
                    System.Console.Write("OK");
                    System.Console.WriteLine();
                }
                else
                {
                    System.Console.Write("FAILED");
                    System.Console.WriteLine();
                }
                System.Console.Write("* COMPLEX TASKS:                         ");
                if( !complexerrors )
                {
                    System.Console.Write("OK");
                    System.Console.WriteLine();
                }
                else
                {
                    System.Console.Write("FAILED");
                    System.Console.WriteLine();
                }
                System.Console.Write("* OTHER PROPERTIES:                      ");
                if( !othererrors )
                {
                    System.Console.Write("OK");
                    System.Console.WriteLine();
                }
                else
                {
                    System.Console.Write("FAILED");
                    System.Console.WriteLine();
                }
                System.Console.Write("* RESTARTS PROPERTIES:                   ");
                if( !restartserrors )
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
                    System.Console.Write("TEST SUMMARY: FAILED");
                    System.Console.WriteLine();
                }
                else
                {
                    System.Console.Write("TEST SUMMARY: PASSED");
                    System.Console.WriteLine();
                }
                System.Console.WriteLine();
                System.Console.WriteLine();
            }
            result = !waserrors;
            return result;
        }


        /*************************************************************************
        Simple test 1: ellipsoid in NF-dimensional space.
        compare k-means centers with random centers
        *************************************************************************/
        private static void simpletest1(int nvars,
            int nc,
            int passcount,
            ref bool converrors,
            ref bool othererrors,
            ref bool simpleerrors)
        {
            int npoints = 0;
            int majoraxis = 0;
            double[,] xy = new double[0,0];
            double[] tmp = new double[0];
            double v = 0;
            int i = 0;
            int j = 0;
            int info = 0;
            double[,] c = new double[0,0];
            int[] xyc = new int[0];
            int pass = 0;
            int restarts = 0;
            double ekmeans = 0;
            double erandom = 0;
            double dclosest = 0;
            int cclosest = 0;
            int i_ = 0;

            npoints = nc*100;
            restarts = 5;
            passcount = 10;
            tmp = new double[nvars-1+1];
            for(pass=1; pass<=passcount; pass++)
            {
                
                //
                // Fill
                //
                xy = new double[npoints-1+1, nvars-1+1];
                majoraxis = AP.Math.RandomInteger(nvars);
                for(i=0; i<=npoints-1; i++)
                {
                    rsphere(ref xy, nvars, i);
                    xy[i,majoraxis] = nc*xy[i,majoraxis];
                }
                
                //
                // Test
                //
                kmeans.kmeansgenerate(ref xy, npoints, nvars, nc, restarts, ref info, ref c, ref xyc);
                if( info<0 )
                {
                    converrors = true;
                    return;
                }
                
                //
                // Test that XYC is correct mapping to cluster centers
                //
                for(i=0; i<=npoints-1; i++)
                {
                    cclosest = -1;
                    dclosest = AP.Math.MaxRealNumber;
                    for(j=0; j<=nc-1; j++)
                    {
                        for(i_=0; i_<=nvars-1;i_++)
                        {
                            tmp[i_] = xy[i,i_];
                        }
                        for(i_=0; i_<=nvars-1;i_++)
                        {
                            tmp[i_] = tmp[i_] - c[i_,j];
                        }
                        v = 0.0;
                        for(i_=0; i_<=nvars-1;i_++)
                        {
                            v += tmp[i_]*tmp[i_];
                        }
                        if( (double)(v)<(double)(dclosest) )
                        {
                            cclosest = j;
                            dclosest = v;
                        }
                    }
                    if( cclosest!=xyc[i] )
                    {
                        othererrors = true;
                        return;
                    }
                }
                
                //
                // Use first NC rows of XY as random centers
                // (XY is totally random, so it is as good as any other choice).
                //
                // Compare potential functions.
                //
                ekmeans = 0;
                for(i=0; i<=npoints-1; i++)
                {
                    for(i_=0; i_<=nvars-1;i_++)
                    {
                        tmp[i_] = xy[i,i_];
                    }
                    for(i_=0; i_<=nvars-1;i_++)
                    {
                        tmp[i_] = tmp[i_] - c[i_,xyc[i]];
                    }
                    v = 0.0;
                    for(i_=0; i_<=nvars-1;i_++)
                    {
                        v += tmp[i_]*tmp[i_];
                    }
                    ekmeans = ekmeans+v;
                }
                erandom = 0;
                for(i=0; i<=npoints-1; i++)
                {
                    dclosest = AP.Math.MaxRealNumber;
                    for(j=0; j<=nc-1; j++)
                    {
                        for(i_=0; i_<=nvars-1;i_++)
                        {
                            tmp[i_] = xy[i,i_];
                        }
                        for(i_=0; i_<=nvars-1;i_++)
                        {
                            tmp[i_] = tmp[i_] - xy[j,i_];
                        }
                        v = 0.0;
                        for(i_=0; i_<=nvars-1;i_++)
                        {
                            v += tmp[i_]*tmp[i_];
                        }
                        if( (double)(v)<(double)(dclosest) )
                        {
                            dclosest = v;
                        }
                    }
                    erandom = erandom+v;
                }
                if( (double)(erandom)<(double)(ekmeans) )
                {
                    simpleerrors = true;
                    return;
                }
            }
        }


        /*************************************************************************
        This non-deterministic test checks that Restarts>1 significantly  improves
        quality of results.

        Subroutine generates random task 3 unit balls in 2D, each with 20  points,
        separated by 5 units wide gaps, and solves it  with  Restarts=1  and  with
        Restarts=5. Potential functions are compared,  outcome  of  the  trial  is
        either 0 or 1 (depending on what is better).

        Sequence of 1000 such tasks is  solved.  If  Restarts>1  actually  improve
        quality of solution, sum of outcome will be non-binomial.  If  it  doesn't
        matter, it will be binomially distributed.

        P.S. This test was added after report from Gianluca  Borello  who  noticed
        error in the handling of multiple restarts.
        *************************************************************************/
        private static void restartstest(ref bool converrors,
            ref bool restartserrors)
        {
            int npoints = 0;
            int nvars = 0;
            int nclusters = 0;
            int clustersize = 0;
            int restarts = 0;
            int passcount = 0;
            double sigmathreshold = 0;
            double p = 0;
            double s = 0;
            double[,] xy = new double[0,0];
            double[,] ca = new double[0,0];
            double[,] cb = new double[0,0];
            int[] xyca = new int[0];
            int[] xycb = new int[0];
            double[] tmp = new double[0];
            int i = 0;
            int j = 0;
            int info = 0;
            int pass = 0;
            double ea = 0;
            double eb = 0;
            double v = 0;
            int i_ = 0;

            restarts = 5;
            passcount = 1000;
            clustersize = 20;
            nclusters = 3;
            nvars = 2;
            npoints = nclusters*clustersize;
            sigmathreshold = 5;
            xy = new double[npoints, nvars];
            tmp = new double[nvars];
            p = 0;
            for(pass=1; pass<=passcount; pass++)
            {
                
                //
                // Fill
                //
                for(i=0; i<=npoints-1; i++)
                {
                    rsphere(ref xy, nvars, i);
                    for(j=0; j<=nvars-1; j++)
                    {
                        xy[i,j] = xy[i,j]+(double)(i)/(double)(clustersize)*5;
                    }
                }
                
                //
                // Test: Restarts=1
                //
                kmeans.kmeansgenerate(ref xy, npoints, nvars, nclusters, 1, ref info, ref ca, ref xyca);
                if( info<0 )
                {
                    converrors = true;
                    return;
                }
                ea = 0;
                for(i=0; i<=npoints-1; i++)
                {
                    for(i_=0; i_<=nvars-1;i_++)
                    {
                        tmp[i_] = xy[i,i_];
                    }
                    for(i_=0; i_<=nvars-1;i_++)
                    {
                        tmp[i_] = tmp[i_] - ca[i_,xyca[i]];
                    }
                    v = 0.0;
                    for(i_=0; i_<=nvars-1;i_++)
                    {
                        v += tmp[i_]*tmp[i_];
                    }
                    ea = ea+v;
                }
                
                //
                // Test: Restarts>1
                //
                kmeans.kmeansgenerate(ref xy, npoints, nvars, nclusters, restarts, ref info, ref cb, ref xycb);
                if( info<0 )
                {
                    converrors = true;
                    return;
                }
                eb = 0;
                for(i=0; i<=npoints-1; i++)
                {
                    for(i_=0; i_<=nvars-1;i_++)
                    {
                        tmp[i_] = xy[i,i_];
                    }
                    for(i_=0; i_<=nvars-1;i_++)
                    {
                        tmp[i_] = tmp[i_] - cb[i_,xycb[i]];
                    }
                    v = 0.0;
                    for(i_=0; i_<=nvars-1;i_++)
                    {
                        v += tmp[i_]*tmp[i_];
                    }
                    eb = eb+v;
                }
                
                //
                // Calculate statistic.
                //
                if( (double)(ea)<(double)(eb) )
                {
                    p = p+1;
                }
                if( (double)(ea)==(double)(eb) )
                {
                    p = p+0.5;
                }
            }
            
            //
            // If Restarts doesn't influence quality of centers found, P must be
            // binomially distributed random value with mean 0.5*PassCount and
            // standard deviation Sqrt(PassCount/4).
            //
            // If Restarts do influence quality of solution, P must be significantly
            // lower than 0.5*PassCount.
            //
            s = (p-0.5*passcount)/Math.Sqrt((double)(passcount)/(double)(4));
            restartserrors = restartserrors | (double)(s)>(double)(-sigmathreshold);
        }


        /*************************************************************************
        Random normal number
        *************************************************************************/
        private static double rnormal()
        {
            double result = 0;
            double u = 0;
            double v = 0;
            double s = 0;
            double x1 = 0;
            double x2 = 0;

            while( true )
            {
                u = 2*AP.Math.RandomReal()-1;
                v = 2*AP.Math.RandomReal()-1;
                s = AP.Math.Sqr(u)+AP.Math.Sqr(v);
                if( (double)(s)>(double)(0) & (double)(s)<(double)(1) )
                {
                    s = Math.Sqrt(-(2*Math.Log(s)/s));
                    x1 = u*s;
                    x2 = v*s;
                    break;
                }
            }
            result = x1;
            return result;
        }


        /*************************************************************************
        Random point from sphere
        *************************************************************************/
        private static double rsphere(ref double[,] xy,
            int n,
            int i)
        {
            double result = 0;
            int j = 0;
            double v = 0;
            int i_ = 0;

            for(j=0; j<=n-1; j++)
            {
                xy[i,j] = rnormal();
            }
            v = 0.0;
            for(i_=0; i_<=n-1;i_++)
            {
                v += xy[i,i_]*xy[i,i_];
            }
            v = AP.Math.RandomReal()/Math.Sqrt(v);
            for(i_=0; i_<=n-1;i_++)
            {
                xy[i,i_] = v*xy[i,i_];
            }
            return result;
        }


        /*************************************************************************
        Silent unit test
        *************************************************************************/
        public static bool testkmeansunit_test_silent()
        {
            bool result = new bool();

            result = testkmeans(true);
            return result;
        }


        /*************************************************************************
        Unit test
        *************************************************************************/
        public static bool testkmeansunit_test()
        {
            bool result = new bool();

            result = testkmeans(false);
            return result;
        }
    }
}
