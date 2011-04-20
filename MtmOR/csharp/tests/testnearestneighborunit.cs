
using System;

namespace alglib
{
    public class testnearestneighborunit
    {
        /*************************************************************************
        Testing Nearest Neighbor Search
        *************************************************************************/
        public static bool testnearestneighbor(bool silent)
        {
            bool result = new bool();
            double[,] xy = new double[0,0];
            int i = 0;
            int j = 0;
            double v = 0;
            int normtype = 0;
            int nx = 0;
            int ny = 0;
            int n = 0;
            int smalln = 0;
            int largen = 0;
            int passcount = 0;
            int pass = 0;
            bool waserrors = new bool();
            bool kdterrors = new bool();

            kdterrors = false;
            passcount = 2;
            smalln = 256;
            largen = 2048;
            ny = 3;
            
            //
            //
            //
            for(pass=1; pass<=passcount; pass++)
            {
                for(normtype=0; normtype<=2; normtype++)
                {
                    for(nx=1; nx<=3; nx++)
                    {
                        
                        //
                        // Test in hypercube
                        //
                        xy = new double[largen, nx+ny];
                        for(i=0; i<=largen-1; i++)
                        {
                            for(j=0; j<=nx+ny-1; j++)
                            {
                                xy[i,j] = 10*AP.Math.RandomReal()-5;
                            }
                        }
                        for(n=1; n<=10; n++)
                        {
                            testkdtuniform(ref xy, n, nx, AP.Math.RandomInteger(ny+1), normtype, ref kdterrors);
                        }
                        testkdtuniform(ref xy, largen, nx, AP.Math.RandomInteger(ny+1), normtype, ref kdterrors);
                        
                        //
                        // Test clustered (2*N points, pairs of equal points)
                        //
                        xy = new double[2*smalln, nx+ny];
                        for(i=0; i<=smalln-1; i++)
                        {
                            for(j=0; j<=nx+ny-1; j++)
                            {
                                xy[2*i+0,j] = 10*AP.Math.RandomReal()-5;
                                xy[2*i+1,j] = xy[2*i+0,j];
                            }
                        }
                        testkdtuniform(ref xy, 2*smalln, nx, AP.Math.RandomInteger(ny+1), normtype, ref kdterrors);
                        
                        //
                        // Test degenerate case: all points are same except for one
                        //
                        xy = new double[smalln, nx+ny];
                        v = AP.Math.RandomReal();
                        for(i=0; i<=smalln-2; i++)
                        {
                            for(j=0; j<=nx+ny-1; j++)
                            {
                                xy[i,j] = v;
                            }
                        }
                        for(j=0; j<=nx+ny-1; j++)
                        {
                            xy[smalln-1,j] = 10*AP.Math.RandomReal()-5;
                        }
                        testkdtuniform(ref xy, smalln, nx, AP.Math.RandomInteger(ny+1), normtype, ref kdterrors);
                    }
                }
            }
            
            //
            // report
            //
            waserrors = kdterrors;
            if( !silent )
            {
                System.Console.Write("TESTING NEAREST NEIGHBOR SEARCH");
                System.Console.WriteLine();
                System.Console.Write("* KD TREES:                              ");
                if( !kdterrors )
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
        Compare results from different queries:
        * X     just X-values
        * XY    X-values and Y-values
        * XT    X-values and tag values
        *************************************************************************/
        private static bool kdtresultsdifferent(ref double[,] refxy,
            int ntotal,
            ref double[,] qx,
            ref double[,] qxy,
            ref int[] qt,
            int n,
            int nx,
            int ny)
        {
            bool result = new bool();
            int i = 0;
            int j = 0;

            result = false;
            for(i=0; i<=n-1; i++)
            {
                if( qt[i]<0 | qt[i]>=ntotal )
                {
                    result = true;
                    return result;
                }
                for(j=0; j<=nx-1; j++)
                {
                    result = result | (double)(qx[i,j])!=(double)(refxy[qt[i],j]);
                    result = result | (double)(qxy[i,j])!=(double)(refxy[qt[i],j]);
                }
                for(j=0; j<=ny-1; j++)
                {
                    result = result | (double)(qxy[i,nx+j])!=(double)(refxy[qt[i],nx+j]);
                }
            }
            return result;
        }


        /*************************************************************************
        Returns norm
        *************************************************************************/
        private static double vnorm(ref double[] x,
            int n,
            int normtype)
        {
            double result = 0;
            int i = 0;

            result = AP.Math.RandomReal();
            if( normtype==0 )
            {
                result = 0;
                for(i=0; i<=n-1; i++)
                {
                    result = Math.Max(result, Math.Abs(x[i]));
                }
                return result;
            }
            if( normtype==1 )
            {
                result = 0;
                for(i=0; i<=n-1; i++)
                {
                    result = result+Math.Abs(x[i]);
                }
                return result;
            }
            if( normtype==2 )
            {
                result = 0;
                for(i=0; i<=n-1; i++)
                {
                    result = result+AP.Math.Sqr(x[i]);
                }
                result = Math.Sqrt(result);
                return result;
            }
            return result;
        }


        /*************************************************************************
        Testing Nearest Neighbor Search on uniformly distributed hypercube

        NormType: 0, 1, 2
        D: space dimension
        N: points count
        *************************************************************************/
        private static void testkdtuniform(ref double[,] xy,
            int n,
            int nx,
            int ny,
            int normtype,
            ref bool kdterrors)
        {
            double errtol = 0;
            int[] tags = new int[0];
            double[] ptx = new double[0];
            double[] tmpx = new double[0];
            bool[] tmpb = new bool[0];
            nearestneighbor.kdtree treex = new nearestneighbor.kdtree();
            nearestneighbor.kdtree treexy = new nearestneighbor.kdtree();
            nearestneighbor.kdtree treext = new nearestneighbor.kdtree();
            double[,] qx = new double[0,0];
            double[,] qxy = new double[0,0];
            int[] qtags = new int[0];
            double[] qr = new double[0];
            int kx = 0;
            int kxy = 0;
            int kt = 0;
            int kr = 0;
            double eps = 0;
            int i = 0;
            int j = 0;
            int k = 0;
            int task = 0;
            bool isequal = new bool();
            double r = 0;
            int q = 0;
            int qcount = 0;
            int i_ = 0;

            qcount = 10;
            
            //
            // Tol - roundoff error tolerance (for '>=' comparisons)
            //
            errtol = 100000*AP.Math.MachineEpsilon;
            
            //
            // fill tags
            //
            tags = new int[n];
            for(i=0; i<=n-1; i++)
            {
                tags[i] = i;
            }
            
            //
            // build trees
            //
            nearestneighbor.kdtreebuild(ref xy, n, nx, 0, normtype, ref treex);
            nearestneighbor.kdtreebuild(ref xy, n, nx, ny, normtype, ref treexy);
            nearestneighbor.kdtreebuildtagged(ref xy, ref tags, n, nx, 0, normtype, ref treext);
            
            //
            // allocate arrays
            //
            tmpx = new double[nx];
            tmpb = new bool[n];
            qx = new double[n, nx];
            qxy = new double[n, nx+ny];
            qtags = new int[n];
            qr = new double[n];
            ptx = new double[nx];
            
            //
            // test general K-NN queries (with self-matches):
            // * compare results from different trees (must be equal) and
            //   check that correct (value,tag) pairs are returned
            // * test results from XT tree - let R be radius of query result.
            //   then all points not in result must be not closer than R.
            //
            for(q=1; q<=qcount; q++)
            {
                
                //
                // Select K: 1..N
                //
                if( (double)(AP.Math.RandomReal())>(double)(0.5) )
                {
                    k = 1+AP.Math.RandomInteger(n);
                }
                else
                {
                    k = 1;
                }
                
                //
                // Select point (either one of the points, or random)
                //
                if( (double)(AP.Math.RandomReal())>(double)(0.5) )
                {
                    i = AP.Math.RandomInteger(n);
                    for(i_=0; i_<=nx-1;i_++)
                    {
                        ptx[i_] = xy[i,i_];
                    }
                }
                else
                {
                    for(i=0; i<=nx-1; i++)
                    {
                        ptx[i] = 2*AP.Math.RandomReal()-1;
                    }
                }
                
                //
                // Test:
                // * consistency of results from different queries
                // * points in query are IN the R-sphere (or at the boundary),
                //   and points not in query are outside of the R-sphere (or at the boundary)
                // * distances are correct and are ordered
                //
                kx = nearestneighbor.kdtreequeryknn(ref treex, ref ptx, k, true);
                kxy = nearestneighbor.kdtreequeryknn(ref treexy, ref ptx, k, true);
                kt = nearestneighbor.kdtreequeryknn(ref treext, ref ptx, k, true);
                if( kx!=k | kxy!=k | kt!=k )
                {
                    kdterrors = true;
                    return;
                }
                kx = 0;
                kxy = 0;
                kt = 0;
                nearestneighbor.kdtreequeryresultsx(ref treex, ref qx, ref kx);
                nearestneighbor.kdtreequeryresultsxy(ref treexy, ref qxy, ref kxy);
                nearestneighbor.kdtreequeryresultstags(ref treext, ref qtags, ref kt);
                nearestneighbor.kdtreequeryresultsdistances(ref treext, ref qr, ref kr);
                if( kx!=k | kxy!=k | kt!=k | kr!=k )
                {
                    kdterrors = true;
                    return;
                }
                kdterrors = kdterrors | kdtresultsdifferent(ref xy, n, ref qx, ref qxy, ref qtags, k, nx, ny);
                for(i=0; i<=n-1; i++)
                {
                    tmpb[i] = true;
                }
                r = 0;
                for(i=0; i<=k-1; i++)
                {
                    tmpb[qtags[i]] = false;
                    for(i_=0; i_<=nx-1;i_++)
                    {
                        tmpx[i_] = ptx[i_];
                    }
                    for(i_=0; i_<=nx-1;i_++)
                    {
                        tmpx[i_] = tmpx[i_] - qx[i,i_];
                    }
                    r = Math.Max(r, vnorm(ref tmpx, nx, normtype));
                }
                for(i=0; i<=n-1; i++)
                {
                    if( tmpb[i] )
                    {
                        for(i_=0; i_<=nx-1;i_++)
                        {
                            tmpx[i_] = ptx[i_];
                        }
                        for(i_=0; i_<=nx-1;i_++)
                        {
                            tmpx[i_] = tmpx[i_] - xy[i,i_];
                        }
                        kdterrors = kdterrors | (double)(vnorm(ref tmpx, nx, normtype))<(double)(r*(1-errtol));
                    }
                }
                for(i=0; i<=k-2; i++)
                {
                    kdterrors = kdterrors | (double)(qr[i])>(double)(qr[i+1]);
                }
                for(i=0; i<=k-1; i++)
                {
                    for(i_=0; i_<=nx-1;i_++)
                    {
                        tmpx[i_] = ptx[i_];
                    }
                    for(i_=0; i_<=nx-1;i_++)
                    {
                        tmpx[i_] = tmpx[i_] - xy[qtags[i],i_];
                    }
                    kdterrors = kdterrors | (double)(Math.Abs(vnorm(ref tmpx, nx, normtype)-qr[i]))>(double)(errtol);
                }
            }
            
            //
            // test general approximate K-NN queries (with self-matches):
            // * compare results from different trees (must be equal) and
            //   check that correct (value,tag) pairs are returned
            // * test results from XT tree - let R be radius of query result.
            //   then all points not in result must be not closer than R/(1+Eps).
            //
            for(q=1; q<=qcount; q++)
            {
                
                //
                // Select K: 1..N
                //
                if( (double)(AP.Math.RandomReal())>(double)(0.5) )
                {
                    k = 1+AP.Math.RandomInteger(n);
                }
                else
                {
                    k = 1;
                }
                
                //
                // Select Eps
                //
                eps = 0.5+AP.Math.RandomReal();
                
                //
                // Select point (either one of the points, or random)
                //
                if( (double)(AP.Math.RandomReal())>(double)(0.5) )
                {
                    i = AP.Math.RandomInteger(n);
                    for(i_=0; i_<=nx-1;i_++)
                    {
                        ptx[i_] = xy[i,i_];
                    }
                }
                else
                {
                    for(i=0; i<=nx-1; i++)
                    {
                        ptx[i] = 2*AP.Math.RandomReal()-1;
                    }
                }
                
                //
                // Test:
                // * consistency of results from different queries
                // * points in query are IN the R-sphere (or at the boundary),
                //   and points not in query are outside of the R-sphere (or at the boundary)
                // * distances are correct and are ordered
                //
                kx = nearestneighbor.kdtreequeryaknn(ref treex, ref ptx, k, true, eps);
                kxy = nearestneighbor.kdtreequeryaknn(ref treexy, ref ptx, k, true, eps);
                kt = nearestneighbor.kdtreequeryaknn(ref treext, ref ptx, k, true, eps);
                if( kx!=k | kxy!=k | kt!=k )
                {
                    kdterrors = true;
                    return;
                }
                kx = 0;
                kxy = 0;
                kt = 0;
                nearestneighbor.kdtreequeryresultsx(ref treex, ref qx, ref kx);
                nearestneighbor.kdtreequeryresultsxy(ref treexy, ref qxy, ref kxy);
                nearestneighbor.kdtreequeryresultstags(ref treext, ref qtags, ref kt);
                nearestneighbor.kdtreequeryresultsdistances(ref treext, ref qr, ref kr);
                if( kx!=k | kxy!=k | kt!=k | kr!=k )
                {
                    kdterrors = true;
                    return;
                }
                kdterrors = kdterrors | kdtresultsdifferent(ref xy, n, ref qx, ref qxy, ref qtags, k, nx, ny);
                for(i=0; i<=n-1; i++)
                {
                    tmpb[i] = true;
                }
                r = 0;
                for(i=0; i<=k-1; i++)
                {
                    tmpb[qtags[i]] = false;
                    for(i_=0; i_<=nx-1;i_++)
                    {
                        tmpx[i_] = ptx[i_];
                    }
                    for(i_=0; i_<=nx-1;i_++)
                    {
                        tmpx[i_] = tmpx[i_] - qx[i,i_];
                    }
                    r = Math.Max(r, vnorm(ref tmpx, nx, normtype));
                }
                for(i=0; i<=n-1; i++)
                {
                    if( tmpb[i] )
                    {
                        for(i_=0; i_<=nx-1;i_++)
                        {
                            tmpx[i_] = ptx[i_];
                        }
                        for(i_=0; i_<=nx-1;i_++)
                        {
                            tmpx[i_] = tmpx[i_] - xy[i,i_];
                        }
                        kdterrors = kdterrors | (double)(vnorm(ref tmpx, nx, normtype))<(double)(r*(1-errtol)/(1+eps));
                    }
                }
                for(i=0; i<=k-2; i++)
                {
                    kdterrors = kdterrors | (double)(qr[i])>(double)(qr[i+1]);
                }
                for(i=0; i<=k-1; i++)
                {
                    for(i_=0; i_<=nx-1;i_++)
                    {
                        tmpx[i_] = ptx[i_];
                    }
                    for(i_=0; i_<=nx-1;i_++)
                    {
                        tmpx[i_] = tmpx[i_] - xy[qtags[i],i_];
                    }
                    kdterrors = kdterrors | (double)(Math.Abs(vnorm(ref tmpx, nx, normtype)-qr[i]))>(double)(errtol);
                }
            }
            
            //
            // test general R-NN queries  (with self-matches):
            // * compare results from different trees (must be equal) and
            //   check that correct (value,tag) pairs are returned
            // * test results from XT tree - let R be radius of query result.
            //   then all points not in result must be not closer than R.
            //
            for(q=1; q<=qcount; q++)
            {
                
                //
                // Select R
                //
                if( (double)(AP.Math.RandomReal())>(double)(0.3) )
                {
                    r = Math.Max(AP.Math.RandomReal(), AP.Math.MachineEpsilon);
                }
                else
                {
                    r = AP.Math.MachineEpsilon;
                }
                
                //
                // Select point (either one of the points, or random)
                //
                if( (double)(AP.Math.RandomReal())>(double)(0.5) )
                {
                    i = AP.Math.RandomInteger(n);
                    for(i_=0; i_<=nx-1;i_++)
                    {
                        ptx[i_] = xy[i,i_];
                    }
                }
                else
                {
                    for(i=0; i<=nx-1; i++)
                    {
                        ptx[i] = 2*AP.Math.RandomReal()-1;
                    }
                }
                
                //
                // Test:
                // * consistency of results from different queries
                // * points in query are IN the R-sphere (or at the boundary),
                //   and points not in query are outside of the R-sphere (or at the boundary)
                // * distances are correct and are ordered
                //
                kx = nearestneighbor.kdtreequeryrnn(ref treex, ref ptx, r, true);
                kxy = nearestneighbor.kdtreequeryrnn(ref treexy, ref ptx, r, true);
                kt = nearestneighbor.kdtreequeryrnn(ref treext, ref ptx, r, true);
                if( kxy!=kx | kt!=kx )
                {
                    kdterrors = true;
                    return;
                }
                kx = 0;
                kxy = 0;
                kt = 0;
                nearestneighbor.kdtreequeryresultsx(ref treex, ref qx, ref kx);
                nearestneighbor.kdtreequeryresultsxy(ref treexy, ref qxy, ref kxy);
                nearestneighbor.kdtreequeryresultstags(ref treext, ref qtags, ref kt);
                nearestneighbor.kdtreequeryresultsdistances(ref treext, ref qr, ref kr);
                if( kxy!=kx | kt!=kx | kr!=kx )
                {
                    kdterrors = true;
                    return;
                }
                kdterrors = kdterrors | kdtresultsdifferent(ref xy, n, ref qx, ref qxy, ref qtags, kx, nx, ny);
                for(i=0; i<=n-1; i++)
                {
                    tmpb[i] = true;
                }
                for(i=0; i<=kx-1; i++)
                {
                    tmpb[qtags[i]] = false;
                }
                for(i=0; i<=n-1; i++)
                {
                    for(i_=0; i_<=nx-1;i_++)
                    {
                        tmpx[i_] = ptx[i_];
                    }
                    for(i_=0; i_<=nx-1;i_++)
                    {
                        tmpx[i_] = tmpx[i_] - xy[i,i_];
                    }
                    if( tmpb[i] )
                    {
                        kdterrors = kdterrors | (double)(vnorm(ref tmpx, nx, normtype))<(double)(r*(1-errtol));
                    }
                    else
                    {
                        kdterrors = kdterrors | (double)(vnorm(ref tmpx, nx, normtype))>(double)(r*(1+errtol));
                    }
                }
                for(i=0; i<=kx-2; i++)
                {
                    kdterrors = kdterrors | (double)(qr[i])>(double)(qr[i+1]);
                }
            }
            
            //
            // Test self-matching:
            // * self-match - nearest neighbor of each point in XY is the point itself
            // * no self-match - nearest neighbor is NOT the point itself
            //
            if( n>1 )
            {
                
                //
                // test for N=1 have non-general form, but it is not really needed
                //
                for(task=0; task<=1; task++)
                {
                    for(i=0; i<=n-1; i++)
                    {
                        for(i_=0; i_<=nx-1;i_++)
                        {
                            ptx[i_] = xy[i,i_];
                        }
                        kx = nearestneighbor.kdtreequeryknn(ref treex, ref ptx, 1, task==0);
                        nearestneighbor.kdtreequeryresultsx(ref treex, ref qx, ref kx);
                        if( kx!=1 )
                        {
                            kdterrors = true;
                            return;
                        }
                        isequal = true;
                        for(j=0; j<=nx-1; j++)
                        {
                            isequal = isequal & (double)(qx[0,j])==(double)(ptx[j]);
                        }
                        if( task==0 )
                        {
                            kdterrors = kdterrors | !isequal;
                        }
                        else
                        {
                            kdterrors = kdterrors | isequal;
                        }
                    }
                }
            }
        }


        /*************************************************************************
        Silent unit test
        *************************************************************************/
        public static bool testnearestneighborunit_test_silent()
        {
            bool result = new bool();

            result = testnearestneighbor(true);
            return result;
        }


        /*************************************************************************
        Unit test
        *************************************************************************/
        public static bool testnearestneighborunit_test()
        {
            bool result = new bool();

            result = testnearestneighbor(false);
            return result;
        }
    }
}
