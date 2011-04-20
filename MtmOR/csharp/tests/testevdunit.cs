
using System;

namespace alglib
{
    public class testevdunit
    {
        /*************************************************************************
        Testing symmetric EVD subroutine
        *************************************************************************/
        public static bool testevd(bool silent)
        {
            bool result = new bool();
            double[,] ra = new double[0,0];
            int n = 0;
            int j = 0;
            int failc = 0;
            int runs = 0;
            double failr = 0;
            double failthreshold = 0;
            double threshold = 0;
            double bithreshold = 0;
            bool waserrors = new bool();
            bool nserrors = new bool();
            bool serrors = new bool();
            bool herrors = new bool();
            bool tderrors = new bool();
            bool sbierrors = new bool();
            bool hbierrors = new bool();
            bool tdbierrors = new bool();
            bool wfailed = new bool();

            failthreshold = 0.005;
            threshold = 100000*AP.Math.MachineEpsilon;
            bithreshold = 1.0E-6;
            nserrors = false;
            serrors = false;
            herrors = false;
            tderrors = false;
            sbierrors = false;
            hbierrors = false;
            tdbierrors = false;
            failc = 0;
            runs = 0;
            
            //
            // Test problems
            //
            for(n=1; n<=ablas.ablasblocksize(ref ra); n++)
            {
                testevdset(n, threshold, bithreshold, ref failc, ref runs, ref nserrors, ref serrors, ref herrors, ref tderrors, ref sbierrors, ref hbierrors, ref tdbierrors);
            }
            for(j=2; j<=3; j++)
            {
                for(n=j*ablas.ablasblocksize(ref ra)-1; n<=j*ablas.ablasblocksize(ref ra)+1; n++)
                {
                    testevdset(n, threshold, bithreshold, ref failc, ref runs, ref nserrors, ref serrors, ref herrors, ref tderrors, ref sbierrors, ref hbierrors, ref tdbierrors);
                }
            }
            
            //
            // report
            //
            wfailed = (double)((double)(failc)/(double)(runs))>(double)(failthreshold);
            waserrors = nserrors | serrors | herrors | tderrors | sbierrors | hbierrors | tdbierrors | wfailed;
            if( !silent )
            {
                System.Console.Write("TESTING EVD UNIT");
                System.Console.WriteLine();
                System.Console.Write("NS ERRORS:                               ");
                if( !nserrors )
                {
                    System.Console.Write("OK");
                    System.Console.WriteLine();
                }
                else
                {
                    System.Console.Write("FAILED");
                    System.Console.WriteLine();
                }
                System.Console.Write("S ERRORS:                                ");
                if( !serrors )
                {
                    System.Console.Write("OK");
                    System.Console.WriteLine();
                }
                else
                {
                    System.Console.Write("FAILED");
                    System.Console.WriteLine();
                }
                System.Console.Write("H ERRORS:                                ");
                if( !herrors )
                {
                    System.Console.Write("OK");
                    System.Console.WriteLine();
                }
                else
                {
                    System.Console.Write("FAILED");
                    System.Console.WriteLine();
                }
                System.Console.Write("TD ERRORS:                               ");
                if( !tderrors )
                {
                    System.Console.Write("OK");
                    System.Console.WriteLine();
                }
                else
                {
                    System.Console.Write("FAILED");
                    System.Console.WriteLine();
                }
                System.Console.Write("SBI ERRORS:                              ");
                if( !sbierrors )
                {
                    System.Console.Write("OK");
                    System.Console.WriteLine();
                }
                else
                {
                    System.Console.Write("FAILED");
                    System.Console.WriteLine();
                }
                System.Console.Write("HBI ERRORS:                              ");
                if( !hbierrors )
                {
                    System.Console.Write("OK");
                    System.Console.WriteLine();
                }
                else
                {
                    System.Console.Write("FAILED");
                    System.Console.WriteLine();
                }
                System.Console.Write("TDBI ERRORS:                             ");
                if( !tdbierrors )
                {
                    System.Console.Write("OK");
                    System.Console.WriteLine();
                }
                else
                {
                    System.Console.Write("FAILED");
                    System.Console.WriteLine();
                }
                System.Console.Write("FAILURE THRESHOLD:                       ");
                if( !wfailed )
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
        Sparse fill
        *************************************************************************/
        private static void rmatrixfillsparsea(ref double[,] a,
            int m,
            int n,
            double sparcity)
        {
            int i = 0;
            int j = 0;

            for(i=0; i<=m-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    if( (double)(AP.Math.RandomReal())>=(double)(sparcity) )
                    {
                        a[i,j] = 2*AP.Math.RandomReal()-1;
                    }
                    else
                    {
                        a[i,j] = 0;
                    }
                }
            }
        }


        /*************************************************************************
        Sparse fill
        *************************************************************************/
        private static void cmatrixfillsparsea(ref AP.Complex[,] a,
            int m,
            int n,
            double sparcity)
        {
            int i = 0;
            int j = 0;

            for(i=0; i<=m-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    if( (double)(AP.Math.RandomReal())>=(double)(sparcity) )
                    {
                        a[i,j].x = 2*AP.Math.RandomReal()-1;
                        a[i,j].y = 2*AP.Math.RandomReal()-1;
                    }
                    else
                    {
                        a[i,j] = 0;
                    }
                }
            }
        }


        /*************************************************************************
        Copies A to AL (lower half) and AU (upper half), filling unused parts by
        random garbage.
        *************************************************************************/
        private static void rmatrixsymmetricsplit(ref double[,] a,
            int n,
            ref double[,] al,
            ref double[,] au)
        {
            int i = 0;
            int j = 0;

            for(i=0; i<=n-1; i++)
            {
                for(j=i+1; j<=n-1; j++)
                {
                    al[i,j] = 2*AP.Math.RandomReal()-1;
                    al[j,i] = a[i,j];
                    au[i,j] = a[i,j];
                    au[j,i] = 2*AP.Math.RandomReal()-1;
                }
                al[i,i] = a[i,i];
                au[i,i] = a[i,i];
            }
        }


        /*************************************************************************
        Copies A to AL (lower half) and AU (upper half), filling unused parts by
        random garbage.
        *************************************************************************/
        private static void cmatrixhermitiansplit(ref AP.Complex[,] a,
            int n,
            ref AP.Complex[,] al,
            ref AP.Complex[,] au)
        {
            int i = 0;
            int j = 0;

            for(i=0; i<=n-1; i++)
            {
                for(j=i+1; j<=n-1; j++)
                {
                    al[i,j] = 2*AP.Math.RandomReal()-1;
                    al[j,i] = AP.Math.Conj(a[i,j]);
                    au[i,j] = a[i,j];
                    au[j,i] = 2*AP.Math.RandomReal()-1;
                }
                al[i,i] = a[i,i];
                au[i,i] = a[i,i];
            }
        }


        /*************************************************************************
        Unsets 2D array.
        *************************************************************************/
        private static void unset2d(ref double[,] a)
        {
            a = new double[0+1, 0+1];
            a[0,0] = 2*AP.Math.RandomReal()-1;
        }


        /*************************************************************************
        Unsets 2D array.
        *************************************************************************/
        private static void cunset2d(ref AP.Complex[,] a)
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
        Unsets 1D array.
        *************************************************************************/
        private static void cunset1d(ref AP.Complex[] a)
        {
            a = new AP.Complex[0+1];
            a[0] = 2*AP.Math.RandomReal()-1;
        }


        /*************************************************************************
        Tests Z*Lambda*Z' against tridiag(D,E).
        Returns relative error.
        *************************************************************************/
        private static double tdtestproduct(ref double[] d,
            ref double[] e,
            int n,
            ref double[,] z,
            ref double[] lambda)
        {
            double result = 0;
            int i = 0;
            int j = 0;
            int k = 0;
            double v = 0;
            double mx = 0;

            result = 0;
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    
                    //
                    // Calculate V = A[i,j], A = Z*Lambda*Z'
                    //
                    v = 0;
                    for(k=0; k<=n-1; k++)
                    {
                        v = v+z[i,k]*lambda[k]*z[j,k];
                    }
                    
                    //
                    // Compare
                    //
                    if( Math.Abs(i-j)==0 )
                    {
                        result = Math.Max(result, Math.Abs(v-d[i]));
                    }
                    if( Math.Abs(i-j)==1 )
                    {
                        result = Math.Max(result, Math.Abs(v-e[Math.Min(i, j)]));
                    }
                    if( Math.Abs(i-j)>1 )
                    {
                        result = Math.Max(result, Math.Abs(v));
                    }
                }
            }
            mx = 0;
            for(i=0; i<=n-1; i++)
            {
                mx = Math.Max(mx, Math.Abs(d[i]));
            }
            for(i=0; i<=n-2; i++)
            {
                mx = Math.Max(mx, Math.Abs(e[i]));
            }
            if( (double)(mx)==(double)(0) )
            {
                mx = 1;
            }
            result = result/mx;
            return result;
        }


        /*************************************************************************
        Tests Z*Lambda*Z' against A
        Returns relative error.
        *************************************************************************/
        private static double testproduct(ref double[,] a,
            int n,
            ref double[,] z,
            ref double[] lambda)
        {
            double result = 0;
            int i = 0;
            int j = 0;
            int k = 0;
            double v = 0;
            double mx = 0;

            result = 0;
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    
                    //
                    // Calculate V = A[i,j], A = Z*Lambda*Z'
                    //
                    v = 0;
                    for(k=0; k<=n-1; k++)
                    {
                        v = v+z[i,k]*lambda[k]*z[j,k];
                    }
                    
                    //
                    // Compare
                    //
                    result = Math.Max(result, Math.Abs(v-a[i,j]));
                }
            }
            mx = 0;
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    mx = Math.Max(mx, Math.Abs(a[i,j]));
                }
            }
            if( (double)(mx)==(double)(0) )
            {
                mx = 1;
            }
            result = result/mx;
            return result;
        }


        /*************************************************************************
        Tests Z*Z' against diag(1...1)
        Returns absolute error.
        *************************************************************************/
        private static double testort(ref double[,] z,
            int n)
        {
            double result = 0;
            int i = 0;
            int j = 0;
            double v = 0;
            int i_ = 0;

            result = 0;
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    v = 0.0;
                    for(i_=0; i_<=n-1;i_++)
                    {
                        v += z[i_,i]*z[i_,j];
                    }
                    if( i==j )
                    {
                        v = v-1;
                    }
                    result = Math.Max(result, Math.Abs(v));
                }
            }
            return result;
        }


        /*************************************************************************
        Tests Z*Lambda*Z' against A
        Returns relative error.
        *************************************************************************/
        private static double testcproduct(ref AP.Complex[,] a,
            int n,
            ref AP.Complex[,] z,
            ref double[] lambda)
        {
            double result = 0;
            int i = 0;
            int j = 0;
            int k = 0;
            AP.Complex v = 0;
            double mx = 0;

            result = 0;
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    
                    //
                    // Calculate V = A[i,j], A = Z*Lambda*Z'
                    //
                    v = 0;
                    for(k=0; k<=n-1; k++)
                    {
                        v = v+z[i,k]*lambda[k]*AP.Math.Conj(z[j,k]);
                    }
                    
                    //
                    // Compare
                    //
                    result = Math.Max(result, AP.Math.AbsComplex(v-a[i,j]));
                }
            }
            mx = 0;
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    mx = Math.Max(mx, AP.Math.AbsComplex(a[i,j]));
                }
            }
            if( (double)(mx)==(double)(0) )
            {
                mx = 1;
            }
            result = result/mx;
            return result;
        }


        /*************************************************************************
        Tests Z*Z' against diag(1...1)
        Returns absolute error.
        *************************************************************************/
        private static double testcort(ref AP.Complex[,] z,
            int n)
        {
            double result = 0;
            int i = 0;
            int j = 0;
            AP.Complex v = 0;
            int i_ = 0;

            result = 0;
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    v = 0.0;
                    for(i_=0; i_<=n-1;i_++)
                    {
                        v += z[i_,i]*AP.Math.Conj(z[i_,j]);
                    }
                    if( i==j )
                    {
                        v = v-1;
                    }
                    result = Math.Max(result, AP.Math.AbsComplex(v));
                }
            }
            return result;
        }


        /*************************************************************************
        Tests SEVD problem
        *************************************************************************/
        private static void testsevdproblem(ref double[,] a,
            ref double[,] al,
            ref double[,] au,
            int n,
            double threshold,
            ref bool serrors,
            ref int failc,
            ref int runs)
        {
            double[] lambda = new double[0];
            double[] lambdaref = new double[0];
            double[,] z = new double[0,0];
            int i = 0;
            int j = 0;
            double v = 0;

            
            //
            // Test simple EVD: values and full vectors, lower A
            //
            unset1d(ref lambdaref);
            unset2d(ref z);
            runs = runs+1;
            if( !evd.smatrixevd(al, n, 1, false, ref lambdaref, ref z) )
            {
                failc = failc+1;
                return;
            }
            serrors = serrors | (double)(testproduct(ref a, n, ref z, ref lambdaref))>(double)(threshold);
            serrors = serrors | (double)(testort(ref z, n))>(double)(threshold);
            for(i=0; i<=n-2; i++)
            {
                if( (double)(lambdaref[i+1])<(double)(lambdaref[i]) )
                {
                    serrors = true;
                    return;
                }
            }
            
            //
            // Test simple EVD: values and full vectors, upper A
            //
            unset1d(ref lambda);
            unset2d(ref z);
            runs = runs+1;
            if( !evd.smatrixevd(au, n, 1, true, ref lambda, ref z) )
            {
                failc = failc+1;
                return;
            }
            serrors = serrors | (double)(testproduct(ref a, n, ref z, ref lambda))>(double)(threshold);
            serrors = serrors | (double)(testort(ref z, n))>(double)(threshold);
            for(i=0; i<=n-2; i++)
            {
                if( (double)(lambda[i+1])<(double)(lambda[i]) )
                {
                    serrors = true;
                    return;
                }
            }
            
            //
            // Test simple EVD: values only, lower A
            //
            unset1d(ref lambda);
            unset2d(ref z);
            runs = runs+1;
            if( !evd.smatrixevd(al, n, 0, false, ref lambda, ref z) )
            {
                failc = failc+1;
                return;
            }
            for(i=0; i<=n-1; i++)
            {
                serrors = serrors | (double)(Math.Abs(lambda[i]-lambdaref[i]))>(double)(threshold);
            }
            
            //
            // Test simple EVD: values only, upper A
            //
            unset1d(ref lambda);
            unset2d(ref z);
            runs = runs+1;
            if( !evd.smatrixevd(au, n, 0, true, ref lambda, ref z) )
            {
                failc = failc+1;
                return;
            }
            for(i=0; i<=n-1; i++)
            {
                serrors = serrors | (double)(Math.Abs(lambda[i]-lambdaref[i]))>(double)(threshold);
            }
        }


        /*************************************************************************
        Tests SEVD problem
        *************************************************************************/
        private static void testhevdproblem(ref AP.Complex[,] a,
            ref AP.Complex[,] al,
            ref AP.Complex[,] au,
            int n,
            double threshold,
            ref bool herrors,
            ref int failc,
            ref int runs)
        {
            double[] lambda = new double[0];
            double[] lambdaref = new double[0];
            AP.Complex[,] z = new AP.Complex[0,0];
            int i = 0;
            int j = 0;
            AP.Complex v = 0;

            
            //
            // Test simple EVD: values and full vectors, lower A
            //
            unset1d(ref lambdaref);
            cunset2d(ref z);
            runs = runs+1;
            if( !evd.hmatrixevd(al, n, 1, false, ref lambdaref, ref z) )
            {
                failc = failc+1;
                return;
            }
            herrors = herrors | (double)(testcproduct(ref a, n, ref z, ref lambdaref))>(double)(threshold);
            herrors = herrors | (double)(testcort(ref z, n))>(double)(threshold);
            for(i=0; i<=n-2; i++)
            {
                if( (double)(lambdaref[i+1])<(double)(lambdaref[i]) )
                {
                    herrors = true;
                    return;
                }
            }
            
            //
            // Test simple EVD: values and full vectors, upper A
            //
            unset1d(ref lambda);
            cunset2d(ref z);
            runs = runs+1;
            if( !evd.hmatrixevd(au, n, 1, true, ref lambda, ref z) )
            {
                failc = failc+1;
                return;
            }
            herrors = herrors | (double)(testcproduct(ref a, n, ref z, ref lambda))>(double)(threshold);
            herrors = herrors | (double)(testcort(ref z, n))>(double)(threshold);
            for(i=0; i<=n-2; i++)
            {
                if( (double)(lambda[i+1])<(double)(lambda[i]) )
                {
                    herrors = true;
                    return;
                }
            }
            
            //
            // Test simple EVD: values only, lower A
            //
            unset1d(ref lambda);
            cunset2d(ref z);
            runs = runs+1;
            if( !evd.hmatrixevd(al, n, 0, false, ref lambda, ref z) )
            {
                failc = failc+1;
                return;
            }
            for(i=0; i<=n-1; i++)
            {
                herrors = herrors | (double)(Math.Abs(lambda[i]-lambdaref[i]))>(double)(threshold);
            }
            
            //
            // Test simple EVD: values only, upper A
            //
            unset1d(ref lambda);
            cunset2d(ref z);
            runs = runs+1;
            if( !evd.hmatrixevd(au, n, 0, true, ref lambda, ref z) )
            {
                failc = failc+1;
                return;
            }
            for(i=0; i<=n-1; i++)
            {
                herrors = herrors | (double)(Math.Abs(lambda[i]-lambdaref[i]))>(double)(threshold);
            }
        }


        /*************************************************************************
        Tests EVD problem

        DistVals    -   is True, when eigenvalues are distinct. Is False, when we
                        are solving sparse task with  lots  of  zero  eigenvalues.
                        In such cases some tests related to the  eigenvectors  are
                        not performed.
        *************************************************************************/
        private static void testsevdbiproblem(ref double[,] afull,
            ref double[,] al,
            ref double[,] au,
            int n,
            bool distvals,
            double threshold,
            ref bool serrors,
            ref int failc,
            ref int runs)
        {
            double[] lambda = new double[0];
            double[] lambdaref = new double[0];
            double[,] z = new double[0,0];
            double[,] zref = new double[0,0];
            double[,] a1 = new double[0,0];
            double[,] a2 = new double[0,0];
            double[,] ar = new double[0,0];
            bool wsucc = new bool();
            int i = 0;
            int j = 0;
            int k = 0;
            int m = 0;
            int i1 = 0;
            int i2 = 0;
            double v = 0;
            double a = 0;
            double b = 0;
            int i_ = 0;

            lambdaref = new double[n-1+1];
            zref = new double[n-1+1, n-1+1];
            a1 = new double[n-1+1, n-1+1];
            a2 = new double[n-1+1, n-1+1];
            
            //
            // Reference EVD
            //
            runs = runs+1;
            if( !evd.smatrixevd(afull, n, 1, true, ref lambdaref, ref zref) )
            {
                failc = failc+1;
                return;
            }
            
            //
            // Select random interval boundaries.
            // If there are non-distinct eigenvalues at the boundaries,
            // we move indexes further until values splits. It is done to
            // avoid situations where we can't get definite answer.
            //
            i1 = AP.Math.RandomInteger(n);
            i2 = i1+AP.Math.RandomInteger(n-i1);
            while( i1>0 )
            {
                if( (double)(Math.Abs(lambdaref[i1-1]-lambdaref[i1]))>(double)(10*threshold) )
                {
                    break;
                }
                i1 = i1-1;
            }
            while( i2<n-1 )
            {
                if( (double)(Math.Abs(lambdaref[i2+1]-lambdaref[i2]))>(double)(10*threshold) )
                {
                    break;
                }
                i2 = i2+1;
            }
            
            //
            // Select A, B
            //
            if( i1>0 )
            {
                a = 0.5*(lambdaref[i1]+lambdaref[i1-1]);
            }
            else
            {
                a = lambdaref[0]-1;
            }
            if( i2<n-1 )
            {
                b = 0.5*(lambdaref[i2]+lambdaref[i2+1]);
            }
            else
            {
                b = lambdaref[n-1]+1;
            }
            
            //
            // Test interval, no vectors, lower A
            //
            unset1d(ref lambda);
            unset2d(ref z);
            runs = runs+1;
            if( !evd.smatrixevdr(al, n, 0, false, a, b, ref m, ref lambda, ref z) )
            {
                failc = failc+1;
                return;
            }
            if( m!=i2-i1+1 )
            {
                failc = failc+1;
                return;
            }
            for(k=0; k<=m-1; k++)
            {
                serrors = serrors | (double)(Math.Abs(lambda[k]-lambdaref[i1+k]))>(double)(threshold);
            }
            
            //
            // Test interval, no vectors, upper A
            //
            unset1d(ref lambda);
            unset2d(ref z);
            runs = runs+1;
            if( !evd.smatrixevdr(au, n, 0, true, a, b, ref m, ref lambda, ref z) )
            {
                failc = failc+1;
                return;
            }
            if( m!=i2-i1+1 )
            {
                failc = failc+1;
                return;
            }
            for(k=0; k<=m-1; k++)
            {
                serrors = serrors | (double)(Math.Abs(lambda[k]-lambdaref[i1+k]))>(double)(threshold);
            }
            
            //
            // Test indexes, no vectors, lower A
            //
            unset1d(ref lambda);
            unset2d(ref z);
            runs = runs+1;
            if( !evd.smatrixevdi(al, n, 0, false, i1, i2, ref lambda, ref z) )
            {
                failc = failc+1;
                return;
            }
            m = i2-i1+1;
            for(k=0; k<=m-1; k++)
            {
                serrors = serrors | (double)(Math.Abs(lambda[k]-lambdaref[i1+k]))>(double)(threshold);
            }
            
            //
            // Test indexes, no vectors, upper A
            //
            unset1d(ref lambda);
            unset2d(ref z);
            runs = runs+1;
            if( !evd.smatrixevdi(au, n, 0, true, i1, i2, ref lambda, ref z) )
            {
                failc = failc+1;
                return;
            }
            m = i2-i1+1;
            for(k=0; k<=m-1; k++)
            {
                serrors = serrors | (double)(Math.Abs(lambda[k]-lambdaref[i1+k]))>(double)(threshold);
            }
            
            //
            // Test interval, vectors, lower A
            //
            unset1d(ref lambda);
            unset2d(ref z);
            runs = runs+1;
            if( !evd.smatrixevdr(al, n, 1, false, a, b, ref m, ref lambda, ref z) )
            {
                failc = failc+1;
                return;
            }
            if( m!=i2-i1+1 )
            {
                failc = failc+1;
                return;
            }
            for(k=0; k<=m-1; k++)
            {
                serrors = serrors | (double)(Math.Abs(lambda[k]-lambdaref[i1+k]))>(double)(threshold);
            }
            if( distvals )
            {
                
                //
                // Distinct eigenvalues, test vectors
                //
                for(j=0; j<=m-1; j++)
                {
                    v = 0.0;
                    for(i_=0; i_<=n-1;i_++)
                    {
                        v += z[i_,j]*zref[i_,i1+j];
                    }
                    if( (double)(v)<(double)(0) )
                    {
                        for(i_=0; i_<=n-1;i_++)
                        {
                            z[i_,j] = -1*z[i_,j];
                        }
                    }
                }
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=m-1; j++)
                    {
                        serrors = serrors | (double)(Math.Abs(z[i,j]-zref[i,i1+j]))>(double)(threshold);
                    }
                }
            }
            
            //
            // Test interval, vectors, upper A
            //
            unset1d(ref lambda);
            unset2d(ref z);
            runs = runs+1;
            if( !evd.smatrixevdr(au, n, 1, true, a, b, ref m, ref lambda, ref z) )
            {
                failc = failc+1;
                return;
            }
            if( m!=i2-i1+1 )
            {
                failc = failc+1;
                return;
            }
            for(k=0; k<=m-1; k++)
            {
                serrors = serrors | (double)(Math.Abs(lambda[k]-lambdaref[i1+k]))>(double)(threshold);
            }
            if( distvals )
            {
                
                //
                // Distinct eigenvalues, test vectors
                //
                for(j=0; j<=m-1; j++)
                {
                    v = 0.0;
                    for(i_=0; i_<=n-1;i_++)
                    {
                        v += z[i_,j]*zref[i_,i1+j];
                    }
                    if( (double)(v)<(double)(0) )
                    {
                        for(i_=0; i_<=n-1;i_++)
                        {
                            z[i_,j] = -1*z[i_,j];
                        }
                    }
                }
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=m-1; j++)
                    {
                        serrors = serrors | (double)(Math.Abs(z[i,j]-zref[i,i1+j]))>(double)(threshold);
                    }
                }
            }
            
            //
            // Test indexes, vectors, lower A
            //
            unset1d(ref lambda);
            unset2d(ref z);
            runs = runs+1;
            if( !evd.smatrixevdi(al, n, 1, false, i1, i2, ref lambda, ref z) )
            {
                failc = failc+1;
                return;
            }
            m = i2-i1+1;
            for(k=0; k<=m-1; k++)
            {
                serrors = serrors | (double)(Math.Abs(lambda[k]-lambdaref[i1+k]))>(double)(threshold);
            }
            if( distvals )
            {
                
                //
                // Distinct eigenvalues, test vectors
                //
                for(j=0; j<=m-1; j++)
                {
                    v = 0.0;
                    for(i_=0; i_<=n-1;i_++)
                    {
                        v += z[i_,j]*zref[i_,i1+j];
                    }
                    if( (double)(v)<(double)(0) )
                    {
                        for(i_=0; i_<=n-1;i_++)
                        {
                            z[i_,j] = -1*z[i_,j];
                        }
                    }
                }
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=m-1; j++)
                    {
                        serrors = serrors | (double)(Math.Abs(z[i,j]-zref[i,i1+j]))>(double)(threshold);
                    }
                }
            }
            
            //
            // Test indexes, vectors, upper A
            //
            unset1d(ref lambda);
            unset2d(ref z);
            runs = runs+1;
            if( !evd.smatrixevdi(au, n, 1, true, i1, i2, ref lambda, ref z) )
            {
                failc = failc+1;
                return;
            }
            m = i2-i1+1;
            for(k=0; k<=m-1; k++)
            {
                serrors = serrors | (double)(Math.Abs(lambda[k]-lambdaref[i1+k]))>(double)(threshold);
            }
            if( distvals )
            {
                
                //
                // Distinct eigenvalues, test vectors
                //
                for(j=0; j<=m-1; j++)
                {
                    v = 0.0;
                    for(i_=0; i_<=n-1;i_++)
                    {
                        v += z[i_,j]*zref[i_,i1+j];
                    }
                    if( (double)(v)<(double)(0) )
                    {
                        for(i_=0; i_<=n-1;i_++)
                        {
                            z[i_,j] = -1*z[i_,j];
                        }
                    }
                }
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=m-1; j++)
                    {
                        serrors = serrors | (double)(Math.Abs(z[i,j]-zref[i,i1+j]))>(double)(threshold);
                    }
                }
            }
        }


        /*************************************************************************
        Tests EVD problem

        DistVals    -   is True, when eigenvalues are distinct. Is False, when we
                        are solving sparse task with  lots  of  zero  eigenvalues.
                        In such cases some tests related to the  eigenvectors  are
                        not performed.
        *************************************************************************/
        private static void testhevdbiproblem(ref AP.Complex[,] afull,
            ref AP.Complex[,] al,
            ref AP.Complex[,] au,
            int n,
            bool distvals,
            double threshold,
            ref bool herrors,
            ref int failc,
            ref int runs)
        {
            double[] lambda = new double[0];
            double[] lambdaref = new double[0];
            AP.Complex[,] z = new AP.Complex[0,0];
            AP.Complex[,] zref = new AP.Complex[0,0];
            AP.Complex[,] a1 = new AP.Complex[0,0];
            AP.Complex[,] a2 = new AP.Complex[0,0];
            AP.Complex[,] ar = new AP.Complex[0,0];
            bool wsucc = new bool();
            int i = 0;
            int j = 0;
            int k = 0;
            int m = 0;
            int i1 = 0;
            int i2 = 0;
            AP.Complex v = 0;
            double a = 0;
            double b = 0;
            int i_ = 0;

            lambdaref = new double[n-1+1];
            zref = new AP.Complex[n-1+1, n-1+1];
            a1 = new AP.Complex[n-1+1, n-1+1];
            a2 = new AP.Complex[n-1+1, n-1+1];
            
            //
            // Reference EVD
            //
            runs = runs+1;
            if( !evd.hmatrixevd(afull, n, 1, true, ref lambdaref, ref zref) )
            {
                failc = failc+1;
                return;
            }
            
            //
            // Select random interval boundaries.
            // If there are non-distinct eigenvalues at the boundaries,
            // we move indexes further until values splits. It is done to
            // avoid situations where we can't get definite answer.
            //
            i1 = AP.Math.RandomInteger(n);
            i2 = i1+AP.Math.RandomInteger(n-i1);
            while( i1>0 )
            {
                if( (double)(Math.Abs(lambdaref[i1-1]-lambdaref[i1]))>(double)(10*threshold) )
                {
                    break;
                }
                i1 = i1-1;
            }
            while( i2<n-1 )
            {
                if( (double)(Math.Abs(lambdaref[i2+1]-lambdaref[i2]))>(double)(10*threshold) )
                {
                    break;
                }
                i2 = i2+1;
            }
            
            //
            // Select A, B
            //
            if( i1>0 )
            {
                a = 0.5*(lambdaref[i1]+lambdaref[i1-1]);
            }
            else
            {
                a = lambdaref[0]-1;
            }
            if( i2<n-1 )
            {
                b = 0.5*(lambdaref[i2]+lambdaref[i2+1]);
            }
            else
            {
                b = lambdaref[n-1]+1;
            }
            
            //
            // Test interval, no vectors, lower A
            //
            unset1d(ref lambda);
            cunset2d(ref z);
            runs = runs+1;
            if( !evd.hmatrixevdr(al, n, 0, false, a, b, ref m, ref lambda, ref z) )
            {
                failc = failc+1;
                return;
            }
            if( m!=i2-i1+1 )
            {
                failc = failc+1;
                return;
            }
            for(k=0; k<=m-1; k++)
            {
                herrors = herrors | (double)(Math.Abs(lambda[k]-lambdaref[i1+k]))>(double)(threshold);
            }
            
            //
            // Test interval, no vectors, upper A
            //
            unset1d(ref lambda);
            cunset2d(ref z);
            runs = runs+1;
            if( !evd.hmatrixevdr(au, n, 0, true, a, b, ref m, ref lambda, ref z) )
            {
                failc = failc+1;
                return;
            }
            if( m!=i2-i1+1 )
            {
                failc = failc+1;
                return;
            }
            for(k=0; k<=m-1; k++)
            {
                herrors = herrors | (double)(Math.Abs(lambda[k]-lambdaref[i1+k]))>(double)(threshold);
            }
            
            //
            // Test indexes, no vectors, lower A
            //
            unset1d(ref lambda);
            cunset2d(ref z);
            runs = runs+1;
            if( !evd.hmatrixevdi(al, n, 0, false, i1, i2, ref lambda, ref z) )
            {
                failc = failc+1;
                return;
            }
            m = i2-i1+1;
            for(k=0; k<=m-1; k++)
            {
                herrors = herrors | (double)(Math.Abs(lambda[k]-lambdaref[i1+k]))>(double)(threshold);
            }
            
            //
            // Test indexes, no vectors, upper A
            //
            unset1d(ref lambda);
            cunset2d(ref z);
            runs = runs+1;
            if( !evd.hmatrixevdi(au, n, 0, true, i1, i2, ref lambda, ref z) )
            {
                failc = failc+1;
                return;
            }
            m = i2-i1+1;
            for(k=0; k<=m-1; k++)
            {
                herrors = herrors | (double)(Math.Abs(lambda[k]-lambdaref[i1+k]))>(double)(threshold);
            }
            
            //
            // Test interval, vectors, lower A
            //
            unset1d(ref lambda);
            cunset2d(ref z);
            runs = runs+1;
            if( !evd.hmatrixevdr(al, n, 1, false, a, b, ref m, ref lambda, ref z) )
            {
                failc = failc+1;
                return;
            }
            if( m!=i2-i1+1 )
            {
                failc = failc+1;
                return;
            }
            for(k=0; k<=m-1; k++)
            {
                herrors = herrors | (double)(Math.Abs(lambda[k]-lambdaref[i1+k]))>(double)(threshold);
            }
            if( distvals )
            {
                
                //
                // Distinct eigenvalues, test vectors
                //
                for(j=0; j<=m-1; j++)
                {
                    v = 0.0;
                    for(i_=0; i_<=n-1;i_++)
                    {
                        v += z[i_,j]*AP.Math.Conj(zref[i_,i1+j]);
                    }
                    v = AP.Math.Conj(v/AP.Math.AbsComplex(v));
                    for(i_=0; i_<=n-1;i_++)
                    {
                        z[i_,j] = v*z[i_,j];
                    }
                }
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=m-1; j++)
                    {
                        herrors = herrors | (double)(AP.Math.AbsComplex(z[i,j]-zref[i,i1+j]))>(double)(threshold);
                    }
                }
            }
            
            //
            // Test interval, vectors, upper A
            //
            unset1d(ref lambda);
            cunset2d(ref z);
            runs = runs+1;
            if( !evd.hmatrixevdr(au, n, 1, true, a, b, ref m, ref lambda, ref z) )
            {
                failc = failc+1;
                return;
            }
            if( m!=i2-i1+1 )
            {
                failc = failc+1;
                return;
            }
            for(k=0; k<=m-1; k++)
            {
                herrors = herrors | (double)(Math.Abs(lambda[k]-lambdaref[i1+k]))>(double)(threshold);
            }
            if( distvals )
            {
                
                //
                // Distinct eigenvalues, test vectors
                //
                for(j=0; j<=m-1; j++)
                {
                    v = 0.0;
                    for(i_=0; i_<=n-1;i_++)
                    {
                        v += z[i_,j]*AP.Math.Conj(zref[i_,i1+j]);
                    }
                    v = AP.Math.Conj(v/AP.Math.AbsComplex(v));
                    for(i_=0; i_<=n-1;i_++)
                    {
                        z[i_,j] = v*z[i_,j];
                    }
                }
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=m-1; j++)
                    {
                        herrors = herrors | (double)(AP.Math.AbsComplex(z[i,j]-zref[i,i1+j]))>(double)(threshold);
                    }
                }
            }
            
            //
            // Test indexes, vectors, lower A
            //
            unset1d(ref lambda);
            cunset2d(ref z);
            runs = runs+1;
            if( !evd.hmatrixevdi(al, n, 1, false, i1, i2, ref lambda, ref z) )
            {
                failc = failc+1;
                return;
            }
            m = i2-i1+1;
            for(k=0; k<=m-1; k++)
            {
                herrors = herrors | (double)(Math.Abs(lambda[k]-lambdaref[i1+k]))>(double)(threshold);
            }
            if( distvals )
            {
                
                //
                // Distinct eigenvalues, test vectors
                //
                for(j=0; j<=m-1; j++)
                {
                    v = 0.0;
                    for(i_=0; i_<=n-1;i_++)
                    {
                        v += z[i_,j]*AP.Math.Conj(zref[i_,i1+j]);
                    }
                    v = AP.Math.Conj(v/AP.Math.AbsComplex(v));
                    for(i_=0; i_<=n-1;i_++)
                    {
                        z[i_,j] = v*z[i_,j];
                    }
                }
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=m-1; j++)
                    {
                        herrors = herrors | (double)(AP.Math.AbsComplex(z[i,j]-zref[i,i1+j]))>(double)(threshold);
                    }
                }
            }
            
            //
            // Test indexes, vectors, upper A
            //
            unset1d(ref lambda);
            cunset2d(ref z);
            runs = runs+1;
            if( !evd.hmatrixevdi(au, n, 1, true, i1, i2, ref lambda, ref z) )
            {
                failc = failc+1;
                return;
            }
            m = i2-i1+1;
            for(k=0; k<=m-1; k++)
            {
                herrors = herrors | (double)(Math.Abs(lambda[k]-lambdaref[i1+k]))>(double)(threshold);
            }
            if( distvals )
            {
                
                //
                // Distinct eigenvalues, test vectors
                //
                for(j=0; j<=m-1; j++)
                {
                    v = 0.0;
                    for(i_=0; i_<=n-1;i_++)
                    {
                        v += z[i_,j]*AP.Math.Conj(zref[i_,i1+j]);
                    }
                    v = AP.Math.Conj(v/AP.Math.AbsComplex(v));
                    for(i_=0; i_<=n-1;i_++)
                    {
                        z[i_,j] = v*z[i_,j];
                    }
                }
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=m-1; j++)
                    {
                        herrors = herrors | (double)(AP.Math.AbsComplex(z[i,j]-zref[i,i1+j]))>(double)(threshold);
                    }
                }
            }
        }


        /*************************************************************************
        Tests EVD problem
        *************************************************************************/
        private static void testtdevdproblem(ref double[] d,
            ref double[] e,
            int n,
            double threshold,
            ref bool tderrors,
            ref int failc,
            ref int runs)
        {
            double[] lambda = new double[0];
            double[] ee = new double[0];
            double[] lambda2 = new double[0];
            double[,] z = new double[0,0];
            double[,] zref = new double[0,0];
            double[,] a1 = new double[0,0];
            double[,] a2 = new double[0,0];
            bool wsucc = new bool();
            int i = 0;
            int j = 0;
            double v = 0;
            int i_ = 0;

            lambda = new double[n-1+1];
            lambda2 = new double[n-1+1];
            zref = new double[n-1+1, n-1+1];
            a1 = new double[n-1+1, n-1+1];
            a2 = new double[n-1+1, n-1+1];
            if( n>1 )
            {
                ee = new double[n-2+1];
            }
            
            //
            // Test simple EVD: values and full vectors
            //
            for(i=0; i<=n-1; i++)
            {
                lambda[i] = d[i];
            }
            for(i=0; i<=n-2; i++)
            {
                ee[i] = e[i];
            }
            runs = runs+1;
            wsucc = evd.smatrixtdevd(ref lambda, ee, n, 2, ref z);
            if( !wsucc )
            {
                failc = failc+1;
                return;
            }
            tderrors = tderrors | (double)(tdtestproduct(ref d, ref e, n, ref z, ref lambda))>(double)(threshold);
            tderrors = tderrors | (double)(testort(ref z, n))>(double)(threshold);
            for(i=0; i<=n-2; i++)
            {
                if( (double)(lambda[i+1])<(double)(lambda[i]) )
                {
                    tderrors = true;
                    return;
                }
            }
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    zref[i,j] = z[i,j];
                }
            }
            
            //
            // Test values only variant
            //
            for(i=0; i<=n-1; i++)
            {
                lambda2[i] = d[i];
            }
            for(i=0; i<=n-2; i++)
            {
                ee[i] = e[i];
            }
            runs = runs+1;
            wsucc = evd.smatrixtdevd(ref lambda2, ee, n, 0, ref z);
            if( !wsucc )
            {
                failc = failc+1;
                return;
            }
            for(i=0; i<=n-1; i++)
            {
                tderrors = tderrors | (double)(Math.Abs(lambda2[i]-lambda[i]))>(double)(threshold);
            }
            
            //
            // Test multiplication variant
            //
            for(i=0; i<=n-1; i++)
            {
                lambda2[i] = d[i];
            }
            for(i=0; i<=n-2; i++)
            {
                ee[i] = e[i];
            }
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    a1[i,j] = 2*AP.Math.RandomReal()-1;
                    a2[i,j] = a1[i,j];
                }
            }
            runs = runs+1;
            wsucc = evd.smatrixtdevd(ref lambda2, ee, n, 1, ref a1);
            if( !wsucc )
            {
                failc = failc+1;
                return;
            }
            for(i=0; i<=n-1; i++)
            {
                tderrors = tderrors | (double)(Math.Abs(lambda2[i]-lambda[i]))>(double)(threshold);
            }
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    v = 0.0;
                    for(i_=0; i_<=n-1;i_++)
                    {
                        v += a2[i,i_]*zref[i_,j];
                    }
                    
                    //
                    // next line is a bit complicated because
                    // depending on algorithm used we can get either
                    // z or -z as eigenvector. so we compare result
                    // with both A*ZRef and -A*ZRef
                    //
                    tderrors = tderrors | (double)(Math.Abs(v-a1[i,j]))>(double)(threshold) & (double)(Math.Abs(v+a1[i,j]))>(double)(threshold);
                }
            }
            
            //
            // Test first row variant
            //
            for(i=0; i<=n-1; i++)
            {
                lambda2[i] = d[i];
            }
            for(i=0; i<=n-2; i++)
            {
                ee[i] = e[i];
            }
            runs = runs+1;
            wsucc = evd.smatrixtdevd(ref lambda2, ee, n, 3, ref z);
            if( !wsucc )
            {
                failc = failc+1;
                return;
            }
            for(i=0; i<=n-1; i++)
            {
                tderrors = tderrors | (double)(Math.Abs(lambda2[i]-lambda[i]))>(double)(threshold);
                
                //
                // next line is a bit complicated because
                // depending on algorithm used we can get either
                // z or -z as eigenvector. so we compare result
                // with both z and -z
                //
                tderrors = tderrors | (double)(Math.Abs(z[0,i]-zref[0,i]))>(double)(threshold) & (double)(Math.Abs(z[0,i]+zref[0,i]))>(double)(threshold);
            }
        }


        /*************************************************************************
        Tests EVD problem

        DistVals    -   is True, when eigenvalues are distinct. Is False, when we
                        are solving sparse task with  lots  of  zero  eigenvalues.
                        In such cases some tests related to the  eigenvectors  are
                        not performed.
        *************************************************************************/
        private static void testtdevdbiproblem(ref double[] d,
            ref double[] e,
            int n,
            bool distvals,
            double threshold,
            ref bool serrors,
            ref int failc,
            ref int runs)
        {
            double[] lambda = new double[0];
            double[] lambdaref = new double[0];
            double[,] z = new double[0,0];
            double[,] zref = new double[0,0];
            double[,] a1 = new double[0,0];
            double[,] a2 = new double[0,0];
            double[,] ar = new double[0,0];
            bool wsucc = new bool();
            int i = 0;
            int j = 0;
            int k = 0;
            int m = 0;
            int i1 = 0;
            int i2 = 0;
            double v = 0;
            double a = 0;
            double b = 0;
            int i_ = 0;

            lambdaref = new double[n-1+1];
            zref = new double[n-1+1, n-1+1];
            a1 = new double[n-1+1, n-1+1];
            a2 = new double[n-1+1, n-1+1];
            
            //
            // Reference EVD
            //
            lambdaref = new double[n];
            for(i_=0; i_<=n-1;i_++)
            {
                lambdaref[i_] = d[i_];
            }
            runs = runs+1;
            if( !evd.smatrixtdevd(ref lambdaref, e, n, 2, ref zref) )
            {
                failc = failc+1;
                return;
            }
            
            //
            // Select random interval boundaries.
            // If there are non-distinct eigenvalues at the boundaries,
            // we move indexes further until values splits. It is done to
            // avoid situations where we can't get definite answer.
            //
            i1 = AP.Math.RandomInteger(n);
            i2 = i1+AP.Math.RandomInteger(n-i1);
            while( i1>0 )
            {
                if( (double)(Math.Abs(lambdaref[i1-1]-lambdaref[i1]))>(double)(10*threshold) )
                {
                    break;
                }
                i1 = i1-1;
            }
            while( i2<n-1 )
            {
                if( (double)(Math.Abs(lambdaref[i2+1]-lambdaref[i2]))>(double)(10*threshold) )
                {
                    break;
                }
                i2 = i2+1;
            }
            
            //
            // Test different combinations
            //
            
            //
            // Select A, B
            //
            if( i1>0 )
            {
                a = 0.5*(lambdaref[i1]+lambdaref[i1-1]);
            }
            else
            {
                a = lambdaref[0]-1;
            }
            if( i2<n-1 )
            {
                b = 0.5*(lambdaref[i2]+lambdaref[i2+1]);
            }
            else
            {
                b = lambdaref[n-1]+1;
            }
            
            //
            // Test interval, no vectors
            //
            lambda = new double[n-1+1];
            for(i=0; i<=n-1; i++)
            {
                lambda[i] = d[i];
            }
            runs = runs+1;
            if( !evd.smatrixtdevdr(ref lambda, ref e, n, 0, a, b, ref m, ref z) )
            {
                failc = failc+1;
                return;
            }
            if( m!=i2-i1+1 )
            {
                failc = failc+1;
                return;
            }
            for(k=0; k<=m-1; k++)
            {
                serrors = serrors | (double)(Math.Abs(lambda[k]-lambdaref[i1+k]))>(double)(threshold);
            }
            
            //
            // Test indexes, no vectors
            //
            lambda = new double[n-1+1];
            for(i=0; i<=n-1; i++)
            {
                lambda[i] = d[i];
            }
            runs = runs+1;
            if( !evd.smatrixtdevdi(ref lambda, ref e, n, 0, i1, i2, ref z) )
            {
                failc = failc+1;
                return;
            }
            m = i2-i1+1;
            for(k=0; k<=m-1; k++)
            {
                serrors = serrors | (double)(Math.Abs(lambda[k]-lambdaref[i1+k]))>(double)(threshold);
            }
            
            //
            // Test interval, transform vectors
            //
            lambda = new double[n-1+1];
            for(i=0; i<=n-1; i++)
            {
                lambda[i] = d[i];
            }
            a1 = new double[n-1+1, n-1+1];
            a2 = new double[n-1+1, n-1+1];
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    a1[i,j] = 2*AP.Math.RandomReal()-1;
                    a2[i,j] = a1[i,j];
                }
            }
            runs = runs+1;
            if( !evd.smatrixtdevdr(ref lambda, ref e, n, 1, a, b, ref m, ref a1) )
            {
                failc = failc+1;
                return;
            }
            if( m!=i2-i1+1 )
            {
                failc = failc+1;
                return;
            }
            for(k=0; k<=m-1; k++)
            {
                serrors = serrors | (double)(Math.Abs(lambda[k]-lambdaref[i1+k]))>(double)(threshold);
            }
            if( distvals )
            {
                ar = new double[n-1+1, m-1+1];
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=m-1; j++)
                    {
                        v = 0.0;
                        for(i_=0; i_<=n-1;i_++)
                        {
                            v += a2[i,i_]*zref[i_,i1+j];
                        }
                        ar[i,j] = v;
                    }
                }
                for(j=0; j<=m-1; j++)
                {
                    v = 0.0;
                    for(i_=0; i_<=n-1;i_++)
                    {
                        v += a1[i_,j]*ar[i_,j];
                    }
                    if( (double)(v)<(double)(0) )
                    {
                        for(i_=0; i_<=n-1;i_++)
                        {
                            ar[i_,j] = -1*ar[i_,j];
                        }
                    }
                }
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=m-1; j++)
                    {
                        serrors = serrors | (double)(Math.Abs(a1[i,j]-ar[i,j]))>(double)(threshold);
                    }
                }
            }
            
            //
            // Test indexes, transform vectors
            //
            lambda = new double[n-1+1];
            for(i=0; i<=n-1; i++)
            {
                lambda[i] = d[i];
            }
            a1 = new double[n-1+1, n-1+1];
            a2 = new double[n-1+1, n-1+1];
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    a1[i,j] = 2*AP.Math.RandomReal()-1;
                    a2[i,j] = a1[i,j];
                }
            }
            runs = runs+1;
            if( !evd.smatrixtdevdi(ref lambda, ref e, n, 1, i1, i2, ref a1) )
            {
                failc = failc+1;
                return;
            }
            m = i2-i1+1;
            for(k=0; k<=m-1; k++)
            {
                serrors = serrors | (double)(Math.Abs(lambda[k]-lambdaref[i1+k]))>(double)(threshold);
            }
            if( distvals )
            {
                ar = new double[n-1+1, m-1+1];
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=m-1; j++)
                    {
                        v = 0.0;
                        for(i_=0; i_<=n-1;i_++)
                        {
                            v += a2[i,i_]*zref[i_,i1+j];
                        }
                        ar[i,j] = v;
                    }
                }
                for(j=0; j<=m-1; j++)
                {
                    v = 0.0;
                    for(i_=0; i_<=n-1;i_++)
                    {
                        v += a1[i_,j]*ar[i_,j];
                    }
                    if( (double)(v)<(double)(0) )
                    {
                        for(i_=0; i_<=n-1;i_++)
                        {
                            ar[i_,j] = -1*ar[i_,j];
                        }
                    }
                }
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=m-1; j++)
                    {
                        serrors = serrors | (double)(Math.Abs(a1[i,j]-ar[i,j]))>(double)(threshold);
                    }
                }
            }
            
            //
            // Test interval, do not transform vectors
            //
            lambda = new double[n-1+1];
            for(i=0; i<=n-1; i++)
            {
                lambda[i] = d[i];
            }
            z = new double[0+1, 0+1];
            runs = runs+1;
            if( !evd.smatrixtdevdr(ref lambda, ref e, n, 2, a, b, ref m, ref z) )
            {
                failc = failc+1;
                return;
            }
            if( m!=i2-i1+1 )
            {
                failc = failc+1;
                return;
            }
            for(k=0; k<=m-1; k++)
            {
                serrors = serrors | (double)(Math.Abs(lambda[k]-lambdaref[i1+k]))>(double)(threshold);
            }
            if( distvals )
            {
                for(j=0; j<=m-1; j++)
                {
                    v = 0.0;
                    for(i_=0; i_<=n-1;i_++)
                    {
                        v += z[i_,j]*zref[i_,i1+j];
                    }
                    if( (double)(v)<(double)(0) )
                    {
                        for(i_=0; i_<=n-1;i_++)
                        {
                            z[i_,j] = -1*z[i_,j];
                        }
                    }
                }
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=m-1; j++)
                    {
                        serrors = serrors | (double)(Math.Abs(z[i,j]-zref[i,i1+j]))>(double)(threshold);
                    }
                }
            }
            
            //
            // Test indexes, do not transform vectors
            //
            lambda = new double[n-1+1];
            for(i=0; i<=n-1; i++)
            {
                lambda[i] = d[i];
            }
            z = new double[0+1, 0+1];
            runs = runs+1;
            if( !evd.smatrixtdevdi(ref lambda, ref e, n, 2, i1, i2, ref z) )
            {
                failc = failc+1;
                return;
            }
            m = i2-i1+1;
            for(k=0; k<=m-1; k++)
            {
                serrors = serrors | (double)(Math.Abs(lambda[k]-lambdaref[i1+k]))>(double)(threshold);
            }
            if( distvals )
            {
                for(j=0; j<=m-1; j++)
                {
                    v = 0.0;
                    for(i_=0; i_<=n-1;i_++)
                    {
                        v += z[i_,j]*zref[i_,i1+j];
                    }
                    if( (double)(v)<(double)(0) )
                    {
                        for(i_=0; i_<=n-1;i_++)
                        {
                            z[i_,j] = -1*z[i_,j];
                        }
                    }
                }
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=m-1; j++)
                    {
                        serrors = serrors | (double)(Math.Abs(z[i,j]-zref[i,i1+j]))>(double)(threshold);
                    }
                }
            }
        }


        /*************************************************************************
        Non-symmetric problem
        *************************************************************************/
        private static void testnsevdproblem(ref double[,] a,
            int n,
            double threshold,
            ref bool nserrors,
            ref int failc,
            ref int runs)
        {
            double mx = 0;
            int i = 0;
            int j = 0;
            int k = 0;
            int vjob = 0;
            bool needl = new bool();
            bool needr = new bool();
            double[] wr0 = new double[0];
            double[] wi0 = new double[0];
            double[] wr1 = new double[0];
            double[] wi1 = new double[0];
            double[] wr0s = new double[0];
            double[] wi0s = new double[0];
            double[] wr1s = new double[0];
            double[] wi1s = new double[0];
            double[,] vl = new double[0,0];
            double[,] vr = new double[0,0];
            double[] vec1r = new double[0];
            double[] vec1i = new double[0];
            double[] vec2r = new double[0];
            double[] vec2i = new double[0];
            double[] vec3r = new double[0];
            double[] vec3i = new double[0];
            double curwr = 0;
            double curwi = 0;
            double vt = 0;
            double tmp = 0;
            int i_ = 0;

            vec1r = new double[n-1+1];
            vec2r = new double[n-1+1];
            vec3r = new double[n-1+1];
            vec1i = new double[n-1+1];
            vec2i = new double[n-1+1];
            vec3i = new double[n-1+1];
            wr0s = new double[n-1+1];
            wr1s = new double[n-1+1];
            wi0s = new double[n-1+1];
            wi1s = new double[n-1+1];
            mx = 0;
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    if( (double)(Math.Abs(a[i,j]))>(double)(mx) )
                    {
                        mx = Math.Abs(a[i,j]);
                    }
                }
            }
            if( (double)(mx)==(double)(0) )
            {
                mx = 1;
            }
            
            //
            // Load values-only
            //
            runs = runs+1;
            if( !evd.rmatrixevd(a, n, 0, ref wr0, ref wi0, ref vl, ref vr) )
            {
                failc = failc+1;
                return;
            }
            
            //
            // Test different jobs
            //
            for(vjob=1; vjob<=3; vjob++)
            {
                needr = vjob==1 | vjob==3;
                needl = vjob==2 | vjob==3;
                runs = runs+1;
                if( !evd.rmatrixevd(a, n, vjob, ref wr1, ref wi1, ref vl, ref vr) )
                {
                    failc = failc+1;
                    return;
                }
                
                //
                // Test values:
                // 1. sort by real part
                // 2. test
                //
                for(i_=0; i_<=n-1;i_++)
                {
                    wr0s[i_] = wr0[i_];
                }
                for(i_=0; i_<=n-1;i_++)
                {
                    wi0s[i_] = wi0[i_];
                }
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=n-2-i; j++)
                    {
                        if( (double)(wr0s[j])>(double)(wr0s[j+1]) )
                        {
                            tmp = wr0s[j];
                            wr0s[j] = wr0s[j+1];
                            wr0s[j+1] = tmp;
                            tmp = wi0s[j];
                            wi0s[j] = wi0s[j+1];
                            wi0s[j+1] = tmp;
                        }
                    }
                }
                for(i_=0; i_<=n-1;i_++)
                {
                    wr1s[i_] = wr1[i_];
                }
                for(i_=0; i_<=n-1;i_++)
                {
                    wi1s[i_] = wi1[i_];
                }
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=n-2-i; j++)
                    {
                        if( (double)(wr1s[j])>(double)(wr1s[j+1]) )
                        {
                            tmp = wr1s[j];
                            wr1s[j] = wr1s[j+1];
                            wr1s[j+1] = tmp;
                            tmp = wi1s[j];
                            wi1s[j] = wi1s[j+1];
                            wi1s[j+1] = tmp;
                        }
                    }
                }
                for(i=0; i<=n-1; i++)
                {
                    nserrors = nserrors | (double)(Math.Abs(wr0s[i]-wr1s[i]))>(double)(threshold);
                    nserrors = nserrors | (double)(Math.Abs(wi0s[i]-wi1s[i]))>(double)(threshold);
                }
                
                //
                // Test right vectors
                //
                if( needr )
                {
                    k = 0;
                    while( k<=n-1 )
                    {
                        if( (double)(wi1[k])==(double)(0) )
                        {
                            for(i_=0; i_<=n-1;i_++)
                            {
                                vec1r[i_] = vr[i_,k];
                            }
                            for(i=0; i<=n-1; i++)
                            {
                                vec1i[i] = 0;
                            }
                            curwr = wr1[k];
                            curwi = 0;
                        }
                        if( (double)(wi1[k])>(double)(0) )
                        {
                            for(i_=0; i_<=n-1;i_++)
                            {
                                vec1r[i_] = vr[i_,k];
                            }
                            for(i_=0; i_<=n-1;i_++)
                            {
                                vec1i[i_] = vr[i_,k+1];
                            }
                            curwr = wr1[k];
                            curwi = wi1[k];
                        }
                        if( (double)(wi1[k])<(double)(0) )
                        {
                            for(i_=0; i_<=n-1;i_++)
                            {
                                vec1r[i_] = vr[i_,k-1];
                            }
                            for(i_=0; i_<=n-1;i_++)
                            {
                                vec1i[i_] = -vr[i_,k];
                            }
                            curwr = wr1[k];
                            curwi = wi1[k];
                        }
                        for(i=0; i<=n-1; i++)
                        {
                            vt = 0.0;
                            for(i_=0; i_<=n-1;i_++)
                            {
                                vt += a[i,i_]*vec1r[i_];
                            }
                            vec2r[i] = vt;
                            vt = 0.0;
                            for(i_=0; i_<=n-1;i_++)
                            {
                                vt += a[i,i_]*vec1i[i_];
                            }
                            vec2i[i] = vt;
                        }
                        for(i_=0; i_<=n-1;i_++)
                        {
                            vec3r[i_] = curwr*vec1r[i_];
                        }
                        for(i_=0; i_<=n-1;i_++)
                        {
                            vec3r[i_] = vec3r[i_] - curwi*vec1i[i_];
                        }
                        for(i_=0; i_<=n-1;i_++)
                        {
                            vec3i[i_] = curwi*vec1r[i_];
                        }
                        for(i_=0; i_<=n-1;i_++)
                        {
                            vec3i[i_] = vec3i[i_] + curwr*vec1i[i_];
                        }
                        for(i=0; i<=n-1; i++)
                        {
                            nserrors = nserrors | (double)(Math.Abs(vec2r[i]-vec3r[i]))>(double)(threshold);
                            nserrors = nserrors | (double)(Math.Abs(vec2i[i]-vec3i[i]))>(double)(threshold);
                        }
                        k = k+1;
                    }
                }
                
                //
                // Test left vectors
                //
                if( needl )
                {
                    k = 0;
                    while( k<=n-1 )
                    {
                        if( (double)(wi1[k])==(double)(0) )
                        {
                            for(i_=0; i_<=n-1;i_++)
                            {
                                vec1r[i_] = vl[i_,k];
                            }
                            for(i=0; i<=n-1; i++)
                            {
                                vec1i[i] = 0;
                            }
                            curwr = wr1[k];
                            curwi = 0;
                        }
                        if( (double)(wi1[k])>(double)(0) )
                        {
                            for(i_=0; i_<=n-1;i_++)
                            {
                                vec1r[i_] = vl[i_,k];
                            }
                            for(i_=0; i_<=n-1;i_++)
                            {
                                vec1i[i_] = vl[i_,k+1];
                            }
                            curwr = wr1[k];
                            curwi = wi1[k];
                        }
                        if( (double)(wi1[k])<(double)(0) )
                        {
                            for(i_=0; i_<=n-1;i_++)
                            {
                                vec1r[i_] = vl[i_,k-1];
                            }
                            for(i_=0; i_<=n-1;i_++)
                            {
                                vec1i[i_] = -vl[i_,k];
                            }
                            curwr = wr1[k];
                            curwi = wi1[k];
                        }
                        for(j=0; j<=n-1; j++)
                        {
                            vt = 0.0;
                            for(i_=0; i_<=n-1;i_++)
                            {
                                vt += vec1r[i_]*a[i_,j];
                            }
                            vec2r[j] = vt;
                            vt = 0.0;
                            for(i_=0; i_<=n-1;i_++)
                            {
                                vt += vec1i[i_]*a[i_,j];
                            }
                            vec2i[j] = -vt;
                        }
                        for(i_=0; i_<=n-1;i_++)
                        {
                            vec3r[i_] = curwr*vec1r[i_];
                        }
                        for(i_=0; i_<=n-1;i_++)
                        {
                            vec3r[i_] = vec3r[i_] + curwi*vec1i[i_];
                        }
                        for(i_=0; i_<=n-1;i_++)
                        {
                            vec3i[i_] = curwi*vec1r[i_];
                        }
                        for(i_=0; i_<=n-1;i_++)
                        {
                            vec3i[i_] = vec3i[i_] - curwr*vec1i[i_];
                        }
                        for(i=0; i<=n-1; i++)
                        {
                            nserrors = nserrors | (double)(Math.Abs(vec2r[i]-vec3r[i]))>(double)(threshold);
                            nserrors = nserrors | (double)(Math.Abs(vec2i[i]-vec3i[i]))>(double)(threshold);
                        }
                        k = k+1;
                    }
                }
            }
        }


        /*************************************************************************
        Testing EVD subroutines for one N

        NOTES:
        * BIThreshold is a threshold for bisection-and-inverse-iteration subroutines.
          special threshold is needed because these subroutines may have much more
          larger error than QR-based algorithms.
        *************************************************************************/
        private static void testevdset(int n,
            double threshold,
            double bithreshold,
            ref int failc,
            ref int runs,
            ref bool nserrors,
            ref bool serrors,
            ref bool herrors,
            ref bool tderrors,
            ref bool sbierrors,
            ref bool hbierrors,
            ref bool tdbierrors)
        {
            double[,] ra = new double[0,0];
            double[,] ral = new double[0,0];
            double[,] rau = new double[0,0];
            AP.Complex[,] ca = new AP.Complex[0,0];
            AP.Complex[,] cal = new AP.Complex[0,0];
            AP.Complex[,] cau = new AP.Complex[0,0];
            double[] d = new double[0];
            double[] e = new double[0];
            int pass = 0;
            int i = 0;
            int j = 0;
            int mkind = 0;

            
            //
            // Test symmetric problems
            //
            
            //
            // Test symmetric problem: zero, random, sparse matrices.
            //
            ra = new double[n, n];
            ral = new double[n, n];
            rau = new double[n, n];
            ca = new AP.Complex[n, n];
            cal = new AP.Complex[n, n];
            cau = new AP.Complex[n, n];
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    ra[i,j] = 0;
                    ca[i,j] = 0;
                }
            }
            rmatrixsymmetricsplit(ref ra, n, ref ral, ref rau);
            cmatrixhermitiansplit(ref ca, n, ref cal, ref cau);
            testsevdproblem(ref ra, ref ral, ref rau, n, threshold, ref serrors, ref failc, ref runs);
            testhevdproblem(ref ca, ref cal, ref cau, n, threshold, ref herrors, ref failc, ref runs);
            testsevdbiproblem(ref ra, ref ral, ref rau, n, false, bithreshold, ref sbierrors, ref failc, ref runs);
            testhevdbiproblem(ref ca, ref cal, ref cau, n, false, bithreshold, ref hbierrors, ref failc, ref runs);
            for(i=0; i<=n-1; i++)
            {
                for(j=i+1; j<=n-1; j++)
                {
                    ra[i,j] = 2*AP.Math.RandomReal()-1;
                    ca[i,j].x = 2*AP.Math.RandomReal()-1;
                    ca[i,j].y = 2*AP.Math.RandomReal()-1;
                    ra[j,i] = ra[i,j];
                    ca[j,i] = AP.Math.Conj(ca[i,j]);
                }
                ra[i,i] = 2*AP.Math.RandomReal()-1;
                ca[i,i] = 2*AP.Math.RandomReal()-1;
            }
            rmatrixsymmetricsplit(ref ra, n, ref ral, ref rau);
            cmatrixhermitiansplit(ref ca, n, ref cal, ref cau);
            testsevdproblem(ref ra, ref ral, ref rau, n, threshold, ref serrors, ref failc, ref runs);
            testhevdproblem(ref ca, ref cal, ref cau, n, threshold, ref herrors, ref failc, ref runs);
            testsevdbiproblem(ref ra, ref ral, ref rau, n, true, bithreshold, ref sbierrors, ref failc, ref runs);
            testhevdbiproblem(ref ca, ref cal, ref cau, n, true, bithreshold, ref hbierrors, ref failc, ref runs);
            rmatrixfillsparsea(ref ra, n, n, 0.995);
            cmatrixfillsparsea(ref ca, n, n, 0.995);
            for(i=0; i<=n-1; i++)
            {
                for(j=i+1; j<=n-1; j++)
                {
                    ra[j,i] = ra[i,j];
                    ca[j,i] = AP.Math.Conj(ca[i,j]);
                }
                ca[i,i].y = 0;
            }
            rmatrixsymmetricsplit(ref ra, n, ref ral, ref rau);
            cmatrixhermitiansplit(ref ca, n, ref cal, ref cau);
            testsevdproblem(ref ra, ref ral, ref rau, n, threshold, ref serrors, ref failc, ref runs);
            testhevdproblem(ref ca, ref cal, ref cau, n, threshold, ref herrors, ref failc, ref runs);
            testsevdbiproblem(ref ra, ref ral, ref rau, n, false, bithreshold, ref sbierrors, ref failc, ref runs);
            testhevdbiproblem(ref ca, ref cal, ref cau, n, false, bithreshold, ref hbierrors, ref failc, ref runs);
            
            //
            // testing tridiagonal problems
            //
            for(mkind=0; mkind<=4; mkind++)
            {
                d = new double[n];
                if( n>1 )
                {
                    e = new double[n-1];
                }
                if( mkind==0 )
                {
                    
                    //
                    // Zero matrix
                    //
                    for(i=0; i<=n-1; i++)
                    {
                        d[i] = 0;
                    }
                    for(i=0; i<=n-2; i++)
                    {
                        e[i] = 0;
                    }
                }
                if( mkind==1 )
                {
                    
                    //
                    // Diagonal matrix
                    //
                    for(i=0; i<=n-1; i++)
                    {
                        d[i] = 2*AP.Math.RandomReal()-1;
                    }
                    for(i=0; i<=n-2; i++)
                    {
                        e[i] = 0;
                    }
                }
                if( mkind==2 )
                {
                    
                    //
                    // Off-diagonal matrix
                    //
                    for(i=0; i<=n-1; i++)
                    {
                        d[i] = 0;
                    }
                    for(i=0; i<=n-2; i++)
                    {
                        e[i] = 2*AP.Math.RandomReal()-1;
                    }
                }
                if( mkind==3 )
                {
                    
                    //
                    // Dense matrix with blocks
                    //
                    for(i=0; i<=n-1; i++)
                    {
                        d[i] = 2*AP.Math.RandomReal()-1;
                    }
                    for(i=0; i<=n-2; i++)
                    {
                        e[i] = 2*AP.Math.RandomReal()-1;
                    }
                    j = 1;
                    i = 2;
                    while( j<=n-2 )
                    {
                        e[j] = 0;
                        j = j+i;
                        i = i+1;
                    }
                }
                if( mkind==4 )
                {
                    
                    //
                    // dense matrix
                    //
                    for(i=0; i<=n-1; i++)
                    {
                        d[i] = 2*AP.Math.RandomReal()-1;
                    }
                    for(i=0; i<=n-2; i++)
                    {
                        e[i] = 2*AP.Math.RandomReal()-1;
                    }
                }
                testtdevdproblem(ref d, ref e, n, threshold, ref tderrors, ref failc, ref runs);
                testtdevdbiproblem(ref d, ref e, n, mkind==1 | mkind==2 | mkind==4, bithreshold, ref tdbierrors, ref failc, ref runs);
            }
            
            //
            // Test non-symmetric problems
            //
            
            //
            // Test non-symmetric problems: zero, random, sparse matrices.
            //
            ra = new double[n, n];
            ca = new AP.Complex[n, n];
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    ra[i,j] = 0;
                    ca[i,j] = 0;
                }
            }
            testnsevdproblem(ref ra, n, threshold, ref nserrors, ref failc, ref runs);
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    ra[i,j] = 2*AP.Math.RandomReal()-1;
                    ca[i,j].x = 2*AP.Math.RandomReal()-1;
                    ca[i,j].y = 2*AP.Math.RandomReal()-1;
                }
            }
            testnsevdproblem(ref ra, n, threshold, ref nserrors, ref failc, ref runs);
            rmatrixfillsparsea(ref ra, n, n, 0.995);
            cmatrixfillsparsea(ref ca, n, n, 0.995);
            testnsevdproblem(ref ra, n, threshold, ref nserrors, ref failc, ref runs);
        }


        /*************************************************************************
        Silent unit test
        *************************************************************************/
        public static bool testevdunit_test_silent()
        {
            bool result = new bool();

            result = testevd(true);
            return result;
        }


        /*************************************************************************
        Unit test
        *************************************************************************/
        public static bool testevdunit_test()
        {
            bool result = new bool();

            result = testevd(false);
            return result;
        }
    }
}
