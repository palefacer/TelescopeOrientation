
using System;

namespace alglib
{
    public class testortfacunit
    {
        /*************************************************************************
        Main unittest subroutine
        *************************************************************************/
        public static bool testortfac(bool silent)
        {
            bool result = new bool();
            int maxmn = 0;
            double threshold = 0;
            int passcount = 0;
            int mx = 0;
            double[,] ra = new double[0,0];
            AP.Complex[,] ca = new AP.Complex[0,0];
            int m = 0;
            int n = 0;
            int pass = 0;
            int i = 0;
            int j = 0;
            bool rqrerrors = new bool();
            bool rlqerrors = new bool();
            bool cqrerrors = new bool();
            bool clqerrors = new bool();
            bool rbderrors = new bool();
            bool rhesserrors = new bool();
            bool rtderrors = new bool();
            bool ctderrors = new bool();
            bool waserrors = new bool();

            waserrors = false;
            rqrerrors = false;
            rlqerrors = false;
            cqrerrors = false;
            clqerrors = false;
            rbderrors = false;
            rhesserrors = false;
            rtderrors = false;
            ctderrors = false;
            maxmn = 3*ablas.ablasblocksize(ref ra)+1;
            passcount = 1;
            threshold = 5*1000*AP.Math.MachineEpsilon;
            
            //
            // Different problems
            //
            for(mx=1; mx<=maxmn; mx++)
            {
                for(pass=1; pass<=passcount; pass++)
                {
                    
                    //
                    // Rectangular factorizations: QR, LQ, bidiagonal
                    // Matrix types: zero, dense, sparse
                    //
                    n = 1+AP.Math.RandomInteger(mx);
                    m = 1+AP.Math.RandomInteger(mx);
                    if( (double)(AP.Math.RandomReal())>(double)(0.5) )
                    {
                        n = mx;
                    }
                    else
                    {
                        m = mx;
                    }
                    ra = new double[m, n];
                    ca = new AP.Complex[m, n];
                    for(i=0; i<=m-1; i++)
                    {
                        for(j=0; j<=n-1; j++)
                        {
                            ra[i,j] = 0;
                            ca[i,j] = 0;
                        }
                    }
                    testrqrproblem(ref ra, m, n, threshold, ref rqrerrors);
                    testrlqproblem(ref ra, m, n, threshold, ref rlqerrors);
                    testcqrproblem(ref ca, m, n, threshold, ref cqrerrors);
                    testclqproblem(ref ca, m, n, threshold, ref clqerrors);
                    testrbdproblem(ref ra, m, n, threshold, ref rbderrors);
                    for(i=0; i<=m-1; i++)
                    {
                        for(j=0; j<=n-1; j++)
                        {
                            ra[i,j] = 2*AP.Math.RandomReal()-1;
                            ca[i,j].x = 2*AP.Math.RandomReal()-1;
                            ca[i,j].y = 2*AP.Math.RandomReal()-1;
                        }
                    }
                    testrqrproblem(ref ra, m, n, threshold, ref rqrerrors);
                    testrlqproblem(ref ra, m, n, threshold, ref rlqerrors);
                    testcqrproblem(ref ca, m, n, threshold, ref cqrerrors);
                    testclqproblem(ref ca, m, n, threshold, ref clqerrors);
                    testrbdproblem(ref ra, m, n, threshold, ref rbderrors);
                    rmatrixfillsparsea(ref ra, m, n, 0.95);
                    cmatrixfillsparsea(ref ca, m, n, 0.95);
                    testrqrproblem(ref ra, m, n, threshold, ref rqrerrors);
                    testrlqproblem(ref ra, m, n, threshold, ref rlqerrors);
                    testcqrproblem(ref ca, m, n, threshold, ref cqrerrors);
                    testclqproblem(ref ca, m, n, threshold, ref clqerrors);
                    testrbdproblem(ref ra, m, n, threshold, ref rbderrors);
                    
                    //
                    // Square factorizations: Hessenberg, tridiagonal
                    // Matrix types: zero, dense, sparse
                    //
                    ra = new double[mx, mx];
                    ca = new AP.Complex[mx, mx];
                    for(i=0; i<=mx-1; i++)
                    {
                        for(j=0; j<=mx-1; j++)
                        {
                            ra[i,j] = 0;
                            ca[i,j] = 0;
                        }
                    }
                    testrhessproblem(ref ra, mx, threshold, ref rhesserrors);
                    for(i=0; i<=mx-1; i++)
                    {
                        for(j=0; j<=mx-1; j++)
                        {
                            ra[i,j] = 2*AP.Math.RandomReal()-1;
                            ca[i,j].x = 2*AP.Math.RandomReal()-1;
                            ca[i,j].y = 2*AP.Math.RandomReal()-1;
                        }
                    }
                    testrhessproblem(ref ra, mx, threshold, ref rhesserrors);
                    rmatrixfillsparsea(ref ra, mx, mx, 0.95);
                    cmatrixfillsparsea(ref ca, mx, mx, 0.95);
                    testrhessproblem(ref ra, mx, threshold, ref rhesserrors);
                    
                    //
                    // Symetric factorizations: tridiagonal
                    // Matrix types: zero, dense, sparse
                    //
                    ra = new double[mx, mx];
                    ca = new AP.Complex[mx, mx];
                    for(i=0; i<=mx-1; i++)
                    {
                        for(j=0; j<=mx-1; j++)
                        {
                            ra[i,j] = 0;
                            ca[i,j] = 0;
                        }
                    }
                    testrtdproblem(ref ra, mx, threshold, ref rtderrors);
                    testctdproblem(ref ca, mx, threshold, ref ctderrors);
                    for(i=0; i<=mx-1; i++)
                    {
                        for(j=i; j<=mx-1; j++)
                        {
                            ra[i,j] = 2*AP.Math.RandomReal()-1;
                            ca[i,j].x = 2*AP.Math.RandomReal()-1;
                            ca[i,j].y = 2*AP.Math.RandomReal()-1;
                            ra[j,i] = ra[i,j];
                            ca[j,i] = AP.Math.Conj(ca[i,j]);
                        }
                    }
                    for(i=0; i<=mx-1; i++)
                    {
                        ca[i,i] = 2*AP.Math.RandomReal()-1;
                    }
                    testrtdproblem(ref ra, mx, threshold, ref rtderrors);
                    testctdproblem(ref ca, mx, threshold, ref ctderrors);
                    rmatrixfillsparsea(ref ra, mx, mx, 0.95);
                    cmatrixfillsparsea(ref ca, mx, mx, 0.95);
                    for(i=0; i<=mx-1; i++)
                    {
                        for(j=i; j<=mx-1; j++)
                        {
                            ra[j,i] = ra[i,j];
                            ca[j,i] = AP.Math.Conj(ca[i,j]);
                        }
                    }
                    for(i=0; i<=mx-1; i++)
                    {
                        ca[i,i] = 2*AP.Math.RandomReal()-1;
                    }
                    testrtdproblem(ref ra, mx, threshold, ref rtderrors);
                    testctdproblem(ref ca, mx, threshold, ref ctderrors);
                }
            }
            
            //
            // report
            //
            waserrors = rqrerrors | rlqerrors | cqrerrors | clqerrors | rbderrors | rhesserrors | rtderrors | ctderrors;
            if( !silent )
            {
                System.Console.Write("TESTING ORTFAC UNIT");
                System.Console.WriteLine();
                System.Console.Write("RQR ERRORS:                              ");
                if( !rqrerrors )
                {
                    System.Console.Write("OK");
                    System.Console.WriteLine();
                }
                else
                {
                    System.Console.Write("FAILED");
                    System.Console.WriteLine();
                }
                System.Console.Write("RLQ ERRORS:                              ");
                if( !rlqerrors )
                {
                    System.Console.Write("OK");
                    System.Console.WriteLine();
                }
                else
                {
                    System.Console.Write("FAILED");
                    System.Console.WriteLine();
                }
                System.Console.Write("CQR ERRORS:                              ");
                if( !cqrerrors )
                {
                    System.Console.Write("OK");
                    System.Console.WriteLine();
                }
                else
                {
                    System.Console.Write("FAILED");
                    System.Console.WriteLine();
                }
                System.Console.Write("CLQ ERRORS:                              ");
                if( !clqerrors )
                {
                    System.Console.Write("OK");
                    System.Console.WriteLine();
                }
                else
                {
                    System.Console.Write("FAILED");
                    System.Console.WriteLine();
                }
                System.Console.Write("RBD ERRORS:                              ");
                if( !rbderrors )
                {
                    System.Console.Write("OK");
                    System.Console.WriteLine();
                }
                else
                {
                    System.Console.Write("FAILED");
                    System.Console.WriteLine();
                }
                System.Console.Write("RHESS ERRORS:                            ");
                if( !rhesserrors )
                {
                    System.Console.Write("OK");
                    System.Console.WriteLine();
                }
                else
                {
                    System.Console.Write("FAILED");
                    System.Console.WriteLine();
                }
                System.Console.Write("RTD ERRORS:                              ");
                if( !rtderrors )
                {
                    System.Console.Write("OK");
                    System.Console.WriteLine();
                }
                else
                {
                    System.Console.Write("FAILED");
                    System.Console.WriteLine();
                }
                System.Console.Write("CTD ERRORS:                              ");
                if( !ctderrors )
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
        Diff
        *************************************************************************/
        private static double rmatrixdiff(ref double[,] a,
            ref double[,] b,
            int m,
            int n)
        {
            double result = 0;
            int i = 0;
            int j = 0;

            result = 0;
            for(i=0; i<=m-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    result = Math.Max(result, Math.Abs(b[i,j]-a[i,j]));
                }
            }
            return result;
        }


        /*************************************************************************
        Copy
        *************************************************************************/
        private static void rmatrixmakeacopy(ref double[,] a,
            int m,
            int n,
            ref double[,] b)
        {
            int i = 0;
            int j = 0;

            b = new double[m-1+1, n-1+1];
            for(i=0; i<=m-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    b[i,j] = a[i,j];
                }
            }
        }


        /*************************************************************************
        Copy
        *************************************************************************/
        private static void cmatrixmakeacopy(ref AP.Complex[,] a,
            int m,
            int n,
            ref AP.Complex[,] b)
        {
            int i = 0;
            int j = 0;

            b = new AP.Complex[m-1+1, n-1+1];
            for(i=0; i<=m-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    b[i,j] = a[i,j];
                }
            }
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
        Matrix multiplication
        *************************************************************************/
        private static void internalmatrixmatrixmultiply(ref double[,] a,
            int ai1,
            int ai2,
            int aj1,
            int aj2,
            bool transa,
            ref double[,] b,
            int bi1,
            int bi2,
            int bj1,
            int bj2,
            bool transb,
            ref double[,] c,
            int ci1,
            int ci2,
            int cj1,
            int cj2)
        {
            int arows = 0;
            int acols = 0;
            int brows = 0;
            int bcols = 0;
            int crows = 0;
            int ccols = 0;
            int i = 0;
            int j = 0;
            int k = 0;
            int l = 0;
            int r = 0;
            double v = 0;
            double[] work = new double[0];
            double beta = 0;
            double alpha = 0;
            int i_ = 0;
            int i1_ = 0;

            
            //
            // Pre-setup
            //
            k = Math.Max(ai2-ai1+1, aj2-aj1+1);
            k = Math.Max(k, bi2-bi1+1);
            k = Math.Max(k, bj2-bj1+1);
            work = new double[k+1];
            beta = 0;
            alpha = 1;
            
            //
            // Setup
            //
            if( !transa )
            {
                arows = ai2-ai1+1;
                acols = aj2-aj1+1;
            }
            else
            {
                arows = aj2-aj1+1;
                acols = ai2-ai1+1;
            }
            if( !transb )
            {
                brows = bi2-bi1+1;
                bcols = bj2-bj1+1;
            }
            else
            {
                brows = bj2-bj1+1;
                bcols = bi2-bi1+1;
            }
            System.Diagnostics.Debug.Assert(acols==brows, "MatrixMatrixMultiply: incorrect matrix sizes!");
            if( arows<=0 | acols<=0 | brows<=0 | bcols<=0 )
            {
                return;
            }
            crows = arows;
            ccols = bcols;
            
            //
            // Test WORK
            //
            i = Math.Max(arows, acols);
            i = Math.Max(brows, i);
            i = Math.Max(i, bcols);
            work[1] = 0;
            work[i] = 0;
            
            //
            // Prepare C
            //
            if( (double)(beta)==(double)(0) )
            {
                for(i=ci1; i<=ci2; i++)
                {
                    for(j=cj1; j<=cj2; j++)
                    {
                        c[i,j] = 0;
                    }
                }
            }
            else
            {
                for(i=ci1; i<=ci2; i++)
                {
                    for(i_=cj1; i_<=cj2;i_++)
                    {
                        c[i,i_] = beta*c[i,i_];
                    }
                }
            }
            
            //
            // A*B
            //
            if( !transa & !transb )
            {
                for(l=ai1; l<=ai2; l++)
                {
                    for(r=bi1; r<=bi2; r++)
                    {
                        v = alpha*a[l,aj1+r-bi1];
                        k = ci1+l-ai1;
                        i1_ = (bj1) - (cj1);
                        for(i_=cj1; i_<=cj2;i_++)
                        {
                            c[k,i_] = c[k,i_] + v*b[r,i_+i1_];
                        }
                    }
                }
                return;
            }
            
            //
            // A*B'
            //
            if( !transa & transb )
            {
                if( arows*acols<brows*bcols )
                {
                    for(r=bi1; r<=bi2; r++)
                    {
                        for(l=ai1; l<=ai2; l++)
                        {
                            i1_ = (bj1)-(aj1);
                            v = 0.0;
                            for(i_=aj1; i_<=aj2;i_++)
                            {
                                v += a[l,i_]*b[r,i_+i1_];
                            }
                            c[ci1+l-ai1,cj1+r-bi1] = c[ci1+l-ai1,cj1+r-bi1]+alpha*v;
                        }
                    }
                    return;
                }
                else
                {
                    for(l=ai1; l<=ai2; l++)
                    {
                        for(r=bi1; r<=bi2; r++)
                        {
                            i1_ = (bj1)-(aj1);
                            v = 0.0;
                            for(i_=aj1; i_<=aj2;i_++)
                            {
                                v += a[l,i_]*b[r,i_+i1_];
                            }
                            c[ci1+l-ai1,cj1+r-bi1] = c[ci1+l-ai1,cj1+r-bi1]+alpha*v;
                        }
                    }
                    return;
                }
            }
            
            //
            // A'*B
            //
            if( transa & !transb )
            {
                for(l=aj1; l<=aj2; l++)
                {
                    for(r=bi1; r<=bi2; r++)
                    {
                        v = alpha*a[ai1+r-bi1,l];
                        k = ci1+l-aj1;
                        i1_ = (bj1) - (cj1);
                        for(i_=cj1; i_<=cj2;i_++)
                        {
                            c[k,i_] = c[k,i_] + v*b[r,i_+i1_];
                        }
                    }
                }
                return;
            }
            
            //
            // A'*B'
            //
            if( transa & transb )
            {
                if( arows*acols<brows*bcols )
                {
                    for(r=bi1; r<=bi2; r++)
                    {
                        for(i=1; i<=crows; i++)
                        {
                            work[i] = 0.0;
                        }
                        for(l=ai1; l<=ai2; l++)
                        {
                            v = alpha*b[r,bj1+l-ai1];
                            k = cj1+r-bi1;
                            i1_ = (aj1) - (1);
                            for(i_=1; i_<=crows;i_++)
                            {
                                work[i_] = work[i_] + v*a[l,i_+i1_];
                            }
                        }
                        i1_ = (1) - (ci1);
                        for(i_=ci1; i_<=ci2;i_++)
                        {
                            c[i_,k] = c[i_,k] + work[i_+i1_];
                        }
                    }
                    return;
                }
                else
                {
                    for(l=aj1; l<=aj2; l++)
                    {
                        k = ai2-ai1+1;
                        i1_ = (ai1) - (1);
                        for(i_=1; i_<=k;i_++)
                        {
                            work[i_] = a[i_+i1_,l];
                        }
                        for(r=bi1; r<=bi2; r++)
                        {
                            i1_ = (bj1)-(1);
                            v = 0.0;
                            for(i_=1; i_<=k;i_++)
                            {
                                v += work[i_]*b[r,i_+i1_];
                            }
                            c[ci1+l-aj1,cj1+r-bi1] = c[ci1+l-aj1,cj1+r-bi1]+alpha*v;
                        }
                    }
                    return;
                }
            }
        }


        /*************************************************************************
        Problem testing
        *************************************************************************/
        private static void testrqrproblem(ref double[,] a,
            int m,
            int n,
            double threshold,
            ref bool qrerrors)
        {
            int i = 0;
            int j = 0;
            int k = 0;
            double[,] b = new double[0,0];
            double[] taub = new double[0];
            double[,] q = new double[0,0];
            double[,] r = new double[0,0];
            double[,] q2 = new double[0,0];
            double v = 0;
            int i_ = 0;

            
            //
            // Test decompose-and-unpack error
            //
            rmatrixmakeacopy(ref a, m, n, ref b);
            ortfac.rmatrixqr(ref b, m, n, ref taub);
            ortfac.rmatrixqrunpackq(ref b, m, n, ref taub, m, ref q);
            ortfac.rmatrixqrunpackr(ref b, m, n, ref r);
            for(i=0; i<=m-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    v = 0.0;
                    for(i_=0; i_<=m-1;i_++)
                    {
                        v += q[i,i_]*r[i_,j];
                    }
                    qrerrors = qrerrors | (double)(Math.Abs(v-a[i,j]))>(double)(threshold);
                }
            }
            for(i=0; i<=m-1; i++)
            {
                for(j=0; j<=Math.Min(i, n-1)-1; j++)
                {
                    qrerrors = qrerrors | (double)(r[i,j])!=(double)(0);
                }
            }
            for(i=0; i<=m-1; i++)
            {
                for(j=0; j<=m-1; j++)
                {
                    v = 0.0;
                    for(i_=0; i_<=m-1;i_++)
                    {
                        v += q[i,i_]*q[j,i_];
                    }
                    if( i==j )
                    {
                        v = v-1;
                    }
                    qrerrors = qrerrors | (double)(Math.Abs(v))>=(double)(threshold);
                }
            }
            
            //
            // Test for other errors
            //
            k = 1+AP.Math.RandomInteger(m);
            ortfac.rmatrixqrunpackq(ref b, m, n, ref taub, k, ref q2);
            for(i=0; i<=m-1; i++)
            {
                for(j=0; j<=k-1; j++)
                {
                    qrerrors = qrerrors | (double)(Math.Abs(q2[i,j]-q[i,j]))>(double)(10*AP.Math.MachineEpsilon);
                }
            }
        }


        /*************************************************************************
        Problem testing
        *************************************************************************/
        private static void testcqrproblem(ref AP.Complex[,] a,
            int m,
            int n,
            double threshold,
            ref bool qrerrors)
        {
            int i = 0;
            int j = 0;
            int k = 0;
            AP.Complex[,] b = new AP.Complex[0,0];
            AP.Complex[] taub = new AP.Complex[0];
            AP.Complex[,] q = new AP.Complex[0,0];
            AP.Complex[,] r = new AP.Complex[0,0];
            AP.Complex[,] q2 = new AP.Complex[0,0];
            AP.Complex v = 0;
            int i_ = 0;

            
            //
            // Test decompose-and-unpack error
            //
            cmatrixmakeacopy(ref a, m, n, ref b);
            ortfac.cmatrixqr(ref b, m, n, ref taub);
            ortfac.cmatrixqrunpackq(ref b, m, n, ref taub, m, ref q);
            ortfac.cmatrixqrunpackr(ref b, m, n, ref r);
            for(i=0; i<=m-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    v = 0.0;
                    for(i_=0; i_<=m-1;i_++)
                    {
                        v += q[i,i_]*r[i_,j];
                    }
                    qrerrors = qrerrors | (double)(AP.Math.AbsComplex(v-a[i,j]))>(double)(threshold);
                }
            }
            for(i=0; i<=m-1; i++)
            {
                for(j=0; j<=Math.Min(i, n-1)-1; j++)
                {
                    qrerrors = qrerrors | r[i,j]!=0;
                }
            }
            for(i=0; i<=m-1; i++)
            {
                for(j=0; j<=m-1; j++)
                {
                    v = 0.0;
                    for(i_=0; i_<=m-1;i_++)
                    {
                        v += q[i,i_]*AP.Math.Conj(q[j,i_]);
                    }
                    if( i==j )
                    {
                        v = v-1;
                    }
                    qrerrors = qrerrors | (double)(AP.Math.AbsComplex(v))>=(double)(threshold);
                }
            }
            
            //
            // Test for other errors
            //
            k = 1+AP.Math.RandomInteger(m);
            ortfac.cmatrixqrunpackq(ref b, m, n, ref taub, k, ref q2);
            for(i=0; i<=m-1; i++)
            {
                for(j=0; j<=k-1; j++)
                {
                    qrerrors = qrerrors | (double)(AP.Math.AbsComplex(q2[i,j]-q[i,j]))>(double)(10*AP.Math.MachineEpsilon);
                }
            }
        }


        /*************************************************************************
        Problem testing
        *************************************************************************/
        private static void testrlqproblem(ref double[,] a,
            int m,
            int n,
            double threshold,
            ref bool lqerrors)
        {
            int i = 0;
            int j = 0;
            int k = 0;
            double[,] b = new double[0,0];
            double[] taub = new double[0];
            double[,] q = new double[0,0];
            double[,] l = new double[0,0];
            double[,] q2 = new double[0,0];
            double v = 0;
            int i_ = 0;

            
            //
            // Test decompose-and-unpack error
            //
            rmatrixmakeacopy(ref a, m, n, ref b);
            ortfac.rmatrixlq(ref b, m, n, ref taub);
            ortfac.rmatrixlqunpackq(ref b, m, n, ref taub, n, ref q);
            ortfac.rmatrixlqunpackl(ref b, m, n, ref l);
            for(i=0; i<=m-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    v = 0.0;
                    for(i_=0; i_<=n-1;i_++)
                    {
                        v += l[i,i_]*q[i_,j];
                    }
                    lqerrors = lqerrors | (double)(Math.Abs(v-a[i,j]))>=(double)(threshold);
                }
            }
            for(i=0; i<=m-1; i++)
            {
                for(j=Math.Min(i, n-1)+1; j<=n-1; j++)
                {
                    lqerrors = lqerrors | (double)(l[i,j])!=(double)(0);
                }
            }
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    v = 0.0;
                    for(i_=0; i_<=n-1;i_++)
                    {
                        v += q[i,i_]*q[j,i_];
                    }
                    if( i==j )
                    {
                        v = v-1;
                    }
                    lqerrors = lqerrors | (double)(Math.Abs(v))>=(double)(threshold);
                }
            }
            
            //
            // Test for other errors
            //
            k = 1+AP.Math.RandomInteger(n);
            ortfac.rmatrixlqunpackq(ref b, m, n, ref taub, k, ref q2);
            for(i=0; i<=k-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    lqerrors = lqerrors | (double)(Math.Abs(q2[i,j]-q[i,j]))>(double)(10*AP.Math.MachineEpsilon);
                }
            }
        }


        /*************************************************************************
        Problem testing
        *************************************************************************/
        private static void testclqproblem(ref AP.Complex[,] a,
            int m,
            int n,
            double threshold,
            ref bool lqerrors)
        {
            int i = 0;
            int j = 0;
            int k = 0;
            AP.Complex[,] b = new AP.Complex[0,0];
            AP.Complex[] taub = new AP.Complex[0];
            AP.Complex[,] q = new AP.Complex[0,0];
            AP.Complex[,] l = new AP.Complex[0,0];
            AP.Complex[,] q2 = new AP.Complex[0,0];
            AP.Complex v = 0;
            int i_ = 0;

            
            //
            // Test decompose-and-unpack error
            //
            cmatrixmakeacopy(ref a, m, n, ref b);
            ortfac.cmatrixlq(ref b, m, n, ref taub);
            ortfac.cmatrixlqunpackq(ref b, m, n, ref taub, n, ref q);
            ortfac.cmatrixlqunpackl(ref b, m, n, ref l);
            for(i=0; i<=m-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    v = 0.0;
                    for(i_=0; i_<=n-1;i_++)
                    {
                        v += l[i,i_]*q[i_,j];
                    }
                    lqerrors = lqerrors | (double)(AP.Math.AbsComplex(v-a[i,j]))>=(double)(threshold);
                }
            }
            for(i=0; i<=m-1; i++)
            {
                for(j=Math.Min(i, n-1)+1; j<=n-1; j++)
                {
                    lqerrors = lqerrors | l[i,j]!=0;
                }
            }
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    v = 0.0;
                    for(i_=0; i_<=n-1;i_++)
                    {
                        v += q[i,i_]*AP.Math.Conj(q[j,i_]);
                    }
                    if( i==j )
                    {
                        v = v-1;
                    }
                    lqerrors = lqerrors | (double)(AP.Math.AbsComplex(v))>=(double)(threshold);
                }
            }
            
            //
            // Test for other errors
            //
            k = 1+AP.Math.RandomInteger(n);
            ortfac.cmatrixlqunpackq(ref b, m, n, ref taub, k, ref q2);
            for(i=0; i<=k-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    lqerrors = lqerrors | (double)(AP.Math.AbsComplex(q2[i,j]-q[i,j]))>(double)(10*AP.Math.MachineEpsilon);
                }
            }
        }


        /*************************************************************************
        Problem testing
        *************************************************************************/
        private static void testrbdproblem(ref double[,] a,
            int m,
            int n,
            double threshold,
            ref bool bderrors)
        {
            int i = 0;
            int j = 0;
            int k = 0;
            double[,] t = new double[0,0];
            double[,] pt = new double[0,0];
            double[,] q = new double[0,0];
            double[,] r = new double[0,0];
            double[,] bd = new double[0,0];
            double[,] x = new double[0,0];
            double[,] r1 = new double[0,0];
            double[,] r2 = new double[0,0];
            double[] taup = new double[0];
            double[] tauq = new double[0];
            double[] d = new double[0];
            double[] e = new double[0];
            bool up = new bool();
            double v = 0;
            int mtsize = 0;
            int i_ = 0;

            
            //
            // Bidiagonal decomposition error
            //
            rmatrixmakeacopy(ref a, m, n, ref t);
            ortfac.rmatrixbd(ref t, m, n, ref tauq, ref taup);
            ortfac.rmatrixbdunpackq(ref t, m, n, ref tauq, m, ref q);
            ortfac.rmatrixbdunpackpt(ref t, m, n, ref taup, n, ref pt);
            ortfac.rmatrixbdunpackdiagonals(ref t, m, n, ref up, ref d, ref e);
            bd = new double[m, n];
            for(i=0; i<=m-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    bd[i,j] = 0;
                }
            }
            for(i=0; i<=Math.Min(m, n)-1; i++)
            {
                bd[i,i] = d[i];
            }
            if( up )
            {
                for(i=0; i<=Math.Min(m, n)-2; i++)
                {
                    bd[i,i+1] = e[i];
                }
            }
            else
            {
                for(i=0; i<=Math.Min(m, n)-2; i++)
                {
                    bd[i+1,i] = e[i];
                }
            }
            r = new double[m, n];
            for(i=0; i<=m-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    v = 0.0;
                    for(i_=0; i_<=m-1;i_++)
                    {
                        v += q[i,i_]*bd[i_,j];
                    }
                    r[i,j] = v;
                }
            }
            for(i=0; i<=m-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    v = 0.0;
                    for(i_=0; i_<=n-1;i_++)
                    {
                        v += r[i,i_]*pt[i_,j];
                    }
                    bderrors = bderrors | (double)(Math.Abs(v-a[i,j]))>(double)(threshold);
                }
            }
            
            //
            // Orthogonality test for Q/PT
            //
            for(i=0; i<=m-1; i++)
            {
                for(j=0; j<=m-1; j++)
                {
                    v = 0.0;
                    for(i_=0; i_<=m-1;i_++)
                    {
                        v += q[i_,i]*q[i_,j];
                    }
                    if( i==j )
                    {
                        bderrors = bderrors | (double)(Math.Abs(v-1))>(double)(threshold);
                    }
                    else
                    {
                        bderrors = bderrors | (double)(Math.Abs(v))>(double)(threshold);
                    }
                }
            }
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    v = 0.0;
                    for(i_=0; i_<=n-1;i_++)
                    {
                        v += pt[i,i_]*pt[j,i_];
                    }
                    if( i==j )
                    {
                        bderrors = bderrors | (double)(Math.Abs(v-1))>(double)(threshold);
                    }
                    else
                    {
                        bderrors = bderrors | (double)(Math.Abs(v))>(double)(threshold);
                    }
                }
            }
            
            //
            // Partial unpacking test
            //
            k = 1+AP.Math.RandomInteger(m);
            ortfac.rmatrixbdunpackq(ref t, m, n, ref tauq, k, ref r);
            for(i=0; i<=m-1; i++)
            {
                for(j=0; j<=k-1; j++)
                {
                    bderrors = bderrors | (double)(Math.Abs(r[i,j]-q[i,j]))>(double)(10*AP.Math.MachineEpsilon);
                }
            }
            k = 1+AP.Math.RandomInteger(n);
            ortfac.rmatrixbdunpackpt(ref t, m, n, ref taup, k, ref r);
            for(i=0; i<=k-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    bderrors = bderrors | (double)(r[i,j]-pt[i,j])!=(double)(0);
                }
            }
            
            //
            // Multiplication test
            //
            x = new double[Math.Max(m, n)-1+1, Math.Max(m, n)-1+1];
            r = new double[Math.Max(m, n)-1+1, Math.Max(m, n)-1+1];
            r1 = new double[Math.Max(m, n)-1+1, Math.Max(m, n)-1+1];
            r2 = new double[Math.Max(m, n)-1+1, Math.Max(m, n)-1+1];
            for(i=0; i<=Math.Max(m, n)-1; i++)
            {
                for(j=0; j<=Math.Max(m, n)-1; j++)
                {
                    x[i,j] = 2*AP.Math.RandomReal()-1;
                }
            }
            mtsize = 1+AP.Math.RandomInteger(Math.Max(m, n));
            rmatrixmakeacopy(ref x, mtsize, m, ref r);
            internalmatrixmatrixmultiply(ref r, 0, mtsize-1, 0, m-1, false, ref q, 0, m-1, 0, m-1, false, ref r1, 0, mtsize-1, 0, m-1);
            rmatrixmakeacopy(ref x, mtsize, m, ref r2);
            ortfac.rmatrixbdmultiplybyq(ref t, m, n, ref tauq, ref r2, mtsize, m, true, false);
            bderrors = bderrors | (double)(rmatrixdiff(ref r1, ref r2, mtsize, m))>(double)(threshold);
            rmatrixmakeacopy(ref x, mtsize, m, ref r);
            internalmatrixmatrixmultiply(ref r, 0, mtsize-1, 0, m-1, false, ref q, 0, m-1, 0, m-1, true, ref r1, 0, mtsize-1, 0, m-1);
            rmatrixmakeacopy(ref x, mtsize, m, ref r2);
            ortfac.rmatrixbdmultiplybyq(ref t, m, n, ref tauq, ref r2, mtsize, m, true, true);
            bderrors = bderrors | (double)(rmatrixdiff(ref r1, ref r2, mtsize, m))>(double)(threshold);
            rmatrixmakeacopy(ref x, m, mtsize, ref r);
            internalmatrixmatrixmultiply(ref q, 0, m-1, 0, m-1, false, ref r, 0, m-1, 0, mtsize-1, false, ref r1, 0, m-1, 0, mtsize-1);
            rmatrixmakeacopy(ref x, m, mtsize, ref r2);
            ortfac.rmatrixbdmultiplybyq(ref t, m, n, ref tauq, ref r2, m, mtsize, false, false);
            bderrors = bderrors | (double)(rmatrixdiff(ref r1, ref r2, m, mtsize))>(double)(threshold);
            rmatrixmakeacopy(ref x, m, mtsize, ref r);
            internalmatrixmatrixmultiply(ref q, 0, m-1, 0, m-1, true, ref r, 0, m-1, 0, mtsize-1, false, ref r1, 0, m-1, 0, mtsize-1);
            rmatrixmakeacopy(ref x, m, mtsize, ref r2);
            ortfac.rmatrixbdmultiplybyq(ref t, m, n, ref tauq, ref r2, m, mtsize, false, true);
            bderrors = bderrors | (double)(rmatrixdiff(ref r1, ref r2, m, mtsize))>(double)(threshold);
            rmatrixmakeacopy(ref x, mtsize, n, ref r);
            internalmatrixmatrixmultiply(ref r, 0, mtsize-1, 0, n-1, false, ref pt, 0, n-1, 0, n-1, true, ref r1, 0, mtsize-1, 0, n-1);
            rmatrixmakeacopy(ref x, mtsize, n, ref r2);
            ortfac.rmatrixbdmultiplybyp(ref t, m, n, ref taup, ref r2, mtsize, n, true, false);
            bderrors = bderrors | (double)(rmatrixdiff(ref r1, ref r2, mtsize, n))>(double)(threshold);
            rmatrixmakeacopy(ref x, mtsize, n, ref r);
            internalmatrixmatrixmultiply(ref r, 0, mtsize-1, 0, n-1, false, ref pt, 0, n-1, 0, n-1, false, ref r1, 0, mtsize-1, 0, n-1);
            rmatrixmakeacopy(ref x, mtsize, n, ref r2);
            ortfac.rmatrixbdmultiplybyp(ref t, m, n, ref taup, ref r2, mtsize, n, true, true);
            bderrors = bderrors | (double)(rmatrixdiff(ref r1, ref r2, mtsize, n))>(double)(threshold);
            rmatrixmakeacopy(ref x, n, mtsize, ref r);
            internalmatrixmatrixmultiply(ref pt, 0, n-1, 0, n-1, true, ref r, 0, n-1, 0, mtsize-1, false, ref r1, 0, n-1, 0, mtsize-1);
            rmatrixmakeacopy(ref x, n, mtsize, ref r2);
            ortfac.rmatrixbdmultiplybyp(ref t, m, n, ref taup, ref r2, n, mtsize, false, false);
            bderrors = bderrors | (double)(rmatrixdiff(ref r1, ref r2, n, mtsize))>(double)(threshold);
            rmatrixmakeacopy(ref x, n, mtsize, ref r);
            internalmatrixmatrixmultiply(ref pt, 0, n-1, 0, n-1, false, ref r, 0, n-1, 0, mtsize-1, false, ref r1, 0, n-1, 0, mtsize-1);
            rmatrixmakeacopy(ref x, n, mtsize, ref r2);
            ortfac.rmatrixbdmultiplybyp(ref t, m, n, ref taup, ref r2, n, mtsize, false, true);
            bderrors = bderrors | (double)(rmatrixdiff(ref r1, ref r2, n, mtsize))>(double)(threshold);
        }


        /*************************************************************************
        Problem testing
        *************************************************************************/
        private static void testrhessproblem(ref double[,] a,
            int n,
            double threshold,
            ref bool hesserrors)
        {
            double[,] b = new double[0,0];
            double[,] h = new double[0,0];
            double[,] q = new double[0,0];
            double[,] t1 = new double[0,0];
            double[,] t2 = new double[0,0];
            double[] tau = new double[0];
            int i = 0;
            int j = 0;
            double v = 0;
            int i_ = 0;

            rmatrixmakeacopy(ref a, n, n, ref b);
            
            //
            // Decomposition
            //
            ortfac.rmatrixhessenberg(ref b, n, ref tau);
            ortfac.rmatrixhessenbergunpackq(ref b, n, ref tau, ref q);
            ortfac.rmatrixhessenbergunpackh(ref b, n, ref h);
            
            //
            // Matrix properties
            //
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    v = 0.0;
                    for(i_=0; i_<=n-1;i_++)
                    {
                        v += q[i_,i]*q[i_,j];
                    }
                    if( i==j )
                    {
                        v = v-1;
                    }
                    hesserrors = hesserrors | (double)(Math.Abs(v))>(double)(threshold);
                }
            }
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=i-2; j++)
                {
                    hesserrors = hesserrors | (double)(h[i,j])!=(double)(0);
                }
            }
            
            //
            // Decomposition error
            //
            t1 = new double[n, n];
            t2 = new double[n, n];
            internalmatrixmatrixmultiply(ref q, 0, n-1, 0, n-1, false, ref h, 0, n-1, 0, n-1, false, ref t1, 0, n-1, 0, n-1);
            internalmatrixmatrixmultiply(ref t1, 0, n-1, 0, n-1, false, ref q, 0, n-1, 0, n-1, true, ref t2, 0, n-1, 0, n-1);
            hesserrors = hesserrors | (double)(rmatrixdiff(ref t2, ref a, n, n))>(double)(threshold);
        }


        /*************************************************************************
        Tridiagonal tester
        *************************************************************************/
        private static void testrtdproblem(ref double[,] a,
            int n,
            double threshold,
            ref bool tderrors)
        {
            int i = 0;
            int j = 0;
            double[,] ua = new double[0,0];
            double[,] la = new double[0,0];
            double[,] t = new double[0,0];
            double[,] q = new double[0,0];
            double[,] t2 = new double[0,0];
            double[,] t3 = new double[0,0];
            double[] tau = new double[0];
            double[] d = new double[0];
            double[] e = new double[0];
            double v = 0;
            int i_ = 0;

            ua = new double[n-1+1, n-1+1];
            la = new double[n-1+1, n-1+1];
            t = new double[n-1+1, n-1+1];
            q = new double[n-1+1, n-1+1];
            t2 = new double[n-1+1, n-1+1];
            t3 = new double[n-1+1, n-1+1];
            
            //
            // fill
            //
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    ua[i,j] = 0;
                }
            }
            for(i=0; i<=n-1; i++)
            {
                for(j=i; j<=n-1; j++)
                {
                    ua[i,j] = a[i,j];
                }
            }
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    la[i,j] = 0;
                }
            }
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=i; j++)
                {
                    la[i,j] = a[i,j];
                }
            }
            
            //
            // Test 2tridiagonal: upper
            //
            ortfac.smatrixtd(ref ua, n, true, ref tau, ref d, ref e);
            ortfac.smatrixtdunpackq(ref ua, n, true, ref tau, ref q);
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    t[i,j] = 0;
                }
            }
            for(i=0; i<=n-1; i++)
            {
                t[i,i] = d[i];
            }
            for(i=0; i<=n-2; i++)
            {
                t[i,i+1] = e[i];
                t[i+1,i] = e[i];
            }
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    v = 0.0;
                    for(i_=0; i_<=n-1;i_++)
                    {
                        v += q[i_,i]*a[i_,j];
                    }
                    t2[i,j] = v;
                }
            }
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    v = 0.0;
                    for(i_=0; i_<=n-1;i_++)
                    {
                        v += t2[i,i_]*q[i_,j];
                    }
                    t3[i,j] = v;
                }
            }
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    tderrors = tderrors | (double)(Math.Abs(t3[i,j]-t[i,j]))>(double)(threshold);
                }
            }
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    v = 0.0;
                    for(i_=0; i_<=n-1;i_++)
                    {
                        v += q[i,i_]*q[j,i_];
                    }
                    if( i==j )
                    {
                        v = v-1;
                    }
                    tderrors = tderrors | (double)(Math.Abs(v))>(double)(threshold);
                }
            }
            
            //
            // Test 2tridiagonal: lower
            //
            ortfac.smatrixtd(ref la, n, false, ref tau, ref d, ref e);
            ortfac.smatrixtdunpackq(ref la, n, false, ref tau, ref q);
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    t[i,j] = 0;
                }
            }
            for(i=0; i<=n-1; i++)
            {
                t[i,i] = d[i];
            }
            for(i=0; i<=n-2; i++)
            {
                t[i,i+1] = e[i];
                t[i+1,i] = e[i];
            }
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    v = 0.0;
                    for(i_=0; i_<=n-1;i_++)
                    {
                        v += q[i_,i]*a[i_,j];
                    }
                    t2[i,j] = v;
                }
            }
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    v = 0.0;
                    for(i_=0; i_<=n-1;i_++)
                    {
                        v += t2[i,i_]*q[i_,j];
                    }
                    t3[i,j] = v;
                }
            }
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    tderrors = tderrors | (double)(Math.Abs(t3[i,j]-t[i,j]))>(double)(threshold);
                }
            }
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    v = 0.0;
                    for(i_=0; i_<=n-1;i_++)
                    {
                        v += q[i,i_]*q[j,i_];
                    }
                    if( i==j )
                    {
                        v = v-1;
                    }
                    tderrors = tderrors | (double)(Math.Abs(v))>(double)(threshold);
                }
            }
        }


        /*************************************************************************
        Hermitian problem tester
        *************************************************************************/
        private static void testctdproblem(ref AP.Complex[,] a,
            int n,
            double threshold,
            ref bool tderrors)
        {
            int i = 0;
            int j = 0;
            AP.Complex[,] ua = new AP.Complex[0,0];
            AP.Complex[,] la = new AP.Complex[0,0];
            AP.Complex[,] t = new AP.Complex[0,0];
            AP.Complex[,] q = new AP.Complex[0,0];
            AP.Complex[,] t2 = new AP.Complex[0,0];
            AP.Complex[,] t3 = new AP.Complex[0,0];
            AP.Complex[] tau = new AP.Complex[0];
            double[] d = new double[0];
            double[] e = new double[0];
            AP.Complex v = 0;
            int i_ = 0;

            ua = new AP.Complex[n-1+1, n-1+1];
            la = new AP.Complex[n-1+1, n-1+1];
            t = new AP.Complex[n-1+1, n-1+1];
            q = new AP.Complex[n-1+1, n-1+1];
            t2 = new AP.Complex[n-1+1, n-1+1];
            t3 = new AP.Complex[n-1+1, n-1+1];
            
            //
            // fill
            //
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    ua[i,j] = 0;
                }
            }
            for(i=0; i<=n-1; i++)
            {
                for(j=i; j<=n-1; j++)
                {
                    ua[i,j] = a[i,j];
                }
            }
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    la[i,j] = 0;
                }
            }
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=i; j++)
                {
                    la[i,j] = a[i,j];
                }
            }
            
            //
            // Test 2tridiagonal: upper
            //
            ortfac.hmatrixtd(ref ua, n, true, ref tau, ref d, ref e);
            ortfac.hmatrixtdunpackq(ref ua, n, true, ref tau, ref q);
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    t[i,j] = 0;
                }
            }
            for(i=0; i<=n-1; i++)
            {
                t[i,i] = d[i];
            }
            for(i=0; i<=n-2; i++)
            {
                t[i,i+1] = e[i];
                t[i+1,i] = e[i];
            }
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    v = 0.0;
                    for(i_=0; i_<=n-1;i_++)
                    {
                        v += AP.Math.Conj(q[i_,i])*a[i_,j];
                    }
                    t2[i,j] = v;
                }
            }
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    v = 0.0;
                    for(i_=0; i_<=n-1;i_++)
                    {
                        v += t2[i,i_]*q[i_,j];
                    }
                    t3[i,j] = v;
                }
            }
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    tderrors = tderrors | (double)(AP.Math.AbsComplex(t3[i,j]-t[i,j]))>(double)(threshold);
                }
            }
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    v = 0.0;
                    for(i_=0; i_<=n-1;i_++)
                    {
                        v += q[i,i_]*AP.Math.Conj(q[j,i_]);
                    }
                    if( i==j )
                    {
                        v = v-1;
                    }
                    tderrors = tderrors | (double)(AP.Math.AbsComplex(v))>(double)(threshold);
                }
            }
            
            //
            // Test 2tridiagonal: lower
            //
            ortfac.hmatrixtd(ref la, n, false, ref tau, ref d, ref e);
            ortfac.hmatrixtdunpackq(ref la, n, false, ref tau, ref q);
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    t[i,j] = 0;
                }
            }
            for(i=0; i<=n-1; i++)
            {
                t[i,i] = d[i];
            }
            for(i=0; i<=n-2; i++)
            {
                t[i,i+1] = e[i];
                t[i+1,i] = e[i];
            }
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    v = 0.0;
                    for(i_=0; i_<=n-1;i_++)
                    {
                        v += AP.Math.Conj(q[i_,i])*a[i_,j];
                    }
                    t2[i,j] = v;
                }
            }
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    v = 0.0;
                    for(i_=0; i_<=n-1;i_++)
                    {
                        v += t2[i,i_]*q[i_,j];
                    }
                    t3[i,j] = v;
                }
            }
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    tderrors = tderrors | (double)(AP.Math.AbsComplex(t3[i,j]-t[i,j]))>(double)(threshold);
                }
            }
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    v = 0.0;
                    for(i_=0; i_<=n-1;i_++)
                    {
                        v += q[i,i_]*AP.Math.Conj(q[j,i_]);
                    }
                    if( i==j )
                    {
                        v = v-1;
                    }
                    tderrors = tderrors | (double)(AP.Math.AbsComplex(v))>(double)(threshold);
                }
            }
        }


        /*************************************************************************
        Silent unit test
        *************************************************************************/
        public static bool testortfacunit_test_silent()
        {
            bool result = new bool();

            result = testortfac(true);
            return result;
        }


        /*************************************************************************
        Unit test
        *************************************************************************/
        public static bool testortfacunit_test()
        {
            bool result = new bool();

            result = testortfac(false);
            return result;
        }
    }
}
