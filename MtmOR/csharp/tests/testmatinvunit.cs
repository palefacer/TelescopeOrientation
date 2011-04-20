
using System;

namespace alglib
{
    public class testmatinvunit
    {
        /*************************************************************************
        Test
        *************************************************************************/
        public static bool testmatinv(bool silent)
        {
            bool result = new bool();
            int maxrn = 0;
            int maxcn = 0;
            int passcount = 0;
            double threshold = 0;
            double rcondtol = 0;
            bool rtrerrors = new bool();
            bool ctrerrors = new bool();
            bool rerrors = new bool();
            bool cerrors = new bool();
            bool spderrors = new bool();
            bool hpderrors = new bool();
            bool waserrors = new bool();
            double[,] emptyra = new double[0,0];
            double[,] emptyca = new double[0,0];

            maxrn = 3*ablas.ablasblocksize(ref emptyra)+1;
            maxcn = 3*ablas.ablasblocksize(ref emptyca)+1;
            passcount = 1;
            threshold = 10000*AP.Math.MachineEpsilon;
            rcondtol = 0.01;
            rtrerrors = false;
            ctrerrors = false;
            rerrors = false;
            cerrors = false;
            spderrors = false;
            hpderrors = false;
            testrtrinv(maxrn, passcount, threshold, ref rtrerrors);
            testctrinv(maxcn, passcount, threshold, ref ctrerrors);
            testrinv(maxrn, passcount, threshold, ref rerrors);
            testspdinv(maxrn, passcount, threshold, ref spderrors);
            testcinv(maxcn, passcount, threshold, ref cerrors);
            testhpdinv(maxcn, passcount, threshold, ref hpderrors);
            waserrors = rtrerrors | ctrerrors | rerrors | cerrors | spderrors | hpderrors;
            if( !silent )
            {
                System.Console.Write("TESTING MATINV");
                System.Console.WriteLine();
                System.Console.Write("* REAL TRIANGULAR:                        ");
                if( rtrerrors )
                {
                    System.Console.Write("FAILED");
                    System.Console.WriteLine();
                }
                else
                {
                    System.Console.Write("OK");
                    System.Console.WriteLine();
                }
                System.Console.Write("* COMPLEX TRIANGULAR:                     ");
                if( ctrerrors )
                {
                    System.Console.Write("FAILED");
                    System.Console.WriteLine();
                }
                else
                {
                    System.Console.Write("OK");
                    System.Console.WriteLine();
                }
                System.Console.Write("* REAL:                                   ");
                if( rerrors )
                {
                    System.Console.Write("FAILED");
                    System.Console.WriteLine();
                }
                else
                {
                    System.Console.Write("OK");
                    System.Console.WriteLine();
                }
                System.Console.Write("* COMPLEX:                                ");
                if( cerrors )
                {
                    System.Console.Write("FAILED");
                    System.Console.WriteLine();
                }
                else
                {
                    System.Console.Write("OK");
                    System.Console.WriteLine();
                }
                System.Console.Write("* SPD:                                    ");
                if( spderrors )
                {
                    System.Console.Write("FAILED");
                    System.Console.WriteLine();
                }
                else
                {
                    System.Console.Write("OK");
                    System.Console.WriteLine();
                }
                System.Console.Write("* HPD:                                    ");
                if( hpderrors )
                {
                    System.Console.Write("FAILED");
                    System.Console.WriteLine();
                }
                else
                {
                    System.Console.Write("OK");
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
            }
            result = !waserrors;
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
        Checks whether inverse is correct
        Returns True on success.
        *************************************************************************/
        private static bool rmatrixcheckinverse(ref double[,] a,
            ref double[,] inva,
            int n,
            double threshold,
            int info,
            ref matinv.matinvreport rep)
        {
            bool result = new bool();
            int i = 0;
            int j = 0;
            double v = 0;
            int i_ = 0;

            result = true;
            if( info<=0 )
            {
                result = false;
            }
            else
            {
                result = result & !((double)(rep.r1)<(double)(100*AP.Math.MachineEpsilon) | (double)(rep.r1)>(double)(1+1000*AP.Math.MachineEpsilon));
                result = result & !((double)(rep.rinf)<(double)(100*AP.Math.MachineEpsilon) | (double)(rep.rinf)>(double)(1+1000*AP.Math.MachineEpsilon));
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=n-1; j++)
                    {
                        v = 0.0;
                        for(i_=0; i_<=n-1;i_++)
                        {
                            v += a[i,i_]*inva[i_,j];
                        }
                        if( i==j )
                        {
                            v = v-1;
                        }
                        result = result & (double)(Math.Abs(v))<=(double)(threshold);
                    }
                }
            }
            return result;
        }


        /*************************************************************************
        Checks whether inverse is correct
        Returns True on success.
        *************************************************************************/
        private static bool spdmatrixcheckinverse(double[,] a,
            double[,] inva,
            bool isupper,
            int n,
            double threshold,
            int info,
            ref matinv.matinvreport rep)
        {
            bool result = new bool();
            int i = 0;
            int j = 0;
            double v = 0;
            int i_ = 0;

            a = (double[,])a.Clone();
            inva = (double[,])inva.Clone();

            for(i=0; i<=n-2; i++)
            {
                if( isupper )
                {
                    for(i_=i+1; i_<=n-1;i_++)
                    {
                        a[i_,i] = a[i,i_];
                    }
                    for(i_=i+1; i_<=n-1;i_++)
                    {
                        inva[i_,i] = inva[i,i_];
                    }
                }
                else
                {
                    for(i_=i+1; i_<=n-1;i_++)
                    {
                        a[i,i_] = a[i_,i];
                    }
                    for(i_=i+1; i_<=n-1;i_++)
                    {
                        inva[i,i_] = inva[i_,i];
                    }
                }
            }
            result = true;
            if( info<=0 )
            {
                result = false;
            }
            else
            {
                result = result & !((double)(rep.r1)<(double)(100*AP.Math.MachineEpsilon) | (double)(rep.r1)>(double)(1+1000*AP.Math.MachineEpsilon));
                result = result & !((double)(rep.rinf)<(double)(100*AP.Math.MachineEpsilon) | (double)(rep.rinf)>(double)(1+1000*AP.Math.MachineEpsilon));
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=n-1; j++)
                    {
                        v = 0.0;
                        for(i_=0; i_<=n-1;i_++)
                        {
                            v += a[i,i_]*inva[i_,j];
                        }
                        if( i==j )
                        {
                            v = v-1;
                        }
                        result = result & (double)(Math.Abs(v))<=(double)(threshold);
                    }
                }
            }
            return result;
        }


        /*************************************************************************
        Checks whether inverse is correct
        Returns True on success.
        *************************************************************************/
        private static bool hpdmatrixcheckinverse(AP.Complex[,] a,
            AP.Complex[,] inva,
            bool isupper,
            int n,
            double threshold,
            int info,
            ref matinv.matinvreport rep)
        {
            bool result = new bool();
            int i = 0;
            int j = 0;
            AP.Complex v = 0;
            int i_ = 0;

            a = (AP.Complex[,])a.Clone();
            inva = (AP.Complex[,])inva.Clone();

            for(i=0; i<=n-2; i++)
            {
                if( isupper )
                {
                    for(i_=i+1; i_<=n-1;i_++)
                    {
                        a[i_,i] = AP.Math.Conj(a[i,i_]);
                    }
                    for(i_=i+1; i_<=n-1;i_++)
                    {
                        inva[i_,i] = AP.Math.Conj(inva[i,i_]);
                    }
                }
                else
                {
                    for(i_=i+1; i_<=n-1;i_++)
                    {
                        a[i,i_] = AP.Math.Conj(a[i_,i]);
                    }
                    for(i_=i+1; i_<=n-1;i_++)
                    {
                        inva[i,i_] = AP.Math.Conj(inva[i_,i]);
                    }
                }
            }
            result = true;
            if( info<=0 )
            {
                result = false;
            }
            else
            {
                result = result & !((double)(rep.r1)<(double)(100*AP.Math.MachineEpsilon) | (double)(rep.r1)>(double)(1+1000*AP.Math.MachineEpsilon));
                result = result & !((double)(rep.rinf)<(double)(100*AP.Math.MachineEpsilon) | (double)(rep.rinf)>(double)(1+1000*AP.Math.MachineEpsilon));
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=n-1; j++)
                    {
                        v = 0.0;
                        for(i_=0; i_<=n-1;i_++)
                        {
                            v += a[i,i_]*inva[i_,j];
                        }
                        if( i==j )
                        {
                            v = v-1;
                        }
                        result = result & (double)(AP.Math.AbsComplex(v))<=(double)(threshold);
                    }
                }
            }
            return result;
        }


        /*************************************************************************
        Checks whether inversion result indicate singular matrix
        Returns True on success.
        *************************************************************************/
        private static bool rmatrixcheckinversesingular(ref double[,] inva,
            int n,
            double threshold,
            int info,
            ref matinv.matinvreport rep)
        {
            bool result = new bool();
            int i = 0;
            int j = 0;

            result = true;
            if( info!=-3 & info!=1 )
            {
                result = false;
            }
            else
            {
                result = result & !((double)(rep.r1)<(double)(0) | (double)(rep.r1)>(double)(1000*AP.Math.MachineEpsilon));
                result = result & !((double)(rep.rinf)<(double)(0) | (double)(rep.rinf)>(double)(1000*AP.Math.MachineEpsilon));
                if( info==-3 )
                {
                    for(i=0; i<=n-1; i++)
                    {
                        for(j=0; j<=n-1; j++)
                        {
                            result = result & (double)(inva[i,j])==(double)(0);
                        }
                    }
                }
            }
            return result;
        }


        /*************************************************************************
        Checks whether inverse is correct
        Returns True on success.
        *************************************************************************/
        private static bool cmatrixcheckinverse(ref AP.Complex[,] a,
            ref AP.Complex[,] inva,
            int n,
            double threshold,
            int info,
            ref matinv.matinvreport rep)
        {
            bool result = new bool();
            int i = 0;
            int j = 0;
            AP.Complex v = 0;
            int i_ = 0;

            result = true;
            if( info<=0 )
            {
                result = false;
            }
            else
            {
                result = result & !((double)(rep.r1)<(double)(100*AP.Math.MachineEpsilon) | (double)(rep.r1)>(double)(1+1000*AP.Math.MachineEpsilon));
                result = result & !((double)(rep.rinf)<(double)(100*AP.Math.MachineEpsilon) | (double)(rep.rinf)>(double)(1+1000*AP.Math.MachineEpsilon));
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=n-1; j++)
                    {
                        v = 0.0;
                        for(i_=0; i_<=n-1;i_++)
                        {
                            v += a[i,i_]*inva[i_,j];
                        }
                        if( i==j )
                        {
                            v = v-1;
                        }
                        result = result & (double)(AP.Math.AbsComplex(v))<=(double)(threshold);
                    }
                }
            }
            return result;
        }


        /*************************************************************************
        Checks whether inversion result indicate singular matrix
        Returns True on success.
        *************************************************************************/
        private static bool cmatrixcheckinversesingular(ref AP.Complex[,] inva,
            int n,
            double threshold,
            int info,
            ref matinv.matinvreport rep)
        {
            bool result = new bool();
            int i = 0;
            int j = 0;

            result = true;
            if( info!=-3 & info!=1 )
            {
                result = false;
            }
            else
            {
                result = result & !((double)(rep.r1)<(double)(0) | (double)(rep.r1)>(double)(1000*AP.Math.MachineEpsilon));
                result = result & !((double)(rep.rinf)<(double)(0) | (double)(rep.rinf)>(double)(1000*AP.Math.MachineEpsilon));
                if( info==-3 )
                {
                    for(i=0; i<=n-1; i++)
                    {
                        for(j=0; j<=n-1; j++)
                        {
                            result = result & inva[i,j]==0;
                        }
                    }
                }
            }
            return result;
        }


        /*************************************************************************
        Drops upper or lower half of the matrix - fills it by special pattern
        which may be used later to ensure that this part wasn't changed
        *************************************************************************/
        private static void rmatrixdrophalf(ref double[,] a,
            int n,
            bool droplower)
        {
            int i = 0;
            int j = 0;

            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    if( droplower & i>j | !droplower & i<j )
                    {
                        a[i,j] = 1+2*i+3*j;
                    }
                }
            }
        }


        /*************************************************************************
        Drops upper or lower half of the matrix - fills it by special pattern
        which may be used later to ensure that this part wasn't changed
        *************************************************************************/
        private static void cmatrixdrophalf(ref AP.Complex[,] a,
            int n,
            bool droplower)
        {
            int i = 0;
            int j = 0;

            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    if( droplower & i>j | !droplower & i<j )
                    {
                        a[i,j] = 1+2*i+3*j;
                    }
                }
            }
        }


        /*************************************************************************
        Real TR inverse
        *************************************************************************/
        private static void testrtrinv(int maxn,
            int passcount,
            double threshold,
            ref bool rtrerrors)
        {
            double[,] a = new double[0,0];
            double[,] b = new double[0,0];
            int n = 0;
            int pass = 0;
            int i = 0;
            int j = 0;
            int task = 0;
            bool isupper = new bool();
            bool isunit = new bool();
            double v = 0;
            bool waserrors = new bool();
            int info = 0;
            matinv.matinvreport rep = new matinv.matinvreport();
            int i_ = 0;

            waserrors = false;
            
            //
            // Test
            //
            for(n=1; n<=maxn; n++)
            {
                a = new double[n, n];
                b = new double[n, n];
                for(task=0; task<=3; task++)
                {
                    for(pass=1; pass<=passcount; pass++)
                    {
                        
                        //
                        // Determine task
                        //
                        isupper = task%2==0;
                        isunit = task/2%2==0;
                        
                        //
                        // Generate matrix
                        //
                        for(i=0; i<=n-1; i++)
                        {
                            for(j=0; j<=n-1; j++)
                            {
                                if( i==j )
                                {
                                    a[i,i] = 1+AP.Math.RandomReal();
                                }
                                else
                                {
                                    a[i,j] = 0.2*AP.Math.RandomReal()-0.1;
                                }
                                b[i,j] = a[i,j];
                            }
                        }
                        
                        //
                        // Inverse
                        //
                        matinv.rmatrixtrinverse(ref b, n, isupper, isunit, ref info, ref rep);
                        if( info<=0 )
                        {
                            rtrerrors = true;
                            return;
                        }
                        
                        //
                        // Structural test
                        //
                        if( isunit )
                        {
                            for(i=0; i<=n-1; i++)
                            {
                                rtrerrors = rtrerrors | (double)(a[i,i])!=(double)(b[i,i]);
                            }
                        }
                        if( isupper )
                        {
                            for(i=0; i<=n-1; i++)
                            {
                                for(j=0; j<=i-1; j++)
                                {
                                    rtrerrors = rtrerrors | (double)(a[i,j])!=(double)(b[i,j]);
                                }
                            }
                        }
                        else
                        {
                            for(i=0; i<=n-1; i++)
                            {
                                for(j=i+1; j<=n-1; j++)
                                {
                                    rtrerrors = rtrerrors | (double)(a[i,j])!=(double)(b[i,j]);
                                }
                            }
                        }
                        
                        //
                        // Inverse test
                        //
                        for(i=0; i<=n-1; i++)
                        {
                            for(j=0; j<=n-1; j++)
                            {
                                if( j<i & isupper | j>i & !isupper )
                                {
                                    a[i,j] = 0;
                                    b[i,j] = 0;
                                }
                            }
                        }
                        if( isunit )
                        {
                            for(i=0; i<=n-1; i++)
                            {
                                a[i,i] = 1;
                                b[i,i] = 1;
                            }
                        }
                        for(i=0; i<=n-1; i++)
                        {
                            for(j=0; j<=n-1; j++)
                            {
                                v = 0.0;
                                for(i_=0; i_<=n-1;i_++)
                                {
                                    v += a[i,i_]*b[i_,j];
                                }
                                if( j!=i )
                                {
                                    rtrerrors = rtrerrors | (double)(Math.Abs(v))>(double)(threshold);
                                }
                                else
                                {
                                    rtrerrors = rtrerrors | (double)(Math.Abs(v-1))>(double)(threshold);
                                }
                            }
                        }
                    }
                }
            }
        }


        /*************************************************************************
        Complex TR inverse
        *************************************************************************/
        private static void testctrinv(int maxn,
            int passcount,
            double threshold,
            ref bool ctrerrors)
        {
            AP.Complex[,] a = new AP.Complex[0,0];
            AP.Complex[,] b = new AP.Complex[0,0];
            int n = 0;
            int pass = 0;
            int i = 0;
            int j = 0;
            int task = 0;
            bool isupper = new bool();
            bool isunit = new bool();
            AP.Complex v = 0;
            bool waserrors = new bool();
            int info = 0;
            matinv.matinvreport rep = new matinv.matinvreport();
            int i_ = 0;

            waserrors = false;
            
            //
            // Test
            //
            for(n=1; n<=maxn; n++)
            {
                a = new AP.Complex[n, n];
                b = new AP.Complex[n, n];
                for(task=0; task<=3; task++)
                {
                    for(pass=1; pass<=passcount; pass++)
                    {
                        
                        //
                        // Determine task
                        //
                        isupper = task%2==0;
                        isunit = task/2%2==0;
                        
                        //
                        // Generate matrix
                        //
                        for(i=0; i<=n-1; i++)
                        {
                            for(j=0; j<=n-1; j++)
                            {
                                if( i==j )
                                {
                                    a[i,i].x = 1+AP.Math.RandomReal();
                                    a[i,i].y = 1+AP.Math.RandomReal();
                                }
                                else
                                {
                                    a[i,j].x = 0.2*AP.Math.RandomReal()-0.1;
                                    a[i,j].y = 0.2*AP.Math.RandomReal()-0.1;
                                }
                                b[i,j] = a[i,j];
                            }
                        }
                        
                        //
                        // Inverse
                        //
                        matinv.cmatrixtrinverse(ref b, n, isupper, isunit, ref info, ref rep);
                        if( info<=0 )
                        {
                            ctrerrors = true;
                            return;
                        }
                        
                        //
                        // Structural test
                        //
                        if( isunit )
                        {
                            for(i=0; i<=n-1; i++)
                            {
                                ctrerrors = ctrerrors | a[i,i]!=b[i,i];
                            }
                        }
                        if( isupper )
                        {
                            for(i=0; i<=n-1; i++)
                            {
                                for(j=0; j<=i-1; j++)
                                {
                                    ctrerrors = ctrerrors | a[i,j]!=b[i,j];
                                }
                            }
                        }
                        else
                        {
                            for(i=0; i<=n-1; i++)
                            {
                                for(j=i+1; j<=n-1; j++)
                                {
                                    ctrerrors = ctrerrors | a[i,j]!=b[i,j];
                                }
                            }
                        }
                        
                        //
                        // Inverse test
                        //
                        for(i=0; i<=n-1; i++)
                        {
                            for(j=0; j<=n-1; j++)
                            {
                                if( j<i & isupper | j>i & !isupper )
                                {
                                    a[i,j] = 0;
                                    b[i,j] = 0;
                                }
                            }
                        }
                        if( isunit )
                        {
                            for(i=0; i<=n-1; i++)
                            {
                                a[i,i] = 1;
                                b[i,i] = 1;
                            }
                        }
                        for(i=0; i<=n-1; i++)
                        {
                            for(j=0; j<=n-1; j++)
                            {
                                v = 0.0;
                                for(i_=0; i_<=n-1;i_++)
                                {
                                    v += a[i,i_]*b[i_,j];
                                }
                                if( j!=i )
                                {
                                    ctrerrors = ctrerrors | (double)(AP.Math.AbsComplex(v))>(double)(threshold);
                                }
                                else
                                {
                                    ctrerrors = ctrerrors | (double)(AP.Math.AbsComplex(v-1))>(double)(threshold);
                                }
                            }
                        }
                    }
                }
            }
        }


        /*************************************************************************
        Real test
        *************************************************************************/
        private static void testrinv(int maxn,
            int passcount,
            double threshold,
            ref bool rerrors)
        {
            double[,] a = new double[0,0];
            double[,] lua = new double[0,0];
            double[,] inva = new double[0,0];
            double[,] invlua = new double[0,0];
            int[] p = new int[0];
            int i = 0;
            int j = 0;
            int k = 0;
            int n = 0;
            int pass = 0;
            int taskkind = 0;
            double v = 0;
            int info = 0;
            matinv.matinvreport rep = new matinv.matinvreport();
            int i_ = 0;

            
            //
            // General square matrices:
            // * test general solvers
            // * test least squares solver
            //
            for(pass=1; pass<=passcount; pass++)
            {
                for(n=1; n<=maxn; n++)
                {
                    
                    //
                    // ********************************************************
                    // WELL CONDITIONED TASKS
                    // ability to find correct solution is tested
                    // ********************************************************
                    //
                    // 1. generate random well conditioned matrix A.
                    // 2. generate random solution vector xe
                    // 3. generate right part b=A*xe
                    // 4. test different methods on original A
                    //
                    matgen.rmatrixrndcond(n, 1000, ref a);
                    rmatrixmakeacopy(ref a, n, n, ref lua);
                    trfac.rmatrixlu(ref lua, n, n, ref p);
                    rmatrixmakeacopy(ref a, n, n, ref inva);
                    rmatrixmakeacopy(ref lua, n, n, ref invlua);
                    info = 0;
                    unsetrep(ref rep);
                    matinv.rmatrixinverse(ref inva, n, ref info, ref rep);
                    rerrors = rerrors | !rmatrixcheckinverse(ref a, ref inva, n, threshold, info, ref rep);
                    info = 0;
                    unsetrep(ref rep);
                    matinv.rmatrixluinverse(ref invlua, ref p, n, ref info, ref rep);
                    rerrors = rerrors | !rmatrixcheckinverse(ref a, ref invlua, n, threshold, info, ref rep);
                    
                    //
                    // ********************************************************
                    // EXACTLY SINGULAR MATRICES
                    // ability to detect singularity is tested
                    // ********************************************************
                    //
                    // 1. generate different types of singular matrices:
                    //    * zero
                    //    * with zero columns
                    //    * with zero rows
                    //    * with equal rows/columns
                    // 2. test different methods
                    //
                    for(taskkind=0; taskkind<=4; taskkind++)
                    {
                        unset2d(ref a);
                        if( taskkind==0 )
                        {
                            
                            //
                            // all zeros
                            //
                            a = new double[n, n];
                            for(i=0; i<=n-1; i++)
                            {
                                for(j=0; j<=n-1; j++)
                                {
                                    a[i,j] = 0;
                                }
                            }
                        }
                        if( taskkind==1 )
                        {
                            
                            //
                            // there is zero column
                            //
                            a = new double[n, n];
                            for(i=0; i<=n-1; i++)
                            {
                                for(j=0; j<=n-1; j++)
                                {
                                    a[i,j] = 2*AP.Math.RandomReal()-1;
                                }
                            }
                            k = AP.Math.RandomInteger(n);
                            for(i_=0; i_<=n-1;i_++)
                            {
                                a[i_,k] = 0*a[i_,k];
                            }
                        }
                        if( taskkind==2 )
                        {
                            
                            //
                            // there is zero row
                            //
                            a = new double[n, n];
                            for(i=0; i<=n-1; i++)
                            {
                                for(j=0; j<=n-1; j++)
                                {
                                    a[i,j] = 2*AP.Math.RandomReal()-1;
                                }
                            }
                            k = AP.Math.RandomInteger(n);
                            for(i_=0; i_<=n-1;i_++)
                            {
                                a[k,i_] = 0*a[k,i_];
                            }
                        }
                        if( taskkind==3 )
                        {
                            
                            //
                            // equal columns
                            //
                            if( n<2 )
                            {
                                continue;
                            }
                            a = new double[n, n];
                            for(i=0; i<=n-1; i++)
                            {
                                for(j=0; j<=n-1; j++)
                                {
                                    a[i,j] = 2*AP.Math.RandomReal()-1;
                                }
                            }
                            k = 1+AP.Math.RandomInteger(n-1);
                            for(i_=0; i_<=n-1;i_++)
                            {
                                a[i_,0] = a[i_,k];
                            }
                        }
                        if( taskkind==4 )
                        {
                            
                            //
                            // equal rows
                            //
                            if( n<2 )
                            {
                                continue;
                            }
                            a = new double[n, n];
                            for(i=0; i<=n-1; i++)
                            {
                                for(j=0; j<=n-1; j++)
                                {
                                    a[i,j] = 2*AP.Math.RandomReal()-1;
                                }
                            }
                            k = 1+AP.Math.RandomInteger(n-1);
                            for(i_=0; i_<=n-1;i_++)
                            {
                                a[0,i_] = a[k,i_];
                            }
                        }
                        rmatrixmakeacopy(ref a, n, n, ref lua);
                        trfac.rmatrixlu(ref lua, n, n, ref p);
                        info = 0;
                        unsetrep(ref rep);
                        matinv.rmatrixinverse(ref a, n, ref info, ref rep);
                        rerrors = rerrors | !rmatrixcheckinversesingular(ref a, n, threshold, info, ref rep);
                        info = 0;
                        unsetrep(ref rep);
                        matinv.rmatrixluinverse(ref lua, ref p, n, ref info, ref rep);
                        rerrors = rerrors | !rmatrixcheckinversesingular(ref lua, n, threshold, info, ref rep);
                    }
                }
            }
        }


        /*************************************************************************
        Complex test
        *************************************************************************/
        private static void testcinv(int maxn,
            int passcount,
            double threshold,
            ref bool cerrors)
        {
            AP.Complex[,] a = new AP.Complex[0,0];
            AP.Complex[,] lua = new AP.Complex[0,0];
            AP.Complex[,] inva = new AP.Complex[0,0];
            AP.Complex[,] invlua = new AP.Complex[0,0];
            int[] p = new int[0];
            int i = 0;
            int j = 0;
            int k = 0;
            int n = 0;
            int pass = 0;
            int taskkind = 0;
            double v = 0;
            int info = 0;
            matinv.matinvreport rep = new matinv.matinvreport();
            int i_ = 0;

            
            //
            // General square matrices:
            // * test general solvers
            // * test least squares solver
            //
            for(pass=1; pass<=passcount; pass++)
            {
                for(n=1; n<=maxn; n++)
                {
                    
                    //
                    // ********************************************************
                    // WELL CONDITIONED TASKS
                    // ability to find correct solution is tested
                    // ********************************************************
                    //
                    // 1. generate random well conditioned matrix A.
                    // 2. generate random solution vector xe
                    // 3. generate right part b=A*xe
                    // 4. test different methods on original A
                    //
                    matgen.cmatrixrndcond(n, 1000, ref a);
                    cmatrixmakeacopy(ref a, n, n, ref lua);
                    trfac.cmatrixlu(ref lua, n, n, ref p);
                    cmatrixmakeacopy(ref a, n, n, ref inva);
                    cmatrixmakeacopy(ref lua, n, n, ref invlua);
                    info = 0;
                    unsetrep(ref rep);
                    matinv.cmatrixinverse(ref inva, n, ref info, ref rep);
                    cerrors = cerrors | !cmatrixcheckinverse(ref a, ref inva, n, threshold, info, ref rep);
                    info = 0;
                    unsetrep(ref rep);
                    matinv.cmatrixluinverse(ref invlua, ref p, n, ref info, ref rep);
                    cerrors = cerrors | !cmatrixcheckinverse(ref a, ref invlua, n, threshold, info, ref rep);
                    
                    //
                    // ********************************************************
                    // EXACTLY SINGULAR MATRICES
                    // ability to detect singularity is tested
                    // ********************************************************
                    //
                    // 1. generate different types of singular matrices:
                    //    * zero
                    //    * with zero columns
                    //    * with zero rows
                    //    * with equal rows/columns
                    // 2. test different methods
                    //
                    for(taskkind=0; taskkind<=4; taskkind++)
                    {
                        cunset2d(ref a);
                        if( taskkind==0 )
                        {
                            
                            //
                            // all zeros
                            //
                            a = new AP.Complex[n, n];
                            for(i=0; i<=n-1; i++)
                            {
                                for(j=0; j<=n-1; j++)
                                {
                                    a[i,j] = 0;
                                }
                            }
                        }
                        if( taskkind==1 )
                        {
                            
                            //
                            // there is zero column
                            //
                            a = new AP.Complex[n, n];
                            for(i=0; i<=n-1; i++)
                            {
                                for(j=0; j<=n-1; j++)
                                {
                                    a[i,j].x = 2*AP.Math.RandomReal()-1;
                                    a[i,j].y = 2*AP.Math.RandomReal()-1;
                                }
                            }
                            k = AP.Math.RandomInteger(n);
                            for(i_=0; i_<=n-1;i_++)
                            {
                                a[i_,k] = 0*a[i_,k];
                            }
                        }
                        if( taskkind==2 )
                        {
                            
                            //
                            // there is zero row
                            //
                            a = new AP.Complex[n, n];
                            for(i=0; i<=n-1; i++)
                            {
                                for(j=0; j<=n-1; j++)
                                {
                                    a[i,j].x = 2*AP.Math.RandomReal()-1;
                                    a[i,j].y = 2*AP.Math.RandomReal()-1;
                                }
                            }
                            k = AP.Math.RandomInteger(n);
                            for(i_=0; i_<=n-1;i_++)
                            {
                                a[k,i_] = 0*a[k,i_];
                            }
                        }
                        if( taskkind==3 )
                        {
                            
                            //
                            // equal columns
                            //
                            if( n<2 )
                            {
                                continue;
                            }
                            a = new AP.Complex[n, n];
                            for(i=0; i<=n-1; i++)
                            {
                                for(j=0; j<=n-1; j++)
                                {
                                    a[i,j].x = 2*AP.Math.RandomReal()-1;
                                    a[i,j].y = 2*AP.Math.RandomReal()-1;
                                }
                            }
                            k = 1+AP.Math.RandomInteger(n-1);
                            for(i_=0; i_<=n-1;i_++)
                            {
                                a[i_,0] = a[i_,k];
                            }
                        }
                        if( taskkind==4 )
                        {
                            
                            //
                            // equal rows
                            //
                            if( n<2 )
                            {
                                continue;
                            }
                            a = new AP.Complex[n, n];
                            for(i=0; i<=n-1; i++)
                            {
                                for(j=0; j<=n-1; j++)
                                {
                                    a[i,j].x = 2*AP.Math.RandomReal()-1;
                                    a[i,j].y = 2*AP.Math.RandomReal()-1;
                                }
                            }
                            k = 1+AP.Math.RandomInteger(n-1);
                            for(i_=0; i_<=n-1;i_++)
                            {
                                a[0,i_] = a[k,i_];
                            }
                        }
                        cmatrixmakeacopy(ref a, n, n, ref lua);
                        trfac.cmatrixlu(ref lua, n, n, ref p);
                        info = 0;
                        unsetrep(ref rep);
                        matinv.cmatrixinverse(ref a, n, ref info, ref rep);
                        cerrors = cerrors | !cmatrixcheckinversesingular(ref a, n, threshold, info, ref rep);
                        info = 0;
                        unsetrep(ref rep);
                        matinv.cmatrixluinverse(ref lua, ref p, n, ref info, ref rep);
                        cerrors = cerrors | !cmatrixcheckinversesingular(ref lua, n, threshold, info, ref rep);
                    }
                }
            }
        }


        /*************************************************************************
        SPD test
        *************************************************************************/
        private static void testspdinv(int maxn,
            int passcount,
            double threshold,
            ref bool spderrors)
        {
            double[,] a = new double[0,0];
            double[,] cha = new double[0,0];
            double[,] inva = new double[0,0];
            double[,] invcha = new double[0,0];
            bool isupper = new bool();
            int i = 0;
            int j = 0;
            int k = 0;
            int n = 0;
            int pass = 0;
            int taskkind = 0;
            double v = 0;
            int info = 0;
            matinv.matinvreport rep = new matinv.matinvreport();
            int i_ = 0;

            
            //
            // General square matrices:
            // * test general solvers
            // * test least squares solver
            //
            for(pass=1; pass<=passcount; pass++)
            {
                for(n=1; n<=maxn; n++)
                {
                    isupper = (double)(AP.Math.RandomReal())>(double)(0.5);
                    
                    //
                    // ********************************************************
                    // WELL CONDITIONED TASKS
                    // ability to find correct solution is tested
                    // ********************************************************
                    //
                    // 1. generate random well conditioned matrix A.
                    // 2. generate random solution vector xe
                    // 3. generate right part b=A*xe
                    // 4. test different methods on original A
                    //
                    matgen.spdmatrixrndcond(n, 1000, ref a);
                    rmatrixdrophalf(ref a, n, isupper);
                    rmatrixmakeacopy(ref a, n, n, ref cha);
                    if( !trfac.spdmatrixcholesky(ref cha, n, isupper) )
                    {
                        continue;
                    }
                    rmatrixmakeacopy(ref a, n, n, ref inva);
                    rmatrixmakeacopy(ref cha, n, n, ref invcha);
                    info = 0;
                    unsetrep(ref rep);
                    matinv.spdmatrixinverse(ref inva, n, isupper, ref info, ref rep);
                    spderrors = spderrors | !spdmatrixcheckinverse(a, inva, isupper, n, threshold, info, ref rep);
                    info = 0;
                    unsetrep(ref rep);
                    matinv.spdmatrixcholeskyinverse(ref invcha, n, isupper, ref info, ref rep);
                    spderrors = spderrors | !spdmatrixcheckinverse(a, invcha, isupper, n, threshold, info, ref rep);
                    
                    //
                    // ********************************************************
                    // EXACTLY SINGULAR MATRICES
                    // ability to detect singularity is tested
                    // ********************************************************
                    //
                    // 1. generate different types of singular matrices:
                    //    * zero
                    //    * with zero columns
                    //    * with zero rows
                    // 2. test different methods
                    //
                    for(taskkind=0; taskkind<=2; taskkind++)
                    {
                        unset2d(ref a);
                        if( taskkind==0 )
                        {
                            
                            //
                            // all zeros
                            //
                            a = new double[n, n];
                            for(i=0; i<=n-1; i++)
                            {
                                for(j=0; j<=n-1; j++)
                                {
                                    a[i,j] = 0;
                                }
                            }
                        }
                        if( taskkind==1 )
                        {
                            
                            //
                            // there is zero column
                            //
                            a = new double[n, n];
                            for(i=0; i<=n-1; i++)
                            {
                                for(j=0; j<=n-1; j++)
                                {
                                    a[i,j] = 2*AP.Math.RandomReal()-1;
                                }
                            }
                            k = AP.Math.RandomInteger(n);
                            for(i_=0; i_<=n-1;i_++)
                            {
                                a[i_,k] = 0*a[i_,k];
                            }
                        }
                        if( taskkind==2 )
                        {
                            
                            //
                            // there is zero row
                            //
                            a = new double[n, n];
                            for(i=0; i<=n-1; i++)
                            {
                                for(j=0; j<=n-1; j++)
                                {
                                    a[i,j] = 2*AP.Math.RandomReal()-1;
                                }
                            }
                            k = AP.Math.RandomInteger(n);
                            for(i_=0; i_<=n-1;i_++)
                            {
                                a[k,i_] = 0*a[k,i_];
                            }
                        }
                        info = 0;
                        unsetrep(ref rep);
                        matinv.spdmatrixcholeskyinverse(ref a, n, isupper, ref info, ref rep);
                        if( info!=-3 & info!=1 )
                        {
                            spderrors = true;
                        }
                        else
                        {
                            spderrors = spderrors | (double)(rep.r1)<(double)(0) | (double)(rep.r1)>(double)(1000*AP.Math.MachineEpsilon);
                            spderrors = spderrors | (double)(rep.rinf)<(double)(0) | (double)(rep.rinf)>(double)(1000*AP.Math.MachineEpsilon);
                        }
                    }
                }
            }
        }


        /*************************************************************************
        HPD test
        *************************************************************************/
        private static void testhpdinv(int maxn,
            int passcount,
            double threshold,
            ref bool hpderrors)
        {
            AP.Complex[,] a = new AP.Complex[0,0];
            AP.Complex[,] cha = new AP.Complex[0,0];
            AP.Complex[,] inva = new AP.Complex[0,0];
            AP.Complex[,] invcha = new AP.Complex[0,0];
            bool isupper = new bool();
            int i = 0;
            int j = 0;
            int k = 0;
            int n = 0;
            int pass = 0;
            int taskkind = 0;
            AP.Complex v = 0;
            int info = 0;
            matinv.matinvreport rep = new matinv.matinvreport();
            int i_ = 0;

            
            //
            // General square matrices:
            // * test general solvers
            // * test least squares solver
            //
            for(pass=1; pass<=passcount; pass++)
            {
                for(n=1; n<=maxn; n++)
                {
                    isupper = (double)(AP.Math.RandomReal())>(double)(0.5);
                    
                    //
                    // ********************************************************
                    // WELL CONDITIONED TASKS
                    // ability to find correct solution is tested
                    // ********************************************************
                    //
                    // 1. generate random well conditioned matrix A.
                    // 2. generate random solution vector xe
                    // 3. generate right part b=A*xe
                    // 4. test different methods on original A
                    //
                    matgen.hpdmatrixrndcond(n, 1000, ref a);
                    cmatrixdrophalf(ref a, n, isupper);
                    cmatrixmakeacopy(ref a, n, n, ref cha);
                    if( !trfac.hpdmatrixcholesky(ref cha, n, isupper) )
                    {
                        continue;
                    }
                    cmatrixmakeacopy(ref a, n, n, ref inva);
                    cmatrixmakeacopy(ref cha, n, n, ref invcha);
                    info = 0;
                    unsetrep(ref rep);
                    matinv.hpdmatrixinverse(ref inva, n, isupper, ref info, ref rep);
                    hpderrors = hpderrors | !hpdmatrixcheckinverse(a, inva, isupper, n, threshold, info, ref rep);
                    info = 0;
                    unsetrep(ref rep);
                    matinv.hpdmatrixcholeskyinverse(ref invcha, n, isupper, ref info, ref rep);
                    hpderrors = hpderrors | !hpdmatrixcheckinverse(a, invcha, isupper, n, threshold, info, ref rep);
                    
                    //
                    // ********************************************************
                    // EXACTLY SINGULAR MATRICES
                    // ability to detect singularity is tested
                    // ********************************************************
                    //
                    // 1. generate different types of singular matrices:
                    //    * zero
                    //    * with zero columns
                    //    * with zero rows
                    // 2. test different methods
                    //
                    for(taskkind=0; taskkind<=2; taskkind++)
                    {
                        cunset2d(ref a);
                        if( taskkind==0 )
                        {
                            
                            //
                            // all zeros
                            //
                            a = new AP.Complex[n, n];
                            for(i=0; i<=n-1; i++)
                            {
                                for(j=0; j<=n-1; j++)
                                {
                                    a[i,j] = 0;
                                }
                            }
                        }
                        if( taskkind==1 )
                        {
                            
                            //
                            // there is zero column
                            //
                            a = new AP.Complex[n, n];
                            for(i=0; i<=n-1; i++)
                            {
                                for(j=0; j<=n-1; j++)
                                {
                                    a[i,j].x = 2*AP.Math.RandomReal()-1;
                                    a[i,j].y = 2*AP.Math.RandomReal()-1;
                                }
                            }
                            k = AP.Math.RandomInteger(n);
                            for(i_=0; i_<=n-1;i_++)
                            {
                                a[i_,k] = 0*a[i_,k];
                            }
                            for(i_=0; i_<=n-1;i_++)
                            {
                                a[k,i_] = 0*a[k,i_];
                            }
                        }
                        if( taskkind==2 )
                        {
                            
                            //
                            // there is zero row
                            //
                            a = new AP.Complex[n, n];
                            for(i=0; i<=n-1; i++)
                            {
                                for(j=0; j<=n-1; j++)
                                {
                                    a[i,j].x = 2*AP.Math.RandomReal()-1;
                                    a[i,j].y = 2*AP.Math.RandomReal()-1;
                                }
                            }
                            k = AP.Math.RandomInteger(n);
                            for(i_=0; i_<=n-1;i_++)
                            {
                                a[k,i_] = 0*a[k,i_];
                            }
                            for(i_=0; i_<=n-1;i_++)
                            {
                                a[i_,k] = 0*a[i_,k];
                            }
                        }
                        info = 0;
                        unsetrep(ref rep);
                        matinv.hpdmatrixcholeskyinverse(ref a, n, isupper, ref info, ref rep);
                        if( info!=-3 & info!=1 )
                        {
                            hpderrors = true;
                        }
                        else
                        {
                            hpderrors = hpderrors | (double)(rep.r1)<(double)(0) | (double)(rep.r1)>(double)(1000*AP.Math.MachineEpsilon);
                            hpderrors = hpderrors | (double)(rep.rinf)<(double)(0) | (double)(rep.rinf)>(double)(1000*AP.Math.MachineEpsilon);
                        }
                    }
                }
            }
        }


        /*************************************************************************
        Unsets real matrix
        *************************************************************************/
        private static void unset2d(ref double[,] x)
        {
            x = new double[1, 1];
            x[0,0] = 2*AP.Math.RandomReal()-1;
        }


        /*************************************************************************
        Unsets real vector
        *************************************************************************/
        private static void unset1d(ref double[] x)
        {
            x = new double[1];
            x[0] = 2*AP.Math.RandomReal()-1;
        }


        /*************************************************************************
        Unsets real matrix
        *************************************************************************/
        private static void cunset2d(ref AP.Complex[,] x)
        {
            x = new AP.Complex[1, 1];
            x[0,0] = 2*AP.Math.RandomReal()-1;
        }


        /*************************************************************************
        Unsets real vector
        *************************************************************************/
        private static void cunset1d(ref AP.Complex[] x)
        {
            x = new AP.Complex[1];
            x[0] = 2*AP.Math.RandomReal()-1;
        }


        /*************************************************************************
        Unsets report
        *************************************************************************/
        private static void unsetrep(ref matinv.matinvreport r)
        {
            r.r1 = -1;
            r.rinf = -1;
        }


        /*************************************************************************
        Silent unit test
        *************************************************************************/
        public static bool testmatinvunit_test_silent()
        {
            bool result = new bool();

            result = testmatinv(true);
            return result;
        }


        /*************************************************************************
        Unit test
        *************************************************************************/
        public static bool testmatinvunit_test()
        {
            bool result = new bool();

            result = testmatinv(false);
            return result;
        }
    }
}
