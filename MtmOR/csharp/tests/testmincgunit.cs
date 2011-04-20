
using System;

namespace alglib
{
    public class testmincgunit
    {
        public static bool testmincg(bool silent)
        {
            bool result = new bool();
            bool waserrors = new bool();
            bool referror = new bool();
            bool eqerror = new bool();
            bool linerror1 = new bool();
            bool linerror2 = new bool();
            bool converror = new bool();
            bool othererrors = new bool();
            int n = 0;
            double[] x = new double[0];
            double[] xe = new double[0];
            double[] b = new double[0];
            double[] xlast = new double[0];
            double fprev = 0;
            double xprev = 0;
            double stpmax = 0;
            int i = 0;
            int j = 0;
            double v = 0;
            double[,] a = new double[0,0];
            mincg.mincgstate state = new mincg.mincgstate();
            mincg.mincgreport rep = new mincg.mincgreport();
            int cgtype = 0;
            int i_ = 0;

            waserrors = false;
            referror = false;
            linerror1 = false;
            linerror2 = false;
            eqerror = false;
            converror = false;
            othererrors = false;
            for(cgtype=0; cgtype<=1; cgtype++)
            {
                
                //
                // Reference problem
                //
                x = new double[2+1];
                n = 3;
                x[0] = 100*AP.Math.RandomReal()-50;
                x[1] = 100*AP.Math.RandomReal()-50;
                x[2] = 100*AP.Math.RandomReal()-50;
                mincg.mincgcreate(n, ref x, ref state);
                mincg.mincgsetcgtype(ref state, cgtype);
                while( mincg.mincgiteration(ref state) )
                {
                    state.f = AP.Math.Sqr(state.x[0]-2)+AP.Math.Sqr(state.x[1])+AP.Math.Sqr(state.x[2]-state.x[0]);
                    state.g[0] = 2*(state.x[0]-2)+2*(state.x[0]-state.x[2]);
                    state.g[1] = 2*state.x[1];
                    state.g[2] = 2*(state.x[2]-state.x[0]);
                }
                mincg.mincgresults(ref state, ref x, ref rep);
                referror = referror | rep.terminationtype<=0 | (double)(Math.Abs(x[0]-2))>(double)(0.001) | (double)(Math.Abs(x[1]))>(double)(0.001) | (double)(Math.Abs(x[2]-2))>(double)(0.001);
                
                //
                // 1D problem #1
                //
                x = new double[0+1];
                n = 1;
                x[0] = 100*AP.Math.RandomReal()-50;
                mincg.mincgcreate(n, ref x, ref state);
                mincg.mincgsetcgtype(ref state, cgtype);
                while( mincg.mincgiteration(ref state) )
                {
                    state.f = -Math.Cos(state.x[0]);
                    state.g[0] = Math.Sin(state.x[0]);
                }
                mincg.mincgresults(ref state, ref x, ref rep);
                linerror1 = linerror1 | rep.terminationtype<=0 | (double)(Math.Abs(x[0]/Math.PI-(int)Math.Round(x[0]/Math.PI)))>(double)(0.001);
                
                //
                // 1D problem #2
                //
                x = new double[0+1];
                n = 1;
                x[0] = 100*AP.Math.RandomReal()-50;
                mincg.mincgcreate(n, ref x, ref state);
                mincg.mincgsetcgtype(ref state, cgtype);
                while( mincg.mincgiteration(ref state) )
                {
                    state.f = AP.Math.Sqr(state.x[0])/(1+AP.Math.Sqr(state.x[0]));
                    state.g[0] = (2*state.x[0]*(1+AP.Math.Sqr(state.x[0]))-AP.Math.Sqr(state.x[0])*2*state.x[0])/AP.Math.Sqr(1+AP.Math.Sqr(state.x[0]));
                }
                mincg.mincgresults(ref state, ref x, ref rep);
                linerror2 = linerror2 | rep.terminationtype<=0 | (double)(Math.Abs(x[0]))>(double)(0.001);
                
                //
                // Linear equations
                //
                for(n=1; n<=10; n++)
                {
                    
                    //
                    // Prepare task
                    //
                    a = new double[n-1+1, n-1+1];
                    x = new double[n-1+1];
                    xe = new double[n-1+1];
                    b = new double[n-1+1];
                    for(i=0; i<=n-1; i++)
                    {
                        xe[i] = 2*AP.Math.RandomReal()-1;
                    }
                    for(i=0; i<=n-1; i++)
                    {
                        for(j=0; j<=n-1; j++)
                        {
                            a[i,j] = 2*AP.Math.RandomReal()-1;
                        }
                        a[i,i] = a[i,i]+3*Math.Sign(a[i,i]);
                    }
                    for(i=0; i<=n-1; i++)
                    {
                        v = 0.0;
                        for(i_=0; i_<=n-1;i_++)
                        {
                            v += a[i,i_]*xe[i_];
                        }
                        b[i] = v;
                    }
                    
                    //
                    // Solve task
                    //
                    for(i=0; i<=n-1; i++)
                    {
                        x[i] = 2*AP.Math.RandomReal()-1;
                    }
                    mincg.mincgcreate(n, ref x, ref state);
                    mincg.mincgsetcgtype(ref state, cgtype);
                    while( mincg.mincgiteration(ref state) )
                    {
                        state.f = 0;
                        for(i=0; i<=n-1; i++)
                        {
                            state.g[i] = 0;
                        }
                        for(i=0; i<=n-1; i++)
                        {
                            v = 0.0;
                            for(i_=0; i_<=n-1;i_++)
                            {
                                v += a[i,i_]*state.x[i_];
                            }
                            state.f = state.f+AP.Math.Sqr(v-b[i]);
                            for(j=0; j<=n-1; j++)
                            {
                                state.g[j] = state.g[j]+2*(v-b[i])*a[i,j];
                            }
                        }
                    }
                    mincg.mincgresults(ref state, ref x, ref rep);
                    eqerror = eqerror | rep.terminationtype<=0;
                    for(i=0; i<=n-1; i++)
                    {
                        eqerror = eqerror | (double)(Math.Abs(x[i]-xe[i]))>(double)(0.001);
                    }
                }
                
                //
                // Testing convergence properties
                //
                x = new double[2+1];
                n = 3;
                for(i=0; i<=2; i++)
                {
                    x[i] = 6*AP.Math.RandomReal()-3;
                }
                mincg.mincgcreate(n, ref x, ref state);
                mincg.mincgsetcond(ref state, 0.001, 0.0, 0.0, 0);
                mincg.mincgsetcgtype(ref state, cgtype);
                while( mincg.mincgiteration(ref state) )
                {
                    testfunc3(ref state);
                }
                mincg.mincgresults(ref state, ref x, ref rep);
                converror = converror | rep.terminationtype!=4;
                for(i=0; i<=2; i++)
                {
                    x[i] = 6*AP.Math.RandomReal()-3;
                }
                mincg.mincgcreate(n, ref x, ref state);
                mincg.mincgsetcond(ref state, 0.0, 0.001, 0.0, 0);
                mincg.mincgsetcgtype(ref state, cgtype);
                while( mincg.mincgiteration(ref state) )
                {
                    testfunc3(ref state);
                }
                mincg.mincgresults(ref state, ref x, ref rep);
                converror = converror | rep.terminationtype!=1;
                for(i=0; i<=2; i++)
                {
                    x[i] = 6*AP.Math.RandomReal()-3;
                }
                mincg.mincgcreate(n, ref x, ref state);
                mincg.mincgsetcond(ref state, 0.0, 0.0, 0.001, 0);
                mincg.mincgsetcgtype(ref state, cgtype);
                while( mincg.mincgiteration(ref state) )
                {
                    testfunc3(ref state);
                }
                mincg.mincgresults(ref state, ref x, ref rep);
                converror = converror | rep.terminationtype!=2;
                for(i=0; i<=2; i++)
                {
                    x[i] = 2*AP.Math.RandomReal()-1;
                }
                mincg.mincgcreate(n, ref x, ref state);
                mincg.mincgsetcond(ref state, 0.0, 0.0, 0.0, 10);
                mincg.mincgsetcgtype(ref state, cgtype);
                while( mincg.mincgiteration(ref state) )
                {
                    testfunc3(ref state);
                }
                mincg.mincgresults(ref state, ref x, ref rep);
                converror = converror | !(rep.terminationtype==5 & rep.iterationscount==10 | rep.terminationtype==7);
                
                //
                // Other properties:
                // 1. test reports (F should form monotone sequence)
                // 2. test maximum step
                //
                n = 50;
                x = new double[n];
                xlast = new double[n];
                for(i=0; i<=n-1; i++)
                {
                    x[i] = 1;
                }
                mincg.mincgcreate(n, ref x, ref state);
                mincg.mincgsetcond(ref state, 0, 0, 0, 100);
                mincg.mincgsetxrep(ref state, true);
                fprev = AP.Math.MaxRealNumber;
                while( mincg.mincgiteration(ref state) )
                {
                    if( state.needfg )
                    {
                        state.f = 0;
                        for(i=0; i<=n-1; i++)
                        {
                            state.f = state.f+AP.Math.Sqr((1+i)*state.x[i]);
                            state.g[i] = 2*(1+i)*state.x[i];
                        }
                    }
                    if( state.xupdated )
                    {
                        othererrors = othererrors | (double)(state.f)>(double)(fprev);
                        if( (double)(fprev)==(double)(AP.Math.MaxRealNumber) )
                        {
                            for(i=0; i<=n-1; i++)
                            {
                                othererrors = othererrors | (double)(state.x[i])!=(double)(x[i]);
                            }
                        }
                        fprev = state.f;
                        for(i_=0; i_<=n-1;i_++)
                        {
                            xlast[i_] = state.x[i_];
                        }
                    }
                }
                mincg.mincgresults(ref state, ref x, ref rep);
                for(i=0; i<=n-1; i++)
                {
                    othererrors = othererrors | (double)(x[i])!=(double)(xlast[i]);
                }
                n = 1;
                x = new double[n];
                x[0] = 100;
                stpmax = 0.05+0.05*AP.Math.RandomReal();
                mincg.mincgcreate(n, ref x, ref state);
                mincg.mincgsetcond(ref state, 1.0E-9, 0, 0, 0);
                mincg.mincgsetstpmax(ref state, stpmax);
                mincg.mincgsetxrep(ref state, true);
                xprev = x[0];
                while( mincg.mincgiteration(ref state) )
                {
                    if( state.needfg )
                    {
                        state.f = Math.Exp(state.x[0])+Math.Exp(-state.x[0]);
                        state.g[0] = Math.Exp(state.x[0])-Math.Exp(-state.x[0]);
                        othererrors = othererrors | (double)(Math.Abs(state.x[0]-xprev))>(double)((1+Math.Sqrt(AP.Math.MachineEpsilon))*stpmax);
                    }
                    if( state.xupdated )
                    {
                        othererrors = othererrors | (double)(Math.Abs(state.x[0]-xprev))>(double)((1+Math.Sqrt(AP.Math.MachineEpsilon))*stpmax);
                        xprev = state.x[0];
                    }
                }
            }
            
            //
            // end
            //
            waserrors = referror | eqerror | linerror1 | linerror2 | converror | othererrors;
            if( !silent )
            {
                System.Console.Write("TESTING CG OPTIMIZATION");
                System.Console.WriteLine();
                System.Console.Write("REFERENCE PROBLEM:                        ");
                if( referror )
                {
                    System.Console.Write("FAILED");
                    System.Console.WriteLine();
                }
                else
                {
                    System.Console.Write("OK");
                    System.Console.WriteLine();
                }
                System.Console.Write("LIN-1 PROBLEM:                            ");
                if( linerror1 )
                {
                    System.Console.Write("FAILED");
                    System.Console.WriteLine();
                }
                else
                {
                    System.Console.Write("OK");
                    System.Console.WriteLine();
                }
                System.Console.Write("LIN-2 PROBLEM:                            ");
                if( linerror2 )
                {
                    System.Console.Write("FAILED");
                    System.Console.WriteLine();
                }
                else
                {
                    System.Console.Write("OK");
                    System.Console.WriteLine();
                }
                System.Console.Write("LINEAR EQUATIONS:                         ");
                if( eqerror )
                {
                    System.Console.Write("FAILED");
                    System.Console.WriteLine();
                }
                else
                {
                    System.Console.Write("OK");
                    System.Console.WriteLine();
                }
                System.Console.Write("CONVERGENCE PROPERTIES:                   ");
                if( converror )
                {
                    System.Console.Write("FAILED");
                    System.Console.WriteLine();
                }
                else
                {
                    System.Console.Write("OK");
                    System.Console.WriteLine();
                }
                System.Console.Write("OTHER PROPERTIES:                         ");
                if( othererrors )
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
                System.Console.WriteLine();
                System.Console.WriteLine();
            }
            result = !waserrors;
            return result;
        }


        /*************************************************************************
        Calculate test function #1
        *************************************************************************/
        private static void testfunc1(ref mincg.mincgstate state)
        {
            if( (double)(state.x[0])<(double)(100) )
            {
                state.f = AP.Math.Sqr(Math.Exp(state.x[0])-2)+AP.Math.Sqr(state.x[1])+AP.Math.Sqr(state.x[2]-state.x[0]);
                state.g[0] = 2*(Math.Exp(state.x[0])-2)*Math.Exp(state.x[0])+2*(state.x[0]-state.x[2]);
                state.g[1] = 2*state.x[1];
                state.g[2] = 2*(state.x[2]-state.x[0]);
            }
            else
            {
                state.f = Math.Sqrt(AP.Math.MaxRealNumber);
                state.g[0] = Math.Sqrt(AP.Math.MaxRealNumber);
                state.g[1] = 0;
                state.g[2] = 0;
            }
        }


        /*************************************************************************
        Calculate test function #2

        Simple variation of #1, much more nonlinear, which makes unlikely premature
        convergence of algorithm .
        *************************************************************************/
        private static void testfunc2(ref mincg.mincgstate state)
        {
            if( (double)(state.x[0])<(double)(100) )
            {
                state.f = AP.Math.Sqr(Math.Exp(state.x[0])-2)+AP.Math.Sqr(AP.Math.Sqr(state.x[1]))+AP.Math.Sqr(state.x[2]-state.x[0]);
                state.g[0] = 2*(Math.Exp(state.x[0])-2)*Math.Exp(state.x[0])+2*(state.x[0]-state.x[2]);
                state.g[1] = 4*state.x[1]*AP.Math.Sqr(state.x[1]);
                state.g[2] = 2*(state.x[2]-state.x[0]);
            }
            else
            {
                state.f = Math.Sqrt(AP.Math.MaxRealNumber);
                state.g[0] = Math.Sqrt(AP.Math.MaxRealNumber);
                state.g[1] = 0;
                state.g[2] = 0;
            }
        }


        /*************************************************************************
        Calculate test function #3

        Simple variation of #1, much more nonlinear, with non-zero value at minimum.
        It achieve two goals:
        * makes unlikely premature convergence of algorithm .
        * solves some issues with EpsF stopping condition which arise when
          F(minimum) is zero

        *************************************************************************/
        private static void testfunc3(ref mincg.mincgstate state)
        {
            double s = 0;

            s = 0.001;
            if( (double)(state.x[0])<(double)(100) )
            {
                state.f = AP.Math.Sqr(Math.Exp(state.x[0])-2)+AP.Math.Sqr(AP.Math.Sqr(state.x[1])+s)+AP.Math.Sqr(state.x[2]-state.x[0]);
                state.g[0] = 2*(Math.Exp(state.x[0])-2)*Math.Exp(state.x[0])+2*(state.x[0]-state.x[2]);
                state.g[1] = 2*(AP.Math.Sqr(state.x[1])+s)*2*state.x[1];
                state.g[2] = 2*(state.x[2]-state.x[0]);
            }
            else
            {
                state.f = Math.Sqrt(AP.Math.MaxRealNumber);
                state.g[0] = Math.Sqrt(AP.Math.MaxRealNumber);
                state.g[1] = 0;
                state.g[2] = 0;
            }
        }


        /*************************************************************************
        Silent unit test
        *************************************************************************/
        public static bool testmincgunit_test_silent()
        {
            bool result = new bool();

            result = testmincg(true);
            return result;
        }


        /*************************************************************************
        Unit test
        *************************************************************************/
        public static bool testmincgunit_test()
        {
            bool result = new bool();

            result = testmincg(false);
            return result;
        }
    }
}
