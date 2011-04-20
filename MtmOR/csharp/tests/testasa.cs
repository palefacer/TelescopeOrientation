
using System;

namespace alglib
{
    public class testasa
    {
        public static bool testminasa(bool silent)
        {
            bool result = new bool();
            bool waserrors = new bool();
            bool referror = new bool();
            bool converror = new bool();
            bool othererrors = new bool();
            int n = 0;
            double[] x = new double[0];
            double[] xe = new double[0];
            double[] c = new double[0];
            double[] bndl = new double[0];
            double[] bndu = new double[0];
            double[] xlast = new double[0];
            double fprev = 0;
            double xprev = 0;
            double stpmax = 0;
            int i = 0;
            int j = 0;
            double v = 0;
            double s = 0;
            double tol = 0;
            int algotype = 0;
            double[,] a = new double[0,0];
            minasa.minasastate state = new minasa.minasastate();
            minasa.minasareport rep = new minasa.minasareport();
            int i_ = 0;

            waserrors = false;
            referror = false;
            converror = false;
            othererrors = false;
            
            //
            // Different algorithms
            //
            for(algotype=-1; algotype<=1; algotype++)
            {
                
                //
                // reference problem, simple convex optimization
                //
                for(n=1; n<=5; n++)
                {
                    
                    //
                    // min(x'*diag(c)*x) on a random box
                    //
                    x = new double[n];
                    xe = new double[n];
                    c = new double[n];
                    bndl = new double[n];
                    bndu = new double[n];
                    for(i=0; i<=n-1; i++)
                    {
                        c[i] = 1+AP.Math.RandomReal();
                        xe[i] = 4*AP.Math.RandomReal()-2;
                        bndl[i] = -Math.Max(AP.Math.RandomReal(), 0.2);
                        bndu[i] = +Math.Max(AP.Math.RandomReal(), 0.2);
                        x[i] = 0.5*(bndl[i]+bndu[i]);
                    }
                    tol = 0.001;
                    minasa.minasacreate(n, ref x, ref bndl, ref bndu, ref state);
                    minasa.minasasetcond(ref state, tol, 0.0, 0.0, 0);
                    minasa.minasasetalgorithm(ref state, algotype);
                    while( minasa.minasaiteration(ref state) )
                    {
                        checkbounds(ref state.x, ref bndl, ref bndu, n, ref othererrors);
                        state.f = 0;
                        for(i=0; i<=n-1; i++)
                        {
                            state.f = state.f+c[i]*AP.Math.Sqr(state.x[i]-xe[i]);
                            state.g[i] = 2*c[i]*(state.x[i]-xe[i]);
                        }
                    }
                    minasa.minasaresults(ref state, ref x, ref rep);
                    referror = referror | rep.terminationtype<=0;
                    for(i=0; i<=n-1; i++)
                    {
                        referror = referror | (double)(Math.Abs(asaboundval(xe[i], bndl[i], bndu[i])-x[i]))>(double)(0.01);
                    }
                }
                
                //
                // reference problem 2: non-convex optimization on [-2,2] x [1,2]
                //
                // A saddle function is minimized:
                // * stationary point [0,0] (non-feasible)
                // * constrained minimum [-2,2].
                // * starting point [+2,2]
                //
                // Path from start to end may be very complex, with multiple changes
                // in active constraints, so it is interesting task for our method.
                //
                // Scale parameter is used to make optimization more interesting
                // during GPA runs.
                //
                x = new double[2];
                bndl = new double[2];
                bndu = new double[2];
                bndl[0] = -2;
                bndu[0] = 2;
                x[0] = 2;
                bndl[1] = 1;
                bndu[1] = 2;
                x[1] = 2;
                tol = 0.001;
                s = 0.01;
                minasa.minasacreate(2, ref x, ref bndl, ref bndu, ref state);
                minasa.minasasetcond(ref state, tol, 0.0, 0.0, 0);
                minasa.minasasetalgorithm(ref state, algotype);
                while( minasa.minasaiteration(ref state) )
                {
                    checkbounds(ref state.x, ref bndl, ref bndu, 2, ref othererrors);
                    state.f = s*(AP.Math.Sqr(state.x[0]+state.x[1])-AP.Math.Sqr(state.x[0]-state.x[1]));
                    state.g[0] = s*(2*(state.x[0]+state.x[1])-2*(state.x[0]-state.x[1]));
                    state.g[1] = s*(2*(state.x[0]+state.x[1])+2*(state.x[0]-state.x[1]));
                }
                minasa.minasaresults(ref state, ref x, ref rep);
                referror = referror | rep.terminationtype<=0 | (double)(Math.Abs(state.x[0]+2))>(double)(0.01) | (double)(Math.Abs(state.x[1]-2))>(double)(0.01);
                
                //
                // function #1 with 'x[0]>=ln(2)' constraint.
                // may show very interesting behavior.
                //
                x = new double[3];
                bndl = new double[3];
                bndu = new double[3];
                n = 3;
                for(i=0; i<=2; i++)
                {
                    bndl[i] = -10000;
                    bndu[i] = +10000;
                }
                bndl[0] = Math.Log(2);
                for(i=0; i<=2; i++)
                {
                    x[i] = 3*AP.Math.RandomReal()+3;
                }
                minasa.minasacreate(n, ref x, ref bndl, ref bndu, ref state);
                minasa.minasasetcond(ref state, 0.0000001, 0.0, 0.0, 0);
                minasa.minasasetalgorithm(ref state, algotype);
                while( minasa.minasaiteration(ref state) )
                {
                    checkbounds(ref state.x, ref bndl, ref bndu, n, ref othererrors);
                    testfunc1(ref state);
                }
                minasa.minasaresults(ref state, ref x, ref rep);
                referror = referror | rep.terminationtype<=0;
                referror = referror | (double)(Math.Abs(x[0]-Math.Log(2)))>(double)(0.05);
                referror = referror | (double)(Math.Abs(x[1]))>(double)(0.05);
                referror = referror | (double)(Math.Abs(x[2]-Math.Log(2)))>(double)(0.05);
                
                //
                // Testing convergence properties
                //
                x = new double[3];
                bndl = new double[3];
                bndu = new double[3];
                n = 3;
                for(i=0; i<=2; i++)
                {
                    bndl[i] = -10000;
                    bndu[i] = +10000;
                }
                bndl[0] = Math.Log(2);
                for(i=0; i<=2; i++)
                {
                    x[i] = 3*AP.Math.RandomReal()+3;
                }
                minasa.minasacreate(n, ref x, ref bndl, ref bndu, ref state);
                minasa.minasasetcond(ref state, 0.001, 0.0, 0.0, 0);
                minasa.minasasetalgorithm(ref state, algotype);
                while( minasa.minasaiteration(ref state) )
                {
                    checkbounds(ref state.x, ref bndl, ref bndu, n, ref othererrors);
                    testfunc3(ref state);
                }
                minasa.minasaresults(ref state, ref x, ref rep);
                converror = converror | rep.terminationtype!=4;
                for(i=0; i<=2; i++)
                {
                    x[i] = 3*AP.Math.RandomReal()+3;
                }
                minasa.minasacreate(n, ref x, ref bndl, ref bndu, ref state);
                minasa.minasasetcond(ref state, 0.0, 0.001, 0.0, 0);
                minasa.minasasetalgorithm(ref state, algotype);
                while( minasa.minasaiteration(ref state) )
                {
                    checkbounds(ref state.x, ref bndl, ref bndu, n, ref othererrors);
                    testfunc3(ref state);
                }
                minasa.minasaresults(ref state, ref x, ref rep);
                converror = converror | rep.terminationtype!=1;
                for(i=0; i<=2; i++)
                {
                    x[i] = 3*AP.Math.RandomReal()+3;
                }
                minasa.minasacreate(n, ref x, ref bndl, ref bndu, ref state);
                minasa.minasasetcond(ref state, 0.0, 0.0, 0.001, 0);
                minasa.minasasetalgorithm(ref state, algotype);
                while( minasa.minasaiteration(ref state) )
                {
                    checkbounds(ref state.x, ref bndl, ref bndu, n, ref othererrors);
                    testfunc3(ref state);
                }
                minasa.minasaresults(ref state, ref x, ref rep);
                converror = converror | rep.terminationtype!=2;
                for(i=0; i<=2; i++)
                {
                    x[i] = 3*AP.Math.RandomReal()+3;
                }
                minasa.minasacreate(n, ref x, ref bndl, ref bndu, ref state);
                minasa.minasasetcond(ref state, 0.0, 0.0, 0.0, 3);
                minasa.minasasetalgorithm(ref state, algotype);
                while( minasa.minasaiteration(ref state) )
                {
                    checkbounds(ref state.x, ref bndl, ref bndu, n, ref othererrors);
                    testfunc3(ref state);
                }
                minasa.minasaresults(ref state, ref x, ref rep);
                converror = converror | !(rep.terminationtype==5 & rep.iterationscount==3 | rep.terminationtype==7);
                
                //
                // Other properties
                //
                //
                // Other properties:
                // 1. test reports (F should form monotone sequence)
                // 2. test maximum step
                //
                n = 50;
                x = new double[n];
                xlast = new double[n];
                bndl = new double[n];
                bndu = new double[n];
                for(i=0; i<=n-1; i++)
                {
                    x[i] = 1;
                    xlast[i] = AP.Math.RandomReal();
                    bndl[i] = -100000;
                    bndu[i] = +100000;
                }
                minasa.minasacreate(n, ref x, ref bndl, ref bndu, ref state);
                minasa.minasasetcond(ref state, 0, 0, 0, 100);
                minasa.minasasetxrep(ref state, true);
                fprev = AP.Math.MaxRealNumber;
                while( minasa.minasaiteration(ref state) )
                {
                    checkbounds(ref state.x, ref bndl, ref bndu, n, ref othererrors);
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
                minasa.minasaresults(ref state, ref x, ref rep);
                for(i=0; i<=n-1; i++)
                {
                    othererrors = othererrors | (double)(x[i])!=(double)(xlast[i]);
                }
                n = 1;
                x = new double[n];
                bndl = new double[n];
                bndu = new double[n];
                x[0] = 100;
                bndl[0] = -1000000;
                bndu[0] = +1000000;
                stpmax = 0.05+0.05*AP.Math.RandomReal();
                minasa.minasacreate(n, ref x, ref bndl, ref bndu, ref state);
                minasa.minasasetcond(ref state, 1.0E-9, 0, 0, 0);
                minasa.minasasetstpmax(ref state, stpmax);
                minasa.minasasetxrep(ref state, true);
                xprev = x[0];
                while( minasa.minasaiteration(ref state) )
                {
                    checkbounds(ref state.x, ref bndl, ref bndu, n, ref othererrors);
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
            waserrors = referror | converror | othererrors;
            if( !silent )
            {
                System.Console.Write("TESTING ASA OPTIMIZATION");
                System.Console.WriteLine();
                System.Console.Write("REFERENCE PROBLEMS:                       ");
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

        It may show very interesting behavior when optimized with 'x[0]>=ln(2)'
        constraint.
        *************************************************************************/
        private static void testfunc1(ref minasa.minasastate state)
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
        private static void testfunc2(ref minasa.minasastate state)
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
        private static void testfunc3(ref minasa.minasastate state)
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
        Checks that X is bounded with respect to BndL/BndU.

        If it is not, True is assigned to the Err variable (which is not changed
        otherwise).
        *************************************************************************/
        private static void checkbounds(ref double[] x,
            ref double[] bndl,
            ref double[] bndu,
            int n,
            ref bool err)
        {
            int i = 0;

            for(i=0; i<=n-1; i++)
            {
                if( (double)(x[i])<(double)(bndl[i]) | (double)(x[i])>(double)(bndu[i]) )
                {
                    err = true;
                }
            }
        }


        /*************************************************************************
        'bound' value: map X to [B1,B2]
        *************************************************************************/
        private static double asaboundval(double x,
            double b1,
            double b2)
        {
            double result = 0;

            if( (double)(x)<=(double)(b1) )
            {
                result = b1;
                return result;
            }
            if( (double)(x)>=(double)(b2) )
            {
                result = b2;
                return result;
            }
            result = x;
            return result;
        }


        /*************************************************************************
        Silent unit test
        *************************************************************************/
        public static bool testasa_test_silent()
        {
            bool result = new bool();

            result = testminasa(true);
            return result;
        }


        /*************************************************************************
        Unit test
        *************************************************************************/
        public static bool testasa_test()
        {
            bool result = new bool();

            result = testminasa(false);
            return result;
        }
    }
}
