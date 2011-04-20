using System;
using alglib;
public class _Demo
{
    public static int Main(string[] args)
    {
        int n = 0;
        int i = 0;
        minasa.minasastate state = new minasa.minasastate();
        minasa.minasareport rep = new minasa.minasareport();
        double[] s = new double[0];
        double[] bndl = new double[0];
        double[] bndu = new double[0];
        double x = 0;
        double y = 0;
        double z = 0;

        
        //
        // Function being minimized:
        //     F = x+4y+9z subject to 0<=x<=1, 0<=y<=1, 0<=z<=1.
        //
        // Take a look at MinASASetStpMax() - it restricts step length by
        // a small value, so we can see the current point traveling through
        // a feasible set, sticking to its bounds.
        //
        n = 3;
        s = new double[n];
        bndl = new double[n];
        bndu = new double[n];
        for(i=0; i<=n-1; i++)
        {
            s[i] = 1;
            bndl[i] = 0;
            bndu[i] = 1;
        }
        minasa.minasacreate(n, ref s, ref bndl, ref bndu, ref state);
        minasa.minasasetcond(ref state, 0.0, 0.0, 0.00001, 0);
        minasa.minasasetxrep(ref state, true);
        minasa.minasasetstpmax(ref state, 0.2);
        System.Console.WriteLine();
        System.Console.WriteLine();
        System.Console.Write("F = x+4y+9z subject to 0<=x<=1, 0<=y<=1, 0<=z<=1");
        System.Console.WriteLine();
        System.Console.Write("OPTIMIZATION STARTED");
        System.Console.WriteLine();
        while( minasa.minasaiteration(ref state) )
        {
            if( state.needfg )
            {
                x = state.x[0];
                y = state.x[1];
                z = state.x[2];
                state.f = x+4*y+9*z;
                state.g[0] = 1;
                state.g[1] = 4;
                state.g[2] = 9;
            }
            if( state.xupdated )
            {
                System.Console.Write("    F(");
                System.Console.Write("{0,4:F2}",state.x[0]);
                System.Console.Write(", ");
                System.Console.Write("{0,4:F2}",state.x[1]);
                System.Console.Write(", ");
                System.Console.Write("{0,4:F2}",state.x[2]);
                System.Console.Write(") = ");
                System.Console.Write("{0,0:F3}",state.f);
                System.Console.WriteLine();
            }
        }
        System.Console.Write("OPTIMIZATION STOPPED");
        System.Console.WriteLine();
        minasa.minasaresults(ref state, ref s, ref rep);
        
        //
        // output results
        //
        System.Console.Write("X = ");
        System.Console.Write("{0,4:F2}",s[0]);
        System.Console.Write(" (should be 0.00)");
        System.Console.WriteLine();
        System.Console.Write("Y = ");
        System.Console.Write("{0,4:F2}",s[1]);
        System.Console.Write(" (should be 0.00)");
        System.Console.WriteLine();
        System.Console.Write("Z = ");
        System.Console.Write("{0,4:F2}",s[2]);
        System.Console.Write(" (should be 0.00)");
        System.Console.WriteLine();
        System.Console.WriteLine();
        System.Console.WriteLine();
        return 0;
    }
}