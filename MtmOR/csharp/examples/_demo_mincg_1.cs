using System;
using alglib;
public class _Demo
{
    public static int Main(string[] args)
    {
        int n = 0;
        mincg.mincgstate state = new mincg.mincgstate();
        mincg.mincgreport rep = new mincg.mincgreport();
        double[] s = new double[0];
        double x = 0;
        double y = 0;

        
        //
        // Function minimized:
        //     F = (x-1)^4 + (y-x)^2
        // N = 2 - task dimension.
        //
        n = 2;
        s = new double[2];
        s[0] = 10;
        s[1] = 11;
        mincg.mincgcreate(n, ref s, ref state);
        mincg.mincgsetcond(ref state, 0.0, 0.0, 0.00001, 0);
        mincg.mincgsetxrep(ref state, true);
        System.Console.WriteLine();
        System.Console.WriteLine();
        System.Console.Write("F = (x-1)^4 + (y-x)^2");
        System.Console.WriteLine();
        System.Console.Write("OPTIMIZATION STARTED");
        System.Console.WriteLine();
        while( mincg.mincgiteration(ref state) )
        {
            if( state.needfg )
            {
                x = state.x[0];
                y = state.x[1];
                state.f = AP.Math.Sqr(AP.Math.Sqr(x-1))+AP.Math.Sqr(y-x);
                state.g[0] = 4*AP.Math.Sqr(x-1)*(x-1)+2*(x-y);
                state.g[1] = 2*(y-x);
            }
            if( state.xupdated )
            {
                System.Console.Write("    F(");
                System.Console.Write("{0,8:F5}",state.x[0]);
                System.Console.Write(",");
                System.Console.Write("{0,8:F5}",state.x[1]);
                System.Console.Write(")=");
                System.Console.Write("{0,0:F5}",state.f);
                System.Console.WriteLine();
            }
        }
        System.Console.Write("OPTIMIZATION STOPPED");
        System.Console.WriteLine();
        mincg.mincgresults(ref state, ref s, ref rep);
        
        //
        // output results
        //
        System.Console.Write("X = ");
        System.Console.Write("{0,4:F2}",s[0]);
        System.Console.Write(" (should be 1.00)");
        System.Console.WriteLine();
        System.Console.Write("Y = ");
        System.Console.Write("{0,4:F2}",s[1]);
        System.Console.Write(" (should be 1.00)");
        System.Console.WriteLine();
        System.Console.WriteLine();
        System.Console.WriteLine();
        return 0;
    }
}