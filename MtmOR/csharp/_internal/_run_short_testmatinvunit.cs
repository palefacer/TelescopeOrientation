using System;
public class _Test
{
	public static int Main(string[] args)
	{
        int seed;
        if( args.Length==1 )
            seed = Convert.ToInt32(args[0]);
        else
            seed = System.DateTime.Now.Millisecond + 1000*System.DateTime.Now.Second + 60*1000*System.DateTime.Now.Minute;
        AP.Math.RndObject = new System.Random(seed);
        try
        {
            if(!alglib.testmatinvunit.testmatinvunit_test_silent())
                throw new Exception("");
        }
        catch(Exception e)
        {
            System.Console.Write("matinv                           FAILED(seed=");
            System.Console.Write("{0,0:d}",seed);
            System.Console.WriteLine(")");
            return 1;
        }
        System.Console.WriteLine("matinv                           OK");
        return 0;
	}
}