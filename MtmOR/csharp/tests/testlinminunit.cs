
using System;

namespace alglib
{
    public class testlinminunit
    {
        public static bool testlinmin(bool silent)
        {
            bool result = new bool();
            bool waserrors = new bool();

            waserrors = false;
            if( !silent )
            {
                System.Console.Write("TESTING LINMIN");
                System.Console.WriteLine();
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
        Silent unit test
        *************************************************************************/
        public static bool testlinminunit_test_silent()
        {
            bool result = new bool();

            result = testlinmin(true);
            return result;
        }


        /*************************************************************************
        Unit test
        *************************************************************************/
        public static bool testlinminunit_test()
        {
            bool result = new bool();

            result = testlinmin(false);
            return result;
        }
    }
}
