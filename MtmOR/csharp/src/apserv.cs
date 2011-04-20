/*************************************************************************
Copyright (c) 2009, Sergey Bochkanov (ALGLIB project).

>>> SOURCE LICENSE >>>
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation (www.fsf.org); either version 2 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License is available at
http://www.fsf.org/licensing/licenses

>>> END OF LICENSE >>>
*************************************************************************/

using System;

namespace alglib
{
    public class apserv
    {
        /*************************************************************************
        This  function  generates  1-dimensional  general  interpolation task with
        moderate Lipshitz constant (close to 1.0)

        If N=1 then suborutine generates only one point at the middle of [A,B]

          -- ALGLIB --
             Copyright 02.12.2009 by Bochkanov Sergey
        *************************************************************************/
        public static void taskgenint1d(double a,
            double b,
            int n,
            ref double[] x,
            ref double[] y)
        {
            int i = 0;
            double h = 0;

            System.Diagnostics.Debug.Assert(n>=1, "TaskGenInterpolationEqdist1D: N<1!");
            x = new double[n];
            y = new double[n];
            if( n>1 )
            {
                x[0] = a;
                y[0] = 2*AP.Math.RandomReal()-1;
                h = (b-a)/(n-1);
                for(i=1; i<=n-1; i++)
                {
                    if( i!=n-1 )
                    {
                        x[i] = a+(i+0.2*(2*AP.Math.RandomReal()-1))*h;
                    }
                    else
                    {
                        x[i] = b;
                    }
                    y[i] = y[i-1]+(2*AP.Math.RandomReal()-1)*(x[i]-x[i-1]);
                }
            }
            else
            {
                x[0] = 0.5*(a+b);
                y[0] = 2*AP.Math.RandomReal()-1;
            }
        }


        /*************************************************************************
        This function generates  1-dimensional equidistant interpolation task with
        moderate Lipshitz constant (close to 1.0)

        If N=1 then suborutine generates only one point at the middle of [A,B]

          -- ALGLIB --
             Copyright 02.12.2009 by Bochkanov Sergey
        *************************************************************************/
        public static void taskgenint1dequidist(double a,
            double b,
            int n,
            ref double[] x,
            ref double[] y)
        {
            int i = 0;
            double h = 0;

            System.Diagnostics.Debug.Assert(n>=1, "TaskGenInterpolationEqdist1D: N<1!");
            x = new double[n];
            y = new double[n];
            if( n>1 )
            {
                x[0] = a;
                y[0] = 2*AP.Math.RandomReal()-1;
                h = (b-a)/(n-1);
                for(i=1; i<=n-1; i++)
                {
                    x[i] = a+i*h;
                    y[i] = y[i-1]+(2*AP.Math.RandomReal()-1)*h;
                }
            }
            else
            {
                x[0] = 0.5*(a+b);
                y[0] = 2*AP.Math.RandomReal()-1;
            }
        }


        /*************************************************************************
        This function generates  1-dimensional Chebyshev-1 interpolation task with
        moderate Lipshitz constant (close to 1.0)

        If N=1 then suborutine generates only one point at the middle of [A,B]

          -- ALGLIB --
             Copyright 02.12.2009 by Bochkanov Sergey
        *************************************************************************/
        public static void taskgenint1dcheb1(double a,
            double b,
            int n,
            ref double[] x,
            ref double[] y)
        {
            int i = 0;

            System.Diagnostics.Debug.Assert(n>=1, "TaskGenInterpolation1DCheb1: N<1!");
            x = new double[n];
            y = new double[n];
            if( n>1 )
            {
                for(i=0; i<=n-1; i++)
                {
                    x[i] = 0.5*(b+a)+0.5*(b-a)*Math.Cos(Math.PI*(2*i+1)/(2*n));
                    if( i==0 )
                    {
                        y[i] = 2*AP.Math.RandomReal()-1;
                    }
                    else
                    {
                        y[i] = y[i-1]+(2*AP.Math.RandomReal()-1)*(x[i]-x[i-1]);
                    }
                }
            }
            else
            {
                x[0] = 0.5*(a+b);
                y[0] = 2*AP.Math.RandomReal()-1;
            }
        }


        /*************************************************************************
        This function generates  1-dimensional Chebyshev-2 interpolation task with
        moderate Lipshitz constant (close to 1.0)

        If N=1 then suborutine generates only one point at the middle of [A,B]

          -- ALGLIB --
             Copyright 02.12.2009 by Bochkanov Sergey
        *************************************************************************/
        public static void taskgenint1dcheb2(double a,
            double b,
            int n,
            ref double[] x,
            ref double[] y)
        {
            int i = 0;

            System.Diagnostics.Debug.Assert(n>=1, "TaskGenInterpolation1DCheb2: N<1!");
            x = new double[n];
            y = new double[n];
            if( n>1 )
            {
                for(i=0; i<=n-1; i++)
                {
                    x[i] = 0.5*(b+a)+0.5*(b-a)*Math.Cos(Math.PI*i/(n-1));
                    if( i==0 )
                    {
                        y[i] = 2*AP.Math.RandomReal()-1;
                    }
                    else
                    {
                        y[i] = y[i-1]+(2*AP.Math.RandomReal()-1)*(x[i]-x[i-1]);
                    }
                }
            }
            else
            {
                x[0] = 0.5*(a+b);
                y[0] = 2*AP.Math.RandomReal()-1;
            }
        }


        /*************************************************************************
        This function checks that all values from X[] are distinct. It does more
        than just usual floating point comparison:
        * first, it calculates max(X) and min(X)
        * second, it maps X[] from [min,max] to [1,2]
        * only at this stage actual comparison is done

        The meaning of such check is to ensure that all values are "distinct enough"
        and will not cause interpolation subroutine to fail.

        NOTE:
            X[] must be sorted by ascending (subroutine ASSERT's it)

          -- ALGLIB --
             Copyright 02.12.2009 by Bochkanov Sergey
        *************************************************************************/
        public static bool apservaredistinct(double[] x,
            int n)
        {
            bool result = new bool();
            bool issorted = new bool();
            double a = 0;
            double b = 0;
            int i = 0;

            x = (double[])x.Clone();

            System.Diagnostics.Debug.Assert(n>=1, "APSERVIsDistinct: internal error!");
            if( n==1 )
            {
                
                //
                // everything is alright, it is up to caller to decide whether it
                // can interpolate something with just one point
                //
                result = true;
                return result;
            }
            a = x[0];
            b = x[0];
            for(i=1; i<=n-1; i++)
            {
                a = Math.Min(a, x[i]);
                b = Math.Max(b, x[i]);
                System.Diagnostics.Debug.Assert((double)(x[i])>=(double)(x[i-1]), "APSERVIsDistinct: internal error!");
            }
            for(i=0; i<=n-1; i++)
            {
                x[i] = (x[i]-a)/(b-a)+1;
            }
            for(i=1; i<=n-1; i++)
            {
                if( (double)(x[i])==(double)(x[i-1]) )
                {
                    result = false;
                    return result;
                }
            }
            result = true;
            return result;
        }


        /*************************************************************************
        Safe sqrt(x^2+y^2)

          -- ALGLIB --
             Copyright by Bochkanov Sergey
        *************************************************************************/
        public static double safepythag2(double x,
            double y)
        {
            double result = 0;
            double w = 0;
            double xabs = 0;
            double yabs = 0;
            double z = 0;

            xabs = Math.Abs(x);
            yabs = Math.Abs(y);
            w = Math.Max(xabs, yabs);
            z = Math.Min(xabs, yabs);
            if( (double)(z)==(double)(0) )
            {
                result = w;
            }
            else
            {
                result = w*Math.Sqrt(1+AP.Math.Sqr(z/w));
            }
            return result;
        }


        /*************************************************************************
        Safe sqrt(x^2+y^2)

          -- ALGLIB --
             Copyright by Bochkanov Sergey
        *************************************************************************/
        public static double safepythag3(double x,
            double y,
            double z)
        {
            double result = 0;
            double w = 0;

            w = Math.Max(Math.Abs(x), Math.Max(Math.Abs(y), Math.Abs(z)));
            if( (double)(w)==(double)(0) )
            {
                result = 0;
                return result;
            }
            x = x/w;
            y = y/w;
            z = z/w;
            result = w*Math.Sqrt(AP.Math.Sqr(x)+AP.Math.Sqr(y)+AP.Math.Sqr(z));
            return result;
        }


        /*************************************************************************
        This function makes periodic mapping of X to [A,B].

        It accepts X, A, B (A>B). It returns T which lies in  [A,B] and integer K,
        such that X = T + K*(B-A).

        NOTES:
        * K is represented as real value, although actually it is integer
        * T is guaranteed to be in [A,B]
        * T replaces X

          -- ALGLIB --
             Copyright by Bochkanov Sergey
        *************************************************************************/
        public static void apperiodicmap(ref double x,
            double a,
            double b,
            ref double k)
        {
            System.Diagnostics.Debug.Assert((double)(a)<(double)(b), "APPeriodicMap: internal error!");
            k = (int)Math.Floor((x-a)/(b-a));
            x = x-k*(b-a);
            while( (double)(x)<(double)(a) )
            {
                x = x+(b-a);
                k = k-1;
            }
            while( (double)(x)>(double)(b) )
            {
                x = x-(b-a);
                k = k+1;
            }
            x = Math.Max(x, a);
            x = Math.Min(x, b);
        }
    }
}
