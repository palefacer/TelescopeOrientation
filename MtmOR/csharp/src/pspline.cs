/*************************************************************************
Copyright (c) 2006-2010, Sergey Bochkanov (ALGLIB project).

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
    public class pspline
    {
        /*************************************************************************
        Parametric spline inteprolant: 2-dimensional curve.

        You should not try to access its members directly - use PSpline2XXXXXXXX()
        functions instead.
        *************************************************************************/
        public struct pspline2interpolant
        {
            public int n;
            public bool periodic;
            public double[] p;
            public spline1d.spline1dinterpolant x;
            public spline1d.spline1dinterpolant y;
        };


        /*************************************************************************
        Parametric spline inteprolant: 3-dimensional curve.

        You should not try to access its members directly - use PSpline3XXXXXXXX()
        functions instead.
        *************************************************************************/
        public struct pspline3interpolant
        {
            public int n;
            public bool periodic;
            public double[] p;
            public spline1d.spline1dinterpolant x;
            public spline1d.spline1dinterpolant y;
            public spline1d.spline1dinterpolant z;
        };




        /*************************************************************************
        This function  builds  non-periodic 2-dimensional parametric spline  which
        starts at (X[0],Y[0]) and ends at (X[N-1],Y[N-1]).

        INPUT PARAMETERS:
            XY  -   points, array[0..N-1,0..1].
                    XY[I,0:1] corresponds to the Ith point.
                    Order of points is important!
            N   -   points count, N>=5 for Akima splines, N>=2 for other types  of
                    splines.
            ST  -   spline type:
                    * 0     Akima spline
                    * 1     parabolically terminated Catmull-Rom spline (Tension=0)
                    * 2     parabolically terminated cubic spline
            PT  -   parameterization type:
                    * 0     uniform
                    * 1     chord length
                    * 2     centripetal

        OUTPUT PARAMETERS:
            P   -   parametric spline interpolant


        NOTES:
        * this function  assumes  that  there all consequent points  are distinct.
          I.e. (x0,y0)<>(x1,y1),  (x1,y1)<>(x2,y2),  (x2,y2)<>(x3,y3)  and  so on.
          However, non-consequent points may coincide, i.e. we can  have  (x0,y0)=
          =(x2,y2).

          -- ALGLIB PROJECT --
             Copyright 28.05.2010 by Bochkanov Sergey
        *************************************************************************/
        public static void pspline2build(double[,] xy,
            int n,
            int st,
            int pt,
            ref pspline2interpolant p)
        {
            double[] tmp = new double[0];
            double v = 0;
            int i = 0;
            int i_ = 0;

            xy = (double[,])xy.Clone();

            System.Diagnostics.Debug.Assert(st>=0 & st<=2, "PSpline2Build: incorrect spline type!");
            System.Diagnostics.Debug.Assert(pt>=0 & pt<=2, "PSpline2Build: incorrect parameterization type!");
            if( st==0 )
            {
                System.Diagnostics.Debug.Assert(n>=5, "PSpline2Build: N<5 (minimum value for Akima splines)!");
            }
            else
            {
                System.Diagnostics.Debug.Assert(n>=2, "PSpline2Build: N<2!");
            }
            
            //
            // Prepare
            //
            p.n = n;
            p.periodic = false;
            tmp = new double[n];
            
            //
            // Build parameterization, check that all parameters are distinct
            //
            pspline2par(ref xy, n, pt, ref p.p);
            System.Diagnostics.Debug.Assert(apserv.apservaredistinct(p.p, n), "PSpline2Build: consequent points are too close!");
            
            //
            // Build splines
            //
            if( st==0 )
            {
                for(i_=0; i_<=n-1;i_++)
                {
                    tmp[i_] = xy[i_,0];
                }
                spline1d.spline1dbuildakima(p.p, tmp, n, ref p.x);
                for(i_=0; i_<=n-1;i_++)
                {
                    tmp[i_] = xy[i_,1];
                }
                spline1d.spline1dbuildakima(p.p, tmp, n, ref p.y);
            }
            if( st==1 )
            {
                for(i_=0; i_<=n-1;i_++)
                {
                    tmp[i_] = xy[i_,0];
                }
                spline1d.spline1dbuildcatmullrom(p.p, tmp, n, 0, 0.0, ref p.x);
                for(i_=0; i_<=n-1;i_++)
                {
                    tmp[i_] = xy[i_,1];
                }
                spline1d.spline1dbuildcatmullrom(p.p, tmp, n, 0, 0.0, ref p.y);
            }
            if( st==2 )
            {
                for(i_=0; i_<=n-1;i_++)
                {
                    tmp[i_] = xy[i_,0];
                }
                spline1d.spline1dbuildcubic(p.p, tmp, n, 0, 0.0, 0, 0.0, ref p.x);
                for(i_=0; i_<=n-1;i_++)
                {
                    tmp[i_] = xy[i_,1];
                }
                spline1d.spline1dbuildcubic(p.p, tmp, n, 0, 0.0, 0, 0.0, ref p.y);
            }
        }


        /*************************************************************************
        This function  builds  non-periodic 3-dimensional parametric spline  which
        starts at (X[0],Y[0],Z[0]) and ends at (X[N-1],Y[N-1],Z[N-1]).

        Same as PSpline2Build() function, but for 3D, so we  won't  duplicate  its
        description here.

          -- ALGLIB PROJECT --
             Copyright 28.05.2010 by Bochkanov Sergey
        *************************************************************************/
        public static void pspline3build(double[,] xy,
            int n,
            int st,
            int pt,
            ref pspline3interpolant p)
        {
            double[] tmp = new double[0];
            double v = 0;
            int i = 0;
            int i_ = 0;

            xy = (double[,])xy.Clone();

            System.Diagnostics.Debug.Assert(st>=0 & st<=2, "PSpline3Build: incorrect spline type!");
            System.Diagnostics.Debug.Assert(pt>=0 & pt<=2, "PSpline3Build: incorrect parameterization type!");
            if( st==0 )
            {
                System.Diagnostics.Debug.Assert(n>=5, "PSpline3Build: N<5 (minimum value for Akima splines)!");
            }
            else
            {
                System.Diagnostics.Debug.Assert(n>=2, "PSpline3Build: N<2!");
            }
            
            //
            // Prepare
            //
            p.n = n;
            p.periodic = false;
            tmp = new double[n];
            
            //
            // Build parameterization, check that all parameters are distinct
            //
            pspline3par(ref xy, n, pt, ref p.p);
            System.Diagnostics.Debug.Assert(apserv.apservaredistinct(p.p, n), "PSpline3Build: consequent points are too close!");
            
            //
            // Build splines
            //
            if( st==0 )
            {
                for(i_=0; i_<=n-1;i_++)
                {
                    tmp[i_] = xy[i_,0];
                }
                spline1d.spline1dbuildakima(p.p, tmp, n, ref p.x);
                for(i_=0; i_<=n-1;i_++)
                {
                    tmp[i_] = xy[i_,1];
                }
                spline1d.spline1dbuildakima(p.p, tmp, n, ref p.y);
                for(i_=0; i_<=n-1;i_++)
                {
                    tmp[i_] = xy[i_,2];
                }
                spline1d.spline1dbuildakima(p.p, tmp, n, ref p.z);
            }
            if( st==1 )
            {
                for(i_=0; i_<=n-1;i_++)
                {
                    tmp[i_] = xy[i_,0];
                }
                spline1d.spline1dbuildcatmullrom(p.p, tmp, n, 0, 0.0, ref p.x);
                for(i_=0; i_<=n-1;i_++)
                {
                    tmp[i_] = xy[i_,1];
                }
                spline1d.spline1dbuildcatmullrom(p.p, tmp, n, 0, 0.0, ref p.y);
                for(i_=0; i_<=n-1;i_++)
                {
                    tmp[i_] = xy[i_,2];
                }
                spline1d.spline1dbuildcatmullrom(p.p, tmp, n, 0, 0.0, ref p.z);
            }
            if( st==2 )
            {
                for(i_=0; i_<=n-1;i_++)
                {
                    tmp[i_] = xy[i_,0];
                }
                spline1d.spline1dbuildcubic(p.p, tmp, n, 0, 0.0, 0, 0.0, ref p.x);
                for(i_=0; i_<=n-1;i_++)
                {
                    tmp[i_] = xy[i_,1];
                }
                spline1d.spline1dbuildcubic(p.p, tmp, n, 0, 0.0, 0, 0.0, ref p.y);
                for(i_=0; i_<=n-1;i_++)
                {
                    tmp[i_] = xy[i_,2];
                }
                spline1d.spline1dbuildcubic(p.p, tmp, n, 0, 0.0, 0, 0.0, ref p.z);
            }
        }


        /*************************************************************************
        This  function  builds  periodic  2-dimensional  parametric  spline  which
        starts at (X[0],Y[0]), goes through all points to (X[N-1],Y[N-1]) and then
        back to (X[0],Y[0]).

        INPUT PARAMETERS:
            XY  -   points, array[0..N-1,0..1].
                    XY[I,0:1] corresponds to the Ith point.
                    XY[N-1,0:1] must be different from XY[0,0:1].
                    Order of points is important!
            N   -   points count, N>=3 for other types of splines.
            ST  -   spline type:
                    * 1     Catmull-Rom spline (Tension=0) with cyclic boundary conditions
                    * 2     cubic spline with cyclic boundary conditions
            PT  -   parameterization type:
                    * 0     uniform
                    * 1     chord length
                    * 2     centripetal

        OUTPUT PARAMETERS:
            P   -   parametric spline interpolant


        NOTES:
        * this function  assumes  that there all consequent points  are  distinct.
          I.e. (x0,y0)<>(x1,y1), (x1,y1)<>(x2,y2),  (x2,y2)<>(x3,y3)  and  so  on.
          However, non-consequent points may coincide, i.e. we can  have  (x0,y0)=
          =(x2,y2).
        * last point of sequence is NOT equal to the first  point.  You  shouldn't
          make curve "explicitly periodic" by making them equal.

          -- ALGLIB PROJECT --
             Copyright 28.05.2010 by Bochkanov Sergey
        *************************************************************************/
        public static void pspline2buildperiodic(double[,] xy,
            int n,
            int st,
            int pt,
            ref pspline2interpolant p)
        {
            double[,] xyp = new double[0,0];
            double[] tmp = new double[0];
            double v = 0;
            int i = 0;
            int i_ = 0;

            xy = (double[,])xy.Clone();

            System.Diagnostics.Debug.Assert(st>=1 & st<=2, "PSpline2BuildPeriodic: incorrect spline type!");
            System.Diagnostics.Debug.Assert(pt>=0 & pt<=2, "PSpline2BuildPeriodic: incorrect parameterization type!");
            System.Diagnostics.Debug.Assert(n>=3, "PSpline2BuildPeriodic: N<3!");
            
            //
            // Prepare
            //
            p.n = n;
            p.periodic = true;
            tmp = new double[n+1];
            xyp = new double[n+1, 2];
            for(i_=0; i_<=n-1;i_++)
            {
                xyp[i_,0] = xy[i_,0];
            }
            for(i_=0; i_<=n-1;i_++)
            {
                xyp[i_,1] = xy[i_,1];
            }
            for(i_=0; i_<=1;i_++)
            {
                xyp[n,i_] = xy[0,i_];
            }
            
            //
            // Build parameterization, check that all parameters are distinct
            //
            pspline2par(ref xyp, n+1, pt, ref p.p);
            System.Diagnostics.Debug.Assert(apserv.apservaredistinct(p.p, n+1), "PSpline2BuildPeriodic: consequent (or first and last) points are too close!");
            
            //
            // Build splines
            //
            if( st==1 )
            {
                for(i_=0; i_<=n;i_++)
                {
                    tmp[i_] = xyp[i_,0];
                }
                spline1d.spline1dbuildcatmullrom(p.p, tmp, n+1, -1, 0.0, ref p.x);
                for(i_=0; i_<=n;i_++)
                {
                    tmp[i_] = xyp[i_,1];
                }
                spline1d.spline1dbuildcatmullrom(p.p, tmp, n+1, -1, 0.0, ref p.y);
            }
            if( st==2 )
            {
                for(i_=0; i_<=n;i_++)
                {
                    tmp[i_] = xyp[i_,0];
                }
                spline1d.spline1dbuildcubic(p.p, tmp, n+1, -1, 0.0, -1, 0.0, ref p.x);
                for(i_=0; i_<=n;i_++)
                {
                    tmp[i_] = xyp[i_,1];
                }
                spline1d.spline1dbuildcubic(p.p, tmp, n+1, -1, 0.0, -1, 0.0, ref p.y);
            }
        }


        /*************************************************************************
        This  function  builds  periodic  3-dimensional  parametric  spline  which
        starts at (X[0],Y[0],Z[0]), goes through all points to (X[N-1],Y[N-1],Z[N-1])
        and then back to (X[0],Y[0],Z[0]).

        Same as PSpline2Build() function, but for 3D, so we  won't  duplicate  its
        description here.

          -- ALGLIB PROJECT --
             Copyright 28.05.2010 by Bochkanov Sergey
        *************************************************************************/
        public static void pspline3buildperiodic(double[,] xy,
            int n,
            int st,
            int pt,
            ref pspline3interpolant p)
        {
            double[,] xyp = new double[0,0];
            double[] tmp = new double[0];
            double v = 0;
            int i = 0;
            int i_ = 0;

            xy = (double[,])xy.Clone();

            System.Diagnostics.Debug.Assert(st>=1 & st<=2, "PSpline3BuildPeriodic: incorrect spline type!");
            System.Diagnostics.Debug.Assert(pt>=0 & pt<=2, "PSpline3BuildPeriodic: incorrect parameterization type!");
            System.Diagnostics.Debug.Assert(n>=3, "PSpline3BuildPeriodic: N<3!");
            
            //
            // Prepare
            //
            p.n = n;
            p.periodic = true;
            tmp = new double[n+1];
            xyp = new double[n+1, 3];
            for(i_=0; i_<=n-1;i_++)
            {
                xyp[i_,0] = xy[i_,0];
            }
            for(i_=0; i_<=n-1;i_++)
            {
                xyp[i_,1] = xy[i_,1];
            }
            for(i_=0; i_<=n-1;i_++)
            {
                xyp[i_,2] = xy[i_,2];
            }
            for(i_=0; i_<=2;i_++)
            {
                xyp[n,i_] = xy[0,i_];
            }
            
            //
            // Build parameterization, check that all parameters are distinct
            //
            pspline3par(ref xyp, n+1, pt, ref p.p);
            System.Diagnostics.Debug.Assert(apserv.apservaredistinct(p.p, n+1), "PSplineBuild2Periodic: consequent (or first and last) points are too close!");
            
            //
            // Build splines
            //
            if( st==1 )
            {
                for(i_=0; i_<=n;i_++)
                {
                    tmp[i_] = xyp[i_,0];
                }
                spline1d.spline1dbuildcatmullrom(p.p, tmp, n+1, -1, 0.0, ref p.x);
                for(i_=0; i_<=n;i_++)
                {
                    tmp[i_] = xyp[i_,1];
                }
                spline1d.spline1dbuildcatmullrom(p.p, tmp, n+1, -1, 0.0, ref p.y);
                for(i_=0; i_<=n;i_++)
                {
                    tmp[i_] = xyp[i_,2];
                }
                spline1d.spline1dbuildcatmullrom(p.p, tmp, n+1, -1, 0.0, ref p.z);
            }
            if( st==2 )
            {
                for(i_=0; i_<=n;i_++)
                {
                    tmp[i_] = xyp[i_,0];
                }
                spline1d.spline1dbuildcubic(p.p, tmp, n+1, -1, 0.0, -1, 0.0, ref p.x);
                for(i_=0; i_<=n;i_++)
                {
                    tmp[i_] = xyp[i_,1];
                }
                spline1d.spline1dbuildcubic(p.p, tmp, n+1, -1, 0.0, -1, 0.0, ref p.y);
                for(i_=0; i_<=n;i_++)
                {
                    tmp[i_] = xyp[i_,2];
                }
                spline1d.spline1dbuildcubic(p.p, tmp, n+1, -1, 0.0, -1, 0.0, ref p.z);
            }
        }


        /*************************************************************************
        This function returns vector of parameter values correspoding to points.

        I.e. for P created from (X[0],Y[0])...(X[N-1],Y[N-1]) and U=TValues(P)  we
        have
            (X[0],Y[0]) = PSpline2Calc(P,U[0]),
            (X[1],Y[1]) = PSpline2Calc(P,U[1]),
            (X[2],Y[2]) = PSpline2Calc(P,U[2]),
            ...

        INPUT PARAMETERS:
            P   -   parametric spline interpolant

        OUTPUT PARAMETERS:
            N   -   array size
            T   -   array[0..N-1]


        NOTES:
        * for non-periodic splines U[0]=0, U[0]<U[1]<...<U[N-1], U[N-1]=1
        * for periodic splines     U[0]=0, U[0]<U[1]<...<U[N-1], U[N-1]<1

          -- ALGLIB PROJECT --
             Copyright 28.05.2010 by Bochkanov Sergey
        *************************************************************************/
        public static void pspline2parametervalues(ref pspline2interpolant p,
            ref int n,
            ref double[] t)
        {
            int i_ = 0;

            System.Diagnostics.Debug.Assert(p.n>=2, "PSpline2ParameterValues: internal error!");
            n = p.n;
            t = new double[n];
            for(i_=0; i_<=n-1;i_++)
            {
                t[i_] = p.p[i_];
            }
            t[0] = 0;
            if( !p.periodic )
            {
                t[n-1] = 1;
            }
        }


        /*************************************************************************
        This function returns vector of parameter values correspoding to points.

        Same as PSpline2ParameterValues(), but for 3D.

          -- ALGLIB PROJECT --
             Copyright 28.05.2010 by Bochkanov Sergey
        *************************************************************************/
        public static void pspline3parametervalues(ref pspline3interpolant p,
            ref int n,
            ref double[] t)
        {
            int i_ = 0;

            System.Diagnostics.Debug.Assert(p.n>=2, "PSpline3ParameterValues: internal error!");
            n = p.n;
            t = new double[n];
            for(i_=0; i_<=n-1;i_++)
            {
                t[i_] = p.p[i_];
            }
            t[0] = 0;
            if( !p.periodic )
            {
                t[n-1] = 1;
            }
        }


        /*************************************************************************
        This function  calculates  the value of the parametric spline for a  given
        value of parameter T

        INPUT PARAMETERS:
            P   -   parametric spline interpolant
            T   -   point:
                    * T in [0,1] corresponds to interval spanned by points
                    * for non-periodic splines T<0 (or T>1) correspond to parts of
                      the curve before the first (after the last) point
                    * for periodic splines T<0 (or T>1) are projected  into  [0,1]
                      by making T=T-floor(T).

        OUTPUT PARAMETERS:
            X   -   X-position
            Y   -   Y-position


          -- ALGLIB PROJECT --
             Copyright 28.05.2010 by Bochkanov Sergey
        *************************************************************************/
        public static void pspline2calc(ref pspline2interpolant p,
            double t,
            ref double x,
            ref double y)
        {
            if( p.periodic )
            {
                t = t-(int)Math.Floor(t);
            }
            x = spline1d.spline1dcalc(ref p.x, t);
            y = spline1d.spline1dcalc(ref p.y, t);
        }


        /*************************************************************************
        This function  calculates  the value of the parametric spline for a  given
        value of parameter T.

        INPUT PARAMETERS:
            P   -   parametric spline interpolant
            T   -   point:
                    * T in [0,1] corresponds to interval spanned by points
                    * for non-periodic splines T<0 (or T>1) correspond to parts of
                      the curve before the first (after the last) point
                    * for periodic splines T<0 (or T>1) are projected  into  [0,1]
                      by making T=T-floor(T).

        OUTPUT PARAMETERS:
            X   -   X-position
            Y   -   Y-position
            Z   -   Z-position


          -- ALGLIB PROJECT --
             Copyright 28.05.2010 by Bochkanov Sergey
        *************************************************************************/
        public static void pspline3calc(ref pspline3interpolant p,
            double t,
            ref double x,
            ref double y,
            ref double z)
        {
            if( p.periodic )
            {
                t = t-(int)Math.Floor(t);
            }
            x = spline1d.spline1dcalc(ref p.x, t);
            y = spline1d.spline1dcalc(ref p.y, t);
            z = spline1d.spline1dcalc(ref p.z, t);
        }


        /*************************************************************************
        This function  calculates  tangent vector for a given value of parameter T

        INPUT PARAMETERS:
            P   -   parametric spline interpolant
            T   -   point:
                    * T in [0,1] corresponds to interval spanned by points
                    * for non-periodic splines T<0 (or T>1) correspond to parts of
                      the curve before the first (after the last) point
                    * for periodic splines T<0 (or T>1) are projected  into  [0,1]
                      by making T=T-floor(T).

        OUTPUT PARAMETERS:
            X    -   X-component of tangent vector (normalized)
            Y    -   Y-component of tangent vector (normalized)
            
        NOTE:
            X^2+Y^2 is either 1 (for non-zero tangent vector) or 0.


          -- ALGLIB PROJECT --
             Copyright 28.05.2010 by Bochkanov Sergey
        *************************************************************************/
        public static void pspline2tangent(ref pspline2interpolant p,
            double t,
            ref double x,
            ref double y)
        {
            double v = 0;
            double v0 = 0;
            double v1 = 0;

            if( p.periodic )
            {
                t = t-(int)Math.Floor(t);
            }
            pspline2diff(ref p, t, ref v0, ref x, ref v1, ref y);
            if( (double)(x)!=(double)(0) | (double)(y)!=(double)(0) )
            {
                
                //
                // this code is a bit more complex than X^2+Y^2 to avoid
                // overflow for large values of X and Y.
                //
                v = apserv.safepythag2(x, y);
                x = x/v;
                y = y/v;
            }
        }


        /*************************************************************************
        This function  calculates  tangent vector for a given value of parameter T

        INPUT PARAMETERS:
            P   -   parametric spline interpolant
            T   -   point:
                    * T in [0,1] corresponds to interval spanned by points
                    * for non-periodic splines T<0 (or T>1) correspond to parts of
                      the curve before the first (after the last) point
                    * for periodic splines T<0 (or T>1) are projected  into  [0,1]
                      by making T=T-floor(T).

        OUTPUT PARAMETERS:
            X    -   X-component of tangent vector (normalized)
            Y    -   Y-component of tangent vector (normalized)
            Z    -   Z-component of tangent vector (normalized)

        NOTE:
            X^2+Y^2+Z^2 is either 1 (for non-zero tangent vector) or 0.


          -- ALGLIB PROJECT --
             Copyright 28.05.2010 by Bochkanov Sergey
        *************************************************************************/
        public static void pspline3tangent(ref pspline3interpolant p,
            double t,
            ref double x,
            ref double y,
            ref double z)
        {
            double v = 0;
            double v0 = 0;
            double v1 = 0;
            double v2 = 0;

            if( p.periodic )
            {
                t = t-(int)Math.Floor(t);
            }
            pspline3diff(ref p, t, ref v0, ref x, ref v1, ref y, ref v2, ref z);
            if( (double)(x)!=(double)(0) | (double)(y)!=(double)(0) | (double)(z)!=(double)(0) )
            {
                v = apserv.safepythag3(x, y, z);
                x = x/v;
                y = y/v;
                z = z/v;
            }
        }


        /*************************************************************************
        This function calculates derivative, i.e. it returns (dX/dT,dY/dT).

        INPUT PARAMETERS:
            P   -   parametric spline interpolant
            T   -   point:
                    * T in [0,1] corresponds to interval spanned by points
                    * for non-periodic splines T<0 (or T>1) correspond to parts of
                      the curve before the first (after the last) point
                    * for periodic splines T<0 (or T>1) are projected  into  [0,1]
                      by making T=T-floor(T).

        OUTPUT PARAMETERS:
            X   -   X-value
            DX  -   X-derivative
            Y   -   Y-value
            DY  -   Y-derivative


          -- ALGLIB PROJECT --
             Copyright 28.05.2010 by Bochkanov Sergey
        *************************************************************************/
        public static void pspline2diff(ref pspline2interpolant p,
            double t,
            ref double x,
            ref double dx,
            ref double y,
            ref double dy)
        {
            double d2s = 0;

            if( p.periodic )
            {
                t = t-(int)Math.Floor(t);
            }
            spline1d.spline1ddiff(ref p.x, t, ref x, ref dx, ref d2s);
            spline1d.spline1ddiff(ref p.y, t, ref y, ref dy, ref d2s);
        }


        /*************************************************************************
        This function calculates derivative, i.e. it returns (dX/dT,dY/dT,dZ/dT).

        INPUT PARAMETERS:
            P   -   parametric spline interpolant
            T   -   point:
                    * T in [0,1] corresponds to interval spanned by points
                    * for non-periodic splines T<0 (or T>1) correspond to parts of
                      the curve before the first (after the last) point
                    * for periodic splines T<0 (or T>1) are projected  into  [0,1]
                      by making T=T-floor(T).

        OUTPUT PARAMETERS:
            X   -   X-value
            DX  -   X-derivative
            Y   -   Y-value
            DY  -   Y-derivative
            Z   -   Z-value
            DZ  -   Z-derivative


          -- ALGLIB PROJECT --
             Copyright 28.05.2010 by Bochkanov Sergey
        *************************************************************************/
        public static void pspline3diff(ref pspline3interpolant p,
            double t,
            ref double x,
            ref double dx,
            ref double y,
            ref double dy,
            ref double z,
            ref double dz)
        {
            double d2s = 0;

            if( p.periodic )
            {
                t = t-(int)Math.Floor(t);
            }
            spline1d.spline1ddiff(ref p.x, t, ref x, ref dx, ref d2s);
            spline1d.spline1ddiff(ref p.y, t, ref y, ref dy, ref d2s);
            spline1d.spline1ddiff(ref p.z, t, ref z, ref dz, ref d2s);
        }


        /*************************************************************************
        This function calculates first and second derivative with respect to T.

        INPUT PARAMETERS:
            P   -   parametric spline interpolant
            T   -   point:
                    * T in [0,1] corresponds to interval spanned by points
                    * for non-periodic splines T<0 (or T>1) correspond to parts of
                      the curve before the first (after the last) point
                    * for periodic splines T<0 (or T>1) are projected  into  [0,1]
                      by making T=T-floor(T).

        OUTPUT PARAMETERS:
            X   -   X-value
            DX  -   derivative
            D2X -   second derivative
            Y   -   Y-value
            DY  -   derivative
            D2Y -   second derivative


          -- ALGLIB PROJECT --
             Copyright 28.05.2010 by Bochkanov Sergey
        *************************************************************************/
        public static void pspline2diff2(ref pspline2interpolant p,
            double t,
            ref double x,
            ref double dx,
            ref double d2x,
            ref double y,
            ref double dy,
            ref double d2y)
        {
            if( p.periodic )
            {
                t = t-(int)Math.Floor(t);
            }
            spline1d.spline1ddiff(ref p.x, t, ref x, ref dx, ref d2x);
            spline1d.spline1ddiff(ref p.y, t, ref y, ref dy, ref d2y);
        }


        /*************************************************************************
        This function calculates first and second derivative with respect to T.

        INPUT PARAMETERS:
            P   -   parametric spline interpolant
            T   -   point:
                    * T in [0,1] corresponds to interval spanned by points
                    * for non-periodic splines T<0 (or T>1) correspond to parts of
                      the curve before the first (after the last) point
                    * for periodic splines T<0 (or T>1) are projected  into  [0,1]
                      by making T=T-floor(T).

        OUTPUT PARAMETERS:
            X   -   X-value
            DX  -   derivative
            D2X -   second derivative
            Y   -   Y-value
            DY  -   derivative
            D2Y -   second derivative
            Z   -   Z-value
            DZ  -   derivative
            D2Z -   second derivative


          -- ALGLIB PROJECT --
             Copyright 28.05.2010 by Bochkanov Sergey
        *************************************************************************/
        public static void pspline3diff2(ref pspline3interpolant p,
            double t,
            ref double x,
            ref double dx,
            ref double d2x,
            ref double y,
            ref double dy,
            ref double d2y,
            ref double z,
            ref double dz,
            ref double d2z)
        {
            if( p.periodic )
            {
                t = t-(int)Math.Floor(t);
            }
            spline1d.spline1ddiff(ref p.x, t, ref x, ref dx, ref d2x);
            spline1d.spline1ddiff(ref p.y, t, ref y, ref dy, ref d2y);
            spline1d.spline1ddiff(ref p.z, t, ref z, ref dz, ref d2z);
        }


        /*************************************************************************
        This function  calculates  arc length, i.e. length of  curve  between  t=a
        and t=b.

        INPUT PARAMETERS:
            P   -   parametric spline interpolant
            A,B -   parameter values corresponding to arc ends:
                    * B>A will result in positive length returned
                    * B<A will result in negative length returned

        RESULT:
            length of arc starting at T=A and ending at T=B.


          -- ALGLIB PROJECT --
             Copyright 30.05.2010 by Bochkanov Sergey
        *************************************************************************/
        public static double pspline2arclength(ref pspline2interpolant p,
            double a,
            double b)
        {
            double result = 0;
            autogk.autogkstate state = new autogk.autogkstate();
            autogk.autogkreport rep = new autogk.autogkreport();
            double sx = 0;
            double dsx = 0;
            double d2sx = 0;
            double sy = 0;
            double dsy = 0;
            double d2sy = 0;

            autogk.autogksmooth(a, b, ref state);
            while( autogk.autogkiteration(ref state) )
            {
                spline1d.spline1ddiff(ref p.x, state.x, ref sx, ref dsx, ref d2sx);
                spline1d.spline1ddiff(ref p.y, state.x, ref sy, ref dsy, ref d2sy);
                state.f = apserv.safepythag2(dsx, dsy);
            }
            autogk.autogkresults(ref state, ref result, ref rep);
            System.Diagnostics.Debug.Assert(rep.terminationtype>0, "PSpline2ArcLength: internal error!");
            return result;
        }


        /*************************************************************************
        This function  calculates  arc length, i.e. length of  curve  between  t=a
        and t=b.

        INPUT PARAMETERS:
            P   -   parametric spline interpolant
            A,B -   parameter values corresponding to arc ends:
                    * B>A will result in positive length returned
                    * B<A will result in negative length returned

        RESULT:
            length of arc starting at T=A and ending at T=B.


          -- ALGLIB PROJECT --
             Copyright 30.05.2010 by Bochkanov Sergey
        *************************************************************************/
        public static double pspline3arclength(ref pspline3interpolant p,
            double a,
            double b)
        {
            double result = 0;
            autogk.autogkstate state = new autogk.autogkstate();
            autogk.autogkreport rep = new autogk.autogkreport();
            double sx = 0;
            double dsx = 0;
            double d2sx = 0;
            double sy = 0;
            double dsy = 0;
            double d2sy = 0;
            double sz = 0;
            double dsz = 0;
            double d2sz = 0;

            autogk.autogksmooth(a, b, ref state);
            while( autogk.autogkiteration(ref state) )
            {
                spline1d.spline1ddiff(ref p.x, state.x, ref sx, ref dsx, ref d2sx);
                spline1d.spline1ddiff(ref p.y, state.x, ref sy, ref dsy, ref d2sy);
                spline1d.spline1ddiff(ref p.z, state.x, ref sz, ref dsz, ref d2sz);
                state.f = apserv.safepythag3(dsx, dsy, dsz);
            }
            autogk.autogkresults(ref state, ref result, ref rep);
            System.Diagnostics.Debug.Assert(rep.terminationtype>0, "PSpline3ArcLength: internal error!");
            return result;
        }


        /*************************************************************************
        Builds non-periodic parameterization for 2-dimensional spline
        *************************************************************************/
        private static void pspline2par(ref double[,] xy,
            int n,
            int pt,
            ref double[] p)
        {
            double v = 0;
            int i = 0;
            int i_ = 0;

            System.Diagnostics.Debug.Assert(pt>=0 & pt<=2, "PSpline2Par: internal error!");
            
            //
            // Build parameterization:
            // * fill by non-normalized values
            // * normalize them so we have P[0]=0, P[N-1]=1.
            //
            p = new double[n];
            if( pt==0 )
            {
                for(i=0; i<=n-1; i++)
                {
                    p[i] = i;
                }
            }
            if( pt==1 )
            {
                p[0] = 0;
                for(i=1; i<=n-1; i++)
                {
                    p[i] = p[i-1]+apserv.safepythag2(xy[i,0]-xy[i-1,0], xy[i,1]-xy[i-1,1]);
                }
            }
            if( pt==2 )
            {
                p[0] = 0;
                for(i=1; i<=n-1; i++)
                {
                    p[i] = p[i-1]+Math.Sqrt(apserv.safepythag2(xy[i,0]-xy[i-1,0], xy[i,1]-xy[i-1,1]));
                }
            }
            v = 1/p[n-1];
            for(i_=0; i_<=n-1;i_++)
            {
                p[i_] = v*p[i_];
            }
        }


        /*************************************************************************
        Builds non-periodic parameterization for 3-dimensional spline
        *************************************************************************/
        private static void pspline3par(ref double[,] xy,
            int n,
            int pt,
            ref double[] p)
        {
            double v = 0;
            int i = 0;
            int i_ = 0;

            System.Diagnostics.Debug.Assert(pt>=0 & pt<=2, "PSpline3Par: internal error!");
            
            //
            // Build parameterization:
            // * fill by non-normalized values
            // * normalize them so we have P[0]=0, P[N-1]=1.
            //
            p = new double[n];
            if( pt==0 )
            {
                for(i=0; i<=n-1; i++)
                {
                    p[i] = i;
                }
            }
            if( pt==1 )
            {
                p[0] = 0;
                for(i=1; i<=n-1; i++)
                {
                    p[i] = p[i-1]+apserv.safepythag3(xy[i,0]-xy[i-1,0], xy[i,1]-xy[i-1,1], xy[i,2]-xy[i-1,2]);
                }
            }
            if( pt==2 )
            {
                p[0] = 0;
                for(i=1; i<=n-1; i++)
                {
                    p[i] = p[i-1]+Math.Sqrt(apserv.safepythag3(xy[i,0]-xy[i-1,0], xy[i,1]-xy[i-1,1], xy[i,2]-xy[i-1,2]));
                }
            }
            v = 1/p[n-1];
            for(i_=0; i_<=n-1;i_++)
            {
                p[i_] = v*p[i_];
            }
        }
    }
}
