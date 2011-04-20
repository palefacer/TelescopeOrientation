
using System;

namespace alglib
{
    public class testpsplineunit
    {
        public static bool testpsplineinterpolation(bool silent)
        {
            bool result = new bool();
            bool waserrors = new bool();
            bool p2errors = new bool();
            bool p3errors = new bool();
            double nonstrictthreshold = 0;
            double threshold = 0;
            int passcount = 0;
            double lstep = 0;
            double h = 0;
            int maxn = 0;
            int periodicity = 0;
            int skind = 0;
            int pkind = 0;
            bool periodic = new bool();
            double a = 0;
            double b = 0;
            int n = 0;
            int tmpn = 0;
            int i = 0;
            int k = 0;
            double vx = 0;
            double vy = 0;
            double vz = 0;
            double vx2 = 0;
            double vy2 = 0;
            double vz2 = 0;
            double vdx = 0;
            double vdy = 0;
            double vdz = 0;
            double vdx2 = 0;
            double vdy2 = 0;
            double vdz2 = 0;
            double vd2x = 0;
            double vd2y = 0;
            double vd2z = 0;
            double vd2x2 = 0;
            double vd2y2 = 0;
            double vd2z2 = 0;
            double v0 = 0;
            double v1 = 0;
            double[] x = new double[0];
            double[] y = new double[0];
            double[] z = new double[0];
            double[] t = new double[0];
            double[] t2 = new double[0];
            double[] t3 = new double[0];
            double[,] xy = new double[0,0];
            double[,] xyz = new double[0,0];
            pspline.pspline2interpolant p2 = new pspline.pspline2interpolant();
            pspline.pspline3interpolant p3 = new pspline.pspline3interpolant();
            spline1d.spline1dinterpolant s = new spline1d.spline1dinterpolant();
            int i_ = 0;

            waserrors = false;
            passcount = 20;
            lstep = 0.005;
            h = 0.00001;
            maxn = 10;
            threshold = 10000*AP.Math.MachineEpsilon;
            nonstrictthreshold = 0.00001;
            p2errors = false;
            p3errors = false;
            
            //
            // Test basic properties of 2- and 3-dimensional splines:
            // * PSpline2ParameterValues() properties
            // * values at nodes
            // * for periodic splines - periodicity properties
            //
            // Variables used:
            // * N              points count
            // * SKind          spline
            // * PKind          parameterization
            // * Periodicity    whether we have periodic spline or not
            //
            for(n=2; n<=maxn; n++)
            {
                for(skind=0; skind<=2; skind++)
                {
                    for(pkind=0; pkind<=2; pkind++)
                    {
                        for(periodicity=0; periodicity<=1; periodicity++)
                        {
                            periodic = periodicity==1;
                            
                            //
                            // skip unsupported combinations of parameters
                            //
                            if( periodic & n<3 )
                            {
                                continue;
                            }
                            if( periodic & skind==0 )
                            {
                                continue;
                            }
                            if( n<5 & skind==0 )
                            {
                                continue;
                            }
                            
                            //
                            // init
                            //
                            xy = new double[n, 2];
                            xyz = new double[n, 3];
                            apserv.taskgenint1dequidist(-1, +1, n, ref t2, ref x);
                            for(i_=0; i_<=n-1;i_++)
                            {
                                xy[i_,0] = x[i_];
                            }
                            for(i_=0; i_<=n-1;i_++)
                            {
                                xyz[i_,0] = x[i_];
                            }
                            apserv.taskgenint1dequidist(-1, +1, n, ref t2, ref y);
                            for(i_=0; i_<=n-1;i_++)
                            {
                                xy[i_,1] = y[i_];
                            }
                            for(i_=0; i_<=n-1;i_++)
                            {
                                xyz[i_,1] = y[i_];
                            }
                            apserv.taskgenint1dequidist(-1, +1, n, ref t2, ref z);
                            for(i_=0; i_<=n-1;i_++)
                            {
                                xyz[i_,2] = z[i_];
                            }
                            unsetp2(ref p2);
                            unsetp3(ref p3);
                            if( periodic )
                            {
                                pspline.pspline2buildperiodic(xy, n, skind, pkind, ref p2);
                                pspline.pspline3buildperiodic(xyz, n, skind, pkind, ref p3);
                            }
                            else
                            {
                                pspline.pspline2build(xy, n, skind, pkind, ref p2);
                                pspline.pspline3build(xyz, n, skind, pkind, ref p3);
                            }
                            
                            //
                            // PSpline2ParameterValues() properties
                            //
                            pspline.pspline2parametervalues(ref p2, ref tmpn, ref t2);
                            if( tmpn!=n )
                            {
                                p2errors = true;
                                continue;
                            }
                            pspline.pspline3parametervalues(ref p3, ref tmpn, ref t3);
                            if( tmpn!=n )
                            {
                                p3errors = true;
                                continue;
                            }
                            p2errors = p2errors | (double)(t2[0])!=(double)(0);
                            p3errors = p3errors | (double)(t3[0])!=(double)(0);
                            for(i=1; i<=n-1; i++)
                            {
                                p2errors = p2errors | (double)(t2[i])<=(double)(t2[i-1]);
                                p3errors = p3errors | (double)(t3[i])<=(double)(t3[i-1]);
                            }
                            if( periodic )
                            {
                                p2errors = p2errors | (double)(t2[n-1])>=(double)(1);
                                p3errors = p3errors | (double)(t3[n-1])>=(double)(1);
                            }
                            else
                            {
                                p2errors = p2errors | (double)(t2[n-1])!=(double)(1);
                                p3errors = p3errors | (double)(t3[n-1])!=(double)(1);
                            }
                            
                            //
                            // Now we have parameter values stored at T,
                            // and want to test whether the actully correspond to
                            // points
                            //
                            for(i=0; i<=n-1; i++)
                            {
                                
                                //
                                // 2-dimensional test
                                //
                                pspline.pspline2calc(ref p2, t2[i], ref vx, ref vy);
                                p2errors = p2errors | (double)(Math.Abs(vx-x[i]))>(double)(threshold);
                                p2errors = p2errors | (double)(Math.Abs(vy-y[i]))>(double)(threshold);
                                
                                //
                                // 3-dimensional test
                                //
                                pspline.pspline3calc(ref p3, t3[i], ref vx, ref vy, ref vz);
                                p3errors = p3errors | (double)(Math.Abs(vx-x[i]))>(double)(threshold);
                                p3errors = p3errors | (double)(Math.Abs(vy-y[i]))>(double)(threshold);
                                p3errors = p3errors | (double)(Math.Abs(vz-z[i]))>(double)(threshold);
                            }
                            
                            //
                            // Test periodicity (if needed)
                            //
                            if( periodic )
                            {
                                
                                //
                                // periodicity at nodes
                                //
                                for(i=0; i<=n-1; i++)
                                {
                                    
                                    //
                                    // 2-dimensional test
                                    //
                                    pspline.pspline2calc(ref p2, t2[i]+AP.Math.RandomInteger(10)-5, ref vx, ref vy);
                                    p2errors = p2errors | (double)(Math.Abs(vx-x[i]))>(double)(threshold);
                                    p2errors = p2errors | (double)(Math.Abs(vy-y[i]))>(double)(threshold);
                                    pspline.pspline2diff(ref p2, t2[i]+AP.Math.RandomInteger(10)-5, ref vx, ref vdx, ref vy, ref vdy);
                                    p2errors = p2errors | (double)(Math.Abs(vx-x[i]))>(double)(threshold);
                                    p2errors = p2errors | (double)(Math.Abs(vy-y[i]))>(double)(threshold);
                                    pspline.pspline2diff2(ref p2, t2[i]+AP.Math.RandomInteger(10)-5, ref vx, ref vdx, ref vd2x, ref vy, ref vdy, ref vd2y);
                                    p2errors = p2errors | (double)(Math.Abs(vx-x[i]))>(double)(threshold);
                                    p2errors = p2errors | (double)(Math.Abs(vy-y[i]))>(double)(threshold);
                                    
                                    //
                                    // 3-dimensional test
                                    //
                                    pspline.pspline3calc(ref p3, t3[i]+AP.Math.RandomInteger(10)-5, ref vx, ref vy, ref vz);
                                    p3errors = p3errors | (double)(Math.Abs(vx-x[i]))>(double)(threshold);
                                    p3errors = p3errors | (double)(Math.Abs(vy-y[i]))>(double)(threshold);
                                    p3errors = p3errors | (double)(Math.Abs(vz-z[i]))>(double)(threshold);
                                    pspline.pspline3diff(ref p3, t3[i]+AP.Math.RandomInteger(10)-5, ref vx, ref vdx, ref vy, ref vdy, ref vz, ref vdz);
                                    p3errors = p3errors | (double)(Math.Abs(vx-x[i]))>(double)(threshold);
                                    p3errors = p3errors | (double)(Math.Abs(vy-y[i]))>(double)(threshold);
                                    p3errors = p3errors | (double)(Math.Abs(vz-z[i]))>(double)(threshold);
                                    pspline.pspline3diff2(ref p3, t3[i]+AP.Math.RandomInteger(10)-5, ref vx, ref vdx, ref vd2x, ref vy, ref vdy, ref vd2y, ref vz, ref vdz, ref vd2z);
                                    p3errors = p3errors | (double)(Math.Abs(vx-x[i]))>(double)(threshold);
                                    p3errors = p3errors | (double)(Math.Abs(vy-y[i]))>(double)(threshold);
                                    p3errors = p3errors | (double)(Math.Abs(vz-z[i]))>(double)(threshold);
                                }
                                
                                //
                                // periodicity between nodes
                                //
                                v0 = AP.Math.RandomReal();
                                pspline.pspline2calc(ref p2, v0, ref vx, ref vy);
                                pspline.pspline2calc(ref p2, v0+AP.Math.RandomInteger(10)-5, ref vx2, ref vy2);
                                p2errors = p2errors | (double)(Math.Abs(vx-vx2))>(double)(threshold);
                                p2errors = p2errors | (double)(Math.Abs(vy-vy2))>(double)(threshold);
                                pspline.pspline3calc(ref p3, v0, ref vx, ref vy, ref vz);
                                pspline.pspline3calc(ref p3, v0+AP.Math.RandomInteger(10)-5, ref vx2, ref vy2, ref vz2);
                                p3errors = p3errors | (double)(Math.Abs(vx-vx2))>(double)(threshold);
                                p3errors = p3errors | (double)(Math.Abs(vy-vy2))>(double)(threshold);
                                p3errors = p3errors | (double)(Math.Abs(vz-vz2))>(double)(threshold);
                                
                                //
                                // near-boundary test for continuity of function values and derivatives:
                                // 2-dimensional curve
                                //
                                System.Diagnostics.Debug.Assert(skind==1 | skind==2, "TEST: unexpected spline type!");
                                v0 = 100*AP.Math.MachineEpsilon;
                                v1 = 1-v0;
                                pspline.pspline2calc(ref p2, v0, ref vx, ref vy);
                                pspline.pspline2calc(ref p2, v1, ref vx2, ref vy2);
                                p2errors = p2errors | (double)(Math.Abs(vx-vx2))>(double)(threshold);
                                p2errors = p2errors | (double)(Math.Abs(vy-vy2))>(double)(threshold);
                                pspline.pspline2diff(ref p2, v0, ref vx, ref vdx, ref vy, ref vdy);
                                pspline.pspline2diff(ref p2, v1, ref vx2, ref vdx2, ref vy2, ref vdy2);
                                p2errors = p2errors | (double)(Math.Abs(vx-vx2))>(double)(threshold);
                                p2errors = p2errors | (double)(Math.Abs(vy-vy2))>(double)(threshold);
                                p2errors = p2errors | (double)(Math.Abs(vdx-vdx2))>(double)(nonstrictthreshold);
                                p2errors = p2errors | (double)(Math.Abs(vdy-vdy2))>(double)(nonstrictthreshold);
                                pspline.pspline2diff2(ref p2, v0, ref vx, ref vdx, ref vd2x, ref vy, ref vdy, ref vd2y);
                                pspline.pspline2diff2(ref p2, v1, ref vx2, ref vdx2, ref vd2x2, ref vy2, ref vdy2, ref vd2y2);
                                p2errors = p2errors | (double)(Math.Abs(vx-vx2))>(double)(threshold);
                                p2errors = p2errors | (double)(Math.Abs(vy-vy2))>(double)(threshold);
                                p2errors = p2errors | (double)(Math.Abs(vdx-vdx2))>(double)(nonstrictthreshold);
                                p2errors = p2errors | (double)(Math.Abs(vdy-vdy2))>(double)(nonstrictthreshold);
                                if( skind==2 )
                                {
                                    
                                    //
                                    // second derivative test only for cubic splines
                                    //
                                    p2errors = p2errors | (double)(Math.Abs(vd2x-vd2x2))>(double)(nonstrictthreshold);
                                    p2errors = p2errors | (double)(Math.Abs(vd2y-vd2y2))>(double)(nonstrictthreshold);
                                }
                                
                                //
                                // near-boundary test for continuity of function values and derivatives:
                                // 3-dimensional curve
                                //
                                System.Diagnostics.Debug.Assert(skind==1 | skind==2, "TEST: unexpected spline type!");
                                v0 = 100*AP.Math.MachineEpsilon;
                                v1 = 1-v0;
                                pspline.pspline3calc(ref p3, v0, ref vx, ref vy, ref vz);
                                pspline.pspline3calc(ref p3, v1, ref vx2, ref vy2, ref vz2);
                                p3errors = p3errors | (double)(Math.Abs(vx-vx2))>(double)(threshold);
                                p3errors = p3errors | (double)(Math.Abs(vy-vy2))>(double)(threshold);
                                p3errors = p3errors | (double)(Math.Abs(vz-vz2))>(double)(threshold);
                                pspline.pspline3diff(ref p3, v0, ref vx, ref vdx, ref vy, ref vdy, ref vz, ref vdz);
                                pspline.pspline3diff(ref p3, v1, ref vx2, ref vdx2, ref vy2, ref vdy2, ref vz2, ref vdz2);
                                p3errors = p3errors | (double)(Math.Abs(vx-vx2))>(double)(threshold);
                                p3errors = p3errors | (double)(Math.Abs(vy-vy2))>(double)(threshold);
                                p3errors = p3errors | (double)(Math.Abs(vz-vz2))>(double)(threshold);
                                p3errors = p3errors | (double)(Math.Abs(vdx-vdx2))>(double)(nonstrictthreshold);
                                p3errors = p3errors | (double)(Math.Abs(vdy-vdy2))>(double)(nonstrictthreshold);
                                p3errors = p3errors | (double)(Math.Abs(vdz-vdz2))>(double)(nonstrictthreshold);
                                pspline.pspline3diff2(ref p3, v0, ref vx, ref vdx, ref vd2x, ref vy, ref vdy, ref vd2y, ref vz, ref vdz, ref vd2z);
                                pspline.pspline3diff2(ref p3, v1, ref vx2, ref vdx2, ref vd2x2, ref vy2, ref vdy2, ref vd2y2, ref vz2, ref vdz2, ref vd2z2);
                                p3errors = p3errors | (double)(Math.Abs(vx-vx2))>(double)(threshold);
                                p3errors = p3errors | (double)(Math.Abs(vy-vy2))>(double)(threshold);
                                p3errors = p3errors | (double)(Math.Abs(vz-vz2))>(double)(threshold);
                                p3errors = p3errors | (double)(Math.Abs(vdx-vdx2))>(double)(nonstrictthreshold);
                                p3errors = p3errors | (double)(Math.Abs(vdy-vdy2))>(double)(nonstrictthreshold);
                                p3errors = p3errors | (double)(Math.Abs(vdz-vdz2))>(double)(nonstrictthreshold);
                                if( skind==2 )
                                {
                                    
                                    //
                                    // second derivative test only for cubic splines
                                    //
                                    p3errors = p3errors | (double)(Math.Abs(vd2x-vd2x2))>(double)(nonstrictthreshold);
                                    p3errors = p3errors | (double)(Math.Abs(vd2y-vd2y2))>(double)(nonstrictthreshold);
                                    p3errors = p3errors | (double)(Math.Abs(vd2z-vd2z2))>(double)(nonstrictthreshold);
                                }
                            }
                        }
                    }
                }
            }
            
            //
            // Test differentiation, tangents, calculation between nodes.
            //
            // Because differentiation is done in parameterization/spline/periodicity
            // oblivious manner, we don't have to test all possible combinations
            // of spline types and parameterizations.
            //
            // Actually we test special combination with properties which allow us
            // to easily solve this problem:
            // * 2 (3) variables
            // * first variable is sampled from equidistant grid on [0,1]
            // * other variables are random
            // * uniform parameterization is used
            // * periodicity - none
            // * spline type - any (we use cubic splines)
            // Same problem allows us to test calculation BETWEEN nodes.
            //
            for(n=2; n<=maxn; n++)
            {
                
                //
                // init
                //
                xy = new double[n, 2];
                xyz = new double[n, 3];
                apserv.taskgenint1dequidist(0, +1, n, ref t, ref x);
                for(i_=0; i_<=n-1;i_++)
                {
                    xy[i_,0] = x[i_];
                }
                for(i_=0; i_<=n-1;i_++)
                {
                    xyz[i_,0] = x[i_];
                }
                apserv.taskgenint1dequidist(0, +1, n, ref t, ref y);
                for(i_=0; i_<=n-1;i_++)
                {
                    xy[i_,1] = y[i_];
                }
                for(i_=0; i_<=n-1;i_++)
                {
                    xyz[i_,1] = y[i_];
                }
                apserv.taskgenint1dequidist(0, +1, n, ref t, ref z);
                for(i_=0; i_<=n-1;i_++)
                {
                    xyz[i_,2] = z[i_];
                }
                unsetp2(ref p2);
                unsetp3(ref p3);
                pspline.pspline2build(xy, n, 2, 0, ref p2);
                pspline.pspline3build(xyz, n, 2, 0, ref p3);
                
                //
                // Test 2D/3D spline:
                // * build non-parametric cubic spline from T and X/Y
                // * calculate its value and derivatives at V0
                // * compare with Spline2Calc/Spline2Diff/Spline2Diff2
                // Because of task properties both variants should
                // return same answer.
                //
                v0 = AP.Math.RandomReal();
                spline1d.spline1dbuildcubic(t, x, n, 0, 0.0, 0, 0.0, ref s);
                spline1d.spline1ddiff(ref s, v0, ref vx2, ref vdx2, ref vd2x2);
                spline1d.spline1dbuildcubic(t, y, n, 0, 0.0, 0, 0.0, ref s);
                spline1d.spline1ddiff(ref s, v0, ref vy2, ref vdy2, ref vd2y2);
                spline1d.spline1dbuildcubic(t, z, n, 0, 0.0, 0, 0.0, ref s);
                spline1d.spline1ddiff(ref s, v0, ref vz2, ref vdz2, ref vd2z2);
                
                //
                // 2D test
                //
                pspline.pspline2calc(ref p2, v0, ref vx, ref vy);
                p2errors = p2errors | (double)(Math.Abs(vx-vx2))>(double)(threshold);
                p2errors = p2errors | (double)(Math.Abs(vy-vy2))>(double)(threshold);
                pspline.pspline2diff(ref p2, v0, ref vx, ref vdx, ref vy, ref vdy);
                p2errors = p2errors | (double)(Math.Abs(vx-vx2))>(double)(threshold);
                p2errors = p2errors | (double)(Math.Abs(vy-vy2))>(double)(threshold);
                p2errors = p2errors | (double)(Math.Abs(vdx-vdx2))>(double)(threshold);
                p2errors = p2errors | (double)(Math.Abs(vdy-vdy2))>(double)(threshold);
                pspline.pspline2diff2(ref p2, v0, ref vx, ref vdx, ref vd2x, ref vy, ref vdy, ref vd2y);
                p2errors = p2errors | (double)(Math.Abs(vx-vx2))>(double)(threshold);
                p2errors = p2errors | (double)(Math.Abs(vy-vy2))>(double)(threshold);
                p2errors = p2errors | (double)(Math.Abs(vdx-vdx2))>(double)(threshold);
                p2errors = p2errors | (double)(Math.Abs(vdy-vdy2))>(double)(threshold);
                p2errors = p2errors | (double)(Math.Abs(vd2x-vd2x2))>(double)(threshold);
                p2errors = p2errors | (double)(Math.Abs(vd2y-vd2y2))>(double)(threshold);
                
                //
                // 3D test
                //
                pspline.pspline3calc(ref p3, v0, ref vx, ref vy, ref vz);
                p3errors = p3errors | (double)(Math.Abs(vx-vx2))>(double)(threshold);
                p3errors = p3errors | (double)(Math.Abs(vy-vy2))>(double)(threshold);
                p3errors = p3errors | (double)(Math.Abs(vz-vz2))>(double)(threshold);
                pspline.pspline3diff(ref p3, v0, ref vx, ref vdx, ref vy, ref vdy, ref vz, ref vdz);
                p3errors = p3errors | (double)(Math.Abs(vx-vx2))>(double)(threshold);
                p3errors = p3errors | (double)(Math.Abs(vy-vy2))>(double)(threshold);
                p3errors = p3errors | (double)(Math.Abs(vz-vz2))>(double)(threshold);
                p3errors = p3errors | (double)(Math.Abs(vdx-vdx2))>(double)(threshold);
                p3errors = p3errors | (double)(Math.Abs(vdy-vdy2))>(double)(threshold);
                p3errors = p3errors | (double)(Math.Abs(vdz-vdz2))>(double)(threshold);
                pspline.pspline3diff2(ref p3, v0, ref vx, ref vdx, ref vd2x, ref vy, ref vdy, ref vd2y, ref vz, ref vdz, ref vd2z);
                p3errors = p3errors | (double)(Math.Abs(vx-vx2))>(double)(threshold);
                p3errors = p3errors | (double)(Math.Abs(vy-vy2))>(double)(threshold);
                p3errors = p3errors | (double)(Math.Abs(vz-vz2))>(double)(threshold);
                p3errors = p3errors | (double)(Math.Abs(vdx-vdx2))>(double)(threshold);
                p3errors = p3errors | (double)(Math.Abs(vdy-vdy2))>(double)(threshold);
                p3errors = p3errors | (double)(Math.Abs(vdz-vdz2))>(double)(threshold);
                p3errors = p3errors | (double)(Math.Abs(vd2x-vd2x2))>(double)(threshold);
                p3errors = p3errors | (double)(Math.Abs(vd2y-vd2y2))>(double)(threshold);
                p3errors = p3errors | (double)(Math.Abs(vd2z-vd2z2))>(double)(threshold);
                
                //
                // Test tangents for 2D/3D
                //
                pspline.pspline2tangent(ref p2, v0, ref vx, ref vy);
                p2errors = p2errors | (double)(Math.Abs(vx-vdx2/apserv.safepythag2(vdx2, vdy2)))>(double)(threshold);
                p2errors = p2errors | (double)(Math.Abs(vy-vdy2/apserv.safepythag2(vdx2, vdy2)))>(double)(threshold);
                pspline.pspline3tangent(ref p3, v0, ref vx, ref vy, ref vz);
                p3errors = p3errors | (double)(Math.Abs(vx-vdx2/apserv.safepythag3(vdx2, vdy2, vdz2)))>(double)(threshold);
                p3errors = p3errors | (double)(Math.Abs(vy-vdy2/apserv.safepythag3(vdx2, vdy2, vdz2)))>(double)(threshold);
                p3errors = p3errors | (double)(Math.Abs(vz-vdz2/apserv.safepythag3(vdx2, vdy2, vdz2)))>(double)(threshold);
            }
            
            //
            // Arc length test.
            //
            // Simple problem with easy solution (points on a straight line with
            // uniform parameterization).
            //
            for(n=2; n<=maxn; n++)
            {
                xy = new double[n, 2];
                xyz = new double[n, 3];
                for(i=0; i<=n-1; i++)
                {
                    xy[i,0] = i;
                    xy[i,1] = i;
                    xyz[i,0] = i;
                    xyz[i,1] = i;
                    xyz[i,2] = i;
                }
                pspline.pspline2build(xy, n, 1, 0, ref p2);
                pspline.pspline3build(xyz, n, 1, 0, ref p3);
                a = AP.Math.RandomReal();
                b = AP.Math.RandomReal();
                p2errors = p2errors | (double)(Math.Abs(pspline.pspline2arclength(ref p2, a, b)-(b-a)*Math.Sqrt(2)*(n-1)))>(double)(nonstrictthreshold);
                p3errors = p3errors | (double)(Math.Abs(pspline.pspline3arclength(ref p3, a, b)-(b-a)*Math.Sqrt(3)*(n-1)))>(double)(nonstrictthreshold);
            }
            
            //
            // report
            //
            waserrors = p2errors | p3errors;
            if( !silent )
            {
                System.Console.Write("TESTING SPLINE INTERPOLATION");
                System.Console.WriteLine();
                
                //
                // Normal tests
                //
                System.Console.Write("2D TEST:                                 ");
                if( p2errors )
                {
                    System.Console.Write("FAILED");
                    System.Console.WriteLine();
                }
                else
                {
                    System.Console.Write("OK");
                    System.Console.WriteLine();
                }
                System.Console.Write("3D TEST:                                 ");
                if( p3errors )
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
            
            //
            // end
            //
            result = !waserrors;
            return result;
        }


        /*************************************************************************
        Unset spline, i.e. initialize it with random garbage
        *************************************************************************/
        private static void unsetp2(ref pspline.pspline2interpolant p)
        {
            double[,] xy = new double[0,0];

            xy = new double[2, 2];
            xy[0,0] = -1;
            xy[0,1] = -1;
            xy[1,0] = +1;
            xy[1,1] = +1;
            pspline.pspline2build(xy, 2, 1, 0, ref p);
        }


        /*************************************************************************
        Unset spline, i.e. initialize it with random garbage
        *************************************************************************/
        private static void unsetp3(ref pspline.pspline3interpolant p)
        {
            double[,] xy = new double[0,0];

            xy = new double[2, 3];
            xy[0,0] = -1;
            xy[0,1] = -1;
            xy[0,2] = -1;
            xy[1,0] = +1;
            xy[1,1] = +1;
            xy[1,2] = +1;
            pspline.pspline3build(xy, 2, 1, 0, ref p);
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
        Silent unit test
        *************************************************************************/
        public static bool testpsplineunit_test_silent()
        {
            bool result = new bool();

            result = testpsplineinterpolation(true);
            return result;
        }


        /*************************************************************************
        Unit test
        *************************************************************************/
        public static bool testpsplineunit_test()
        {
            bool result = new bool();

            result = testpsplineinterpolation(false);
            return result;
        }
    }
}
