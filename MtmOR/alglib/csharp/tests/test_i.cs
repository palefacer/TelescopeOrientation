
#pragma warning disable 219
#pragma warning disable 162
using System;
public class MainTest
{

    public static bool doc_test_bool(bool v, bool t)
    { return (v && t) || (!v && !t); }

    public static bool doc_test_int(int v, int t)
    { return v==t; }

    public static bool doc_test_real(double v, double t, double threshold)
    { return Math.Abs(v-t)<=threshold; }

    public static bool doc_test_complex(alglib.complex v, alglib.complex t, double threshold)
    { return alglib.math.abscomplex(v-t)<=threshold; }


    public static bool doc_test_bool_vector(bool[] v, bool[] t)
    {
        int i;
        if( alglib.ap.len(v)!=alglib.ap.len(t) )
            return false;
        for(i=0; i<alglib.ap.len(v); i++)
            if( v[i]!=t[i] )
                return false;
        return true;
    }

    public static bool doc_test_bool_matrix(bool[,] v, bool[,] t)
    {
        int i, j;
        if( alglib.ap.rows(v)!=alglib.ap.rows(t) )
            return false;
        if( alglib.ap.cols(v)!=alglib.ap.cols(t) )
            return false;
        for(i=0; i<alglib.ap.rows(v); i++)
            for(j=0; j<alglib.ap.cols(v); j++)
                if( v[i,j]!=t[i,j] )
                    return false;
        return true;
    }

    public static bool doc_test_int_vector(int[] v, int[] t)
    {
        int i;
        if( alglib.ap.len(v)!=alglib.ap.len(t) )
            return false;
        for(i=0; i<alglib.ap.len(v); i++)
            if( v[i]!=t[i] )
                return false;
        return true;
    }

    public static bool doc_test_int_matrix(int[,] v, int[,] t)
    {
        int i, j;
        if( alglib.ap.rows(v)!=alglib.ap.rows(t) )
            return false;
        if( alglib.ap.cols(v)!=alglib.ap.cols(t) )
            return false;
        for(i=0; i<alglib.ap.rows(v); i++)
            for(j=0; j<alglib.ap.cols(v); j++)
                if( v[i,j]!=t[i,j] )
                    return false;
        return true;
    }

    public static bool doc_test_real_vector(double[] v, double[] t, double threshold)
    {
        int i;
        if( alglib.ap.len(v)!=alglib.ap.len(t) )
            return false;
        for(i=0; i<alglib.ap.len(v); i++)
            if( Math.Abs(v[i]-t[i])>threshold )
                return false;
        return true;
    }

    public static bool doc_test_real_matrix(double[,] v, double[,] t, double threshold)
    {
        int i, j;
        if( alglib.ap.rows(v)!=alglib.ap.rows(t) )
            return false;
        if( alglib.ap.cols(v)!=alglib.ap.cols(t) )
            return false;
        for(i=0; i<alglib.ap.rows(v); i++)
            for(j=0; j<alglib.ap.cols(v); j++)
                if( Math.Abs(v[i,j]-t[i,j])>threshold )
                    return false;
        return true;
    }

    public static bool doc_test_complex_vector(alglib.complex[] v, alglib.complex[] t, double threshold)
    {
        int i;
        if( alglib.ap.len(v)!=alglib.ap.len(t) )
            return false;
        for(i=0; i<alglib.ap.len(v); i++)
            if( alglib.math.abscomplex(v[i]-t[i])>threshold )
                return false;
        return true;
    }

    public static bool doc_test_complex_matrix(alglib.complex[,] v, alglib.complex[,] t, double threshold)
    {
        int i, j;
        if( alglib.ap.rows(v)!=alglib.ap.rows(t) )
            return false;
        if( alglib.ap.cols(v)!=alglib.ap.cols(t) )
            return false;
        for(i=0; i<alglib.ap.rows(v); i++)
            for(j=0; j<alglib.ap.cols(v); j++)
                if( alglib.math.abscomplex(v[i,j]-t[i,j])>threshold )
                    return false;
        return true;
    }

    public static void spoil_vector_by_adding_element<T>(ref T[] x, T val) where T : new()
    {
        int i;
        T[] y = x;
        x = new T[y.Length+1];
        for(i=0; i<y.Length; i++)
            x[i] = y[i];
        x[y.Length] = val;
    }

    public static void spoil_vector_by_deleting_element<T>(ref T[] x) where T : new()
    {
        int i;
        T[] y = x;
        x = new T[y.Length-1];
        for(i=0; i<y.Length-1; i++)
            x[i] = y[i];
    }

    public static void spoil_matrix_by_adding_row<T>(ref T[,] x, T val) where T : new()
    {
        int i, j;
        T[,] y = x;
        x = new T[y.GetLength(0)+1,y.GetLength(1)];
        for(i=0; i<y.GetLength(0); i++)
            for(j=0; j<y.GetLength(1); j++)
                x[i,j] = y[i,j];
        for(j=0; j<y.GetLength(1); j++)
            x[y.GetLength(0),j] = val;
    }

    public static void spoil_matrix_by_deleting_row<T>(ref T[,] x) where T : new()
    {
        int i, j;
        T[,] y = x;
        x = new T[y.GetLength(0)-1,y.GetLength(1)];
        for(i=0; i<y.GetLength(0)-1; i++)
            for(j=0; j<y.GetLength(1); j++)
                x[i,j] = y[i,j];
    }

    public static void spoil_matrix_by_adding_col<T>(ref T[,] x, T val) where T : new()
    {
        int i, j;
        T[,] y = x;
        x = new T[y.GetLength(0), y.GetLength(1)+1];
        for(i=0; i<y.GetLength(0); i++)
            for(j=0; j<y.GetLength(1); j++)
                x[i,j] = y[i,j];
        for(i=0; i<y.GetLength(0); i++)
            x[i,y.GetLength(1)] = val;
    }

    public static void spoil_matrix_by_deleting_col<T>(ref T[,] x) where T : new()
    {
        int i, j;
        T[,] y = x;
        x = new T[y.GetLength(0), y.GetLength(1)-1];
        for(i=0; i<y.GetLength(0); i++)
            for(j=0; j<y.GetLength(1)-1; j++)
                x[i,j] = y[i,j];
    }

    public static void spoil_vector_by_value<T>(ref T[] x, T val)
    {
        if( x.Length!=0 )
            x[alglib.math.randominteger(x.Length)] = val;
    }
    public static void spoil_matrix_by_value<T>(ref T[,] x, T val)
    {
        if( x.GetLength(0)!=0 && x.GetLength(1)!=0 )
            x[alglib.math.randominteger(x.GetLength(0)),alglib.math.randominteger(x.GetLength(1))] = val;
    }

    public static void function1_func(double[] x, ref double func, object obj)
    {
        // this callback calculates f(x0,x1) = 100*(x0+3)^4 + (x1-3)^4
        func = 100*System.Math.Pow(x[0]+3,4) + System.Math.Pow(x[1]-3,4);
    }
    public static void function1_grad(double[] x, ref double func, double[] grad, object obj)
    {
        // this callback calculates f(x0,x1) = 100*(x0+3)^4 + (x1-3)^4
        // and its derivatives df/d0 and df/dx1
        func = 100*System.Math.Pow(x[0]+3,4) + System.Math.Pow(x[1]-3,4);
        grad[0] = 400*System.Math.Pow(x[0]+3,3);
        grad[1] = 4*System.Math.Pow(x[1]-3,3);
    }
    public static void function1_hess(double[] x, ref double func, double[] grad, double[,] hess, object obj)
    {
        // this callback calculates f(x0,x1) = 100*(x0+3)^4 + (x1-3)^4
        // its derivatives df/d0 and df/dx1
        // and its Hessian.
        func = 100*System.Math.Pow(x[0]+3,4) + System.Math.Pow(x[1]-3,4);
        grad[0] = 400*System.Math.Pow(x[0]+3,3);
        grad[1] = 4*System.Math.Pow(x[1]-3,3);
        hess[0,0] = 1200*System.Math.Pow(x[0]+3,2);
        hess[0,1] = 0;
        hess[1,0] = 0;
        hess[1,1] = 12*System.Math.Pow(x[1]-3,2);
    }
    public static void  function1_jac(double[] x, double[] fi, double[,] jac, object obj)
    {
        // this callback calculates
        // f0(x0,x1) = 100*(x0+3)^4,
        // f1(x0,x1) = (x1-3)^4
        // and Jacobian matrix J = [dfi/dxj]
        fi[0] = 10*System.Math.Pow(x[0]+3,2);
        fi[1] = System.Math.Pow(x[1]-3,2);
        jac[0,0] = 20*(x[0]+3);
        jac[0,1] = 0;
        jac[1,0] = 0;
        jac[1,1] = 2*(x[1]-3);
    }
    public static void function_cx_1_func(double[] c, double[] x, ref double func, object obj)
    {
        // this callback calculates f(c,x)=exp(-c0*sqr(x0))
        // where x is a position on X-axis and c is adjustable parameter
        func = System.Math.Exp(-c[0]*x[0]*x[0]);
    }
    public static void function_cx_1_grad(double[] c, double[] x, ref double func, double[] grad, object obj)
    {
        // this callback calculates f(c,x)=exp(-c0*sqr(x0)) and gradient G={df/dc[i]}
        // where x is a position on X-axis and c is adjustable parameter.
        // IMPORTANT: gradient is calculated with respect to C, not to X
        func = System.Math.Exp(-c[0]*System.Math.Pow(x[0],2));
        grad[0] = -System.Math.Pow(x[0],2)*func;
    }
    public static void function_cx_1_hess(double[] c, double[] x, ref double func, double[] grad, double[,] hess, object obj)
    {
        // this callback calculates f(c,x)=exp(-c0*sqr(x0)), gradient G={df/dc[i]} and Hessian H={d2f/(dc[i]*dc[j])}
        // where x is a position on X-axis and c is adjustable parameter.
        // IMPORTANT: gradient/Hessian are calculated with respect to C, not to X
        func = System.Math.Exp(-c[0]*System.Math.Pow(x[0],2));
        grad[0] = -System.Math.Pow(x[0],2)*func;
        hess[0,0] = System.Math.Pow(x[0],4)*func;
    }
    public static void ode_function_1_diff(double[] y, double x, double[] dy, object obj)
    {
        // this callback calculates f(y[],x)=-y[0]
        dy[0] = -y[0];
    }
    public static void int_function_1_func(double x, double xminusa, double bminusx, ref double y, object obj)
    {
        // this callback calculates f(x)=exp(x)
        y = Math.Exp(x);
    }

    public static int Main(string[] args)
    {
        bool _TotalResult = true;
        bool _TestResult;
        int _spoil_scenario;
        System.Console.WriteLine("C# interface tests. Please wait...");
        try
        {
            //
            // TEST matinv_d_r1
            //      Real matrix inverse
            //
            System.Console.WriteLine("0/66");
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<7; _spoil_scenario++)
            {
                try
                {
                    double[,] a = new double[,]{{1,-1},{1,1}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref a, Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref a, Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref a, Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_matrix_by_adding_row(ref a, 0);
                    if( _spoil_scenario==4 )
                        spoil_matrix_by_adding_col(ref a, 0);
                    if( _spoil_scenario==5 )
                        spoil_matrix_by_deleting_row(ref a);
                    if( _spoil_scenario==6 )
                        spoil_matrix_by_deleting_col(ref a);
                    int info;
                    alglib.matinvreport rep;
                    alglib.rmatrixinverse(ref a, out info, out rep);
                    _TestResult = _TestResult && doc_test_int(info, 1);
                    _TestResult = _TestResult && doc_test_real_matrix(a, new double[,]{{0.5,0.5},{-0.5,0.5}}, 0.00005);
                    _TestResult = _TestResult && doc_test_real(rep.r1, 0.5, 0.00005);
                    _TestResult = _TestResult && doc_test_real(rep.rinf, 0.5, 0.00005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "matinv_d_r1");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST matinv_d_c1
            //      Complex matrix inverse
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<7; _spoil_scenario++)
            {
                try
                {
                    alglib.complex[,] a = new alglib.complex[,]{{new alglib.complex(0,1),-1},{new alglib.complex(0,1),1}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref a, Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref a, Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref a, Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_matrix_by_adding_row(ref a, 0);
                    if( _spoil_scenario==4 )
                        spoil_matrix_by_adding_col(ref a, 0);
                    if( _spoil_scenario==5 )
                        spoil_matrix_by_deleting_row(ref a);
                    if( _spoil_scenario==6 )
                        spoil_matrix_by_deleting_col(ref a);
                    int info;
                    alglib.matinvreport rep;
                    alglib.cmatrixinverse(ref a, out info, out rep);
                    _TestResult = _TestResult && doc_test_int(info, 1);
                    _TestResult = _TestResult && doc_test_complex_matrix(a, new alglib.complex[,]{{new alglib.complex(0,-0.5),new alglib.complex(0,-0.5)},{-0.5,0.5}}, 0.00005);
                    _TestResult = _TestResult && doc_test_real(rep.r1, 0.5, 0.00005);
                    _TestResult = _TestResult && doc_test_real(rep.rinf, 0.5, 0.00005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "matinv_d_c1");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST matinv_d_spd1
            //      SPD matrix inverse
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<7; _spoil_scenario++)
            {
                try
                {
                    double[,] a = new double[,]{{2,1},{1,2}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref a, Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref a, Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref a, Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_matrix_by_adding_row(ref a, 0);
                    if( _spoil_scenario==4 )
                        spoil_matrix_by_adding_col(ref a, 0);
                    if( _spoil_scenario==5 )
                        spoil_matrix_by_deleting_row(ref a);
                    if( _spoil_scenario==6 )
                        spoil_matrix_by_deleting_col(ref a);
                    int info;
                    alglib.matinvreport rep;
                    alglib.spdmatrixinverse(ref a, out info, out rep);
                    _TestResult = _TestResult && doc_test_int(info, 1);
                    _TestResult = _TestResult && doc_test_real_matrix(a, new double[,]{{0.666666,-0.333333},{-0.333333,0.666666}}, 0.00005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "matinv_d_spd1");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST matinv_d_hpd1
            //      HPD matrix inverse
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<7; _spoil_scenario++)
            {
                try
                {
                    alglib.complex[,] a = new alglib.complex[,]{{2,1},{1,2}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref a, Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref a, Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref a, Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_matrix_by_adding_row(ref a, 0);
                    if( _spoil_scenario==4 )
                        spoil_matrix_by_adding_col(ref a, 0);
                    if( _spoil_scenario==5 )
                        spoil_matrix_by_deleting_row(ref a);
                    if( _spoil_scenario==6 )
                        spoil_matrix_by_deleting_col(ref a);
                    int info;
                    alglib.matinvreport rep;
                    alglib.hpdmatrixinverse(ref a, out info, out rep);
                    _TestResult = _TestResult && doc_test_int(info, 1);
                    _TestResult = _TestResult && doc_test_complex_matrix(a, new alglib.complex[,]{{0.666666,-0.333333},{-0.333333,0.666666}}, 0.00005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "matinv_d_hpd1");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST matinv_t_r1
            //      Real matrix inverse: singular matrix
            //
            _TestResult = true;
            try
            {
                double[,] a = new double[,]{{1,-1},{-2,2}};
                int info;
                alglib.matinvreport rep;
                alglib.rmatrixinverse(ref a, out info, out rep);
                _TestResult = _TestResult && doc_test_int(info, -3);
                _TestResult = _TestResult && doc_test_real(rep.r1, 0.0, 0.00005);
                _TestResult = _TestResult && doc_test_real(rep.rinf, 0.0, 0.00005);
            }
            catch(alglib.alglibexception)
            { _TestResult = false; }
            catch
            { throw; }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "matinv_t_r1");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST matinv_t_c1
            //      Complex matrix inverse: singular matrix
            //
            _TestResult = true;
            try
            {
                alglib.complex[,] a = new alglib.complex[,]{{new alglib.complex(0,1),new alglib.complex(0,-1)},{-2,2}};
                int info;
                alglib.matinvreport rep;
                alglib.cmatrixinverse(ref a, out info, out rep);
                _TestResult = _TestResult && doc_test_int(info, -3);
                _TestResult = _TestResult && doc_test_real(rep.r1, 0.0, 0.00005);
                _TestResult = _TestResult && doc_test_real(rep.rinf, 0.0, 0.00005);
            }
            catch(alglib.alglibexception)
            { _TestResult = false; }
            catch
            { throw; }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "matinv_t_c1");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST matinv_e_spd1
            //      Attempt to use SPD function on nonsymmetrix matrix
            //
            _TestResult = true;
            try
            {
                double[,] a = new double[,]{{1,0},{1,1}};
                int info;
                alglib.matinvreport rep;
                alglib.spdmatrixinverse(ref a, out info, out rep);
                _TestResult = false;
            }
            catch(alglib.alglibexception)
            {}
            catch
            { throw; }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "matinv_e_spd1");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST matinv_e_hpd1
            //      Attempt to use SPD function on nonsymmetrix matrix
            //
            _TestResult = true;
            try
            {
                alglib.complex[,] a = new alglib.complex[,]{{1,0},{1,1}};
                int info;
                alglib.matinvreport rep;
                alglib.hpdmatrixinverse(ref a, out info, out rep);
                _TestResult = false;
            }
            catch(alglib.alglibexception)
            {}
            catch
            { throw; }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "matinv_e_hpd1");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST minlbfgs_d_1
            //      Nonlinear optimization, F(x,y) = 100*(x+3)^4+(y-3)^4
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<12; _spoil_scenario++)
            {
                try
                {
                    //
                    // This example demonstrates minimization of f(x,y) = 100*(x+3)^4+(y-3)^4
                    // using LBFGS method.
                    //
                    double[] x = new double[]{0,0};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, Double.NegativeInfinity);
                    double epsg = 0.0000000001;
                    if( _spoil_scenario==3 )
                        epsg = Double.NaN;
                    if( _spoil_scenario==4 )
                        epsg = Double.PositiveInfinity;
                    if( _spoil_scenario==5 )
                        epsg = Double.NegativeInfinity;
                    double epsf = 0;
                    if( _spoil_scenario==6 )
                        epsf = Double.NaN;
                    if( _spoil_scenario==7 )
                        epsf = Double.PositiveInfinity;
                    if( _spoil_scenario==8 )
                        epsf = Double.NegativeInfinity;
                    double epsx = 0;
                    if( _spoil_scenario==9 )
                        epsx = Double.NaN;
                    if( _spoil_scenario==10 )
                        epsx = Double.PositiveInfinity;
                    if( _spoil_scenario==11 )
                        epsx = Double.NegativeInfinity;
                    int maxits = 0;
                    alglib.minlbfgsstate state;
                    alglib.minlbfgsreport rep;

                    alglib.minlbfgscreate(1, x, out state);
                    alglib.minlbfgssetcond(state, epsg, epsf, epsx, maxits);
                    alglib.minlbfgsoptimize(state, function1_grad, null, null);
                    alglib.minlbfgsresults(state, out x, out rep);

                    _TestResult = _TestResult && doc_test_int(rep.terminationtype, 4);
                    _TestResult = _TestResult && doc_test_real_vector(x, new double[]{-3,3}, 0.005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "minlbfgs_d_1");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST minlbfgs_d_2
            //      Nonlinear optimization with additional settings and restarts, F(x,y) = 100*(x+3)^4+(y-3)^4
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<18; _spoil_scenario++)
            {
                try
                {
                    //
                    // This example demonstrates minimization of f(x,y) = 100*(x+3)^4+(y-3)^4
                    // using LBFGS method.
                    //
                    // Several advanced techniques are demonstrated:
                    // * upper limit on step size
                    // * restart from new point
                    //
                    double[] x = new double[]{0,0};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, Double.NegativeInfinity);
                    double epsg = 0.0000000001;
                    if( _spoil_scenario==3 )
                        epsg = Double.NaN;
                    if( _spoil_scenario==4 )
                        epsg = Double.PositiveInfinity;
                    if( _spoil_scenario==5 )
                        epsg = Double.NegativeInfinity;
                    double epsf = 0;
                    if( _spoil_scenario==6 )
                        epsf = Double.NaN;
                    if( _spoil_scenario==7 )
                        epsf = Double.PositiveInfinity;
                    if( _spoil_scenario==8 )
                        epsf = Double.NegativeInfinity;
                    double epsx = 0;
                    if( _spoil_scenario==9 )
                        epsx = Double.NaN;
                    if( _spoil_scenario==10 )
                        epsx = Double.PositiveInfinity;
                    if( _spoil_scenario==11 )
                        epsx = Double.NegativeInfinity;
                    double stpmax = 0.1;
                    if( _spoil_scenario==12 )
                        stpmax = Double.NaN;
                    if( _spoil_scenario==13 )
                        stpmax = Double.PositiveInfinity;
                    if( _spoil_scenario==14 )
                        stpmax = Double.NegativeInfinity;
                    int maxits = 0;
                    alglib.minlbfgsstate state;
                    alglib.minlbfgsreport rep;

                    // first run
                    alglib.minlbfgscreate(1, x, out state);
                    alglib.minlbfgssetcond(state, epsg, epsf, epsx, maxits);
                    alglib.minlbfgssetstpmax(state, stpmax);
                    alglib.minlbfgsoptimize(state, function1_grad, null, null);
                    alglib.minlbfgsresults(state, out x, out rep);

                    _TestResult = _TestResult && doc_test_real_vector(x, new double[]{-3,3}, 0.005);

                    // second run - algorithm is restarted
                    x = new double[]{10,10};
                    if( _spoil_scenario==15 )
                        spoil_vector_by_value(ref x, Double.NaN);
                    if( _spoil_scenario==16 )
                        spoil_vector_by_value(ref x, Double.PositiveInfinity);
                    if( _spoil_scenario==17 )
                        spoil_vector_by_value(ref x, Double.NegativeInfinity);
                    alglib.minlbfgsrestartfrom(state, x);
                    alglib.minlbfgsoptimize(state, function1_grad, null, null);
                    alglib.minlbfgsresults(state, out x, out rep);

                    _TestResult = _TestResult && doc_test_int(rep.terminationtype, 4);
                    _TestResult = _TestResult && doc_test_real_vector(x, new double[]{-3,3}, 0.005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "minlbfgs_d_2");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST odesolver_d1
            //      Solving y'=-y with ODE solver
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<13; _spoil_scenario++)
            {
                try
                {
                    double[] y = new double[]{1};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref y, Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref y, Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref y, Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_vector_by_deleting_element(ref y);
                    double[] x = new double[]{0,1,2,3};
                    if( _spoil_scenario==4 )
                        spoil_vector_by_value(ref x, Double.NaN);
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref x, Double.PositiveInfinity);
                    if( _spoil_scenario==6 )
                        spoil_vector_by_value(ref x, Double.NegativeInfinity);
                    double eps = 0.00001;
                    if( _spoil_scenario==7 )
                        eps = Double.NaN;
                    if( _spoil_scenario==8 )
                        eps = Double.PositiveInfinity;
                    if( _spoil_scenario==9 )
                        eps = Double.NegativeInfinity;
                    double h = 0;
                    if( _spoil_scenario==10 )
                        h = Double.NaN;
                    if( _spoil_scenario==11 )
                        h = Double.PositiveInfinity;
                    if( _spoil_scenario==12 )
                        h = Double.NegativeInfinity;
                    alglib.odesolverstate s;
                    int m;
                    double[] xtbl;
                    double[,] ytbl;
                    alglib.odesolverreport rep;
                    alglib.odesolverrkck(y, x, eps, h, out s);
                    alglib.odesolversolve(s, ode_function_1_diff, null);
                    alglib.odesolverresults(s, out m, out xtbl, out ytbl, out rep);
                    _TestResult = _TestResult && doc_test_int(m, 4);
                    _TestResult = _TestResult && doc_test_real_vector(xtbl, new double[]{0,1,2,3}, 0.005);
                    _TestResult = _TestResult && doc_test_real_matrix(ytbl, new double[,]{{1},{0.367},{0.135},{0.050}}, 0.005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "odesolver_d1");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST fft_complex_d1
            //      Complex FFT: simple example
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
            {
                try
                {
                    //
                    // first we demonstrate forward FFT:
                    // [1i,1i,1i,1i] is converted to [4i, 0, 0, 0]
                    //
                    alglib.complex[] z = new alglib.complex[]{new alglib.complex(0,1),new alglib.complex(0,1),new alglib.complex(0,1),new alglib.complex(0,1)};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref z, Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref z, Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref z, Double.NegativeInfinity);
                    alglib.fftc1d(ref z);
                    _TestResult = _TestResult && doc_test_complex_vector(z, new alglib.complex[]{new alglib.complex(0,4),0,0,0}, 0.0001);

                    //
                    // now we convert [4i, 0, 0, 0] back to [1i,1i,1i,1i]
                    // with backward FFT
                    //
                    alglib.fftc1dinv(ref z);
                    _TestResult = _TestResult && doc_test_complex_vector(z, new alglib.complex[]{new alglib.complex(0,1),new alglib.complex(0,1),new alglib.complex(0,1),new alglib.complex(0,1)}, 0.0001);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "fft_complex_d1");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST fft_complex_d2
            //      Complex FFT: advanced example
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
            {
                try
                {
                    //
                    // first we demonstrate forward FFT:
                    // [0,1,0,1i] is converted to [1+1i, -1-1i, -1-1i, 1+1i]
                    //
                    alglib.complex[] z = new alglib.complex[]{0,1,0,new alglib.complex(0,1)};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref z, Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref z, Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref z, Double.NegativeInfinity);
                    alglib.fftc1d(ref z);
                    _TestResult = _TestResult && doc_test_complex_vector(z, new alglib.complex[]{new alglib.complex(1,+1),new alglib.complex(-1,-1),new alglib.complex(-1,-1),new alglib.complex(1,+1)}, 0.0001);

                    //
                    // now we convert result back with backward FFT
                    //
                    alglib.fftc1dinv(ref z);
                    _TestResult = _TestResult && doc_test_complex_vector(z, new alglib.complex[]{0,1,0,new alglib.complex(0,1)}, 0.0001);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "fft_complex_d2");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST fft_real_d1
            //      Real FFT: simple example
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
            {
                try
                {
                    //
                    // first we demonstrate forward FFT:
                    // [1,1,1,1] is converted to [4, 0, 0, 0]
                    //
                    double[] x = new double[]{1,1,1,1};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, Double.NegativeInfinity);
                    alglib.complex[] f;
                    double[] x2;
                    alglib.fftr1d(x, out f);
                    _TestResult = _TestResult && doc_test_complex_vector(f, new alglib.complex[]{4,0,0,0}, 0.0001);

                    //
                    // now we convert [4, 0, 0, 0] back to [1,1,1,1]
                    // with backward FFT
                    //
                    alglib.fftr1dinv(f, out x2);
                    _TestResult = _TestResult && doc_test_real_vector(x2, new double[]{1,1,1,1}, 0.0001);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "fft_real_d1");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST fft_real_d2
            //      Real FFT: advanced example
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
            {
                try
                {
                    //
                    // first we demonstrate forward FFT:
                    // [1,2,3,4] is converted to [10, -2+2i, -2, -2-2i]
                    //
                    // note that output array is self-adjoint:
                    // * f[0] = conj(f[0])
                    // * f[1] = conj(f[3])
                    // * f[2] = conj(f[2])
                    //
                    double[] x = new double[]{1,2,3,4};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, Double.NegativeInfinity);
                    alglib.complex[] f;
                    double[] x2;
                    alglib.fftr1d(x, out f);
                    _TestResult = _TestResult && doc_test_complex_vector(f, new alglib.complex[]{10,new alglib.complex(-2,+2),-2,new alglib.complex(-2,-2)}, 0.0001);

                    //
                    // now we convert [10, -2+2i, -2, -2-2i] back to [1,2,3,4]
                    //
                    alglib.fftr1dinv(f, out x2);
                    _TestResult = _TestResult && doc_test_real_vector(x2, new double[]{1,2,3,4}, 0.0001);

                    //
                    // remember that F is self-adjoint? It means that we can pass just half
                    // (slightly larger than half) of F to inverse real FFT and still get our result.
                    //
                    // I.e. instead [10, -2+2i, -2, -2-2i] we pass just [10, -2+2i, -2] and everything works!
                    //
                    // NOTE: in this case we should explicitly pass array length (which is 4) to ALGLIB;
                    // if not, it will automatically use array length to determine FFT size and
                    // will erroneously make half-length FFT.
                    //
                    f = new alglib.complex[]{10,new alglib.complex(-2,+2),-2};
                    alglib.fftr1dinv(f, 4, out x2);
                    _TestResult = _TestResult && doc_test_real_vector(x2, new double[]{1,2,3,4}, 0.0001);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "fft_real_d2");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST fft_complex_e1
            //      error detection in backward FFT
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
            {
                try
                {
                    alglib.complex[] z = new alglib.complex[]{0,2,0,-2};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref z, Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref z, Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref z, Double.NegativeInfinity);
                    alglib.fftc1dinv(ref z);
                    _TestResult = _TestResult && doc_test_complex_vector(z, new alglib.complex[]{0,new alglib.complex(0,1),0,new alglib.complex(0,-1)}, 0.0001);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "fft_complex_e1");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST autogk_d1
            //      Integrating f=exp(x) by adaptive integrator
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<6; _spoil_scenario++)
            {
                try
                {
                    //
                    // This example demonstrates integration of f=exp(x) on [0,1]:
                    // * first, autogkstate is initialized
                    // * then we call integration function
                    // * and finally we obtain results with autogkresults() call
                    //
                    double a = 0;
                    if( _spoil_scenario==0 )
                        a = Double.NaN;
                    if( _spoil_scenario==1 )
                        a = Double.PositiveInfinity;
                    if( _spoil_scenario==2 )
                        a = Double.NegativeInfinity;
                    double b = 1;
                    if( _spoil_scenario==3 )
                        b = Double.NaN;
                    if( _spoil_scenario==4 )
                        b = Double.PositiveInfinity;
                    if( _spoil_scenario==5 )
                        b = Double.NegativeInfinity;
                    alglib.autogkstate s;
                    double v;
                    alglib.autogkreport rep;

                    alglib.autogksmooth(a, b, out s);
                    alglib.autogkintegrate(s, int_function_1_func, null);
                    alglib.autogkresults(s, out v, out rep);

                    _TestResult = _TestResult && doc_test_real(v, 1.7182, 0.005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "autogk_d1");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST lsfit_d_1
            //      Nonlinear fitting by f(x,c) = exp(-c*x^2) (FG scheme)
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<19; _spoil_scenario++)
            {
                try
                {
                    double[,] x = new double[,]{{-1},{-0.8},{-0.6},{-0.4},{-0.2},{0},{0.2},{0.4},{0.6},{0.8},{1.0}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref x, Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref x, Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref x, Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_matrix_by_deleting_row(ref x);
                    if( _spoil_scenario==4 )
                        spoil_matrix_by_deleting_col(ref x);
                    double[] y = new double[]{0.223130,0.382893,0.582748,0.786628,0.941765,1.000000,0.941765,0.786628,0.582748,0.382893,0.223130};
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref y, Double.NaN);
                    if( _spoil_scenario==6 )
                        spoil_vector_by_value(ref y, Double.PositiveInfinity);
                    if( _spoil_scenario==7 )
                        spoil_vector_by_value(ref y, Double.NegativeInfinity);
                    if( _spoil_scenario==8 )
                        spoil_vector_by_adding_element(ref y, 0);
                    if( _spoil_scenario==9 )
                        spoil_vector_by_deleting_element(ref y);
                    double[] c = new double[]{0.3};
                    if( _spoil_scenario==10 )
                        spoil_vector_by_value(ref c, Double.NaN);
                    if( _spoil_scenario==11 )
                        spoil_vector_by_value(ref c, Double.PositiveInfinity);
                    if( _spoil_scenario==12 )
                        spoil_vector_by_value(ref c, Double.NegativeInfinity);
                    double epsf = 0;
                    if( _spoil_scenario==13 )
                        epsf = Double.NaN;
                    if( _spoil_scenario==14 )
                        epsf = Double.PositiveInfinity;
                    if( _spoil_scenario==15 )
                        epsf = Double.NegativeInfinity;
                    double epsx = 0.000001;
                    if( _spoil_scenario==16 )
                        epsx = Double.NaN;
                    if( _spoil_scenario==17 )
                        epsx = Double.PositiveInfinity;
                    if( _spoil_scenario==18 )
                        epsx = Double.NegativeInfinity;
                    int maxits = 0;
                    int info;
                    alglib.lsfitstate state;
                    alglib.lsfitreport rep;
                    alglib.lsfitcreatefg(x, y, c, true, out state);
                    alglib.lsfitsetcond(state, epsf, epsx, maxits);
                    alglib.lsfitfit(state, function_cx_1_func, function_cx_1_grad, null, null);
                    alglib.lsfitresults(state, out info, out c, out rep);
                    _TestResult = _TestResult && doc_test_int(info, 2);
                    _TestResult = _TestResult && doc_test_real_vector(c, new double[]{1.5}, 0.05);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "lsfit_d_1");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST lsfit_d_2
            //      Nonlinear fitting by f(x,c) = exp(-c*x^2) (FGH scheme)
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<19; _spoil_scenario++)
            {
                try
                {
                    double[,] x = new double[,]{{-1},{-0.8},{-0.6},{-0.4},{-0.2},{0},{0.2},{0.4},{0.6},{0.8},{1.0}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref x, Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref x, Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref x, Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_matrix_by_deleting_row(ref x);
                    if( _spoil_scenario==4 )
                        spoil_matrix_by_deleting_col(ref x);
                    double[] y = new double[]{0.223130,0.382893,0.582748,0.786628,0.941765,1.000000,0.941765,0.786628,0.582748,0.382893,0.223130};
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref y, Double.NaN);
                    if( _spoil_scenario==6 )
                        spoil_vector_by_value(ref y, Double.PositiveInfinity);
                    if( _spoil_scenario==7 )
                        spoil_vector_by_value(ref y, Double.NegativeInfinity);
                    if( _spoil_scenario==8 )
                        spoil_vector_by_adding_element(ref y, 0);
                    if( _spoil_scenario==9 )
                        spoil_vector_by_deleting_element(ref y);
                    double[] c = new double[]{0.3};
                    if( _spoil_scenario==10 )
                        spoil_vector_by_value(ref c, Double.NaN);
                    if( _spoil_scenario==11 )
                        spoil_vector_by_value(ref c, Double.PositiveInfinity);
                    if( _spoil_scenario==12 )
                        spoil_vector_by_value(ref c, Double.NegativeInfinity);
                    double epsf = 0;
                    if( _spoil_scenario==13 )
                        epsf = Double.NaN;
                    if( _spoil_scenario==14 )
                        epsf = Double.PositiveInfinity;
                    if( _spoil_scenario==15 )
                        epsf = Double.NegativeInfinity;
                    double epsx = 0.000001;
                    if( _spoil_scenario==16 )
                        epsx = Double.NaN;
                    if( _spoil_scenario==17 )
                        epsx = Double.PositiveInfinity;
                    if( _spoil_scenario==18 )
                        epsx = Double.NegativeInfinity;
                    int maxits = 0;
                    int info;
                    alglib.lsfitstate state;
                    alglib.lsfitreport rep;
                    alglib.lsfitcreatefgh(x, y, c, out state);
                    alglib.lsfitsetcond(state, epsf, epsx, maxits);
                    alglib.lsfitfit(state, function_cx_1_func, function_cx_1_grad, function_cx_1_hess, null, null);
                    alglib.lsfitresults(state, out info, out c, out rep);
                    _TestResult = _TestResult && doc_test_int(info, 2);
                    _TestResult = _TestResult && doc_test_real_vector(c, new double[]{1.5}, 0.05);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "lsfit_d_2");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST lsfit_d_3
            //      Nonlinear fitting by f(x,c) = exp(-c*x^2) (weighted FG scheme)
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<24; _spoil_scenario++)
            {
                try
                {
                    double[,] x = new double[,]{{-1},{-0.8},{-0.6},{-0.4},{-0.2},{0},{0.2},{0.4},{0.6},{0.8},{1.0}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref x, Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref x, Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref x, Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_matrix_by_deleting_row(ref x);
                    if( _spoil_scenario==4 )
                        spoil_matrix_by_deleting_col(ref x);
                    double[] y = new double[]{0.223130,0.382893,0.582748,0.786628,0.941765,1.000000,0.941765,0.786628,0.582748,0.382893,0.223130};
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref y, Double.NaN);
                    if( _spoil_scenario==6 )
                        spoil_vector_by_value(ref y, Double.PositiveInfinity);
                    if( _spoil_scenario==7 )
                        spoil_vector_by_value(ref y, Double.NegativeInfinity);
                    if( _spoil_scenario==8 )
                        spoil_vector_by_adding_element(ref y, 0);
                    if( _spoil_scenario==9 )
                        spoil_vector_by_deleting_element(ref y);
                    double[] w = new double[]{1,1,1,1,1,1,1,1,1,1,1};
                    if( _spoil_scenario==10 )
                        spoil_vector_by_value(ref w, Double.NaN);
                    if( _spoil_scenario==11 )
                        spoil_vector_by_value(ref w, Double.PositiveInfinity);
                    if( _spoil_scenario==12 )
                        spoil_vector_by_value(ref w, Double.NegativeInfinity);
                    if( _spoil_scenario==13 )
                        spoil_vector_by_adding_element(ref w, 0);
                    if( _spoil_scenario==14 )
                        spoil_vector_by_deleting_element(ref w);
                    double[] c = new double[]{0.3};
                    if( _spoil_scenario==15 )
                        spoil_vector_by_value(ref c, Double.NaN);
                    if( _spoil_scenario==16 )
                        spoil_vector_by_value(ref c, Double.PositiveInfinity);
                    if( _spoil_scenario==17 )
                        spoil_vector_by_value(ref c, Double.NegativeInfinity);
                    double epsf = 0;
                    if( _spoil_scenario==18 )
                        epsf = Double.NaN;
                    if( _spoil_scenario==19 )
                        epsf = Double.PositiveInfinity;
                    if( _spoil_scenario==20 )
                        epsf = Double.NegativeInfinity;
                    double epsx = 0.000001;
                    if( _spoil_scenario==21 )
                        epsx = Double.NaN;
                    if( _spoil_scenario==22 )
                        epsx = Double.PositiveInfinity;
                    if( _spoil_scenario==23 )
                        epsx = Double.NegativeInfinity;
                    int maxits = 0;
                    int info;
                    alglib.lsfitstate state;
                    alglib.lsfitreport rep;
                    alglib.lsfitcreatewfg(x, y, w, c, true, out state);
                    alglib.lsfitsetcond(state, epsf, epsx, maxits);
                    alglib.lsfitfit(state, function_cx_1_func, function_cx_1_grad, null, null);
                    alglib.lsfitresults(state, out info, out c, out rep);
                    _TestResult = _TestResult && doc_test_int(info, 2);
                    _TestResult = _TestResult && doc_test_real_vector(c, new double[]{1.5}, 0.05);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "lsfit_d_3");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST lsfit_d_4
            //      Nonlinear fitting by f(x,c) = exp(-c*x^2) (weighted FGH scheme)
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<24; _spoil_scenario++)
            {
                try
                {
                    double[,] x = new double[,]{{-1},{-0.8},{-0.6},{-0.4},{-0.2},{0},{0.2},{0.4},{0.6},{0.8},{1.0}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref x, Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref x, Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref x, Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_matrix_by_deleting_row(ref x);
                    if( _spoil_scenario==4 )
                        spoil_matrix_by_deleting_col(ref x);
                    double[] y = new double[]{0.223130,0.382893,0.582748,0.786628,0.941765,1.000000,0.941765,0.786628,0.582748,0.382893,0.223130};
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref y, Double.NaN);
                    if( _spoil_scenario==6 )
                        spoil_vector_by_value(ref y, Double.PositiveInfinity);
                    if( _spoil_scenario==7 )
                        spoil_vector_by_value(ref y, Double.NegativeInfinity);
                    if( _spoil_scenario==8 )
                        spoil_vector_by_adding_element(ref y, 0);
                    if( _spoil_scenario==9 )
                        spoil_vector_by_deleting_element(ref y);
                    double[] w = new double[]{1,1,1,1,1,1,1,1,1,1,1};
                    if( _spoil_scenario==10 )
                        spoil_vector_by_value(ref w, Double.NaN);
                    if( _spoil_scenario==11 )
                        spoil_vector_by_value(ref w, Double.PositiveInfinity);
                    if( _spoil_scenario==12 )
                        spoil_vector_by_value(ref w, Double.NegativeInfinity);
                    if( _spoil_scenario==13 )
                        spoil_vector_by_adding_element(ref w, 0);
                    if( _spoil_scenario==14 )
                        spoil_vector_by_deleting_element(ref w);
                    double[] c = new double[]{0.3};
                    if( _spoil_scenario==15 )
                        spoil_vector_by_value(ref c, Double.NaN);
                    if( _spoil_scenario==16 )
                        spoil_vector_by_value(ref c, Double.PositiveInfinity);
                    if( _spoil_scenario==17 )
                        spoil_vector_by_value(ref c, Double.NegativeInfinity);
                    double epsf = 0;
                    if( _spoil_scenario==18 )
                        epsf = Double.NaN;
                    if( _spoil_scenario==19 )
                        epsf = Double.PositiveInfinity;
                    if( _spoil_scenario==20 )
                        epsf = Double.NegativeInfinity;
                    double epsx = 0.000001;
                    if( _spoil_scenario==21 )
                        epsx = Double.NaN;
                    if( _spoil_scenario==22 )
                        epsx = Double.PositiveInfinity;
                    if( _spoil_scenario==23 )
                        epsx = Double.NegativeInfinity;
                    int maxits = 0;
                    int info;
                    alglib.lsfitstate state;
                    alglib.lsfitreport rep;
                    alglib.lsfitcreatewfgh(x, y, w, c, out state);
                    alglib.lsfitsetcond(state, epsf, epsx, maxits);
                    alglib.lsfitfit(state, function_cx_1_func, function_cx_1_grad, function_cx_1_hess, null, null);
                    alglib.lsfitresults(state, out info, out c, out rep);
                    _TestResult = _TestResult && doc_test_int(info, 2);
                    _TestResult = _TestResult && doc_test_real_vector(c, new double[]{1.5}, 0.05);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "lsfit_d_4");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST lsfit_d_5
            //      Linear fitting by f(x,a) = a*exp(0.5*x)
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<8; _spoil_scenario++)
            {
                try
                {
                    double[,] fmatrix = new double[,]{{0.606531},{0.670320},{0.740818},{0.818731},{0.904837},{1.000000},{1.105171},{1.221403},{1.349859},{1.491825},{1.648721}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref fmatrix, Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref fmatrix, Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref fmatrix, Double.NegativeInfinity);
                    double[] y = new double[]{1.133719,1.306522,1.504604,1.554663,1.884638,2.072436,2.257285,2.534068,2.622017,2.897713,3.219371};
                    if( _spoil_scenario==3 )
                        spoil_vector_by_value(ref y, Double.NaN);
                    if( _spoil_scenario==4 )
                        spoil_vector_by_value(ref y, Double.PositiveInfinity);
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref y, Double.NegativeInfinity);
                    if( _spoil_scenario==6 )
                        spoil_vector_by_adding_element(ref y, 0);
                    if( _spoil_scenario==7 )
                        spoil_vector_by_deleting_element(ref y);
                    int info;
                    double[] c;
                    alglib.lsfitreport rep;
                    alglib.lsfitlinear(y, fmatrix, out info, out c, out rep);
                    _TestResult = _TestResult && doc_test_int(info, 1);
                    _TestResult = _TestResult && doc_test_real_vector(c, new double[]{1.98650}, 0.00005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "lsfit_d_5");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST lsfit_d_6
            //      Linear fitting by f(x,a) = a*exp(0.5*x), nonuniform weights
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<13; _spoil_scenario++)
            {
                try
                {
                    double[,] fmatrix = new double[,]{{0.606531},{0.670320},{0.740818},{0.818731},{0.904837},{1.000000},{1.105171},{1.221403},{1.349859},{1.491825},{1.648721}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref fmatrix, Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref fmatrix, Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref fmatrix, Double.NegativeInfinity);
                    double[] y = new double[]{1.133719,1.306522,1.504604,1.554663,1.884638,2.072436,2.257285,2.534068,2.622017,2.897713,3.219371};
                    if( _spoil_scenario==3 )
                        spoil_vector_by_value(ref y, Double.NaN);
                    if( _spoil_scenario==4 )
                        spoil_vector_by_value(ref y, Double.PositiveInfinity);
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref y, Double.NegativeInfinity);
                    if( _spoil_scenario==6 )
                        spoil_vector_by_adding_element(ref y, 0);
                    if( _spoil_scenario==7 )
                        spoil_vector_by_deleting_element(ref y);
                    double[] w = new double[]{1.414213,1,1,1,1,1,1,1,1,1,1};
                    if( _spoil_scenario==8 )
                        spoil_vector_by_value(ref w, Double.NaN);
                    if( _spoil_scenario==9 )
                        spoil_vector_by_value(ref w, Double.PositiveInfinity);
                    if( _spoil_scenario==10 )
                        spoil_vector_by_value(ref w, Double.NegativeInfinity);
                    if( _spoil_scenario==11 )
                        spoil_vector_by_adding_element(ref w, 0);
                    if( _spoil_scenario==12 )
                        spoil_vector_by_deleting_element(ref w);
                    int info;
                    double[] c;
                    alglib.lsfitreport rep;
                    alglib.lsfitlinearw(y, w, fmatrix, out info, out c, out rep);
                    _TestResult = _TestResult && doc_test_int(info, 1);
                    _TestResult = _TestResult && doc_test_real_vector(c, new double[]{1.983354}, 0.00005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "lsfit_d_6");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST lsfit_d_7
            //      Linear fitting by f(x,a,b) = a*x+b, f(0)=0
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<15; _spoil_scenario++)
            {
                try
                {
                    double[] y = new double[]{0.072436,0.246944,0.491263,0.522300,0.714064,0.921929};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref y, Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref y, Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref y, Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_vector_by_adding_element(ref y, 0);
                    if( _spoil_scenario==4 )
                        spoil_vector_by_deleting_element(ref y);
                    double[,] fmatrix = new double[,]{{1,0.0},{1,0.2},{1,0.4},{1,0.6},{1,0.8},{1,1.0}};
                    if( _spoil_scenario==5 )
                        spoil_matrix_by_value(ref fmatrix, Double.NaN);
                    if( _spoil_scenario==6 )
                        spoil_matrix_by_value(ref fmatrix, Double.PositiveInfinity);
                    if( _spoil_scenario==7 )
                        spoil_matrix_by_value(ref fmatrix, Double.NegativeInfinity);
                    if( _spoil_scenario==8 )
                        spoil_matrix_by_adding_row(ref fmatrix, 0);
                    if( _spoil_scenario==9 )
                        spoil_matrix_by_adding_col(ref fmatrix, 0);
                    if( _spoil_scenario==10 )
                        spoil_matrix_by_deleting_row(ref fmatrix);
                    if( _spoil_scenario==11 )
                        spoil_matrix_by_deleting_col(ref fmatrix);
                    double[,] cmatrix = new double[,]{{1,0,0}};
                    if( _spoil_scenario==12 )
                        spoil_matrix_by_value(ref cmatrix, Double.NaN);
                    if( _spoil_scenario==13 )
                        spoil_matrix_by_value(ref cmatrix, Double.PositiveInfinity);
                    if( _spoil_scenario==14 )
                        spoil_matrix_by_value(ref cmatrix, Double.NegativeInfinity);
                    int info;
                    double[] c;
                    alglib.lsfitreport rep;
                    alglib.lsfitlinearc(y, fmatrix, cmatrix, out info, out c, out rep);
                    _TestResult = _TestResult && doc_test_int(info, 1);
                    _TestResult = _TestResult && doc_test_real_vector(c, new double[]{0,0.932933}, 0.0005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "lsfit_d_7");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST lsfit_d_8
            //      Linear fitting by f(x,a,b) = a*x+b, f(0)=0, nonuniform weights
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<20; _spoil_scenario++)
            {
                try
                {
                    double[] y = new double[]{0.072436,0.246944,0.491263,0.522300,0.714064,0.921929};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref y, Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref y, Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref y, Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_vector_by_adding_element(ref y, 0);
                    if( _spoil_scenario==4 )
                        spoil_vector_by_deleting_element(ref y);
                    double[] w = new double[]{1,1.414213,1,1,1,1};
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref w, Double.NaN);
                    if( _spoil_scenario==6 )
                        spoil_vector_by_value(ref w, Double.PositiveInfinity);
                    if( _spoil_scenario==7 )
                        spoil_vector_by_value(ref w, Double.NegativeInfinity);
                    if( _spoil_scenario==8 )
                        spoil_vector_by_adding_element(ref w, 0);
                    if( _spoil_scenario==9 )
                        spoil_vector_by_deleting_element(ref w);
                    double[,] fmatrix = new double[,]{{1,0.0},{1,0.2},{1,0.4},{1,0.6},{1,0.8},{1,1.0}};
                    if( _spoil_scenario==10 )
                        spoil_matrix_by_value(ref fmatrix, Double.NaN);
                    if( _spoil_scenario==11 )
                        spoil_matrix_by_value(ref fmatrix, Double.PositiveInfinity);
                    if( _spoil_scenario==12 )
                        spoil_matrix_by_value(ref fmatrix, Double.NegativeInfinity);
                    if( _spoil_scenario==13 )
                        spoil_matrix_by_adding_row(ref fmatrix, 0);
                    if( _spoil_scenario==14 )
                        spoil_matrix_by_adding_col(ref fmatrix, 0);
                    if( _spoil_scenario==15 )
                        spoil_matrix_by_deleting_row(ref fmatrix);
                    if( _spoil_scenario==16 )
                        spoil_matrix_by_deleting_col(ref fmatrix);
                    double[,] cmatrix = new double[,]{{1,0,0}};
                    if( _spoil_scenario==17 )
                        spoil_matrix_by_value(ref cmatrix, Double.NaN);
                    if( _spoil_scenario==18 )
                        spoil_matrix_by_value(ref cmatrix, Double.PositiveInfinity);
                    if( _spoil_scenario==19 )
                        spoil_matrix_by_value(ref cmatrix, Double.NegativeInfinity);
                    int info;
                    double[] c;
                    alglib.lsfitreport rep;
                    alglib.lsfitlinearwc(y, w, fmatrix, cmatrix, out info, out c, out rep);
                    _TestResult = _TestResult && doc_test_int(info, 1);
                    _TestResult = _TestResult && doc_test_real_vector(c, new double[]{0,0.938322}, 0.0005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "lsfit_d_8");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST minlm_d_1
            //      Nonlinear optimization, F(x0,x1) = 100*(x0+3)^4+(x1-3)^4, posed as Levenberg-Marquardt problem (FJ scheme)
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<12; _spoil_scenario++)
            {
                try
                {
                    //
                    // This example demonstrates minimization of f(x0,x1) = 100*(x0+3)^4+(x1-3)^4
                    // using classic Levenberg-Marquardt problem (FJ scheme).
                    //
                    // F is treated like a sum of squares of two functions,
                    // f0=10*(x0+3)^2, and f1=(x1-3)^2.
                    //
                    // Optimization algorithm uses:
                    // * function value F(x0,x1)=f0^2+f1^2
                    // * Jacobian matrix J={dfi/dxj}.
                    //
                    double[] x = new double[]{0,0};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, Double.NegativeInfinity);
                    double epsg = 0.0000000001;
                    if( _spoil_scenario==3 )
                        epsg = Double.NaN;
                    if( _spoil_scenario==4 )
                        epsg = Double.PositiveInfinity;
                    if( _spoil_scenario==5 )
                        epsg = Double.NegativeInfinity;
                    double epsf = 0;
                    if( _spoil_scenario==6 )
                        epsf = Double.NaN;
                    if( _spoil_scenario==7 )
                        epsf = Double.PositiveInfinity;
                    if( _spoil_scenario==8 )
                        epsf = Double.NegativeInfinity;
                    double epsx = 0;
                    if( _spoil_scenario==9 )
                        epsx = Double.NaN;
                    if( _spoil_scenario==10 )
                        epsx = Double.PositiveInfinity;
                    if( _spoil_scenario==11 )
                        epsx = Double.NegativeInfinity;
                    int maxits = 0;
                    alglib.minlmstate state;
                    alglib.minlmreport rep;

                    alglib.minlmcreatefj(2, x, out state);
                    alglib.minlmsetcond(state, epsg, epsf, epsx, maxits);
                    alglib.minlmoptimize(state, function1_func, function1_jac, null, null);
                    alglib.minlmresults(state, out x, out rep);

                    _TestResult = _TestResult && doc_test_int(rep.terminationtype, 4);
                    _TestResult = _TestResult && doc_test_real_vector(x, new double[]{-3,+3}, 0.005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "minlm_d_1");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST minlm_d_2
            //      Nonlinear optimization, F(x,y) = 100*(x+3)^4+(y-3)^4, posed as Levenberg-Marquardt problem (FGJ scheme)
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<12; _spoil_scenario++)
            {
                try
                {
                    //
                    // This example demonstrates minimization of f(x0,x1) = 100*(x0+3)^4+(x1-3)^4
                    // using improved Levenberg-Marquardt problem (FGJ scheme).
                    //
                    // F is treated like a sum of squares of two functions,
                    // f0=10*(x0+3)^2, and f1=(x1-3)^2.
                    //
                    // Optimization algorithm uses:
                    // * function value F(x0,x1)=f0^2+f1^2
                    // * Jacobian matrix J={dfi/dxj}
                    // * gradient G={dF/dxi}
                    //
                    double[] x = new double[]{0,0};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, Double.NegativeInfinity);
                    double epsg = 0.0000000001;
                    if( _spoil_scenario==3 )
                        epsg = Double.NaN;
                    if( _spoil_scenario==4 )
                        epsg = Double.PositiveInfinity;
                    if( _spoil_scenario==5 )
                        epsg = Double.NegativeInfinity;
                    double epsf = 0;
                    if( _spoil_scenario==6 )
                        epsf = Double.NaN;
                    if( _spoil_scenario==7 )
                        epsf = Double.PositiveInfinity;
                    if( _spoil_scenario==8 )
                        epsf = Double.NegativeInfinity;
                    double epsx = 0;
                    if( _spoil_scenario==9 )
                        epsx = Double.NaN;
                    if( _spoil_scenario==10 )
                        epsx = Double.PositiveInfinity;
                    if( _spoil_scenario==11 )
                        epsx = Double.NegativeInfinity;
                    int maxits = 0;
                    alglib.minlmstate state;
                    alglib.minlmreport rep;

                    alglib.minlmcreatefgj(2, x, out state);
                    alglib.minlmsetcond(state, epsg, epsf, epsx, maxits);
                    alglib.minlmoptimize(state, function1_func, function1_grad, function1_jac, null, null);
                    alglib.minlmresults(state, out x, out rep);

                    _TestResult = _TestResult && doc_test_int(rep.terminationtype, 4);
                    _TestResult = _TestResult && doc_test_real_vector(x, new double[]{-3,+3}, 0.005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "minlm_d_2");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST minlm_d_3
            //      Nonlinear optimization, F(x,y) = 100*(x+3)^4+(y-3)^4, posed as Levenberg-Marquardt problem (FGH scheme)
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<12; _spoil_scenario++)
            {
                try
                {
                    //
                    // This example demonstrates minimization of f(x0,x1) = 100*(x0+3)^4+(x1-3)^4
                    // using improved Levenberg-Marquardt problem (FGH scheme).
                    //
                    // F is treated like a monolitic function without internal structure.
                    // I.e. we do NOT represent it as a sum of squares.
                    //
                    // Optimization algorithm uses:
                    // * function value F(x0,x1)=f0^2+f1^2
                    // * gradient G={dF/dxi}
                    // * Hessian H={d2F/(dxi*dxj)}
                    //
                    double[] x = new double[]{0,0};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, Double.NegativeInfinity);
                    double epsg = 0.0000000001;
                    if( _spoil_scenario==3 )
                        epsg = Double.NaN;
                    if( _spoil_scenario==4 )
                        epsg = Double.PositiveInfinity;
                    if( _spoil_scenario==5 )
                        epsg = Double.NegativeInfinity;
                    double epsf = 0;
                    if( _spoil_scenario==6 )
                        epsf = Double.NaN;
                    if( _spoil_scenario==7 )
                        epsf = Double.PositiveInfinity;
                    if( _spoil_scenario==8 )
                        epsf = Double.NegativeInfinity;
                    double epsx = 0;
                    if( _spoil_scenario==9 )
                        epsx = Double.NaN;
                    if( _spoil_scenario==10 )
                        epsx = Double.PositiveInfinity;
                    if( _spoil_scenario==11 )
                        epsx = Double.NegativeInfinity;
                    int maxits = 0;
                    alglib.minlmstate state;
                    alglib.minlmreport rep;

                    alglib.minlmcreatefgh(x, out state);
                    alglib.minlmsetcond(state, epsg, epsf, epsx, maxits);
                    alglib.minlmoptimize(state, function1_func, function1_grad, function1_hess, null, null);
                    alglib.minlmresults(state, out x, out rep);

                    _TestResult = _TestResult && doc_test_int(rep.terminationtype, 4);
                    _TestResult = _TestResult && doc_test_real_vector(x, new double[]{-3,+3}, 0.005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "minlm_d_3");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST polint_d_1a
            //      Polynomial interpolation: y=x^2-x, general grid, barycentric form
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<12; _spoil_scenario++)
            {
                try
                {
                    double[] x = new double[]{0,1,2};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_vector_by_adding_element(ref x, 0);
                    if( _spoil_scenario==4 )
                        spoil_vector_by_deleting_element(ref x);
                    double[] y = new double[]{0,0,2};
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref y, Double.NaN);
                    if( _spoil_scenario==6 )
                        spoil_vector_by_value(ref y, Double.PositiveInfinity);
                    if( _spoil_scenario==7 )
                        spoil_vector_by_value(ref y, Double.NegativeInfinity);
                    if( _spoil_scenario==8 )
                        spoil_vector_by_adding_element(ref y, 0);
                    if( _spoil_scenario==9 )
                        spoil_vector_by_deleting_element(ref y);
                    double t = -1;
                    if( _spoil_scenario==10 )
                        t = Double.PositiveInfinity;
                    if( _spoil_scenario==11 )
                        t = Double.NegativeInfinity;
                    double v;
                    alglib.barycentricinterpolant p;
                    alglib.polynomialbuild(x, y, out p);
                    v = alglib.barycentriccalc(p, t);
                    _TestResult = _TestResult && doc_test_real(v, 2.0, 0.00005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "polint_d_1a");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST polint_d_1b
            //      Polynomial differentiation: y=x^2-x, general grid, barycentric form
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<12; _spoil_scenario++)
            {
                try
                {
                    double[] x = new double[]{0,1,2};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_vector_by_adding_element(ref x, 0);
                    if( _spoil_scenario==4 )
                        spoil_vector_by_deleting_element(ref x);
                    double[] y = new double[]{0,0,2};
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref y, Double.NaN);
                    if( _spoil_scenario==6 )
                        spoil_vector_by_value(ref y, Double.PositiveInfinity);
                    if( _spoil_scenario==7 )
                        spoil_vector_by_value(ref y, Double.NegativeInfinity);
                    if( _spoil_scenario==8 )
                        spoil_vector_by_adding_element(ref y, 0);
                    if( _spoil_scenario==9 )
                        spoil_vector_by_deleting_element(ref y);
                    double t = -1;
                    if( _spoil_scenario==10 )
                        t = Double.PositiveInfinity;
                    if( _spoil_scenario==11 )
                        t = Double.NegativeInfinity;
                    double v;
                    double dv;
                    alglib.barycentricinterpolant p;
                    alglib.polynomialbuild(x, y, out p);
                    alglib.barycentricdiff1(p, t, out v, out dv);
                    _TestResult = _TestResult && doc_test_real(v, 2.0, 0.00005);
                    _TestResult = _TestResult && doc_test_real(dv, -3.0, 0.00005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "polint_d_1b");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST polint_d_1c
            //      Polynomial differentiation (2): y=x^2-x, general grid, barycentric form
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<12; _spoil_scenario++)
            {
                try
                {
                    double[] x = new double[]{0,1,2};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_vector_by_adding_element(ref x, 0);
                    if( _spoil_scenario==4 )
                        spoil_vector_by_deleting_element(ref x);
                    double[] y = new double[]{0,0,2};
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref y, Double.NaN);
                    if( _spoil_scenario==6 )
                        spoil_vector_by_value(ref y, Double.PositiveInfinity);
                    if( _spoil_scenario==7 )
                        spoil_vector_by_value(ref y, Double.NegativeInfinity);
                    if( _spoil_scenario==8 )
                        spoil_vector_by_adding_element(ref y, 0);
                    if( _spoil_scenario==9 )
                        spoil_vector_by_deleting_element(ref y);
                    double t = -1;
                    if( _spoil_scenario==10 )
                        t = Double.PositiveInfinity;
                    if( _spoil_scenario==11 )
                        t = Double.NegativeInfinity;
                    double v;
                    double dv;
                    double d2v;
                    alglib.barycentricinterpolant p;
                    alglib.polynomialbuild(x, y, out p);
                    alglib.barycentricdiff2(p, t, out v, out dv, out d2v);
                    _TestResult = _TestResult && doc_test_real(v, 2.0, 0.00005);
                    _TestResult = _TestResult && doc_test_real(dv, -3.0, 0.00005);
                    _TestResult = _TestResult && doc_test_real(d2v, 2.0, 0.00005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "polint_d_1c");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST polint_d_2
            //      Polynomial interpolation: y=x^2-x, equidistant grid, barycentric form
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<5; _spoil_scenario++)
            {
                try
                {
                    double[] y = new double[]{0,0,2};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref y, Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref y, Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref y, Double.NegativeInfinity);
                    double t = -1;
                    if( _spoil_scenario==3 )
                        t = Double.PositiveInfinity;
                    if( _spoil_scenario==4 )
                        t = Double.NegativeInfinity;
                    alglib.barycentricinterpolant p;
                    double v;
                    alglib.polynomialbuildeqdist(0.0, 2.0, y, out p);
                    v = alglib.barycentriccalc(p, t);
                    _TestResult = _TestResult && doc_test_real(v, 2.0, 0.00005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "polint_d_2");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST polint_d_3
            //      Polynomial interpolation: y=x^2-x, Chebyshev grid (first kind), barycentric form
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<11; _spoil_scenario++)
            {
                try
                {
                    double[] y = new double[]{-0.116025,0.000000,1.616025};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref y, Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref y, Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref y, Double.NegativeInfinity);
                    double t = -1;
                    if( _spoil_scenario==3 )
                        t = Double.PositiveInfinity;
                    if( _spoil_scenario==4 )
                        t = Double.NegativeInfinity;
                    double a = -1;
                    if( _spoil_scenario==5 )
                        a = Double.NaN;
                    if( _spoil_scenario==6 )
                        a = Double.PositiveInfinity;
                    if( _spoil_scenario==7 )
                        a = Double.NegativeInfinity;
                    double b = +1;
                    if( _spoil_scenario==8 )
                        b = Double.NaN;
                    if( _spoil_scenario==9 )
                        b = Double.PositiveInfinity;
                    if( _spoil_scenario==10 )
                        b = Double.NegativeInfinity;
                    alglib.barycentricinterpolant p;
                    double v;
                    alglib.polynomialbuildcheb1(a, b, y, out p);
                    v = alglib.barycentriccalc(p, t);
                    _TestResult = _TestResult && doc_test_real(v, 2.0, 0.00005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "polint_d_3");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST polint_d_4
            //      Polynomial interpolation: y=x^2-x, Chebyshev grid (second kind), barycentric form
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<11; _spoil_scenario++)
            {
                try
                {
                    double[] y = new double[]{0,0,2};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref y, Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref y, Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref y, Double.NegativeInfinity);
                    double t = -2;
                    if( _spoil_scenario==3 )
                        t = Double.PositiveInfinity;
                    if( _spoil_scenario==4 )
                        t = Double.NegativeInfinity;
                    double a = -1;
                    if( _spoil_scenario==5 )
                        a = Double.NaN;
                    if( _spoil_scenario==6 )
                        a = Double.PositiveInfinity;
                    if( _spoil_scenario==7 )
                        a = Double.NegativeInfinity;
                    double b = +1;
                    if( _spoil_scenario==8 )
                        b = Double.NaN;
                    if( _spoil_scenario==9 )
                        b = Double.PositiveInfinity;
                    if( _spoil_scenario==10 )
                        b = Double.NegativeInfinity;
                    alglib.barycentricinterpolant p;
                    double v;
                    alglib.polynomialbuildcheb2(a, b, y, out p);
                    v = alglib.barycentriccalc(p, t);
                    _TestResult = _TestResult && doc_test_real(v, 6.0, 0.00005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "polint_d_4");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST polint_d_5
            //      Polynomial interpolation: y=x^2-x, equidistant grid
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<5; _spoil_scenario++)
            {
                try
                {
                    double[] y = new double[]{0,0,2};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref y, Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref y, Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref y, Double.NegativeInfinity);
                    double t = -1;
                    if( _spoil_scenario==3 )
                        t = Double.PositiveInfinity;
                    if( _spoil_scenario==4 )
                        t = Double.NegativeInfinity;
                    double v;
                    v = alglib.polynomialcalceqdist(0.0, 2.0, y, t);
                    _TestResult = _TestResult && doc_test_real(v, 2.0, 0.00005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "polint_d_5");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST polint_d_6
            //      Polynomial interpolation: y=x^2-x, Chebyshev grid (first kind)
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<11; _spoil_scenario++)
            {
                try
                {
                    double[] y = new double[]{-0.116025,0.000000,1.616025};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref y, Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref y, Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref y, Double.NegativeInfinity);
                    double t = -1;
                    if( _spoil_scenario==3 )
                        t = Double.PositiveInfinity;
                    if( _spoil_scenario==4 )
                        t = Double.NegativeInfinity;
                    double a = -1;
                    if( _spoil_scenario==5 )
                        a = Double.NaN;
                    if( _spoil_scenario==6 )
                        a = Double.PositiveInfinity;
                    if( _spoil_scenario==7 )
                        a = Double.NegativeInfinity;
                    double b = +1;
                    if( _spoil_scenario==8 )
                        b = Double.NaN;
                    if( _spoil_scenario==9 )
                        b = Double.PositiveInfinity;
                    if( _spoil_scenario==10 )
                        b = Double.NegativeInfinity;
                    double v;
                    v = alglib.polynomialcalccheb1(a, b, y, t);
                    _TestResult = _TestResult && doc_test_real(v, 2.0, 0.00005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "polint_d_6");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST polint_d_7
            //      Polynomial interpolation: y=x^2-x, Chebyshev grid (second kind)
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<11; _spoil_scenario++)
            {
                try
                {
                    double[] y = new double[]{0,0,2};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref y, Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref y, Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref y, Double.NegativeInfinity);
                    double t = -2;
                    if( _spoil_scenario==3 )
                        t = Double.PositiveInfinity;
                    if( _spoil_scenario==4 )
                        t = Double.NegativeInfinity;
                    double a = -1;
                    if( _spoil_scenario==5 )
                        a = Double.NaN;
                    if( _spoil_scenario==6 )
                        a = Double.PositiveInfinity;
                    if( _spoil_scenario==7 )
                        a = Double.NegativeInfinity;
                    double b = +1;
                    if( _spoil_scenario==8 )
                        b = Double.NaN;
                    if( _spoil_scenario==9 )
                        b = Double.PositiveInfinity;
                    if( _spoil_scenario==10 )
                        b = Double.NegativeInfinity;
                    double v;
                    v = alglib.polynomialcalccheb2(a, b, y, t);
                    _TestResult = _TestResult && doc_test_real(v, 6.0, 0.00005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "polint_d_7");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST polint_d_8
            //      Polynomial fitting, simple form. Linear fit is used (M=2).
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<12; _spoil_scenario++)
            {
                try
                {
                    double[] x = new double[]{0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_vector_by_adding_element(ref x, 0);
                    if( _spoil_scenario==4 )
                        spoil_vector_by_deleting_element(ref x);
                    double[] y = new double[]{0.00,0.05,0.26,0.32,0.33,0.43,0.60,0.60,0.77,0.98,1.02};
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref y, Double.NaN);
                    if( _spoil_scenario==6 )
                        spoil_vector_by_value(ref y, Double.PositiveInfinity);
                    if( _spoil_scenario==7 )
                        spoil_vector_by_value(ref y, Double.NegativeInfinity);
                    if( _spoil_scenario==8 )
                        spoil_vector_by_adding_element(ref y, 0);
                    if( _spoil_scenario==9 )
                        spoil_vector_by_deleting_element(ref y);
                    int m = 2;
                    double t = 2;
                    if( _spoil_scenario==10 )
                        t = Double.PositiveInfinity;
                    if( _spoil_scenario==11 )
                        t = Double.NegativeInfinity;
                    int info;
                    alglib.barycentricinterpolant p;
                    alglib.polynomialfitreport rep;
                    double v;
                    alglib.polynomialfit(x, y, m, out info, out p, out rep);
                    v = alglib.barycentriccalc(p, t);
                    _TestResult = _TestResult && doc_test_real(v, 2.011, 0.002);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "polint_d_8");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST polint_d_9
            //      Polynomial fitting, weighted form. Linear fit is used: M=2.
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<20; _spoil_scenario++)
            {
                try
                {
                    double[] x = new double[]{0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_vector_by_adding_element(ref x, 0);
                    if( _spoil_scenario==4 )
                        spoil_vector_by_deleting_element(ref x);
                    double[] y = new double[]{0.00,0.05,0.26,0.32,0.33,0.43,0.60,0.60,0.77,0.98,1.02};
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref y, Double.NaN);
                    if( _spoil_scenario==6 )
                        spoil_vector_by_value(ref y, Double.PositiveInfinity);
                    if( _spoil_scenario==7 )
                        spoil_vector_by_value(ref y, Double.NegativeInfinity);
                    if( _spoil_scenario==8 )
                        spoil_vector_by_adding_element(ref y, 0);
                    if( _spoil_scenario==9 )
                        spoil_vector_by_deleting_element(ref y);
                    double[] w = new double[]{1,1.414213562,1,1,1,1,1,1,1,1,1};
                    if( _spoil_scenario==10 )
                        spoil_vector_by_value(ref w, Double.NaN);
                    if( _spoil_scenario==11 )
                        spoil_vector_by_value(ref w, Double.PositiveInfinity);
                    if( _spoil_scenario==12 )
                        spoil_vector_by_value(ref w, Double.NegativeInfinity);
                    if( _spoil_scenario==13 )
                        spoil_vector_by_adding_element(ref w, 0);
                    if( _spoil_scenario==14 )
                        spoil_vector_by_deleting_element(ref w);
                    double[] xc = new double[]{};
                    if( _spoil_scenario==15 )
                        spoil_vector_by_adding_element(ref xc, 0);
                    double[] yc = new double[]{};
                    if( _spoil_scenario==16 )
                        spoil_vector_by_adding_element(ref yc, 0);
                    int[] dc = new int[]{};
                    if( _spoil_scenario==17 )
                        spoil_vector_by_adding_element(ref dc, 0);
                    int m = 2;
                    double t = 2;
                    if( _spoil_scenario==18 )
                        t = Double.PositiveInfinity;
                    if( _spoil_scenario==19 )
                        t = Double.NegativeInfinity;
                    int info;
                    alglib.barycentricinterpolant p;
                    alglib.polynomialfitreport rep;
                    double v;
                    alglib.polynomialfitwc(x, y, w, xc, yc, dc, m, out info, out p, out rep);
                    v = alglib.barycentriccalc(p, t);
                    _TestResult = _TestResult && doc_test_real(v, 2.023, 0.002);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "polint_d_9");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST polint_d_10
            //      Polynomial fitting with constrants. Linear fit is used: M=2.
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<29; _spoil_scenario++)
            {
                try
                {
                    double[] x = new double[]{1.0,1.0};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_vector_by_adding_element(ref x, 0);
                    if( _spoil_scenario==4 )
                        spoil_vector_by_deleting_element(ref x);
                    double[] y = new double[]{0.9,1.1};
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref y, Double.NaN);
                    if( _spoil_scenario==6 )
                        spoil_vector_by_value(ref y, Double.PositiveInfinity);
                    if( _spoil_scenario==7 )
                        spoil_vector_by_value(ref y, Double.NegativeInfinity);
                    if( _spoil_scenario==8 )
                        spoil_vector_by_adding_element(ref y, 0);
                    if( _spoil_scenario==9 )
                        spoil_vector_by_deleting_element(ref y);
                    double[] w = new double[]{1,1};
                    if( _spoil_scenario==10 )
                        spoil_vector_by_value(ref w, Double.NaN);
                    if( _spoil_scenario==11 )
                        spoil_vector_by_value(ref w, Double.PositiveInfinity);
                    if( _spoil_scenario==12 )
                        spoil_vector_by_value(ref w, Double.NegativeInfinity);
                    if( _spoil_scenario==13 )
                        spoil_vector_by_adding_element(ref w, 0);
                    if( _spoil_scenario==14 )
                        spoil_vector_by_deleting_element(ref w);
                    double[] xc = new double[]{0};
                    if( _spoil_scenario==15 )
                        spoil_vector_by_value(ref xc, Double.NaN);
                    if( _spoil_scenario==16 )
                        spoil_vector_by_value(ref xc, Double.PositiveInfinity);
                    if( _spoil_scenario==17 )
                        spoil_vector_by_value(ref xc, Double.NegativeInfinity);
                    if( _spoil_scenario==18 )
                        spoil_vector_by_adding_element(ref xc, 0);
                    if( _spoil_scenario==19 )
                        spoil_vector_by_deleting_element(ref xc);
                    double[] yc = new double[]{0};
                    if( _spoil_scenario==20 )
                        spoil_vector_by_value(ref yc, Double.NaN);
                    if( _spoil_scenario==21 )
                        spoil_vector_by_value(ref yc, Double.PositiveInfinity);
                    if( _spoil_scenario==22 )
                        spoil_vector_by_value(ref yc, Double.NegativeInfinity);
                    if( _spoil_scenario==23 )
                        spoil_vector_by_adding_element(ref yc, 0);
                    if( _spoil_scenario==24 )
                        spoil_vector_by_deleting_element(ref yc);
                    int[] dc = new int[]{0};
                    if( _spoil_scenario==25 )
                        spoil_vector_by_adding_element(ref dc, 0);
                    if( _spoil_scenario==26 )
                        spoil_vector_by_deleting_element(ref dc);
                    double t = 2;
                    if( _spoil_scenario==27 )
                        t = Double.PositiveInfinity;
                    if( _spoil_scenario==28 )
                        t = Double.NegativeInfinity;
                    int m = 2;
                    int info;
                    alglib.barycentricinterpolant p;
                    alglib.polynomialfitreport rep;
                    double v;
                    alglib.polynomialfitwc(x, y, w, xc, yc, dc, m, out info, out p, out rep);
                    v = alglib.barycentriccalc(p, t);
                    _TestResult = _TestResult && doc_test_real(v, 2.000, 0.001);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "polint_d_10");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST polint_t_1
            //      Polynomial interpolation, full list of parameters.
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<10; _spoil_scenario++)
            {
                try
                {
                    double[] x = new double[]{0,1,2};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_vector_by_deleting_element(ref x);
                    double[] y = new double[]{0,0,2};
                    if( _spoil_scenario==4 )
                        spoil_vector_by_value(ref y, Double.NaN);
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref y, Double.PositiveInfinity);
                    if( _spoil_scenario==6 )
                        spoil_vector_by_value(ref y, Double.NegativeInfinity);
                    if( _spoil_scenario==7 )
                        spoil_vector_by_deleting_element(ref y);
                    double t = -1;
                    if( _spoil_scenario==8 )
                        t = Double.PositiveInfinity;
                    if( _spoil_scenario==9 )
                        t = Double.NegativeInfinity;
                    alglib.barycentricinterpolant p;
                    double v;
                    alglib.polynomialbuild(x, y, 3, out p);
                    v = alglib.barycentriccalc(p, t);
                    _TestResult = _TestResult && doc_test_real(v, 2.0, 0.00005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "polint_t_1");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST polint_t_2
            //      Polynomial interpolation, full list of parameters.
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<6; _spoil_scenario++)
            {
                try
                {
                    double[] y = new double[]{0,0,2};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref y, Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref y, Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref y, Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_vector_by_deleting_element(ref y);
                    double t = -1;
                    if( _spoil_scenario==4 )
                        t = Double.PositiveInfinity;
                    if( _spoil_scenario==5 )
                        t = Double.NegativeInfinity;
                    alglib.barycentricinterpolant p;
                    double v;
                    alglib.polynomialbuildeqdist(0.0, 2.0, y, 3, out p);
                    v = alglib.barycentriccalc(p, t);
                    _TestResult = _TestResult && doc_test_real(v, 2.0, 0.00005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "polint_t_2");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST polint_t_3
            //      Polynomial interpolation, full list of parameters.
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<6; _spoil_scenario++)
            {
                try
                {
                    double[] y = new double[]{-0.116025,0.000000,1.616025};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref y, Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref y, Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref y, Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_vector_by_deleting_element(ref y);
                    double t = -1;
                    if( _spoil_scenario==4 )
                        t = Double.PositiveInfinity;
                    if( _spoil_scenario==5 )
                        t = Double.NegativeInfinity;
                    alglib.barycentricinterpolant p;
                    double v;
                    alglib.polynomialbuildcheb1(-1.0, +1.0, y, 3, out p);
                    v = alglib.barycentriccalc(p, t);
                    _TestResult = _TestResult && doc_test_real(v, 2.0, 0.00005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "polint_t_3");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST polint_t_4
            //      Polynomial interpolation, full list of parameters.
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<12; _spoil_scenario++)
            {
                try
                {
                    double[] y = new double[]{0,0,2};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref y, Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref y, Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref y, Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_vector_by_deleting_element(ref y);
                    double t = -2;
                    if( _spoil_scenario==4 )
                        t = Double.PositiveInfinity;
                    if( _spoil_scenario==5 )
                        t = Double.NegativeInfinity;
                    double a = -1;
                    if( _spoil_scenario==6 )
                        a = Double.NaN;
                    if( _spoil_scenario==7 )
                        a = Double.PositiveInfinity;
                    if( _spoil_scenario==8 )
                        a = Double.NegativeInfinity;
                    double b = +1;
                    if( _spoil_scenario==9 )
                        b = Double.NaN;
                    if( _spoil_scenario==10 )
                        b = Double.PositiveInfinity;
                    if( _spoil_scenario==11 )
                        b = Double.NegativeInfinity;
                    alglib.barycentricinterpolant p;
                    double v;
                    alglib.polynomialbuildcheb2(a, b, y, 3, out p);
                    v = alglib.barycentriccalc(p, t);
                    _TestResult = _TestResult && doc_test_real(v, 6.0, 0.00005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "polint_t_4");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST polint_t_5
            //      Polynomial interpolation, full list of parameters.
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<6; _spoil_scenario++)
            {
                try
                {
                    double[] y = new double[]{0,0,2};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref y, Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref y, Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref y, Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_vector_by_deleting_element(ref y);
                    double t = -1;
                    if( _spoil_scenario==4 )
                        t = Double.PositiveInfinity;
                    if( _spoil_scenario==5 )
                        t = Double.NegativeInfinity;
                    double v;
                    v = alglib.polynomialcalceqdist(0.0, 2.0, y, 3, t);
                    _TestResult = _TestResult && doc_test_real(v, 2.0, 0.00005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "polint_t_5");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST polint_t_6
            //      Polynomial interpolation, full list of parameters.
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<12; _spoil_scenario++)
            {
                try
                {
                    double[] y = new double[]{-0.116025,0.000000,1.616025};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref y, Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref y, Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref y, Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_vector_by_deleting_element(ref y);
                    double t = -1;
                    if( _spoil_scenario==4 )
                        t = Double.PositiveInfinity;
                    if( _spoil_scenario==5 )
                        t = Double.NegativeInfinity;
                    double a = -1;
                    if( _spoil_scenario==6 )
                        a = Double.NaN;
                    if( _spoil_scenario==7 )
                        a = Double.PositiveInfinity;
                    if( _spoil_scenario==8 )
                        a = Double.NegativeInfinity;
                    double b = +1;
                    if( _spoil_scenario==9 )
                        b = Double.NaN;
                    if( _spoil_scenario==10 )
                        b = Double.PositiveInfinity;
                    if( _spoil_scenario==11 )
                        b = Double.NegativeInfinity;
                    double v;
                    v = alglib.polynomialcalccheb1(a, b, y, 3, t);
                    _TestResult = _TestResult && doc_test_real(v, 2.0, 0.00005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "polint_t_6");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST polint_t_7
            //      Polynomial interpolation, full list of parameters.
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<12; _spoil_scenario++)
            {
                try
                {
                    double[] y = new double[]{0,0,2};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref y, Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref y, Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref y, Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_vector_by_deleting_element(ref y);
                    double t = -2;
                    if( _spoil_scenario==4 )
                        t = Double.PositiveInfinity;
                    if( _spoil_scenario==5 )
                        t = Double.NegativeInfinity;
                    double a = -1;
                    if( _spoil_scenario==6 )
                        a = Double.NaN;
                    if( _spoil_scenario==7 )
                        a = Double.PositiveInfinity;
                    if( _spoil_scenario==8 )
                        a = Double.NegativeInfinity;
                    double b = +1;
                    if( _spoil_scenario==9 )
                        b = Double.NaN;
                    if( _spoil_scenario==10 )
                        b = Double.PositiveInfinity;
                    if( _spoil_scenario==11 )
                        b = Double.NegativeInfinity;
                    double v;
                    v = alglib.polynomialcalccheb2(a, b, y, 3, t);
                    _TestResult = _TestResult && doc_test_real(v, 6.0, 0.00005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "polint_t_7");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST polint_t_8
            //      Polynomial fitting, full list of parameters.
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<10; _spoil_scenario++)
            {
                try
                {
                    double[] x = new double[]{0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_vector_by_deleting_element(ref x);
                    double[] y = new double[]{0.00,0.05,0.26,0.32,0.33,0.43,0.60,0.60,0.77,0.98,1.02};
                    if( _spoil_scenario==4 )
                        spoil_vector_by_value(ref y, Double.NaN);
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref y, Double.PositiveInfinity);
                    if( _spoil_scenario==6 )
                        spoil_vector_by_value(ref y, Double.NegativeInfinity);
                    if( _spoil_scenario==7 )
                        spoil_vector_by_deleting_element(ref y);
                    int m = 2;
                    double t = 2;
                    if( _spoil_scenario==8 )
                        t = Double.PositiveInfinity;
                    if( _spoil_scenario==9 )
                        t = Double.NegativeInfinity;
                    int info;
                    alglib.barycentricinterpolant p;
                    alglib.polynomialfitreport rep;
                    double v;
                    alglib.polynomialfit(x, y, 11, m, out info, out p, out rep);
                    v = alglib.barycentriccalc(p, t);
                    _TestResult = _TestResult && doc_test_real(v, 2.011, 0.002);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "polint_t_8");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST polint_t_9
            //      Polynomial fitting, full list of parameters.
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<14; _spoil_scenario++)
            {
                try
                {
                    double[] x = new double[]{0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_vector_by_deleting_element(ref x);
                    double[] y = new double[]{0.00,0.05,0.26,0.32,0.33,0.43,0.60,0.60,0.77,0.98,1.02};
                    if( _spoil_scenario==4 )
                        spoil_vector_by_value(ref y, Double.NaN);
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref y, Double.PositiveInfinity);
                    if( _spoil_scenario==6 )
                        spoil_vector_by_value(ref y, Double.NegativeInfinity);
                    if( _spoil_scenario==7 )
                        spoil_vector_by_deleting_element(ref y);
                    double[] w = new double[]{1,1.414213562,1,1,1,1,1,1,1,1,1};
                    if( _spoil_scenario==8 )
                        spoil_vector_by_value(ref w, Double.NaN);
                    if( _spoil_scenario==9 )
                        spoil_vector_by_value(ref w, Double.PositiveInfinity);
                    if( _spoil_scenario==10 )
                        spoil_vector_by_value(ref w, Double.NegativeInfinity);
                    if( _spoil_scenario==11 )
                        spoil_vector_by_deleting_element(ref w);
                    double[] xc = new double[]{};
                    double[] yc = new double[]{};
                    int[] dc = new int[]{};
                    int m = 2;
                    double t = 2;
                    if( _spoil_scenario==12 )
                        t = Double.PositiveInfinity;
                    if( _spoil_scenario==13 )
                        t = Double.NegativeInfinity;
                    int info;
                    alglib.barycentricinterpolant p;
                    alglib.polynomialfitreport rep;
                    double v;
                    alglib.polynomialfitwc(x, y, w, 11, xc, yc, dc, 0, m, out info, out p, out rep);
                    v = alglib.barycentriccalc(p, t);
                    _TestResult = _TestResult && doc_test_real(v, 2.023, 0.002);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "polint_t_9");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST polint_t_10
            //      Polynomial fitting, full list of parameters.
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<23; _spoil_scenario++)
            {
                try
                {
                    double[] x = new double[]{1.0,1.0};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_vector_by_deleting_element(ref x);
                    double[] y = new double[]{0.9,1.1};
                    if( _spoil_scenario==4 )
                        spoil_vector_by_value(ref y, Double.NaN);
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref y, Double.PositiveInfinity);
                    if( _spoil_scenario==6 )
                        spoil_vector_by_value(ref y, Double.NegativeInfinity);
                    if( _spoil_scenario==7 )
                        spoil_vector_by_deleting_element(ref y);
                    double[] w = new double[]{1,1};
                    if( _spoil_scenario==8 )
                        spoil_vector_by_value(ref w, Double.NaN);
                    if( _spoil_scenario==9 )
                        spoil_vector_by_value(ref w, Double.PositiveInfinity);
                    if( _spoil_scenario==10 )
                        spoil_vector_by_value(ref w, Double.NegativeInfinity);
                    if( _spoil_scenario==11 )
                        spoil_vector_by_deleting_element(ref w);
                    double[] xc = new double[]{0};
                    if( _spoil_scenario==12 )
                        spoil_vector_by_value(ref xc, Double.NaN);
                    if( _spoil_scenario==13 )
                        spoil_vector_by_value(ref xc, Double.PositiveInfinity);
                    if( _spoil_scenario==14 )
                        spoil_vector_by_value(ref xc, Double.NegativeInfinity);
                    if( _spoil_scenario==15 )
                        spoil_vector_by_deleting_element(ref xc);
                    double[] yc = new double[]{0};
                    if( _spoil_scenario==16 )
                        spoil_vector_by_value(ref yc, Double.NaN);
                    if( _spoil_scenario==17 )
                        spoil_vector_by_value(ref yc, Double.PositiveInfinity);
                    if( _spoil_scenario==18 )
                        spoil_vector_by_value(ref yc, Double.NegativeInfinity);
                    if( _spoil_scenario==19 )
                        spoil_vector_by_deleting_element(ref yc);
                    int[] dc = new int[]{0};
                    if( _spoil_scenario==20 )
                        spoil_vector_by_deleting_element(ref dc);
                    int m = 2;
                    double t = 2;
                    if( _spoil_scenario==21 )
                        t = Double.PositiveInfinity;
                    if( _spoil_scenario==22 )
                        t = Double.NegativeInfinity;
                    int info;
                    alglib.barycentricinterpolant p;
                    alglib.polynomialfitreport rep;
                    double v;
                    alglib.polynomialfitwc(x, y, w, 2, xc, yc, dc, 1, m, out info, out p, out rep);
                    v = alglib.barycentriccalc(p, t);
                    _TestResult = _TestResult && doc_test_real(v, 2.000, 0.001);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "polint_t_10");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST nneighbor_d_1
            //      Nearest neighbor search, KNN queries
            //
            System.Console.WriteLine("50/66");
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
            {
                try
                {
                    double[,] a = new double[,]{{0,0},{0,1},{1,0},{1,1}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref a, Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref a, Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref a, Double.NegativeInfinity);
                    int nx = 2;
                    int ny = 0;
                    int normtype = 2;
                    alglib.kdtree kdt;
                    double[] x;
                    double[,] r = new double[,]{{}};
                    int k;
                    alglib.kdtreebuild(a, nx, ny, normtype, out kdt);
                    x = new double[]{-1,0};
                    k = alglib.kdtreequeryknn(kdt, x, 1);
                    _TestResult = _TestResult && doc_test_int(k, 1);
                    alglib.kdtreequeryresultsx(kdt, ref r);
                    _TestResult = _TestResult && doc_test_real_matrix(r, new double[,]{{0,0}}, 0.05);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "nneighbor_d_1");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST nneighbor_t_2
            //      Subsequent queries; buffered functions must use previously allocated storage (if large enough), so buffer may contain some info from previous call
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
            {
                try
                {
                    double[,] a = new double[,]{{0,0},{0,1},{1,0},{1,1}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref a, Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref a, Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref a, Double.NegativeInfinity);
                    int nx = 2;
                    int ny = 0;
                    int normtype = 2;
                    alglib.kdtree kdt;
                    double[] x;
                    double[,] rx = new double[,]{{}};
                    int k;
                    alglib.kdtreebuild(a, nx, ny, normtype, out kdt);
                    x = new double[]{+2,0};
                    k = alglib.kdtreequeryknn(kdt, x, 2, true);
                    _TestResult = _TestResult && doc_test_int(k, 2);
                    alglib.kdtreequeryresultsx(kdt, ref rx);
                    _TestResult = _TestResult && doc_test_real_matrix(rx, new double[,]{{1,0},{1,1}}, 0.05);
                    x = new double[]{-2,0};
                    k = alglib.kdtreequeryknn(kdt, x, 1, true);
                    _TestResult = _TestResult && doc_test_int(k, 1);
                    alglib.kdtreequeryresultsx(kdt, ref rx);
                    _TestResult = _TestResult && doc_test_real_matrix(rx, new double[,]{{0,0},{1,1}}, 0.05);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "nneighbor_t_2");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST matdet_d_1
            //      Determinant calculation, real matrix, short form
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<7; _spoil_scenario++)
            {
                try
                {
                    double[,] b = new double[,]{{1,2},{2,1}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref b, Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref b, Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref b, Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_matrix_by_adding_row(ref b, 0);
                    if( _spoil_scenario==4 )
                        spoil_matrix_by_adding_col(ref b, 0);
                    if( _spoil_scenario==5 )
                        spoil_matrix_by_deleting_row(ref b);
                    if( _spoil_scenario==6 )
                        spoil_matrix_by_deleting_col(ref b);
                    double a;
                    a = alglib.rmatrixdet(b);
                    _TestResult = _TestResult && doc_test_real(a, -3, 0.0001);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "matdet_d_1");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST matdet_d_2
            //      Determinant calculation, real matrix, full form
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<5; _spoil_scenario++)
            {
                try
                {
                    double[,] b = new double[,]{{5,4},{4,5}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref b, Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref b, Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref b, Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_matrix_by_deleting_row(ref b);
                    if( _spoil_scenario==4 )
                        spoil_matrix_by_deleting_col(ref b);
                    double a;
                    a = alglib.rmatrixdet(b, 2);
                    _TestResult = _TestResult && doc_test_real(a, 9, 0.0001);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "matdet_d_2");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST matdet_d_3
            //      Determinant calculation, complex matrix, short form
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<7; _spoil_scenario++)
            {
                try
                {
                    alglib.complex[,] b = new alglib.complex[,]{{new alglib.complex(1,+1),2},{2,new alglib.complex(1,-1)}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref b, Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref b, Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref b, Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_matrix_by_adding_row(ref b, 0);
                    if( _spoil_scenario==4 )
                        spoil_matrix_by_adding_col(ref b, 0);
                    if( _spoil_scenario==5 )
                        spoil_matrix_by_deleting_row(ref b);
                    if( _spoil_scenario==6 )
                        spoil_matrix_by_deleting_col(ref b);
                    alglib.complex a;
                    a = alglib.cmatrixdet(b);
                    _TestResult = _TestResult && doc_test_complex(a, -2, 0.0001);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "matdet_d_3");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST matdet_d_4
            //      Determinant calculation, complex matrix, full form
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<5; _spoil_scenario++)
            {
                try
                {
                    alglib.complex a;
                    alglib.complex[,] b = new alglib.complex[,]{{new alglib.complex(0,5),4},{new alglib.complex(0,4),5}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref b, Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref b, Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref b, Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_matrix_by_deleting_row(ref b);
                    if( _spoil_scenario==4 )
                        spoil_matrix_by_deleting_col(ref b);
                    a = alglib.cmatrixdet(b, 2);
                    _TestResult = _TestResult && doc_test_complex(a, new alglib.complex(0,9), 0.0001);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "matdet_d_4");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST matdet_d_5
            //      Determinant calculation, complex matrix with zero imaginary part, short form
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<7; _spoil_scenario++)
            {
                try
                {
                    alglib.complex a;
                    alglib.complex[,] b = new alglib.complex[,]{{9,1},{2,1}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref b, Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref b, Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref b, Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_matrix_by_adding_row(ref b, 0);
                    if( _spoil_scenario==4 )
                        spoil_matrix_by_adding_col(ref b, 0);
                    if( _spoil_scenario==5 )
                        spoil_matrix_by_deleting_row(ref b);
                    if( _spoil_scenario==6 )
                        spoil_matrix_by_deleting_col(ref b);
                    a = alglib.cmatrixdet(b);
                    _TestResult = _TestResult && doc_test_complex(a, 7, 0.0001);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "matdet_d_5");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST matdet_t_0
            //      Determinant calculation, real matrix, full form
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<5; _spoil_scenario++)
            {
                try
                {
                    double a;
                    double[,] b = new double[,]{{3,4},{-4,3}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref b, Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref b, Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref b, Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_matrix_by_deleting_row(ref b);
                    if( _spoil_scenario==4 )
                        spoil_matrix_by_deleting_col(ref b);
                    a = alglib.rmatrixdet(b, 2);
                    _TestResult = _TestResult && doc_test_real(a, 25, 0.0001);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "matdet_t_0");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST matdet_t_1
            //      Determinant calculation, real matrix, LU, short form
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<9; _spoil_scenario++)
            {
                try
                {
                    double a;
                    double[,] b = new double[,]{{1,2},{2,5}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref b, Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref b, Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref b, Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_matrix_by_adding_row(ref b, 0);
                    if( _spoil_scenario==4 )
                        spoil_matrix_by_adding_col(ref b, 0);
                    if( _spoil_scenario==5 )
                        spoil_matrix_by_deleting_row(ref b);
                    if( _spoil_scenario==6 )
                        spoil_matrix_by_deleting_col(ref b);
                    int[] p = new int[]{1,1};
                    if( _spoil_scenario==7 )
                        spoil_vector_by_adding_element(ref p, 0);
                    if( _spoil_scenario==8 )
                        spoil_vector_by_deleting_element(ref p);
                    a = alglib.rmatrixludet(b, p);
                    _TestResult = _TestResult && doc_test_real(a, -5, 0.0001);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "matdet_t_1");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST matdet_t_2
            //      Determinant calculation, real matrix, LU, full form
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<6; _spoil_scenario++)
            {
                try
                {
                    double a;
                    double[,] b = new double[,]{{5,4},{4,5}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref b, Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref b, Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref b, Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_matrix_by_deleting_row(ref b);
                    if( _spoil_scenario==4 )
                        spoil_matrix_by_deleting_col(ref b);
                    int[] p = new int[]{0,1};
                    if( _spoil_scenario==5 )
                        spoil_vector_by_deleting_element(ref p);
                    a = alglib.rmatrixludet(b, p, 2);
                    _TestResult = _TestResult && doc_test_real(a, 25, 0.0001);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "matdet_t_2");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST matdet_t_3
            //      Determinant calculation, complex matrix, full form
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<5; _spoil_scenario++)
            {
                try
                {
                    alglib.complex a;
                    alglib.complex[,] b = new alglib.complex[,]{{new alglib.complex(0,5),4},{-4,new alglib.complex(0,5)}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref b, Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref b, Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref b, Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_matrix_by_deleting_row(ref b);
                    if( _spoil_scenario==4 )
                        spoil_matrix_by_deleting_col(ref b);
                    a = alglib.cmatrixdet(b, 2);
                    _TestResult = _TestResult && doc_test_complex(a, -9, 0.0001);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "matdet_t_3");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST matdet_t_4
            //      Determinant calculation, complex matrix, LU, short form
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<9; _spoil_scenario++)
            {
                try
                {
                    alglib.complex a;
                    alglib.complex[,] b = new alglib.complex[,]{{1,2},{2,new alglib.complex(0,5)}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref b, Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref b, Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref b, Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_matrix_by_adding_row(ref b, 0);
                    if( _spoil_scenario==4 )
                        spoil_matrix_by_adding_col(ref b, 0);
                    if( _spoil_scenario==5 )
                        spoil_matrix_by_deleting_row(ref b);
                    if( _spoil_scenario==6 )
                        spoil_matrix_by_deleting_col(ref b);
                    int[] p = new int[]{1,1};
                    if( _spoil_scenario==7 )
                        spoil_vector_by_adding_element(ref p, 0);
                    if( _spoil_scenario==8 )
                        spoil_vector_by_deleting_element(ref p);
                    a = alglib.cmatrixludet(b, p);
                    _TestResult = _TestResult && doc_test_complex(a, new alglib.complex(0,-5), 0.0001);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "matdet_t_4");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST matdet_t_5
            //      Determinant calculation, complex matrix, LU, full form
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<6; _spoil_scenario++)
            {
                try
                {
                    alglib.complex a;
                    alglib.complex[,] b = new alglib.complex[,]{{5,new alglib.complex(0,4)},{4,5}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref b, Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref b, Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref b, Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_matrix_by_deleting_row(ref b);
                    if( _spoil_scenario==4 )
                        spoil_matrix_by_deleting_col(ref b);
                    int[] p = new int[]{0,1};
                    if( _spoil_scenario==5 )
                        spoil_vector_by_deleting_element(ref p);
                    a = alglib.cmatrixludet(b, p, 2);
                    _TestResult = _TestResult && doc_test_complex(a, 25, 0.0001);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "matdet_t_5");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST mincg_d_1
            //      Nonlinear optimization, f(x,y) = 100*(x+3)^4+(y-3)^4
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<12; _spoil_scenario++)
            {
                try
                {
                    //
                    // This example demonstrates minimization of f(x,y) = 100*(x+3)^4+(y-3)^4
                    // with nonlinear conjugate gradient method.
                    //
                    double[] x = new double[]{0,0};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, Double.NegativeInfinity);
                    double epsg = 0.0000000001;
                    if( _spoil_scenario==3 )
                        epsg = Double.NaN;
                    if( _spoil_scenario==4 )
                        epsg = Double.PositiveInfinity;
                    if( _spoil_scenario==5 )
                        epsg = Double.NegativeInfinity;
                    double epsf = 0;
                    if( _spoil_scenario==6 )
                        epsf = Double.NaN;
                    if( _spoil_scenario==7 )
                        epsf = Double.PositiveInfinity;
                    if( _spoil_scenario==8 )
                        epsf = Double.NegativeInfinity;
                    double epsx = 0;
                    if( _spoil_scenario==9 )
                        epsx = Double.NaN;
                    if( _spoil_scenario==10 )
                        epsx = Double.PositiveInfinity;
                    if( _spoil_scenario==11 )
                        epsx = Double.NegativeInfinity;
                    int maxits = 0;
                    alglib.mincgstate state;
                    alglib.mincgreport rep;

                    alglib.mincgcreate(x, out state);
                    alglib.mincgsetcond(state, epsg, epsf, epsx, maxits);
                    alglib.mincgoptimize(state, function1_grad, null, null);
                    alglib.mincgresults(state, out x, out rep);

                    _TestResult = _TestResult && doc_test_int(rep.terminationtype, 4);
                    _TestResult = _TestResult && doc_test_real_vector(x, new double[]{-3,3}, 0.005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "mincg_d_1");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST mincg_d_2
            //      Nonlinear optimization with additional settings and restarts, F(x,y) = 100*(x+3)^4+(y-3)^4
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<18; _spoil_scenario++)
            {
                try
                {
                    //
                    // This example demonstrates minimization of f(x,y) = 100*(x+3)^4+(y-3)^4
                    // with nonlinear conjugate gradient method.
                    //
                    // Several advanced techniques are demonstrated:
                    // * upper limit on step size
                    // * restart from new point
                    //
                    double[] x = new double[]{0,0};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, Double.NegativeInfinity);
                    double epsg = 0.0000000001;
                    if( _spoil_scenario==3 )
                        epsg = Double.NaN;
                    if( _spoil_scenario==4 )
                        epsg = Double.PositiveInfinity;
                    if( _spoil_scenario==5 )
                        epsg = Double.NegativeInfinity;
                    double epsf = 0;
                    if( _spoil_scenario==6 )
                        epsf = Double.NaN;
                    if( _spoil_scenario==7 )
                        epsf = Double.PositiveInfinity;
                    if( _spoil_scenario==8 )
                        epsf = Double.NegativeInfinity;
                    double epsx = 0;
                    if( _spoil_scenario==9 )
                        epsx = Double.NaN;
                    if( _spoil_scenario==10 )
                        epsx = Double.PositiveInfinity;
                    if( _spoil_scenario==11 )
                        epsx = Double.NegativeInfinity;
                    double stpmax = 0.1;
                    if( _spoil_scenario==12 )
                        stpmax = Double.NaN;
                    if( _spoil_scenario==13 )
                        stpmax = Double.PositiveInfinity;
                    if( _spoil_scenario==14 )
                        stpmax = Double.NegativeInfinity;
                    int maxits = 0;
                    alglib.mincgstate state;
                    alglib.mincgreport rep;

                    // first run
                    alglib.mincgcreate(x, out state);
                    alglib.mincgsetcond(state, epsg, epsf, epsx, maxits);
                    alglib.mincgsetstpmax(state, stpmax);
                    alglib.mincgoptimize(state, function1_grad, null, null);
                    alglib.mincgresults(state, out x, out rep);

                    _TestResult = _TestResult && doc_test_real_vector(x, new double[]{-3,3}, 0.005);

                    // second run - algorithm is restarted with mincgrestartfrom()
                    x = new double[]{10,10};
                    if( _spoil_scenario==15 )
                        spoil_vector_by_value(ref x, Double.NaN);
                    if( _spoil_scenario==16 )
                        spoil_vector_by_value(ref x, Double.PositiveInfinity);
                    if( _spoil_scenario==17 )
                        spoil_vector_by_value(ref x, Double.NegativeInfinity);
                    alglib.mincgrestartfrom(state, x);
                    alglib.mincgoptimize(state, function1_grad, null, null);
                    alglib.mincgresults(state, out x, out rep);

                    _TestResult = _TestResult && doc_test_int(rep.terminationtype, 4);
                    _TestResult = _TestResult && doc_test_real_vector(x, new double[]{-3,3}, 0.005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "mincg_d_2");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST minasa_d_1
            //      Nonlinear optimization, F(x,y) = 100*(x+3)^4+(y-3)^4, bound constraints
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<24; _spoil_scenario++)
            {
                try
                {
                    //
                    // This example demonstrates minimization of f(x,y) = 100*(x+3)^4+(y-3)^4
                    // subject to bound constraints -1<=x0<=+1, -1<=x1<=+1, using ASA.
                    //
                    double[] x = new double[]{0,0};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_vector_by_adding_element(ref x, 0);
                    if( _spoil_scenario==4 )
                        spoil_vector_by_deleting_element(ref x);
                    double[] bndl = new double[]{-1,-1};
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref bndl, Double.NaN);
                    if( _spoil_scenario==6 )
                        spoil_vector_by_value(ref bndl, Double.PositiveInfinity);
                    if( _spoil_scenario==7 )
                        spoil_vector_by_value(ref bndl, Double.NegativeInfinity);
                    if( _spoil_scenario==8 )
                        spoil_vector_by_adding_element(ref bndl, 0);
                    if( _spoil_scenario==9 )
                        spoil_vector_by_deleting_element(ref bndl);
                    double[] bndu = new double[]{+1,+1};
                    if( _spoil_scenario==10 )
                        spoil_vector_by_value(ref bndu, Double.NaN);
                    if( _spoil_scenario==11 )
                        spoil_vector_by_value(ref bndu, Double.PositiveInfinity);
                    if( _spoil_scenario==12 )
                        spoil_vector_by_value(ref bndu, Double.NegativeInfinity);
                    if( _spoil_scenario==13 )
                        spoil_vector_by_adding_element(ref bndu, 0);
                    if( _spoil_scenario==14 )
                        spoil_vector_by_deleting_element(ref bndu);
                    double epsg = 0.0000000001;
                    if( _spoil_scenario==15 )
                        epsg = Double.NaN;
                    if( _spoil_scenario==16 )
                        epsg = Double.PositiveInfinity;
                    if( _spoil_scenario==17 )
                        epsg = Double.NegativeInfinity;
                    double epsf = 0;
                    if( _spoil_scenario==18 )
                        epsf = Double.NaN;
                    if( _spoil_scenario==19 )
                        epsf = Double.PositiveInfinity;
                    if( _spoil_scenario==20 )
                        epsf = Double.NegativeInfinity;
                    double epsx = 0;
                    if( _spoil_scenario==21 )
                        epsx = Double.NaN;
                    if( _spoil_scenario==22 )
                        epsx = Double.PositiveInfinity;
                    if( _spoil_scenario==23 )
                        epsx = Double.NegativeInfinity;
                    int maxits = 0;
                    alglib.minasastate state;
                    alglib.minasareport rep;

                    alglib.minasacreate(x, bndl, bndu, out state);
                    alglib.minasasetcond(state, epsg, epsf, epsx, maxits);
                    alglib.minasaoptimize(state, function1_grad, null, null);
                    alglib.minasaresults(state, out x, out rep);

                    _TestResult = _TestResult && doc_test_int(rep.terminationtype, 4);
                    _TestResult = _TestResult && doc_test_real_vector(x, new double[]{-1,1}, 0.005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "minasa_d_1");
            _TotalResult = _TotalResult && _TestResult;


            System.Console.WriteLine("66/66");
        }
        catch
        {
            System.Console.WriteLine("Unhandled exception was raised!");
            return 1;
        }
        return _TotalResult ? 0 : 1;
    }
}
