using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
//using System.Linq;
using System.Text;
using System.Windows.Forms;
using System.IO;
using System.Reflection;
using alglib;
//#pragma comment(lib, "libalglib.lib")

namespace MtmOR
{
    public partial class Form1 : Form
    {
        public Form2 Fform2;
        public bool close;

        public Form1()
        {
            InitializeComponent();
            Fform2 = new Form2();
            Fform2.Owner = this;
            
            //Fform2.Top = Top;
            //System.Globalization.CultureInfo.CurrentCulture.NumberFormat.NumberDecimalSeparator
            //System.Globalization.CultureInfo culture = new System.Globalization.CultureInfo("en-US");
        }


        public class Observation
        {
            public string Time;
            public double LimbA;
            public double LimbD;
            public string Object;
            public string File;
            public int Orientation;
            public int Error1;
            public int Error2;
            public string Pos;
            public string Date;
            public string Dir;
            public Observation()
            {
                Time = "";
                LimbA = 500;
                LimbD = 500;
                Object = "";
                File = "";
                Orientation = -1;
                Error1 = 1;
                Error2 = 1;
                Pos = "";
                Date = "";

            }
            public void AddtoObs(ListBox lb)
            {
                string endstring = ""; ;
                string ttime = Time.Replace(':', ' ');
                if (ttime.IndexOf(' ') == 1) ttime = "0" + ttime;
                string tdate = Date.Replace('-', ' ');
                endstring += tdate;
                endstring += "  ";
                endstring += ttime;
                endstring += "  ";
                endstring += Orientation.ToString();
                endstring += "  ";

                //now RA
                int pos1 = Pos.IndexOf(',');
                string stA = Pos.Substring(0, pos1);
                string stD = Pos.Substring(pos1 + 8, Pos.Length - pos1 - 8);
                if (stA[2] == ' ') stA = "0" + stA;
                for (int i = 0; i < stA.Length; i++)
                {
                    if (((stA[i] > 57) || (stA[i] < 48)) && (stA[i] != '.') && (stA[i] != '+') && (stA[i] != '-'))
                    {
                        stA = stA.Remove(i, 1);
                        //	stA=stA.Replace(stA[i],' ');
                    }
                }
                for (int i = 0; i < stD.Length; i++)
                {
                    if (((stD[i] > 57) || (stD[i] < 48)) && (stD[i] != '.') && (stD[i] != '+') && (stD[i] != '-'))
                    {
                        stD = stD.Remove(i, 1);
                        //	stD=stD.Replace(stD[i],' ');
                    }
                }
                if (stD[2] == ' ')
                {
                    string tst = stD[0] + "0" + stD[1];
                    stD = tst + stD.Substring(2, stD.Length - 2);

                }
                endstring += stA;
                endstring += "  ";
                endstring += stD;
                endstring += "  ";
                //LimbT, LimbD
                double grad = Math.Floor(LimbA);
                string stgrad = grad.ToString();
                double min = Math.Floor((LimbA - grad) * 60);
                string stmin = min.ToString();
                double sec = Math.Floor(((LimbA - grad) * 60 - min) * 60);
                string stsec = sec.ToString();
                stgrad = stgrad.PadRight(3);
                stmin = stmin.PadRight(2);
                stsec = stsec.PadRight(2);
                endstring = endstring + stgrad + " " + stmin + " " + stsec + "  ";
                grad = Math.Floor(LimbD);
                stgrad = grad.ToString();
                min = Math.Floor((LimbD - grad) * 60);
                stmin = min.ToString();
                sec = Math.Floor(((LimbD - grad) * 60 - min) * 60);
                stsec = sec.ToString();
                stgrad = stgrad.PadRight(3);
                stmin = stmin.PadRight(2);
                stsec = stsec.PadRight(2);
                endstring = endstring + stgrad + " " + stmin + " " + stsec + "  ";
                int[] coordint= ProclogCoord(File);
                string coord = coordint[0].ToString()+" "+coordint[1].ToString()+"  ";
                
                //    endstring = endstring + "512 512  ";
                
                
                if (coord != "512 512  ")  //512 512
                {
                    //coord = "512 512  ";
                    endstring += coord;
                    endstring += File;
                    if ((Math.Abs(coordint[0] - 512) < 150) && (Math.Abs(coordint[1] - 512) < 150))
                    {
                        lb.Items.Add(endstring);
                    }
                }
            }

            public int[] ProclogCoord(string filename)
            {

                System.Globalization.CultureInfo ci = new System.Globalization.CultureInfo("ru-RU");

                string fullname = Dir+filename + ".proclog";

                if (System.IO.File.Exists(fullname))
                {
                    StreamReader sr = new StreamReader(fullname);
                    int[] coordint = new int[2];
                    coordint[0] = 512;   //512
                    coordint[1] = 512;   //512
                    string str = sr.ReadToEnd();
                    sr.Dispose();
                    int num1 = str.IndexOf(" * ");
                    if (num1 == -1)
                        return coordint;
                    int num2 = str.IndexOf(" * ", num1 + 2);
                    if (num2 == -1)
                        return coordint;

                    int num3 = str.IndexOf("Object: ");
                    if (num3 == -1)
                        return coordint;
                    int num4=str.IndexOf("Catalog",num3);
                    if (num4==-1)
                        return coordint;

                    string line3 = str.Substring(num3 + 7, num4 - num3 - 7);
                    line3=line3.Replace(" ","");
                    line3 = line3.Replace("\r\n", "");
                    if (filename.IndexOf(line3) == -1)
                        return coordint;
                    string line = str.Substring(num2 + 5, 50);
                    line = line.Replace(".", ci.NumberFormat.NumberDecimalSeparator);
                    line = line.Replace(",", ci.NumberFormat.NumberDecimalSeparator);
                    char[] sep = new char[1];
                    sep[0] = ' ';
                    string[] spl = line.Split(sep, System.StringSplitOptions.RemoveEmptyEntries);
                    coordint[0] = (int)Convert.ToDouble(spl[0]);
                    coordint[1] = (int)Convert.ToDouble(spl[1]);

                    
                    return coordint;
                }

                else
                {
                    int[] coordint = new int[2];
                    coordint[0] = 512;
                    coordint[1] = 512;
                    return coordint;
                }
            }
        }




         string LoggerToObs(string loggerfile)
        {
            progressBar1.Value = 0;
            lObserv.Items.Clear();



            //прописываем основные строки
            string refStr1 = "[CHAOS] - Observing object \"";
            string refStr2 = "[CHAOS] - Pointing started: RA =";
            //string refStr2="[Agate/BKU] - Fine pointing complete";
            string refStr3 = "'AGATE LIMB DEC',";
            string refStr4 = "'AGATE LIMB HA',";
            string refStr5 = "'AGATE LIMB TIME',";
            //string refStr6="[CHAOS] -   RA =";
            string refStr7 = "[CHAOS] -   Orientation:";
            string refStr8 = "[logger] - Exposure complete.";
            string refStr9 = "[logger] - 127.0.0.1> grab(1,'";



            StreamReader sr = new StreamReader(loggerfile);
            
             
            String line;
            Observation Obs = new Observation();
            System.Globalization.CultureInfo ci = new System.Globalization.CultureInfo("ru-RU");
            line = sr.ReadLine();
            int count = 1;
            while (line != null)
            {
                if (line.IndexOf(refStr1) != -1)
                {
                    Obs = new Observation();
                    
                    int pos1 = line.IndexOf("\"", 54);
                    string name = line.Substring(54, pos1 - 54);
                    Obs.Object = name;
                    count = 2;
                }


                if ((line.IndexOf(refStr2) != -1) && (count == 2))
                {
                    int pos2 = line.IndexOf("\"", 59);
                    string name = line.Substring(59, pos2 - 59);
                    Obs.Pos = name;


                    Obs.Error1 = 0;
                    count = 3;
                }


                if ((line.IndexOf(refStr3) != -1) && (count == 3))
                {

                    int pos2 = line.IndexOf(")", 72);
                    string name = line.Substring(72, pos2 - 72);
                    name = name.Replace('.', ',');
                    name = name.Replace(".", ci.NumberFormat.NumberDecimalSeparator);
                    name = name.Replace(",", ci.NumberFormat.NumberDecimalSeparator);
                    Obs.LimbD = Convert.ToDouble(name);
                    count = 4;
                }
                if ((line.IndexOf(refStr4) != -1) && (count == 4))
                {
                    int pos2 = line.IndexOf(")", 71);
                    string name = line.Substring(71, pos2 - 71);
                    name = name.Replace('.', ',');
                    name = name.Replace(".", ci.NumberFormat.NumberDecimalSeparator);
                    name = name.Replace(",", ci.NumberFormat.NumberDecimalSeparator);
                    Obs.LimbA = Convert.ToDouble(name);
                    count = 5;
                }
                if ((line.IndexOf(refStr5) != -1) && (count == 5))
                {

                    int pos2 = line.IndexOf(")", 73);
                    string name = line.Substring(73, pos2 - 74);
                    Obs.Time = name;
                    pos2 = line.IndexOf(' ', 0);
                    name = line.Substring(0, pos2);
                    
                    Obs.Date = name;
                    count = 7;
                }

                /*if ((line.IndexOf(refStr6)!=-1)&&(count==6))
                {
                    int pos2=line.IndexOf("\"",43);
                    string name=line.Substring(43,pos2-43);
                    Obs.Pos=name;
                		
                    count=7;
                }*/

                if ((line.IndexOf(refStr7) != -1) && (count == 7))
                {
                    char nname = line[51];
                    if (nname == '0') Obs.Orientation = 2;
                    else Obs.Orientation = 1;
                    count = 9;
                }
                if ((line.IndexOf(refStr8) != -1) && (count == 8))
                {

                    Obs.Error2 = 0;
                    count = 9;
                }

                if ((line.IndexOf(refStr9) != -1) && (count == 9))
                {

                    int pos2 = line.IndexOf("'", 56);
                    string name = line.Substring(56, pos2 - 56);
                    Obs.File = name;

                    //if ((Obs.Error1==0)&&(Obs.Error2==0))
                    //{
                    Obs.Dir = loggerfile.Substring(0,loggerfile.LastIndexOf('\\'))+'\\';
                    Obs.AddtoObs(lObserv);
                    //}
                    count = 1;
                }

                //Obs.AddtoObs(lOrient);	

                //lOrient.Items.Add(line);

                //progressBar1.Value += line.Length;
                line = sr.ReadLine();
            }
            progressBar1.Value = progressBar1.Maximum;
            sr.Close();
            string fobsnew = "dl_" + loggerfile.Substring(loggerfile.LastIndexOf('\\') + 1, loggerfile.Length - loggerfile.LastIndexOf('\\') - 1) + ".obs";
            StreamWriter sr2 = new StreamWriter(fobsnew);
            for (int i = 0; i < lObserv.Items.Count; i++)
            {
                sr2.WriteLine(lObserv.Items[i].ToString());
            }
            sr2.Close();
            return fobsnew;

        }



        /*
         * void Button2Click(object sender, EventArgs e)
		{
			if (saveFileDialog1.ShowDialog()==DialogResult.OK)
			{
				StreamWriter sr = new StreamWriter(saveFileDialog1.FileName);
				for (int i=0;i<lOrient.Items.Count;i++)
				{
					sr.WriteLine(lOrient.Items[i].ToString());
				}
				sr.Close();
			}
			
		}
		
		void BLoadClick(object sender, EventArgs e)
		{
			DialogResult dr = openFileDialog1.ShowDialog();
            if (dr == DialogResult.OK)
            {
            tFile.Text = openFileDialog1.FileName;
            System.IO.FileInfo fi = new System.IO.FileInfo(tFile.Text);
            progressBar1.Maximum =  (int)fi.Length;
            }
			
			
		}*/


         private void button1_Click(object sender, EventArgs e)
         {
             Calculate(false,"");
         }


         private double PopVert(double t)
         {
             double r = 0;
             if (t < 0)
             {
                 r = -1 / 5.0 * t - 1;
             }
             if (t >= 0)
             {
                 r = -1.2 * t - 1;
             }

             return r / 3.7;
         }

         private double PopHor(double t)
         {
             double r = 0;
             if (t < 0)
             {
                 r = -t;
             }
             if ((t >= 0) && (t <= 5))
             {
                 r = 0;
             }
             if (t > 5)
             {
                 r = 3 / 5.0 * t - 3;
             }
             return r / 3.7;
         }

         private double PopVert2(double x, double y)
         {
             double alpha_gr = Convert.ToDouble(tangle.Text);
             double alpha = alpha_gr / 180.0 * Math.PI;
             return y * Math.Sin(alpha) + x * Math.Cos(alpha);
         }

         private double PopHor2(double x, double y)
         {
             double alpha_gr = Convert.ToDouble(tangle.Text);
             double alpha = alpha_gr / 180.0 * Math.PI;
             return x * Math.Sin(alpha) + y * Math.Cos(alpha);
         }
        
        private void Calculate (bool recalc, string filename)
        {
                        
            
            //FILE fin = fopen("obs.dat", "r");
            
            //double Fi=43.7486111;   //{иЁа®в  ¬Ґбв  ­ Ў«о¤Ґ­Ёп ў Ја ¤гб е}
            //double Lamb=2.844444;  //{¤®«Ј®в  ¬Ґбв  ­ Ў«о¤Ґ­Ёп ў з б е}
            //double Fi = 59.4612;
            //double Lamb = 2 + 0.021944444444444444444444444443778;


            if (cbMethod.SelectedIndex == -1)
            {
                MessageBox.Show("Не выбран метод решения");
                return;
            }
            
            System.Globalization.CultureInfo ci = new System.Globalization.CultureInfo("ru-RU");
            int method = cbMethod.SelectedIndex;

            
            lObserv.Items.Clear();
            lOper2.Items.Clear();


            string ftscope = "e:\\1\\dpas\\MTM500M.dat";
            ftscope = tfiletscope.Text;
            if (ftscope == "")
            {
                MessageBox.Show("Выберите файл с настройками телескопа");
                return;
            }
            StreamReader finmtm = new StreamReader(ftscope);
            char[] separator = new char[] { (char)09, (char)32, '/', ':' };
            string line;
            string[] pars;
            line = finmtm.ReadLine();
            pars = line.Split(separator, StringSplitOptions.RemoveEmptyEntries);
            int Dt = Convert.ToInt32(pars[0]);
            line = finmtm.ReadLine();
            pars = line.Split(separator, StringSplitOptions.RemoveEmptyEntries);
            int JD_R = Convert.ToInt32(pars[0]);
            line = finmtm.ReadLine();
            line = line.Replace(".", ci.NumberFormat.NumberDecimalSeparator);
            line = line.Replace(",", ci.NumberFormat.NumberDecimalSeparator);
            pars = line.Split(separator, StringSplitOptions.RemoveEmptyEntries);
            double Temp = Convert.ToDouble(pars[0]);
            line = finmtm.ReadLine();
            line = line.Replace(".", ci.NumberFormat.NumberDecimalSeparator);
            line = line.Replace(",", ci.NumberFormat.NumberDecimalSeparator);
            pars = line.Split(separator, StringSplitOptions.RemoveEmptyEntries);
            double Dav = Convert.ToDouble(pars[0]);
            /*line = finmtm.ReadLine();
            pars = line.Split(separator, StringSplitOptions.RemoveEmptyEntries);
            int God0 = Convert.ToInt32(pars[0]);
            int Mes0 = Convert.ToInt32(pars[1]);
            int Den0 = Convert.ToInt32(pars[2]);
            int Hour0 = Convert.ToInt32(pars[3]);
            int Min0 = Convert.ToInt32(pars[4]);
            line = finmtm.ReadLine();
            line = line.Replace(".", ci.NumberFormat.NumberDecimalSeparator);
            line = line.Replace(",", ci.NumberFormat.NumberDecimalSeparator);
            pars = line.Split(separator, StringSplitOptions.RemoveEmptyEntries);
            double ds0 = Convert.ToDouble(pars[0]);
            double V = Convert.ToDouble(pars[1]);*/
            line = finmtm.ReadLine();
            line = line.Replace(".", ci.NumberFormat.NumberDecimalSeparator);
            line = line.Replace(",", ci.NumberFormat.NumberDecimalSeparator);
            pars = line.Split(separator, StringSplitOptions.RemoveEmptyEntries);
            double Fi = Convert.ToDouble(pars[0]);
            line = finmtm.ReadLine();
            line = line.Replace(".", ci.NumberFormat.NumberDecimalSeparator);
            line = line.Replace(",", ci.NumberFormat.NumberDecimalSeparator);
            pars = line.Split(separator, StringSplitOptions.RemoveEmptyEntries);
            double Lamb = Convert.ToDouble(pars[0]);
            line = finmtm.ReadLine();
            line = line.Replace(".", ci.NumberFormat.NumberDecimalSeparator);
            line = line.Replace(",", ci.NumberFormat.NumberDecimalSeparator);
            pars = line.Split(separator, StringSplitOptions.RemoveEmptyEntries);
            double ex_t = Convert.ToDouble(pars[0]);
            double Eex_t = Convert.ToDouble(pars[1]);
            line = finmtm.ReadLine();
            line = line.Replace(".", ci.NumberFormat.NumberDecimalSeparator);
            line = line.Replace(",", ci.NumberFormat.NumberDecimalSeparator);
            pars = line.Split(separator, StringSplitOptions.RemoveEmptyEntries);
            double ex_d = Convert.ToDouble(pars[0]);
            double Eex_d = Convert.ToDouble(pars[1]);



            string fobs;
            int ii2;
            StreamReader fin2;

            if (recalc == false)
            {
                //StreamReader fin = new StreamReader("obs.dat");
                bool loggerfile = false;

                Fform2.lResult.Clear();
                fobs = tfileobs.Text;
                if (fobs == "")
                {
                    MessageBox.Show("Выберите файл с наблюдениями");
                    return;
                }
                fin2 = new StreamReader(fobs);
                line = "";
                ii2 = 0;
                line = fin2.ReadLine();
                if (line.IndexOf("Log") != -1)
                {
                    loggerfile = true;
                }

                if (loggerfile == true)
                {
                    progressBar1.Visible = true;
                    System.IO.FileInfo fi = new System.IO.FileInfo(fobs);
                    progressBar1.Maximum = (int)fi.Length;
                    progressBar1.Value = 0;
                    string fobsnew = LoggerToObs(fobs);
                    fobs = fobsnew;
                    progressBar1.Visible = false;
                    fin2.Close();
                    fin2 = new StreamReader(fobs);
                    line = "";
                    line = fin2.ReadLine();
                }

            }

            else //recalc==true
            {
                fobs = filename;
                fin2 = new StreamReader(fobs);
                line = "";
                ii2 = 0;
                line = fin2.ReadLine();
            }
            
            
            while (line != null)
            {
                
                ii2++;
                line = fin2.ReadLine();
            }
            
            fin2.Close();


            if (recalc == false)
            {
                FileInfo fn = new FileInfo(fobs);
                fobs = fobs + ".temp";
                fn.CopyTo(fobs, true);
            }
            
            int NN = ii2 * 2;
            int MM=0;
            if (method == 0) MM = 8;
            if (method == 1) // пр. восхождение
            {
                MM = 6;
                NN = NN / 2;
            }
            if (method == 2) // сколонение
            {
                MM = 4;
                NN = NN / 2;
            }

            StreamReader fin = new StreamReader(fobs);
            separator = new char[] { (char)09, (char)32, '/', ':' };
            double[] OC=new double [NN];
            double[] y = new double[NN];
            double[,] Fmatrix=new double[NN,MM];//8,36
            double[] c = new double[NN];//36
            lsfit.lsfitreport rep = new lsfit.lsfitreport();
            int info = 0;
            int i_count = 0;
            double Delta=0;
            lOper.Items.Clear();
            dg1.RowCount = MM;
            dg1.ColumnCount = 3;
            if (recalc == false)
            {
                for (int i = 0; i < dg1.RowCount; i++)
                {
                    dg1[0, i].Value = 1;
                }
            }
            else
            {
                int method_old = Convert.ToInt32((lSummary.Items[0].ToString())[7]) - 48;
                if (method != method_old)
                {
                    for (int i = 0; i < dg1.RowCount; i++)
                    {
                        dg1[0, i].Value = 1;
                    }
                }
                   
            }
            lSummary.Items.Clear();
           
            dg1.Columns[1].Width = 225;
            dg2.RowCount = NN;
            dg2.ColumnCount = 9;
            
            dg2.Columns[0].Name = "Number";
            dg2.Columns[1].Name = "Observation";
            dg2.Columns[2].Name = "Calculation";
            dg2.Columns[3].Name = "O-C";
            dg2.Columns[4].Name = "Alpha";
            dg2.Columns[5].Name = "Delta";
            dg2.Columns[6].Name = "T";
            dg2.Columns[7].Name = "Delta";
            dg2.Columns[8].Name = "Alpha";


            //MessageBox.Show(Convert.ToString(((CheckBox)dg1[0,0].Value).Checked));
               // MessageBox.Show(Convert.ToString(dg1[0,0].Value));
            
            int Alf1,Alf2,Delt1,Delt2;
            double Alf3,Delt3;
            
            while ((line = fin.ReadLine()) != null)
            {
                
                line = line.Replace(".", ci.NumberFormat.NumberDecimalSeparator);
                line = line.Replace(",", ci.NumberFormat.NumberDecimalSeparator);
                lOper2.Items.Add(line);
                pars = line.Split(separator, StringSplitOptions.RemoveEmptyEntries);

                //Разбиваешь на составляющие. Дальше работаешь с массивом строк pars
                //textBox1.Text = Convert.ToString(pars.Length);
                int God=Convert.ToInt32(pars[0]);
                int Mes=Convert.ToInt32(pars[1]);
                int Den=Convert.ToInt32(pars[2]);
                int Hour=Convert.ToInt32(pars[3]);
                int Min=Convert.ToInt32(pars[4]);
                //double ss = Convert.ToDouble("11,1");
                pars[5]=pars[5].Replace(".", ci.NumberFormat.NumberDecimalSeparator);
                pars[5]=pars[5].Replace(",", ci.NumberFormat.NumberDecimalSeparator);
                double Sec=Convert.ToDouble(pars[5]);
                int NWE=Convert.ToInt32(pars[6]);
                Alf1=Convert.ToInt32(pars[7]);
                Alf2=Convert.ToInt32(pars[8]);
                pars[9]=pars[9].Replace(".", ci.NumberFormat.NumberDecimalSeparator);
                pars[9]=pars[9].Replace(",", ci.NumberFormat.NumberDecimalSeparator);
                Alf3=Convert.ToDouble(pars[9]);
                Delt1=Convert.ToInt32(pars[10]);
                Delt2=Convert.ToInt32(pars[11]);
                pars[12]=pars[12].Replace(".", ci.NumberFormat.NumberDecimalSeparator);
                pars[12]=pars[12].Replace(",", ci.NumberFormat.NumberDecimalSeparator);
                Delt3=Convert.ToDouble(pars[12]);
                double t1=Convert.ToDouble(pars[13]);
                double t2=Convert.ToDouble(pars[14]);
                double t3=Convert.ToDouble(pars[15]);
                double d1=Convert.ToDouble(pars[16]);
                double d2=Convert.ToDouble(pars[17]);
                double d3=Convert.ToDouble(pars[18]);
                int dtt=Convert.ToInt32(pars[19]);
                int ddd=Convert.ToInt32(pars[20]);

                double tin=t1+t2/60.0+t3/3600.0;
                double din=d1+d2/60.0+d3/3600.0;
                string WE;
                int Znnn = 1;
                if (NWE == 1)
                {
                    WE = "W";
                    Znnn = -1;
                }
                else
                {
                    WE = "E";
                    Znnn = 1;
                }
                int Zn_Delt = 1;
                if (Delt1 < 0)
                {
                    Zn_Delt = -1;
                    Delt1 = Math.Abs(Delt1);
                }
                if ((Delt1 == 0) && (Delt2 < 0))
                {
                    Zn_Delt = -1;
                    Delt2 = Math.Abs(Delt2);
                }
                if ((Delt1 == 0) && (Delt2 == 0) && (Delt3 < 0))
                {
                    Zn_Delt = -1;
                    Delt3 = Math.Abs(Delt3);
                }
                Delta = Zn_Delt * (Delt1 + Delt2 / 60.0 + Delt3 / 3600.0);
//                Delta = (Delt1 + Delt2 / 60.0 + Delt3 / 3600.0);

                Mathgr Mgr=new Mathgr();
                double dt_CCD=0;
                double dd_CCD = 0;
                if (NWE == 1)
                {
                    dt_CCD = -1.205 * ((dtt - 512.0) / Mgr.cos(Delta)) / 15.0;
                    dd_CCD = -1.205 * (ddd - 512.0);
                }
                if (NWE == 2)
                {
                    dt_CCD = 1.205 * ((dtt - 512.0) / Mgr.cos(Delta)) / 15.0;
                    dd_CCD = 1.205 * (ddd - 512.0);
                }

                Alf3 = Alf3 + dt_CCD;
                Delt3 = Delt3 + dd_CCD;
                
                //Delta = Zn_Delt * (Delt1 + Delt2 / 60.0 + Delt3 / 3600.0);
                
                Telescope Tscope = new Telescope();
                Tscope.Fi = Fi;
                Tscope.Lamb = Lamb;
                Tscope.Dt = Dt;
                Tscope.Temp = Temp;
                Tscope.Dav = Dav;
                Tscope.ex_t = ex_t;
                Tscope.ex_d = ex_d;
                Tscope.Eex_t = Eex_t;
                Tscope.Eex_d = Eex_d;
                Tscope.JD_R = JD_R;
                Tscope.Alf1 = Alf1;
                Tscope.Alf2 = Alf2;
                Tscope.Alf3 = Alf3;
                Tscope.WE = WE;
                Tscope.Zn_Delt = Zn_Delt;
                Tscope.Delt1 = Delt1;
                Tscope.Delt2 = Delt2;
                Tscope.Delt3 = Delt3;
                Tscope.God = God;
                Tscope.Mes = Mes;
                Tscope.Den = Den;
                Tscope.Hour = Hour;
                Tscope.Min = Min;
                Tscope.Sec = Sec;
                Tscope.Ustanovka();
                Type type = Tscope.GetType();
                BindingFlags flags = BindingFlags.Public | BindingFlags.Instance;
                PropertyInfo[] properties = type.GetProperties(flags);

                dg2[4, i_count].Value = Convert.ToString(Alf1) + " " + Convert.ToString(Alf2) + " " + Convert.ToString(Alf3);
                dg2[5, i_count].Value = Convert.ToString(Delt1) + " " + Convert.ToString(Delt2) + " " + Convert.ToString(Delt3);
                dg2[6, i_count].Value = Convert.ToString(Tscope.t);
                dg2[7, i_count].Value = Convert.ToString(Delta);
                double Alff=(Alf1 + Alf2 / 60.0 + Alf3 / 3600.0);
                dg2[8, i_count].Value = Convert.ToString(Alff);
                

               /* foreach (PropertyInfo property in properties)
                {

                    lOper.Items.Add("1");
                    lOper.Items.Add("Name: " + property.Name + ", Value: " + property.GetValue(Tscope, null));
                    //Console.WriteLine();
                }
                */
                
                //Tscope.Ustanovka(Fi,Lamb,Dt,WE,Temp,Dav,ex_t,Eex_t,ex_d,Eex_d,God0,Mes0,Den0,Hour0,Min0,ds0,V,JD_R,Alf1,Alf2,Alf3,Zn_Delt,Delt1,Delt2,Delt3,God,Mes,Den,Hour,Min,Sec);
                


                /*public static void lsfitlinear(ref double[] y,
    ref double[,] fmatrix,
    int n,
    int m,
    ref int info,
    ref double[] c,
    ref lsfitreport rep)*/
                if (method == 0)
                {

                    dg1[1, 0].Value = "Нуль-пункт по склонению";
                    dg1[1, 1].Value = "Установки пол.оси по азимуту";
                    dg1[1, 2].Value = "Установки пол.оси по широте (склонению)";
                    dg1[1, 3].Value = "Нуль-пункт по часовому углу";
                    dg1[1, 4].Value = "Наклонности";
                    dg1[1, 5].Value = "Коллимации";
                    dg1[1, 6].Value = "Гнутия по склонению";
                    dg1[1, 7].Value = "Гнутия по прямому восхождению";

                    Fmatrix[i_count, 0] = -1.0;
                    Fmatrix[i_count, 1] = Znnn * Mgr.cos(Fi) * Mgr.cos(15.0 * Tscope.t);
                    Fmatrix[i_count, 2] = Znnn * Mgr.sin(15.0 * Tscope.t);
                    Fmatrix[i_count, 3] = 0;
                    Fmatrix[i_count, 4] = 0;
                    Fmatrix[i_count, 5] = 0;
                    Fmatrix[i_count, 6] = Znnn * (-Mgr.sin(Fi) * Mgr.cos(Delta) - Mgr.cos(Fi) * Mgr.sin(Delta) * Mgr.cos(15.0 * Tscope.t));
                    Fmatrix[i_count, 7] = 0;

                    
                    // Fmatrix[i_count, 8] = 0;
                    // Fmatrix[i_count, 9] = 0;

                    Fmatrix[i_count + 1, 0] = 0;
                    Fmatrix[i_count + 1, 1] = Mgr.cos(Fi) * Mgr.sin(15.0 * Tscope.t) * Mgr.tn(Delta);
                    Fmatrix[i_count + 1, 2] = -Mgr.cos(15.0 * Tscope.t) * Mgr.tn(Delta);
                    Fmatrix[i_count + 1, 3] = 1.0;
                    Fmatrix[i_count + 1, 4] = (NWE == 1) ? Mgr.tn(Delta) : -Mgr.tn(Delta);
                    Fmatrix[i_count + 1, 5] = (NWE == 1) ? 1.0 / Mgr.cos(Delta) : -1.0 / Mgr.cos(Delta);
                    Fmatrix[i_count + 1, 6] = Mgr.sin(Fi) * Mgr.cos(Fi) * Mgr.sin(Delta) * Mgr.tn(Delta) + Mgr.cos(Fi) * Mgr.cos(15.0 * Tscope.t) * Mgr.cos(Fi) * Mgr.sin(Delta);
                    if (NWE != 1)
                    {
                        Fmatrix[i_count + 1, 6] = -Fmatrix[i_count + 1, 6];
                    }


                    Fmatrix[i_count + 1, 7] = Mgr.cos(Fi) * Mgr.sin(15.0 * Tscope.t) / Mgr.cos(Delta);

                    for (int j = 0; j < dg1.RowCount; j++)
                    {
                        Fmatrix[i_count, j] *= Convert.ToInt32(dg1[0, j].Value);
                        Fmatrix[i_count+1, j] *= Convert.ToInt32(dg1[0, j].Value);
                    }
                    
                    OC[i_count] = din - Tscope.DeltaL2;
                    OC[i_count + 1] = tin - Tscope.tL2;
                    lOper.Items.Add(Convert.ToString(Sec) + " " + Convert.ToString(i_count) + "   " + Convert.ToString(OC[i_count + 1]) + "   " + Convert.ToString(OC[i_count]) + "   " + Convert.ToString(Tscope.t) + "   " + Convert.ToString(Tscope.tL2) + "   " + Convert.ToString(Tscope.DeltaL2));
                    

                    i_count += 2;
                } //if

                if (method == 1)
                {
                    dg1[1, 0].Value = "Нуль пункт по часовому углу";
                    dg1[1, 1].Value = "Установки пол.оси по азимуту";
                    dg1[1, 2].Value = "Установки пол.оси по широте";
                    dg1[1, 3].Value = "Наклонности";
                    dg1[1, 4].Value = "Коллимации";
                    dg1[1, 5].Value = "Гнутия по прямому восхождению";
                    

                    Fmatrix[i_count, 0] = 1.0;
                    Fmatrix[i_count, 1] = Mgr.cos(Fi) * Mgr.sin(15.0 * Tscope.t) * Mgr.tn(Delta);
                    Fmatrix[i_count, 2] = -Mgr.cos(15.0 * Tscope.t) * Mgr.tn(Delta);
                    Fmatrix[i_count, 3] = (NWE == 1) ? Mgr.tn(Delta) : -Mgr.tn(Delta);
                    Fmatrix[i_count, 4] = (NWE == 1) ? 1.0 / Mgr.cos(Delta) : -1 / Mgr.cos(Delta);
                    Fmatrix[i_count, 5] = Mgr.cos(Fi) * Mgr.sin(15.0 * Tscope.t) / Mgr.cos(Delta);

                    for (int j = 0; j < dg1.RowCount; j++)
                    {
                        Fmatrix[i_count, j] *= Convert.ToInt32(dg1[0, j].Value);
                    }


                    OC[i_count] = tin - Tscope.tL2 + PopHor2(PopHor(Tscope.t), PopVert(Tscope.t)) / 60 + PopVert2(PopHor(Tscope.t), PopVert(Tscope.t)) / 60;
                    lOper.Items.Add(Convert.ToString(Sec) + " " + Convert.ToString(i_count) + "   " +  "   " + Convert.ToString(OC[i_count]) + "   " + Convert.ToString(Tscope.t) + "   " + Convert.ToString(Tscope.tL2) + "   " + Convert.ToString(Tscope.DeltaL2));
                    i_count += 1;
                }

                if (method == 2)
                {
                    dg1[1, 0].Value = "Нуль-пункт по склонению";
                    dg1[1, 1].Value = "Установки пол.оси по азимуту";
                    dg1[1, 2].Value = "Установки пол.оси по склонению";
                    dg1[1, 3].Value = "Гнутия по склонению";
                    
                    Fmatrix[i_count, 0] = -1.0;
                    Fmatrix[i_count, 1] = Znnn * Mgr.cos(Fi) * Mgr.cos(15.0 * Tscope.t);
                    Fmatrix[i_count, 2] = Znnn * Mgr.sin(15.0 * Tscope.t);
                    Fmatrix[i_count, 3] = Znnn * (-Mgr.sin(Fi) * Mgr.cos(Delta) - Mgr.cos(Fi) * Mgr.sin(Delta) * Mgr.cos(15.0 * Tscope.t));

                    for (int j = 0; j < dg1.RowCount; j++)
                    {
                        Fmatrix[i_count, j] *= Convert.ToInt32(dg1[0, j].Value);
                    }

                    OC[i_count] = din - Tscope.DeltaL2;
                    lOper.Items.Add(Convert.ToString(Sec) + " " + Convert.ToString(i_count) + "   "+ "   " + Convert.ToString(OC[i_count]) + "   " + Convert.ToString(Tscope.t) + "   " + Convert.ToString(Tscope.tL2) + "   " + Convert.ToString(Tscope.DeltaL2));
                    i_count += 1;
                }

                
            
            } //while

            //Start least squares method
            //MessageBox.Show(Convert.ToString(i_count));
            lsfit.lsfitlinear(ref OC, ref Fmatrix, i_count, MM, ref info, ref c, ref rep);
            //MessageBox.Show(Convert.ToString(info));
            for (int ii = 0; ii < MM; ii++)
            {
                dg1[2, ii].Value = Math.Round(c[ii]*60,5);
                lOper.Items.Add(Convert.ToString(c[ii]));
            }
            lOper.Items.Add(" ");
            lSummary.Items.Add("Method " + Convert.ToString(method)); //7th symbol 
            lSummary.Items.Add("FileName " + fobs);
            lSummary.Items.Add("Average error     " + Convert.ToString(rep.avgerror*60));
            lSummary.Items.Add("Maximum error     "+Convert.ToString(rep.maxerror*60));


            lOper.Items.Add(rep.avgerror);
            lOper.Items.Add(rep.maxerror);
            // zdes' budem probivat' schitat' nevyazku
            for (int ii = 0; ii < NN; ii++)
            {
                double ypol=0;
                for (int jj = 0; jj < MM; jj++)
                {
                    ypol += Fmatrix[ii, jj] * c[jj];
                }
                lOper.Items.Add(OC[ii]-ypol);
                //dg2[0, ii].Value = (int)ii / 2;
                if (method == 0) dg2[0, ii].Value = (int)ii / 2;
                else dg2[0, ii].Value = ii;

                dg2[1, ii].Value = Math.Round(OC[ii]*60,5);
                dg2[2, ii].Value = Math.Round(ypol*60,5);
                dg2[3, ii].Value = Math.Round((OC[ii] - ypol)*60,5);



                /*for (int jj = 0; jj < 4; jj++)
                {
                    if (Math.Abs(OC[ii] - ypol) / 5 > rep.avgerror)
                    {
                        dg2[jj, ii].Style.BackColor = Color.Red;
                    }
                }*/
            }


            Fform2.lResult.Text += "\n";
            
            Fform2.lResult.Text += "---------------------------";
            Fform2.lResult.Text += "\n";
                     

            for (int i = 0; i < dg1.RowCount; i++)
            {
                string str = "";
                for (int j = 0; j < dg1.ColumnCount; j++)
                {
                    str += dg1[j, i].Value + "   ";
                }
                Fform2.lResult.Text += (str+"\n");
            }
            

            //Обработка и раскраска
            double xav = 0;
            for (int ii = 0; ii < dg2.Rows.Count; ii++)
            {
                xav += Convert.ToDouble(dg2[3, ii].Value);
            }
            xav = xav / dg2.Rows.Count;
           
            double ss = 0;
            for (int ii = 0; ii < dg2.Rows.Count; ii++)
            {
                ss += (Convert.ToDouble(dg2[3, ii].Value) - xav) * (Convert.ToDouble(dg2[3, ii].Value) - xav);
            }
            double ss2 = Math.Sqrt(ss / dg2.Rows.Count);

            for (int ii = 0; ii < dg2.Rows.Count; ii++)
            {
                if (Math.Abs(Convert.ToDouble(dg2[3, ii].Value)) > 3*ss2)
                {
                    for (int jj = 0; jj < 6; jj++)
                    {
                        dg2[jj, ii].Style.BackColor = Color.Red;
                    }
                }
            }




            fin.Close();
        }

        private void comboBox1_SelectedIndexChanged(object sender, EventArgs e)
        {

        }

        private void listBox1_SelectedIndexChanged(object sender, EventArgs e)
        {

        }

        private void dg_CellContentClick(object sender, DataGridViewCellEventArgs e)
        {

        }

        private void lOper_SelectedIndexChanged(object sender, EventArgs e)
        {

        }

        private void lOper_DoubleClick(object sender, EventArgs e)
        {
            

        }

        private void bfileobs_Click(object sender, EventArgs e)
        {
            OpenFileDialog ofd = new OpenFileDialog();
            ofd.Title = "Выберите файл с наблюдениями";
            ofd.InitialDirectory = "";
            ofd.Filter = "txt files (*.txt,*.dat,*.log)|*.txt;*.dat;*.log;*.temp;|All files (*.*)|*.*";
            if (ofd.ShowDialog() == DialogResult.OK)
            {
                try
                {
                    tfileobs.Text = ofd.FileName;
                }
                catch (Exception ex)
                {
                    MessageBox.Show("Error: Could not read file from disk. Original error: " + ex.Message);
                }
            }
        }

        private void bfiletscope_Click(object sender, EventArgs e)
        {
            OpenFileDialog ofd = new OpenFileDialog();
            ofd.Title = "Выберите файл с настройками телескопа";
            ofd.Filter = "txt files (*.txt,*.dat)|*.txt;*.dat|All files (*.*)|*.*";
            ofd.InitialDirectory = "";
            if (ofd.ShowDialog() == DialogResult.OK)
            {
                try
                {
                    tfiletscope.Text = ofd.FileName;
                }
                catch (Exception ex)
                {
                    MessageBox.Show("Error: Could not read file from disk. Original error: " + ex.Message);
                }
            }
        }

        private void Form1_Load(object sender, EventArgs e)
        {
            cbMethod.SelectedIndex = 0;
            
        }

        private void Form1_Shown(object sender, EventArgs e)
        {
            cbMethod.SelectedIndex = 0;
            //MessageBox.Show(Convert.ToString(Top));
            
        }

        private void cbMethod_SelectedIndexChanged(object sender, EventArgs e)
        {

        }

        private void bDelsel_Click(object sender, EventArgs e)
        {
            int method = Convert.ToInt32((lSummary.Items[0].ToString())[7])-48 ;
            if (dg2.SelectedRows.Count >= 1)
            {
                   
                        if ((method == 1) || (method == 2))
                        {
                            foreach (DataGridViewRow dgvRow in dg2.SelectedRows)
                            {
                                lOper2.Items.RemoveAt(dgvRow.Index);
                                
                                dg2.Rows.RemoveAt(dgvRow.Index);
                                
                            }
                        }

                        if (method == 0)
                        {
                            foreach (DataGridViewRow dgvRow in dg2.SelectedRows)
                            {
                                if (dgvRow.Index % 2 == 1)
                                {
                                    dg2.Rows[dgvRow.Index - 1].Selected = true;
                                    dg2.Rows[dgvRow.Index].Selected = false;
                                    
                                }
                                
                            }
                            foreach (DataGridViewRow dgvRow in dg2.SelectedRows)
                            {
                                if (dgvRow.Index % 2 == 0)
                                {
                                    lOper2.Items.RemoveAt(dgvRow.Index / 2);
                                    int ind = dgvRow.Index;
                                    dg2.Rows.RemoveAt(ind);
                                    dg2.Rows.RemoveAt(ind);
                                }
                            }
                           
                        }

                       
                                   
            }
        }

        private void bRecalc_Click(object sender, EventArgs e)
        {
            string st=lSummary.Items[1].ToString();
            string fin = st.Substring(9, st.Length - 9);
            //MessageBox.Show(fin);
            StreamWriter sr = new StreamWriter(fin);
            for (int i = 0; i < lOper2.Items.Count; i++)
            {
                sr.WriteLine(lOper2.Items[i].ToString());
            }
            sr.Close();


            for (int ii = 0; ii < dg2.Rows.Count; ii++)
            {
                    for (int jj = 0; jj < 4; jj++)
                    {
                        dg2[jj, ii].Style.BackColor = Color.White;
                    }
            }
            
            
            Calculate(true,fin);
        }

        private void pictureBox1_Click(object sender, EventArgs e)
        {

        }

        private void bSave_Click(object sender, EventArgs e)
        {
            SaveFileDialog sfd = new SaveFileDialog();
            sfd.Filter = "txt files (*.txt,*.dat)|*.txt;*.dat|All files (*.*)|*.*";
            if (sfd.ShowDialog() == DialogResult.OK)
            {
                try
                {
                    StreamWriter resfile = new StreamWriter(sfd.FileName);

                    for (int i = 0; i < dg1.RowCount;i++ )
                    {
                        string str="";
                        for (int j=0;j<dg1.ColumnCount;j++)
                        {
                            str+=dg1[j,i].Value+"   ";
                        }
                        resfile.WriteLine(str);
                    }
                    resfile.Close();
                }
                catch (Exception ex)
                {
                    MessageBox.Show("Error: Could not read file from disk. Original error: " + ex.Message);
                }
            }
        }

        private void bAppend_Click(object sender, EventArgs e)
        {
            SaveFileDialog sfd = new SaveFileDialog();
            sfd.Filter = "txt files (*.txt,*.dat)|*.txt;*.dat|All files (*.*)|*.*";
            sfd.Title = "Дописать в файл";
            if (sfd.ShowDialog() == DialogResult.OK)
            {
                try
                {
                    StreamWriter resfile = new StreamWriter(sfd.FileName,true);
                    resfile.WriteLine();
                    resfile.WriteLine("---------------------------");
                    resfile.WriteLine();

                    for (int i = 0; i < dg1.RowCount; i++)
                    {
                        string str = "";
                        for (int j = 0; j < dg1.ColumnCount; j++)
                        {
                            str += dg1[j, i].Value + "   ";
                        }
                        resfile.WriteLine(str);
                    }
                    resfile.Close();
                }
                catch (Exception ex)
                {
                    MessageBox.Show("Error: Could not read file from disk. Original error: " + ex.Message);
                }
            }
        }

        private void button1_Click_1(object sender, EventArgs e)
        {


            //Fform2.Top = Top;

            if (Fform2.Visible == false)
            {
                Fform2.Show();
                bExt.Text = ">";
            }

            else
            {
                Fform2.Hide();
                bExt.Text = "<";
            }
            
            
        }
    }

    public class Mathgr
    {
        public double rad = 0.0174532925199433;

        public double sin(double x)
        {
            return Math.Sin(x * rad);
        }
        public double cos(double x)
        {
            return Math.Cos(x * rad);
        }
        public double asn(double x)
        {
            return Math.Asin(x) / rad;

        }
        public double acs(double x)
        {
            return Math.Acos(x) / rad;
        }
        public double tn(double x)
        {
            return Math.Tan(x * rad);

        }
        public double atn(double y,double x)
        {
            return Math.Atan2(y, x) / rad;
        }

        public double atn2(double y, double x)
        {
            double rad = 0.0174532925199433;
            double ax, ay, phi;
            if ((x == 0) && (y == 0)) return 0.0;
            else
            {
                ax = Math.Abs(x);
                ay = Math.Abs(y);
                if (ax > ay) phi = Math.Atan(ay / ax) / rad;
                else phi = 90.0 - Math.Atan(ax / ay) / rad;
                if (x < 0) phi = 180.0 - phi;
                if (y < 0) phi = -phi;
                return phi;
            
            }
        }


    }
    public class Calendar
    {
        public int year;
        public int month;
        public double day;
        public double juldate;
        public double st; //startime
        public Calendar(int y,int m,double d)
        {
            month=m;
            year=y;
            day=d;
            JULDAT();
        }
        public void JULDAT() //Юлианская дата для любой даты
        {
            double a = 10000.0 * year + 100.0 * month + day;
            double b;
            if (month <= 2)
            {
                month += 12;
                year--;
            }
            if (a <= 15821004.1)
            {
                b = -2 + Math.Truncate((double)(year + 4716.0) / 4.0) - 1179;
            }
            else
            {
                b = Math.Truncate((double)year / 400) - Math.Truncate((double)year / 100) + Math.Truncate((double)year / 4);
            }
            a = 365.0 * year - 679004.0;
            juldate = 2400000.5 + a + b + Math.Truncate(30.6001 * (month + 1)) + day;
        }
        
        // Определение календарной даты по юлианской
        public void calcdat()
        {
            double b,d,f;
            double jd0,c,e;
            jd0=Math.Floor(juldate+0.5);
            if (jd0<2299161.0)
            {
                b=0;
                c=jd0+1524.0; //julianskiy calendar
            }
            else
            {
                b=Math.Truncate((jd0-1867216.25)/36524.25);
                c=jd0+(b-Math.Truncate(b/4.0))+1525.0; //grigiriansky calendar
            }
            d=Math.Truncate((c-122.1)/365.25);
            e=365.0*d+Math.Truncate(d/4.0);
            f=Math.Truncate((c-e)/30.6001);                 //!!!! Pomenyal int na double
            day=Math.Truncate(c-e+0.5)-Math.Truncate(30.6001*f)+(juldate+0.5-jd0);
            month=(int)(f-1-12*Math.Truncate(f/14.0));
            year=(int)(d-4715-Math.Truncate((7+month)/10.0));

        }
        public double startime(int y,double jdt,double gmt)
        {
            Calendar j1= new Calendar(y,1,0.0);
            double t=(j1.juldate-2415020.0)/36525.0;
            double r=6.6460656+2400.051262*t+0.00002581*t*t;
            double u=r-24*(year-1900);
            double b=24-u;
            double c=Math.Truncate(jdt-j1.juldate)*0.0657098-b; //v originale zdesi ne trunc a int()
            double d=gmt*1.002738;
            d+=c;
            if (d<0) d+=24;
            if (d>24) d-=24;
            st=d;
            return st;
        }




    }

    public class Telescope
    {
        public double alf;
        public double del;

        public double RKEP1(double SA, double EX)
        {
            double eps = 0.000001;
            double SAR = SA * Math.PI / 180;
            double EA=SAR;
            double D, DEA;
            D = EA - EX * Math.Sin(EA) - SAR;
            while (Math.Abs(D) > eps)
            {
                DEA = D / (1 - EX * Math.Cos(EA));
                EA -= DEA;
                D = EA - EX * Math.Sin(EA) - SAR;    //Zdes' izmenyal
            }
            EA = EA / (Math.PI / 180);
            return EA;
        }
        public void ALFDEL1(double lam,double bet,double jd) //NEPONYATNO GDE ISPOLZUETSA
        {
            double tb=(jd-2415020)/36525;
            double ei=23.452294-0.0130125*tb-0.0000016*tb*tb+0.00000028*tb*tb*tb;
            double s1=Math.Sin(lam*Math.PI/180)*Math.Cos(ei*Math.PI/180)-Math.Sin(bet*Math.PI/180)/Math.Cos(bet*Math.PI/180)*Math.Sin(ei*Math.PI/180);
            double s2=Math.Cos(lam*Math.PI/180);
            alf=Math.Atan2(s1,s2);
            if (alf < 0) alf += 360;
            alf = alf / 15;
            double sdel = Math.Sin(bet * Math.PI / 180) * Math.Cos(ei * Math.PI / 180) + Math.Cos(bet * Math.PI / 180) * Math.Sin(ei * Math.PI / 180) * Math.Sin(lam * Math.PI / 180);
            del=Math.Atan2(sdel, Math.Sqrt(1-sdel*sdel));
        }
        public double ESPn(double JD)
        {
            Mathgr Mgr = new Mathgr();
            double C = 360.0;
            double T=(JD-2415020)/36525;
            double L0=279.6966800+36000.7689200*T+0.0003025*T*T;
            while(L0<0) 
            {
               L0=L0+C;
            }
            while(L0>=C)
            {
                L0=L0-C;
            }
            double M0=358.4758300+35999.0497500*T-0.0001500*T*T-0.0000033*T*T*T;
            while (M0<0) 
            {
                M0=M0+C;
            } 
            while(M0>=C) 
            {
                M0=M0-C;
            }
            double E0=0.01675104-0.0000418*T-0.000000126*T*T;
            //double A0=1.0000002;
            double U0=259.1800000-1934.1420000*T;
            double U1=153.23+22518.7541*T;
            double U2=216.57+45037.5082*T;
            double U3=312.69+32964.3577*T;
            double U4=350.74+445267.1142*T-0.00144*T*T;
            double U5=231.19+20.20*T;
            double U6=353.40E0+65928.7155E0*T;
            double EH=23.4522940-0.0130125*T-0.00000164*T*T+0.000000503*T*T*T+0.0025600*Math.Cos(U0*Math.PI/180);
            double DI0=0.0013400*Math.Cos(U1*Math.PI/180)+0.00154*Math.Cos(U2*Math.PI/180)+0.00200*Math.Cos(U3*Math.PI/180)+0.0017900*Math.Sin(U4*Math.PI/180)+
                0.0017800*Math.Sin(U5*Math.PI/180)-0.0056900-0.0047900*Math.Sin(U0*Math.PI/180);
            double DR0=0.00000543*Math.Sin(U1*Math.PI/180)+0.00001575*Math.Sin(U2*Math.PI/180)+
                0.00001627*Math.Sin(U3*Math.PI/180)+0.00003076*Math.Cos(U4*Math.PI/180)+
                0.00000927*Math.Sin(U6*Math.PI/180);

            double EA0=RKEP1(M0,E0);
            double IA0=2*Mgr.atn2((Math.Sqrt((1+E0)/(1-E0))*Math.Sin(EA0*Math.PI/360)),(Math.Cos(EA0*Math.PI/360)));
            double IL0=L0+IA0-M0+DI0;
            return IL0;
        }

        /*Z - з_-ит-R_ р сстRя-и_ в _р дус х,
           T - т_мп_р тур  вRздух  в _р дус х ц_<ьсия,
        P - д в<_-и_ вRздух  в мм ртут-R_R стR<б ,
         R - р_фр кция в с_ку-д х ду_и }*/
        public double Refr(double Z,double T,double P)
        {
            double F=0.0005;
            double D=0;
            double FI=0;
            double H=0;
            Mathgr MathG=new Mathgr();
            double T1=MathG.tn(Z);
            double T2=T1*T1;
            double T3=T2*T1;
            double T4=T3*T1;
            double E=1.0;
            double EM=0.003681;
            double LR0=Math.Log10(T1)+1.756522-5.11E-4*T2+1.09E-6*T4;
            double GAM=Math.Log10((E+EM*15.0)/(E+EM*T));
            double LAM=E-1.2E-3*T1+1.92E-3*T2-(7.5E-5+3.3E-7*T)*T3;
            double B=Math.Log10(P/760.0);
            F=F*1.333224;
            double A=1.0004+1.1E-4*T2;
            double C=(-6.6E-5*F-2.41E-6*F*F)*(E+0.001E0*T2);
            double LR=LR0+LAM*GAM+A*B+C+D+FI+H;
            double R=Math.Pow(10,LR);
            return R;

        }



        public double Fi,Lamb,Dt,Temp,Dav,ex_t,Eex_t,ex_d,Eex_d;
        public int JD_R,Alf1,Alf2;
        public double Alf3;
        public string WE;
        public int Zn_Delt;
        public int Delt1,Delt2;
        public double Delt3;
        public int God, Mes, Den, Hour, Min;
        public double Sec;
        public double t, tL2, DeltaL2;

                
        public void Ustanovka()
        {
            if (Delt1 > 87)
            {
                MessageBox.Show("Delta>87");
                return;
            }

            Calendar JDC = new Calendar(God, Mes, Den + (Hour + Min / 60.0 + Sec / 3600.0) / 24.0);
//            Calendar JD0C = new Calendar(God0, Mes0, Den0 + Hour0 / 24.0);
            
//            double AA=Den+(Hour+Min/60.0+(Sec+ds0+(JDC.juldate-JD0C.juldate)*V)/3600.0)/24.0;
//            JDC=new Calendar(God,Mes,AA);     //Какая-то хренотень
            JDC.calcdat();                              
            JDC.JULDAT();
            double D = JDC.day;
            int G = JDC.year;
            double A0=(D-Math.Truncate(D))*24;
            double H=Math.Truncate(A0);
            double A1=A0-H;
            double Mi=Math.Truncate(A1*60);
            double A2=A1-Mi/60.0;
            double S=A2*60*60;

            double Gmt=H+Mi/60.0+S/3600.0-Dt;
            double STg=JDC.startime(G,JDC.juldate,Gmt);
            double STm=STg+Lamb;
            if (STm>24) STm-=24.0;
            if (STm<0) STm+=24.0;

            double STm1=Math.Truncate(STm);
            A1=STm-STm1;
            double STm2=Math.Truncate(A1*60);
            A2=A1-STm2/60.0;
            double STm3=A2*60.0*60.0;

            double t0=(JDC.juldate-2415020.0)/36525.0;
            double mm=3.07234+0.00186*t0;
            double n=20.0467-0.0085*t0;
            double Alf=Alf1+Alf2/60.0+Alf3/3600.0;
            double Zn;
            if ((Zn_Delt==1)||(1==2/*Zn_Delt==" "*/)) Zn=1;
                    else Zn=-1;
            double Delt=Zn*(Delt1+Delt2/60.0+Delt3/3600.0);
            Calendar JDTC=new Calendar(JD_R,1,0.0);
            Mathgr Mgr=new Mathgr();
            double Alfa=Alf+(JDC.juldate-JDTC.juldate)/365.25*(mm+n/15.0*Mgr.sin(Alf*15.0)*Mgr.tn(Delt))/3600.0;
            double test_1 = (mm + n / 15.0 + Mgr.sin(Alf * 15.0) * Mgr.tn(Delt)) ;
            double Delta=Delt+(JDC.juldate-JDTC.juldate)/365.25*(n*Mgr.cos(Alf*15))/3600.0;

            t = STm - Alfa;
            if (Math.Abs(t) > 12.0) t = -24.0 * t / Math.Abs(t) + t;

            double z = Mgr.acs(Mgr.sin(Delta)*Mgr.sin(Fi)+Mgr.cos(Delta)*Mgr.cos(Fi)*Mgr.cos(t*15.0));
            if (Math.Abs(z) > 87) { MessageBox.Show("Z>87"); return; };
            double Ref = Refr(z, Temp, Dav);
            double Refa, Refd;
            if (z != 0)
            {
                A0 = Mgr.acs((Mgr.sin(Fi) - Mgr.sin(Delta) * Mgr.cos(z)) / (Mgr.cos(Delta) * Mgr.sin(z)));
                if (t != 0)
                    Refa = Ref * Mgr.sin(A0) / Mgr.cos(Delta) * t / Math.Abs(t) / 15.0;
                else
                    Refa = 0;
                Refd = Ref * Mgr.cos(A0);
            }
            else
            {
                Refa = 0;
                Refd = 0;
            }


            double tL0 = STm - (Alfa + Refa / (3600.0 * 15.0));
            if (Math.Abs(tL0) > 12.0) tL0 = -24.0 * tL0 / Math.Abs(tL0) + tL0;
            double tLOL = tL0 * 15.0;
            if (tLOL < 0)
                tLOL = -tLOL;
            else
                tLOL = 360 - tLOL;
            if (WE == "E")
                tLOL = tLOL + 180.0;
            if (tLOL > 360)
                tLOL -= 360;
            if (tLOL < 0)
                tLOL += 360.0;


            double IL0 = ESPn(JDC.juldate);
            double beta = -20.4958 / 3600.0;
            double d_alf = beta * Mgr.cos(IL0) * Mgr.cos(23.5) * Mgr.cos(Alfa) / Mgr.cos(Delta) + beta * Mgr.sin(IL0) * Mgr.sin(Alfa) / Mgr.cos(Delta);
            double d_delt = beta * Mgr.cos(IL0) * Mgr.cos(23.5) * (Mgr.tn(23.5) * Mgr.cos(Delta) - Mgr.sin(Alfa) * Mgr.sin(Delta)) + beta * Mgr.sin(IL0) * Mgr.cos(Alfa) * Mgr.sin(Delta);

            double AlfaL2 = Alfa + d_alf / 15.0 + Refa / (3600 * 15);
            tL2 = STm - AlfaL2;
            if (Math.Abs(tL2) > 12.0)
                tL2 = -24.0 * tL2 / Math.Abs(tL2) + tL2;
            tL2 = tL2 * 15.0;
            if (tL2 < 0)
                tL2 = -tL2;
            else
                tL2 = 360 - tL2;
            if (WE == "E")
                tL2 += 180;
            double D_ex_t = (ex_t * Mgr.sin(tL2 - Eex_t)) / 60.0;
            tL2 = tL2 + D_ex_t;
            if (tL2 > 360)
                tL2 -= 360;
            if (tL2 < 0)
                tL2 += 360;

            DeltaL2 = Delta + d_delt + Refd / 3600.0;
            if (WE == "W")
                DeltaL2 = DeltaL2+0.0;
            if (WE == "E")
                DeltaL2 = 180.0 - DeltaL2;


            double D_ex_d = (ex_d * Mgr.sin(DeltaL2 - Eex_d)) / 60.0;
            DeltaL2 = DeltaL2 + D_ex_d;

            double Mz=1.0/Mgr.cos(z)-40.0*Math.Exp(-2.25*Math.Log(90.0-z));



        }




    }

   
    
    

}
