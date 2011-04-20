using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Text;
using System.Windows.Forms;

namespace MtmOR
{
    public partial class Form2 : Form
    {
        public Form2()
        {
            InitializeComponent();
            
            Form1 Form1 = this.Owner as Form1;
            if (Form1 != null)
            {
                //string s = Form1.textBox1.Text;
                //Form1.textBox1.Text = "OK";
            }


        }

        private void copyToolStripMenuItem_Click(object sender, EventArgs e)
        {
            lResult.Copy();
        }

        private void clearToolStripMenuItem_Click(object sender, EventArgs e)
        {
            lResult.Clear();
        }

        private void pasteToolStripMenuItem_Click(object sender, EventArgs e)
        {
            lResult.Paste();
        }

        private void Form2_FormClosing(object sender, FormClosingEventArgs e)
        {
            if (e.CloseReason == CloseReason.UserClosing)
            {
                Hide();
                e.Cancel = true;
            }
        }
    }
}
