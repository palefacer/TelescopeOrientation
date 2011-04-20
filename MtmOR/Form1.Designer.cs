namespace MtmOR
{
    partial class Form1
    {
        /// <summary>
        /// Требуется переменная конструктора.
        /// </summary>
        private System.ComponentModel.IContainer components = null;

        /// <summary>
        /// Освободить все используемые ресурсы.
        /// </summary>
        /// <param name="disposing">истинно, если управляемый ресурс должен быть удален; иначе ложно.</param>
        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null))
            {
                components.Dispose();
            }
            base.Dispose(disposing);
        }

        #region Код, автоматически созданный конструктором форм Windows

        /// <summary>
        /// Обязательный метод для поддержки конструктора - не изменяйте
        /// содержимое данного метода при помощи редактора кода.
        /// </summary>
        private void InitializeComponent()
        {
            this.bStart = new System.Windows.Forms.Button();
            this.lOper = new System.Windows.Forms.ListBox();
            this.dg1 = new System.Windows.Forms.DataGridView();
            this.boxcol = new System.Windows.Forms.DataGridViewCheckBoxColumn();
            this.dg2 = new System.Windows.Forms.DataGridView();
            this.lSummary = new System.Windows.Forms.ListBox();
            this.cbMethod = new System.Windows.Forms.ComboBox();
            this.tfileobs = new System.Windows.Forms.TextBox();
            this.tfiletscope = new System.Windows.Forms.TextBox();
            this.bfileobs = new System.Windows.Forms.Button();
            this.bfiletscope = new System.Windows.Forms.Button();
            this.textBox1 = new System.Windows.Forms.TextBox();
            this.lObserv = new System.Windows.Forms.ListBox();
            this.progressBar1 = new System.Windows.Forms.ProgressBar();
            this.lOper2 = new System.Windows.Forms.ListBox();
            this.bDelsel = new System.Windows.Forms.Button();
            this.bRecalc = new System.Windows.Forms.Button();
            this.bSave = new System.Windows.Forms.Button();
            this.bAppend = new System.Windows.Forms.Button();
            this.bExt = new System.Windows.Forms.Button();
            this.pictureBox1 = new System.Windows.Forms.PictureBox();
            this.tangle = new System.Windows.Forms.TextBox();
            ((System.ComponentModel.ISupportInitialize)(this.dg1)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.dg2)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.pictureBox1)).BeginInit();
            this.SuspendLayout();
            // 
            // bStart
            // 
            this.bStart.ForeColor = System.Drawing.Color.Red;
            this.bStart.Location = new System.Drawing.Point(29, 12);
            this.bStart.Name = "bStart";
            this.bStart.Size = new System.Drawing.Size(117, 62);
            this.bStart.TabIndex = 0;
            this.bStart.Text = "Загрузить данные и выполнить расчет";
            this.bStart.UseVisualStyleBackColor = true;
            this.bStart.Click += new System.EventHandler(this.button1_Click);
            // 
            // lOper
            // 
            this.lOper.FormattingEnabled = true;
            this.lOper.Location = new System.Drawing.Point(1030, 492);
            this.lOper.Name = "lOper";
            this.lOper.Size = new System.Drawing.Size(99, 82);
            this.lOper.TabIndex = 1;
            this.lOper.Visible = false;
            this.lOper.SelectedIndexChanged += new System.EventHandler(this.lOper_SelectedIndexChanged);
            this.lOper.DoubleClick += new System.EventHandler(this.lOper_DoubleClick);
            // 
            // dg1
            // 
            this.dg1.AllowUserToAddRows = false;
            this.dg1.AllowUserToDeleteRows = false;
            this.dg1.AllowUserToResizeColumns = false;
            this.dg1.AllowUserToResizeRows = false;
            this.dg1.ColumnHeadersHeightSizeMode = System.Windows.Forms.DataGridViewColumnHeadersHeightSizeMode.AutoSize;
            this.dg1.Columns.AddRange(new System.Windows.Forms.DataGridViewColumn[] {
            this.boxcol});
            this.dg1.ImeMode = System.Windows.Forms.ImeMode.NoControl;
            this.dg1.Location = new System.Drawing.Point(29, 124);
            this.dg1.Name = "dg1";
            this.dg1.RowHeadersWidth = 5;
            this.dg1.SelectionMode = System.Windows.Forms.DataGridViewSelectionMode.FullRowSelect;
            this.dg1.Size = new System.Drawing.Size(362, 321);
            this.dg1.TabIndex = 5;
            this.dg1.CellContentClick += new System.Windows.Forms.DataGridViewCellEventHandler(this.dg_CellContentClick);
            // 
            // boxcol
            // 
            this.boxcol.HeaderText = "";
            this.boxcol.Name = "boxcol";
            this.boxcol.Resizable = System.Windows.Forms.DataGridViewTriState.False;
            this.boxcol.Width = 20;
            // 
            // dg2
            // 
            this.dg2.AllowUserToAddRows = false;
            this.dg2.AllowUserToDeleteRows = false;
            this.dg2.AllowUserToResizeColumns = false;
            this.dg2.AllowUserToResizeRows = false;
            this.dg2.ColumnHeadersHeightSizeMode = System.Windows.Forms.DataGridViewColumnHeadersHeightSizeMode.AutoSize;
            this.dg2.Location = new System.Drawing.Point(397, 124);
            this.dg2.Name = "dg2";
            this.dg2.ReadOnly = true;
            this.dg2.RowHeadersWidth = 25;
            this.dg2.SelectionMode = System.Windows.Forms.DataGridViewSelectionMode.FullRowSelect;
            this.dg2.Size = new System.Drawing.Size(453, 337);
            this.dg2.TabIndex = 6;
            // 
            // lSummary
            // 
            this.lSummary.Font = new System.Drawing.Font("Courier New", 8.25F, System.Drawing.FontStyle.Regular, System.Drawing.GraphicsUnit.Point, ((byte)(204)));
            this.lSummary.FormattingEnabled = true;
            this.lSummary.ItemHeight = 14;
            this.lSummary.Location = new System.Drawing.Point(29, 451);
            this.lSummary.Name = "lSummary";
            this.lSummary.Size = new System.Drawing.Size(248, 74);
            this.lSummary.TabIndex = 7;
            // 
            // cbMethod
            // 
            this.cbMethod.DisplayMember = "1";
            this.cbMethod.DropDownStyle = System.Windows.Forms.ComboBoxStyle.DropDownList;
            this.cbMethod.Items.AddRange(new object[] {
            "Совместное",
            "Пр. восхождение",
            "Склонение"});
            this.cbMethod.Location = new System.Drawing.Point(29, 80);
            this.cbMethod.Name = "cbMethod";
            this.cbMethod.Size = new System.Drawing.Size(117, 21);
            this.cbMethod.TabIndex = 8;
            this.cbMethod.SelectedIndexChanged += new System.EventHandler(this.cbMethod_SelectedIndexChanged);
            // 
            // tfileobs
            // 
            this.tfileobs.Enabled = false;
            this.tfileobs.Location = new System.Drawing.Point(167, 12);
            this.tfileobs.Name = "tfileobs";
            this.tfileobs.Size = new System.Drawing.Size(464, 20);
            this.tfileobs.TabIndex = 9;
            // 
            // tfiletscope
            // 
            this.tfiletscope.Enabled = false;
            this.tfiletscope.Location = new System.Drawing.Point(167, 47);
            this.tfiletscope.Name = "tfiletscope";
            this.tfiletscope.Size = new System.Drawing.Size(465, 20);
            this.tfiletscope.TabIndex = 10;
            // 
            // bfileobs
            // 
            this.bfileobs.Location = new System.Drawing.Point(637, 5);
            this.bfileobs.Name = "bfileobs";
            this.bfileobs.Size = new System.Drawing.Size(79, 27);
            this.bfileobs.TabIndex = 11;
            this.bfileobs.Text = "Наблюдения";
            this.bfileobs.UseVisualStyleBackColor = true;
            this.bfileobs.Click += new System.EventHandler(this.bfileobs_Click);
            // 
            // bfiletscope
            // 
            this.bfiletscope.Location = new System.Drawing.Point(638, 43);
            this.bfiletscope.Name = "bfiletscope";
            this.bfiletscope.Size = new System.Drawing.Size(78, 27);
            this.bfiletscope.TabIndex = 12;
            this.bfiletscope.Text = "Телескоп";
            this.bfiletscope.UseVisualStyleBackColor = true;
            this.bfiletscope.Click += new System.EventHandler(this.bfiletscope_Click);
            // 
            // textBox1
            // 
            this.textBox1.Location = new System.Drawing.Point(56, 235);
            this.textBox1.Name = "textBox1";
            this.textBox1.Size = new System.Drawing.Size(117, 20);
            this.textBox1.TabIndex = 2;
            this.textBox1.Visible = false;
            // 
            // lObserv
            // 
            this.lObserv.Font = new System.Drawing.Font("Courier New", 8.25F, System.Drawing.FontStyle.Regular, System.Drawing.GraphicsUnit.Point, ((byte)(204)));
            this.lObserv.FormattingEnabled = true;
            this.lObserv.ItemHeight = 14;
            this.lObserv.Location = new System.Drawing.Point(29, 533);
            this.lObserv.Name = "lObserv";
            this.lObserv.Size = new System.Drawing.Size(821, 60);
            this.lObserv.TabIndex = 13;
            // 
            // progressBar1
            // 
            this.progressBar1.Location = new System.Drawing.Point(167, 80);
            this.progressBar1.Name = "progressBar1";
            this.progressBar1.Size = new System.Drawing.Size(464, 20);
            this.progressBar1.TabIndex = 14;
            this.progressBar1.Visible = false;
            // 
            // lOper2
            // 
            this.lOper2.Font = new System.Drawing.Font("Courier New", 8.25F, System.Drawing.FontStyle.Regular, System.Drawing.GraphicsUnit.Point, ((byte)(204)));
            this.lOper2.FormattingEnabled = true;
            this.lOper2.HorizontalScrollbar = true;
            this.lOper2.ItemHeight = 14;
            this.lOper2.Location = new System.Drawing.Point(852, 10);
            this.lOper2.Name = "lOper2";
            this.lOper2.Size = new System.Drawing.Size(422, 564);
            this.lOper2.TabIndex = 15;
            this.lOper2.Visible = false;
            // 
            // bDelsel
            // 
            this.bDelsel.Location = new System.Drawing.Point(397, 467);
            this.bDelsel.Name = "bDelsel";
            this.bDelsel.Size = new System.Drawing.Size(199, 45);
            this.bDelsel.TabIndex = 16;
            this.bDelsel.Text = "Удалить выбранные ";
            this.bDelsel.UseVisualStyleBackColor = true;
            this.bDelsel.Click += new System.EventHandler(this.bDelsel_Click);
            // 
            // bRecalc
            // 
            this.bRecalc.ForeColor = System.Drawing.Color.Red;
            this.bRecalc.Location = new System.Drawing.Point(651, 467);
            this.bRecalc.Name = "bRecalc";
            this.bRecalc.Size = new System.Drawing.Size(199, 45);
            this.bRecalc.TabIndex = 17;
            this.bRecalc.Text = "Пересчитать";
            this.bRecalc.UseVisualStyleBackColor = true;
            this.bRecalc.Click += new System.EventHandler(this.bRecalc_Click);
            // 
            // bSave
            // 
            this.bSave.Location = new System.Drawing.Point(288, 451);
            this.bSave.Name = "bSave";
            this.bSave.Size = new System.Drawing.Size(103, 35);
            this.bSave.TabIndex = 19;
            this.bSave.Text = "Сохранить в файл";
            this.bSave.UseVisualStyleBackColor = true;
            this.bSave.Click += new System.EventHandler(this.bSave_Click);
            // 
            // bAppend
            // 
            this.bAppend.Location = new System.Drawing.Point(288, 492);
            this.bAppend.Name = "bAppend";
            this.bAppend.Size = new System.Drawing.Size(103, 35);
            this.bAppend.TabIndex = 20;
            this.bAppend.Text = "Дописать в файл";
            this.bAppend.UseVisualStyleBackColor = true;
            this.bAppend.Click += new System.EventHandler(this.bAppend_Click);
            // 
            // bExt
            // 
            this.bExt.Font = new System.Drawing.Font("Microsoft Sans Serif", 18F, System.Drawing.FontStyle.Regular, System.Drawing.GraphicsUnit.Point, ((byte)(204)));
            this.bExt.Location = new System.Drawing.Point(-1, 124);
            this.bExt.Name = "bExt";
            this.bExt.Size = new System.Drawing.Size(24, 321);
            this.bExt.TabIndex = 21;
            this.bExt.Text = "<";
            this.bExt.UseVisualStyleBackColor = true;
            this.bExt.Click += new System.EventHandler(this.button1_Click_1);
            // 
            // pictureBox1
            // 
            this.pictureBox1.BorderStyle = System.Windows.Forms.BorderStyle.Fixed3D;
            this.pictureBox1.Image = global::MtmOR.Properties.Resources.mtm500;
            this.pictureBox1.InitialImage = global::MtmOR.Properties.Resources.mtm500;
            this.pictureBox1.Location = new System.Drawing.Point(748, 8);
            this.pictureBox1.Name = "pictureBox1";
            this.pictureBox1.Size = new System.Drawing.Size(102, 102);
            this.pictureBox1.TabIndex = 18;
            this.pictureBox1.TabStop = false;
            this.pictureBox1.Click += new System.EventHandler(this.pictureBox1_Click);
            // 
            // tangle
            // 
            this.tangle.Location = new System.Drawing.Point(644, 90);
            this.tangle.Name = "tangle";
            this.tangle.Size = new System.Drawing.Size(71, 20);
            this.tangle.TabIndex = 22;
            this.tangle.Text = "83";
            // 
            // Form1
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(867, 600);
            this.Controls.Add(this.tangle);
            this.Controls.Add(this.bExt);
            this.Controls.Add(this.bAppend);
            this.Controls.Add(this.bSave);
            this.Controls.Add(this.pictureBox1);
            this.Controls.Add(this.bRecalc);
            this.Controls.Add(this.bDelsel);
            this.Controls.Add(this.lOper2);
            this.Controls.Add(this.progressBar1);
            this.Controls.Add(this.lObserv);
            this.Controls.Add(this.bfiletscope);
            this.Controls.Add(this.bfileobs);
            this.Controls.Add(this.tfiletscope);
            this.Controls.Add(this.tfileobs);
            this.Controls.Add(this.lOper);
            this.Controls.Add(this.cbMethod);
            this.Controls.Add(this.lSummary);
            this.Controls.Add(this.dg2);
            this.Controls.Add(this.dg1);
            this.Controls.Add(this.textBox1);
            this.Controls.Add(this.bStart);
            this.Name = "Form1";
            this.Text = "Вычисление погрешностей наведения телескопа";
            this.Load += new System.EventHandler(this.Form1_Load);
            this.Shown += new System.EventHandler(this.Form1_Shown);
            ((System.ComponentModel.ISupportInitialize)(this.dg1)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.dg2)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.pictureBox1)).EndInit();
            this.ResumeLayout(false);
            this.PerformLayout();

        }

        #endregion

        private System.Windows.Forms.Button bStart;
        private System.Windows.Forms.ListBox lOper;
        private System.Windows.Forms.DataGridView dg1;
        private System.Windows.Forms.DataGridView dg2;
        private System.Windows.Forms.ListBox lSummary;
        private System.Windows.Forms.ComboBox cbMethod;
        private System.Windows.Forms.TextBox tfileobs;
        private System.Windows.Forms.TextBox tfiletscope;
        private System.Windows.Forms.Button bfileobs;
        private System.Windows.Forms.Button bfiletscope;
        private System.Windows.Forms.TextBox textBox1;
        private System.Windows.Forms.ListBox lObserv;
        private System.Windows.Forms.ProgressBar progressBar1;
        private System.Windows.Forms.ListBox lOper2;
        private System.Windows.Forms.Button bDelsel;
        private System.Windows.Forms.Button bRecalc;
        private System.Windows.Forms.Button bSave;
        private System.Windows.Forms.Button bAppend;
        private System.Windows.Forms.DataGridViewCheckBoxColumn boxcol;
        private System.Windows.Forms.Button bExt;
        private System.Windows.Forms.PictureBox pictureBox1;
        private System.Windows.Forms.TextBox tangle;
    }
}

