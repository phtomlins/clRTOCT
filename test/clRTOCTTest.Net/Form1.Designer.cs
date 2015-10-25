namespace clRTOCTTest
{
    partial class Form1
    {
        /// <summary>
        /// Required designer variable.
        /// </summary>
        private System.ComponentModel.IContainer components = null;

        /// <summary>
        /// Clean up any resources being used.
        /// </summary>
        /// <param name="disposing">true if managed resources should be disposed; otherwise, false.</param>
        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null))
            {
                components.Dispose();
            }
            base.Dispose(disposing);
        }

        #region Windows Form Designer generated code

        /// <summary>
        /// Required method for Designer support - do not modify
        /// the contents of this method with the code editor.
        /// </summary>
        private void InitializeComponent()
        {
            this.btnLoadRawData = new System.Windows.Forms.Button();
            this.dlgOpenFile = new System.Windows.Forms.OpenFileDialog();
            this.lblFileSize = new System.Windows.Forms.Label();
            this.label1 = new System.Windows.Forms.Label();
            this.label2 = new System.Windows.Forms.Label();
            this.label3 = new System.Windows.Forms.Label();
            this.label4 = new System.Windows.Forms.Label();
            this.label5 = new System.Windows.Forms.Label();
            this.label6 = new System.Windows.Forms.Label();
            this.btnSetKernelSource = new System.Windows.Forms.Button();
            this.lblKernelSourcePath = new System.Windows.Forms.Label();
            this.nudNumAScansPerBScan = new System.Windows.Forms.NumericUpDown();
            this.nudAScanAveragingFactor = new System.Windows.Forms.NumericUpDown();
            this.nudBScanAveragingFactor = new System.Windows.Forms.NumericUpDown();
            this.nudNumBScans = new System.Windows.Forms.NumericUpDown();
            this.nudOutputAScanLength = new System.Windows.Forms.NumericUpDown();
            this.nudInputSpectraLength = new System.Windows.Forms.NumericUpDown();
            this.lblExpectedFileLength = new System.Windows.Forms.Label();
            this.btnProcess = new System.Windows.Forms.Button();
            this.nudDeviceIndex = new System.Windows.Forms.NumericUpDown();
            this.label7 = new System.Windows.Forms.Label();
            this.lblResamplingTable = new System.Windows.Forms.Label();
            this.btnLoadResamplingTable = new System.Windows.Forms.Button();
            this.lblReferenceSpectrum = new System.Windows.Forms.Label();
            this.btnLoadReferenceSpectrum = new System.Windows.Forms.Button();
            this.lblReferenceAScan = new System.Windows.Forms.Label();
            this.btnLoadReferenceAScan = new System.Windows.Forms.Button();
            this.txtResamplingTable = new System.Windows.Forms.TextBox();
            this.txtReferenceSpectrum = new System.Windows.Forms.TextBox();
            this.txtReferenceAScan = new System.Windows.Forms.TextBox();
            this.progLoadRaw = new System.Windows.Forms.ProgressBar();
            this.btnLoad = new System.Windows.Forms.Button();
            this.btnBScans = new System.Windows.Forms.Button();
            this.lblStatus = new System.Windows.Forms.Label();
            this.pictureBox1 = new System.Windows.Forms.PictureBox();
            this.label8 = new System.Windows.Forms.Label();
            this.label11 = new System.Windows.Forms.Label();
            this.nudBmpMin = new System.Windows.Forms.NumericUpDown();
            this.nudBmpMax = new System.Windows.Forms.NumericUpDown();
            ((System.ComponentModel.ISupportInitialize)(this.nudNumAScansPerBScan)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.nudAScanAveragingFactor)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.nudBScanAveragingFactor)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.nudNumBScans)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.nudOutputAScanLength)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.nudInputSpectraLength)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.nudDeviceIndex)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.pictureBox1)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.nudBmpMin)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.nudBmpMax)).BeginInit();
            this.SuspendLayout();
            // 
            // btnLoadRawData
            // 
            this.btnLoadRawData.Location = new System.Drawing.Point(12, 12);
            this.btnLoadRawData.Name = "btnLoadRawData";
            this.btnLoadRawData.Size = new System.Drawing.Size(343, 23);
            this.btnLoadRawData.TabIndex = 0;
            this.btnLoadRawData.Text = "Set Raw OCT Data File";
            this.btnLoadRawData.UseVisualStyleBackColor = true;
            this.btnLoadRawData.Click += new System.EventHandler(this.btnLoadRawData_Click);
            // 
            // dlgOpenFile
            // 
            this.dlgOpenFile.FileName = "RawData.bin";
            // 
            // lblFileSize
            // 
            this.lblFileSize.BorderStyle = System.Windows.Forms.BorderStyle.Fixed3D;
            this.lblFileSize.Location = new System.Drawing.Point(361, 12);
            this.lblFileSize.Name = "lblFileSize";
            this.lblFileSize.Size = new System.Drawing.Size(360, 23);
            this.lblFileSize.TabIndex = 1;
            this.lblFileSize.TextAlign = System.Drawing.ContentAlignment.MiddleCenter;
            // 
            // label1
            // 
            this.label1.Location = new System.Drawing.Point(12, 150);
            this.label1.Name = "label1";
            this.label1.Size = new System.Drawing.Size(123, 24);
            this.label1.TabIndex = 2;
            this.label1.Text = "Input Spectra Length";
            this.label1.TextAlign = System.Drawing.ContentAlignment.MiddleLeft;
            // 
            // label2
            // 
            this.label2.Location = new System.Drawing.Point(12, 174);
            this.label2.Name = "label2";
            this.label2.Size = new System.Drawing.Size(123, 24);
            this.label2.TabIndex = 3;
            this.label2.Text = "Output A-Scan Length";
            this.label2.TextAlign = System.Drawing.ContentAlignment.MiddleLeft;
            // 
            // label3
            // 
            this.label3.Location = new System.Drawing.Point(12, 198);
            this.label3.Name = "label3";
            this.label3.Size = new System.Drawing.Size(123, 24);
            this.label3.TabIndex = 4;
            this.label3.Text = "Number of B-Scans";
            this.label3.TextAlign = System.Drawing.ContentAlignment.MiddleLeft;
            // 
            // label4
            // 
            this.label4.Location = new System.Drawing.Point(12, 222);
            this.label4.Name = "label4";
            this.label4.Size = new System.Drawing.Size(161, 24);
            this.label4.TabIndex = 5;
            this.label4.Text = "Number of A-Scans Per B-Scan";
            this.label4.TextAlign = System.Drawing.ContentAlignment.MiddleLeft;
            // 
            // label5
            // 
            this.label5.Location = new System.Drawing.Point(12, 246);
            this.label5.Name = "label5";
            this.label5.Size = new System.Drawing.Size(161, 24);
            this.label5.TabIndex = 6;
            this.label5.Text = "A-Scan Averaging Factor";
            this.label5.TextAlign = System.Drawing.ContentAlignment.MiddleLeft;
            // 
            // label6
            // 
            this.label6.Location = new System.Drawing.Point(12, 270);
            this.label6.Name = "label6";
            this.label6.Size = new System.Drawing.Size(161, 24);
            this.label6.TabIndex = 7;
            this.label6.Text = "B-Scan Averaging Factor";
            this.label6.TextAlign = System.Drawing.ContentAlignment.MiddleLeft;
            // 
            // btnSetKernelSource
            // 
            this.btnSetKernelSource.Location = new System.Drawing.Point(12, 128);
            this.btnSetKernelSource.Name = "btnSetKernelSource";
            this.btnSetKernelSource.Size = new System.Drawing.Size(343, 23);
            this.btnSetKernelSource.TabIndex = 8;
            this.btnSetKernelSource.Text = "Kernel Source";
            this.btnSetKernelSource.UseVisualStyleBackColor = true;
            this.btnSetKernelSource.Click += new System.EventHandler(this.btnSetKernelSource_Click);
            // 
            // lblKernelSourcePath
            // 
            this.lblKernelSourcePath.BorderStyle = System.Windows.Forms.BorderStyle.Fixed3D;
            this.lblKernelSourcePath.Location = new System.Drawing.Point(361, 128);
            this.lblKernelSourcePath.Name = "lblKernelSourcePath";
            this.lblKernelSourcePath.Size = new System.Drawing.Size(360, 190);
            this.lblKernelSourcePath.TabIndex = 9;
            this.lblKernelSourcePath.TextAlign = System.Drawing.ContentAlignment.MiddleCenter;
            // 
            // nudNumAScansPerBScan
            // 
            this.nudNumAScansPerBScan.Location = new System.Drawing.Point(179, 226);
            this.nudNumAScansPerBScan.Maximum = new decimal(new int[] {
            100000,
            0,
            0,
            0});
            this.nudNumAScansPerBScan.Minimum = new decimal(new int[] {
            1,
            0,
            0,
            0});
            this.nudNumAScansPerBScan.Name = "nudNumAScansPerBScan";
            this.nudNumAScansPerBScan.Size = new System.Drawing.Size(176, 20);
            this.nudNumAScansPerBScan.TabIndex = 10;
            this.nudNumAScansPerBScan.Value = new decimal(new int[] {
            500,
            0,
            0,
            0});
            this.nudNumAScansPerBScan.ValueChanged += new System.EventHandler(this.nudNumAScansPerBScan_ValueChanged);
            // 
            // nudAScanAveragingFactor
            // 
            this.nudAScanAveragingFactor.Location = new System.Drawing.Point(179, 250);
            this.nudAScanAveragingFactor.Minimum = new decimal(new int[] {
            1,
            0,
            0,
            0});
            this.nudAScanAveragingFactor.Name = "nudAScanAveragingFactor";
            this.nudAScanAveragingFactor.Size = new System.Drawing.Size(176, 20);
            this.nudAScanAveragingFactor.TabIndex = 11;
            this.nudAScanAveragingFactor.Value = new decimal(new int[] {
            1,
            0,
            0,
            0});
            this.nudAScanAveragingFactor.ValueChanged += new System.EventHandler(this.nudAScanAveragingFactor_ValueChanged);
            // 
            // nudBScanAveragingFactor
            // 
            this.nudBScanAveragingFactor.Location = new System.Drawing.Point(179, 274);
            this.nudBScanAveragingFactor.Minimum = new decimal(new int[] {
            1,
            0,
            0,
            0});
            this.nudBScanAveragingFactor.Name = "nudBScanAveragingFactor";
            this.nudBScanAveragingFactor.Size = new System.Drawing.Size(176, 20);
            this.nudBScanAveragingFactor.TabIndex = 12;
            this.nudBScanAveragingFactor.Value = new decimal(new int[] {
            1,
            0,
            0,
            0});
            this.nudBScanAveragingFactor.ValueChanged += new System.EventHandler(this.nudBScanAveragingFactor_ValueChanged);
            // 
            // nudNumBScans
            // 
            this.nudNumBScans.Location = new System.Drawing.Point(179, 202);
            this.nudNumBScans.Maximum = new decimal(new int[] {
            10000,
            0,
            0,
            0});
            this.nudNumBScans.Minimum = new decimal(new int[] {
            1,
            0,
            0,
            0});
            this.nudNumBScans.Name = "nudNumBScans";
            this.nudNumBScans.Size = new System.Drawing.Size(176, 20);
            this.nudNumBScans.TabIndex = 13;
            this.nudNumBScans.Value = new decimal(new int[] {
            10,
            0,
            0,
            0});
            this.nudNumBScans.ValueChanged += new System.EventHandler(this.nudNumBScans_ValueChanged);
            // 
            // nudOutputAScanLength
            // 
            this.nudOutputAScanLength.Location = new System.Drawing.Point(179, 178);
            this.nudOutputAScanLength.Maximum = new decimal(new int[] {
            4096,
            0,
            0,
            0});
            this.nudOutputAScanLength.Minimum = new decimal(new int[] {
            32,
            0,
            0,
            0});
            this.nudOutputAScanLength.Name = "nudOutputAScanLength";
            this.nudOutputAScanLength.Size = new System.Drawing.Size(176, 20);
            this.nudOutputAScanLength.TabIndex = 14;
            this.nudOutputAScanLength.Value = new decimal(new int[] {
            1024,
            0,
            0,
            0});
            this.nudOutputAScanLength.ValueChanged += new System.EventHandler(this.nudOutputAScanLength_ValueChanged);
            // 
            // nudInputSpectraLength
            // 
            this.nudInputSpectraLength.Location = new System.Drawing.Point(179, 154);
            this.nudInputSpectraLength.Maximum = new decimal(new int[] {
            2048,
            0,
            0,
            0});
            this.nudInputSpectraLength.Minimum = new decimal(new int[] {
            64,
            0,
            0,
            0});
            this.nudInputSpectraLength.Name = "nudInputSpectraLength";
            this.nudInputSpectraLength.Size = new System.Drawing.Size(176, 20);
            this.nudInputSpectraLength.TabIndex = 15;
            this.nudInputSpectraLength.Value = new decimal(new int[] {
            1024,
            0,
            0,
            0});
            this.nudInputSpectraLength.ValueChanged += new System.EventHandler(this.nudInputSpectraLength_ValueChanged);
            // 
            // lblExpectedFileLength
            // 
            this.lblExpectedFileLength.Location = new System.Drawing.Point(15, 323);
            this.lblExpectedFileLength.Name = "lblExpectedFileLength";
            this.lblExpectedFileLength.Size = new System.Drawing.Size(706, 23);
            this.lblExpectedFileLength.TabIndex = 16;
            this.lblExpectedFileLength.Text = "Expected File Length:";
            this.lblExpectedFileLength.TextAlign = System.Drawing.ContentAlignment.MiddleCenter;
            // 
            // btnProcess
            // 
            this.btnProcess.Location = new System.Drawing.Point(12, 387);
            this.btnProcess.Name = "btnProcess";
            this.btnProcess.Size = new System.Drawing.Size(155, 41);
            this.btnProcess.TabIndex = 17;
            this.btnProcess.Text = "Process on GPU";
            this.btnProcess.UseVisualStyleBackColor = true;
            this.btnProcess.Click += new System.EventHandler(this.btnProcess_Click);
            // 
            // nudDeviceIndex
            // 
            this.nudDeviceIndex.Location = new System.Drawing.Point(179, 298);
            this.nudDeviceIndex.Maximum = new decimal(new int[] {
            10,
            0,
            0,
            0});
            this.nudDeviceIndex.Name = "nudDeviceIndex";
            this.nudDeviceIndex.Size = new System.Drawing.Size(176, 20);
            this.nudDeviceIndex.TabIndex = 19;
            // 
            // label7
            // 
            this.label7.Location = new System.Drawing.Point(12, 294);
            this.label7.Name = "label7";
            this.label7.Size = new System.Drawing.Size(161, 24);
            this.label7.TabIndex = 18;
            this.label7.Text = "GPU Device Number";
            this.label7.TextAlign = System.Drawing.ContentAlignment.MiddleLeft;
            // 
            // lblResamplingTable
            // 
            this.lblResamplingTable.BorderStyle = System.Windows.Forms.BorderStyle.Fixed3D;
            this.lblResamplingTable.Location = new System.Drawing.Point(361, 41);
            this.lblResamplingTable.Name = "lblResamplingTable";
            this.lblResamplingTable.Size = new System.Drawing.Size(360, 23);
            this.lblResamplingTable.TabIndex = 21;
            this.lblResamplingTable.TextAlign = System.Drawing.ContentAlignment.MiddleCenter;
            // 
            // btnLoadResamplingTable
            // 
            this.btnLoadResamplingTable.Location = new System.Drawing.Point(12, 41);
            this.btnLoadResamplingTable.Name = "btnLoadResamplingTable";
            this.btnLoadResamplingTable.Size = new System.Drawing.Size(343, 23);
            this.btnLoadResamplingTable.TabIndex = 20;
            this.btnLoadResamplingTable.Text = "Load Resampling Table";
            this.btnLoadResamplingTable.UseVisualStyleBackColor = true;
            this.btnLoadResamplingTable.Click += new System.EventHandler(this.btnLoadResamplingTable_Click);
            // 
            // lblReferenceSpectrum
            // 
            this.lblReferenceSpectrum.BorderStyle = System.Windows.Forms.BorderStyle.Fixed3D;
            this.lblReferenceSpectrum.Location = new System.Drawing.Point(361, 70);
            this.lblReferenceSpectrum.Name = "lblReferenceSpectrum";
            this.lblReferenceSpectrum.Size = new System.Drawing.Size(360, 23);
            this.lblReferenceSpectrum.TabIndex = 23;
            this.lblReferenceSpectrum.TextAlign = System.Drawing.ContentAlignment.MiddleCenter;
            // 
            // btnLoadReferenceSpectrum
            // 
            this.btnLoadReferenceSpectrum.Location = new System.Drawing.Point(12, 70);
            this.btnLoadReferenceSpectrum.Name = "btnLoadReferenceSpectrum";
            this.btnLoadReferenceSpectrum.Size = new System.Drawing.Size(343, 23);
            this.btnLoadReferenceSpectrum.TabIndex = 22;
            this.btnLoadReferenceSpectrum.Text = "Load Reference Spectrum";
            this.btnLoadReferenceSpectrum.UseVisualStyleBackColor = true;
            this.btnLoadReferenceSpectrum.Click += new System.EventHandler(this.btnLoadReferenceSpectrum_Click);
            // 
            // lblReferenceAScan
            // 
            this.lblReferenceAScan.BorderStyle = System.Windows.Forms.BorderStyle.Fixed3D;
            this.lblReferenceAScan.Location = new System.Drawing.Point(361, 99);
            this.lblReferenceAScan.Name = "lblReferenceAScan";
            this.lblReferenceAScan.Size = new System.Drawing.Size(360, 23);
            this.lblReferenceAScan.TabIndex = 25;
            this.lblReferenceAScan.TextAlign = System.Drawing.ContentAlignment.MiddleCenter;
            // 
            // btnLoadReferenceAScan
            // 
            this.btnLoadReferenceAScan.Location = new System.Drawing.Point(12, 99);
            this.btnLoadReferenceAScan.Name = "btnLoadReferenceAScan";
            this.btnLoadReferenceAScan.Size = new System.Drawing.Size(343, 23);
            this.btnLoadReferenceAScan.TabIndex = 24;
            this.btnLoadReferenceAScan.Text = "Load Reference A-Scan";
            this.btnLoadReferenceAScan.UseVisualStyleBackColor = true;
            this.btnLoadReferenceAScan.Click += new System.EventHandler(this.btnLoadReferenceAScan_Click);
            // 
            // txtResamplingTable
            // 
            this.txtResamplingTable.Location = new System.Drawing.Point(727, 12);
            this.txtResamplingTable.Multiline = true;
            this.txtResamplingTable.Name = "txtResamplingTable";
            this.txtResamplingTable.ReadOnly = true;
            this.txtResamplingTable.ScrollBars = System.Windows.Forms.ScrollBars.Both;
            this.txtResamplingTable.Size = new System.Drawing.Size(147, 360);
            this.txtResamplingTable.TabIndex = 26;
            // 
            // txtReferenceSpectrum
            // 
            this.txtReferenceSpectrum.Location = new System.Drawing.Point(880, 12);
            this.txtReferenceSpectrum.Multiline = true;
            this.txtReferenceSpectrum.Name = "txtReferenceSpectrum";
            this.txtReferenceSpectrum.ReadOnly = true;
            this.txtReferenceSpectrum.ScrollBars = System.Windows.Forms.ScrollBars.Both;
            this.txtReferenceSpectrum.Size = new System.Drawing.Size(147, 360);
            this.txtReferenceSpectrum.TabIndex = 27;
            // 
            // txtReferenceAScan
            // 
            this.txtReferenceAScan.Location = new System.Drawing.Point(1033, 12);
            this.txtReferenceAScan.Multiline = true;
            this.txtReferenceAScan.Name = "txtReferenceAScan";
            this.txtReferenceAScan.ReadOnly = true;
            this.txtReferenceAScan.ScrollBars = System.Windows.Forms.ScrollBars.Both;
            this.txtReferenceAScan.Size = new System.Drawing.Size(147, 360);
            this.txtReferenceAScan.TabIndex = 28;
            // 
            // progLoadRaw
            // 
            this.progLoadRaw.Location = new System.Drawing.Point(361, 349);
            this.progLoadRaw.Name = "progLoadRaw";
            this.progLoadRaw.Size = new System.Drawing.Size(360, 23);
            this.progLoadRaw.TabIndex = 29;
            // 
            // btnLoad
            // 
            this.btnLoad.Location = new System.Drawing.Point(12, 349);
            this.btnLoad.Name = "btnLoad";
            this.btnLoad.Size = new System.Drawing.Size(343, 23);
            this.btnLoad.TabIndex = 30;
            this.btnLoad.Text = "Load Raw Spectra";
            this.btnLoad.UseVisualStyleBackColor = true;
            this.btnLoad.Click += new System.EventHandler(this.btnLoad_Click);
            // 
            // btnBScans
            // 
            this.btnBScans.Location = new System.Drawing.Point(179, 387);
            this.btnBScans.Name = "btnBScans";
            this.btnBScans.Size = new System.Drawing.Size(129, 41);
            this.btnBScans.TabIndex = 31;
            this.btnBScans.Text = "Output B-Scans";
            this.btnBScans.UseVisualStyleBackColor = true;
            this.btnBScans.Click += new System.EventHandler(this.btnBScans_Click);
            // 
            // lblStatus
            // 
            this.lblStatus.BorderStyle = System.Windows.Forms.BorderStyle.FixedSingle;
            this.lblStatus.Location = new System.Drawing.Point(361, 387);
            this.lblStatus.Name = "lblStatus";
            this.lblStatus.Size = new System.Drawing.Size(360, 41);
            this.lblStatus.TabIndex = 32;
            this.lblStatus.TextAlign = System.Drawing.ContentAlignment.MiddleCenter;
            // 
            // pictureBox1
            // 
            this.pictureBox1.BorderStyle = System.Windows.Forms.BorderStyle.FixedSingle;
            this.pictureBox1.Location = new System.Drawing.Point(727, 387);
            this.pictureBox1.Name = "pictureBox1";
            this.pictureBox1.Size = new System.Drawing.Size(453, 264);
            this.pictureBox1.SizeMode = System.Windows.Forms.PictureBoxSizeMode.StretchImage;
            this.pictureBox1.TabIndex = 33;
            this.pictureBox1.TabStop = false;
            // 
            // label8
            // 
            this.label8.AutoSize = true;
            this.label8.Location = new System.Drawing.Point(12, 443);
            this.label8.Name = "label8";
            this.label8.Size = new System.Drawing.Size(62, 13);
            this.label8.TabIndex = 34;
            this.label8.Text = "Bitmap Min.";
            // 
            // label11
            // 
            this.label11.AutoSize = true;
            this.label11.Location = new System.Drawing.Point(15, 465);
            this.label11.Name = "label11";
            this.label11.Size = new System.Drawing.Size(65, 13);
            this.label11.TabIndex = 37;
            this.label11.Text = "Bitmap Max.";
            this.label11.TextAlign = System.Drawing.ContentAlignment.MiddleLeft;
            // 
            // nudBmpMin
            // 
            this.nudBmpMin.Location = new System.Drawing.Point(89, 441);
            this.nudBmpMin.Maximum = new decimal(new int[] {
            1000,
            0,
            0,
            0});
            this.nudBmpMin.Minimum = new decimal(new int[] {
            1000,
            0,
            0,
            -2147483648});
            this.nudBmpMin.Name = "nudBmpMin";
            this.nudBmpMin.Size = new System.Drawing.Size(120, 20);
            this.nudBmpMin.TabIndex = 38;
            this.nudBmpMin.Value = new decimal(new int[] {
            155,
            0,
            0,
            -2147483648});
            // 
            // nudBmpMax
            // 
            this.nudBmpMax.Location = new System.Drawing.Point(89, 463);
            this.nudBmpMax.Maximum = new decimal(new int[] {
            1000,
            0,
            0,
            0});
            this.nudBmpMax.Minimum = new decimal(new int[] {
            1000,
            0,
            0,
            -2147483648});
            this.nudBmpMax.Name = "nudBmpMax";
            this.nudBmpMax.Size = new System.Drawing.Size(120, 20);
            this.nudBmpMax.TabIndex = 39;
            this.nudBmpMax.Value = new decimal(new int[] {
            55,
            0,
            0,
            -2147483648});
            // 
            // Form1
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(1208, 663);
            this.Controls.Add(this.nudBmpMax);
            this.Controls.Add(this.nudBmpMin);
            this.Controls.Add(this.label11);
            this.Controls.Add(this.label8);
            this.Controls.Add(this.pictureBox1);
            this.Controls.Add(this.lblStatus);
            this.Controls.Add(this.btnBScans);
            this.Controls.Add(this.btnLoad);
            this.Controls.Add(this.progLoadRaw);
            this.Controls.Add(this.txtReferenceAScan);
            this.Controls.Add(this.txtReferenceSpectrum);
            this.Controls.Add(this.txtResamplingTable);
            this.Controls.Add(this.lblReferenceAScan);
            this.Controls.Add(this.btnLoadReferenceAScan);
            this.Controls.Add(this.lblReferenceSpectrum);
            this.Controls.Add(this.btnLoadReferenceSpectrum);
            this.Controls.Add(this.lblResamplingTable);
            this.Controls.Add(this.btnLoadResamplingTable);
            this.Controls.Add(this.nudDeviceIndex);
            this.Controls.Add(this.label7);
            this.Controls.Add(this.btnProcess);
            this.Controls.Add(this.lblExpectedFileLength);
            this.Controls.Add(this.nudInputSpectraLength);
            this.Controls.Add(this.nudOutputAScanLength);
            this.Controls.Add(this.nudNumBScans);
            this.Controls.Add(this.nudBScanAveragingFactor);
            this.Controls.Add(this.nudAScanAveragingFactor);
            this.Controls.Add(this.nudNumAScansPerBScan);
            this.Controls.Add(this.lblKernelSourcePath);
            this.Controls.Add(this.btnSetKernelSource);
            this.Controls.Add(this.label6);
            this.Controls.Add(this.label5);
            this.Controls.Add(this.label4);
            this.Controls.Add(this.label3);
            this.Controls.Add(this.label2);
            this.Controls.Add(this.label1);
            this.Controls.Add(this.lblFileSize);
            this.Controls.Add(this.btnLoadRawData);
            this.Name = "Form1";
            this.Text = "OpenCL Real-Time OCT Processing Library Offline Testing";
            ((System.ComponentModel.ISupportInitialize)(this.nudNumAScansPerBScan)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.nudAScanAveragingFactor)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.nudBScanAveragingFactor)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.nudNumBScans)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.nudOutputAScanLength)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.nudInputSpectraLength)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.nudDeviceIndex)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.pictureBox1)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.nudBmpMin)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.nudBmpMax)).EndInit();
            this.ResumeLayout(false);
            this.PerformLayout();

        }

        #endregion

        private System.Windows.Forms.Button btnLoadRawData;
        private System.Windows.Forms.OpenFileDialog dlgOpenFile;
        private System.Windows.Forms.Label lblFileSize;
        private System.Windows.Forms.Label label1;
        private System.Windows.Forms.Label label2;
        private System.Windows.Forms.Label label3;
        private System.Windows.Forms.Label label4;
        private System.Windows.Forms.Label label5;
        private System.Windows.Forms.Label label6;
        private System.Windows.Forms.Button btnSetKernelSource;
        private System.Windows.Forms.Label lblKernelSourcePath;
        private System.Windows.Forms.NumericUpDown nudNumAScansPerBScan;
        private System.Windows.Forms.NumericUpDown nudAScanAveragingFactor;
        private System.Windows.Forms.NumericUpDown nudBScanAveragingFactor;
        private System.Windows.Forms.NumericUpDown nudNumBScans;
        private System.Windows.Forms.NumericUpDown nudOutputAScanLength;
        private System.Windows.Forms.NumericUpDown nudInputSpectraLength;
        private System.Windows.Forms.Label lblExpectedFileLength;
        private System.Windows.Forms.Button btnProcess;
        private System.Windows.Forms.NumericUpDown nudDeviceIndex;
        private System.Windows.Forms.Label label7;
        private System.Windows.Forms.Label lblResamplingTable;
        private System.Windows.Forms.Button btnLoadResamplingTable;
        private System.Windows.Forms.Label lblReferenceSpectrum;
        private System.Windows.Forms.Button btnLoadReferenceSpectrum;
        private System.Windows.Forms.Label lblReferenceAScan;
        private System.Windows.Forms.Button btnLoadReferenceAScan;
        private System.Windows.Forms.TextBox txtResamplingTable;
        private System.Windows.Forms.TextBox txtReferenceSpectrum;
        private System.Windows.Forms.TextBox txtReferenceAScan;
        private System.Windows.Forms.ProgressBar progLoadRaw;
        private System.Windows.Forms.Button btnLoad;
        private System.Windows.Forms.Button btnBScans;
        private System.Windows.Forms.Label lblStatus;
        private System.Windows.Forms.PictureBox pictureBox1;
        private System.Windows.Forms.Label label8;
        private System.Windows.Forms.Label label11;
        private System.Windows.Forms.NumericUpDown nudBmpMin;
        private System.Windows.Forms.NumericUpDown nudBmpMax;
    }
}

