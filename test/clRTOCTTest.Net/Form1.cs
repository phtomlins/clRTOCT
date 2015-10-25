using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Windows.Forms;
using System.IO;
using OCTProcessingLibrary;
using System.Diagnostics;



namespace clRTOCTTest
{
    public partial class Form1 : Form
    {

        #region Class Variables

        FileStream pFs;
        BinaryReader pBr;
        long pFileLength;
        string pKernelPath="C:\\Users\\OCT\\Google Drive\\OCT Control Software - octView\\clRTOCT\\devel\\current\\src\\kernels\\octProcessingKernels.cl";
        string pSpectraFileName = "";

        const long BYTES_PER_DATA_POINT  = 2;

        float[] pHostResamplingTable;
        float[] pHostReferenceSpectrum;
        float[] pHostReferenceAScan;

        byte[] pRawSpectra;
        short[,] pHostSpectra;

        float[,] pHostLogEnvBScan;


        #endregion

        public Form1()
        {
            InitializeComponent();
            UpdateKernelPath(pKernelPath);
            UpdateExpectedFileLength();
        }

        private void btnLoadRawData_Click(object sender, EventArgs e)
        {
            if (dlgOpenFile.ShowDialog() == DialogResult.OK)
            {
                pSpectraFileName = dlgOpenFile.FileName;
                pFs = new FileStream(pSpectraFileName,FileMode.Open);
                //
                // Attempt to determine the file length
                //
                pFileLength = pFs.Length / BYTES_PER_DATA_POINT;
                pFs.Close();
             
                this.lblFileSize.Text = "There are " + pFileLength.ToString() + " elements in this file.";
                //
                // Try to load the resampling table, reference spectrum and reference a-scan
                //
                string dir = new FileInfo(dlgOpenFile.FileName).Directory.FullName;
                LoadResamplingTable(dir + "\\resamplingTable.csv");
                LoadReferenceSpectrum(dir + "\\referenceSpectrum.csv");
                LoadReferenceAScan(dir + "\\referenceAScan.csv");
            }
        }

        private void btnSetKernelSource_Click(object sender, EventArgs e)
        {
            if (dlgOpenFile.ShowDialog() == DialogResult.OK)
            {
                UpdateKernelPath(dlgOpenFile.FileName);
            }

            
        }

        private void UpdateKernelPath(string kernelPath)
        {
            pKernelPath = kernelPath;
            this.lblKernelSourcePath.Text = pKernelPath;
        }

        private void UpdateExpectedFileLength()
        {
            long expectedFileLength = (long)this.nudInputSpectraLength.Value *
                (long)this.nudNumAScansPerBScan.Value *
                (long)this.nudNumBScans.Value *
                (long)this.nudAScanAveragingFactor.Value *
                (long)this.nudBScanAveragingFactor.Value;

            this.lblExpectedFileLength.Text = "Expected File Length: " + expectedFileLength.ToString();
        }

        private void nudInputSpectraLength_ValueChanged(object sender, EventArgs e)
        {
            UpdateExpectedFileLength();
        }

        private void nudOutputAScanLength_ValueChanged(object sender, EventArgs e)
        {
            UpdateExpectedFileLength();

        }

        private void nudNumBScans_ValueChanged(object sender, EventArgs e)
        {
            UpdateExpectedFileLength();

        }

        private void nudNumAScansPerBScan_ValueChanged(object sender, EventArgs e)
        {
            UpdateExpectedFileLength();

        }

        private void nudAScanAveragingFactor_ValueChanged(object sender, EventArgs e)
        {
            UpdateExpectedFileLength();

        }

        private void nudBScanAveragingFactor_ValueChanged(object sender, EventArgs e)
        {
            UpdateExpectedFileLength();

        }

        private void btnProcess_Click(object sender, EventArgs e)
        {
            this.Enabled = false;
            //
            // Initialise the GPU
            //
            clRTOCT clrtoct;
            try
            {
                DateTime end;

                lblStatus.Text = "Initialising GPU...";
                Application.DoEvents();

                clrtoct = new clRTOCT(
                    (int)nudDeviceIndex.Value,
                    (int)nudInputSpectraLength.Value,
                    (int)nudOutputAScanLength.Value,
                    (int)nudNumBScans.Value,
                    (int)nudNumAScansPerBScan.Value,
                    (int)nudAScanAveragingFactor.Value,
                    (int)nudBScanAveragingFactor.Value,
                    pHostResamplingTable,
                    pHostReferenceSpectrum,
                    pHostReferenceAScan,
                    pKernelPath,
                    2,
                    2,
                    "-cl-fast-relaxed-math -cl-mad-enable"
                    );
                //
                // Reformat the input spectra into the correct 2D array
                //

                DateTime start = DateTime.Now;

                // Process
                lblStatus.Text = "Pre-Processing...";
                Application.DoEvents();
                clrtoct.PreProcessSpectra(pHostSpectra, clRTOCT.eWINDOW_TYPE.BLACKMAN);

               // float[,] preProcessed = new float[pHostSpectra.GetLength(0),pHostSpectra.GetLength(1)*2];
               // clrtoct.GetPreProcessed(preProcessed);

                lblStatus.Text = "Fourier transform...";
                Application.DoEvents();
                clrtoct.InverseFourierTransform();

                lblStatus.Text = "Post-processing...";
                Application.DoEvents();

                clrtoct.PostProcessOCTSignal((float)nudBmpMin.Value,(float)nudBmpMax.Value);

                lblStatus.Text = "Copying back to host...";
                Application.DoEvents();
                clrtoct.GetLogEnvelope(pHostLogEnvBScan);

                this.pictureBox1.Image = clrtoct.getBScanBitmap(1);
                byte[] bscanBmps = new byte[clrtoct.TotalBScans * clrtoct.BScanBitmapHeight * clrtoct.BScanStride];
                clrtoct.GetBScanBitmaps(bscanBmps);
                end = DateTime.Now;



                lblStatus.Text = String.Format("Done in {0:0.000} s",end.Subtract(start).TotalSeconds);
                Application.DoEvents();

                

                clrtoct.Dispose();// = null;    // Force the clRTOCT object to dispose
                
            }
            catch (Exception ex)
            {

                MessageBox.Show(ex.Message);
                clrtoct = null;

            }

            this.Enabled = true;
        }

        private float[] ReadCSV(string filename, int col)
        {
            //
            // Read a column of numbers from a csv file
            //
            StreamReader sr = new StreamReader(filename);
            string line;
            string[] elements;
            float[] data = new float[1000000];
            int i=0;



            while (!sr.EndOfStream)
            {
                line = sr.ReadLine();
                elements = line.Split(',');
                data[i] = Convert.ToSingle(elements[col]);
                i++;
            }

            Array.Resize(ref data, i);

            sr.Close();
            return data;

        }


        private float[] LoadTable(string filename, ref TextBox textbox, ref Label label)
        {
            try
            {
                string txt = "";
                //this.txtResamplingTable.Text = "";
                float[] col0 = ReadCSV(filename, 0);
                float[] col1 = ReadCSV(filename, 1);

                for (int i = 0; i < col0.Length; i++)
                {
                    txt += String.Format("{0:0.000},{1:0.000}" + Environment.NewLine, col0[i], col1[i]);
                }
                textbox.Text = txt;
                label.Text = "Success.";
                return col1;
            }
            catch (Exception ex)
            {
                label.Text = "Failed!";
                textbox.Text = ex.Message;
                return null;
            }
        }
        

        private void LoadResamplingTable(string filename)
        {
            pHostResamplingTable = this.LoadTable(filename, ref this.txtResamplingTable, ref this.lblResamplingTable);
        }

        private void LoadReferenceSpectrum(string filename)
        {
            pHostReferenceSpectrum = this.LoadTable(filename, ref this.txtReferenceSpectrum, ref this.lblReferenceSpectrum);
        }

        private void LoadReferenceAScan(string filename)
        {
            pHostReferenceAScan = this.LoadTable(filename, ref this.txtReferenceAScan, ref this.lblReferenceAScan);
        }

        private void btnLoadResamplingTable_Click(object sender, EventArgs e)
        {
            if (dlgOpenFile.ShowDialog() == DialogResult.OK)
            {
                LoadResamplingTable(dlgOpenFile.FileName);
            }
        }

        private void btnLoadReferenceSpectrum_Click(object sender, EventArgs e)
        {
            if (dlgOpenFile.ShowDialog() == DialogResult.OK)
            {
                LoadReferenceSpectrum(dlgOpenFile.FileName);
            }

        }

        private void btnLoadReferenceAScan_Click(object sender, EventArgs e)
        {
            if (dlgOpenFile.ShowDialog() == DialogResult.OK)
            {
                LoadReferenceAScan(dlgOpenFile.FileName);
            }

        }

        private void btnLoad_Click(object sender, EventArgs e)
        {

            //
            // Load the whole data file into memory and initialise the GPU
            //
            this.progLoadRaw.Minimum = 0;
            this.progLoadRaw.Value = 0;
            this.progLoadRaw.Maximum = (int)pFileLength;
            //
            // Read in the raw spectra
            //
            pFs = new FileStream(pSpectraFileName, FileMode.Open);
            pFs.Seek(0, SeekOrigin.Begin);  // got to the beginning of the file
            pBr = new BinaryReader(pFs);

            int totalBScans = (int)(nudNumBScans.Value * nudBScanAveragingFactor.Value);
            int totalAScansPerBScan = (int)(nudNumAScansPerBScan.Value * nudAScanAveragingFactor.Value);
            int inputSpectraLength = (int)nudInputSpectraLength.Value;
            int outputAScanLength = (int)nudOutputAScanLength.Value;

            pHostSpectra = new short[totalBScans, inputSpectraLength * totalAScansPerBScan];
            pHostLogEnvBScan = new float[totalBScans, outputAScanLength * totalAScansPerBScan];

            int progress = 0;

            this.Enabled = false;
            for (int b = 0; b < totalBScans; b++)
            {
                for (int a = 0; a < totalAScansPerBScan; a++)
                {
                    for (int p = 0; p < inputSpectraLength; p++)
                    {
                        pHostSpectra[b, a * inputSpectraLength + p] = pBr.ReadInt16();
                        progress++;
                    }
                }
                Debug.WriteLine("Loading B-Scan: " + b);
                this.progLoadRaw.Value = progress;
                Application.DoEvents();
            }
            this.Enabled = true;
            pBr.Close();
            pFs.Close();


            //

        }

        private void btnBScans_Click(object sender, EventArgs e)
        {
            float[,] img = new float[(int)(nudNumAScansPerBScan.Value * nudAScanAveragingFactor.Value), (int)nudOutputAScanLength.Value];
            int totalBScans = (int)(nudNumBScans.Value * nudBScanAveragingFactor.Value);
            int totalAScansPerBScan = (int)(nudNumAScansPerBScan.Value * nudAScanAveragingFactor.Value);
            int outputAScanLength = (int)nudOutputAScanLength.Value;
            int inputSpectraLength = (int)nudInputSpectraLength.Value;

            for (int i = 0; i < totalBScans; i++)
            {
                //
                // Copy the bscan out of the image array
                //
                for (int a = 0; a < totalAScansPerBScan; a++)
                {
                    for (int p = 0; p < outputAScanLength; p++)
                    {
                        img[a, p] = pHostLogEnvBScan[i, a * outputAScanLength + p];
                    }
                }
                clRTOCT.SaveToBitmap("test" + i + ".bmp", img, 512, -1550, -550);
            }
            lblStatus.Text="Bitmaps Saved.";
        }
    }
}
