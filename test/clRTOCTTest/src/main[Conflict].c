
/*
 Offline OCT Processing Example - Calls clRTOCT library
 Dr P Tomlins, School of Medicine and Dentisty, Queen Mary University of London
 
 
 Copyright (c) 2014, Queen Mary University of London
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright
 notice, this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright
 notice, this list of conditions and the following disclaimer in the
 documentation and/or other materials provided with the distribution.
 * Neither the name of Queen Mary University of London nor the
 names of its contributors may be used to endorse or promote products
 derived from this software without specific prior written permission.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL QUEEN MARY UNIVERSITY OF LONDON BE LIABLE FOR ANY
 DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#define MAX_CSV_FILE_LENGTH 1000000
//#define TRUE 1
//#define FALSE 0
#ifdef _WIN32
    #include "timing.h"
    #define GPU_DEVICE_INDEX 0          // GPU device to be used by OpenCL
    #define PLATFORM 0

#elif __APPLE__
    #define GPU_DEVICE_INDEX 1          // GPU device to be used by OpenCL
    #define PLATFORM 1
#elif
    #define PLATFORM -1
#endif

#include "main.h"
//#include <sys/time.h>
//#include "clfft_test.h"

int main(int argc, char **argv)
{

	//
	// Create some input data
	//
	/*
	DLL clOCTInit(
	size_t inputSpectraLength,
	size_t outputAScanLength,
	size_t numBScans,
	size_t numAScansPerBScan,
	size_t ascanAveragingFactor,
	size_t bscanAveragingFactor,
	float* hostResamplingTable,		// Resampling table
	float* hostReferenceSpectrum,	// Reference spectrum to be subtracted (or maybe divided!!) from each spectrum
	float* hostReferenceAScan		// A reference background A-Scan to be subtracted.. maybe log domain!!
	)
	*/
    cl_uint clDeviceIndex = GPU_DEVICE_INDEX;

    
    unsigned int inputLen; // = 1024;
    unsigned int outputLen; // = 1024;
    //
    // Set the specifics of the input data - eventually will be with be either passed from command line
    // or read from parameters file
    //
    unsigned int totalBScans = 100;
	unsigned int numBScansPerBatch = 50;           // Number of BScans to process in one go
    unsigned int numBScanProcessingIteratations = (unsigned int)floor(totalBScans / numBScansPerBatch);
	unsigned int numAScans = 500;
	unsigned int ascanAve = 1;
	unsigned int bscanAve = 1;
    //
    // Bitmap options
    //
    float bmpMinVal = -155;    // Threshold values for storing log images as bitmaps
    float bmpMaxVal = -55;
    //
    // Other options
    //
    unsigned int saveBmp = true;
    //
    // Set the OpenCL build log length - used internally when compiling kernels
    //
	size_t loglen = 1000000;
    cl_uint numAScansPerBScan;
    cl_uint totalAScans; // = numAScans * numBScans * ascanAve * bscanAve;
    cl_uint totalInputLen; // = inputLen * totalAScans;
    cl_uint totalOutputLen; // = outputLen * totalAScans;
    cl_uint totalBatchLength;    // Number of array elements in a batch of bscans

    float* resampTable; // = (float*)malloc(sizeof(float) * inputLen);
    float* refSpec; // = (float*)malloc(sizeof(float) * inputLen);
    float* refAScan; // = (float*)malloc(sizeof(float) * inputLen);

    short* inputSpectra; // = (short*)malloc(sizeof(short) * totalInputLen);
    char* buildLog; // = (char*)malloc(sizeof(char) * loglen);

    float* referenceSpectrumFileData;
    float* referenceAScanFileData;
    float* resampleTableFileData;
    float* tempArray;
    //
    // Time variables
    //
    struct timeval timeBegin, timeEnd;
    struct timeval timePreProcessBegin, timePreProcessEnd;
    struct timeval timePostProcessBegin, timePostProcessEnd;
    struct timeval timeFTProcessBegin, timeFTProcessEnd;
    struct timeval timeAllProcessBegin, timeAllProcessEnd;
    long preProcessElapsedMicroSeconds, elapsedMicroSeconds, ftProcessElapsedMicroSeconds, postProcessElapsedMicroSeconds, allProcessElapsedMicroSeconds;
    double preProcessElapsedMilliSeconds, elapsedSeconds, postProcessElapsedMilliSeconds, ftProcessElapsedMilliSeconds, allProcessElapsedMilliSeconds;
    
    //
    // Set the kernel path... this depends on the platform
    //
    //****************************************************
//	char* kernelPath = "C:\\Users\\oct\\Google Drive\\OCT Control Software - octView\\clOCTKernels\\current\\octProcessingKernels2.cl\0";
    char* kernelPath;// = "../../../../../src/kernels/octProcessingKernels.cl\0";
    //****************************************************
	//
	// Allocate arrays for the output data
	//
    float* preProcessed;// = (float*)malloc(sizeof(float) * totalOutputLen*2);
    float* cmplx;// = (float*)malloc(sizeof(float) * totalOutputLen*2);
//	float* real = (float*)malloc(sizeof(float) * totalOutputLen);
//	float* imag = (float*)malloc(sizeof(float) * totalOutputLen);
    float* linEnv;// = (float*)malloc(sizeof(float) * totalOutputLen);
    float* logEnv;// = (float*)malloc(sizeof(float) * totalOutputLen);
    float* sum;// = (float*)malloc(sizeof(float) * totalAScans);
    float* sam;// = (float*)malloc(sizeof(float) * totalAScans);
    float* ad;// = (float*)malloc(sizeof(float) * totalAScans);
    float* correlationMap;
    
	int res;

//	int repeat = 500;
	int r;
	int i;
	int j;
    




	//
	// Files for IO
	//
    FILE* inputSpecFile;// = fopen("inputSpec.bin","wb");
//    FILE* preProcessedFile;// = fopen("preProcessed.bin","wb");
//	FILE* realFile = fopen("realPart.bin","wb");
//	FILE* imagFile = fopen("imagPart.bin","wb");
//    FILE* cmplxFile;// = fopen("cmplx.bin","wb");
//    FILE* linEnvFile;// = fopen("linEnv.bin","wb");
//    FILE* logEnvFile;// = fopen("logEnv.bin","wb");



    //
    // Use the command line to specify a directory containing the OCT data to be processed
    //
 //   if (argc<4)
 //   {
 //       exit(-1);   // Quite with an error code if the argument count is less than 1
 //   }
    char* resamplingTablePath;//="/Users/phtomlins/Google Drive/OCT Data/Test/resamplingTable.csv";
    char* spectraPath;//="/Users/phtomlins/Google Drive/OCT Data/Test/Spectra.bin";
    char* referenceSpectrumPath;//="/Users/phtomlins/Google Drive/OCT Data/Test/referenceSpectrum.csv";
    char* referenceAScanPath;//="/Users/phtomlins/Google Drive/OCT Data/Test/referenceAScan.csv";
    char* rootPath;//="/Users/phtomlins/Google Drive/OCT Data/Test/";
    char bmpPath[2048];

    if (PLATFORM == 0)  // Windows
    {
        kernelPath = "C:\\Users\\OCT\\Google Drive\\OCT Control Software - octView\\clRTOCT\\devel\\current\\src\\kernels\\octProcessingKernels.cl\0";
        resamplingTablePath="C:\\Users\\OCT\\Google Drive\\OCT Data\\Test\\resamplingTable.csv";
        spectraPath="C:\\Users\\OCT\\Google Drive\\OCT Data\\Test\\Spectra.bin";
        referenceSpectrumPath="C:\\Users\\OCT\\Google Drive\\OCT Data\\Test\\referenceSpectrum.csv";
        referenceAScanPath="C:\\Users\\OCT\\Google Drive\\OCT Data\\Test\\referenceAScan.csv";
        rootPath="C:\\Users\\OCT\\Google Drive\\OCT Data\\Test\\";
    }
    else if (PLATFORM == 1) // MacOS
    {
        kernelPath = "../../../../../src/kernels/octProcessingKernels.cl\0";
        resamplingTablePath="/Users/phtomlins/Google Drive/OCT Data/Test/resamplingTable.csv";
        spectraPath="/Users/phtomlins/Google Drive/OCT Data/Test/Spectra.bin";
        referenceSpectrumPath="/Users/phtomlins/Google Drive/OCT Data/Test/referenceSpectrum.csv";
        referenceAScanPath="/Users/phtomlins/Google Drive/OCT Data/Test/referenceAScan.csv";
        rootPath="/Users/phtomlins/Google Drive/OCT Data/Test/";
    }
    //
    // Try to load the resampling table and reference spectrum from the csv files in the
    // specified folder
    //
    tempArray= (float*)malloc(sizeof(float)*MAX_CSV_FILE_LENGTH);
    resampleTableFileData = (float*)malloc(sizeof(float)*MAX_CSV_FILE_LENGTH);
    referenceSpectrumFileData = (float*)malloc(sizeof(float)*MAX_CSV_FILE_LENGTH);
    referenceAScanFileData = (float*)malloc(sizeof(float)*MAX_CSV_FILE_LENGTH);
    //
    // Try to read the reference spectrum and resampling table
    //
    if (readTwoColumnCSVFile(resamplingTablePath, &inputLen, MAX_CSV_FILE_LENGTH, tempArray, resampleTableFileData) != EXIT_SUCCESS)
        exit(EXIT_FAILURE);
    
    if (readTwoColumnCSVFile(referenceSpectrumPath, &inputLen,MAX_CSV_FILE_LENGTH, tempArray, referenceSpectrumFileData) != EXIT_SUCCESS)
        exit(EXIT_FAILURE);

    if (readTwoColumnCSVFile(referenceAScanPath, &outputLen,MAX_CSV_FILE_LENGTH, tempArray, referenceAScanFileData) != EXIT_SUCCESS)
        exit(EXIT_FAILURE);
    
    
    //
    // Now we have the spectrum length...
    // Allocate arrays and Populate variables
    //
    //outputLen = inputLen;
    numAScansPerBScan = numAScans * ascanAve * bscanAve;
    totalAScans = numBScansPerBatch * numAScansPerBScan;
    totalInputLen = inputLen * totalAScans;
    totalOutputLen = outputLen * totalAScans;
    resampTable = (float*)malloc(sizeof(float) * inputLen);
    refSpec = (float*)malloc(sizeof(float) * inputLen);
    refAScan = (float*)malloc(sizeof(float) * inputLen);
    inputSpectra = (short*)malloc(sizeof(short) * totalInputLen);
    buildLog = (char*)malloc(sizeof(char) * loglen);
    
    //
    // Allocate arrays for the output data
    //
    preProcessed = (float*)malloc(sizeof(float) * totalOutputLen*2);
    cmplx = (float*)malloc(sizeof(float) * totalOutputLen*2);
    //at* real = (float*)malloc(sizeof(float) * totalOutputLen);
    //at* imag = (float*)malloc(sizeof(float) * totalOutputLen);
    linEnv = (float*)malloc(sizeof(float) * totalOutputLen);
    logEnv = (float*)malloc(sizeof(float) * totalOutputLen);
    sum = (float*)malloc(sizeof(float) * totalAScans);
    sam = (float*)malloc(sizeof(float) * totalAScans);
    ad = (float*)malloc(sizeof(float) * totalAScans);
    correlationMap = (float*)malloc(sizeof(float) * totalOutputLen);
    
    //
    // Copy the data that was read from the files into the arrays used for input to clRTOCT
    //
    for (i=0; i<inputLen; i++)
    {
        resampTable[i] = resampleTableFileData[i];
        refSpec[i] = referenceSpectrumFileData[i];
        refAScan[i] = referenceAScanFileData[i];
        
    }
    
    
    //
    // Open binary files for reading and writing
    inputSpecFile = fopen(spectraPath,"rb");
    //
    // Output files
    //
  //  preProcessedFile = fopen("preProcessed.bin","wb");
    //LE* realFile = fopen("realPart.bin","wb");
    //LE* imagFile = fopen("imagPart.bin","wb");
  //  cmplxFile = fopen("cmplx.bin","wb");
  //  linEnvFile = fopen("linEnv.bin","wb");
  //  logEnvFile = fopen("logEnv.bin","wb");
 
    
    
    

//	test();
    allProcessElapsedMicroSeconds = 0;
    allProcessElapsedMilliSeconds = 0.0;
    ftProcessElapsedMicroSeconds = 0;
    ftProcessElapsedMilliSeconds = 0.0;
    postProcessElapsedMicroSeconds = 0;
    postProcessElapsedMilliSeconds = 0.0;
    preProcessElapsedMicroSeconds = 0;
    preProcessElapsedMilliSeconds = 0.0;
    elapsedSeconds = 0.0;
    elapsedMicroSeconds = 0;
    
	printf("Initialising GPU...");

	res = clOCTInit(
		clDeviceIndex,
		inputLen,
		outputLen,
		numBScansPerBatch,
		numAScans,
		ascanAve,
		bscanAve,
		resampTable,
		refSpec,
		refAScan, 
		kernelPath,
		buildLog, 
		&loglen,
		3*numAScans*ascanAve,
        64,
        2
		);

	if (res != CL_SUCCESS)
	{
		printf("\nError %d, Build log output:\n",res);
		printf(buildLog);
		printf("\n");
	}
	else
	{
		//t1=clock();
        gettimeofday(&timeBegin, NULL);
        //
        // Loop through each of the input spectra and process bscans in the specified batch size
        //
        printf("Processing %i B-Scans in batches of %i, each batch comprising %i A-Scans...\n", totalBScans, numBScansPerBatch, totalAScans);
        for (i=0;i<numBScanProcessingIteratations; i++)
        {
            //
            // Read the next set of BScans from file
            //
            printf(".");
            fseek(inputSpecFile, sizeof(short) * totalInputLen * i, SEEK_SET);
            fread(inputSpectra,sizeof(short), totalInputLen, inputSpecFile);
            // Need to add preprocessing step
//            printf("\nPre-Processing...");
            printf(".");
            
            
            //
            // Measure the time taken for the pre-processing step
            //
            res=Wait();
            gettimeofday(&timePreProcessBegin, NULL);
            
            
            res = PreProcess(inputSpectra,WINDOW_TYPE_BLACKMAN);	// This step includes copying the raw spectra from host to device
            
//            res=CopyPreProcessedSpectraToHost(preProcessed);
//            printf("preProcessed[0]=%i, preProcessed[100]=%i",(int)preProcessed[0],(int)preProcessed[1]);
//            printf("Done.  \nPre-processing returned %d.", res);
//            printf("\nComputing Transform...");
            printf(".");
            res = InverseTransform();
//            printf("Done.\n");
            //printf("Transform returned %d.\n", res);
//            printf("Post-processing...");
            printf(".");
            
            
            res = PostProcess();
            
            res=Wait();
            gettimeofday(&timePreProcessEnd, NULL);
            
            preProcessElapsedMicroSeconds = (long)(timePreProcessEnd.tv_sec - timePreProcessBegin.tv_sec)*1000000 + (long)(timePreProcessEnd.tv_usec-timePreProcessBegin.tv_usec);
            preProcessElapsedMilliSeconds += (double)preProcessElapsedMicroSeconds / 1000.0f;
            
            //res = Wait();
            
            //res = CorrelationMap(2, 2, 2, 0, 1, 0);
		
            if (res == CL_SUCCESS)
            {
                
//                printf("Done.\n");
//                printf("Copying data back to host...");
            
//            CopyLogEnvelopeToHost(
//                                  logEnv
//                                  );
            
                res = CopyAllResultsToHost(
                                       preProcessed,
                                       cmplx,
                                       //imag,
                                       linEnv,
                                       logEnv,
                                       sum,
                                       sam,
                                       ad
                                       );
                

                
                res = CopyCorrelationMapToHost(correlationMap);
                
                res=Wait(); // Wait for the GPU to complete
                //
                // Write the bscans to bmp file
                //
                if (saveBmp == true)
                {
                    for (j=0; j<numBScansPerBatch;j++)
                    {
                        sprintf(bmpPath, "%sbscan%0.4i.bmp\0", rootPath, i*numBScansPerBatch + j);
                        saveBitmap(bmpPath,&logEnv[j*numAScansPerBScan*outputLen], (unsigned short)numAScansPerBScan, outputLen, outputLen/2, bmpMinVal, bmpMaxVal);
                        
                        sprintf(bmpPath, "%scorrelationMap%0.4i.bmp\0", rootPath, i*numBScansPerBatch + j);
                        saveBitmap(bmpPath,&correlationMap[j*numAScansPerBScan*outputLen], (unsigned short)numAScansPerBScan, outputLen, outputLen/2, -1, 1);
                    }
                }
                
            
                printf("%f%% complete\n",(float)i/(float)numBScanProcessingIteratations*100.0f);
            
//                printf("Done.  Returned %d.\n", res);
            }
            else
            {
                printf("\nError %d during post-processing.\n",res);
            }
        }
//        printf("Done.");
        gettimeofday(&timeEnd, NULL);
        elapsedMicroSeconds = (long)(timeEnd.tv_sec-timeBegin.tv_sec)*1000000 + (long)(timeEnd.tv_usec - timeBegin.tv_usec);
        elapsedSeconds = elapsedMicroSeconds / 1000000.0f;
        
        //t2 = clock();
        //diff = (((float)t2 - (float)t1) / CLOCKS_PER_SEC ) ;
        printf(" Done.\n");
        printf("Executed in %f s \n",elapsedSeconds);
        printf("Timed processing took %f ms per batch.\n", preProcessElapsedMilliSeconds/(double)numBScanProcessingIteratations);
        
        
        //
        // Compute correlation on entire volume
        //
     //   res = CorrelationMapOnBScanVolume(logEnv, totalBScans, 2, 2, 1, 0, 0, 1, <#float *correlationMapVolume#>);
    }
	

	//
	// Save the output to file
	//
	if (res == CL_SUCCESS)
	{
//		fwrite(preProcessed,sizeof(float),totalAScans*outputLen*2,preProcessedFile);//
//		fwrite(cmplx,sizeof(float),totalAScans*outputLen*2,cmplxFile);
		//fwrite(real,sizeof(float),totalAScans*outputLen,realFile);
		//fwrite(imag,sizeof(float),totalAScans*outputLen,imagFile);
//		fwrite(linEnv,sizeof(float),totalAScans*outputLen,linEnvFile);
//		fwrite(logEnv,sizeof(float),totalAScans*outputLen,logEnvFile);
	}

    fclose(inputSpecFile);
//	fclose(preProcessedFile);
//	fclose(cmplxFile);
	//fclose(realFile);
	//fclose(imagFile);
//	fclose(linEnvFile);
//	fclose(logEnvFile);

	clOCTDispose();

	free(inputSpectra);
	free(resampTable);
	free(refSpec);
	free(refAScan);

	free(preProcessed);
	free(cmplx);
	//free(real);
	//free(imag);
	free(linEnv);
	free(logEnv);
    free(correlationMap);
	free(sum);
	free(sam);
    free(referenceSpectrumFileData);
    free(referenceAScanFileData);
    free(resampleTableFileData);
    free(tempArray);
    free(ad);
    free(buildLog);
}


int readTwoColumnCSVFile(char* sourceFile, unsigned int* len, unsigned int maxRows, float* col1Output, float* col2Output)
{
    FILE* fp;
     ssize_t read;
    size_t lineLen=255;
    char* line = (char*)malloc(sizeof(char)*lineLen);
   // float* output = (float*)malloc(sizeof(float)*2*maxRows);
    float col1;
    float col2;

     *len = 0;
   
    fp = fopen(sourceFile, "r");
    if (fp == NULL)
        exit(EXIT_FAILURE);
    
    
    //r = fscanf(fp, "%f,%f\r\n", &col1, &col2);
    
    while ( (*len <= maxRows) && ((read = getline(&line, &lineLen, fp)) != -1))
    {
        sscanf( line, "%f,%f", &col1, &col2 );
        col1Output[*len]=col1;
        col2Output[*len]=col2;
      //  printf("Retrieved line %i of length %zu : %s\n", *len,read, line);
        //printf("%s", line);
        (*len)++;
    }
    
    fclose(fp);
    if (line)
        free(line);
    
    return EXIT_SUCCESS;

}

int saveBitmap(char* path, float* sourceArray, unsigned short width, unsigned short sourceHeight, unsigned outputHeight, float minVal, float maxVal){
    //
    // Save BScan to RGB 24bpp device independent bitmap
    //
    float bitmapFileSize;
    float actualWidthInBytes;
    float requiredWidthInBytes;
    unsigned int pixelArraySize;
    unsigned short bpp = 24;
    unsigned short bytesPerPixel = bpp/8;
    unsigned int headerLen = 14;    // 14 byte bitmap header
    unsigned int dibHeaderLen = 12; // BITMAPCOREHEADER
    unsigned int offset = headerLen+dibHeaderLen;
    unsigned short planes = 1;
    unsigned short w;
    unsigned short h;
   // unsigned short paddedWidth = (unsigned short)ceil(((double)width * (double)bytesPerPixel)/4.0);
    unsigned short paddingWidth; // = paddedWidth - width;
    
    float srcPixel;
    float bmpPixel;
    float span = maxVal - minVal;
    unsigned short bytePixel;
	FILE* fbmp;
    
    actualWidthInBytes = (float)(bytesPerPixel * width);
    requiredWidthInBytes = (float)ceil((double)(actualWidthInBytes/4.0f)) * 4.0f;
    
    
    pixelArraySize = (unsigned int)(requiredWidthInBytes)*(unsigned int)outputHeight;         // Pixel array size in bytes
    bitmapFileSize = headerLen + dibHeaderLen + pixelArraySize;
    paddingWidth = (unsigned short)(requiredWidthInBytes - actualWidthInBytes);
    //
    // Write the source array to a device independent bitmap file
    //
    fbmp = fopen(path, "wb");
    //
    // Write a 14 byte header
    //
    //value="BM";
    fwrite("BM", 1, 2, fbmp);
    fwrite(&bitmapFileSize, 4, 1, fbmp);
    fwrite("0", 1, 4, fbmp);   // Write 4 zero bytes to the reserved section
    fwrite(&offset, 4, 1, fbmp); // Write the offset to the beginning of the image data
    //
    // Write DIB header
    //
    fwrite(&dibHeaderLen, 4, 1, fbmp);
    fwrite(&width, 2, 1, fbmp);
    fwrite(&outputHeight, 2, 1, fbmp);
    fwrite(&planes, 2, 1, fbmp);
    fwrite(&bpp, 2, 1, fbmp);
    //
    // Write the pixel array
    //
    for (h=0;h<outputHeight; h++)
    {
        for (w=0;w<width; w++)
        {
            srcPixel = sourceArray[w*sourceHeight+outputHeight-h];
            bmpPixel = (srcPixel - minVal)/span * 255;
            if (bmpPixel < 0.0f)
                bmpPixel = 0.0f;
            if (bmpPixel > 255.0f)
                bmpPixel = 255.0f;
            
            bytePixel = (unsigned short)bmpPixel;
            
            fwrite(&bytePixel,1,1,fbmp);  // Red
            fwrite(&bytePixel,1,1,fbmp);     // Green
            fwrite(&bytePixel,1,1,fbmp);     // Blue
            
        }
        
        for (w=0;w<paddingWidth;w++)    // Pad to a 4 byte boundary
        {
            fwrite("0",1,1,fbmp);
        }
    }
    
    fclose(fbmp);
    
    return EXIT_SUCCESS;
    
}
