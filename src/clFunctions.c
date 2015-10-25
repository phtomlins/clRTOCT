
#include "clFunctions.h"

//
// Initialise the open cl platform and device
//
int clInit(
           cl_context* clContext,
           cl_command_queue* clCommandQueue,
           cl_uint deviceIndex,
           char*** deviceNameList,
           cl_uint* numDevicesInList
           )
{
	cl_uint numPlatforms;	// Number of platforms
	cl_uint platformIDSize = sizeof(cl_platform_id);
	
	cl_uint j;
	cl_int err;

	_selectedDeviceIndex = deviceIndex;	// OpenCL device index in the device list
	//_selectedPlatformIndex = 0;
	_windowType = 0;					// Default to no window - but can be set during call to pre-processing kernel
	//
	// Get the number of platforms
	//
	err = clGetPlatformIDs(NULL, NULL, &numPlatforms);
	if (err != CL_SUCCESS) return err;
	//
	// Create an array to hold the platforms list
	//
	_platformIDs = (cl_platform_id*)malloc(sizeof(cl_platform_id) * numPlatforms); 
	err = clGetPlatformIDs(numPlatforms,_platformIDs,NULL);
	if (err != CL_SUCCESS) return err;
	///
	// Assume that we can use the first platform in the list.  Most systems will have only one platform installed, i.e. NVIDIA
	// although it's possible that the Intel SDK could also be installed.
	//
	if (numPlatforms < 1)
		return -1;
	//
	// Get a list of the platform IDs associated with the first platform
	//
	err = clGetDeviceIDs(_platformIDs[SELECTED_PLATFORM], CL_DEVICE_TYPE_ALL, 0, NULL, &_numDevices);	// Get the number of devices
	if (err < 0)
		return err;
	_deviceIDs = (cl_device_id*)malloc(_numDevices*sizeof(cl_device_id));	// Allocate an array for the devices
	_deviceNameList = (char**)malloc(_numDevices * sizeof(char*));		// Allocate a list of pointers to device names
	err = clGetDeviceIDs(_platformIDs[SELECTED_PLATFORM], CL_DEVICE_TYPE_ALL, _numDevices, _deviceIDs, NULL);	// Populate the array with the device IDs
	if (err < 0)
		return err;
	//
	// Create a list of the device names - mainly for debugging at the moment
	//
	for (j=0; j<_numDevices; j++)
	{

		//infoRetSize = 1024;	// A large return size variable
		_deviceNameList[j] = (char*)malloc(DEVICE_INFO_RETURN_SIZE);
		err = clGetDeviceInfo(_deviceIDs[j], CL_DEVICE_NAME, DEVICE_INFO_RETURN_SIZE, _deviceNameList[j], NULL);
		//printf()
		if (err < 0)
			return err;
        
        
	}
    //
    // Get the maximum local work group size
    // using         //CL_DEVICE_MAX_WORK_GROUP_SIZE
    //
    err = clGetDeviceInfo(_deviceIDs[_selectedDeviceIndex], CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(size_t), &_maxWorkGroupSize, NULL);
    //printf()
    if (err < 0)
        return err;
    printf("\nMaximum local work group size: %i \n", _maxWorkGroupSize);
	//
	// Also return a reference to the name list to the user
	//
	*deviceNameList = _deviceNameList;
	*numDevicesInList = _numDevices;
	//
	// Create a context for GPU execution - use the first platform for now
	//
	*clContext = clCreateContext(NULL,1,&_deviceIDs[_selectedDeviceIndex],NULL,NULL,&err);
	if (err < 0)
		return err;
	//
	// Create a command queue on which to execute the OCT processing pipeline
	//
	*clCommandQueue = clCreateCommandQueue(*clContext,_deviceIDs[_selectedDeviceIndex],0, &err);
	if (err < 0)
		return err;

	return 0;

}

//
// Build a kernel from a source file
//
int clBuild(char* sourceFile, char* buildLog, size_t* buildLogLength)
{
	//
	// Test code to check the expected output from openCL library functions
	// Use this to check correct outputs in C# wrapper functions
	//
//	cl_uint i;
	size_t logLen = 0;
	cl_uint j;
	size_t infoRetSize;
	cl_int err;
	//
	// Compile the kernels
	//
	if (buildLog == NULL)
	{
		logLen = 1000000;
		buildLogLength = &logLen;
		buildLog = (char*)malloc(sizeof(char) * logLen);
	}
	//printf("Comiling kernels...");

	err = clOCTCompileKernels(sourceFile, buildLog, buildLogLength);
	//
	// Clear up internally allocated memory
	//
	if (logLen > 0)
		free(buildLog);

	

	return err;
}
//

int clOCTCompileKernels(char* sourceFile, char* build_log, size_t* buildLogLength)
{
	//
	// Compile the Kernels used by the OCT processing pipeline
	//
//	char* sourceFile = "C:\\Users\\oct\\Google Drive\\OCT Control Software - octView\\clOCTKernels\\current\\octKernels.cl\0";//
//	char* sourceFile = "/Users/phtomlins/Google Drive/OCT Control Software - octView/clOCTKernels/20140520/octKernels.cl\0";
	//char* build_log;
	size_t log_size;
	unsigned int i;

	unsigned int sourceLength;
	cl_int err = 0;


	char* kernelSource = clOCTLoadKernelSourceFromFile(sourceFile, &sourceLength);
			
	
	
	//		return err;		//DEBUGGING!!

	if (sourceLength < 0)
		return sourceLength;
	//
	// Use the openCL library to create a program from the source code
	//
	_clOCTProgram = clCreateProgramWithSource(_context, 1, (const char**)&kernelSource,NULL, &err);
	if (err != CL_SUCCESS)
		return err;

	//
	// Build the kernels that reside within the source
	//
	//err = clBuildProgram(_clOCTProgram, _numDevices, _deviceIDs, NULL, NULL, NULL); 
 //   const char options[] = "-Werror -cl-std=CL1.1";
//    error = clBuildProgram(program, 1, &device, options, NULL, NULL);
	err = clBuildProgram(_clOCTProgram, NUM_DEVICES_TO_BUILD_FOR, &_deviceIDs[_selectedDeviceIndex], _compilerOptions, NULL, NULL);
	if (err != CL_SUCCESS)
	{
		//
		// If there was an error, get the build log
		//
		clGetProgramBuildInfo(_clOCTProgram, _deviceIDs[_selectedDeviceIndex], CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size);	// Get the size of the build log
		//build_log = (char*)malloc(sizeof(char)*log_size + 1);		// Allocate memory for the log
		if (log_size > *buildLogLength)
			build_log = "Insufficient space to store build log.\0";
		else
			clGetProgramBuildInfo(_clOCTProgram, _deviceIDs[_selectedDeviceIndex], CL_PROGRAM_BUILD_LOG, log_size, build_log, NULL);
		//
		// Terminate with null character
		//
		*buildLogLength = log_size;
		build_log[log_size-1] = '\0';
		build_log[log_size] = '\0';
		//
		return err;
	}
	else
	{
		*buildLogLength = 0;	
		build_log[0] = '\0';
		build_log[1] = '\0';
	}
	//
	// Create specific kernels
	//
	_preProcessingKernel = clCreateKernel(_clOCTProgram, "octPreProcessingKernel", &err);
	if (err != CL_SUCCESS)
		return err;
	_postProcessingKernel = clCreateKernel(_clOCTProgram, "octPostProcessingKernel", &err);
	if (err != CL_SUCCESS)
		return err;
    _octCorrelationKernel = clCreateKernel(_clOCTProgram, "octCorrelationKernel", &err);
    if (err != CL_SUCCESS)
        return err;
//	_imageKernel = clCreateKernel(_clOCTProgram, "octImageKernel", &err);
//	if (err != CL_SUCCESS)
//		return err;


	//	free(kernelSource);

	return 0;
}



char* clOCTLoadKernelSourceFromFile(char* sourceFile, unsigned int* len)
{
	//
	// Load a source cl file into a string.  Len is the length of the returned source file in characters, or negative if an error occurred
	// Return the loaded file as a null terminated string, or NULL if an error ocurred
	//
	int c;
	unsigned int i;
	unsigned int index = 0;
	char* tmpFileContent = (char*)malloc(sizeof(char) * MAX_PROGRAM_FILE_LENGTH);	// Specify maximum file length in character
	char* fileContent;
	FILE *file;
	file = fopen(sourceFile, "r");

	if (file) 
	{
		while (((c = getc(file)) != EOF) && (index < MAX_PROGRAM_FILE_LENGTH))	// c is integer because EOF can be signed
		{
			tmpFileContent[index] = (char)c;
			index++;
		}
		fclose(file);


		if (index == MAX_PROGRAM_FILE_LENGTH)
		{
			*len = -2;
			return NULL;
		}

		// Copy to a new array of the correct size

		fileContent = (char*)malloc(sizeof(char)*index);
		for (i = 0; i< index; i++)
		{
			fileContent[i] = tmpFileContent[i];
		}
		fileContent[index] = '\0';	// Terminate string with null character

		free(tmpFileContent);
	}
	else
	{
		*len = -1;
		return NULL;
	}
	*len = index-1;
	return fileContent;
}


int clAlloc(
		cl_context context,
		cl_uint inputSpectraLength,		// Specify the actual length of input spectra
		cl_uint outputAScanLength,		// Specify the length of output (outputLength >= inputLength).  If inputlength < outputlength then the input spectra will be zero padded
		cl_uint numBScans,               // Number of B-Scans processed in one go
		cl_uint numAScansPerBScan,
		cl_uint ascanAveragingFactor,
		cl_uint bscanAveragingFactor,
		float* hostResamplingTable,		// Resampling table
		float* hostInterpolationMatrix,		// Pre-computed interpolation matrix
		float* hostReferenceSpectrum,	// Reference spectrum to be subtracted (or maybe divided!!) from each spectrum
		float* hostReferenceAScan,		// A reference background A-Scan to be subtracted.. maybe log domain!!
		unsigned int imageFormatStride
		)
{
	//
	// Try to allocate space on the GPU in global memory, this avoids repeating the allocation step for each DFT call
	// Copy constant arrays, i.e. resampling table and spectral reference
	//
	cl_int err = 0;
	cl_uint i;
    //
    
	//
	// Copy values to global variables
	//
	_inputSpectrumLength = inputSpectraLength;
	_outputAScanLength = outputAScanLength;
	_outputImageHeight = _outputAScanLength/2;
	_numBScans = numBScans;
	_numAScansPerBScan = numAScansPerBScan;
	_ascanAveragingFactor=ascanAveragingFactor;
	_bscanAveragingFactor=bscanAveragingFactor;
    _imageFormatStride = imageFormatStride;
    //
    // Derived quantities used for allocating memory on the GPU
    //
	_totalAScansPerBScan =  (size_t)numAScansPerBScan * (size_t)ascanAveragingFactor;
    _totalBScans = (size_t)_numBScans * (size_t)_bscanAveragingFactor;
    _totalAScans = _totalBScans * _totalAScansPerBScan;
    _totalInputSpectraLength = _totalAScans * (size_t)_inputSpectrumLength;
    _totalOutputAScanLength = _totalAScans * (size_t)_outputAScanLength;
    _totalPreProcessedSpectraLength = _totalOutputAScanLength * (size_t)2;
    _bitmapBScanSize = (size_t)_outputImageHeight * (size_t)_imageFormatStride;
    _bitmapBScanVolumeSize = (size_t)_numBScans * (size_t)_bscanAveragingFactor * _bitmapBScanSize;
    
    //
    _correlationUsesLogBScans = TRUE;       // default is for correlation mapping to use logarithmic bscans
	//
	// Create a buffer on the device for the input spectra to be copied to
	//
	deviceSpectra = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_ALLOC_HOST_PTR, sizeof(short) * _totalInputSpectraLength, NULL, &err);
	if (err != CL_SUCCESS)
		return err;
	//
	// Copy the resampling table, reference spectrum and reference a-scan immediately (no need to call enqueuewritebuffer)
	// using the CL_MEM_COPY_HOST_PTR flag
	//
	deviceResamplingTable = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(float) * (size_t)_inputSpectrumLength, hostResamplingTable, &err);
	if (err != CL_SUCCESS)
		return err;
	deviceInterpolationMatrix = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(float) * (size_t)_inputSpectrumLength * (size_t)(INTERPOLATION_ORDER + 1) * (size_t)(INTERPOLATION_ORDER + 1), hostInterpolationMatrix, &err);
	if (err != CL_SUCCESS)
		return err;
	deviceReferenceSpectrum = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(float) * (size_t)_inputSpectrumLength, hostReferenceSpectrum, &err);
	if (err != CL_SUCCESS)
		return err;
	deviceReferenceAScan = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(float) * (size_t)_outputAScanLength, hostReferenceAScan, &err);
	if (err != CL_SUCCESS)
		return err;
	//
	// Allocate device memory for the outputs
	//
	// clFFT runs on an array of cl_mem objects, so we'll create the resampled spectra as such an array
	//
	devicePreProcessedSpectra = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(float) * _totalPreProcessedSpectraLength, NULL, &err);	// Complex interleaved pre-processed spectra
	if (err != CL_SUCCESS)
		return err;
	//
	// Likewise create an array of output buffers to hold the transform
	//
	deviceFourierTransform = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(float) * _totalPreProcessedSpectraLength, NULL, &err);	// Complex interleaved
	if (err != CL_SUCCESS)
		return err;
	
	//
	// Create buffers for each of the outputs
	//
	deviceEnvBScan = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(float)*_totalOutputAScanLength, NULL, &err);	// Envelope of the fft
	if (err != CL_SUCCESS)
		return err;
	deviceLogEnvBScan = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(float)*_totalOutputAScanLength, NULL, &err);	// Logarithmic Envelope of the fft
	if (err != CL_SUCCESS)
		return err;
    deviceCorrelationMap = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(float)*_totalOutputAScanLength, NULL, &err);	// Correlation map
    if (err != CL_SUCCESS)
        return err;
	deviceSum = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(float) * _totalAScans, NULL, &err);	// Sum along each a-scan
	if (err != CL_SUCCESS)
		return err;
	deviceSAM = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(float) * _totalAScans, NULL, &err);	// attenuation measured along each a-scan
	if (err != CL_SUCCESS)
		return err;
	deviceAttenuationDepth = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(float) * _totalAScans, NULL, &err);	// attenuation depth along each a-scan
	if (err != CL_SUCCESS)
		return err;
	deviceBScanBmp = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(unsigned char) * _bitmapBScanVolumeSize, NULL, &err);	// BScan bitmap images
	if (err != CL_SUCCESS)
		return err;
		
	return err;
}

//
// Function releases the cl_mem object buffers allocated on the device
//
int clRelease()
{
	cl_uint i;

	// Release inputs
	clReleaseMemObject(deviceSpectra);
	clReleaseMemObject(deviceResamplingTable);
	clReleaseMemObject(deviceInterpolationMatrix);
	clReleaseMemObject(deviceReferenceSpectrum);
	clReleaseMemObject(deviceReferenceAScan);
	//
	// Release outputs
	//
	//for (i=0; i<_totalAScans; i++)
	//{
		clReleaseMemObject(devicePreProcessedSpectra);
		clReleaseMemObject(deviceFourierTransform);
//		clReleaseMemObject(devicePreProcessedSpectra[1]);
//		clReleaseMemObject(deviceFourierTransform[1]);
	//}
	//free(devicePreProcessedSpectra);
	//free(deviceFourierTransform);

	//clReleaseMemObject(deviceRealBScan);
	//clReleaseMemObject(deviceImagBScan);
	clReleaseMemObject(deviceEnvBScan);
	clReleaseMemObject(deviceLogEnvBScan);
    clReleaseMemObject(deviceCorrelationMap);
	clReleaseMemObject(deviceSum);
	clReleaseMemObject(deviceSAM);
	clReleaseMemObject(deviceAttenuationDepth);
	clReleaseMemObject(deviceBScanBmp);
	//
	// Release other resources
	//
	for (i=0; i<_numDevices; i++)
	{
		free(_deviceNameList[i]);
	}
//	clReleaseKernel(dftKernel);
////	clReleaseKernel(testKernel);
	clReleaseKernel(_preProcessingKernel);
	clReleaseKernel(_postProcessingKernel);
	//clReleaseKernel(_imageKernel);
    clReleaseKernel(_octCorrelationKernel);

	clReleaseProgram(_clOCTProgram);


	clReleaseCommandQueue(_commandQueue);
	clReleaseContext(_context);
//	
	free(_deviceNameList);
	free(_deviceIDs);
	free(_platformIDs);


	return 0;
}
/*
int SetImageKernelParameters(
                             cl_mem floatArray,
                             unsigned int inputHeight,
                             unsigned int inputWidth,
                             
                             
                             )
{
    //
    // Set some of the parameters for the image conversion kernel
    //
 
     __kernel void octImageKernel(
     __global float* floatImage,       // 0. Input array of floats
     const unsigned int inputHeight,   // 1. Height of image as in input array
     const unsigned int inputWidth,    // 2. Width of image in input array
     const unsigned int imageFormatStride,    // 3. Stride for bitmap image
     const unsigned int startIndex,    // 4. Index of first image pixel
     const float minVal,               // 5. Min value
     const float maxVal,               // 6. Max value
     const unsigned int outputTop,     // 7. Top pixel of ouptput image
     const unsigned int outputHeight,  // 8. Height of the output image
     __global unsigned char* byteImage // 9. Raw output bitmap
     )
 
    cl_int err = 0;
    err = clSetKernelArg(_imageKernel, 0, sizeof(cl_mem), &deviceLogEnvBScan);	// Takes the log b-scan and converts it to and image array
    if (err != CL_SUCCESS) return err;
    err = clSetKernelArg(_imageKernel, 1, sizeof(unsigned int), &_totalAScansPerBScan);
    if (err != CL_SUCCESS) return err;
    err = clSetKernelArg(_imageKernel, 2, sizeof(unsigned int), &_outputImageHeight);
    if (err != CL_SUCCESS) return err;
    err = clSetKernelArg(_imageKernel, 3, sizeof(unsigned int), &_imageFormatStride);
    if (err != CL_SUCCESS) return err;
    err = clSetKernelArg(_imageKernel, 7, sizeof(cl_mem), &deviceBScanBmp);
    if (err != CL_SUCCESS) return err;
    
    return err;
    
    
}

*/

int SetPreProcessingKernelParameters()
{
	/*
	__kernel void octPreProcessingKernel(
						0 __global const short* spectra,					// Input arrays - raw spectra					
						1 __global const float* resamplingTable,			// resampling table
						2 __global const float* interpolationMatrix,
						3 __global const float* referenceSpectrum,		// reference spectrum for deconvolution
						4 const cl_uint inputSpectraLength,				
						5 const cl_uint outputAScanLength,				
						6 const cl_uint numBScans,						
						7 const cl_uint numAScansPerBScan,				
						8 const cl_uint ascanAveragingFactor,
						9 const cl_uint bscanAveragingFactor,
						10 const cl_uint windowType,
						11 __global float* preProcessedCmplxSpectra				// Output is the pre-processed spectra, ready for FFT
						)
	*/
	cl_int err = 0;
	err = clSetKernelArg(_preProcessingKernel, 0, sizeof(cl_mem), &deviceSpectra);
	if (err != CL_SUCCESS) return err;
	err = clSetKernelArg(_preProcessingKernel, 1, sizeof(cl_mem), &deviceResamplingTable);
	if (err != CL_SUCCESS) return err;
	err = clSetKernelArg(_preProcessingKernel, 2, sizeof(cl_mem), &deviceInterpolationMatrix);
	if (err != CL_SUCCESS) return err;
	err = clSetKernelArg(_preProcessingKernel, 3, sizeof(cl_mem), &deviceReferenceSpectrum);
	if (err != CL_SUCCESS) return err;
	err = clSetKernelArg(_preProcessingKernel, 4, sizeof(cl_uint), &_inputSpectrumLength);
	if (err != CL_SUCCESS) return err;
	err = clSetKernelArg(_preProcessingKernel, 5, sizeof(cl_uint), &_outputAScanLength);
	if (err != CL_SUCCESS) return err;
	err = clSetKernelArg(_preProcessingKernel, 6, sizeof(cl_uint), &_numBScans);
	if (err != CL_SUCCESS) return err;
	err = clSetKernelArg(_preProcessingKernel, 7, sizeof(cl_uint), &_numAScansPerBScan);
	if (err != CL_SUCCESS) return err;
	err = clSetKernelArg(_preProcessingKernel, 8, sizeof(cl_uint), &_ascanAveragingFactor);
	if (err != CL_SUCCESS) return err;
	err = clSetKernelArg(_preProcessingKernel, 9, sizeof(cl_uint), &_bscanAveragingFactor);
	if (err != CL_SUCCESS) return err;
	err = clSetKernelArg(_preProcessingKernel, 10, sizeof(cl_uint), &_windowType);
	if (err != CL_SUCCESS) return err;
	err = clSetKernelArg(_preProcessingKernel, 11, sizeof(cl_mem), &devicePreProcessedSpectra);	// real part
	if (err != CL_SUCCESS) return err;
	//err = clSetKernelArg(_preProcessingKernel, 11, sizeof(cl_mem), &devicePreProcessedSpectra[1]);	// imaginary part
	//if (err != CL_SUCCESS) return err;

	return err;


}

int SetPostProcessingKernelParameters()
{
	/*
__kernel void octPostProcessingKernel(
     __global const float* FourierTransformCmplx,		// 0. Input, the Fourier transform of the pre-processed spectra
     __global const float* referenceAScan,				// 1. Reference a-scan subtracted from envelope of FT
     const unsigned int outputAScanLength,				// 2. Length of the 1D Fourier transform
     const unsigned int numBScans,						// 3. Number of BScans loaded into memory
     const unsigned int numAScansPerBScan,               // 4. Number of A-Scans per BScan
     const unsigned int ascanAveragingFactor,            // 5. AScan averaging factor
     const unsigned int bscanAveragingFactor,            // 6. BScan Averaging factor
     const unsigned int bscanBmpFormatStride,            // 7. Bitmap image row length for bscans
     const float minVal,                                 // 8. Minimum image value for bitmaps
     const float maxVal,                                 // 9. Maximum image value for bitmaps
     __global float* envBScan,							// 10. Envelope of the bscan
     __global float* logEnvBScan,						// 11. Log of envelope
     __global float* Sum,                                // 12. Integral along a-scan
     __global float* SAM,                                // 13. SAM image produced by measuring slope (not implemented yet)
     __global float* AttenuationDepth,                   // 14. Attenuation depth measurement (not yet implemented)
     __global char* bscanBmp                             // 15. BScan Bitmap Array
						)	*/
	cl_int err = 0;
	err = clSetKernelArg(_postProcessingKernel, 0, sizeof(cl_mem), &deviceFourierTransform);
	if (err != CL_SUCCESS) return err;
	err = clSetKernelArg(_postProcessingKernel, 1, sizeof(cl_mem), &deviceReferenceAScan);
	if (err != CL_SUCCESS) return err;
	err = clSetKernelArg(_postProcessingKernel, 2, sizeof(cl_uint), &_outputAScanLength);
	if (err != CL_SUCCESS) return err;
	err = clSetKernelArg(_postProcessingKernel, 3, sizeof(cl_uint), &_numBScans);
	if (err != CL_SUCCESS) return err;
	err = clSetKernelArg(_postProcessingKernel, 4, sizeof(cl_uint), &_numAScansPerBScan);
	if (err != CL_SUCCESS) return err;
	err = clSetKernelArg(_postProcessingKernel, 5, sizeof(cl_uint), &_ascanAveragingFactor);
	if (err != CL_SUCCESS) return err;
	err = clSetKernelArg(_postProcessingKernel, 6, sizeof(cl_uint), &_bscanAveragingFactor);
	if (err != CL_SUCCESS) return err;
    err = clSetKernelArg(_postProcessingKernel, 7, sizeof(cl_uint), &_imageFormatStride);	// Envelope of BScan
    if (err != CL_SUCCESS) return err;
	err = clSetKernelArg(_postProcessingKernel, 10, sizeof(cl_mem), &deviceEnvBScan);	// Envelope of BScan
	if (err != CL_SUCCESS) return err;
	err = clSetKernelArg(_postProcessingKernel, 11, sizeof(cl_mem), &deviceLogEnvBScan);	// Log envelope - usually used for B-Scan images
	if (err != CL_SUCCESS) return err;
	err = clSetKernelArg(_postProcessingKernel, 12, sizeof(cl_mem), &deviceSum);	// Log envelope - usually used for B-Scan images
	if (err != CL_SUCCESS) return err;
	err = clSetKernelArg(_postProcessingKernel, 13, sizeof(cl_mem), &deviceSAM);	// Log envelope - usually used for B-Scan images
	if (err != CL_SUCCESS) return err;
	err = clSetKernelArg(_postProcessingKernel, 14, sizeof(cl_mem), &deviceAttenuationDepth);	// Log envelope - usually used for B-Scan images
	if (err != CL_SUCCESS) return err;
    err = clSetKernelArg(_postProcessingKernel, 15, sizeof(cl_mem), &deviceBScanBmp);	// Log envelope - usually used for B-Scan images
    if (err != CL_SUCCESS) return err;

	return err;


}


int SetCorrelationKernelParameters()
{
    /*
     __kernel void octCorrelationKernel(
   0  __global float* bscans,					// Input bscans
   1  const unsigned int AScanLength,
   2  const unsigned int numBScans,
   3  const unsigned int numAScansPerBScan,
   4  const unsigned int ascanAveragingFactor,
   5  const unsigned int bscanAveragingFactor,
   6  const unsigned int corrSizeX,  // Correlation window size across b-scan
   7  const unsigned int corrSizeY,  // Correlation window size across multiple b-scans (<=numBScans)
   8  const unsigned int corrSizeZ,  // Correlation window size axially
   9  const unsigned int offsetX,      // Offset between the centre of regions being correlated
   10  const unsigned int offsetY,      // Should allow elastographic analysis
   11  const unsigned int offsetZ,      //
   12  __global float* correlationMap				// Correlation coefficients
     )
     )	*/
    cl_int err = 0;
    if (_correlationUsesLogBScans)
    {
        err = clSetKernelArg(_octCorrelationKernel, 0, sizeof(cl_mem), &deviceLogEnvBScan);
        if (err != CL_SUCCESS) return err;
    }
    else
    {
        err = clSetKernelArg(_octCorrelationKernel, 0, sizeof(cl_mem), &deviceEnvBScan);
        if (err != CL_SUCCESS) return err;
    }
    err = clSetKernelArg(_octCorrelationKernel, 1, sizeof(unsigned int), &_outputAScanLength);
    if (err != CL_SUCCESS) return err;
    err = clSetKernelArg(_octCorrelationKernel, 2, sizeof(unsigned int), &_numBScans);
    if (err != CL_SUCCESS) return err;
    err = clSetKernelArg(_octCorrelationKernel, 3, sizeof(unsigned int), &_numAScansPerBScan);
    if (err != CL_SUCCESS) return err;
    err = clSetKernelArg(_octCorrelationKernel, 4, sizeof(unsigned int), &_ascanAveragingFactor);
    if (err != CL_SUCCESS) return err;
    err = clSetKernelArg(_octCorrelationKernel, 5, sizeof(unsigned int), &_bscanAveragingFactor);
    if (err != CL_SUCCESS) return err;
    err = clSetKernelArg(_octCorrelationKernel, 6, sizeof(unsigned int), &_corrSizeX);
    if (err != CL_SUCCESS) return err;
    err = clSetKernelArg(_octCorrelationKernel, 7, sizeof(unsigned int), &_corrSizeY);
    if (err != CL_SUCCESS) return err;
    err = clSetKernelArg(_octCorrelationKernel, 8, sizeof(unsigned int), &_corrSizeZ);
    if (err != CL_SUCCESS) return err;
    err = clSetKernelArg(_octCorrelationKernel, 9, sizeof(unsigned int), &_offsetX);
    if (err != CL_SUCCESS) return err;
    err = clSetKernelArg(_octCorrelationKernel, 10, sizeof(unsigned int), &_offsetY);
    if (err != CL_SUCCESS) return err;
    err = clSetKernelArg(_octCorrelationKernel, 11, sizeof(unsigned int), &_offsetZ);
    if (err != CL_SUCCESS) return err;
    err = clSetKernelArg(_octCorrelationKernel, 12, sizeof(cl_mem), &deviceCorrelationMap);
    if (err != CL_SUCCESS) return err;
    

    return err;
    
}


