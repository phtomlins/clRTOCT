//
#include "clRTOCT.h"
//
// File will contain the main library functions for using clFFT and OpenCL to process OCT data in real-time
//
//
// Exported Functions
//
//
// Initialise the clOCT library
//
DLL clOCTInit(
    cl_uint clDeviceIndex,			// Index in the device list of the OpenCL capable device to use
    cl_uint inputSpectraLength,
    cl_uint outputAScanLength,
    cl_uint numBScans,
    cl_uint numAScansPerBScan,
    cl_uint ascanAveragingFactor,
    cl_uint bscanAveragingFactor,
	float* hostResamplingTable,		// Resampling table
	float* hostReferenceSpectrum,	// Reference spectrum to be subtracted (or maybe divided!!) from each spectrum
	float* hostReferenceAScan,		// A reference background A-Scan to be subtracted.. maybe log domain!!
	char* kernelPath,
	char* clBuildLog, 
	size_t* clBuildLogLength,
	unsigned int imageFormatStride,
    size_t preProcessingkernelWorkgroupSize,
    size_t postProcessingkernelWorkgroupSize,
    char* compilerOptions
	)
{
	clfftStatus status;
	int clErr;
	size_t strideX;
	size_t dist;
	//size_t batchSize;
	char** deviceList;
	cl_uint numDevicesInList;
	//int i;

	//
	// Debugging output
	//
	/*
	printf("\n");
	printf("Resampe Table, Ref. Spectrum, Ref. AScan\n");
	for (i=0; i<10; i++)
	{
		printf("%f8,  %f8,  %f8\n",hostResamplingTable[i],hostReferenceSpectrum[i],hostReferenceAScan[i]);
	}
	*/
	//
	// Set global variable values
	//
	_kernelPath = kernelPath;//"C:\\Users\\oct\\Google Drive\\OCT Control Software - octView\\clOCTKernels\\current\\octProcessingKernels2.cl\0";
	_fftLength = (size_t)outputAScanLength;	// Global FFT length for clFFT to expect
	_numCommandQueues = 1;
    
    _preProcessKernelWorkGroupSize = preProcessingkernelWorkgroupSize;
    _postProcessKernelWorkGroupSize = postProcessingkernelWorkgroupSize;
    _compilerOptions = compilerOptions;

	strideX = 1; 
	dist = _fftLength;
	//batchSize = numBScans*numAScansPerBScan*ascanAveragingFactor*bscanAveragingFactor;  // total number of ascans passed to GPU
	//
	//
	// Initalise the GPU
	//
	clErr = clInit(&_context, &_commandQueue,clDeviceIndex, &deviceList, &numDevicesInList);
	if (clErr != CL_SUCCESS) return(clErr);

	printf("\nCL Selected Device: ");
	printf(deviceList[clDeviceIndex]);
	printf("\n");
	//
	// Pre-compute the interpolation coefficients
	//
	_interpolationMatrix = (float*)malloc(sizeof(float) * (size_t)inputSpectraLength * (INTERPOLATION_ORDER + 1) * (INTERPOLATION_ORDER + 1));
	PreComputeInterpolationCoefficients(hostResamplingTable, _interpolationMatrix, (int)inputSpectraLength);
	//
	//
	// Allocate memory on the GPU device
	//
	clErr = clAlloc(
		_context,
		inputSpectraLength,		// Specify the actual length of input spectra
		outputAScanLength,		// Specify the length of output (outputLength >= inputLength).  If inputlength < outputlength then the input spectra will be zero padded
		numBScans,
		numAScansPerBScan,
		ascanAveragingFactor,
		bscanAveragingFactor,
		hostResamplingTable,		// Resampling table
		_interpolationMatrix,
		hostReferenceSpectrum,	// Reference spectrum to be subtracted (or maybe divided!!) from each spectrum
		hostReferenceAScan,		// A reference background A-Scan to be subtracted.. maybe log domain!!
		imageFormatStride
		);
	if (clErr != CL_SUCCESS) return(clErr);
	//
	// Build the OpenCL kernels
	//
	clErr = clBuild(_kernelPath,clBuildLog,clBuildLogLength);
	if (clErr != CL_SUCCESS) return clErr;
	//
	// Set the pre- and post- processing kernel parameters
	//
	clErr= SetPreProcessingKernelParameters();
	if (clErr != CL_SUCCESS) return clErr;
	clErr= SetPostProcessingKernelParameters();
	if (clErr != CL_SUCCESS) return clErr;
	//clErr = SetImageKernelParameters();
	//if (clErr != CL_SUCCESS) return clErr;
    clErr = SetCorrelationKernelParameters();
    if (clErr != CL_SUCCESS) return clErr;
    //
    // Determine the preferred workgroup size for each kernel
    // using CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE
    /*
    clErr = clGetKernelWorkGroupInfo (_preProcessingKernel,
                                     _deviceIDs[_selectedDeviceIndex],
                                     CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE,
                                     sizeof(size_t),
                                     &_preProcessKernelWorkGroupSize,
                                     NULL);
    */
    //printf("\nPreferred preProcessing work group size: %i\n", _preProcessKernelWorkGroupSize);
	//
	//
	// Configure the setup data
	//
    // NOTE:  clfftInitSetupData is defined __inline in the clFFT source
    // This breaks gcc when compiled via xcode
    // Update clfft.h to declare
    // static __inline
    //	static __inline clfftStatus clfftInitSetupData( clfftSetupData* setupData )
    //
	status = clfftInitSetupData(&_clfftSetupData);
	if (status != CLFFT_SUCCESS)
		return((int)status);

	//
	// Setup clFFT
	//
	status = clfftSetup(&_clfftSetupData);
	if (status != CLFFT_SUCCESS)
		return((int)status);

	//
	// Create an FFT plan using the CL context that was initialised earlier
	//
	status = clfftCreateDefaultPlan(&_clfftPlanHandle, _context, CLFFT_1D, &_fftLength);
	if (status != CLFFT_SUCCESS)
		return((int)status);
	//
	// Set the plan parameters (see example at http://dournac.org/info/fft_gpu)
	//
	status = clfftSetResultLocation(_clfftPlanHandle, CLFFT_OUTOFPLACE);
	status = clfftSetLayout(_clfftPlanHandle, CLFFT_COMPLEX_INTERLEAVED , CLFFT_COMPLEX_INTERLEAVED );	// Currently complex to complex transform... 
	status = clfftSetPlanBatchSize(_clfftPlanHandle,_totalAScans);
	status = clfftSetPlanPrecision(_clfftPlanHandle, CLFFT_SINGLE);

	//status = clfftSetPlanInStride(_clfftPlanHandle, CLFFT_1D, &strideX);
	//status = clfftSetPlanOutStride(_clfftPlanHandle, CLFFT_1D, &strideX);
	//status = clfftSetPlanDistance(_clfftPlanHandle, _fftLength, _fftLength);
	//
	// Finalise the plan (apparently we can specify it a bit more... will look into this later)
	//
	status = clfftBakePlan(_clfftPlanHandle, _numCommandQueues, &_commandQueue,NULL,NULL);
	if (status != CLFFT_SUCCESS)
		return((int)status);
	//
	// Create a temporary buffer for the plan
	//
	status = clfftGetTmpBufSize(_clfftPlanHandle, &_tempBufferSize);
	if (status != CLFFT_SUCCESS) return status;

	if (_tempBufferSize > 0)
	{
		_deviceTemporaryBuffer = clCreateBuffer(_context, CL_MEM_READ_WRITE, _tempBufferSize, 0, &clErr);
		if (clErr != CL_SUCCESS) return clErr;
	}

	return 0;
}

//
// Function used to properly clear up memory objects
//
DLL clOCTDispose()
{
	clfftStatus status;
	//
	// Destroy the plan
	//
	status = clfftDestroyPlan(&_clfftPlanHandle);

	status = clfftTeardown();
	if (status != CLFFT_SUCCESS)
		return((int)status);
	//
	// Release the temporary buffer
	//
	clReleaseMemObject(_deviceTemporaryBuffer);
	//
	// Release the openCL mem objects
	//
	clRelease();
	//
	// Free the interpolation matrix
	//
	free(_interpolationMatrix);

	return 0;
}

//
// Pre compute the interpolation coefficients - this function is used internally
//
void PreComputeInterpolationCoefficients(float* resamplingTable, float* coefs, int ResamplingTableLength)
{
	//
	// Use a quadratic interpolation
	// Therefore, a curve is defined by three points, requiring the calculation of 9 coefficients.
	//
	int i,j;
	int interpOrder = INTERPOLATION_ORDER;		// 2nd order interpolation
	int numCoefs = (interpOrder + 1) * (interpOrder + 1);
	int offset;
	float x0;
	float x1;
	float x2;
	float x3;
	//
	for (i=0; i<ResamplingTableLength; i++)
	{
		x0 = resamplingTable[i];
		x1 = floor(x0);
		x2 = ceil(x0);
		x3 = x2 + 1;
		offset = i*numCoefs;

		if (
				(abs(x1-x2) < 0.5)	// If the first two points are at the same position, then don't interpolate
				||
				(i>=ResamplingTableLength-interpOrder)	// Also don't inerpolate the last interpOrder number of points
			)
		{
			coefs[offset] = 1.0f;
			for (j=1; j< numCoefs; j++)
				coefs[offset+j] = 0.0f;		// All of the other coefficients are zero
		}
		else	// Otherwise, compute the matrix
		{
			// First row
			coefs[offset] = x2*x3 / ( (x1-x2)*(x1-x3) );
			coefs[offset+1] = -x1*x3 / ( (x1-x2)*(x2-x3) );
			coefs[offset+2] = x1*x2 / ( (x1-x3)*(x2-x3) );
			// Second row
			coefs[offset+3] = -(x2+x3) / ( (x1-x2)*(x1-x3) );
			coefs[offset+4] = (x1+x3) / ( (x1-x2)*(x2-x3) );
			coefs[offset+5] = -(x1+x2) / ( (x1-x3)*(x2-x3) );
			// Third row
			coefs[offset+6] = 1.0f / ( (x1-x2)*(x1-x3) );
			coefs[offset+7] = -1.0f / ( (x1-x2)*(x2-x3) );
			coefs[offset+8] = 1.0f / ( (x1-x3)*(x2-x3) );

/*			if (
				_isnan(coefs[offset]) || !_finite(coefs[offset]) ||
				_isnan(coefs[offset+1]) || !_finite(coefs[offset+1]) ||
				_isnan(coefs[offset+2]) || !_finite(coefs[offset+2]) ||
				_isnan(coefs[offset+3]) || !_finite(coefs[offset+3]) ||
				_isnan(coefs[offset+4]) || !_finite(coefs[offset+4]) ||
				_isnan(coefs[offset+5]) || !_finite(coefs[offset+5]) ||
				_isnan(coefs[offset+6]) || !_finite(coefs[offset+6]) ||
				_isnan(coefs[offset+7]) || !_finite(coefs[offset+7]) ||
				_isnan(coefs[offset+8]) || !_finite(coefs[offset+8]) 
				)
			{
		printf("------------------------------------------------\r\n");
		printf("| %f\t%f\t%f |\r\n",coefs[offset],coefs[offset+1],coefs[offset+2]);
		printf("| %f\t%f\t%f |\r\n",coefs[offset+3],coefs[offset+4],coefs[offset+5]);
		printf("| %f\t%f\t%f |\r\n",coefs[offset+6],coefs[offset+7],coefs[offset+8]);
		printf("------------------------------------------------\r\n");
				
			}
*/
		}
		/*
		printf("------------------------------------------------\r\n");
		printf("| %f\t%f\t%f |\r\n",coefs[offset],coefs[offset+1],coefs[offset+2]);
		printf("| %f\t%f\t%f |\r\n",coefs[offset+3],coefs[offset+4],coefs[offset+5]);
		printf("| %f\t%f\t%f |\r\n",coefs[offset+6],coefs[offset+7],coefs[offset+8]);
		printf("------------------------------------------------\r\n");
		*/
	}
}

//
// Function to PreProcess the data - this function handles copying to the GPU device
// resampling, windowing, zero-padding and any other pre-processing that is required
// This is handled by a pre-process kernel
//
DLL PreProcess(short* hostSpectra, int windowType)
{
	//int i;
	cl_int err;
    size_t numWorkItemsPerGroup = _preProcessKernelWorkGroupSize;
//    size_t numWorkGroups;
	size_t totalWorkItems = (size_t)_totalAScans;
	size_t spectraSize = _totalInputSpectraLength * sizeof(short);
    
//    numWorkGroups = (size_t)ceil((double)totalWorkItems / (double)numWorkItemsPerGroup);
	//
	_windowType = windowType;	// Set the window type before calling the kernel...
	//
	// Update the windowtype in the kernel argument
	//
	err = clSetKernelArg(_preProcessingKernel, 10, sizeof(unsigned int), &_windowType);
	if (err != CL_SUCCESS) return err;
/*
	printf("\n");
	for (i=0; i<10; i++)
	{
		printf("%i\n",hostSpectra[i]);
	}
	*/
	//
	// Copy the raw spectra from the host to the device
	//
	err = clEnqueueWriteBuffer(_commandQueue,
								deviceSpectra,
								CL_FALSE,
								0,
								spectraSize,
								hostSpectra,
								0,
								NULL,
								NULL);
	if (err != CL_SUCCESS)
		return err;
    
    //err = clFinish(_commandQueue);
    //if (err != CL_SUCCESS)	return err;
	//
	// Enqueue the pre processing kernel on the device
	//
	err = clEnqueueNDRangeKernel(_commandQueue, _preProcessingKernel, 1, NULL, &totalWorkItems, &numWorkItemsPerGroup, 0, NULL, NULL);
	if (err != CL_SUCCESS)
		return err;

	//err = clFinish(_commandQueue);
	//if (err != CL_SUCCESS)	return err;

	return 0;
}

//
// Computes the Fourier transform on the pre-processed data
//
DLL InverseTransform()
{
	clfftStatus status;
	cl_int err;

//
// For debugging
//
	//err = clEnqueueWriteBuffer(_commandQueue,
	//							deviceSpectra,
	//							CL_TRUE,
	//							0,
	//							sizeof(short) * _totalAScans * _inputSpectrumLength,
	//							hostSpectra,
	//							0,
	//							NULL,
	//							NULL);
	//if (err != CL_SUCCESS)
	//	return err;



	status = clfftEnqueueTransform(_clfftPlanHandle, CLFFT_BACKWARD, _numCommandQueues, &_commandQueue, 0, NULL, NULL, &devicePreProcessedSpectra, &deviceFourierTransform, _deviceTemporaryBuffer);
	if (status != CLFFT_SUCCESS)
		return((int)status);


	//err = clFinish(_commandQueue);	// Wait for the Transform to complete
	//if (err != CL_SUCCESS) return err;

	err = 0;
	return err;
}

//
// Post process and produce B-Scan bitmap images with the specified minimum and maximum thresholds
//
DLL PostProcess(float minVal, float maxVal)
{
	cl_int err;
	size_t numWorkItemsPerGroup = _postProcessKernelWorkGroupSize;
	size_t totalWorkItems = _totalAScans;
    //
    // Set the kernel min and max value parameters
    //
    err = clSetKernelArg(_postProcessingKernel, 8, sizeof(cl_float), &minVal);	// min threshold
    if (err != CL_SUCCESS) return err;
    err = clSetKernelArg(_postProcessingKernel, 9, sizeof(cl_float), &maxVal);	// max threshold
    if (err != CL_SUCCESS) return err;
    //
	// Enqueue the post processing kernel on the device - all of the required data, i.e. the FFT result, is already on the GPU - so nothing to copy
	//
	err = clEnqueueNDRangeKernel(_commandQueue, _postProcessingKernel, 1, NULL, &totalWorkItems, &numWorkItemsPerGroup, 0, NULL, NULL);
	if (err != CL_SUCCESS)
		return err;

	//err = clFinish(_commandQueue);
	//if (err != CL_SUCCESS)	return err;


	return 0;


}

//
// Function to compute correlation map on whatever B-Scans have already been loaded on to the GPU
// This version also uses the linear envelope to compute the correlation.
// It might be worth experimenting with other approaches
//
DLL CorrelationMap(
                   unsigned int corrSizeX,
                   unsigned int corrSizeY,
                   unsigned int corrSizeZ,
                   unsigned int offsetX,
                   unsigned int offsetY,
                   unsigned int offsetZ
                   )
{
    cl_int err;
    size_t numWorkItemsPerGroup = 1;
    size_t numKernels = _totalAScans * _outputAScanLength;
    _corrSizeX = corrSizeX;	// Set the window type before calling the kernel...
    _corrSizeY = corrSizeY;
    _corrSizeZ = corrSizeZ;
    _offsetX = offsetX;
    _offsetY = offsetY;
    _offsetZ = offsetZ;
    //
    // Update the kernel arguments
    //
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
    //
    // Enqueue the correlation kernel on the device.  Assign one thread per pixel in the bscan data - each kernel computes one correlation coefficient
    //
    err = clEnqueueNDRangeKernel(_commandQueue, _octCorrelationKernel, 1, NULL, &numKernels, &numWorkItemsPerGroup, 0, NULL, NULL);
    if (err != CL_SUCCESS)
        return err;
    
    //err = clFinish(_commandQueue);
    //if (err != CL_SUCCESS)	return err;
    
    return 0;
}


//
// Function to compute correlation map from B-Scan volume stored in host variable 'volume'
//
// It loads the volume into all available space on the GPU, using multiple batches if necessary.
//
// This version also uses the linear envelope to compute the correlation.
// It might be worth experimenting with other approaches
//
DLL CorrelationMapOnBScanVolume(
                           float* bscanVolume,
                           int numBScansInVolume,   // Ignoring extras for averaging... assume global _bscanAveragingFactor
                           unsigned int corrSizeX,
                           unsigned int corrSizeY,
                           unsigned int corrSizeZ,
                           unsigned int offsetX,
                           unsigned int offsetY,
                           unsigned int offsetZ,
                           float* correlationMapVolume    // Output correlation map volume - must be pre-allocated by calling function
                   )
{
    cl_int err;
    size_t numKernels = _totalAScans * _outputAScanLength;
    size_t numBatches;
    size_t maxBScansPerBatch;
    size_t numBScansInLastBatch;
    size_t batch;
    long startCopyIndex = 0;
    long endCopyIndex;
    
    cl_mem deviceTarget;
    
    int midY;
    
    //
    // Set the target on the device - this is either log or linear bscans
    //
    if (_correlationUsesLogBScans)
        deviceTarget = deviceLogEnvBScan;
    else
        deviceTarget = deviceEnvBScan;
    
    midY=(int)floor((float)corrSizeY/2.0f); // The offset due to performing the correlation over surrounding bscans
    //
    // Work out how many batches are required to compute a complete correlation map
    //
    maxBScansPerBatch = _totalBScans;   // Maximum batch size is the totalBScans value that has already been allocated on the GPU
    numBatches = (size_t)ceil((float)numBScansInVolume/(float)maxBScansPerBatch);
    //
    // Check to see if this can be completed in a single batch, i.e. has enough memory been
    // allocated on the GPU to hold the whole volume in one go?
    //
    if (numBatches > 1)
    {
        // If more than one batch is required, then we need to re-compute the number of batches based on the
        // fact that to form a continuous correlation map, bscans from previous batches need to be incorporated
        // This is done by making sure that the start point of each batch overlaps with the end point
        // from the previous batch.
        //
        // Do this by using an effective batch size, i.e. the maximumBatch size - correlation kernel radius
        //
        numBatches = (size_t)ceil((float)numBScansInVolume/(float)(maxBScansPerBatch-midY));
        
    }
 //   numBScansInLastBatch = maxBScansPerBatch * numBatches - numBScansInVolume;  // CURRENTLY INCORRECT VALUE - TODO: WORKOUT HOW TO DO THIS!!
    //
    // Loop over all batches
    //
    startCopyIndex = 0;
    for (batch=0; batch<numBatches; batch++)
    {
        //
        // For each batch, copy the maximum number of b-scans onto the GPU device
        //
        // Make sure that the start index of each batch overlaps with the previous one to compensate for
        // required overlap in the correlation kernel
        //
        endCopyIndex=startCopyIndex + _totalAScans * _outputAScanLength;
        //
        if (startCopyIndex<0)
            startCopyIndex = 0;
        if (endCopyIndex >= numBScansInVolume * _totalAScansPerBScan * _bscanAveragingFactor)
            endCopyIndex = numBScansInVolume * _totalAScansPerBScan * _bscanAveragingFactor - 1;
        //
        // Copy the batch on to the GPU
        //
        err = clEnqueueWriteBuffer(_commandQueue,
                                   deviceTarget,
                                   CL_TRUE,
                                   0,
                                   sizeof(float) * _totalAScans * _outputAScanLength,
                                   &bscanVolume[startCopyIndex],
                                   0,
                                   NULL,
                                   NULL);
        if (err != CL_SUCCESS)
            return err;
        //
        // Once the volume has been copied, then we can run the correlation kernel
        //
        err = CorrelationMap(corrSizeX, corrSizeY, corrSizeZ, offsetX, offsetY, offsetZ);
        if (err !=CL_SUCCESS)
            return err;
        //
        // Now copy back the correlation map to host memory.
        //
        if (!((batch!=0) && (startCopyIndex==0))) // If batch!=0 and startCopyIndex==0 don't copy because this will be over-written - wasted transfer
        {
            err=CopyCorrelationMapToHost(&correlationMapVolume[startCopyIndex]);
            if (err !=CL_SUCCESS)
                return err;
        }
        //
        startCopyIndex = endCopyIndex - midY * _numAScansPerBScan * _outputAScanLength; // Overlap batches
    }
    
    return 0;
}

/*
//
// Convert B-Scan log data into an image
//
DLL ConvertBScanToBitmap(
                         unsigned char* bmp,
                         float minVal,
                         float maxVal,
                         unsigned int bscanIndex)
{
	int err = 0;
	size_t numWorkItemsPerGroup = 1;

	err = clSetKernelArg(_imageKernel, 4, sizeof(unsigned int), &minVal);
	if (err != CL_SUCCESS) return err;
	err = clSetKernelArg(_imageKernel, 5, sizeof(unsigned int), &maxVal);
	if (err != CL_SUCCESS) return err;
	err = clSetKernelArg(_imageKernel, 6, sizeof(unsigned int), &bscanIndex);
	if (err != CL_SUCCESS) return err;
	//
	// Call the kernel
	//
	err = clEnqueueNDRangeKernel(_commandQueue, _imageKernel, 1, NULL, &_totalAScansPerBScan, &numWorkItemsPerGroup, 0, NULL, NULL);
	if (err != CL_SUCCESS)
		return err;


	//
	// Copy the data back to the host
	//
	err = clEnqueueReadBuffer(_commandQueue, deviceBScanBmp, CL_TRUE, 0, sizeof(unsigned char) * (_totalAScansPerBScan * 3 + _outputImageHeight * _imageFormatStride), bmp, 0, NULL, NULL);
	if (err != CL_SUCCESS) return err;

    err = clFinish(_commandQueue);
    if (err != CL_SUCCESS)	return err;


	return err;
}
*/
//
// Copy the results back to the host from the GPU
//
DLL CopyPreProcessedSpectraToHost(
	float* preProcessedSpectra
	)
{
	cl_int err;
	//
	// Don't block - use Wait function to wait for the buffer to finish copying
	//
	err = clEnqueueReadBuffer(_commandQueue, devicePreProcessedSpectra, CL_TRUE, 0, sizeof(float)*_totalAScans * _inputSpectrumLength*2, preProcessedSpectra, 0, NULL, NULL);
	if (err != CL_SUCCESS) return err;

	return err;
}

DLL CopyComplexFourierTransformToHost(
	float* cmplx
//	float* imagPart
	)
{
	cl_int err;
	//
	// Don't block - use Wait function to wait for the buffer to finish copying
	//
	err = clEnqueueReadBuffer(_commandQueue, deviceFourierTransform, CL_TRUE, 0, sizeof(float)*_totalAScans * _outputAScanLength*2, cmplx, 0, NULL, NULL);
	if (err != CL_SUCCESS) return err;
//	err = clEnqueueReadBuffer(_commandQueue, deviceFourierTransform[1], CL_TRUE, 0, sizeof(float)*_totalAScans * _outputAScanLength, imagPart, 0, NULL, NULL);
//	if (err != CL_SUCCESS) return err;

	return err;

}

DLL CopyLinearEnvelopeToHost(
	float* linearEnvelope
	)
{
	cl_int err;
	/*
cl_mem devicePreProcessedSpectra[2];			// Resampled, windowed, etc, spectra
cl_mem deviceFourierTransform[2];
//cl_mem deviceRealBScan;
//cl_mem deviceImagBScan;
cl_mem deviceEnvBScan;
cl_mem deviceLogEnvBScan;
cl_mem deviceResamplingTable;
cl_mem deviceReferenceSpectrum;
cl_mem deviceReferenceAScan;
cl_mem deviceSum;
cl_mem deviceSAM;
cl_mem deviceAttenuationDepth;


	*/
	//
	// Don't block - use Wait function to wait for the buffer to finish copying
	//
	err = clEnqueueReadBuffer(_commandQueue, deviceEnvBScan, CL_TRUE, 0, sizeof(float)*_totalAScans * _outputAScanLength, linearEnvelope, 0, NULL, NULL);
	if (err != CL_SUCCESS) return err;

	return err;

}

DLL CopyLogEnvelopeToHost(
	float* logEnvelope
	)
{
	cl_int err;
	int i;
	/*
cl_mem devicePreProcessedSpectra[2];			// Resampled, windowed, etc, spectra
cl_mem deviceFourierTransform[2];
//cl_mem deviceRealBScan;
//cl_mem deviceImagBScan;
cl_mem deviceEnvBScan;
cl_mem deviceLogEnvBScan;
cl_mem deviceSum;
cl_mem deviceSAM;
cl_mem deviceAttenuationDepth;


	*/
	//
	// Don't block - use Wait function to wait for the buffer to finish copying
	//
	err = clEnqueueReadBuffer(_commandQueue, deviceLogEnvBScan, CL_TRUE, 0, sizeof(float)*_totalAScans * _outputAScanLength, logEnvelope, 0, NULL, NULL);
	if (err != CL_SUCCESS) return err;

	printf("\n");
	printf("Log Env.\n");
	for (i=0; i<10; i++)
	{
		printf("%f\n",logEnvelope[i]);
	}

	return err;

}

DLL CopyCorrelationMapToHost(
                          float* correlationMap
                          )
{
    cl_int err;
    /*
     cl_mem devicePreProcessedSpectra[2];			// Resampled, windowed, etc, spectra
     cl_mem deviceFourierTransform[2];
     //cl_mem deviceRealBScan;
     //cl_mem deviceImagBScan;
     cl_mem deviceEnvBScan;
     cl_mem deviceLogEnvBScan;
     cl_mem deviceSum;
     cl_mem deviceSAM;
     cl_mem deviceAttenuationDepth;
     
     
     */
    //
    // Don't block - use Wait function to wait for the buffer to finish copying
    //
    err = clEnqueueReadBuffer(_commandQueue, deviceCorrelationMap, CL_TRUE, 0, sizeof(float)*_totalAScans * _outputAScanLength, correlationMap, 0, NULL, NULL);
    if (err != CL_SUCCESS) return err;
    
    return err;
    
}

DLL CopySumToHost(
	float* sum
	)
{
	cl_int err;
	/*
cl_mem devicePreProcessedSpectra[2];			// Resampled, windowed, etc, spectra
cl_mem deviceFourierTransform[2];
//cl_mem deviceRealBScan;
//cl_mem deviceImagBScan;
cl_mem deviceEnvBScan;
cl_mem deviceLogEnvBScan;
cl_mem deviceSum;
cl_mem deviceSAM;
cl_mem deviceAttenuationDepth;


	*/
	//
	// Don't block - use Wait function to wait for the buffer to finish copying
	//
	err = clEnqueueReadBuffer(_commandQueue, deviceSum, CL_TRUE, 0, sizeof(float)*_totalAScans, sum, 0, NULL, NULL);
	if (err != CL_SUCCESS) return err;

	return err;

}

DLL CopySAMToHost(
	float* sam
	)
{
	cl_int err;
	/*
cl_mem devicePreProcessedSpectra[2];			// Resampled, windowed, etc, spectra
cl_mem deviceFourierTransform[2];
//cl_mem deviceRealBScan;
//cl_mem deviceImagBScan;
cl_mem deviceEnvBScan;
cl_mem deviceLogEnvBScan;
cl_mem deviceSum;
cl_mem deviceSAM;
cl_mem deviceAttenuationDepth;


	*/
	//
	// Don't block - use Wait function to wait for the buffer to finish copying
	//
	err = clEnqueueReadBuffer(_commandQueue, deviceSAM, CL_TRUE, 0, sizeof(float)*_totalAScans, sam, 0, NULL, NULL);
	if (err != CL_SUCCESS) return err;

	return err;

}

DLL CopyAttenuationDepthToHost(
	float* attenuationDepth
	)
{
	cl_int err;
	/*
cl_mem devicePreProcessedSpectra[2];			// Resampled, windowed, etc, spectra
cl_mem deviceFourierTransform[2];
//cl_mem deviceRealBScan;
//cl_mem deviceImagBScan;
cl_mem deviceEnvBScan;
cl_mem deviceLogEnvBScan;
cl_mem deviceSum;
cl_mem deviceSAM;
cl_mem deviceAttenuationDepth;


	*/
	//
	// Don't block - use Wait function to wait for the buffer to finish copying
	//
	err = clEnqueueReadBuffer(_commandQueue, deviceAttenuationDepth, CL_TRUE, 0, sizeof(float)*_totalAScans, attenuationDepth, 0, NULL, NULL);
	if (err != CL_SUCCESS) return err;

	return err;

}


//
// Copy B-Scan Bitmaps back to the host
//
DLL CopyBScanBitmapsToHost(
                               unsigned char* bscanBmp
                               )
{
    cl_int err;
    //
    // Don't block - use Wait function to wait for the buffer to finish copying
    //
    err = clEnqueueReadBuffer(_commandQueue, deviceBScanBmp, CL_TRUE, 0, sizeof(unsigned char)*_bitmapBScanVolumeSize, bscanBmp, 0, NULL, NULL);
    if (err != CL_SUCCESS) return err;
    
    return err;
    
}

//
// Copy a single B-Scan Bitmap back to the host
//
DLL CopyBScanBitmapToHost(
								unsigned int bscanIndex,
								unsigned char* bscanBmp
                           )
{
    cl_int err;
    //
    // Don't block - use Wait function to wait for the buffer to finish copying
    //
    err = clEnqueueReadBuffer(_commandQueue, deviceBScanBmp, CL_TRUE, (size_t)bscanIndex * sizeof(unsigned char)*_bitmapBScanSize, sizeof(unsigned char)*_bitmapBScanSize, bscanBmp, 0, NULL, NULL);
    if (err != CL_SUCCESS) return err;
    
    return err;
    
}


//
// Copy all results back to the host - blocking call, automatically calls the clFinish via the Wait function
//
DLL CopyAllResultsToHost(
	float* preProcessedSpectra,
	float* cmplx,
	//float* imagPart,
	float* linearEnvelope,
	float* logEnvelope,
	float* sum,
	float* sam,
	float* attenuationDepth,
    unsigned char* bscanBmp
	)
{
	cl_int err;
	err = CopyPreProcessedSpectraToHost(preProcessedSpectra);
	if (err !=CL_SUCCESS) return err;
	err = CopyComplexFourierTransformToHost(cmplx);
	if (err !=CL_SUCCESS) return err;
	err = CopyLinearEnvelopeToHost(linearEnvelope);
	if (err !=CL_SUCCESS) return err;
	err = CopyLogEnvelopeToHost(logEnvelope);
	if (err !=CL_SUCCESS) return err;
	err = CopySumToHost(sum);
	if (err !=CL_SUCCESS) return err;
	err = CopySAMToHost(sam);
	if (err !=CL_SUCCESS) return err;
	err = CopyAttenuationDepthToHost(attenuationDepth);
	if (err !=CL_SUCCESS) return err;
    err = CopyBScanBitmapsToHost(bscanBmp);
    if (err !=CL_SUCCESS) return err;
	//Wait();

	return err;
}

//
// Allow the user to wait for the command queue to complete execution
//
DLL Wait()
{
	cl_int err=0;
	err = clFinish(_commandQueue);
	return err;
}


DLL SaveBitmap(
               char* path,
               unsigned char* pixelArray,       // In bitmap format (including padding)
               unsigned int width,
               unsigned int height,
               unsigned int pixelArraySize      // length of the above array
               )
{
    //
    // Save BScan to RGB 24bpp device independent bitmap
    //
    float bitmapFileSize;
    unsigned int bpp = 24;
    unsigned int bytesPerPixel = bpp/8;
    unsigned int headerLen = 14;    // 14 byte bitmap header
    unsigned int dibHeaderLen = 12; // BITMAPCOREHEADER
    unsigned int offset = headerLen+dibHeaderLen;
    unsigned int planes = 1;

    unsigned int bytePixel;
    unsigned int i;
	unsigned int ii;
	FILE* fbmp;
    
    bitmapFileSize = headerLen + dibHeaderLen + pixelArraySize;
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
    fwrite(&height, 2, 1, fbmp);
    fwrite(&planes, 2, 1, fbmp);
    fwrite(&bpp, 2, 1, fbmp);
    //
    // Write the pixel array
    //
    for (i=0;i<pixelArraySize; i++)
    {
			ii = pixelArraySize - i - 1;
            fwrite(&pixelArray[ii],1,1,fbmp);  // Red
    }
    
    fclose(fbmp);
    
    return EXIT_SUCCESS;
    
}

//
// Helper function to calculate the image stride
//
DLL CalculateBitmapStride(unsigned int bpp, unsigned int bmpWidth, unsigned int* stride)
{
	unsigned int bytesPerPixel = bpp/8;
	float actualWidthInBytes = (float)(bytesPerPixel * bmpWidth);
	float fstride = (float)ceil((double)(actualWidthInBytes/4.0f)) * 4.0f;
	*stride = (unsigned int)fstride;

	return EXIT_SUCCESS;
}

DLL SaveBitmapFromFloatArray(
               char* path,
               float* sourceArray,
               unsigned short width,
               unsigned short sourceHeight,
               unsigned outputHeight,
               float minVal,
               float maxVal
               )
{
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



