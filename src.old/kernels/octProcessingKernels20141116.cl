//
#define NULL 0
//
// Define the window type constants
//
#define WINDOW_TYPE_NONE	0
#define WINDOW_TYPE_HANN	1
#define WINDOW_TYPE_BLACKMAN 2
#define WINDOW_TYPE_GAUSSIAN 3
#define MAX_FFT_LEN 1024

void GPUResample(
	__global float* spectra, 
	__global const float* resamplingTable, 
	__global const float* interpolationMatrix, 
	unsigned int inputSpectraLength, 
	__global float* resampledSpectra);
void preProcess(
	__global const short* spectrum,  
	__global const float* resamplingTable,
	__global const float* interpolationMatrix,
	__global const float* referenceSpectrum, 
	__global float* preProcessedCmplxSpectrum, 
	unsigned int inputSpectraLength,
	unsigned int outputAScanLength,
	unsigned int windowType
	);


//
// Kernel to convert an array of floats into a 256 greyscale byte array image.  The idea is that we pass in the
// data part of a managed bitmap object
//
__kernel void ConvertFloatArrayToGreyScaleImage(
											__global const float* floatImage,			// 0.  2D float array to convert into an image
											const unsigned int imageWidth,				// 1.  Width of image
											const unsigned int imageHeight,				// 2.  Height of the image
											const unsigned int imageFormatStride,		// 3.  Number of bytes required to store a single row in the image
											const float minVal,							// 4. Minimum image intensity - used for scaling
											const float maxVal,							// 5. Maximum image intensity for scaling
											__global unsigned char* byteImage			// 6. Output image array of bytes
											)
{
	//
	//
	//
	const unsigned int x = get_global_id(0);		// A-Scan for each instance of the kernel to operate on.
	unsigned int y;
	float span = maxVal - minVal;
	float val;
	unsigned char bVal;
    // Go through the draw area and set the pixels as they should be
	if (x<imageWidth)
	{
		for (int y = 0; y < imageHeight; y++)
		{
			val = (floatImage[x*imageHeight*2 + y] - minVal)/span * 255.0f;    // Transposed image as source - remember image height of source data is double actual image height (assuming FFT as source)
			bVal = (unsigned char)val;
			if (val > 255.0f)
				bVal = 255;
			if (val < 0.0f)
				bVal = 0;
                            
                            

			// layer.GetBitmap().SetPixel(x, y, m_colour);
			byteImage[(x * 3) + y * imageFormatStride] = bVal;
			byteImage[(x * 3) + y * imageFormatStride + 1] = bVal;
			byteImage[(x * 3) + y * imageFormatStride + 2] = bVal;
		}
	}

}


//
// Post processing kernel that produces the OCT output
//
__kernel void octPostProcessingKernel(
						__global const float* FourierTransformCmplx,				// 1. Input, the Fourier transform of the pre-processed spectra	
					//	__global const float* FourierTransformImag,				// 1. Input, the Fourier transform of the pre-processed spectra	
						__global const float* referenceAScan,				// 2. Reference a-scan subtracted from envelope of FT
						const unsigned int outputAScanLength,				// 3. Length of the 1D Fourier transform			
						const unsigned int numBScans,						// 4.
						const unsigned int numAScansPerBScan,				
						const unsigned int ascanAveragingFactor,
						const unsigned int bscanAveragingFactor,
						__global float* envBScan,							
						__global float* logEnvBScan,						
						__global float* Sum,
						__global float* SAM,		// Currently unused
						__global float* AttenuationDepth	// Currently unused
						)
{
	//
	// Reconstruct the envelope
	//
	float re2;
	float im2;
	float sum = 0.0f;
	int i,count;
	int topOffset = 15;
    const unsigned int ascanIndex = get_global_id(0);	// numAScansPerBScan * ascanAveragingFactor * numBScans * bscanAveragingFactor  work items should be exectuting this kernel,
	// i.e. one work item per a-scan.

	const unsigned int currentAScanOffset = ascanIndex * outputAScanLength;	// Index of beginning of current A-Scan in spectra array


	for (i=0; i<outputAScanLength; i++)
	{
		re2 = FourierTransformCmplx[2*(currentAScanOffset+i)] * FourierTransformCmplx[2*(currentAScanOffset+i)];
		im2 = FourierTransformCmplx[2*(currentAScanOffset+i)  + 1] * FourierTransformCmplx[2*(currentAScanOffset+i) + 1];
		envBScan[currentAScanOffset+i] =  re2 + im2 ;
		logEnvBScan[currentAScanOffset + i] = 20.0f*(float)log10(envBScan[currentAScanOffset+i] );
		//sum +=envBScan[currentAScanOffset+i];
	}
	
	count = 0;
	for (i=topOffset; i<outputAScanLength/8; i++)
	{
		sum += envBScan[currentAScanOffset+i];
		count++;
	}

	Sum[ascanIndex] = sum / (float)count;
	//

}



//
// Kernel for calculating correlation mapped images on the gpu
//
__kernel void octCorrelationKernel(
                                     __global float* bscans,					// Input bscans
                                     const unsigned int singleAScanLength,
                                     const unsigned int numBScans,  // Number of individual bscans (not counting duplicates for averaging)
                                     const unsigned int numAScansPerBScan,  // Ignoring extras for averaging
                                     const unsigned int ascanAveragingFactor,
                                     const unsigned int bscanAveragingFactor,
                                     const unsigned int corrSizeX,  // Correlation window size across b-scan
                                     const unsigned int corrSizeY,  // Correlation window size across multiple b-scans (<=numBScans)
                                     const unsigned int corrSizeZ,  // Correlation window size axially
                                     const unsigned int offsetX,      // Offset between the centre of regions being correlated
                                     const unsigned int offsetY,      // Should allow elastographic analysis
                                     const unsigned int offsetZ,      //
                                     __global float* correlationMap				// Correlation coefficients
                                     )

{
    
    //
    // Find the Pearson correlation coefficient between adjacent OCT B-Scans
    //
    // Each kernel calculates the correlation coeficient of a single region and stores its result
    // in the associated pixel of the correlation map array
    //
    // Kernels are launched using a 1D indexing scheme such that sequential indices (work items)
    // operate on axially sequential correlation coefficients in the bscan correlation map, i.e.
    // kernelIndex=0 -> correlationMap[0]
    // kernelIndex=1 -> correlationMap[1]
    const int kernelIndex = get_global_id(0);	// index of current work item
    //
    // Convert the kernelIndex into ascanIndex, bscansIndex and axialIndex values
    //
   // unsigned int totalBScans = numBScans * bscanAveragingFactor;
    int singleBScanLength = singleAScanLength * numAScansPerBScan;
    //
    int totalAScanLength = singleAScanLength * ascanAveragingFactor;
    int totalBScanLength = singleBScanLength * ascanAveragingFactor * bscanAveragingFactor;
    int totalPixels = totalBScanLength * numBScans;
    int totalBScans = numBScans * bscanAveragingFactor;
    int totalAScansPerBScan = numAScansPerBScan * ascanAveragingFactor;
    int ascanIndex;
    int bscanIndex;
    int axialIndex;
    
    int x;
    int y;
    int z;
    int x1;
    int y1;
    int z1;
    int x2;
    int y2;
    int z2;
    int srcIndex1=0; // Source bscan array index
    int srcIndex2=0;
    int numPoints=0; // target array index
    int i;
    
    int midX=(int)floor((float)corrSizeX/2.0f);
    int midY=(int)floor((float)corrSizeY/2.0f);
    int midZ=(int)floor((float)corrSizeZ/2.0f);
    
    float corr=0.0f;
    float meanVal1=0.0f;
    float meanVal2=0.0f;
    float sumVal=0.0f;
    
    float diff1=0.0f;
    float diff2= 0.0f;
    float sumDiff1Diff2 = 0.0f;
    float sumDiff1Diff1 = 0.0f;
    float sumDiff2Diff2 = 0.0f;

    //
    // Determine the ascan, bscan and axial indices pointed to by the kernel index
    //
    bscanIndex = (int)floor((float)kernelIndex/(float)singleBScanLength);
    ascanIndex = (int)floor((float)(kernelIndex - bscanIndex * singleBScanLength)/(float)singleAScanLength);
    axialIndex = kernelIndex - bscanIndex * singleBScanLength - ascanIndex * singleAScanLength;
    //
    //
    if ((bscanIndex-midY >= 0) && (bscanIndex-midY+corrSizeY+offsetY-1 < totalBScans))    // Make sure that we don't try to use non-existent b-scans
    {
        for (y=0; y<corrSizeY;y++)  // Loop over b-scans
        {
            y1=bscanIndex-midY + y;
            y2=y1 + offsetY;
            //
            for (x=0;x<corrSizeX;x++)   // Loop over intra-bscan lateral position
            {
                x1=ascanIndex-midX + x;
                x2=x1 + offsetX;
                
                for (z=0; z<corrSizeZ;z++)  // Loop over axial position
                {
                    z1=axialIndex-midZ + z;
                    z2=z1+offsetZ;
                    srcIndex1 = y1*singleBScanLength + x1*singleAScanLength + z1;
                    srcIndex2 = y2*singleBScanLength + x2*singleAScanLength + z2;
                    //
                    // Check that the index is within a valid bscan and ascan
                    //
                    if (
                            ((0<=z1) && (z2<singleAScanLength)) &&
                            ((0<=y1) && (y2<totalBScans)) &&
                            ((0<=x1) && (x2<totalAScansPerBScan)) &&
                            ((0<=srcIndex1) && (srcIndex1<totalPixels)) &&
                            ((0<=srcIndex2) && (srcIndex2<totalPixels))
                        )
                    {
                        meanVal1 += bscans[srcIndex1];
                        meanVal2 += bscans[srcIndex2];
                        numPoints++;
                    }
                }
            }
        }
        meanVal1 = meanVal1/(float)numPoints;
        meanVal2 = meanVal2/(float)numPoints;
        //
        // Loop back through the ROIs to compute the components of the correlation coefficient
        //
        for (y=0; y<corrSizeY;y++)  // Loop over b-scans
        {
            y1=bscanIndex-midY + y;
            y2=y1 + offsetY;
            //
            for (x=0;x<corrSizeX;x++)   // Loop over intra-bscan lateral position
            {
                x1=ascanIndex-midX + x;
                x2=x1 + offsetX;
                
                for (z=0; z<corrSizeZ;z++)  // Loop over axial position
                {
                    z1=axialIndex-midZ + z;
                    z2=z1+offsetZ;
                    srcIndex1 = y1*singleBScanLength + x1*singleAScanLength + z1;
                    srcIndex2 = y2*singleBScanLength + x2*singleAScanLength + z2;
                    //
                    // Check that the index is within a valid bscan and ascan
                    //
                    if (
                        ((0<=z1) && (z2<singleAScanLength)) &&
                        ((0<=y1) && (y2<totalBScans)) &&
                        ((0<=x1) && (x2<totalAScansPerBScan)) &&
                        ((0<=srcIndex1) && (srcIndex1<totalPixels)) &&
                        ((0<=srcIndex2) && (srcIndex2<totalPixels))
                        )
                    {
                        //
                        // Calculate the difference between each point and the mean
                        //
                        diff1 = bscans[srcIndex1]-meanVal1;
                        diff2 = bscans[srcIndex2]-meanVal2;
                        //
                        // Calculate the sum of the product terms
                        //
                        sumDiff1Diff2 += diff1 * diff2;
                        sumDiff1Diff1 += diff1 * diff1;
                        sumDiff2Diff2 += diff2 * diff2;
                    }
                }
            }
        }
        //
        corr = sumDiff1Diff2/(sqrt(sumDiff1Diff1)*sqrt(sumDiff2Diff2));
        
        correlationMap[kernelIndex] = corr;
        
    
    }
    else
    {
        correlationMap[kernelIndex]=-2.0;   // i.e. correlation not computed (invalid correlation value)
    }

}


// raw spectra
// An array that will hold the resampled spectra
// The real part of the FFT
// The imaginary part of the FFT
// The envelope of the BScan
// The log of the envelope
// Resampling table indicating the index at which to resample each point in the spectra
// Length of a single input spectrum
// Output length of a single A-Scan - should be the length of realBScan, imagBScan and logEnvBScan
// Number of BScans in spectra
// Number of AScans per BScan	(ignoring averaging)
// Number of AScans to group together for averaging

__kernel void octPreProcessingKernel(
						__global const short* spectra,					// Input arrays - raw spectra					
						__global const float* resamplingTable,			// resampling table
						__global const float* interpolationMatrix,
						__global const float* referenceSpectrum,		// reference spectrum for deconvolution
						const unsigned int inputSpectraLength,				
						const unsigned int outputAScanLength,				
						const unsigned int numBScans,						
						const unsigned int numAScansPerBScan,				
						const unsigned int ascanAveragingFactor,
						const unsigned int bscanAveragingFactor,
						const unsigned int windowType,
						__global float* preProcessedCmplxSpectra				// Output is the pre-processed spectra, ready for FFT
						)

{
	
	//
	// PreProcess OCT sepctra ready for input to clFFT
	//
	// Get the global id to specify which ascan is being transformed - operate on B-Scans only
	//
    const unsigned int ascanIndex = get_global_id(0);	// numAScansPerBScan * averagingFactor work items should be exectuting this kernel,
	// i.e. one work item per a-scan.  These are executed in averagingFactor work-groups
	// Use a 2D addressing system such that bscans are indexed on the second dimension

	const unsigned int currentSpectraOffset = ascanIndex * inputSpectraLength;	// Index of beginning of current A-Scan in spectra array
	const unsigned int currentAScanOffset = ascanIndex * outputAScanLength;	// Index of beginning of current A-Scan in spectra array
	int unsigned i = 0 ;
	int unsigned j = 0;

//	float M = (float)outputAScanLength;
//	float N = (float)inputSpectraLength;
	//
	// Apply pre-processing
	//
	preProcess(
		&spectra[currentSpectraOffset],		// To spectrum starting at currentSpectraOffset
		resamplingTable,					// Use the resampling table
		interpolationMatrix,			// Pre-computed coeffients for faster interpolation
		referenceSpectrum,					// Divide by the reference spectrum
		&preProcessedCmplxSpectra[2*currentAScanOffset],	// Output the pre-processing at the current AScan offset, to accomodate any zero padding
		inputSpectraLength,
		outputAScanLength,
		windowType);

}



//
// Polynomial resampling using linear interpolation
//
void GPUResample(
	__global float* spectrum, 
	__global const float* resamplingTable, 
	__global const float* interpolationMatrix, 
	unsigned int inputSpectraLength, 
	__global float* resampledSpectra
	)
{
    //
    // Resample the spectra using the GPU
    //
    int i,r,c;
    int x;
    int o;
	int offset;
	int interpOrder = 2;	// Assume quadratic interpolation (2nd order)
	int numInterpCoefs = (interpOrder+1)*(interpOrder+1);
	float p[3];				// Polynomial coefficients
	float Y[3];				// Used for interpolation
	float x0;
    float x3;    // Index at which to resample the spectrum
    float x1;
    float x2;
	float I0;
    float I1;
    float I2;
    float I3;
    //
    for (x = 0; x < inputSpectraLength; x++)    // Resample pixel by pixel
    {
        x0 = resamplingTable[x];
        //
        // Round to find the closest pixel indices
        //
        x1 = (float)floor(x0);
        x2 = (float)ceil(x0);
		x3 = x2+1.0f;
        //
        // Loop through the a-scans
        //
        //
        // Check that x1 and x3 are not out of the array bounds
        //
        if ((int)x1 < 0)								// x1 must not be less than zero
            I0 = (float)spectrum[0];
        else if ((int)x2 > inputSpectraLength - 1)		// and x2 must not be greater than the spectrum length
            I0 = (float)spectrum[2*(inputSpectraLength - 1)];
        else if (fabs(x2 - x1) < 0.5)					// Also check that x1 != x2.  If there are equal, then don't interplation
            I0 = (float)spectrum[2*(int)x3];
        else
        {
            //
            // Interpolate to estimate I0 at x0 using pre-computed interpolation coefficient matrix
            // First compute the polynomial coefficients a,b,c etc...
			// 
            I1 = (float)spectrum[2*(int)x1];
            I2 = (float)spectrum[2*(int)x2];
			I3 = (float)spectrum[2*(int)x3];
			Y[0] = I1;
			Y[1] = I2;
			Y[2] = I3;
			//
			// Multiple the coefficient matrix with the polynomial coefficient vector
			//
			offset = x*numInterpCoefs;
			//
			i=0;
			for (r=0; r<=interpOrder; r++)
			{
				p[r] = 0.0f;
				for (c=0; c <= interpOrder; c++)
				{
					p[r] += interpolationMatrix[offset+i] * Y[c];
					i++;
				}
			}
			//
			// Now compute the interpolated value of I0 at x0
			//
			I0 = p[0] + p[1] * x0 + p[2] * x0 * x0;
        }
        resampledSpectra[2*x] = I0;
    }
}

void preProcess(
	__global const short* spectrum,  
	__global const float* resamplingTable,
	__global const float* interpolationMatrix,
	__global const float* referenceSpectrum, 
	__global float* preProcessedCmplxSpectrum, 
	unsigned int inputSpectraLength,
	unsigned int outputAScanLength,
	unsigned int windowType
	)
{
	//
	// Divide reference spectrum by the measured spectra
	//
	unsigned int i;
	float N = (float)inputSpectraLength;
	float a0;
	float a1;
	float a2;
	float a3;

	for (i=0; i<inputSpectraLength; i++)
	{
		//
		// Set the imaginary part to zero
		//
		preProcessedCmplxSpectrum[2*i+1] = 0.0f;
		//
		// Deconvolve the reference spectrum from the real part
		//
		if (referenceSpectrum != NULL)
		{
			if ( referenceSpectrum[i] > 0.001 )
			{
				preProcessedCmplxSpectrum[2*i] = (float)spectrum[i] / referenceSpectrum[i] - 1.0f;
			}
			else
			{
				preProcessedCmplxSpectrum[2*i] = 0.0f;
			}
		}
		else
		{
			//
			// If no reference spectrum, then ignore deconvolution
			//
			preProcessedCmplxSpectrum[2*i] = (float)spectrum[i];
		}
	}

	//
	// Apply a window function
	//
	if (windowType == WINDOW_TYPE_HANN)
	{
		for (i=0; i<inputSpectraLength; i++)
		{
			preProcessedCmplxSpectrum[2*i] *= 0.5f*(1-cos(2*M_PI_F*i/(inputSpectraLength-1)));
		}
	}
	else if (windowType == WINDOW_TYPE_BLACKMAN)
	{
		a0 = 0.355768;
		a1 = 0.487396;
		a2 = 0.144232;
		a3 = 0.012604;
		for (i=0; i<inputSpectraLength; i++)
		{
			preProcessedCmplxSpectrum[2*i] *= a0 - a1 * cos(2*M_PI_F*i/(N-1.0f)) + a2 * cos(4*M_PI_F*i/(N-1.0f)) - a3 * cos(6*M_PI_F*i/(N-1.0f));
		}
	
	}
	else if (windowType == WINDOW_TYPE_GAUSSIAN)
	{
		a0 = N/2.0f;
		a1 = N/4.0f;
		a2 = a1*a1;
		for (i=0; i<inputSpectraLength; i++)
		{
			a3 = (float)i-a0;
			preProcessedCmplxSpectrum[2*i] *= exp(-(a3*a3)/a2);
		}
	
	}
	//
	// Use the imaginary part as a workspace
	//
	for (i=0; i<inputSpectraLength; i++)
	{
		preProcessedCmplxSpectrum[2*i+1] = preProcessedCmplxSpectrum[2*i];
	}
	//
	// Resample
	//
	//
	// If a resampling table has been provided, then resample the spectra - only the real part
	//
	if (resamplingTable != NULL)
	{
		 //GPUResample(preProcessedCmplxSpectrum, resamplingTable, inputSpectraLength, preProcessedCmplxSpectrum);
		 GPUResample(&preProcessedCmplxSpectrum[1], resamplingTable, interpolationMatrix, inputSpectraLength, preProcessedCmplxSpectrum);
	}
	//
	// Set imaginary part to zero
	//
	for (i=0; i<inputSpectraLength; i++)
	{
		preProcessedCmplxSpectrum[2*i+1] = 0.0f;
	}

	//
	// Finally, zero pad if necessary
	//
	if (outputAScanLength > inputSpectraLength)
	{
		for (i=inputSpectraLength-1; i<outputAScanLength; i++)
		{
			preProcessedCmplxSpectrum[2*i] = 0.0f;
			preProcessedCmplxSpectrum[2*i+1] = 0.0f;

		}
	}


}
