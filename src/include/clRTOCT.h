
#define DLL extern __declspec(dllexport) int __stdcall
//
// See http://stackoverflow.com/questions/538134/exporting-functions-from-a-dll-with-dllexport
//
#ifdef _WIN32
#   ifdef LIBRARY_EXPORTS
#       define LIBRARY_API __declspec(dllexport)
#   else
#       define LIBRARY_API __declspec(dllimport)
#   endif// Define the window type constants
#elif __APPLE__
#    define LIBRARY_API
#endif
//
#define WINDOW_TYPE_NONE	0
#define WINDOW_TYPE_HANN	1
#define WINDOW_TYPE_BLACKMAN 2
#define WINDOW_TYPE_GAUSSIAN 3

//#include <float.h>
#include <math.h>
#include <clFFT.h>
#include "clFunctions.h"
//
// Create some global variables
//
// clFFT
//
clfftSetupData _clfftSetupData;
clfftPlanHandle _clfftPlanHandle;
size_t _fftLength;
char* _kernelPath;

size_t _tempBufferSize;
cl_mem _deviceTemporaryBuffer;

//
// OCT Variables
//
float* _interpolationMatrix;

//
// Exported Function prototypes
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
);
DLL clOCTDispose();

void PreComputeInterpolationCoefficients(float* resamplingTable, float* coefs, int ResamplingTableLength);	// Internal function

DLL PreProcess(short* hostSpectra, int windowType);
DLL InverseTransform();
DLL PostProcess(float minVal, float maxVal);
//DLL ConvertBScanToBitmap(unsigned char* bmp, float minVal, float maxVal, unsigned int bscanIndex);
DLL CorrelationMap(
                   unsigned int corrSizeX,
                   unsigned int corrSizeY,
                   unsigned int corrSizeZ,
                   unsigned int offsetX,
                   unsigned int offsetY,
                   unsigned int offsetZ
                   );

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
);
DLL CopyPreProcessedSpectraToHost(
	float* preProcessedSpectra
	);
DLL CopyComplexFourierTransformToHost(
	float* cmplx
	//float* imagPart
	);
DLL CopyLinearEnvelopeToHost(
	float* linearEnvelope
	);
DLL CopyLogEnvelopeToHost(
	float* logEnvelope
	);
DLL CopyCorrelationMapToHost(
                          float* correlationMap
                          );
DLL CopySumToHost(
	float* sum
	);
DLL CopySAMToHost(
	float* sam
	);
DLL CopyAttenuationDepthToHost(
	float* attenuationDepth
	);

DLL CopyBScanBitmapsToHost(
                           unsigned char* bscanBmp
                           );

DLL CopyBScanBitmapToHost(
								unsigned int bscanIndex,
								unsigned char* bscanBmp
                           );

DLL CopyAllResultsToHost(
	float* preProcessedSpectra,
	float* cmplx,
//	float* imagPart,
	float* linearEnvelope,
	float* logEnvelope,
	float* sum,
	float* sam,
	float* attenuationDepth,
    unsigned char* bscanBmp
	);
DLL Wait();
DLL SaveBitmap(
               char* path,
               unsigned char* pixelArray,       // In bitmap format (including padding)
               unsigned int width,
               unsigned int height,
               unsigned int pixelArraySize      // length of the above array
                );
DLL CalculateBitmapStride(unsigned int bpp, unsigned int bmpWidth, unsigned int* stride);

DLL SaveBitmapFromFloatArray(
                             char* path,
                             float* sourceArray,
                             unsigned short width,
                             unsigned short sourceHeight,
                             unsigned outputHeight,
                             float minVal,
                             float maxVal
                             );

