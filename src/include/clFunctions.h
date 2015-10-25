
#define NUM_PLATFORMS 1
#define NUM_DEVICES 2
#define SELECTED_PLATFORM 0
#define MAX_PROGRAM_FILE_LENGTH 10240000
#define DEVICE_INFO_RETURN_SIZE 1024
#define NUM_DEVICES_TO_BUILD_FOR 1	// Number of devices to build a program for - should be 1, just build for the main GPU
#define INTERPOLATION_ORDER 2		// Use second order interpolation
//
// Define the window type constants
//
#define WINDOW_TYPE_NONE	0
#define WINDOW_TYPE_HANN	1
#define WINDOW_TYPE_BLACKMAN 2
#define WINDOW_TYPE_GAUSSIAN 3
//
#include <stdio.h>
#include <stdlib.h>
//#include <CL/cl.h>		// I think that this includes cl.h, so no need to explicitly add
#ifdef _WIN32
#	include <CL/cl.h>
#elif __APPLE__
	#include <OpenCL/OpenCL.h>
#endif
//
// Make sure that we have a boolean type
//
typedef int                 BOOL;

#ifndef FALSE
#define FALSE               0
#endif

#ifndef TRUE
#define TRUE                1
#endif
//
// OpenCL Variables
//
cl_platform_id*		_platformIDs;
cl_device_id*		_deviceIDs;
cl_context			_context;
char**				_deviceNameList;	// Array of devices names
cl_uint				_numDevices;
cl_command_queue	_commandQueue;
cl_uint				_numCommandQueues;
cl_uint				_selectedDeviceIndex;		// Device selected - default it 0... only use 1 device forthe moment
size_t              _maxWorkGroupSize;
size_t              _preProcessKernelWorkGroupSize;
size_t              _postProcessKernelWorkGroupSize;

char*               _compilerOptions;
//size_t				_deviceInfoRetSize;
//cl_uint				_selectedPlatformIndex;		// Platform index that will be used - default is zero
//
// OCT variables
//
cl_uint _inputSpectrumLength;
cl_uint _outputAScanLength;
cl_uint _outputImageHeight;
cl_uint _numBScans;
cl_uint _numAScansPerBScan;
cl_uint _ascanAveragingFactor;
cl_uint _bscanAveragingFactor;
cl_uint _windowType;
cl_uint _imageFormatStride;
cl_uint _corrSizeX;
cl_uint _corrSizeY;
cl_uint _corrSizeZ;
cl_uint _offsetX;
cl_uint _offsetY;
cl_uint _offsetZ;
//
size_t _totalBScans;
size_t _totalAScans;
size_t _totalInputSpectraLength;
size_t _totalOutputAScanLength;
size_t _totalPreProcessedSpectraLength;
size_t _totalAScansPerBScan;
size_t _bitmapBScanVolumeSize;
size_t _bitmapBScanSize;
//
BOOL _correlationUsesLogBScans;
//
//size_t _totalAScans;
//
// Program and kernel variables, each one relates to a different program
//
cl_program			_clOCTProgram;
cl_kernel _preProcessingKernel;
cl_kernel _postProcessingKernel;
cl_kernel _imageKernel;
cl_kernel _octCorrelationKernel;
//cl_kernel dftKernel;
//
//
//
cl_mem deviceSpectra;					// Input array of spectra that make up a single B-Scan
cl_mem devicePreProcessedSpectra;			// Resampled, windowed, etc, spectra
cl_mem deviceFourierTransform;
//cl_mem deviceRealBScan;
//cl_mem deviceImagBScan;
cl_mem deviceEnvBScan;
cl_mem deviceCorrelationMap;        // Stored B-Scan correlation maps
cl_mem deviceLogEnvBScan;
cl_mem deviceResamplingTable;
cl_mem deviceInterpolationMatrix;	// Pre-computed matrix of values used for interpolation
cl_mem deviceReferenceSpectrum;
cl_mem deviceReferenceAScan;
cl_mem deviceSum;
cl_mem deviceSAM;
cl_mem deviceAttenuationDepth;
cl_mem deviceBScanBmp;
//cl_mem deviceSumBmp;




//
// Function prototypes
//

int clInit(
           cl_context* clContext,
           cl_command_queue* clCommandQueue,
           cl_uint deviceIndex,
           char*** deviceNameList,
           cl_uint* numDevicesInList
           );
int clAlloc(
		cl_context context,
		cl_uint inputSpectraLength,		// Specify the actual length of input spectra
		cl_uint outputAScanLength,		// Specify the length of output (outputLength >= inputLength).  If inputlength < outputlength then the input spectra will be zero padded
		cl_uint numBScans,
		cl_uint numAScansPerBScan,
		cl_uint ascanAveragingFactor,
		cl_uint bscanAveragingFactor,
		float* hostResamplingTable,		// Resampling table
		float* hostInterpolationMatrix,		// Pre-computed interpolation matrix
		float* hostReferenceSpectrum,	// Reference spectrum to be subtracted (or maybe divided!!) from each spectrum
		float* hostReferenceAScan,		// A reference background A-Scan to be subtracted.. maybe log domain!!
		unsigned int imageFormatStride
		);
int clRelease();

int clBuild(char* sourceFile, char* buildLog, size_t* buildLogLength);
int clOCTCompileKernels(char* sourceFile, char* build_log, size_t* buildLogLength);
char* clOCTLoadKernelSourceFromFile(char* sourceFile, unsigned int* len);
//int SetImageKernelParameters();
int SetPreProcessingKernelParameters();
int SetPostProcessingKernelParameters();
int SetCorrelationKernelParameters();