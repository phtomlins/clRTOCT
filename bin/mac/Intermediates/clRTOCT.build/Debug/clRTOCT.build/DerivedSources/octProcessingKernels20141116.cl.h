/***** GCL Generated File *********************/
/* Automatically generated file, do not edit! */
/**********************************************/

#include <OpenCL/opencl.h>
extern void (^ConvertFloatArrayToGreyScaleImage_kernel)(const cl_ndrange *ndrange, cl_float* floatImage, cl_uint imageWidth, cl_uint imageHeight, cl_uint imageFormatStride, cl_float minVal, cl_float maxVal, cl_uchar* byteImage);
extern void (^octPostProcessingKernel_kernel)(const cl_ndrange *ndrange, cl_float* FourierTransformCmplx, cl_float* referenceAScan, cl_uint outputAScanLength, cl_uint numBScans, cl_uint numAScansPerBScan, cl_uint ascanAveragingFactor, cl_uint bscanAveragingFactor, cl_float* envBScan, cl_float* logEnvBScan, cl_float* Sum, cl_float* SAM, cl_float* AttenuationDepth);
extern void (^octCorrelationKernel_kernel)(const cl_ndrange *ndrange, cl_float* bscans, cl_uint singleAScanLength, cl_uint numBScans, cl_uint numAScansPerBScan, cl_uint ascanAveragingFactor, cl_uint bscanAveragingFactor, cl_uint corrSizeX, cl_uint corrSizeY, cl_uint corrSizeZ, cl_uint offsetX, cl_uint offsetY, cl_uint offsetZ, cl_float* correlationMap);
extern void (^octPreProcessingKernel_kernel)(const cl_ndrange *ndrange, cl_short* spectra, cl_float* resamplingTable, cl_float* interpolationMatrix, cl_float* referenceSpectrum, cl_uint inputSpectraLength, cl_uint outputAScanLength, cl_uint numBScans, cl_uint numAScansPerBScan, cl_uint ascanAveragingFactor, cl_uint bscanAveragingFactor, cl_uint windowType, cl_float* preProcessedCmplxSpectra);
