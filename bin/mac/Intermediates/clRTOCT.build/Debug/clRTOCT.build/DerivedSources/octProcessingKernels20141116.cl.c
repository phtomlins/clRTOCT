/***** GCL Generated File *********************/
/* Automatically generated file, do not edit! */
/**********************************************/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <dispatch/dispatch.h>
#include <OpenCL/opencl.h>
#include <OpenCL/gcl_priv.h>
#include "octProcessingKernels20141116.cl.h"

static void initBlocks(void);

// Initialize static data structures
static block_kernel_pair pair_map[4] = {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL }
};

static block_kernel_map bmap = { 0, 4, initBlocks, pair_map };

// Block function
void (^ConvertFloatArrayToGreyScaleImage_kernel)(const cl_ndrange *ndrange, cl_float* floatImage, cl_uint imageWidth, cl_uint imageHeight, cl_uint imageFormatStride, cl_float minVal, cl_float maxVal, cl_uchar* byteImage) =
^(const cl_ndrange *ndrange, cl_float* floatImage, cl_uint imageWidth, cl_uint imageHeight, cl_uint imageFormatStride, cl_float minVal, cl_float maxVal, cl_uchar* byteImage) {
  int err = 0;
  cl_kernel k = bmap.map[0].kernel;
  if (!k) {
    initBlocks();
    k = bmap.map[0].kernel;
  }
  if (!k)
    gcl_log_fatal("kernel ConvertFloatArrayToGreyScaleImage does not exist for device");
  kargs_struct kargs;
  gclCreateArgsAPPLE(k, &kargs);
  err |= gclSetKernelArgMemAPPLE(k, 0, floatImage, &kargs);
  err |= gclSetKernelArgAPPLE(k, 1, sizeof(imageWidth), &imageWidth, &kargs);
  err |= gclSetKernelArgAPPLE(k, 2, sizeof(imageHeight), &imageHeight, &kargs);
  err |= gclSetKernelArgAPPLE(k, 3, sizeof(imageFormatStride), &imageFormatStride, &kargs);
  err |= gclSetKernelArgAPPLE(k, 4, sizeof(minVal), &minVal, &kargs);
  err |= gclSetKernelArgAPPLE(k, 5, sizeof(maxVal), &maxVal, &kargs);
  err |= gclSetKernelArgMemAPPLE(k, 6, byteImage, &kargs);
  gcl_log_cl_fatal(err, "setting argument for ConvertFloatArrayToGreyScaleImage failed");
  err = gclExecKernelAPPLE(k, ndrange, &kargs);
  gcl_log_cl_fatal(err, "Executing ConvertFloatArrayToGreyScaleImage failed");
  gclDeleteArgsAPPLE(k, &kargs);
};

void (^octPostProcessingKernel_kernel)(const cl_ndrange *ndrange, cl_float* FourierTransformCmplx, cl_float* referenceAScan, cl_uint outputAScanLength, cl_uint numBScans, cl_uint numAScansPerBScan, cl_uint ascanAveragingFactor, cl_uint bscanAveragingFactor, cl_float* envBScan, cl_float* logEnvBScan, cl_float* Sum, cl_float* SAM, cl_float* AttenuationDepth) =
^(const cl_ndrange *ndrange, cl_float* FourierTransformCmplx, cl_float* referenceAScan, cl_uint outputAScanLength, cl_uint numBScans, cl_uint numAScansPerBScan, cl_uint ascanAveragingFactor, cl_uint bscanAveragingFactor, cl_float* envBScan, cl_float* logEnvBScan, cl_float* Sum, cl_float* SAM, cl_float* AttenuationDepth) {
  int err = 0;
  cl_kernel k = bmap.map[1].kernel;
  if (!k) {
    initBlocks();
    k = bmap.map[1].kernel;
  }
  if (!k)
    gcl_log_fatal("kernel octPostProcessingKernel does not exist for device");
  kargs_struct kargs;
  gclCreateArgsAPPLE(k, &kargs);
  err |= gclSetKernelArgMemAPPLE(k, 0, FourierTransformCmplx, &kargs);
  err |= gclSetKernelArgMemAPPLE(k, 1, referenceAScan, &kargs);
  err |= gclSetKernelArgAPPLE(k, 2, sizeof(outputAScanLength), &outputAScanLength, &kargs);
  err |= gclSetKernelArgAPPLE(k, 3, sizeof(numBScans), &numBScans, &kargs);
  err |= gclSetKernelArgAPPLE(k, 4, sizeof(numAScansPerBScan), &numAScansPerBScan, &kargs);
  err |= gclSetKernelArgAPPLE(k, 5, sizeof(ascanAveragingFactor), &ascanAveragingFactor, &kargs);
  err |= gclSetKernelArgAPPLE(k, 6, sizeof(bscanAveragingFactor), &bscanAveragingFactor, &kargs);
  err |= gclSetKernelArgMemAPPLE(k, 7, envBScan, &kargs);
  err |= gclSetKernelArgMemAPPLE(k, 8, logEnvBScan, &kargs);
  err |= gclSetKernelArgMemAPPLE(k, 9, Sum, &kargs);
  err |= gclSetKernelArgMemAPPLE(k, 10, SAM, &kargs);
  err |= gclSetKernelArgMemAPPLE(k, 11, AttenuationDepth, &kargs);
  gcl_log_cl_fatal(err, "setting argument for octPostProcessingKernel failed");
  err = gclExecKernelAPPLE(k, ndrange, &kargs);
  gcl_log_cl_fatal(err, "Executing octPostProcessingKernel failed");
  gclDeleteArgsAPPLE(k, &kargs);
};

void (^octCorrelationKernel_kernel)(const cl_ndrange *ndrange, cl_float* bscans, cl_uint singleAScanLength, cl_uint numBScans, cl_uint numAScansPerBScan, cl_uint ascanAveragingFactor, cl_uint bscanAveragingFactor, cl_uint corrSizeX, cl_uint corrSizeY, cl_uint corrSizeZ, cl_uint offsetX, cl_uint offsetY, cl_uint offsetZ, cl_float* correlationMap) =
^(const cl_ndrange *ndrange, cl_float* bscans, cl_uint singleAScanLength, cl_uint numBScans, cl_uint numAScansPerBScan, cl_uint ascanAveragingFactor, cl_uint bscanAveragingFactor, cl_uint corrSizeX, cl_uint corrSizeY, cl_uint corrSizeZ, cl_uint offsetX, cl_uint offsetY, cl_uint offsetZ, cl_float* correlationMap) {
  int err = 0;
  cl_kernel k = bmap.map[2].kernel;
  if (!k) {
    initBlocks();
    k = bmap.map[2].kernel;
  }
  if (!k)
    gcl_log_fatal("kernel octCorrelationKernel does not exist for device");
  kargs_struct kargs;
  gclCreateArgsAPPLE(k, &kargs);
  err |= gclSetKernelArgMemAPPLE(k, 0, bscans, &kargs);
  err |= gclSetKernelArgAPPLE(k, 1, sizeof(singleAScanLength), &singleAScanLength, &kargs);
  err |= gclSetKernelArgAPPLE(k, 2, sizeof(numBScans), &numBScans, &kargs);
  err |= gclSetKernelArgAPPLE(k, 3, sizeof(numAScansPerBScan), &numAScansPerBScan, &kargs);
  err |= gclSetKernelArgAPPLE(k, 4, sizeof(ascanAveragingFactor), &ascanAveragingFactor, &kargs);
  err |= gclSetKernelArgAPPLE(k, 5, sizeof(bscanAveragingFactor), &bscanAveragingFactor, &kargs);
  err |= gclSetKernelArgAPPLE(k, 6, sizeof(corrSizeX), &corrSizeX, &kargs);
  err |= gclSetKernelArgAPPLE(k, 7, sizeof(corrSizeY), &corrSizeY, &kargs);
  err |= gclSetKernelArgAPPLE(k, 8, sizeof(corrSizeZ), &corrSizeZ, &kargs);
  err |= gclSetKernelArgAPPLE(k, 9, sizeof(offsetX), &offsetX, &kargs);
  err |= gclSetKernelArgAPPLE(k, 10, sizeof(offsetY), &offsetY, &kargs);
  err |= gclSetKernelArgAPPLE(k, 11, sizeof(offsetZ), &offsetZ, &kargs);
  err |= gclSetKernelArgMemAPPLE(k, 12, correlationMap, &kargs);
  gcl_log_cl_fatal(err, "setting argument for octCorrelationKernel failed");
  err = gclExecKernelAPPLE(k, ndrange, &kargs);
  gcl_log_cl_fatal(err, "Executing octCorrelationKernel failed");
  gclDeleteArgsAPPLE(k, &kargs);
};

void (^octPreProcessingKernel_kernel)(const cl_ndrange *ndrange, cl_short* spectra, cl_float* resamplingTable, cl_float* interpolationMatrix, cl_float* referenceSpectrum, cl_uint inputSpectraLength, cl_uint outputAScanLength, cl_uint numBScans, cl_uint numAScansPerBScan, cl_uint ascanAveragingFactor, cl_uint bscanAveragingFactor, cl_uint windowType, cl_float* preProcessedCmplxSpectra) =
^(const cl_ndrange *ndrange, cl_short* spectra, cl_float* resamplingTable, cl_float* interpolationMatrix, cl_float* referenceSpectrum, cl_uint inputSpectraLength, cl_uint outputAScanLength, cl_uint numBScans, cl_uint numAScansPerBScan, cl_uint ascanAveragingFactor, cl_uint bscanAveragingFactor, cl_uint windowType, cl_float* preProcessedCmplxSpectra) {
  int err = 0;
  cl_kernel k = bmap.map[3].kernel;
  if (!k) {
    initBlocks();
    k = bmap.map[3].kernel;
  }
  if (!k)
    gcl_log_fatal("kernel octPreProcessingKernel does not exist for device");
  kargs_struct kargs;
  gclCreateArgsAPPLE(k, &kargs);
  err |= gclSetKernelArgMemAPPLE(k, 0, spectra, &kargs);
  err |= gclSetKernelArgMemAPPLE(k, 1, resamplingTable, &kargs);
  err |= gclSetKernelArgMemAPPLE(k, 2, interpolationMatrix, &kargs);
  err |= gclSetKernelArgMemAPPLE(k, 3, referenceSpectrum, &kargs);
  err |= gclSetKernelArgAPPLE(k, 4, sizeof(inputSpectraLength), &inputSpectraLength, &kargs);
  err |= gclSetKernelArgAPPLE(k, 5, sizeof(outputAScanLength), &outputAScanLength, &kargs);
  err |= gclSetKernelArgAPPLE(k, 6, sizeof(numBScans), &numBScans, &kargs);
  err |= gclSetKernelArgAPPLE(k, 7, sizeof(numAScansPerBScan), &numAScansPerBScan, &kargs);
  err |= gclSetKernelArgAPPLE(k, 8, sizeof(ascanAveragingFactor), &ascanAveragingFactor, &kargs);
  err |= gclSetKernelArgAPPLE(k, 9, sizeof(bscanAveragingFactor), &bscanAveragingFactor, &kargs);
  err |= gclSetKernelArgAPPLE(k, 10, sizeof(windowType), &windowType, &kargs);
  err |= gclSetKernelArgMemAPPLE(k, 11, preProcessedCmplxSpectra, &kargs);
  gcl_log_cl_fatal(err, "setting argument for octPreProcessingKernel failed");
  err = gclExecKernelAPPLE(k, ndrange, &kargs);
  gcl_log_cl_fatal(err, "Executing octPreProcessingKernel failed");
  gclDeleteArgsAPPLE(k, &kargs);
};

// Initialization functions
static void initBlocks(void) {
  const char* build_opts = "";
  static dispatch_once_t once;
  dispatch_once(&once,
    ^{ int err = gclBuildProgramBinaryAPPLE("OpenCL/octProcessingKernels20141116.cl", "", &bmap, build_opts);
       if (!err) {
          assert(bmap.map[0].block_ptr == ConvertFloatArrayToGreyScaleImage_kernel && "mismatch block");
          bmap.map[0].kernel = clCreateKernel(bmap.program, "ConvertFloatArrayToGreyScaleImage", &err);
          assert(bmap.map[1].block_ptr == octPostProcessingKernel_kernel && "mismatch block");
          bmap.map[1].kernel = clCreateKernel(bmap.program, "octPostProcessingKernel", &err);
          assert(bmap.map[2].block_ptr == octCorrelationKernel_kernel && "mismatch block");
          bmap.map[2].kernel = clCreateKernel(bmap.program, "octCorrelationKernel", &err);
          assert(bmap.map[3].block_ptr == octPreProcessingKernel_kernel && "mismatch block");
          bmap.map[3].kernel = clCreateKernel(bmap.program, "octPreProcessingKernel", &err);
       }
     });
}

__attribute__((constructor))
static void RegisterMap(void) {
  gclRegisterBlockKernelMap(&bmap);
  bmap.map[0].block_ptr = ConvertFloatArrayToGreyScaleImage_kernel;
  bmap.map[1].block_ptr = octPostProcessingKernel_kernel;
  bmap.map[2].block_ptr = octCorrelationKernel_kernel;
  bmap.map[3].block_ptr = octPreProcessingKernel_kernel;
}

