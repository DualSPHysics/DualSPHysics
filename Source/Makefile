#DualSPHysics v4.2.023 Full 28-11-2017

#=============== CUDA selection ===============
CUDAVER=75
#=============== Cluster selection ===============
CLUSTER=EPH
### EPH -> Uses g++ compiler
### EPH -> To use Cuda 4.0: module load cuda/7.5
### CSF -> Uses g++ compiler
### CSF -> To use Cuda 4.0: module load libs/cuda/4.0.17
### BSC -> Uses icpc compiler 
### BSC -> To use Cuda 4.0: module unload cuda/4.1 && module load cuda/4.0
### EME -> To use icpc: module load intel/12.1
### EME -> To use Cuda 4.0: module load cuda/4.0.17

#=============== Files to compile ===============
OBJXML=JXml.o tinystr.o tinyxml.o tinyxmlerror.o tinyxmlparser.o
OBJSPHMOTION=JMotion.o JMotionList.o JMotionMov.o JMotionObj.o JMotionPos.o JSphMotion.o
OBCOMMON=Functions.o FunctionsMath.o JBinaryData.o JException.o JLog2.o JMeanValues.o JObject.o JRadixSort.o JRangeFilter.o JReadDatafile.o JSaveCsv.o JTimeControl.o randomc.o
OBCOMMONDSPH=JDsphConfig.o JPartDataBi4.o JPartFloatBi4.o JPartOutBi4Save.o JSpaceCtes.o JSpaceEParms.o JSpaceParts.o JSpaceProperties.o
OBSPH=JArraysCpu.o JCellDivCpu.o JCfgRun.o JDamping.o JGauge.o JPartsOut.o JSaveDt.o JSph.o JSphAccInput.o JSphCpu.o JSphDtFixed.o JSphVisco.o JTimeOut.o main.o
OBSPHSINGLE=JCellDivCpuSingle.o JPartsLoad4.o JSphCpuSingle.o
OBCOMMONGPU=FunctionsCuda.o JObjectGpu.o 
OBSPHGPU=JArraysGpu.o JBlockSizeAuto.o JCellDivGpu.o JSphGpu.o 
OBSPHSINGLEGPU=JCellDivGpuSingle.o JSphGpuSingle.o
OBCUDA=JCellDivGpu_ker.o JCellDivGpuSingle_ker.o JGauge_ker.o JReduSum_ker.o JSphGpu_ker.o JWaveOrder2_ker.o

OBJECTS=$(OBJXML) $(OBJSPHMOTION) $(OBCOMMON) $(OBCOMMONDSPH) $(OBSPH) $(OBSPHSINGLE)
OBJECTS:=$(OBJECTS) $(OBCOMMONGPU) $(OBSPHGPU) $(OBSPHSINGLEGPU) $(OBCUDA)

#=============== To select GPU architecture ===============
#GENCODE:=$(GENCODE) -gencode=arch=compute_10,code=\"sm_10,compute_10\"
#GENCODE:=$(GENCODE) -gencode=arch=compute_12,code=\"sm_12,compute_12\"
#GENCODE:=$(GENCODE) -gencode=arch=compute_13,code=\"sm_13,compute_13\"
ifeq ($(CUDAVER),40)
  GENCODE:=$(GENCODE) -gencode=arch=compute_20,code=\"sm_20,compute_20\"
endif
ifeq ($(CUDAVER),75)
  GENCODE:=$(GENCODE) -gencode=arch=compute_20,code=\"sm_20,compute_20\"
  GENCODE:=$(GENCODE) -gencode=arch=compute_30,code=\"sm_30,compute_30\"
  GENCODE:=$(GENCODE) -gencode=arch=compute_35,code=\"sm_35,compute_35\"
endif
ifeq ($(CUDAVER),80)
  # GENCODE:=$(GENCODE) -gencode=arch=compute_20,code=\"sm_20,compute_20\"
  # GENCODE:=$(GENCODE) -gencode=arch=compute_30,code=\"sm_30,compute_30\"
  GENCODE:=$(GENCODE) -gencode=arch=compute_35,code=\"sm_35,compute_35\"
endif

#=============== DualSPHysics libs to include ===============
JLIBS=-L./ -ljformatfiles2_64 -ljwavegengpu_64

#=============== To compile CPU code ===============
CC=g++
CCFLAGS=-c -O3 -fopenmp -ffast-math -D_WITHGPU
CCLINKFLAGS=-fopenmp -lgomp
ifneq (,$(filter $(CLUSTER),BSC EME))
  CC=icpc
  CCFLAGS=-c -O3 -openmp -D_WITHGPU
  CCLINKFLAGS=-openmp
endif
CCFLAGS := $(CCFLAGS) $(CCMOREFLAGS)

#=============== To compile GPU code ===============
DIRTOOLKIT=
ifeq ($(CLUSTER),EPH)
  ifeq ($(CUDAVER),40)
    DIRTOOLKIT=/exports/opt/NVIDIA/cuda-4.0
  endif
  ifeq ($(CUDAVER),75)
    DIRTOOLKIT=/exports/opt/NVIDIA/cuda-7.5
  endif
  ifeq ($(CUDAVER),80)
    DIRTOOLKIT=/exports/opt/NVIDIA/cuda-8.0
  endif
endif
ifeq ($(CLUSTER),CSF)
  DIRTOOLKIT=/opt/gridware/libs/nvidia-cuda/toolkit/4.0.17
endif
ifeq ($(CLUSTER),BSC)
  DIRTOOLKIT=/opt/cuda/4.0
endif
ifeq ($(CLUSTER),EME)
  DIRTOOLKIT=/apps/cuda/4.0.17/cuda
endif
CCFLAGS := $(CCFLAGS) -I./ -I$(DIRTOOLKIT)/include
CCLINKFLAGS := $(CCLINKFLAGS) -L$(DIRTOOLKIT)/lib64 -lcudart
NCC=nvcc
NCCFLAGS=-c $(GENCODE) -use_fast_math -O3


all:dualsphysics4gpu 
	rm -rf *.o
	@echo "  --- Compiled GPU version for the cluster $(CLUSTER) ---"

dualsphysics4gpu:  $(OBJECTS)
	$(CC) $(CCLINKFLAGS) $(OBJECTS) -o $@ $(JLIBS)

.cpp.o: 
	$(CC) $(CCFLAGS) $< 

JSphGpu_ker.o: JSphGpu_ker.cu
	$(NCC) $(NCCFLAGS) JSphGpu_ker.cu

JCellDivGpu_ker.o: JCellDivGpu_ker.cu
	$(NCC) $(NCCFLAGS) JCellDivGpu_ker.cu

JCellDivGpuSingle_ker.o: JCellDivGpuSingle_ker.cu
	$(NCC) $(NCCFLAGS) JCellDivGpuSingle_ker.cu

JGauge_ker.o: JGauge_ker.cu
	$(NCC) $(NCCFLAGS) JGauge_ker.cu

JWaveOrder2_ker.o: JWaveOrder2_ker.cu
	$(NCC) $(NCCFLAGS) JWaveOrder2_ker.cu

JReduSum_ker.o: JReduSum_ker.cu
	$(NCC) $(NCCFLAGS) JReduSum_ker.cu

clean:
	rm -rf *.o dualsphysics4gpu

help:
	@echo "  make CLUSTER=[EPH,CSF,BSC,EME]"
	@echo "    EPH: Ephyslab - Universidad de Vigo"
	@echo "    CSF: University of Manchester"
	@echo "    BSC: Barcelona Supercomputing Center"
	@echo "    EME: Emerald - e-Infrastructure South"

