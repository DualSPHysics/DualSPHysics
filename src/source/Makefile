#DualSPHysics GPU/CPU v5.4.351 07-12-2024

#=============== Compilation Options (YES/NO) ===============
USE_DEBUG=NO
USE_FAST_MATH=YES
USE_NATIVE_CPU_OPTIMIZATIONS=NO
COMPILE_CHRONO=YES
COMPILE_WAVEGEN=YES
COMPILE_MOORDYNPLUS=YES

LIBS_DIRECTORIES=-L./
ifeq ($(GCC), 7)
  LIBS_DIRECTORIES:=$(LIBS_DIRECTORIES) -L../lib/linux_gcc7
else
  LIBS_DIRECTORIES:=$(LIBS_DIRECTORIES) -L../lib/linux_gcc
endif

EXECNAME=DualSPHysics5.4_linux64
EXECS_DIRECTORY=../../bin/linux

ifeq ($(GCC), 7)
  EXECS_DIRECTORY=../../bin/linux/DSGcc7
endif

# -std=c++0x ---> Used to avoid errors for calls to enums
ifeq ($(USE_DEBUG), YES)
  CCFLAGS=-c -O0 -g -Wall -fopenmp -D_WITHGPU -std=c++0x -D_WITHMR 
else
  CCFLAGS=-c -O3 -fopenmp -D_WITHGPU -std=c++0x -D_WITHMR
  ifeq ($(USE_FAST_MATH), YES)
    CCFLAGS+= -ffast-math
  endif
  ifeq ($(USE_NATIVE_CPU_OPTIMIZATIONS), YES)
    CCFLAGS+= -march=native
  endif
endif
CC=g++
CCLINKFLAGS=-fopenmp -lgomp

ifeq ($(COMPILE_CHRONO), NO)
  CCFLAGS:=$(CCFLAGS) -DDISABLE_CHRONO
endif
ifeq ($(COMPILE_WAVEGEN), NO)
  CCFLAGS:=$(CCFLAGS) -DDISABLE_WAVEGEN
endif
ifeq ($(COMPILE_MOORDYNPLUS), NO)
  CCFLAGS:=$(CCFLAGS) -DDISABLE_MOORDYNPLUS
endif

#=============== CUDA selection ===============
CUDAVER=12

ifneq ($(CUDA),)
  CUDAVER=$(CUDA)
endif

#=============== CUDA toolkit directory (make appropriate for local CUDA installation) ===============
ifeq ($(CUDAVER),00)
  DIRTOOLKIT=/usr/local/cuda
endif
ifeq ($(CUDAVER),75)
  DIRTOOLKIT=/exports/opt/NVIDIA/cuda-7.5
endif
ifeq ($(CUDAVER),91)
  DIRTOOLKIT=/exports/opt/NVIDIA/cuda-9.1
endif
ifeq ($(CUDAVER),92)
  DIRTOOLKIT=/exports/opt/NVIDIA/cuda-9.2
endif
ifeq ($(CUDAVER),11)
  DIRTOOLKIT=/exports/opt/NVIDIA/cuda-11.7
endif
ifeq ($(CUDAVER),12)
  DIRTOOLKIT=/exports/opt/NVIDIA/cuda-12.8
endif

#=============== Select GPU architectures ===============
ifeq ($(CUDAVER),00)
  GENCODE:=$(GENCODE) -gencode=arch=compute_20,code=\"sm_20,compute_20\"
  GENCODE:=$(GENCODE) -gencode=arch=compute_30,code=\"sm_30,compute_30\"
endif
ifeq ($(CUDAVER),75)
  GENCODE:=$(GENCODE) -gencode=arch=compute_20,code=\"sm_20,compute_20\"
  GENCODE:=$(GENCODE) -gencode=arch=compute_30,code=\"sm_30,compute_30\"
  GENCODE:=$(GENCODE) -gencode=arch=compute_35,code=\"sm_35,compute_35\"
  GENCODE:=$(GENCODE) -gencode=arch=compute_37,code=\"sm_37,compute_37\"
  GENCODE:=$(GENCODE) -gencode=arch=compute_50,code=\"sm_50,compute_50\"
  GENCODE:=$(GENCODE) -gencode=arch=compute_52,code=\"sm_52,compute_52\"
endif
ifeq ($(CUDAVER),91)
  GENCODE:=$(GENCODE) -gencode=arch=compute_30,code=\"sm_30,compute_30\"
  GENCODE:=$(GENCODE) -gencode=arch=compute_35,code=\"sm_35,compute_35\"
  GENCODE:=$(GENCODE) -gencode=arch=compute_50,code=\"sm_50,compute_50\"
  GENCODE:=$(GENCODE) -gencode=arch=compute_52,code=\"sm_52,compute_52\"
  GENCODE:=$(GENCODE) -gencode=arch=compute_61,code=\"sm_61,compute_61\"
  GENCODE:=$(GENCODE) -gencode=arch=compute_70,code=\"sm_70,compute_70\"
endif
ifeq ($(CUDAVER),92)
  # module load cuda/9.2
  GENCODE:=$(GENCODE) -gencode=arch=compute_30,code=\"sm_30,compute_30\"
  GENCODE:=$(GENCODE) -gencode=arch=compute_35,code=\"sm_35,compute_35\"
  GENCODE:=$(GENCODE) -gencode=arch=compute_50,code=\"sm_50,compute_50\"
  GENCODE:=$(GENCODE) -gencode=arch=compute_52,code=\"sm_52,compute_52\"
  GENCODE:=$(GENCODE) -gencode=arch=compute_61,code=\"sm_61,compute_61\"
  GENCODE:=$(GENCODE) -gencode=arch=compute_70,code=\"sm_70,compute_70\"
endif
ifeq ($(CUDAVER),11)
  # module load cuda/11.7
  # GENCODE:=$(GENCODE) -gencode=arch=compute_30,code=\"sm_30,compute_30\"
  GENCODE:=$(GENCODE) -gencode=arch=compute_35,code=\"sm_35,compute_35\"
  GENCODE:=$(GENCODE) -gencode=arch=compute_50,code=\"sm_50,compute_50\"
  GENCODE:=$(GENCODE) -gencode=arch=compute_52,code=\"sm_52,compute_52\"
  GENCODE:=$(GENCODE) -gencode=arch=compute_61,code=\"sm_61,compute_61\"
  GENCODE:=$(GENCODE) -gencode=arch=compute_70,code=\"sm_70,compute_70\"
  GENCODE:=$(GENCODE) -gencode=arch=compute_86,code=\"sm_86,compute_86\"
endif
ifeq ($(CUDAVER),12)
  # module load cuda/11.7
  # GENCODE:=$(GENCODE) -gencode=arch=compute_30,code=\"sm_30,compute_30\"
  # GENCODE:=$(GENCODE) -gencode=arch=compute_35,code=\"sm_35,compute_35\"
  # GENCODE:=$(GENCODE) -gencode=arch=compute_50,code=\"sm_50,compute_50\"
  # GENCODE:=$(GENCODE) -gencode=arch=compute_52,code=\"sm_52,compute_52\"
  GENCODE:=$(GENCODE) -gencode=arch=compute_61,code=\"sm_61,compute_61\"
  GENCODE:=$(GENCODE) -gencode=arch=compute_70,code=\"sm_70,compute_70\"
  GENCODE:=$(GENCODE) -gencode=arch=compute_86,code=\"sm_86,compute_86\"
endif

#=============== Files to compile ===============
OBJXML=JXml.o tinystr.o tinyxml.o tinyxmlerror.o tinyxmlparser.o
OBJSPHMOTION=JMotion.o JMotionData.o JMotionList.o JMotionMov.o JMotionObj.o JMotionPos.o JDsMotion.o
OBCOMMON=Functions.o FunGeo3d.o FunSphKernelsCfg.o JAppInfo.o JBinaryData.o JCfgRunBase.o JDataArrays.o
OBCOMMON:=$(OBCOMMON) JException.o JLinearValue.o JLog2.o JObject.o JOutputCsv.o JRadixSort.o
OBCOMMON:=$(OBCOMMON) JRangeFilter.o JReadDatafile.o JSaveCsv2.o JSpVtkData.o JSpVtkShape.o
OBCOMMON:=$(OBCOMMON) JTimeControl.o JTrianglesMesh.o randomc.o
OBCOMMONDSPH=JDsphConfig.o JDsPips.o JPartDataBi4.o JPartDataHead.o JDsPartFloatSave.o JDsPartMotionSave.o
OBCOMMONDSPH:=$(OBCOMMONDSPH) JPartFloatInfoBi4.o JPartMotRefBi4Save.o JPartOutBi4Save.o JCaseCtes.o
OBCOMMONDSPH:=$(OBCOMMONDSPH) JCaseEParms.o JCaseParts.o JCaseProperties.o
OBSPH=JArraysCpu.o JCellDivCpu.o JDebugSphCpu.o JSphCfgRun.o JComputeMotionRef.o JDsDcell.o JDsDamping.o
OBSPH:=$(OBSPH) JDsExtraData.o JDsGaugeItem.o JDsGaugeSystem.o JDsPartsOut.o JDsSaveDt.o JSphShifting.o
OBSPH:=$(OBSPH) JSphShiftingAdv.o JSph.o JDsAccInput.o JSphCpu.o JSphCpu_mdbc.o JSphCpu_preloop.o
OBSPH:=$(OBSPH) JDsInitialize.o JSphMk.o JDsOutputParts.o JDsPartsInit.o JDsFixedDt.o JDsViscoInput.o
OBSPH:=$(OBSPH) JDsOutputTime.o JDsTimers.o JWaveSpectrumGpu.o main.o
OBSPHSINGLE=JCellDivCpuSingle.o JPartsLoad4.o JSphCpuSingle.o
OBCOMMONGPU=FunctionsCuda.o JObjectGpu.o 
OBSPHGPU=JArraysGpu.o JDebugSphGpu.o JCellDivGpu.o JSphGpu.o JDsGpuInfo.o 
OBSPHSINGLEGPU=JCellDivGpuSingle.o JSphGpuSingle.o
OBCUDA=JCellDivGpu_ker.o JCellDivGpuSingle_ker.o JDsPips_ker.o JDsGauge_ker.o JReduSum_ker.o JSphShifting_ker.o
OBCUDA:=$(OBCUDA) JDsAccInput_ker.o JSphGpu_ker.o JSphGpuSimple_ker.o JWaveOrder2_ker.o 
OBWAVERZ=JMLPistonsGpu.o JRelaxZonesGpu.o
OBWAVERZCUDA=JRelaxZone_ker.o
OBCHRONO=JChronoGeoObj.o JChronoObjects.o
OBMOORDYNPLUS=JDsMooredFloatings.o JDsFtForcePoints.o
OBINOUT=JSphCpu_InOut.o JSphCpuSingle_InOut.o JSphInOut.o JSphInOutZone.o JSphInOutGridData.o
OBINOUT:=$(OBINOUT) JSphInOutPoints.o JSphInOutVel.o JSphInOutVelAwas.o JSphInOutZsurf.o JSimpleNeigs.o
OBINOUTGPU=JSphGpuSingle_InOut.o
OBMESH=JMeshData.o JMeshTDatas.o JMeshTDatasLoad.o JMeshTDatasSave.o JMeshTDatasDsVel.o JMeshTDatasXml.o
OBVRES=JSphVRes.o JSphVResZone.o JSphCpuSingle_VRes.o JSphGpuSingle_VRes.o JSphCpu_VRes.o JCaseVRes.o JVResDataLimits.o JDsVresData.o JBoxDef.o
OBFLEXSTRUC=JSphFlexStruc.o

OBJECTS=$(OBJXML) $(OBJSPHMOTION) $(OBCOMMON) $(OBCOMMONDSPH) $(OBSPH) $(OBSPHSINGLE)
OBJECTS:=$(OBJECTS) $(OBCOMMONGPU) $(OBSPHGPU) $(OBSPHSINGLEGPU) $(OBCUDA)
OBJECTS:=$(OBJECTS) $(OBWAVERZ) $(OBWAVERZCUDA) $(OBCHRONO) $(OBMOORDYNPLUS)
OBJECTS:=$(OBJECTS) $(OBINOUT) $(OBINOUTGPU) $(OBMESH) $(OBVRES) $(OBFLEXSTRUC)

#=============== DualSPHysics libs to be included ===============
JLIBS=${LIBS_DIRECTORIES}
ifeq ($(COMPILE_CHRONO), YES)
  JLIBS:=$(JLIBS) -ldsphchrono -lChronoEngine 
endif
ifeq ($(COMPILE_WAVEGEN), YES)
  JLIBS:=$(JLIBS) -ljwavegen_64
endif
ifeq ($(COMPILE_MOORDYNPLUS), YES)
  JLIBS:=$(JLIBS) -ldsphmoordynplus_64
endif

#=============== GPU Code Compilation ===============
CCFLAGS := $(CCFLAGS) -I./ -I$(DIRTOOLKIT)/include
CCLINKFLAGS := $(CCLINKFLAGS) -L$(DIRTOOLKIT)/lib64 -lcudart_static -ldl -lrt
NCC=nvcc
ifeq ($(USE_DEBUG), NO)
  NCCFLAGS+=-c $(GENCODE) -O3 -ccbin $(CC) -D_WITHMR -D_WITHGPU -DCUB_IGNORE_DEPRECATED_CPP_DIALECT -DTHRUST_IGNORE_DEPRECATED_CPP_DIALECT
else
  NCCFLAGS+=-c $(GENCODE) -g -G -ccbin $(CC)  -D_WITHMR -D_WITHGPU -DCUB_IGNORE_DEPRECATED_CPP_DIALECT -DTHRUST_IGNORE_DEPRECATED_CPP_DIALECT
endif
ifeq ($(USE_FAST_MATH), YES)
  NCCFLAGS+= -use_fast_math
endif
ifeq ($(CUDAVER),11)	
  NCCFLAGS+= -Wno-deprecated-gpu-targets
endif

all:compile_libs $(EXECS_DIRECTORY)/$(EXECNAME)
	rm -rf *.o
ifeq ($(USE_DEBUG), NO)
	@echo "  --- Compiled Release GPU version ---"
else
	@echo "  --- Compiled Debug GPU version ---"
	mv $(EXECS_DIRECTORY)/$(EXECNAME) $(EXECNAME)_debug
endif

compile_libs:
ifeq ($(COMPILE_MOORDYNPLUS), YES)
	$(MAKE) -C Source_DSphMoorDynPlus -f Makefile_MoorDynPlus $(MAKEFLAGS)
endif

$(EXECS_DIRECTORY)/$(EXECNAME):  $(OBJECTS)
	$(CC) $(OBJECTS) $(CCLINKFLAGS) -o $@ $(JLIBS)

.cpp.o: 
	$(CC) $(CCFLAGS) $< 

JSphGpu_ker.o: JSphGpu_ker.cu
	$(NCC) $(NCCFLAGS) JSphGpu_ker.cu

JSphGpuSimple_ker.o: JSphGpuSimple_ker.cu
	$(NCC) $(NCCFLAGS) JSphGpuSimple_ker.cu

JCellDivGpu_ker.o: JCellDivGpu_ker.cu
	$(NCC) $(NCCFLAGS) JCellDivGpu_ker.cu

JCellDivGpuSingle_ker.o: JCellDivGpuSingle_ker.cu
	$(NCC) $(NCCFLAGS) JCellDivGpuSingle_ker.cu

JSphShifting_ker.o: JSphShifting_ker.cu
	$(NCC) $(NCCFLAGS) JSphShifting_ker.cu

JDsAccInput_ker.o: JDsAccInput_ker.cu
	$(NCC) $(NCCFLAGS) JDsAccInput_ker.cu

JDsPips_ker.o: JDsPips_ker.cu
	$(NCC) $(NCCFLAGS) JDsPips_ker.cu

JDsGauge_ker.o: JDsGauge_ker.cu
	$(NCC) $(NCCFLAGS) JDsGauge_ker.cu

JWaveOrder2_ker.o: JWaveOrder2_ker.cu
	$(NCC) $(NCCFLAGS) JWaveOrder2_ker.cu

JReduSum_ker.o: JReduSum_ker.cu
	$(NCC) $(NCCFLAGS) JReduSum_ker.cu

JRelaxZone_ker.o: JRelaxZone_ker.cu
	$(NCC) $(NCCFLAGS) JRelaxZone_ker.cu

clean:
ifeq ($(COMPILE_MOORDYNPLUS), YES)
	$(MAKE) -C Source_DSphMoorDynPlus -f Makefile_MoorDynPlus clean $(MAKEFLAGS)
endif
	rm -rf *.o $(EXECS_DIRECTORY)/$(EXECNAME) $(EXECS_DIRECTORY)/$(EXECNAME)_debug
