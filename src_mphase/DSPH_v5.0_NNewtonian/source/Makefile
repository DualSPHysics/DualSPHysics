#DualSPHysics NNewtonian GPU/CPU v5.0.164 21-11-2020

#=============== Compilation Options (YES/NO) ===============
USE_GCC5=NO
USE_DEBUG=NO
USE_FAST_MATH=YES
USE_NATIVE_CPU_OPTIMIZATIONS=NO
COMPILE_VTKLIB=YES
COMPILE_NUMEXLIB=YES
COMPILE_CHRONO=YES
COMPILE_WAVEGEN=YES
COMPILE_MOORDYN=YES

LIBS_DIRECTORIES=-L./
LIBS_DIRECTORIES:=$(LIBS_DIRECTORIES) -L../lib/linux_gcc

EXECNAME=DualSPHysics5.0_NNewtonian_linux64
EXECS_DIRECTORY=../../../bin/linux/DSNNewtonian

# -std=c++0x ---> Used to avoid errors for calls to enums
ifeq ($(USE_DEBUG), YES)
  CCFLAGS=-c -O0 -g -Wall -fopenmp -D_WITHGPU -std=c++0x
else
  CCFLAGS=-c -O3 -fopenmp -D_WITHGPU -std=c++0x
  ifeq ($(USE_FAST_MATH), YES)
    CCFLAGS+= -ffast-math
  endif
  ifeq ($(USE_NATIVE_CPU_OPTIMIZATIONS), YES)
    CCFLAGS+= -march=native
  endif
endif
CC=g++
CCLINKFLAGS=-fopenmp -lgomp

#Required for GCC versions >=5.0
ifeq ($(USE_GCC5), YES)
  CCFLAGS+=-D_GLIBCXX_USE_CXX11_ABI=0
  CCLINKFLAGS+=-D_GLIBCXX_USE_CXX11_ABI=0
  NCCFLAGS+=-D_GLIBCXX_USE_CXX11_ABI=0
endif

ifeq ($(COMPILE_VTKLIB), NO)
  CCFLAGS:=$(CCFLAGS) -DDISABLE_VTKLIB
endif
ifeq ($(COMPILE_NUMEXLIB), NO)
  CCMOREFLAGS:=$(CCMOREFLAGS) -DDISABLE_NUMEXLIB
endif
ifeq ($(COMPILE_CHRONO), NO)
  CCFLAGS:=$(CCFLAGS) -DDISABLE_CHRONO
endif
ifeq ($(COMPILE_WAVEGEN), NO)
  CCFLAGS:=$(CCFLAGS) -DDISABLE_WAVEGEN
endif
ifeq ($(COMPILE_MOORDYN), NO)
  CCFLAGS:=$(CCFLAGS) -DDISABLE_MOORDYN
endif

#=============== CUDA selection ===============
CUDAVER=92

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


#=============== Files to compile ===============
OBJXML=JXml.o tinystr.o tinyxml.o tinyxmlerror.o tinyxmlparser.o
OBJSPHMOTION=JMotion.o JMotionList.o JMotionMov.o JMotionObj.o JMotionPos.o JDsMotion.o
OBCOMMON=Functions.o FunGeo3d.o FunSphKernelsCfg.o JAppInfo.o JBinaryData.o JCfgRunBase.o JDataArrays.o JException.o JLinearValue.o JLog2.o JObject.o JOutputCsv.o JRadixSort.o JRangeFilter.o JReadDatafile.o JSaveCsv2.o JTimeControl.o randomc.o
OBCOMMONDSPH=JDsphConfig.o JDsPips.o JPartDataBi4.o JPartDataHead.o JPartFloatBi4.o JPartOutBi4Save.o JCaseCtes.o JCaseEParms.o JCaseParts.o JCaseProperties.o JCaseUserVars.o JCaseVtkOut.o
OBSPH=JArraysCpu.o JCellDivCpu.o JSphCfgRun.o JComputeMotionRef.o JDsDamping.o JDsGaugeItem.o JDsGaugeSystem.o JDsPartsOut.o JDsSaveDt.o JSphShifting.o JSph.o JDsAccInput.o JSphCpu.o JDsInitialize.o JFtMotionSave.o JSphMk.o JDsPartsInit.o JDsFixedDt.o JDsViscoInput.o JDsOutputTime.o JWaveAwasZsurf.o JWaveSpectrumGpu.o main.o
OBSPHSINGLE=JCellDivCpuSingle.o JPartsLoad4.o JSphCpuSingle.o
OBCOMMONGPU=FunctionsCuda.o JObjectGpu.o 
OBSPHGPU=JArraysGpu.o JDebugSphGpu.o JCellDivGpu.o JSphGpu.o 
OBSPHSINGLEGPU=JCellDivGpuSingle.o JSphGpuSingle.o
OBCUDA=JCellDivGpu_ker.o JCellDivGpuSingle_ker.o JDsPips_ker.o JDsGauge_ker.o JReduSum_ker.o JSphShifting_ker.o JDsAccInput_ker.o JSphGpu_ker.o JSphGpuSimple_ker.o JWaveOrder2_ker.o

OBWAVERZ=JMLPistonsGpu.o JRelaxZonesGpu.o
OBWAVERZCUDA=JRelaxZone_ker.o
OBCHRONO=JChronoObjects.o
OBMOORDYN=JDsMooredFloatings.o JDsFtForcePoints.o
OBINOUT=JSphCpu_InOut.o JSphCpuSingle_InOut.o JSphBoundCorr.o JSphInOut.o JSphInOutZone.o JSphInOutGridData.o JSphInOutPoints.o JSphInOutVel.o JSphInOutVelAwas.o JSphInOutZsurf.o JSimpleNeigs.o
OBINOUTGPU=JSphGpuSingle_InOut.o
OBMDBC=JPartNormalData.o JNormalsMarrone.o
OBNN=JSphCpu_Tensors.o

OBJECTS=$(OBJXML) $(OBJSPHMOTION) $(OBCOMMON) $(OBCOMMONDSPH) $(OBSPH) $(OBSPHSINGLE)
OBJECTS:=$(OBJECTS) $(OBCOMMONGPU) $(OBSPHGPU) $(OBSPHSINGLEGPU) $(OBCUDA)
OBJECTS:=$(OBJECTS) $(OBWAVERZ) $(OBWAVERZCUDA) $(OBCHRONO) $(OBMOORDYN) $(OBINOUT) $(OBINOUTGPU) $(OBMDBC)
OBJECTS:=$(OBJECTS) $(OBNN)

#=============== DualSPHysics libs to be included ===============
JLIBS=${LIBS_DIRECTORIES}
ifeq ($(COMPILE_VTKLIB), YES)
  JLIBS:=$(JLIBS) -ljvtklib_64
endif
ifeq ($(COMPILE_NUMEXLIB), YES)
  JLIBS:=$(JLIBS) -ljnumexlib_64
endif
ifeq ($(COMPILE_CHRONO), YES)
  JLIBS:=$(JLIBS) -ldsphchrono -lChronoEngine 
endif
ifeq ($(COMPILE_WAVEGEN), YES)
  JLIBS:=$(JLIBS) -ljwavegen_64
endif
ifeq ($(COMPILE_MOORDYN), YES)
  JLIBS:=$(JLIBS) -ldsphmoordyn_64
endif

#=============== GPU Code Compilation ===============
CCFLAGS := $(CCFLAGS) -I./ -I$(DIRTOOLKIT)/include
CCLINKFLAGS := $(CCLINKFLAGS) -L$(DIRTOOLKIT)/lib64 -lcudart_static -ldl -lrt
NCC=nvcc
ifeq ($(USE_DEBUG), NO)
  NCCFLAGS+=-c $(GENCODE) -O3 -ccbin $(CC)
else
  NCCFLAGS+=-c $(GENCODE) -O0 -ccbin $(CC) -g
endif
ifeq ($(USE_FAST_MATH), YES)
  NCCFLAGS+= -use_fast_math
endif

all:$(EXECS_DIRECTORY)/$(EXECNAME)
	rm -rf *.o
ifeq ($(USE_DEBUG), NO)
	@echo "  --- Compiled Release GPU version ---"
else
	@echo "  --- Compiled Debug GPU version ---"
	mv $(EXECS_DIRECTORY)/$(EXECNAME) $(EXECNAME)_debug
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
	rm -rf *.o $(EXECNAME) $(EXECNAME)_debug

