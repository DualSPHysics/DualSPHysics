#==================================================================
# DualSPHysics GPU/CPU v5.0.164 21-11-2020
#==================================================================

cmake_minimum_required(VERSION 2.8)

PROJECT(DualSPHysics)

#------------------------------------------------------------------
# Coupling options
#------------------------------------------------------------------
option(ENABLE_MOORDYN "Enable the MoorDyn+ library" ON)
option(ENABLE_CHRONO "Enable the Chrono Engine library" ON)

#------------------------------------------------------------------
# Source files
#------------------------------------------------------------------

# CPU Objects
set(OBJXML JXml.cpp tinystr.cpp tinyxml.cpp tinyxmlerror.cpp tinyxmlparser.cpp)
set(OBJSPHMOTION JMotion.cpp JMotionList.cpp JMotionMov.cpp JMotionObj.cpp JMotionPos.cpp JDsMotion.cpp)
set(OBCOMMON Functions.cpp FunGeo3d.cpp FunSphKernelsCfg.cpp JAppInfo.cpp JBinaryData.cpp JCfgRunBase.cpp JDataArrays.cpp JException.cpp JLinearValue.cpp JLog2.cpp JObject.cpp JOutputCsv.cpp JRadixSort.cpp JRangeFilter.cpp JReadDatafile.cpp JSaveCsv2.cpp JTimeControl.cpp randomc.cpp)
set(OBCOMMONDSPH JDsphConfig.cpp JDsPips.cpp JPartDataBi4.cpp JPartDataHead.cpp JPartFloatBi4.cpp JPartOutBi4Save.cpp JCaseCtes.cpp JCaseEParms.cpp JCaseParts.cpp JCaseProperties.cpp JCaseUserVars.cpp JCaseVtkOut.cpp)
set(OBSPH JArraysCpu.cpp JCellDivCpu.cpp JSphCfgRun.cpp JComputeMotionRef.cpp JDsDamping.cpp JDsGaugeItem.cpp JDsGaugeSystem.cpp JDsPartsOut.cpp JDsSaveDt.cpp JSphShifting.cpp JSph.cpp JDsAccInput.cpp JSphCpu.cpp JDsInitialize.cpp JSphMk.cpp JDsPartsInit.cpp JDsFixedDt.cpp JDsViscoInput.cpp JDsOutputTime.cpp JFtMotionSave.cpp JWaveAwasZsurf.cpp JWaveSpectrumGpu.cpp main.cpp)
set(OBSPHSINGLE JCellDivCpuSingle.cpp JPartsLoad4.cpp JSphCpuSingle.cpp)

# GPU Objects
set(OBCOMMONGPU FunctionsCuda.cpp JObjectGpu.cpp)
set(OBSPHGPU JArraysGpu.cpp JDebugSphGpu.cpp JCellDivGpu.cpp JSphGpu.cpp)
set(OBSPHSINGLEGPU JCellDivGpuSingle.cpp JSphGpuSingle.cpp)
set(OBCUDA JCellDivGpu_ker.cu JCellDivGpuSingle_ker.cu JDsPips_ker.cu JDsGauge_ker.cu JReduSum_ker.cu JSphShifting_ker.cu JDsAccInput_ker.cu JSphGpu_ker.cu JSphGpuSimple_ker.cu JWaveOrder2_ker.cu)

# Other Objects
set(OBWAVERZ JMLPistonsGpu.cpp JRelaxZonesGpu.cpp)
set(OBWAVERZCUDA JRelaxZone_ker.cu)
set(OBCHRONO JChronoObjects.cpp)
set(OBMOORDYN JDsMooredFloatings.cpp JDsFtForcePoints.cpp)
set(OBINOUT JSphCpu_InOut.cpp JSphCpuSingle_InOut.cpp JSphBoundCorr.cpp JSphInOut.cpp JSphInOutZone.cpp JSphInOutGridData.cpp JSphInOutPoints.cpp JSphInOutVel.cpp JSphInOutVelAwas.cpp JSphInOutZsurf.cpp JSimpleNeigs.cpp)
set(OBINOUTGPU JSphGpuSingle_InOut.cpp)
set(OBMDBC JPartNormalData.cpp JNormalsMarrone.cpp)
set(OBNNMULTIPHSE JSphCpu_Tensors.cpp)

#------------------------------------------------------------------
# NVCC Flags
#------------------------------------------------------------------
find_package(CUDA QUIET)

if(CUDA_FOUND)
  if(CUDA_VERSION VERSION_LESS "7.5")
    message("Using cuda version <7.5")
    list(APPEND CUDA_NVCC_FLAGS "-use_fast_math -O3 -gencode=arch=compute_20,code=\"sm_20,compute_20\"")
  elseif(CUDA_VERSION VERSION_GREATER "7.4" AND CUDA_VERSION VERSION_LESS "9.1")
    message("Using cuda version >=7.5 and <9.1")
    list(APPEND CUDA_NVCC_FLAGS "-use_fast_math -O3  -gencode=arch=compute_20,code=\"sm_20,compute_20\" -gencode=arch=compute_30,code=\"sm_30,compute_30\" -gencode=arch=compute_35,code=\"sm_35,compute_35\" -gencode=arch=compute_37,code=\"sm_37,compute_37\" -gencode=arch=compute_50,code=\"sm_50,compute_50\" -gencode=arch=compute_52,code=\"sm_52,compute_52\"")
  else()
    message("Using cuda version >=9.1")
    list(APPEND CUDA_NVCC_FLAGS "-use_fast_math -O3  -gencode=arch=compute_30,code=\"sm_30,compute_30\" -gencode=arch=compute_35,code=\"sm_35,compute_35\" -gencode=arch=compute_37,code=\"sm_37,compute_37\" -gencode=arch=compute_50,code=\"sm_50,compute_50\" -gencode=arch=compute_52,code=\"sm_52,compute_52\" -gencode=arch=compute_61,code=\"sm_61,compute_61\" -gencode=arch=compute_70,code=\"sm_70,compute_70\"")
  endif()
else()
  message("CUDA Libraries were not found.")
endif()

#------------------------------------------------------------------
# Static libraries linker path
#------------------------------------------------------------------
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set(PLATFORM "linux")
  message("Using libraries for gcc")
  LINK_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../lib/linux_gcc)
elseif(MSVC)
  message("Windows version")
  set(PLATFORM "windows")
  #if(MSVC_VERSION VERSION_EQUAL 1900)
    LINK_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../lib/vs2015)
  #endif()
endif()

#--------------------------------------------
# OpenMP flags
#---------------------------------------------
find_package(OpenMP)

if (OPENMP_FOUND)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

#------------------------------------------------------------------
# Don't use C++11 ABI for compatibility between GCC v4.9 and >=5
#------------------------------------------------------------------
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=0 -fPIC")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=0 -fPIC")

#------------------------------------------------------------------
# Binaries
#------------------------------------------------------------------
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../../../bin/${PLATFORM}/DSNNewtonian)# Generic output directory
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})         # Release output directory
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})           # Debug output directory

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  add_executable(DualSPHysics5.0_NNewtonianCPU_linux64 ${OBJXML} ${OBJSPHMOTION} ${OBCOMMON} ${OBCOMMONDSPH} ${OBSPH} ${OBSPHSINGLE} ${OBWAVERZ} ${OBCHRONO} ${OBMOORDYN} ${OBINOUT} ${OBMDBC} ${OBNNMULTIPHSE})
  install(TARGETS	DualSPHysics5.0_NNewtonianCPU_linux64 DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
  if (CUDA_FOUND)
    cuda_add_executable(DualSPHysics5.0_NNewtonian_linux64 ${OBJXML} ${OBJSPHMOTION} ${OBCOMMON} ${OBCOMMONDSPH} ${OBSPH} ${OBSPHSINGLE} ${OBCOMMONGPU} ${OBSPHGPU} ${OBSPHSINGLEGPU} ${OBCUDA} ${OBWAVERZ} ${OBWAVERZCUDA} ${OBCHRONO} ${OBMOORDYN} ${OBINOUT} ${OBINOUTGPU} ${OBMDBC} ${OBNNMULTIPHSE})
    install(TARGETS DualSPHysics5.0_NNewtonian_linux64 DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
  endif(CUDA_FOUND)
elseif(MSVC) 
  add_executable(DualSPHysics5.0_NNewtonianCPU_win64 ${OBJXML} ${OBJSPHMOTION} ${OBCOMMON} ${OBCOMMONDSPH} ${OBSPH} ${OBSPHSINGLE} ${OBWAVERZ} ${OBCHRONO} ${OBMOORDYN} ${OBINOUT} ${OBMDBC} ${OBNNMULTIPHSE})
  install(TARGETS	DualSPHysics5.0_NNewtonianCPU_win64 DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
  if (CUDA_FOUND)
    cuda_add_executable(DualSPHysics5.0_NNewtonian_win64 ${OBJXML} ${OBJSPHMOTION} ${OBCOMMON} ${OBCOMMONDSPH} ${OBSPH} ${OBSPHSINGLE} ${OBCOMMONGPU} ${OBSPHGPU} ${OBSPHSINGLEGPU} ${OBCUDA} ${OBWAVERZ} ${OBWAVERZCUDA} ${OBCHRONO} ${OBMOORDYN} ${OBINOUT} ${OBINOUTGPU} ${OBMDBC} ${OBNNMULTIPHSE})
    install(TARGETS DualSPHysics5.0_NNewtonian_win64 DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
  endif(CUDA_FOUND)
endif()

foreach(OUTPUTCONFIG ${CMAKE_CONFIGURATION_TYPES} )
  string(TOUPPER ${OUTPUTCONFIG} OUTPUTCONFIG )
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} )
endforeach(OUTPUTCONFIG CMAKE_CONFIGURATION_TYPES )

#------------------------------------------------------------------
# General linker flags
#------------------------------------------------------------------
set(LINKER_FLAGS "")
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set(LINKER_FLAGS ${LINKER_FLAGS} jvtklib_64 jwavegen_64 jnumexlib_64)
elseif(MSVC) 
  set(LINKER_FLAGS ${LINKER_FLAGS} LibJVtkLib_x64_v140_Release LibJWaveGen_x64_v140_Release LibJNumexLib_x64_v140_Release)
endif()

#------------------------------------------------------------------
# Coupled libraries
#------------------------------------------------------------------
# MoorDyn+
if(ENABLE_MOORDYN)
  message(STATUS "Using MoorDyn+")
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(LINKER_FLAGS ${LINKER_FLAGS} dsphmoordyn_64)
  elseif(MSVC) 
    set(LINKER_FLAGS ${LINKER_FLAGS} LibDSphMoorDyn_x64_v140_Release)
  endif()
else()
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DDISABLE_MOORDYN")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDISABLE_MOORDYN")
endif(ENABLE_MOORDYN)

# Chrono Engine
if(ENABLE_CHRONO)
  message(STATUS "Using ChronoEngine Module")
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(LINKER_FLAGS ${LINKER_FLAGS} ChronoEngine dsphchrono)
  elseif(MSVC) 
    set(LINKER_FLAGS ${LINKER_FLAGS} dsphchrono)
  endif()
else()
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDISABLE_CHRONO")
endif(ENABLE_CHRONO)

#------------------------------------------------------------------
# Linker flags
#------------------------------------------------------------------
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  target_link_libraries(DualSPHysics5.0_NNewtonianCPU_linux64 ${LINKER_FLAGS})
  set_target_properties(DualSPHysics5.0_NNewtonianCPU_linux64 PROPERTIES COMPILE_FLAGS "-use_fast_math -O3 -D_GLIBCXX_USE_CXX11_ABI=0 -fPIC -std=c++0x")
  
  if (CUDA_FOUND)
    target_link_libraries(DualSPHysics5.0_NNewtonian_linux64 ${LINKER_FLAGS})
    set_target_properties(DualSPHysics5.0_NNewtonian_linux64 PROPERTIES COMPILE_FLAGS "-use_fast_math -O3 -D_WITHGPU -D_GLIBCXX_USE_CXX11_ABI=0 -fPIC -std=c++0x")
  endif(CUDA_FOUND)
  
elseif(MSVC)
  set_target_properties(DualSPHysics5.0_NNewtonian_win64 PROPERTIES COMPILE_FLAGS "/D _WITHGPU")

  #if(MSVC_VERSION VERSION_EQUAL 1900)
    # MSVC 2015
    target_link_libraries(DualSPHysics5.0_NNewtonianCPU_win64  ${LINKER_FLAGS})
    target_link_libraries(DualSPHysics5.0_NNewtonian_win64 ${LINKER_FLAGS})
  #endif()
  
  SET(CUDA_PROPAGATE_HOST_FLAGS OFF CACHE BOOL "Propagate C/CXX Flags and friends to the host compiler in NVCC via -Xompile  " FORCE)
  
  foreach(CPP_FLAGS CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO)
    if(${CPP_FLAGS} MATCHES "/MD")
      string(REGEX REPLACE "/MD" "/MT" ${CPP_FLAGS} "${${CPP_FLAGS}}")
    endif(${CPP_FLAGS} MATCHES "/MD")
  endforeach(CPP_FLAGS)
endif()
