#==================================================================
# DualSPHysics GPU/CPU v5.4.316 21-02-2024
#==================================================================

cmake_minimum_required(VERSION 3.0)

PROJECT(DualSPHysics)

#------------------------------------------------------------------
# Coupling options
#------------------------------------------------------------------
option(ENABLE_CHRONO "Enable the Chrono Engine library" ON)
option(ENABLE_WAVEGEN "Enable the WaveGen library" ON)
option(ENABLE_MOORDYNPLUS "Enable the MoorDynPlus library" ON)
option(ENABLE_CUDA "Enable the use of CUDA for compilation" ON)

#------------------------------------------------------------------
# Source files
#------------------------------------------------------------------

# CPU Objects
set(OBJXML JXml.cpp tinystr.cpp tinyxml.cpp tinyxmlerror.cpp tinyxmlparser.cpp)
set(OBJSPHMOTION JMotion.cpp JMotionData.cpp JMotionList.cpp JMotionMov.cpp JMotionObj.cpp JMotionPos.cpp JDsMotion.cpp)
set(OBCOMMON Functions.cpp FunGeo3d.cpp FunSphKernelsCfg.cpp JAppInfo.cpp JBinaryData.cpp JCfgRunBase.cpp JDataArrays.cpp JException.cpp JLinearValue.cpp JLog2.cpp JObject.cpp JOutputCsv.cpp JRadixSort.cpp JRangeFilter.cpp JReadDatafile.cpp JSaveCsv2.cpp JSpVtkData.cpp JSpVtkShape.cpp JTimeControl.cpp JTrianglesMesh.cpp randomc.cpp)
set(OBCOMMONDSPH JDsphConfig.cpp JDsPips.cpp JPartDataBi4.cpp JPartDataHead.cpp JDsPartFloatSave.cpp JDsPartMotionSave.cpp JPartFloatInfoBi4.cpp JPartMotRefBi4Save.cpp JPartOutBi4Save.cpp JCaseCtes.cpp JCaseEParms.cpp JCaseParts.cpp JCaseProperties.cpp)
set(OBSPH JArraysCpu.cpp JCellDivCpu.cpp JDebugSphCpu.cpp JSphCfgRun.cpp JComputeMotionRef.cpp JDsDcell.cpp JDsDamping.cpp JDsExtraData.cpp JDsGaugeItem.cpp JDsGaugeSystem.cpp JDsPartsOut.cpp JDsSaveDt.cpp JSphShifting.cpp JSphShiftingAdv.cpp JSph.cpp JDsAccInput.cpp JSphCpu.cpp JSphCpu_mdbc.cpp JSphCpu_preloop.cpp JDsInitialize.cpp JSphMk.cpp JDsOutputParts.cpp JDsPartsInit.cpp JDsFixedDt.cpp JDsViscoInput.cpp JDsOutputTime.cpp JDsTimers.cpp JWaveSpectrumGpu.cpp main.cpp)
set(OBSPHSINGLE JCellDivCpuSingle.cpp JPartsLoad4.cpp JSphCpuSingle.cpp)

# GPU Objects
set(OBCOMMONGPU FunctionsCuda.cpp JObjectGpu.cpp)
set(OBSPHGPU JArraysGpu.cpp JDebugSphGpu.cpp JCellDivGpu.cpp JSphGpu.cpp JDsGpuInfo.cpp)
set(OBSPHSINGLEGPU JCellDivGpuSingle.cpp JSphGpuSingle.cpp)
set(OBCUDA JCellDivGpu_ker.cu JCellDivGpuSingle_ker.cu JDsPips_ker.cu JDsGauge_ker.cu JReduSum_ker.cu JSphShifting_ker.cu JDsAccInput_ker.cu JSphGpu_ker.cu JSphGpuSimple_ker.cu JWaveOrder2_ker.cu)
set(OBVRESGPU JSphGpuSingle_VRes.cpp)
set(OBWAVERZGPU JMLPistonsGpu.cpp JRelaxZonesGpu.cpp)
set(OBINOUTGPU JSphGpuSingle_InOut.cpp)

# Other Objects
set(OBWAVERZCUDA JRelaxZone_ker.cu)
set(OBCHRONO JChronoGeoObj.cpp JChronoObjects.cpp)
set(OBMOORDYNPLUS JDsMooredFloatings.cpp JDsFtForcePoints.cpp)
set(OBINOUT JSphCpu_InOut.cpp JSphCpuSingle_InOut.cpp JSphInOut.cpp JSphInOutZone.cpp JSphInOutGridData.cpp JSphInOutPoints.cpp JSphInOutVel.cpp JSphInOutVelAwas.cpp JSphInOutZsurf.cpp JSimpleNeigs.cpp)
set(OBMESH JMeshData.cpp JMeshTDatas.cpp JMeshTDatasLoad.cpp JMeshTDatasSave.cpp JMeshTDatasDsVel.cpp JMeshTDatasXml.cpp)
set(OBVRES JSphVRes.cpp JSphVResZone.cpp JSphCpuSingle_VRes.cpp JSphCpu_VRes.cpp JCaseVRes.cpp JVResDataLimits.cpp JDsVresData.cpp JBoxDef.cpp)
set(OBFLEXSTRUC JSphFlexStruc.cpp)

#------------------------------------------------------------------
# NVCC Flags
#------------------------------------------------------------------
if(ENABLE_CUDA)
  find_package(CUDA QUIET)
  if(CUDA_FOUND)
    if(CUDA_VERSION VERSION_LESS "7.5")
      message("Using CUDA version <7.5")
      list(APPEND CUDA_NVCC_FLAGS "-use_fast_math -O3 -D_WITHMR -D_WITHGPU -gencode=arch=compute_20,code=\"sm_20,compute_20\"")
    elseif(CUDA_VERSION VERSION_GREATER "7.4" AND CUDA_VERSION VERSION_LESS "9.1")
      message("Using CUDA version >=7.5 and <9.1")
      list(APPEND CUDA_NVCC_FLAGS "-use_fast_math -O3 -D_WITHMR -D_WITHGPU -gencode=arch=compute_20,code=\"sm_20,compute_20\" -gencode=arch=compute_30,code=\"sm_30,compute_30\" -gencode=arch=compute_35,code=\"sm_35,compute_35\" -gencode=arch=compute_37,code=\"sm_37,compute_37\" -gencode=arch=compute_50,code=\"sm_50,compute_50\" -gencode=arch=compute_52,code=\"sm_52,compute_52\"")
    elseif(CUDA_VERSION VERSION_GREATER "9.1" AND CUDA_VERSION VERSION_LESS "11.0")
      message("Using CUDA version >=9.1 and <11.0")
      list(APPEND CUDA_NVCC_FLAGS "-use_fast_math -O3 -D_WITHMR -D_WITHGPU -gencode=arch=compute_30,code=\"sm_30,compute_30\" -gencode=arch=compute_35,code=\"sm_35,compute_35\" -gencode=arch=compute_50,code=\"sm_50,compute_50\" -gencode=arch=compute_52,code=\"sm_52,compute_52\" -gencode=arch=compute_61,code=\"sm_61,compute_61\" -gencode=arch=compute_70,code=\"sm_70,compute_70\"")
    elseif(CUDA_VERSION VERSION_GREATER "10.9" AND CUDA_VERSION VERSION_LESS "12.0")
      message("Using CUDA version >=11.0 and <12.0")
      list(APPEND CUDA_NVCC_FLAGS "-use_fast_math -O3 -D_WITHMR -D_WITHGPU -gencode=arch=compute_35,code=\"sm_35,compute_35\" -gencode=arch=compute_50,code=\"sm_50,compute_50\" -gencode=arch=compute_52,code=\"sm_52,compute_52\" -gencode=arch=compute_61,code=\"sm_61,compute_61\" -gencode=arch=compute_70,code=\"sm_70,compute_70\" -gencode=arch=compute_75,code=\"sm_75,compute_75\" -gencode=arch=compute_80,code=\"sm_80,compute_80\"")
    else()
      message("Using CUDA version >=12.0")
      list(APPEND CUDA_NVCC_FLAGS "-use_fast_math -O3 -D_WITHMR -D_WITHGPU -Wno-deprecated-gpu-targets -DCUB_IGNORE_DEPRECATED_CPP_DIALECT -DTHRUST_IGNORE_DEPRECATED_CPP_DIALECT -gencode=arch=compute_61,code=\"sm_61,compute_61\" -gencode=arch=compute_70,code=\"sm_70,compute_70\" -gencode=arch=compute_75,code=\"sm_75,compute_75\" -gencode=arch=compute_80,code=\"sm_80,compute_80\" -gencode=arch=compute_86,code=\"sm_86,compute_86\" -gencode=arch=compute_89,code=\"sm_89,compute_89\" -gencode=arch=compute_90,code=\"sm_90,compute_90\"")
    endif()
  else()
    message("CUDA Libraries were not found.")
  endif(CUDA_FOUND)
endif(ENABLE_CUDA)

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
    LINK_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../lib/vs2022)
  #endif()
endif()

#------------------------------------------------------------------
# Compile MoorDynPlus library using existing Makefile
#------------------------------------------------------------------
# include(ExternalProject)
# if(ENABLE_MOORDYNPLUS)
  # if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  #   ExternalProject_Add(LibDSphMoorDynPlus
  #     SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Source_DSphMoorDynPlus
  #     CONFIGURE_COMMAND ""
  #     BUILD_COMMAND make -f ${CMAKE_CURRENT_SOURCE_DIR}/Source_DSphMoorDynPlus/Makefile_MoorDynPlus
  #     BUILD_IN_SOURCE 1
  #     INSTALL_COMMAND /bin/sh -c "rm -rf ${CMAKE_CURRENT_SOURCE_DIR}/*.o ${CMAKE_CURRENT_SOURCE_DIR}/Source_DSphMoorDynPlus/*.o"
  #     BUILD_BYPRODUCTS ${CMAKE_CURRENT_SOURCE_DIR}/../lib/linux_gcc/libdsphmoordynplus_64.a
  #   )
  # endif()
# endif(ENABLE_MOORDYNPLUS)
#------------------------------------------------------------------
# Compile MoorDynPlus library using existing CMake
#------------------------------------------------------------------
if(ENABLE_MOORDYNPLUS)
  set(MOORDYNPLUS_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Source_DSphMoorDynPlus)
  set(MOORDYNPLUS_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/Build_DSphMoorDynPlus)
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      message(STATUS "Building MoorDynPlus using CMake (GNU)")
      add_subdirectory(${MOORDYNPLUS_SOURCE_DIR} ${MOORDYNPLUS_BINARY_DIR})
  elseif(MSVC)
      message(STATUS "Building MoorDynPlus using CMake (MSVC)")
      add_subdirectory(${MOORDYNPLUS_SOURCE_DIR} ${MOORDYNPLUS_BINARY_DIR})
  endif()
endif(ENABLE_MOORDYNPLUS)

#--------------------------------------------
# OpenMP flags
#---------------------------------------------
find_package(OpenMP)

if (OPENMP_FOUND)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

#------------------------------------------------------------------
# Binaries
#------------------------------------------------------------------
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../../bin/${PLATFORM}) # Generic output directory
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})         # Release output directory
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})           # Debug output directory

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  add_executable(DualSPHysics5.4CPU_linux64 ${OBJXML} ${OBJSPHMOTION} ${OBCOMMON} ${OBCOMMONDSPH} ${OBSPH} ${OBSPHSINGLE} ${OBWAVERZ} ${OBWAVERZGPU}  ${OBCHRONO} ${OBMOORDYNPLUS} ${OBINOUT} ${OBMESH} ${OBVRES} ${OBFLEXSTRUC})
  add_dependencies(DualSPHysics5.4CPU_linux64 LibDSphMoorDynPlus)
  install(TARGETS	DualSPHysics5.4CPU_linux64 DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
  if (CUDA_FOUND)
    cuda_add_executable(DualSPHysics5.4_linux64 ${OBJXML} ${OBJSPHMOTION} ${OBCOMMON} ${OBCOMMONDSPH} ${OBSPH} ${OBSPHSINGLE} ${OBCOMMONGPU} ${OBSPHGPU} ${OBSPHSINGLEGPU} ${OBCUDA} ${OBWAVERZ} ${OBWAVERZCUDA} ${OBCHRONO} ${OBMOORDYNPLUS} ${OBINOUT} ${OBINOUTGPU} ${OBVRESGPU} ${OBWAVERZGPU} ${OBINOUTGPU} ${OBMESH} ${OBVRES} ${OBFLEXSTRUC})
    add_dependencies(DualSPHysics5.4_linux64 LibDSphMoorDynPlus)
    install(TARGETS DualSPHysics5.4_linux64 DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
  endif(CUDA_FOUND)
elseif(MSVC) 
  add_executable(DualSPHysics5.4CPU_win64 ${OBJXML} ${OBJSPHMOTION} ${OBCOMMON} ${OBCOMMONDSPH} ${OBSPH} ${OBSPHSINGLE} ${OBWAVERZ} ${OBWAVERZGPU}  ${OBCHRONO} ${OBMOORDYNPLUS} ${OBINOUT} ${OBMESH} ${OBVRES} ${OBFLEXSTRUC})
  add_dependencies(DualSPHysics5.4CPU_win64 LibDSphMoorDynPlus)
  install(TARGETS	DualSPHysics5.4CPU_win64 DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
  if (CUDA_FOUND)
    cuda_add_executable(DualSPHysics5.4_win64 ${OBJXML} ${OBJSPHMOTION} ${OBCOMMON} ${OBCOMMONDSPH} ${OBSPH} ${OBSPHSINGLE} ${OBCOMMONGPU} ${OBSPHGPU} ${OBSPHSINGLEGPU} ${OBCUDA} ${OBWAVERZ} ${OBWAVERZCUDA} ${OBCHRONO} ${OBMOORDYNPLUS} ${OBINOUT} ${OBINOUTGPU} ${OBVRESGPU} ${OBWAVERZGPU} ${OBINOUTGPU} ${OBMESH} ${OBVRES} ${OBFLEXSTRUC})
    add_dependencies(DualSPHysics5.4_win64 LibDSphMoorDynPlus)
    install(TARGETS DualSPHysics5.4_win64 DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
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

#------------------------------------------------------------------
# Coupled libraries
#------------------------------------------------------------------
# Chrono Engine
if(ENABLE_CHRONO)
  message(STATUS "Using ChronoEngine Module")
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(LINKER_FLAGS ${LINKER_FLAGS} ChronoEngine dsphchrono)
  elseif(MSVC) 
    set(LINKER_FLAGS ${LINKER_FLAGS} dsphchrono)
  endif()
else()
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DDISABLE_CHRONO")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDISABLE_CHRONO")
endif(ENABLE_CHRONO)

# WaveGen
if(ENABLE_WAVEGEN)
  message(STATUS "Using WaveGen module")
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(LINKER_FLAGS ${LINKER_FLAGS} jwavegen_64)
  elseif(MSVC) 
    set(LINKER_FLAGS ${LINKER_FLAGS} LibJWaveGen_x64_v143_Release)
  endif()
else()
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DDISABLE_WAVEGEN")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDISABLE_WAVEGEN")
endif(ENABLE_WAVEGEN)

# MoorDynPlus
if(ENABLE_MOORDYNPLUS)
  message(STATUS "Using MoorDynPlus module")
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(LINKER_FLAGS ${LINKER_FLAGS} dsphmoordynplus_64)
  elseif(MSVC) 
    set(LINKER_FLAGS ${LINKER_FLAGS} LibDSphMoorDynPlus_x64_v143_Release)
  endif()
else()
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DDISABLE_MOORDYNPLUS")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDISABLE_MOORDYNPLUS")
endif(ENABLE_MOORDYNPLUS)



#------------------------------------------------------------------
# Linker flags
#------------------------------------------------------------------
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  target_link_libraries(DualSPHysics5.4CPU_linux64 ${LINKER_FLAGS})
  set_target_properties(DualSPHysics5.4CPU_linux64 PROPERTIES COMPILE_FLAGS "-use_fast_math -O3 -D_WITHMR -fPIC -std=c++0x")
  
  if (CUDA_FOUND)
    target_link_libraries(DualSPHysics5.4_linux64 ${LINKER_FLAGS})
    set_target_properties(DualSPHysics5.4_linux64 PROPERTIES COMPILE_FLAGS "-use_fast_math -O3 -D_WITHMR -D_WITHGPU -fPIC -std=c++0x")
  endif(CUDA_FOUND)
  
elseif(MSVC)
  # MSVC 2022
  set_target_properties(DualSPHysics5.4CPU_win64 PROPERTIES COMPILE_FLAGS "/D _WITHMR")
  target_link_libraries(DualSPHysics5.4CPU_win64  ${LINKER_FLAGS})
  if (CUDA_FOUND)
    set_target_properties(DualSPHysics5.4_win64 PROPERTIES COMPILE_FLAGS "/D _WITHMR /D _WITHGPU")
    target_link_libraries(DualSPHysics5.4_win64 ${LINKER_FLAGS})
  endif(CUDA_FOUND)

  SET(CUDA_PROPAGATE_HOST_FLAGS OFF CACHE BOOL "Propagate C/CXX Flags and friends to the host compiler in NVCC via -Xompile  " FORCE)
  
  foreach(CPP_FLAGS CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO)
    if(${CPP_FLAGS} MATCHES "/MD")
      string(REGEX REPLACE "/MD" "/MT" ${CPP_FLAGS} "${${CPP_FLAGS}}")
    endif(${CPP_FLAGS} MATCHES "/MD")
  endforeach(CPP_FLAGS)
endif()
