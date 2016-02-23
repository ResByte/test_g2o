cmake_minimum_required(VERSION 2.8.3)

project(ba_test)

LIST(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake_modules/")



find_package(cmake_modules REQUIRED)
find_package(Eigen REQUIRED)
find_package(SuiteSparse REQUIRED)
#check for OpenMP
find_package(OpenMP)
if(OPENMP_FOUND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  if(MSVC90 OR MSVC10)
    if(MSVC90)
      set(OPENMP_DLL VCOMP90)
    elseif(MSVC10)
      set(OPENMP_DLL VCOMP100)
    endif(MSVC90)
    set(CMAKE_SHARED_LINKER_FLAGS_DEBUG "${CMAKE_SHARED_LINKER_FLAGS_DEBUG} /DELAYLOAD:${OPENMP_DLL}D.dll")
    set(CMAKE_SHARED_LINKER_FLAGS_RELEASE "${CMAKE_SHARED_LINKER_FLAGS_RELEASE} /DELAYLOAD:${OPENMP_DLL}.dll")
  endif(MSVC)
else(OPENMP_FOUND)
  message (STATUS "OpenMP not found")
endif()

if (MSVC)
  Set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc /O2 ${SSE_FLAGS}")
else (MSVC)
  set(CMAKE_CXX_FLAGS "-O3 ${CMAKE_CXX_FLAGS} ${SSE_FLAGS}")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -O0 -g ${SSE_FLAGS}")
  set(CMAKE_CXX_FLAGS_RELEASE " ${CMAKE_CXX_FLAGS} -O3 ${SSE_FLAGS}")
endif (MSVC)


find_package(G2O REQUIRED)

link_directories(${G2O_LIBRARY_DIRS})
link_libraries(${G2O_LIBRARIES})


include_directories(
  ${EIGEN_INCLUDE_DIRS}
  ${G2O_INLCUDE_DIRS}
#  ${CHOLMOD_INLCUDE_DIR}
  ${CSPARSE_INCLUDE_DIR}
)

add_executable(ba_demo ba_test.cpp)

SET(G2O_LIBS g2o_cli g2o_ext_freeglut_minimal g2o_simulator g2o_solver_slam2d_linear g2o_types_icp g2o_types_slam2d g2o_core g2o_interface g2o_solver_csparse g2o_solver_structure_only g2o_types_sba g2o_types_slam3d g2o_csparse_extension g2o_opengl_helper g2o_solver_dense g2o_stuff g2o_types_sclam2d  g2o_ext_csparse g2o_parser g2o_solver_pcg g2o_types_data g2o_types_sim3 cxsparse )

target_link_libraries(ba_demo

#  ${CHOLMOD_LIBRARIES}	
  ${CSPARSE_LIBRARY}
  cholmod 
  cxsparse
  ${G2O_LIBS}
  ${catkin_LIBRARIES}	
)
