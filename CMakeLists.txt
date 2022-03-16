project( gmshr )

cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

if(DEFINED CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose build type")
else()
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose build type")
endif()

# Compiladores MPI
set(CMAKE_C_COMPILER /usr/bin/mpicc)
set(CMAKE_CXX_COMPILER /usr/bin/mpic++)

# gmsh

find_library(BOOST_FS boost_filesystem)
find_library(GMSH_LIB gmsh)
find_path(GMSH_INC gmsh.h)
include_directories(${GMSH_INC})


# OpenMP
option (USE_OpenMP "Use OpenMP" ON)
IF(USE_OpenMP)
  find_package(OpenMP)
  IF(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  ENDIF()
ENDIF()

# Compilar arquivos-fonte dentro da biblioteca estática libModules.a
ADD_LIBRARY(Modules
  src/Mesh.cpp 
)

# Indicar bibliotecas compartilhadas no sistema
target_link_libraries(Modules ${GMSH_LIB})
target_link_libraries(Modules ${BOOST_FS})


# Compilar e lincar executável final
add_executable( gmshr src/main.cpp  )
target_link_libraries( gmshr Modules )
