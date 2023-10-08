
# [ref](https://www.hiroom2.com/2016/09/07/makefile%E3%82%92cmakelists-txt%E3%81%AB%E7%BD%AE%E3%81%8D%E6%8F%9B%E3%81%88%E3%82%8B/)
# [musl-gcc](https://github.com/actions-rs/toolchain/issues/102#issuecomment-1132620977)

cmake_minimum_required(VERSION 3.10)

#set(CMAKE_CXX_STANDARD 11)


#project("cmake example for Rust" LANGUAGES CXX C)

# before project() #[ref](https://stackoverflow.com/questions/45933732/how-to-specify-a-compiler-in-cmake)
#set(CMAKE_C_COMPILER "gcc")
#set(CMAKE_CXX_COMPILER "g++")
#set(CMAKE_CXX_STANDARD 11)

#project("cmake-example")
#project("cmake-example" LANGUAGES CXX)
#project("cmake-example" CXX)
project("pgenlib")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11") # no effect

# INC           =  -I./simde/ -I./include/
# -l /simde/ too
#include_directories(${CMAKE_CURRENT_SOURCE_DIR})
#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/simde)



# what for?
AUX_SOURCE_DIRECTORY(${CMAKE_CURRENT_SOURCE_DIR} source)
AUX_SOURCE_DIRECTORY(${CMAKE_CURRENT_SOURCE_DIR}/include source)
AUX_SOURCE_DIRECTORY(${CMAKE_CURRENT_SOURCE_DIR}/simde source)

set(TARGET ${CMAKE_PROJECT_NAME})
#set(TARGET cmake-example)

# [ref](https://qiita.com/iwatake2222/items/5800fda029019ce8a276)
# For OpenMP
find_package(OpenMP REQUIRED)
if(OpenMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()
# [ref](https://stackoverflow.com/questions/12399422/how-to-set-linker-flags-for-openmp-in-cmakes-try-compile-function)
# error
# target_link_libraries(${TARGET} PRIVATE OpenMP::OpenMP_CXX)



add_library(${TARGET} STATIC ${source})




install (TARGETS ${TARGET} DESTINATION .)



## for main.cpp
## https://qiita.com/shohirose/items/45fb49c6b429e8b204ac
##add_library(${TARGET} STATIC pgenlibr_wrap.cpp)
#add_executable(a.out main.cpp)
#
#target_include_directories(a.out
#  PUBLIC
#	${CMAKE_CURRENT_SOURCE_DIR} 
#	${CMAKE_CURRENT_SOURCE_DIR}/include
#	${CMAKE_CURRENT_SOURCE_DIR}/simde
#)
