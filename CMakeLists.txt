cmake_minimum_required(VERSION 3.13)
project(baseline)
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -Wall -g -ggdb")
#[[SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")]]
set(CMAKE_CXX_STANDARD 17)
#设置可执行文件的执行目录
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/build)
#设置库文件的输出目录
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/build/lib)

add_subdirectory(${PROJECT_SOURCE_DIR}/src)
add_compile_options(/std:c++latest)
target_link_libraries(csm pthread)


