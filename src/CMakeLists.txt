cmake_minimum_required(VERSION 3.6.3)
set(CMAKE_CXX_STANDARD 17)
set(SRCS main.cpp)
FIND_PACKAGE( OpenMP REQUIRED)
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
#增加当前文件夹下的子文件夹，model1,model2，model3文件夹为子文件夹，进入该文件夹执行其目录下的CMakeLists.txt文件
add_subdirectory(graph)
add_subdirectory(utils)
add_subdirectory(matching)

#设置可执行程序的依赖
add_executable(csm ${SRCS})
#设置可执行文件csm的动态库依赖，设置依赖时，cmake会自动导入对应库的头文件
target_link_libraries(csm graph util match)
#设置可执行文件csm的动态库依赖运行时的加载目录，使用本条命令后，csm每次都会从本文件统计目录的lib文件夹下加载so文件
set_target_properties(csm PROPERTIES INSTALL_RPATH "$ORIGIN/lib")
#add_compile_definitions(PRINT_DEBUG )
SET(CMAKE_BUILD_TYPE "Debug")
add_compile_definitions(COMPUTE_TIME)
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -Wall -g -ggdb")
