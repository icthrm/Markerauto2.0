project(markerauto-pb)

cmake_minimum_required(VERSION 2.6)


#the store of data in clapack or cblas is column main sequence

set(ANN_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/lib/ann_1.2_char/include)
set(CLAPACK_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/lib/clapack_3.1/INCLUDE)
set(MRC_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/lib/mrc)
set(BLASWRAP_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/lib/clapack_3.1/BLAS/WRAP)
set(CBLAS_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/lib/cblas/INCLUDE)
set(SBA_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/lib/sba_1.6)
set(MATRIX_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/lib/matrix)
 
# set(OPENCV_INCLUDE_DIR "usr/local/opencv45/include/opencv4")
set(OpenCV_DIR "usr/local/opencv45/include/opencv4")
set(CMAKE_THREAD_LIBS_INIT "-lpthread")
find_package(Ceres REQUIRED)

set(CERES_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/lib/ceres)
# set(OpenCV_LIBS "usr/local/opencv45/lib")
# 寻找OpenCV库
find_package( OpenCV 4.5.5 REQUIRED )
# 添加头文件
# include_directories( ${OpenCV_INCLUDE_DIRS} )
# # 寻找OpenCV库
# find_package( OpenCV 4 REQUIRED )
# set(OPENCV_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/lib/opencv_2.4.1 
#     ${PROJECT_SOURCE_DIR}/lib/opencv_2.4.1/opencv ${PROJECT_SOURCE_DIR}/lib/opencv_2.4.1/opencv2 )
# 可执行程序
# 链接OpenCV库
set(LEVMAR_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/lib/levmar_2.5)


include_directories(${PROJECT_SOURCE_DIR}/lib ${ANN_INCLUDE_DIR} ${PROJECT_SOURCE_DIR}/src 
    ${CLAPACK_INCLUDE_DIR} ${BLASWRAP_INCLUDE_DIR} ${CBLAS_INCLUDE_DIR} ${MRC_INCLUDE_DIR} ${SBA_INCLUDE_DIR}
    ${MATRIX_INCLUDE_DIR} ${OPENCV_INCLUDE_DIR} ${LEVMAR_INCLUDE_DIR})
link_directories (${PROJECT_BINARY_DIR}/lib)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

add_subdirectory(lib)
add_subdirectory(src)
