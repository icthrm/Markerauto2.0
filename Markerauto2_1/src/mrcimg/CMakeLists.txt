
#set(KEYMATCH_SRCS match_core.cpp)
aux_source_directory(. MRC_IMG_SRCS)

MESSAGE(STATUS "building mrcimg")

add_library (mrcimg STATIC ${MRC_IMG_SRCS})

set(OpenCV_DIR /home/xzh-vm/opencv-4.5.5/build)
find_package(OpenCV 4.5.5 REQUIRED)

MESSAGE(STATUS "libraries:${OpenCV_LIBS}")
# MESSAGE(STATUS "libraries:${OpenCV_INCLUDE_DIRS}")

# 
# target_link_libraries(mrcimg /usr/local/opencv45/lib/libopencv_calib3d.so /usr/local/opencv45/lib/libopencv_core.so /usr/local/opencv45/lib/libopencv_features2d.so /usr/local/opencv45/lib/libopencv_flann.so /usr/local/opencv45/lib/libopencv_highgui.so /usr/local/opencv45/lib/libopencv_imgproc.so /usr/local/opencv45/lib/libopencv_ml.so /usr/local/opencv45/lib/libopencv_objdetect.so /usr/local/opencv45/lib/libopencv_ocl.so /usr/local/opencv45/lib/libopencv_photo.so /usr/local/opencv45/lib/libopencv_stitching.so /usr/local/opencv45/lib/libopencv_superres.so /usr/local/opencv45/lib/libopencv_ts.a /usr/local/opencv45/lib/libopencv_video.so /usr/local/opencv45/lib/libopencv_videostab.so -lrt -lpthread -lm -ldl)
# libopencv_con
target_link_libraries(mrcimg ${OpenCV_LIBS})#opencv_core opencv_imgproc opencv_calib3d opencv_video opencv_features2d 
# 	opencv_ml opencv_highgui opencv_objdetect opencv_contrib opencv_legacy)

set_target_properties(mrcimg PROPERTIES OUTPUT_NAME "mrcimg")
set_target_properties(mrcimg PROPERTIES CLEAN_DIRECT_OUTPUT 1)
