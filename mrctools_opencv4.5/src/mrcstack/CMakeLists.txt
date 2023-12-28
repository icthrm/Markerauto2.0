include_directories(${CMAKE_CURRENT_BINARY_DIR})

add_executable(mrcstack main.cpp)

find_package(OpenCV REQUIRED)

target_link_libraries(mrcstack mrcimg ${OpenCV_LIBS})
