include_directories(${CMAKE_CURRENT_BINARY_DIR})

add_executable(mrc2pgm main.cpp)

find_package(OpenCV REQUIRED)

target_link_libraries(mrc2pgm mrcimg ${OpenCV_LIBS})
