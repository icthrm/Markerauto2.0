

# find_package(Opencv 4.5 REQUIRED)
# include_directories( ${OpenCV_INCLUDE_DIRS} )
MESSAGE(STATUS "building main")

include_directories(${PROJECT_SOURCE_DIR} ${EIGEN3_INCLUDE_DIR})

add_executable(markerauto main.cpp) #cpd_test.cpp)

# list(APPEND CMAKE_CXX_FLAGS "-std=c++0x")

target_link_libraries(markerauto detector util dataf modelmatch detector nbundle)
