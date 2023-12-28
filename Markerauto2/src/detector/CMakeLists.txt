
#set(DETECTOR_SRCS )
aux_source_directory(. DETECTOR_SRCS)

MESSAGE(STATUS "building detector")

include_directories(${PROJECT_SOURCE_DIR})

add_library (detector STATIC ${DETECTOR_SRCS})

target_link_libraries(detector levmar clapack cblas ${OpenCV_LIBS})

set_target_properties(detector PROPERTIES OUTPUT_NAME "detector")
set_target_properties(detector PROPERTIES CLEAN_DIRECT_OUTPUT 1)
