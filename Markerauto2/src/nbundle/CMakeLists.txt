aux_source_directory(. BUNDLE_SRCS)

MESSAGE(STATUS "building newbundle")
# find_package(Ceres REQUIRED)
# set(CMAKE_THREAD_LIBS_INIT "-lpthread")
# find_package(Ceres REQUIRED)
# include_directories(${CERES_INCLUDE_DIRS} )

include_directories(${PROJECT_SOURCE_DIR})
add_library (nbundle STATIC ${BUNDLE_SRCS})

# target_link_libraries(nbundle sba dataf util matrix m)
target_link_libraries(nbundle sba dataf util matrix m ${CERES_LIBRARIES})
# target_link_libraries(nbundle ${CERES_LIBRARIES})

set_target_properties(nbundle PROPERTIES OUTPUT_NAME "bundle")
set_target_properties(nbundle PROPERTIES CLEAN_DIRECT_OUTPUT 1)
