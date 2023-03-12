cmake_minimum_required(VERSION 3.0.0)
project(lets_be_rational VERSION 0.1.0)

include(CTest)
enable_testing()

set(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
set(HEADER_FILES ${DIRECTORY}/rationalcubic.h ${DIRECTORY}/normaldistribution.h ${DIRECTORY}/lets_be_rational.h)
set(SOURCE_FILES ${DIRECTORY}/erf_cody.cpp ${DIRECTORY}/rationalcubic.cpp ${DIRECTORY}/normaldistribution.cpp ${DIRECTORY}/lets_be_rational.cpp)
add_library(lets_be_rational STATIC ${HEADER_FILES} ${SOURCE_FILES})
target_include_directories(lets_be_rational PRIVATE ${DIRECTORY})
add_definitions(-fPIC -DNDEBUG -Ofast)


set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
