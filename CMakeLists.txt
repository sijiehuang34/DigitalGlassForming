cmake_minimum_required(VERSION 3.5)
project(helloworld)

include(CTest)
enable_testing()

set(OpenCV_DIR C:/Users/Cindy Huang/projects/opencv/build)
find_package(OpenCV REQUIRED)
include_directories( ${OpenCV_INCLUDE_DIRS} )

add_executable(${PROJECT_NAME} helloworld.cpp)

target_link_libraries(${PROJECT_NAME} ${OpenCV_LIBS})

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})

include(CPack)