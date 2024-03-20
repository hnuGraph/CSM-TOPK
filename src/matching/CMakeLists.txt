set(MATCHING matching.cpp CSMTopk.cpp SingleCandidate.cpp)
#FIND_PACKAGE( OpenMP REQUIRED)
#[[if(OPENMP_FOUND)
    message("OPENMP FOUND")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()]]
#设置编译match.so文件
add_library(match SHARED ${MATCHING})
#设置将当前文件夹作为依赖match库的头文件的搜索目录
target_include_directories(match PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
SET(CMAKE_BUILD_TYPE "Debug")
#PRINT_DEBUG LOG_TRACK RESULT_TRACK COMPUTE_TRACK
add_compile_definitions()

