set(GRAPH graph.cpp MatchRecord.cpp Neighbor.cpp Neighbor.h StarGraph.cpp ForwardNeighbor.cpp ForwardNeighbor.h Subgraph.cpp)
#设置编译graph.so文件
add_library(graph SHARED ${GRAPH})
#设置将当前文件夹作为依赖graph库的头文件的搜索目录
target_include_directories(graph PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
#add_compile_definitions(PRINT_DEBUG)
SET(CMAKE_BUILD_TYPE "Debug")
