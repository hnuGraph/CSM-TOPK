set(UTILS globals.cpp Log.cpp Timer.cpp)
add_library(util SHARED ${UTILS})
target_include_directories(util PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})