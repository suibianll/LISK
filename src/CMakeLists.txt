file(GLOB LISK_SRC "*.cc" "*.h" "*.c")

add_library(lisk SHARED ${LISK_SRC})
target_link_libraries(lisk glog)
# target_include_directories(lisk PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(lisk PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)
INSTALL(TARGETS lisk 
        ARCHIVE DESTINATION ${CMAKE_SOURCE_DIR}
        )