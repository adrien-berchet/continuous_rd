# Generate the executable
ADD_EXECUTABLE(
    continuous_rd.exe
    continuous_rd.cpp
)

IF( GPU )
  TARGET_COMPILE_DEFINITIONS(continuous_rd.exe PRIVATE -DWITH_GPU)
ENDIF( GPU )

IF( PROFILING )
    TARGET_COMPILE_DEFINITIONS(continuous_rd.exe PRIVATE -DPROFILING)
ENDIF( PROFILING )

TARGET_LINK_LIBRARIES(continuous_rd.exe
                      m
                      stdc++
                      boost_system
                      boost_chrono
                      boost_timer
                      boost_filesystem
                      )

IF(THREAD)
    TARGET_LINK_LIBRARIES(continuous_rd.exe
                          pthread
                          )
ENDIF(THREAD)

INSTALL (TARGETS continuous_rd.exe DESTINATION bin)
