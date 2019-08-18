# Generate the executable
# CUDA_ADD_EXECUTABLE(
#     continuous_rd.exe
#     continuous_rd.cpp
#     continuous_rd_device.cu
# )
ADD_EXECUTABLE(
    continuous_rd.exe
    continuous_rd.cu
    # continuous_rd.cpp
    # continuous_rd_device.cu
)

SET_PROPERTY(TARGET continuous_rd.exe PROPERTY CUDA_STANDARD 14)

IF( GPU )
  TARGET_COMPILE_DEFINITIONS(continuous_rd.exe PRIVATE -DWITH_GPU)
ENDIF( GPU )

IF( PROFILING )
    TARGET_COMPILE_DEFINITIONS(continuous_rd.exe PRIVATE -DPROFILING)
ENDIF( PROFILING )

TARGET_LINK_LIBRARIES(continuous_rd.exe
                      m
                      stdc++
                      boost_chrono
                      boost_filesystem
                      boost_log
                      boost_log_setup
                      boost_system
                      boost_thread
                      boost_timer
                      pthread
                      )

INSTALL (TARGETS continuous_rd.exe DESTINATION bin)
