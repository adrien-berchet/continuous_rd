# Generate the executable
ADD_EXECUTABLE(
    continuous_rd.exe
    continuous_rd.cu
)

SET_PROPERTY(TARGET continuous_rd.exe PROPERTY CUDA_STANDARD 14)
SET_PROPERTY(TARGET continuous_rd.exe PROPERTY CUDA_SEPARABLE_COMPILATION ON)

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
