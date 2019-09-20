set(PROJECT_NAME "blintsov_s_conj_gradient_task2_omp")

add_executable( ${PROJECT_NAME} main.cpp )

#set(MYCOMPILE_FLAGS "-fopenmp")
#set(CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} ${MYCOMPILE_FLAGS}")
#set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${MYCOMPILE_FLAGS}")
#set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

include_directories("/usr/local/include" "/usr/local/opt/llvm/include")
link_directories("/usr/local/lib" "/usr/local/opt/llvm/lib")

if ( USE_OMP )
    project( ${PROJECT_NAME} )
    message( STATUS "-- " ${PROJECT_NAME} )
else( USE_OMP )
    message( STATUS "-- ${PROJECT_NAME} - NOT BUILD!"  )
endif( USE_OMP )