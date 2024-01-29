cmake_minimum_required(VERSION 3.10)



project(IsingWolff VERSION 1.0)




add_subdirectory(Lib)


set(APPS
    test 
    main
)


foreach(_APP ${APPS})
    add_executable(${_APP} ${_APP}.cpp)
    target_link_libraries(${_APP} PUBLIC Ising)
    target_include_directories(${_APP} PUBLIC
        "${PROJECT_BINARY_DIR}"
        "${PROJECT_SOURCE_DIR}/Lib"
    )
    
endforeach()