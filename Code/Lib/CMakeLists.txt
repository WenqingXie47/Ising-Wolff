
set(LIB_NAME Ising)

set(LIB_SRCS
    util.cpp
    gridState.cpp 
    isingModel.cpp 
    monteCarlo.cpp
)

set(LIB_HDRS
    util.h
    gridState.h
    isingModel.h 
    monteCarlo.h
)

add_library(${LIB_NAME} ${LIB_HDRS} ${LIB_SRCS})