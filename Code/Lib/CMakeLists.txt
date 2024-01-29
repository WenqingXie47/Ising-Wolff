
set(LIB_NAME Ising)

set(LIB_SRCS
    gridState.cpp 
    isingState.cpp 
    monteCarlo.cpp
)

set(LIB_HDRS
    gridState.h
    isingState.h 
    monteCarlo.h
)

add_library(${LIB_NAME} ${LIB_HDRS} ${LIB_SRCS})