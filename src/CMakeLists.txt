file(GLOB dbfe_SRC
    "*.h"
    "*.cpp"
)
message(${RCPPLOCATION})
SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${RCPPLOCATION}" )
#include_directories(${RCPPLOCATION})
add_executable(dbfexe ${dbfe_SRC})
