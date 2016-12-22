set(MATERIALS
        isolin.cpp
        )

foreach(matc ${MATERIALS})
    string(REGEX REPLACE "\\.cpp$" "" matso ${matc})
    add_library(${matso} SHARED ${matc})
    set_target_properties(${matso} PROPERTIES PREFIX "")
endforeach(matc)
