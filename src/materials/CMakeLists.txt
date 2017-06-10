set(MATERIALS
        isolin.cpp
        g_local_mu2.cpp
        viscofluid.cpp
        slurry_viscofluid.cpp
        )

foreach(matc ${MATERIALS})
    string(REGEX REPLACE "\\.cpp$" "" matso ${matc})
    add_library(${matso} SHARED ${matc})
    set_target_properties(${matso} PROPERTIES PREFIX "")
endforeach(matc)