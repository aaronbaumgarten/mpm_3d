set(MATERIALS
        isolin.cpp
        g_local_mu2.cpp
        isolin_wheel.cpp
        viscoelastic_polymer_liquid.cpp
        viscofluid.cpp
        isolin_wheel_cst_omega.cpp
        gravitational_body.cpp
        g_local_mu2_full.cpp
        )

foreach(matc ${MATERIALS})
    string(REGEX REPLACE "\\.cpp$" "" matso ${matc})
    add_library(${matso} SHARED ${matc})
    set_target_properties(${matso} PROPERTIES PREFIX "")
endforeach(matc)
