set(MATERIALS
        isolin.cpp
        g_local_mu2.cpp
        viscofluid.cpp
        slurry_viscofluid.cpp
        slurry_mu_iv.cpp
        slurry_mu_iv_int.cpp
        slurry_viscofluid_avg.cpp
        slurry_mu_im_int.cpp
        )

foreach(matc ${MATERIALS})
    string(REGEX REPLACE "\\.cpp$" "" matso ${matc})
    add_library(${matso} SHARED ${matc})
    set_target_properties(${matso} PROPERTIES PREFIX "")
endforeach(matc)