set(SOLVERS
        explicit_usl.cpp
        explicit_usl_es.cpp
        explicit_usl_axisym2D.cpp
        explicit_usl_parallel.cpp
        newton_cg.cpp
        newton_bicgstab.cpp
        newton_bicgstabl.cpp
        newton_direct.cpp
        explicit_usl_cubic.cpp
        explicit_usl_delta_correctedmass.cpp
        explicit_usl_apic.cpp
        )

foreach(solc ${SOLVERS})
    string(REGEX REPLACE "\\.cpp$" "" solso ${solc})
    add_library(${solso} SHARED ${solc})
    set_target_properties(${solso} PROPERTIES PREFIX "")
endforeach(solc)