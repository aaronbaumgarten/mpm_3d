set(SOLVERS
        explicit_usl.cpp
        explicit_usl_es.cpp
        newton_cg.cpp
        newton_bicgstab.cpp
        )

foreach(solc ${SOLVERS})
    string(REGEX REPLACE "\\.cpp$" "" solso ${solc})
    add_library(${solso} SHARED ${solc})
    set_target_properties(${solso} PROPERTIES PREFIX "")
endforeach(solc)