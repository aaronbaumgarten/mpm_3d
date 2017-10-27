set(CONTACTS
        contact_huang.cpp
        mixture_slurry.cpp
        contact_rigid.cpp
        contact_bonded.cpp
        mixture_slurry_explicit.cpp
        )

foreach(conc ${CONTACTS})
    string(REGEX REPLACE "\\.cpp$" "" conso ${conc})
    add_library(${conso} SHARED ${conc})
    set_target_properties(${conso} PROPERTIES PREFIX "")
endforeach(conc)