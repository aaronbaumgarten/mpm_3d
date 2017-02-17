set(CONTACTS
        huang2d.cpp
        huang.cpp
        nocontact.cpp
        slurry.cpp
        )

foreach(conc ${CONTACTS})
    string(REGEX REPLACE "\\.cpp$" "" conso ${conc})
    add_library(${conso} SHARED ${conc})
    set_target_properties(${conso} PROPERTIES PREFIX "")
endforeach(conc)
