set(CONTACTS
        contact_huang.cpp
        )

foreach(conc ${CONTACTS})
    string(REGEX REPLACE "\\.cpp$" "" conso ${conc})
    add_library(${conso} SHARED ${conc})
    set_target_properties(${conso} PROPERTIES PREFIX "")
endforeach(conc)