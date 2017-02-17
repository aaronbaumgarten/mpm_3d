set(BOUNDARIES
        boxBC.cpp
        boxBC2d.cpp
        floorBC.cpp
        floorBC2d.cpp
        syringeBC2d.cpp
        periodicBC2d.cpp
        periodicBC.cpp
        pipeBC2d.cpp
    )

foreach(bcc ${BOUNDARIES})
    string(REGEX REPLACE "\\.cpp$" "" bcso ${bcc})
    add_library(${bcso} SHARED ${bcc})
    set_target_properties(${bcso} PROPERTIES PREFIX "")
endforeach(bcc)
