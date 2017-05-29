set(SERIALIZERS
        default_vtk.cpp
        )

foreach(serc ${SERIALIZERS})
    string(REGEX REPLACE "\\.cpp$" "" serso ${serc})
    add_library(${serso} SHARED ${serc})
    set_target_properties(${serso} PROPERTIES PREFIX "")
endforeach(serc)