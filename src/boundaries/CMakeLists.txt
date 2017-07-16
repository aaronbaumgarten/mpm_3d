set(BOUNDARIES
        cartesian_box.cpp
        cartesian_smoothbox.cpp
        rigid_linear_motion.cpp
        cartesian_periodic_floor.cpp
        cartesian_floor_freewall.cpp
        cartesian_periodic_driven.cpp
        )

foreach(bcc ${BOUNDARIES})
    string(REGEX REPLACE "\\.cpp$" "" bcso ${bcc})
    add_library(${bcso} SHARED ${bcc})
    set_target_properties(${bcso} PROPERTIES PREFIX "")
endforeach(bcc)