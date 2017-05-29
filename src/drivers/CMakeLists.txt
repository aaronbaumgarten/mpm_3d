set(DRIVERS
        default_driver.cpp
        )

foreach(dric ${DRIVERS})
    string(REGEX REPLACE "\\.cpp$" "" driso ${dric})
    add_library(${driso} SHARED ${dric})
    set_target_properties(${driso} PROPERTIES PREFIX "")
endforeach(dric)