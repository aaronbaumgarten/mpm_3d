set(DRIVERS
        default_driver.cpp
        ballistic_impactor.cpp
        shear_driver.cpp
        user_def_grav_driver.cpp
        shear_traction_driver.cpp
        launched_projectile_driver.cpp
        shear_traction_driver_subsample.cpp
        )

foreach(dric ${DRIVERS})
    string(REGEX REPLACE "\\.cpp$" "" driso ${dric})
    add_library(${driso} SHARED ${dric})
    set_target_properties(${driso} PROPERTIES PREFIX "")
endforeach(dric)