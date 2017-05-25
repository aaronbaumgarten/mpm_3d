set(GRIDS
        cartesian.cpp
        )

foreach(gridc ${GRIDS})
    string(REGEX REPLACE "\\.cpp$" "" gridso ${gridc})
    add_library(${gridso} SHARED ${gridc})
    set_target_properties(${gridso} PROPERTIES PREFIX "")
endforeach(gridc)