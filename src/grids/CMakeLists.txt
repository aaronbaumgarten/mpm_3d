set(GRIDS
        cartesian.cpp
        cartesian_periodic.cpp
        gmsh_tri2D.cpp
        cartesian_custom.cpp
        cartesian_cubic.cpp
        cartesian_cubic_custom.cpp
        )

foreach(gridc ${GRIDS})
    string(REGEX REPLACE "\\.cpp$" "" gridso ${gridc})
    add_library(${gridso} SHARED ${gridc})
    set_target_properties(${gridso} PROPERTIES PREFIX "")
endforeach(gridc)