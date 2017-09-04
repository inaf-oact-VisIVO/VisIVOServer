if(SPLVISIVO)
SET(PRJ_SRC

pipe.cpp


pointspipe.cpp
pointssmoothpipe.cpp
splotchpipe.cpp
splotchpipecamera.cpp
volumepipe.cpp
slicerpipe.cpp
isosurfacepipe.cpp
vectorpipe.cpp
vtkimagepipe.cpp
)
else(SPLVISIVO)
SET(PRJ_SRC

pipe.cpp


pointspipe.cpp
pointssmoothpipe.cpp
volumepipe.cpp
slicerpipe.cpp
isosurfacepipe.cpp
vectorpipe.cpp
vtkimagepipe.cpp
)

endif(SPLVISIVO)
