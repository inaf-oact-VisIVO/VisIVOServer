if(GLITE)
SET(PRJ_SRC

startFilter.cpp

visivoutils.cpp

commandline.cpp
optionssetter.cpp

luteditor.cpp
color.cpp
extendedglyph3d.cpp

glitegw.cpp


)
else(GLITE)
SET(PRJ_SRC


startFilter.cpp
visivoutils.cpp

commandline.cpp
optionssetter.cpp

luteditor.cpp
color.cpp
extendedglyph3d.cpp
glitegw.cpp

)
endif(GLITE)