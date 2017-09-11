if(LIGHT)
SET(PRJ_SRC

abstractsource.cpp

fitstablesource.cpp
fitsimagesource.cpp


vosourcenew.cpp

hdf5source.cpp

asciisource.cpp
csvsource.cpp

binsource.cpp
flysource.cpp

gadgetsource.cpp

rawgridsource.cpp
rawpointssource.cpp
muportalsource.cpp
ramsessource.cpp

)

else(LIGHT)

SET(PRJ_SRC

abstractsource.cpp

fitstablesource.cpp
fitsimagesource.cpp

vosource.cpp
xmlsource.cpp

vosourcenew.cpp

hdf5source.cpp

asciisource.cpp
csvsource.cpp

binsource.cpp
flysource.cpp

gadgetsource.cpp

rawgridsource.cpp
rawpointssource.cpp
muportalsource.cpp
ramsessource.cpp


)
endif(LIGHT)

