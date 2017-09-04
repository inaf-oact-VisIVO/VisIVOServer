#!/bin/sh

LOCALDIR=$(pwd)

echo "**********************************************************"
echo "** Starting on $HOSTNAME at: "$(date +%Y%m%d-%H%M%S) 
echo "**********************************************************"
echo ""
echo "Content of Current Dir:"
echo "-----------------------"
ls -lrt $LOCALDIR
echo ""
echo "Compiling source code..."
echo "------------------------"
/usr/bin/make
echo ""
echo "Content of Current Dir:"
echo "-----------------------"
ls -lrt $LOCALDIR
echo ""
./myprogram
echo "**********************************************************"
echo "** Done on $HOSTNAME at: "$(date +%Y%m%d-%H%M%S) 
echo "**********************************************************"
