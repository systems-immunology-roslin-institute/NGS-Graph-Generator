#! /bin/bash

TAB_FILE=$1
BASE_NAME=$(basename "$TAB_FILE")
NO_EXT="${BASE_NAME%.*}"

awk {'print ">"$3"\n"$4'} $TAB_FILE | sed -e '1,2d'
