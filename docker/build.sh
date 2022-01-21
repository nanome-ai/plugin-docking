#!/bin/bash

algorithm=smina
cachebust=0
while [ $# -gt 0 ]; do
  case $1 in
    --algorithm )
      shift
      if ! [[ "$1" =~ ^(smina|autodock4)$ ]]; then
        echo "Algorithm must be either smina or autodock4"
        exit 1
      fi
      algorithm=$1
      ;;
    -u | --update ) cachebust=1 ;;
  esac
  shift
done

if [ ! -f ".cachebust" ] || (($cachebust)); then
  date +%s > .cachebust
fi

cachebust=`cat .cachebust`
docker build -f Dockerfile --build-arg CACHEBUST=$cachebust --build-arg ALGORITHM=$algorithm -t docking-$algorithm:latest ..
