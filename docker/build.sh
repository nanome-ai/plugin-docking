#!/bin/bash

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$parent_path"

build_args=""
cachebust=0
while [ $# -gt 0 ]; do
  case $1 in
    -u | --update ) cachebust=1 ;;
    *) build_args=$build_args" $1" ;;
  esac
  shift
done

if [ ! -f ".cachebust" ] || (($cachebust)); then
  date +%s > .cachebust
fi

cachebust=`cat .cachebust`
docker build -f Dockerfile --build-arg CACHEBUST=$cachebust $build_args -t docking:latest ..
