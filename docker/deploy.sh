#!/bin/bash

echo "./deploy.sh $*" > redeploy.sh
chmod +x redeploy.sh

ARGS="$*"
algorithm=smina
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
  esac
  shift
done

existing=$(docker ps -aq -f name=docking-$algorithm)
if [ -n "$existing" ]; then
    echo "removing existing container"
    docker rm -f $existing
fi

docker run -d \
--name docking-$algorithm \
--restart unless-stopped \
-e ARGS="${ARGS[*]}" \
docking-$algorithm
