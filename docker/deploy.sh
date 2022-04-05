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

container_name="docking-$algorithm"
existing=$(docker ps -aq -f name=$container_name)
if [ -n "$existing" ]; then
    echo "removing existing container"
    docker rm -f $existing
fi

docker run -d \
--name $container_name \
--restart unless-stopped \
-h $(hostname)-$container_name \
-e ARGS="${ARGS[*]}" \
$container_name
