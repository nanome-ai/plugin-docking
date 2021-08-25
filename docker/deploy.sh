#!/bin/bash

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$parent_path"

echo "./deploy.sh $*" > redeploy.sh
chmod +x redeploy.sh

existing=$(docker ps -aq -f name=docking)
if [ -n "$existing" ]; then
    echo "removing existing container"
    docker rm -f $existing
fi

ARGS=("smina" $*)

docker run -d \
--name docking \
--restart unless-stopped \
-e ARGS="${ARGS[*]}" \
docking
