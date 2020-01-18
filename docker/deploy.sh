#!/bin/bash

USER=user
SERVER=localhost
SSH=$USER@$SERVER
TAG=latest

IMAGE=docking
REPO=441665557124.dkr.ecr.us-west-1.amazonaws.com

echo "" > remote-actions.txt

echo "sudo docker pull $REPO/$IMAGE:$TAG" >> remote-actions.txt
echo "sudo docker stop $IMAGE" >> remote-actions.txt
echo "sudo docker rm -f $IMAGE" >> remote-actions.txt
echo "sudo docker run \
    -d --restart unless-stopped \
    --name=$IMAGE  \
    -p 8888:8888 \
    $REPO/$IMAGE:$TAG" >> remote-actions.txt

ssh $SSH 'bash -s' < remote-actions.txt
rm remote-actions.txt