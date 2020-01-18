#!/bin/bash

ENV=production
IMAGE=docking
REPO=441665557124.dkr.ecr.us-west-1.amazonaws.com

# Sign in to AWS
aws ecr get-login --no-include-email --region us-west-1 | sh

# Build docker
docker build -f $IMAGE.Dockerfile -t $IMAGE:$ENV .

# Tag docker
docker tag $IMAGE:$ENV $REPO/$IMAGE:$ENV

# Upload docker to secured repo
docker push $REPO/$IMAGE:$ENV