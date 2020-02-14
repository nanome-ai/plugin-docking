if [ "$(docker ps -aq -f name=docking)" != "" ]; then
        # cleanup
        echo "removing exited container"
        docker rm -f docking
fi

docker run -d \
--restart always \
-e ARGS="$*" \
--name docking docking
