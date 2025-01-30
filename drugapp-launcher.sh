#!/bin/bash

print_help() {
  echo "Usage: $0 [--build [--no-cache]|--remove|--up|--down|--parameters|--help]"
  echo
  echo "Commands:"
  echo "  --build        Build the Docker image"
  echo "  --build --no-cache"
  echo "                 Build the Docker image without using the cache"
  echo "  --up           Start the services in detached mode"
  echo "  --down         Stop the Docker containers"
  echo "  --remove       Stop and remove the Docker containers, image, and volumes"
  echo "  --logs         Returns logs from the currently running instance"
  echo "  --help         Display this help message"
  echo "  --parameters \"(MONARCH,DGIDB,DRUGSIM,NETWORKMODEL,NEGSAMPLES,MAIN)\""
  echo "                    Set the number of instances for each service"
  echo "                    Default is (1,2,2,2,2,4) if not specified"
  echo "                    Example: --parameters \"(2,1,1,1,1,1)\""
}

update_parameters() {
    local params="${1:-"(1,2,2,2,2,4)"}"

    # Remove parentheses and split
    params_clean="${params//[()]/}"
    IFS=',' read -r -a PARAM_ARRAY <<< "$params_clean"

    # Path to Docker_parameters.txt
    PARAMS_FILE="./main/app/Docker_parameters.txt"

    # Create or overwrite the file
    cat > "$PARAMS_FILE" << EOF
MONARCH_INSTANCES=${PARAM_ARRAY[0]:-1}
DGIDB_INSTANCES=${PARAM_ARRAY[1]:-2}
DRUGSIMILARITY_INSTANCES=${PARAM_ARRAY[2]:-2}
NETWORKMODEL_INSTANCES=${PARAM_ARRAY[3]:-2}
NEGSAMPLE_INSTANCES=${PARAM_ARRAY[4]:-2}
MAIN_INSTANCES=${PARAM_ARRAY[5]:-4}
EOF

    echo "Updated Docker parameters:"
    cat "$PARAMS_FILE"
}

# Check for parameters first
if [[ "$1" == "--parameters" ]]; then
    update_parameters "$2"
    shift 2
fi

case "$1" in
  --build)
    if [ "$2" == "--no-cache" ]; then
      docker-compose build --no-cache
    else
      docker-compose build
    fi
    docker-compose up -d
    ;;
  --up)
    docker-compose up -d
    ;;
  --down)
    docker-compose down
    ;;
  --remove)
    docker-compose down
    docker image rm drugrepurposing-main:latest
    docker volume prune
    docker system prune --all --force --volumes
    ;;
  --logs)
    docker-compose logs -f app
    ;;
  --help)
    print_help
    ;;
  *)
    print_help
    exit 1
    ;;
esac
