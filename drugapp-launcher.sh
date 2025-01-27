#!/bin/bash

print_help() {
  echo "Usage: $0 [--build [--no-cache]|--remove|--up|--down|--help]"
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
}

case "$1" in
  --build)
    if [ "$2" == "--no-cache" ]; then
      docker-compose build --no-cache
    else
      docker-compose build
    fi
    docker-compose up -d  # Add -d flag for detached mode
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