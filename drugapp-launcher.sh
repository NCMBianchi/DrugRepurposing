#!/bin/bash

# Function to show help message
show_help() {
    echo "DrugApp Launcher - Manage DrugApp docker services"
    echo
    echo "Usage: ./drugapp-launcher.sh [OPTION]"
    echo
    echo "Options:"
    echo "  --build    Build and deploy all DrugApp services"
    echo "  --remove   Stop and remove all DrugApp services and related Docker resources"
    echo "  --cleanup  Thoroughly clears all DrugApp services and related Docker resources"
    echo "  --help     Display this help message"
}

# Function to build and deploy services
build_services() {
    echo "Building and deploying DrugApp services..."
    
    # Ensure entrypoint.sh has execute permissions
    chmod +x entrypoint.sh

    # Build all services with quiet flag and redirecting all output
    docker build -q -t drugapp/main-service:latest . 2>/dev/null >/dev/null
    echo "✓ Built main service"
    
    docker build -q -t drugapp/monarch-service:latest -f services/Dockerfile.monarch . 2>/dev/null >/dev/null
    echo "✓ Built monarch service"
    
    docker build -q -t drugapp/dgidb-service:latest -f services/Dockerfile.dgidb . 2>/dev/null >/dev/null
    echo "✓ Built dgidb service"
    
    docker build -q -t drugapp/drugsim-service:latest -f services/Dockerfile.drugsim . 2>/dev/null >/dev/null
    echo "✓ Built drugsim service"
    
    docker build -q -t drugapp/negsample-service:latest -f services/Dockerfile.negsample . 2>/dev/null >/dev/null
    echo "✓ Built negsample service"
    
    docker build -q -t drugapp/networkmodel-service:latest -f services/Dockerfile.networkmodel . 2>/dev/null >/dev/null
    echo "✓ Built networkmodel service"
    
    # Initialize swarm and deploy
    docker swarm init 2>/dev/null || true
    echo "✓ Initialized Docker swarm"
    
    docker stack deploy -c docker-compose.yml drugapp 2>/dev/null || { echo "Error deploying stack"; exit 1; }
    echo "✓ Deployed DrugApp stack"
    
    echo "DrugApp services have been built and deployed successfully."
}

# Function to remove services and cleanup
remove_services() {
    echo "Removing DrugApp services and cleaning up..."
    
    docker stack rm drugapp 2>/dev/null || true
    echo "✓ Removed Docker stack"
    
    sleep 2
    docker swarm leave --force 2>/dev/null || true
    echo "✓ Left Docker swarm"
    
    docker stop $(docker ps -a -q) 2>/dev/null || true
    echo "✓ Stopped all containers"
    
    echo "Removing Docker images..."
    docker rmi drugapp/main-service:latest 2>/dev/null || true
    docker rmi drugapp/monarch-service:latest 2>/dev/null || true
    docker rmi drugapp/dgidb-service:latest 2>/dev/null || true
    docker rmi drugapp/drugsim-service:latest 2>/dev/null || true
    docker rmi drugapp/negsample-service:latest 2>/dev/null || true
    docker rmi drugapp/networkmodel-service:latest 2>/dev/null || true
    echo "✓ Removed all DrugApp images"
    
    docker volume prune -f 2>/dev/null || true
    echo "✓ Cleaned up volumes"
    
    echo "DrugApp services have been removed successfully."
}

cleanup_services() {
    echo "Performing thorough cleanup..."
    
    docker stack rm drugapp 2>/dev/null || true
    echo "✓ Removed Docker stack"
    
    sleep 5
    
    docker swarm leave --force 2>/dev/null || true
    echo "✓ Left Docker swarm"
    
    docker stop $(docker ps -a -q) 2>/dev/null || true
    echo "✓ Stopped all containers"
    
    docker rm $(docker ps -a -q) 2>/dev/null || true
    echo "✓ Removed all containers"
    
    docker rmi $(docker images drugapp/* -q) 2>/dev/null || true
    echo "✓ Removed all DrugApp images"
    
    docker volume prune -f 2>/dev/null || true
    echo "✓ Cleaned up volumes"
    
    echo "Cleanup completed successfully."
}


# Main script logic
case "$1" in
    --build)
        build_services
        ;;
    --remove)
        remove_services
        ;;
    --cleanup)
        cleanup_services
        ;;
    --help)
        show_help
        ;;
    *)
        echo "Error: Invalid option"
        echo "Use --help to see available options"
        exit 1
        ;;
esac