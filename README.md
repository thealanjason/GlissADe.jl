# Differentiable Free Surface Flow Solver

This is a self contained repository with code and documentation for the Research Internship work carried out by Tanish Jain at the Chair of Methods for Model-based Development in Computational Engineering, RWTH Aachen University.

## Repository Structure
This repositor is a combination repository that is meant to host both, the computational code and documentation. It contains the following two folders:
- src : The software assets are stored here
- docs : The documentation website assets are stored here
- examples: Some examples to start using the software

## Getting Started with Code
<!-- NOTE: Currently only Intel Processor based Linux is supported. If running Windows, please Windows Subsystem for Linux (WSL). -->

1. Install Julia using instructions from [here](https://julialang.org/downloads/)
2. Clone this repository using `git clone https://git.rwth-aachen.de/mbd_teaching/internship-tj.git`
3. Change directory to cloned folder
4. Setup the Environment using `julia setup.jl`

## Getting Started with Documentation Website

1. Install Quarto
    ```
    export QUARTO_VERSION=1.5.39
    sudo mkdir -p /opt/quarto/${QUARTO_VERSION}
    sudo curl -o quarto.tar.gz -L "https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-amd64.tar.gz"
    sudo tar -zxvf quarto.tar.gz -C "/opt/quarto/${QUARTO_VERSION}" --strip-components=1
    sudo ln -s /opt/quarto/${QUARTO_VERSION}/bin/quarto /usr/local/bin/quarto
    quarto install tinytex
    quarto check
    ```
    
2. Locally render documentation website
    ```
    cd scripts
    ./build.sh
    ```

## Useful Links
- [https://modernjuliaworkflows.github.io/](https://modernjuliaworkflows.github.io/) A non-exhaustive guide with best practices and useful tips for Julia-lang development
