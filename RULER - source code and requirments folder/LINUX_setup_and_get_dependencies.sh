#!/bin/bash

set -e

# Check if Python3 is installed
if ! command -v python3 &> /dev/null
then
    echo "Python3 is not installed. Please install it first."
    exit 1
fi

# Check if pip is installed, install if missing
if ! command -v pip3 &> /dev/null
then
    echo "pip3 not found. Installing pip..."
    sudo apt update
    sudo apt install -y python3-pip
fi

# Create virtual environment
python3 -m venv venv

# Activate virtual environment
source venv/bin/activate

# Upgrade pip inside virtual environment
pip install --upgrade pip

# Install required libraries
pip install opencv-python numpy pandas scikit-image scipy

echo "All dependencies installed successfully in the virtual environment 'venv'."
echo "Activate it with: source venv/bin/activate"
