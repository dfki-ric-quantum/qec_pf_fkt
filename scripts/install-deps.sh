#!/bin/bash
set -e

# System dependencies
apt update
apt install -y \
  build-essential \
  libgmp-dev \
  libgmpxx4ldbl \
  libgsl-dev \
  libhdf5-dev \
  python3-dev \
  python3-pip \
  python3-venv

# Create and activate virtual environment
python3 -m venv .venv
source .venv/bin/activate

# Python requirements
pip install --upgrade pip
pip install -r requirements.txt

echo "Dependencies installed successfully!"
echo "Run this to activate the virtual environment:"
echo "source .venv/bin/activate"