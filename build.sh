#!/usr/bin/env bash
set -o errexit

pip install --upgrade pip
pip install --no-cache-dir numpy scipy pandas scikit-learn
pip install --no-cache-dir -r requirements.txt