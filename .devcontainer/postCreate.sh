#!/bin/bash

conda env create -f .devcontainer/environment.yaml

touch $HOME/.bash_profile
conda init && echo "conda activate nm" >> $HOME/.bash_profile

