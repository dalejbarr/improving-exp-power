#!/bin/bash -eu

cd /var/scripts && Rscript do-all-simulations.R $1 /output
