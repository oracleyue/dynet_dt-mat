#!/bin/bash

# This script is a wraper to run matlab codes.
# Last modified on 23 Nov 2018


matlab=/usr/local/bin/matlab

$matlab -nosplash -nodesktop -r 'sARXMultiSNR; delete(gcp); quit'
$matlab -nosplash -nodesktop -r 'sARXMultiNodes; delete(gcp); quit'
