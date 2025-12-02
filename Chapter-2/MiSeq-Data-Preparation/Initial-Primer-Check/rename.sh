#!/bin/bash

for file in *.collapsed; do mv "$file" "${file%.txt}.fastq"; 

done
