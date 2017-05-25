#!/bin/bash

for i in {1..10000}
do
time ls | wc -l
date
done > times.txt