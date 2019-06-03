#!/bin/bash

git init

git add .

git commit -m “second”

git remote rm origin

git remote add origin git@github.com:zhangtianchi/CHOPPERS.git
