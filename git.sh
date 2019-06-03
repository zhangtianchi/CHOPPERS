#!/bin/bash

git init

git add .

git commit -m “first”

git remote rm origint

git remote add origin git@github.com:zhangtianchi/CHOPPERS.git

git push -f origin master

