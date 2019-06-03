#!/bin/bash

ssh-agent bash

ssh-add ~/.ssh/id_test_rsa

ssh -T git@github.com
