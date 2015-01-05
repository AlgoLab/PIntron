#!/bin/bash
set -e
source /build/buildconfig
set -x

export LC_ALL=C
export DEBIAN_FRONTEND=noninteractive

apt-get update
apt-get install -y git-core build-essential python3 nginx ssh
