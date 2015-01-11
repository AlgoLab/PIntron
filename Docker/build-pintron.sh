#!/bin/bash -x
####
# Copyright (C) 2015  Gianluca Della Vedova
#
# Distributed under the terms of the GNU Affero General Public License (AGPL)
#
#
# This file is part of PIntron.
#
# PIntron is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# PIntron is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with PIntron.  If not, see <http://www.gnu.org/licenses/>.
#
####

# This script can also be used to update pintron

VERSION=${pintron_version:master}

cd /home/pintron
test -d pintron/.git || git clone https://github.com/AlgoLab/PIntron.git pintron
cd pintron
git checkout $VERSION
git pull -r
make dist
