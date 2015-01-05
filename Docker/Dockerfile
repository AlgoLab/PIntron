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
#
# Some parts are originally from baseimage-docker (https://github.com/phusion/baseimage-docker)
# Copyright (c) 2013-2014 Phusion
#
####
FROM phusion/baseimage

# Set correct environment variables.
ENV HOME /root

# Use baseimage-docker's init process.
CMD ["/sbin/my_init"]

# Update and install packages
COPY apt-get.sh /build/
RUN /build/apt-get.sh

# PIntron requisites
COPY pre-pintron.sh /build/
COPY config/* /root/config/
RUN /build/pre-pintron.sh

# Enable nginx
USER root
RUN mkdir /etc/service/nginx
COPY runit-nginx /etc/service/nginx/run
# Add configuration to nginx:
RUN rm /etc/nginx/sites-enabled/default
COPY config/nginx-pintron.conf /etc/nginx/sites-enabled/pintron.conf

# Periodically check if we can update pintron
RUN echo "27 * * * *	pintron test -x ~/pintron/Docker/build-pintron.sh &&  ~/pintron/Docker/build-pintron.sh" >> /etc/crontab


# Get and build PIntron
COPY build-pintron.sh /build/
USER pintron
ENV HOME /home/pintron
WORKDIR /home/pintron
RUN /build/build-pintron.sh

# Clean up APT when done.
USER root
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
EXPOSE 80 22
