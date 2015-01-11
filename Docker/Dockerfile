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
USER root
# Set correct environment variables.
ENV HOME /root
ENV LC_ALL C
ENV DEBIAN_FRONTEND noninteractive

# Use baseimage-docker's init process.
CMD ["/sbin/my_init"]

#PIntron requisites
COPY config/root_key.pub /root/.ssh/authorized_keys
COPY config/* /root/config/
RUN apt-get update
RUN apt-get -y dist-upgrade
RUN apt-get install -y git-core build-essential python3 ssh gengetopt

# Create the PIntron user
RUN adduser --disabled-password --gecos "PIntron" pintron
RUN usermod -L pintron
RUN mkdir -p /home/pintron/.ssh
RUN chmod 700 /home/pintron/.ssh
RUN mkdir -p /home/pintron/input
RUN mkdir -p /home/pintron/results
# Get the web servers public keys
RUN test -f /root/config/web_keys.pub && cp /root/config/web_keys.pub /home/pintron/.ssh/authorized_keys || echo "No frontends configured"
RUN test -f config/web_keys_url && curl -L $(cat config/web_keys_url) >>  /home/pintron/.ssh/authorized_keys || echo "No frontends configured"
RUN chown --recursive pintron:pintron /home/pintron/


# Periodically check if we can update pintron
RUN echo "27 * * * *	pintron test -x ~/pintron/Docker/build-pintron.sh &&  ~/pintron/Docker/build-pintron.sh" >> /etc/crontab

RUN echo "unattended-upgrades	unattended-upgrades/enable_auto_updates	boolean	true" > /tmp/debconf.tmp
RUN debconf-set-selections /tmp/debconf.tmp
RUN rm -f /tmp/debconf.tmp

# Get and build PIntron
COPY build-pintron.sh /build/
USER pintron
ENV HOME /home/pintron
WORKDIR /home/pintron
RUN /build/build-pintron.sh

USER root
ENV HOME /root
EXPOSE 22
