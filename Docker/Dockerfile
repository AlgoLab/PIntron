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

# Update the system
RUN apt-get update
RUN apt-get -y dist-upgrade
RUN apt-get install -y git-core build-essential python3 ssh gengetopt
RUN echo "unattended-upgrades	unattended-upgrades/enable_auto_updates	boolean	true" > /tmp/debconf.tmp
RUN debconf-set-selections /tmp/debconf.tmp
RUN rm -f /tmp/debconf.tmp

# Set the ssh key for root
COPY config/root_key.pub /root/.ssh/authorized_keys

# Create the PIntron user
RUN adduser --disabled-password --gecos "PIntron" pintron
RUN usermod -L pintron
RUN mkdir /config
RUN chmod 777 /config
RUN chown pintron:pintron /config

COPY build-pintron.sh /build/

USER pintron
ENV HOME /home/pintron
WORKDIR /home/pintron
RUN ssh-keygen -q -f /home/pintron/.ssh/id_rsa -N ""
RUN mkdir -p /home/pintron/input
RUN mkdir -p /home/pintron/results
# Build pintron and get the web servers public keys
RUN /build/build-pintron.sh
# Periodically check if we can update pintron
RUN echo 'MAILTO=""' | crontab -
RUN echo '27 * * * * test -x ~/pintron/Docker/build-pintron.sh && ~/pintron/Docker/build-pintron.sh' | crontab -
RUN crontab -l
COPY config/web_keys.pub /config/authorized_keys
RUN ln -s /config/web_keys.pub /home/pintron/.ssh/authorized_keys

USER root
ENV HOME /root
EXPOSE 22
