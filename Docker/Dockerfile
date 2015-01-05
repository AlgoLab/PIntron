FROM phusion/baseimage

# Set correct environment variables.
ENV HOME /root

# Use baseimage-docker's init process.
CMD ["/sbin/my_init"]

# PIntron requisites
COPY pre-pintron.sh /build/
RUN /build/pre-pintron.sh
# PIntron optimized build


# Enable nginx
RUN rm -f /etc/service/nginx/down

WORKDIR /home/app/pintron

# Get and build PIntron
COPY build-pintron.sh /build/
USER pintron
RUN /build/build-pintron.sh
WORKDIR /home/app/pintron

# Add configuration to nginx:
USER root
RUN rm /etc/nginx/sites-enabled/default
COPY config/nginx-pintron.conf /etc/nginx/sites-enabled/pintron.conf
# TODO: better configuration

# Clean up APT when done.
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

EXPOSE 80 22
