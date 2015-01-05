# Design

The Dockerfile builds an image that contains:
* ssh server, to feed input data to pintron
* nginx, to get the results

Therefore ports 80 and 22 must be open.

PIntron is built inside the /home/app/pintron directory.
Input files are placed into the /home/app/input dir, while the results are
stored into /home/app/results.

# Customization
