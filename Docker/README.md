# Deployment
Simply build the image with
```cd Docker && docker build -t pintron .```
and run a container (the equivalent of starting a daemon with
```docker run -d -P --name "pintron-1" pintron```
# Design

The Dockerfile builds an image that contains:
* ssh server, to feed input data to pintron
* nginx, to get the results

Therefore ports 80 and 22 must be open.

PIntron is built inside the /home/app/pintron directory.
Input files are placed into the /home/app/input dir, while the results are
stored into /home/app/results.

If you want to map those directories to some specific locations on your server,
for instance because the input data are large, you can mount them as external
shares when running the image, with something like
```docker run -d -P --name "pintron-1" -v /media/HUGE_HARD_DISK:/home/app/input
pintron```. See
[Managing Data in Containers](http://docs.docker.com/userguide/dockervolumes/#volume-def)
for more information.

# Customization

If you want to track a version of PIntron different from stable (master), put
the desired version into ```config/pintron_version```.

Put the ssh public key(s) of the root user into ```config/root_key.pub```.

Put the URL containing the public keys of the webservers into
```config/web_keys_url``` or put the keys directly into
```config/web_keys.pub```.
