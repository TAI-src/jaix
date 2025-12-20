# docker notes

## constant images

There are some expensive images to build that can be assumed to be stable. These are not rebuilt with CI/CD to save on testing time.
These are:
- tabrepo
- kim

If they do need to be rebuilt
```
DOCKER_BUILDKIT=1 docker build --shm-size=16g \
  -t taisrc/tabrepo \
  -f ../docker/Dockerfile_tabrepo \
  ../deps/
```

```
DOCKER_BUILDKIT=1 docker build \
  -t taisrc/kim \
  -f ../docker/Dockerfile_kim \
  ../deps/
```
