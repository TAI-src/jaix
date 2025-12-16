# docker notes

tabrepo building is currently not done automatically to prevent long runtimes during tests as it is assumed to be stable.

If it does need to be rebuilt:
```
DOCKER_BUILDKIT=1 docker build --shm-size=16g \
  -t taisrc/tabrepo \
  -f ../docker/Dockerfile_tabrepo \
  ../deps/
``
