# Docker images for CI

CI can be achieved with Docker containers.
This has the advantage to be self-contained and independent on the session where the CI is launched.

Two `Dockerfile`s are available for compilation with `gcc` and `clang`.

The images can be built and pushed to the GitLab server with the following commands (the example applies to clang6):

```bash
docker build -t gitlab.chab.ethz.ch:4567/scine/scine/clang6 .
docker push gitlab.chab.ethz.ch:4567/scine/scine/clang6
```

If you want to test the image before pushing it:
```bash
docker run -it gitlab.chab.ethz.ch:4567/scine/scine/clang6
```

Note that the docker service must be running for this to work.
You may need to execute the commands as `root`.
