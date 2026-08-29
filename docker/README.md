# docker notes

## constant images

Dependencies are built into constant docker images that are not updatedthrough CI/CD. This is to save on unnecessary build time (some of them are expensive), but also to make sure that the builds are reproducible and stable. The images are built from this repository and pushed to DockerHub.

To modify the images, you can build them locally and push to DockerHub.

```
docker compose -f compose.build-deps.yml build
docker compose -f compose.build-deps.yml push
```

These dependency images are currently:

- kim, containing kimpy and the kim model for jaix_ase
- tabrepo, containing the tabrepo with preloaded data for jaix_hpo
- cobi, containing cobi package (since it is not on pypi) for jaix_cobi
