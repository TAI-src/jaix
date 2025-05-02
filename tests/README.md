# Tests

Tests are pytest-based and can be run using the following command:

```bash
cd tests
python -m pytest -v
```
Or inside a Docker container:

```bash
cd docker
docker compose run jaix bash
python -m pytest -v
```

The separate folders contain unit tests for each of the modules. Integration tests are in the top level folder.


