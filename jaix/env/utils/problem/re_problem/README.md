# RE problems by Ryoji Tanabe and Hisao Ishibushi

We implement an adapter to be able to run the real-world problems (RE problems) by Ryoji Tanabe and Hisao Ishibushi. The original implementation is available at [https://github.com/ryojitanabe/reproblems].

To make it easier to install, we copy the python implementation and data from the original repository, that is also linked as a [submodule](/deps/reproblems/README.md) in this repository. The original implementation is licensed under the MIT license, and we also use the same license for our copy.

So the duplicated files are vendored in and should never be manually changed. Instead, if there are changes to the submodule, the files can be updated by running the following commands in the root of this repository:

```bash
git submodule update --remote --merge deps/reproblems # to update the submodule
python scripts/update_reproblems.py
```

To avoid missing out on bug fixes on the original repository, there is also a [GitHub action](../../../../../.github/workflows/re_problems.yml) that periodically runs the above commands and creates a pull request if there are any changes.

Note that the version in this repo does go through auto-linting and is thus not the exact same version as in the submodule.
