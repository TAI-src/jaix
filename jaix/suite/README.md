# Suite

In this framework, the suite is responsible for providing the environments /problems. It should inherit from the [`jaix.suite.Suite`](jaix/suite.py) class and implement the following methods:
 * `_get_env`: Returns the environment for a given function and instance.

There are currently three main ways of configuring a suite. They are listed below and described in more detail in subsections.
* [`jaix.suite.ECSuite`](jaix/suite/ec_suite.py): Use if you want a suite comprised of [static problems](jaix/env/utils/problem/static_problem.py)
)
* `jaix.suite.Suite`: Use if you want a suite of [singular environments](jaix/env/singular/README.md).
* Custom: Use if you want complete freedom to implement your own suite. In this case, you should implement a class that inherits from the [`jaix.suite.Suite`](jaix/suite.py) class and implements the `_get_env` method. For an example, see the [`COCO Suite`](jaix/suite/coco_suite.py). The configuration will depend on your implementation. See the [config module](https://github.com/TAI-src/ttex/tree/main/ttex/config) for general instructions.

It is important to note that, independent of the way chosen, only provides (lists of) singular environments. How these environments are combined is described in the [`CompositeEnvironment`](jaix/env/composite/README.md) class and configured by the [`CompositeEnvironmentConfig`](jaix/env/composite/README.md) class.

## Configuring an ECSuite

To configure an `ECSuite` inside an [`EnvironmentConfig`](/experiments/config.md#environment-config), you need to specify the `suite_class` as `jaix.suite.ECSuite` and the `suite_config` as a dictionary with the following keys:
```json
{
          "jaix.suite.ECSuiteConfig": {
            "func_classes": [
              "full.path.to.StaticProblem"
            ],
            "func_configs": [
              {
                "full.path.to.StaticProblemConfig": {
                  "property1": "val1", # Depends on problem
                }
              }
            ],
            "env_config": {
              "jaix.env.singular.ECEnvironmentConfig": {
                "budget_multiplier": 1000
              }
            },
            "agg_instances": 0
          }
        }
        }
```


