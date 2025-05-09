# Static Problems

A static problem is a function that maps an input to an output, and is not influenced by time or the state. For an example of existing static problems, see the [`Sphere`](/jaix/env/utils/problems/sphere.py).


## Adding a Static Problem
To add a new static problem (standard in EC research), create a class that inherits from the [`StaticProblem`](/jaix/env/utils/problem/static_problem.py) class. 

### Implementation
This class should implement the following methods:
* `_eval`: This method should evaluate the problem and return a tuple of (fitness value, reward). For single_objective problems, the reward can just be the same as the fitness value. For multi-objective problem, the reward could be e.g. the hyper-volume. You can also just return a non-meaningful value if your algorithm (as most EC algorithms do) do not use rewards.
* `__str__`: A string representation of the problem. This is used for logging and debugging purposes.

Experiments are then runnable by using the [`ECEnvironment`](/jaix/env/singular/ec_env.py). This class is a wrapper around the static problem that offers several convenience functionalities by automatically determining the number of objectives, the number of dimensions, and the bounds of the problem. It also implements the required methods that make it compatible with the OpenAI gym framework.

### Config
A suitable environment configuration for `MyProblem` would look as follows:
```json
   "env_config": {
      "jaix.EnvironmentConfig": {
        "suite_class": "jaix.suite.ECSuite",
        "suite_config": {
          "jaix.suite.ECSuiteConfig": {
            "func_classes": ["full.path.to.MyProblem"],
            "func_configs": [
              {
                "full.path.to.MyProblem.Config": {
                  "property1": "val",
                  "property2": 0
                }
              }
            ],
            "env_config": {
              "jaix.env.singular.ECEnvironmentConfig": {
                "budget_multiplier": 1000
              }
            },
          }
        },
      }
    }
```
For an example of a full config with a static problem, see [`rbf`](/experiments/rbf/README.md), and for more information on the config format, see [config instructions](/experiments/config.md).

