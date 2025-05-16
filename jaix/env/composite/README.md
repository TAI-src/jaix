# Composite Environments

A composite environment is comprised of multiple singular ones and should only be responsible for initialising the corresponding singular environments and switching between them. An example for a composite environment is the [`SwitchingEnvironment`](jaix/env/composite/switching_environment.py).
 

## Adding a Composite Environment

To add a new composite environment, create a class that inherits from the [`CompositeEnvironment`](/jaix/env/composite/composite_environment.py) class.

### Implementation 

`CompositeEnvironment` inherits from [`gymnasium.Env`](https://gymnasium.farama.org/api/env/) and thus should implement the same methods:
* `step`: Updates an environment with actions returning the next agent observation, the reward for taking that actions, if the environment has terminated or truncated due to the latest action and information from the environment about the step, i.e. metrics, debug info.
* `reset`: Resets the environment to an initial state, required before calling step. Returns the first agent observation for an episode and information, i.e. metrics, debug info.
* `render`: Renders the environments to help visualise what the agent see, examples modes are “human”, “rgb_array”, “ansi” for text.
* `close`: Closes the environment, important when external software is used, i.e. pygame for rendering, databases

### Config

The configuration for a composite environment `MyCompositeEnv` would look something like this:
```json
    "env_config": {
      "jaix.EnvironmentConfig": {
        "suite_class": "jaix.suite.CompositeSuite",
        "suite_config": {
          "jaix.suite.CompositeSuiteConfig": {
            "env_classes": [
              "full.path.to.MySingularEnv1",
              "full.path.to.MySingularEnv2"
            ],
            "env_configs": [
              {
                "full.path.to.MySingularEnv1Config": {
                  "property1": "val1",
                  "property2": 0,
                }
              },
              {
                "full.path.to.MySingularEnv2Config": {
                  "property1": "val2",
                  "property2": 0,
                }
              }
            ],
            "agg_instances": 0,
            "functions": [
              1
            ]
          }
        },
        "env_wrappers": [
          [
            "jaix.env.wrapper.AnyFitWrapper",
            {}
          ]
        ]
      }
    },
```

For an example of a full config with a composite suite, see [`hpo`](/experiments/hpo/binary_comp.json), and for more information on the config format, see the [module documentation](https://github.com/TAI-src/ttex/tree/main/ttex/config).
