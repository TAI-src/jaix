# Singular Environments

A singular environment implements the logic and values that are actually passed to the optimisation algorithm. Examples for singular environments are [`ECEnvironment`](/jaix/env/singular/ec_env.py) and [`MastermindEnvironment`](/jaix/env/singular/mastermind_env.py).

## Adding a Singular Environment
To add a new singular environment, create a class that inherits from the [SingularEnvironment](/jaix/env/singular/singular_environment.py) class. 

### Implementation
`SingularEnvironment` inherits from [`gymnasium.Env`](https://gymnasium.farama.org/api/env/) and thus should implement the same methods:
* `step`: Updates an environment with actions returning the next agent observation, the reward for taking that actions, if the environment has terminated or truncated due to the latest action and information from the environment about the step, i.e. metrics, debug info.
* `reset`: Resets the environment to an initial state, required before calling step. Returns the first agent observation for an episode and information, i.e. metrics, debug info.
* `render`: Renders the environments to help visualise what the agent see, examples modes are “human”, “rgb_array”, “ansi” for text.
* `close`: Closes the environment, important when external software is used, i.e. pygame for rendering, databases

### Config
The configuration for a singular environment `MySingularEnv` would look something like this:
```json
    "env_config": {
      "jaix.EnvironmentConfig": {
        "suite_class": "jaix.suite.Suite",
        "suite_config": {
          "jaix.suite.SuiteConfig": {
            "env_class": "full.path.to.MySingularEnv",
            "env_config": {
              "full.path.to.MySingularEnvConfig": {
                "property1": "val1",
                "property2": 0,
              }
            },
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

