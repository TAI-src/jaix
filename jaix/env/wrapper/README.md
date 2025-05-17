# Environment Wrappers

In order to ensure compatibility across different use cases, the OpenAI gym framework has introduced the concept of [environment wrappers](https://gymnasium.farama.org/api/wrappers/), which can alter a problem environment without needing to modify the underlying code. We release several convenience wrappers with this benchmark (e.g. to ensure early stopping and automatic logging).

 More importantly, wrappers offer a way to explore different problem models and evaluation settings without changing the underlying implementation (for example, by adding additional constraints to the search space, changing how solutions are represented, adding an additional objective, adding context information).


## Configuring a Wrapper

To configure a wrapper inside an [`EnvironmentConfig`](/experiments/config.md#environment-config), you need to specify the `env_wrappers` as a list of tuples, where each tuple contains the class of the wrapper and the configuration of the wrapper (as a dictionary). The class of the wrapper is a subclass of [`jaix.env.wrapper.Wrapper`](jaix/env/wrapper.py).

For example, to specify a [`LoggingWrapper`](/jaix/env/wrapper/logging_wrapper.py) and a [`MaxEvalWrapper`](/jaix/env/wrapper/max_eval_wrapper.py), you would do the following:
```json
        "env_wrappers": [
          [
            "jaix.env.wrapper.LoggingWrapper",
            {
              "logger_name": "logger name",
            }
          ],
          [
            "jaix.env.wrapper.MaxEvalWrapper",
            {
              "max_evals": 100,
            }
          ]
        ]
```

## Adding a Wrapper

To add a new environment wrapper, create a class that inherits from the [`PassthroughWrapper`](jaix/env/wrapper/passthrough_wrapper.py) class. The `PassthroughWrapper` guarantees that all added methods are still available in the wrapped environment (unless they are overwritten). This is the same behaviour that you would expect from inheritance.

### Implementation
`PassthroughWrapper` inherits directly from [`gymnasium.Wrapper`](https://gymnasium.farama.org/api/wrappers/) and can thus implement the same as a [`gymasium.Env`](https://gymnasium.farama.org/api/env/) class. In contrast to an environment, however, it does not need to implement all of them as it simply overwrites them.

### Config

The configuration for a wrapper `MyWrapper` would look something like this:
```json
    "env_config": {
      "jaix.EnvironmentConfig": {
        "suite_class": "jaix.suite.Suite", #depends on suite
        "suite_config": {
          "jaix.suite.SuiteConfig": {
          ...
        },
        "env_wrappers": [
          [
            "full.path.to.MyWrapper",
            {
              "property1": "val1",
              "property2": 0,
            }
          ]
        ]
      }
    },
```
For an example of a full config with a singular environment, see [`hpo`](/experiments/hpo/binary.json), and for more information on the config format, see the [module documentation](https://github.com/TAI-src/ttex/tree/main/ttex/config).

