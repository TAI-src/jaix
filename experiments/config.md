# Experiment Config


This repository is using the [`config`](https://github.com/TAI-src/ttex/tree/main/ttex/config) module provided by the ['tai-ttex`](https://pypi.org/project/tai-ttex/) package. For mode details on their usage, please check out the documentation in that repository.

Here, we detail a specific class, [`ExperimentConfig`](/jaix/experiment.py), that fully specifies a complete experiment. It is made up of several nested configurations, that specify the separate modules. For our reasoning of ensuring configurability in this repository, see [below](#motivation-configuration-driven-development).

## ExperimentConfig class

As defined in [`jaix/experiment.py`](/jaix/experiment.py), the `ExperimentConfig` class requires the following components, which are all described in subsections below:

* [`env_config`](#environment-config): `EnvironmentConfig`. Contains the full specification of the environment, i.e. the problem.
* `runner_class`: `Type[Runner]`. The class of the runner to be used. This is a subclass of [`jaix.runner.Runner`](jaix/runner.py).
* `runner_config`: `RunnerConfig`. The configuration of the runner, i.e. the algorithm paradigm.
* `opt_class`: `Type[Optimiser]`. The class of the optimiser to be used. This is a subclass of [`jaix.optimiser.Optimiser`](jaix/optimiser.py).
* `opt_config`: `OptimiserConfig`. The configuration of the optimiser, i.e. the optimisation algorithm.
* `logging_config`: `LoggingConfig`. The configuration of the logging, i.e. the logging settings.

### Environment Config

As defined in [`jaix/environment_factory.py`](/jaix/environment_factory.py), the `EnvironmentConfig` class requires the following components, which are all described below:

* `suite_class`: `Type[Suite]`. The class of the suite to be used. This is a subclass of [`jaix.suite.Suite`](jaix/suite.py).
* `suite_config`: `SuiteConfig`. The configuration of the suite, i.e. the problem generator. [Suite Configuration details](/jaix/suite/README.md).
* `env_wrappers`: `List[Tuple[Type[Wrapper], Dict]]`. The wrappers to be used. Each wrapper is a tuple of the class of the wrapper and the configuration of the wrapper. The class of the wrapper is a subclass of [`jaix.env.wrapper.Wrapper`](jaix/env/wrapper.py). [Wrapper Configuration details](/jaix/env/wrapper/README.md#wrapper-configuration).
* `comp_config`: `CompositeEnvironmentConfig`. The configuration of the composite environment, i.e. the environment that is used to switch between different environments. This is a subclass of [`jaix.env.composite.CompositeEnvironment`](jaix/env/composite/composite_environment.py). The required keys depend on the environment config. Which instances are passed to the `CompositeEnvironment` is determined by the `instances` and `agg_instances` keys in the [suite configuration](/jaix/suite/README.md).
* `seed`: `int`. The seed to be used for the environment. This is used to ensure that the experiment is reproducible.

### Runner Config

The configuration of the runner completely depends on the runner class. The runner class is a subclass of [`jaix.runner.Runner`](jaix/runner.py). The configuration is passed to the constructor of the runner class. The runner class is responsible for defining the interface between the optimiser and the environment. It is responsible for running the optimisation loop and calling the optimiser to optimise the environment.

An example of a runner class is [`ATRunner`](/jaix/runner/ask_tell/ask_tell_runner.py), which implements the ask-tell interface. Full configuration examples are available in the [experiments](/experiments/README.md) folder.

### Optimiser Config

As with the runner, the configuration of the optimiser completely depends on the optimiser class. The optimiser class is a subclass of [`jaix.runner.Optimiser`](jaix/runner/optimiser.py). The configuration is passed to the constructor of the optimiser class. The optimiser class is responsible for implementing the optimisation logic within the interface defined by the runner. 

An example of an optimiser class is [`ATOptimers`](/jaix/runner/ask_tell/at_optimiser.py). Full configuration examples are available in the [experiments](/experiments/README.md) folder.

### Logging Config

Basic logging capabilities are always added to the [`experiment`](/jaix/experiment.py) to ensure basic tracking. However, this can be modified by using the [python logging facility](https://docs.python.org/3/library/logging.html), and specifically the [logging configuration](https://docs.python.org/3/library/logging.config.html) this framework offers.

So, to add a `logging.Logger`, specify the following in the `jaix.LoggingConfig`:

* `log_level`: `int`. The logging level to be used. This is used to filter the log messages. The default is `logging.INFO`.
* `logger_name`: `str`. The name of the logger to be used. This is used to identify the logger in the log messages. The default is `DefaultLogger`. The default name can be changed by specifying the global variable `LOGGER_NAME`.
* `disable_existing`: `bool`. If `True`, all existing loggers are disabled. This is used to ensure that only the loggers defined in the configuration are used. The default is `True`.
* `dict_config`: `Dict`. The configuration of the logger. This is used to define the logging format, handlers, and other settings. More details can be found in the [python logging documentation](https://docs.python.org/3/library/logging.config.html). If nothing is specified, the [default logging configuration](https://docs.python.org/3/library/logging.config.html) from the `tai-ttex` package is used.

For example, the [experiment launcher](/jaix/utils/launch_experiment.py) adds a specific `logging.Handler` that handles the integration to wandb (an experiment tracking framework). This handler is implemented in [`ttex.log.handler`](https://github.com/TAI-src/ttex/blob/main/ttex/log/handler/wandb_handler.py). This is done simply by adding the appropriate information to the logging config, as well as an additional `LoggingWrapper`.


## Motivation Configuration-driven development

Configuration-driven development is a software development approach that emphasizes the use of configuration files to define the behavior and settings of an application. This approach allows developers to separate the application's logic from its configuration, making it easier to manage and modify settings without changing the codebase. Configuration-driven development is particularly useful in scenarios where applications need to be deployed in different environments or require frequent updates to their settings.

In the context of benchmarking, configuration-driven development can be applied to define the parameters and settings for benchmarking tests. By using configuration files, developers can easily adjust the algorithm and problem hyper-parameters. This allows for an easy overview and tracking of previous experiments and ensures complete repeatability, as no variables are hard-coded in the codebase.



