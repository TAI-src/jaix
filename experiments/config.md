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
* `suite_config`: `SuiteConfig`. The configuration of the suite, i.e. the problem generator.
* `env_wrappers`: `List[Tuple[Type[Wrapper], Dict]]`. The wrappers to be used. Each wrapper is a tuple of the class of the wrapper and the configuration of the wrapper. The class of the wrapper is a subclass of [`jaix.env.wrapper.Wrapper`](jaix/env/wrapper.py).
* `comp_config`: `CompositeEnvironmentConfig`. The configuration of the composite environment, i.e. the environment that is used to switch between different environments. This is a subclass of [`jaix.env.composite.CompositeEnvironment`](jaix/env/composite/composite_environment.py).
* `seed`: `int`. The seed to be used for the environment. This is used to ensure that the experiment is reproducible.




## Motivation Configuration-driven development

Configuration-driven development is a software development approach that emphasizes the use of configuration files to define the behavior and settings of an application. This approach allows developers to separate the application's logic from its configuration, making it easier to manage and modify settings without changing the codebase. Configuration-driven development is particularly useful in scenarios where applications need to be deployed in different environments or require frequent updates to their settings.

In the context of benchmarking, configuration-driven development can be applied to define the parameters and settings for benchmarking tests. By using configuration files, developers can easily adjust the algorithm and problem hyper-parameters. This allows for an easy overview and tracking of previous experiments and ensures complete repeatability, as no variables are hard-coded in the codebase.



