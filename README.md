# jaix Framework for Jacked-up Artificial Intelligence eXperiments

The jaix framework is a toolkit for running optimisation experiments based on the OpenAI gym framework. It's main goal is versatility, in terms of the possible experimental setups as well as the applicable algorithmic approaches. Published on [PyPI](https://pypi.org/project/tai-jaix/).

## Running an experiment

Experiments are fully described as a single [experiment configuration](/experiments/config.md) file. Examples, requirements and more details can be found in in the [experiments](experiments/README.md) folder. The required setup and instructions are detailed there, including different options (using Docker, local python, and launching via the wandb.ai web UI). All essentially boil down to a single command that starts the [experiment launcher](jaix/utiils/launch_experiment.py) with the desired config file.

```
pip install -e .
python jaix/utils/launch_experiment.py --config_file experiments/<path/to/config_file>"
```

You can either run one of the existing configurations in the [experiments](experiments/README.md) folder, or create your own. For full instructions, follow [experiments](experiments/config.md).

## Extending the framework

![modules](https://github.com/user-attachments/assets/aa328c45-9557-4462-aef1-0f7e3e1ac13b)

### Making additions configurable

This framework is using [configuration-driven development](/experiments/config.md#motivation-configuration-driven-development). This means that all experiments should be fully configurable in a json file. This allows easier comparison, overview, repeatability and tracking. This is implemented using the [`config`](https://github.com/TAI-src/ttex/tree/main/ttex/config) module provided by the ['tai-ttex`](https://pypi.org/project/tai-ttex/) package.

To maintain this ability, all new classes that have configurable properties should inherit from [`ConfigurableObject`](https://github.com/TAI-src/ttex/blob/main/ttex/config/configurable_object.py) and expect a config object that inherits from [`Config`](https://github.com/TAI-src/ttex/blob/main/ttex/config/config.py). See the module documentation for more details and examples.

### Static Problem (EC)

To add a new static problem (standard in EC research), create a class that inherits from the [`StaticProblem`](jaix/env/utils/problem/static_problem.py) class. A static problem is a function that maps an input to an output, and is not influenced by time or the state. For detailed instructions and examples, see [Static Problems](jaix/env/utils/problem/README.md).

Note that a static problem can also be implemented as [`SingularEnvironment`](jaix/env/singular/singular_environment.py) (see below). The [`StaticProblem`](jaix/env/utils/problem/static_problem.py) class and [`ECEnvironment`](jaix/env/singular/ec_env.py) class are introduced to allow the addition of new static problems without the need to implement the full environment logic. At the same time, this ensures compatibility with the OpenAI gym framework while being able to maintain standard coding interfaces for functions in EC research.

### Problem Environment (RL/EC)

First decide if you are implementing a [`CompositeEnvironment`](jaix/env/composite/composite_environment.py), [SingularEnvironment](jaix/env/singular/singular_environment.py), or potentially both.

* A composite environment is comprised of multiple singular once and should only be responsible for initialising the corresponding singular environments and switching between them.
* A singular environment implements the logic and values that are actually passed to the optimisation algorithm.

To add a new singular environment, create a class that inherits from the [`SingularEnvironment`](jaix/env/singular/singular_environment.py) class. For detailed instructions and examples, see [Singular Environments](jaix/env/singular/README.md). To add a new composite environment, create a class that inherits from the [`CompositeEnvironment`](jaix/env/composite/composite_environment.py) class. For detailed instructions and examples, see [Composite Environments](jaix/env/composite/README.md).

### Environment Wrapper

In order to ensure compatibility across different use cases, the OpenAI gym framework has introduced the concept of [environment wrappers](https://gymnasium.farama.org/api/wrappers/), which can alter a problem environment without needing to modify the underlying code. More importantly, wrappers offer a way to explore different problem models and evaluation settings without changing the underlying implementation (for example, by adding additional constraints to the search space, changing how solutions are represented, adding an additional objective, adding context information).

To add a new environment wrapper, create a class that inherits from the [`PassthroughWrapper`](jaix/env/wrapper/passthrough_wrapper.py) class. For detailed instructions and examples, see [Environment Wrappers](jaix/env/wrapper/README.md).

## Troubleshooting

Make sure all sub-modules for required dependencies in [`deps`](/deps) are pulled. To do that, navigate to each folder and do:

```
git submodule init
git submodule update
```

### Initialisation Order Issues

In case the wrappers need to pass information to each other, it is important to ensure that the initialisation order is correct. This can be done by specifying the order in which the wrappers are applied in the experiment configuration file. Then wrappers are applied in order, that means that the first wrapper in the list is the innermost wrapper.

The easiest way to pass information is via global variables defined in `jaix.utils.globals`. However, this is not ideal and should be avoided if possible.

It can also be helpful to consider the order in how the experiment setup happens:

1. If passed as dict, the [`ExperimentConfig`](jaix/experiment.py) object is created. That means that all nested Config objects are initialised (depth-first recursive strategy) as defined using [`ttex.config.ConfigFactory`](https://github.com/TAI-src/ttex/blob/main/ttex/config/config.py). This happens either in [launch_experiment.py](jaix/utils/launch_experiment.py) or externally.
2. The experiment is started with the passed config object by calling [`Experiment.run()`](jaix/experiment.py).
3. Setup is called on the `ExperimentConfig` object, which means that the following is called
    a. [`logging_config.setup()`](jaix/experiment.py) which creates the main logger and populates `jaix.utils.globals.LOGGER_NAME`
    b. [`env_config.setup()`](jaix/environment_factory.py) which also calls `setup()` on all wrappers in order.
    c. [`runner_config.setup()`](jaix/runner/runner.py)
    d. [`opt_config.setup()`](jaix/runner/optimiser.py)
    e. Iff a [`WandWrapper`](jaix/env/wrapper/wandb_wrapper.py) was passed, wandb is initialised using [ttex.log](https://github.com/TAI-src/ttex/blob/main/ttex/log/utils/wandb_logging_setup.py). This also means that a wandb run object is created.
4. The experiment id is identified (passed parameter, wandb run id, or newly created uuid) and set as a global variable.
5. The objects configured by the config are initialised (order depending on runner)
