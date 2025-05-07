# jaix Framework for Jacked-up Artificial Intelligence eXperiments

The jaix framework is a toolkit for running optimisation experiments based on the OpenAI gym framework. It's main goal is versatility, in terms of the possible experimental setups as well as the applicable algorithmic approaches.

## Running an experiment


Experiments are fully described as a configuration file. Examples, requirements and more details can be found in in the [experiments](experiments/README.md) folder. The required setup and instructions are detailed there, including different options (using Docker, local python, via wandb.ai). All essentially boil down to a single command that starts the [experiment launcher](jaix/utiils/launch_experiment.py) with the desired config file.

```
pip install -e .
python jaix/utils/launch_experiment.py --config_file experiments/<path/to/config_file>"
```

You can either run one of the existing configurations in the [experiments](experiments/README.md) folder, or create your own. For full instructions, follow [experiments](experiments/config.md).

`



## Troubleshooting

Make sure all submodules for required dependencies in `deps` are pulled. To do that, navigate to each folder and do:
```
git submodule init
git submodule update
```

See 'test_launch_experiment' for example config
