# Experiments

This folder contains the full configuration to several different experiments run for the paper, as well as a few additional ones. Experiments using different singular environments are divided into different subfolders.

* [brachytherapy](rbf): Experiments using the brachytherapy environment
* [game tree](mmind): Experiments using the game tree environment
* [ensemble selection](hpo): Experiments using the ensemble selection environment
* [coco](coco): Experiments using the COCO environment

## Run experiment locally (docker)

Requirements:
* docker
* docker-compose

```bash
git clone git@github.com:TAI-src/jaix.git
cd docker
docker-compose pull
```
Create an `.env` file in the `docker` folder. If you want to use wandb for logging, you need to add the following information to the file.
```
WANDB_API_KEY=<your_wandb_api_key> #from https://wandb.ai/authorize
WANDB_ENTITY=<your_wandb_entity>
```
Start experiment:
```bash
docker run -d jaix bash -c "pip install tai_jaix; python /jaix/utils/launch_experiment.py --config_file /experiments/<path/to/config_file>"
```
This will start a docker container in the background that runs the experiment. The config file needs to be available inside the docker container, so it is easiest to use the experiments folder, or otherwise mount the folder containing the config file into the docker container.

Here is an overview of the available command line options
* `--config_file` (required): path to the config file
* `--project` (optional): name of the wandb project to log to. No logging is done if this is not provided.
* `--group` (optional): name of the wandb group to log to.
* `--repeat` (optional): number of times to repeat the experiment. Default is 1.
* `--sweep_keys` (optional): space-separated list of keys that specify which specific value in the configuration should be modified. See example below.
* `--sweep_values` (optional): space-separated list of values that should be used for the keys specified in `--sweep_keys`. See example below.

### Sweeping

A sweep allows you to run an experiment based on one config file, but vary the values for exactly one parameter. If you run the experiment with wandb active (you passed a project), and do not pass a group, the group name will be automatically generated based on sweep values.

Let's assume the following config:
```json
{
  "config_class": {
    "key1": "value1",
    "key2": "value2",
    "config_class2": "config2",
    "config2": {
      "key3": "value3"
    }
  }
}
```
To try this config with values `a,b,c` for key3, you would specify:
* `--sweep_keys config_class config2 key3`
* `--sweep_values a b c`

The runs would then be logged under three different groups called `key3 a`, `key3 b` and `key3 c`.

[!NOTE]
It is currently only possible to sweep over a single key. This key cannot be part of a list. The sweep values must be numeric

## Run experiment locally (without docker)
This is also possible, but untested and requires manual installation of all required dependencies. Check the Dockerfiles in the [docker folder](../docker) for the required dependencies.

## Run experiments via wandb

Requirements:
* docker
* docker-compose
* wandb


### wandb Configuration
* Set up a queue in wandb with name `$queue` ([wandb_docs](https://docs.wandb.ai/guides/launch/walkthrough#create-a-queue))
* Create a job to start your experiment. You can technically use any type of job (git, code artifact or image), To create, run the following command in an environment with wandb installed (e.g. ttex docker container)
```bash
wandb job create -p $project -e $entity -n $job_name -a $job_alias image $image
```
### Start listener
```bash
git clone git@github.com:TAI-src/jaix.git
cd docker
docker-compose pull
```
Create an `.env` file in the `docker` folder with the following content.```
```
WANDB_API_KEY=<your_wandb_api_key> #from https://wandb.ai/authorize
WANDB_ENTITY=$entity
WANDB_Q=$queue
DOCKER_USER_NAME=<your_docker_user_name>
DOCKER_PWD=<your_docker_password>
```

To start the listener, run the following command:
```bash
docker compose up -d --remove-orphans wandb_launcher
```
### Launch a new job

Now the job can be launched from either the webpage or command line. The full config needs to be passed via the key `run_config`

#### Via the webpage

Navigate to the $project page and click on "Jobs". Find the job with the selected $name and click "Launch". Launch the job with defaults or edit any of the suggested values.

#### Via the command line
[!NOTE]
:warning: This does currently not actually work, there seems to be something missing in wandb (passing the alias)

```
wandb launch -j $job_name -q $queue -e $entity -c path/to/config.json
```

## Reproduce experiment from wandb

To reproduce the data from a specific run on wandb, you just need to identify the run config and run it as described above, either [locally](#run-experiment-locally-docker) or via [wandb](#run-experiments-via-wandb). The full config for a given run can be found under `files/config.yaml`
