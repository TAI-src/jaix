# Game tree exploration

Corresponding singular environment [MastermindEnvironment](../../jaix/env/singular/mastermind_env.py) and composite environment [SwitchingEnvironment](../../jaix/env/composite/switching_environment.py) where relevant.

## Telltale

For this example, decisions have to be made in the correct order (see paper). Values are based on [Telltale's The Walking Dead Season 1](https://venturebeat.com/games/the-walking-dead-season-one-plot-graph/).

### Configuration Details
* [singular config](./telltale.json)
* [composite config](./telltale_comp.json)

Important details:
* `num_slots_range`: The number of decisions range between 25 and 35, depending on the num_slots_range
* `num_colours_range`: The number of colours (decision options) is always 2, based on the inspiration example.
* `update_opts.s`: Is the success ratio that we vary during the exeriments

## Mastermind

This is not part of the paper, but we ran additional examples based on the game [Mastermind](https://en.wikipedia.org/wiki/Mastermind_(board_game). In this case, the goal is to play the game as a player. The values are based on the standard version of the game.

### Configuration Details
* [singular config](./mmind.json)
* [composite config](./mmind_comp.json)

Important details:
* `num_slots_range`: The number of slots (decisions) is always 4 
* `num_colours_range`: The number of colours (decision options) is always 6
* `update_opts.s`: Is the success ratio that we vary during the exeriments
