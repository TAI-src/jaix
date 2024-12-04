import wandb
# 1. Start a W&B Run
wandb.init()

# 2. Save mode inputs and hyperparameters
wandb.config.learning_rate = 0.01

# 3. Log metrics over time to visualize performance
wandb.log({"loss": 5})

