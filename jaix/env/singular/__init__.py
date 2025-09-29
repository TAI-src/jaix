# only import if tabrepo extra options
try:
    from jaix.env.singular.hpo_env import HPOEnvironmentConfig, HPOEnvironment
except ImportError:
    # If the import fails, we set HPOEnvironmentConfig and HPOEnvironment to None
    HPOEnvironmentConfig = None  # type: ignore[assignment,misc]
    HPOEnvironment = None  # type: ignore[assignment,misc]
# only import if ase extra options
try:
    from jaix.env.singular.ljclust_env import (
        LJClustEnvironmentConfig,
        LJClustEnvironment,
    )
except ImportError:
    # If the import fails, we set LJClustEnvironmentConfig and LJClustEnvironment to None
    LJClustEnvironmentConfig = None  # type: ignore[assignment,misc]
    LJClustEnvironment = None  # type: ignore[assignment,misc]
