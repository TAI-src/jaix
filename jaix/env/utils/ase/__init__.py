try:
    from jaix.env.utils.ase.ljclust_adapter import LJClustAdapter, LJClustAdapterConfig
except ImportError:
    # If the import fails, we set LJClustAdapter and LJClustAdapterConfig to None
    LJClustAdapter = None  # type: ignore[assignment,misc]
    LJClustAdapterConfig = None  # type: ignore[assignment,misc]
