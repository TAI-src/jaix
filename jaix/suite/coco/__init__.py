try:
    from jaix.suite.coco.coco_problem import COCOProblem
except ImportError:
    # If the import fails, we set COCOProblem to None
    COCOProblem = None  # type: ignore[assignment,misc]
try:
    from jaix.suite.coco.coco_suite import COCOSuiteConfig, COCOSuite
except ImportError:
    # If the import fails, we set COCOSuiteConfig and COCOSuite to None
    COCOSuiteConfig = None  # type: ignore[assignment,misc]
    COCOSuite = None  # type: ignore[assignment,misc]
