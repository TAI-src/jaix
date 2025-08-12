try:
    from jaix.suite.coco.coco_problem import COCOProblem
except ImportError:
    # If the import fails, we set COCOProblem to None
    COCOProblem = None
try:
    from jaix.suite.coco.coco_suite import COCOSuiteConfig, COCOSuite
except ImportError:
    # If the import fails, we set COCOSuiteConfig and COCOSuite to None
    COCOSuiteConfig = None
    COCOSuite = None
