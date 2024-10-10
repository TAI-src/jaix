import logging
import logging.config
import copy


# This is a test handler so it is easy to capture what was logged
class TestHandler(logging.Handler):
    def __init__(self, level=logging.NOTSET):
        self.record_log = []
        super().__init__(level)

    def emit(self, record):
        msg = self.format(record)
        self.last_record = record
        self.record_log.append(copy.copy(record))
        self.last_record.format_msg = msg


logging_config = {
    "version": 1,
    "disable_existing_loggers": True,
    "formatters": {
        "default": {
            "format": "%(asctime)s - %(name)s - %(filename)s - %(funcName)s - %(threadName)s - %(levelname)s - %(message)s"
        },
    },
    "handlers": {
        "console": {
            "()": TestHandler,
            "level": "DEBUG",
            "formatter": "default",
        },
        "root_console": {
            "level": "DEBUG",
            "class": "logging.StreamHandler",
            "formatter": "default",
        },
    },
    "loggers": {
        "DefaultLogger": {
            "level": "DEBUG",
            "handlers": ["console"],
        },
    },
    "root": {
        "level": "DEBUG",
        "handlers": ["root_console"],
    },
}

logging.config.dictConfig(logging_config)

logger = logging.getLogger("DefaultLogger")
test_handler = logger.handlers[0]
