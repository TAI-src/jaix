from pymoo.core.callback import Callback
from jaix.env.utils.archive.mo_archive import MOArchive


class ArchiveStatsCallback(Callback):
    def __init__(self, archive: MOArchive):
        super().__init__()
        self.archive = archive
        self.data["archive_stats"] = []

    def notify(self, algorithm):
        self.data["archive_stats"].append(self.archive.get_archive_stats())
