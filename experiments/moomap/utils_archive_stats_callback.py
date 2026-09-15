from pymoo.core.callback import Callback
from jaix.env.utils.archive.mo_archive import MOArchive


class ArchiveStatsCallback(Callback):
    def __init__(self, archive: MOArchive):
        super().__init__()
        self.archive = archive
        self.data["archive_stats"] = []

    def notify(self, algorithm):
        archive_stats = self.archive.get_archive_stats()
        archive_stats["generation"] = len(self.data["archive_stats"]) + 1
        self.data["archive_stats"].append(archive_stats)
