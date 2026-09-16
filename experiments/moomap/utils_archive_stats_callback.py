from jaix.env.utils.archive.mo_archive import MOArchive
from pymoo.core.callback import Callback


class ArchiveStatsCallback(Callback):
    def __init__(self, archive: MOArchive, pop_size: int):
        super().__init__()
        self.archive = archive
        self.pop_size = pop_size
        self.data["archive_stats"] = []

    def notify(self, algorithm):
        archive_stats = self.archive.get_archive_stats()
        gen = len(self.data["archive_stats"]) + 1
        archive_stats["generation"] = gen
        archive_stats["fevals"] = gen * self.pop_size
        self.data["archive_stats"].append(archive_stats)
