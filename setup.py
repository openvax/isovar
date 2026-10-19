"""Ensure incremental wheel builds copy only current Isovar sources."""

from pathlib import Path
from shutil import rmtree

from setuptools import setup
from setuptools.command.build_py import build_py


class CleanBuildPy(build_py):
    """Discard staged package files before copying current sources."""

    def run(self):
        built_package = Path(self.build_lib) / "isovar"
        if built_package.is_symlink():
            built_package.unlink()
        elif built_package.exists():
            rmtree(built_package)
        super().run()


setup(cmdclass={"build_py": CleanBuildPy})
