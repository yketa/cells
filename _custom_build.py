"""
This module overrides setuptools build.
https://stackoverflow.com/questions/74409966/how-to-replace-setup-py-with-a-pyproject-toml-for-a-native-c-build-dependency
"""

from distutils.command.build import build as _build
from subprocess import check_call

class build(_build):
    def run(self):
        check_call("make", cwd="src/cells") # compile library
        super().run()                       # running build

