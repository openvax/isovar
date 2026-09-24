# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import subprocess
import sys


def test_import_isovar_does_not_use_pkg_resources():
    result = subprocess.run(
        [
            sys.executable,
            "-W",
            "error:pkg_resources is deprecated as an API:UserWarning",
            "-c",
            (
                "import sys\n"
                "import isovar\n"
                "assert 'pkg_resources' not in sys.modules\n"
            ),
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr


def test_import_isovar_leaves_host_logging_unchanged():
    result = subprocess.run(
        [
            sys.executable,
            "-c",
            (
                "import logging\n"
                "host = logging.getLogger('host.app')\n"
                "import isovar, isovar.cli\n"
                "assert not host.disabled\n"
                "assert not logging.getLogger().handlers\n"
                "assert not logging.getLogger('isovar').handlers\n"
            ),
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr


def test_cli_logging_goes_to_stderr_once_and_keeps_existing_loggers():
    result = subprocess.run(
        [
            sys.executable,
            "-c",
            (
                "import logging\n"
                "host = logging.getLogger('host.app')\n"
                "from isovar.logging import configure_cli_logging\n"
                "configure_cli_logging()\n"
                "configure_cli_logging()\n"
                "assert not host.disabled\n"
                "assert len(logging.getLogger('isovar').handlers) == 1\n"
                "logging.getLogger('isovar.example').info('progress')\n"
                "logging.getLogger('isovar.example').debug('detail')\n"
            ),
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout == ""
    assert result.stderr.count("progress") == 1
    assert "detail" not in result.stderr
