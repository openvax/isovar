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

"""
Library modules only create loggers. Commands configure output explicitly,
so importing Isovar never changes an embedding application's logging.
"""

import logging
import logging.config

CLI_LOGGER_NAMES = ("isovar", "varcode", "pyensembl", "datacache")


def get_logger(name):
    return logging.getLogger(name)


def configure_cli_logging(level="INFO"):
    """Send command progress from Isovar and its annotation libraries to stderr.

    The root logger and any logger not named in ``CLI_LOGGER_NAMES`` are left
    unchanged. Repeated calls replace, rather than duplicate, the handler.
    """
    logging.config.dictConfig({
        "version": 1,
        "disable_existing_loggers": False,
        "formatters": {"isovar": {
            "format": "%(asctime)s - %(name)s:%(lineno)s - %(levelname)s - %(message)s"}},
        "handlers": {"stderr": {"class": "logging.StreamHandler", "formatter": "isovar"}},
        "loggers": {name: {"level": level, "handlers": ["stderr"]} for name in CLI_LOGGER_NAMES},
    })
