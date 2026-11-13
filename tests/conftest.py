import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import pytest

from app import main as main_module
from app.jobs import JobManager
from app.trim_jobs import TrimManager


@pytest.fixture(scope="module")
def client(tmp_path_factory):
    jobs_root = tmp_path_factory.mktemp("jobs")
    trim_root = tmp_path_factory.mktemp("trimming")
    # Reload the module-level manager paths so no test touches /data.
    main_module.JOBS_ROOT = jobs_root
    main_module.TRIM_ROOT = trim_root
    main_module.manager = JobManager(jobs_root)
    main_module.manager.wipe()
    main_module.trim_manager = TrimManager(trim_root, main_module.manager)
    main_module.trim_manager.wipe()
    from fastapi.testclient import TestClient

    with TestClient(main_module.app) as c:
        yield c
