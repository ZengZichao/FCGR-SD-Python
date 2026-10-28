"""pytest 共享配置。"""

from __future__ import annotations

import shutil
import tempfile
from pathlib import Path

import pytest


@pytest.fixture
def workdir():
    """在工作区内创建临时目录，退出时删除。

    部分受限环境（沙箱）无法写系统临时目录，pytest 自带的 ``tmp_path``
    会因 PermissionError 失败；把临时目录放在 ``tests/_tmp`` 下则任何环境可用。
    """
    base = Path(__file__).resolve().parent / "_tmp"
    base.mkdir(parents=True, exist_ok=True)
    d = Path(tempfile.mkdtemp(dir=base))
    try:
        yield d
    finally:
        shutil.rmtree(d, ignore_errors=True)
