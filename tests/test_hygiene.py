from pathlib import Path

from tests.helpers import get_repo_root


FORBIDDEN_PATHS = [
    "/" + "home" + "/shanav",
    "/" + "import" + "/home4" + "/shanav",
]
TEXT_SUFFIXES = {
    ".py",
    ".md",
    ".yaml",
    ".yml",
    ".toml",
    ".json",
    ".txt",
    ".tsv",
    ".csv",
    ".sh",
    ".R",
    ".template",
}


def test_no_hardcoded_private_paths_in_public_tree() -> None:
    repo_root = get_repo_root()
    offenders = []

    for path in repo_root.rglob("*"):
        if not path.is_file():
            continue
        if "__pycache__" in path.parts or ".pytest_cache" in path.parts:
            continue
        if path.suffix not in TEXT_SUFFIXES:
            continue

        text = path.read_text()
        if any(forbidden in text for forbidden in FORBIDDEN_PATHS):
            offenders.append(str(path.relative_to(repo_root)))

    assert offenders == []


def test_no_parallel_scripts_shipped() -> None:
    repo_root = get_repo_root()
    assert list(repo_root.glob("parallel*.sh")) == []
