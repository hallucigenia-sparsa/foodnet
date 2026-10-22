"""The gate's merge-marker check: a conflict left by a merge once reached the changelog unseen."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "checks"))
import gate  # noqa: E402


def test_conflict_markers_are_found_and_a_heading_underline_is_not(monkeypatch):
    files = {"CHANGELOG.md": "# Changes\n<<<<<<< HEAD\n- ours\n=======\n- theirs\n>>>>>>> origin/branch\n",
             "README.txt": "foodnet\n=======\n\nStart it\n--------\n"}
    monkeypatch.setattr(gate, "_read", files.get)
    assert gate.check_merge_markers(list(files)) == [
        "conflict marker left from a merge in CHANGELOG.md, line 2",
        "conflict marker left from a merge in CHANGELOG.md, line 6"]


def test_the_declared_dependency_is_allowed_only_where_it_is_declared(monkeypatch):
    # a review: truststore is a macOS-only dependency, so an import of it anywhere but tls.py would break other
    # systems
    files = {"src/foodnet/tls.py": "def f():\n    import truststore\n",
             "src/foodnet/other.py": "def g():\n    import truststore\n"}
    monkeypatch.setattr(gate, "_read", files.get)
    assert gate.check_self_contained(list(files)) == [
        "src/foodnet/other.py: imports 'truststore' (not stdlib or foodnet; declare it or drop it)"]
