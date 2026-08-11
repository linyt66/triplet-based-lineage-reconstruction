import json
import re
import unittest
from pathlib import Path

import _path  # noqa: F401

PROJECT_ROOT = Path(__file__).resolve().parents[1]
NOTEBOOKS = sorted((PROJECT_ROOT / "notebooks").glob("*.ipynb"))
TEXT_PATTERNS = ("*.py", "*.md", "*.toml", "*.txt", "*.ipynb")
MOJIBAKE_TOKENS = tuple(chr(codepoint) for codepoint in (0x0432, 0x043F, 0x0458, 0x0413))
CJK_PATTERN = re.compile(r"[\u4e00-\u9fff]")


class RepositoryHygieneTests(unittest.TestCase):
    def test_notebooks_have_no_saved_outputs(self):
        offenders = []
        for notebook_path in NOTEBOOKS:
            notebook = json.loads(notebook_path.read_text(encoding="utf-8"))
            for index, cell in enumerate(notebook.get("cells", []), 1):
                if cell.get("cell_type") != "code":
                    continue
                if cell.get("outputs"):
                    offenders.append(f"{notebook_path}:{index}:outputs")
                if cell.get("execution_count") is not None:
                    offenders.append(f"{notebook_path}:{index}:execution_count")
        self.assertEqual(offenders, [])

    def test_text_files_do_not_contain_cjk_or_mojibake(self):
        offenders = []
        for pattern in TEXT_PATTERNS:
            for path in PROJECT_ROOT.glob(f"**/{pattern}"):
                if ".git" in path.parts or ".ipynb_checkpoints" in path.parts:
                    continue
                text = path.read_text(encoding="utf-8", errors="ignore")
                if CJK_PATTERN.search(text) or any(token in text for token in MOJIBAKE_TOKENS):
                    offenders.append(str(path.relative_to(PROJECT_ROOT)))
        self.assertEqual(offenders, [])

    def test_generated_cache_directories_are_not_tracked(self):
        import subprocess

        completed = subprocess.run(
            ["git", "ls-files"],
            cwd=PROJECT_ROOT,
            check=True,
            capture_output=True,
            text=True,
        )
        tracked_generated = [
            path
            for path in completed.stdout.splitlines()
            if "__pycache__/" in path or ".egg-info/" in path
        ]
        self.assertEqual(tracked_generated, [])

    def test_requirements_does_not_pin_python_as_package(self):
        requirements = (PROJECT_ROOT / "requirements.txt").read_text(encoding="utf-8")
        self.assertNotIn("python>=", requirements.lower())


if __name__ == "__main__":
    unittest.main()
