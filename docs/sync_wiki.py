#!/usr/bin/env python3

import argparse
import re
import subprocess
import sys
from pathlib import Path
from urllib.parse import urlparse


ROOT = Path(__file__).resolve().parent.parent
DEFAULT_DOCS = ROOT / "docs"
DEFAULT_WIKI = DEFAULT_DOCS / "residualanalysis.wiki"

SOURCE_FILES = {
    "WikiHome.md": "home.md",
    "L2Residual.md": "L2Residual.md",
    "JER.md": "JER.md",
    "L3Residual.md": "L3Residual.md",
    "Systematics.md": "Systematics.md",
    "Batch.md": "Batch.md",
}

PAGE_ALIASES = {
    "WikiHome": "home",
    "WikiHome.md": "home",
    "residualanalysis": "home",
    "residualanalysis.md": "home",
    "L2Residual.md": "L2Residual",
    "JER.md": "JER",
    "L3Residual.md": "L3Residual",
    "Systematics.md": "Systematics",
    "Batch.md": "Batch",
}


def remote_to_project_url(remote: str) -> str:
    remote = remote.strip()
    if remote.startswith("git@"):
        host_and_path = remote[4:]
        host, path = host_and_path.split(":", 1)
        path = path.removesuffix(".git").removesuffix(".wiki")
        return f"https://{host}/{path}"

    if remote.startswith("ssh://"):
        parsed = urlparse(remote)
        path = parsed.path.lstrip("/").removesuffix(".git").removesuffix(".wiki")
        return f"https://{parsed.hostname}/{path}"

    parsed = urlparse(remote)
    if parsed.scheme in {"http", "https"}:
        path = parsed.path.lstrip("/").removesuffix(".git").removesuffix(".wiki")
        return f"https://{parsed.netloc}/{path}"

    raise ValueError(f"Unsupported remote URL format: {remote}")


def infer_project_url(root: Path) -> str:
    result = subprocess.run(
        ["git", "-C", str(root), "remote", "get-url", "origin"],
        check=True,
        capture_output=True,
        text=True,
    )
    return remote_to_project_url(result.stdout.strip())


def split_target(target: str) -> tuple[str, str]:
    if "#" in target:
        base, anchor = target.split("#", 1)
        return base, f"#{anchor}"
    return target, ""


def wiki_slug_for_target(target: str) -> str:
    base, anchor = split_target(target)
    page = PAGE_ALIASES.get(base, Path(base).stem)
    return page + anchor


def repo_blob_url(project_url: str, docs_dir: Path, src_file: Path, target: str) -> str:
    base, anchor = split_target(target)
    resolved = (src_file.parent / base).resolve()
    rel_path = resolved.relative_to(ROOT).as_posix()
    return f"{project_url}/-/blob/HEAD/{rel_path}{anchor}"


def rewrite_for_wiki(text: str, src_file: Path, docs_dir: Path, project_url: str) -> str:
    def replace_link(match: re.Match[str]) -> str:
        label = match.group(1)
        target = match.group(2)

        if target.startswith(("http://", "https://", "mailto:", "#")):
            return match.group(0)

        base, _anchor = split_target(target)
        if base.endswith(".md"):
            return f"[{label}]({wiki_slug_for_target(target)})"

        if target.startswith(("../", "./")):
            return f"[{label}]({repo_blob_url(project_url, docs_dir, src_file, target)})"

        return match.group(0)

    body = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", replace_link, text)
    return "<!-- Generated from docs/. Edit docs/*.md and run docs/sync_wiki.py. -->\n\n" + body


def write_redirects(wiki_dir: Path) -> None:
    gitlab_dir = wiki_dir / ".gitlab"
    gitlab_dir.mkdir(parents=True, exist_ok=True)
    redirects = {
        "residualanalysis": "home",
        "residualanalysis.md": "home",
        "WikiHome": "home",
        "WikiHome.md": "home",
        "L2Residual.md": "L2Residual",
        "JER.md": "JER",
        "L3Residual.md": "L3Residual",
        "Systematics.md": "Systematics",
        "Batch.md": "Batch",
    }
    lines = [f'"{key}": "{value}"' for key, value in redirects.items()]
    (gitlab_dir / "redirects.yml").write_text("---\n" + "\n".join(lines) + "\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate the GitLab wiki pages from docs/*.md. "
            "Use --wiki-dir to point at an external clone of <project>.wiki.git."
        )
    )
    parser.add_argument("--docs-dir", default=str(DEFAULT_DOCS), help="Path to the canonical docs directory")
    parser.add_argument("--wiki-dir", default=str(DEFAULT_WIKI), help="Path to the GitLab wiki repository checkout")
    parser.add_argument(
        "--project-url",
        default=None,
        help="Project web URL, for example https://gitlab.cern.ch/group/project. If omitted, inferred from git origin.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    docs_dir = Path(args.docs_dir).resolve()
    wiki_dir = Path(args.wiki_dir).resolve()
    project_url = args.project_url or infer_project_url(ROOT)
    project_url = project_url.rstrip("/")

    wiki_dir.mkdir(parents=True, exist_ok=True)
    for src_name, dst_name in SOURCE_FILES.items():
        src = docs_dir / src_name
        dst = wiki_dir / dst_name
        content = rewrite_for_wiki(src.read_text(), src, docs_dir, project_url)
        if dst_name == "home.md":
            content = "---\ntitle: Home\n---\n" + content
        dst.write_text(content)

    write_redirects(wiki_dir)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"sync_wiki.py failed: {exc}", file=sys.stderr)
        raise