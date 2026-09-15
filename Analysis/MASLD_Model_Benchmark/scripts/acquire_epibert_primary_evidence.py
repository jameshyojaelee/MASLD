#!/usr/bin/env python3
"""Freeze EpiBERT repository history and the complete PMC OA paper package."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path, PurePosixPath
import shutil
import subprocess
import tarfile
import tempfile
from urllib.parse import urlparse
from xml.etree import ElementTree


class EvidenceAcquisitionError(RuntimeError):
    """Raised when a primary-source acquisition differs from its requirement."""


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def run(*command: str, cwd: Path | None = None) -> str:
    result = subprocess.run(
        command,
        cwd=cwd,
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout


def safe_member_name(value: str) -> PurePosixPath:
    if not value or value.startswith("/") or "\\" in value or "\x00" in value:
        raise EvidenceAcquisitionError("PMC archive member has an unsafe name")
    path = PurePosixPath(value)
    if any(part in {"", ".", ".."} for part in path.parts):
        raise EvidenceAcquisitionError("PMC archive member traverses its destination")
    return path


def oa_package_url(payload: bytes, pmc_id: str) -> str:
    root = ElementTree.fromstring(payload)
    records = root.findall(".//record")
    if len(records) != 1 or records[0].get("id") != pmc_id:
        raise EvidenceAcquisitionError("PMC OA record identity differs")
    links = [
        row.get("href", "")
        for row in records[0].findall("link")
        if row.get("format") == "tgz"
    ]
    if len(links) != 1:
        raise EvidenceAcquisitionError("PMC OA tgz link roster differs")
    parsed = urlparse(links[0])
    if parsed.scheme not in {"ftp", "https"} or parsed.netloc != "ftp.ncbi.nlm.nih.gov":
        raise EvidenceAcquisitionError("PMC OA tgz origin differs")
    return links[0]


def download(url: str, output: Path) -> None:
    run(
        "/usr/bin/curl",
        "--fail",
        "--silent",
        "--show-error",
        "--location",
        "--proto",
        "=https",
        "--tlsv1.2",
        url,
        "--output",
        str(output),
    )
    if not output.is_file() or output.stat().st_size == 0:
        raise EvidenceAcquisitionError(f"download is empty: {url}")


def download_ncbi_oa_package(source_url: str, output: Path) -> dict[str, object]:
    """Download the exact OA-API object, retaining every failed transport attempt."""

    parsed = urlparse(source_url)
    if parsed.scheme not in {"ftp", "https"} or parsed.netloc != "ftp.ncbi.nlm.nih.gov":
        raise EvidenceAcquisitionError("PMC OA package origin differs")
    candidates = [source_url]
    if parsed.scheme == "ftp":
        candidates.append(f"https://ftp.ncbi.nlm.nih.gov{parsed.path}")
    attempts = []
    for index, candidate in enumerate(candidates, start=1):
        candidate_parsed = urlparse(candidate)
        protocol = candidate_parsed.scheme
        temporary = output.with_name(f"{output.name}.attempt-{index}")
        command = [
            "/usr/bin/curl",
            "--fail",
            "--silent",
            "--show-error",
            "--location",
            "--retry",
            "3",
            "--retry-delay",
            "2",
            "--user-agent",
            "EpiBERT-primary-evidence-audit/1.0",
            "--proto",
            f"={protocol}",
        ]
        if protocol == "ftp":
            command.append("--ftp-pasv")
        command.extend([candidate, "--output", str(temporary)])
        try:
            run(*command)
        except subprocess.CalledProcessError as error:
            attempts.append(
                {
                    "url": candidate,
                    "protocol": protocol,
                    "status": "failed",
                    "curl_exit_code": error.returncode,
                    "retained_partial": temporary.name if temporary.exists() else None,
                }
            )
            continue
        if not temporary.is_file() or temporary.stat().st_size == 0:
            raise EvidenceAcquisitionError("PMC OA package download is empty")
        temporary.rename(output)
        attempts.append(
            {
                "url": candidate,
                "protocol": protocol,
                "status": "passed",
                "curl_exit_code": 0,
                "retained_partial": None,
            }
        )
        return {
            "status": "passed",
            "source_url": source_url,
            "download_url": candidate,
            "attempts": attempts,
        }
    return {
        "status": "failed",
        "source_url": source_url,
        "download_url": None,
        "attempts": attempts,
    }


def jats_supplement_hrefs(payload: bytes) -> list[str]:
    root = ElementTree.fromstring(payload)
    xlink = "{http://www.w3.org/1999/xlink}href"
    values = []
    for supplement in root.findall(".//supplementary-material"):
        for element in supplement.iter():
            href = element.get(xlink)
            if href and href not in values:
                values.append(href)
    if not values:
        raise EvidenceAcquisitionError("PMC JATS has no supplementary-material href")
    for value in values:
        parsed = urlparse(value)
        if parsed.scheme or parsed.netloc or PurePosixPath(parsed.path).name != parsed.path:
            raise EvidenceAcquisitionError("PMC JATS supplementary href is not a basename")
        if PurePosixPath(value).suffix.lower() not in {".pdf", ".xlsx", ".xls", ".docx", ".zip"}:
            raise EvidenceAcquisitionError("PMC JATS supplementary type differs")
    return values


def jats_supplement_inventory(payload: bytes) -> list[dict[str, str]]:
    root = ElementTree.fromstring(payload)
    xlink = "{http://www.w3.org/1999/xlink}href"
    rows = []
    for supplement in root.findall(".//supplementary-material"):
        title_node = supplement.find("caption/title")
        title = (
            " ".join("".join(title_node.itertext()).split())
            if title_node is not None
            else ""
        )
        hrefs = [
            element.get(xlink)
            for element in supplement.iter()
            if element.get(xlink)
        ]
        if len(hrefs) != 1:
            raise EvidenceAcquisitionError("PMC JATS supplement member roster differs")
        rows.append(
            {
                "id": str(supplement.get("id", "")),
                "title": title,
                "href": str(hrefs[0]),
            }
        )
    if [row["href"] for row in rows] != jats_supplement_hrefs(payload):
        raise EvidenceAcquisitionError("PMC JATS supplement ordering differs")
    return rows


def download_ncbi_article_evidence(pmc_id: str, paper: Path) -> dict[str, object]:
    efetch_url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        f"?db=pmc&id={pmc_id.removeprefix('PMC')}"
    )
    fulltext = paper / f"{pmc_id}.efetch.xml"
    download(efetch_url, fulltext)
    inventory = jats_supplement_inventory(fulltext.read_bytes())
    article_url = f"https://pmc.ncbi.nlm.nih.gov/articles/{pmc_id}/"
    article_html = paper / f"{pmc_id}.article.html"
    download(article_url, article_html)
    html = article_html.read_text(encoding="utf-8")
    instance_id = pmc_id.removeprefix("PMC")
    rows = []
    for row in inventory:
        relative = f"/articles/instance/{instance_id}/bin/{row['href']}"
        if relative not in html:
            raise EvidenceAcquisitionError("PMC article supplement binding differs")
        rows.append(
            {
                **row,
                "url": f"https://pmc.ncbi.nlm.nih.gov{relative}",
                "binary_retrieval_status": "not_retrieved_protected_article_member",
            }
        )
    return {
        "efetch_url": efetch_url,
        "efetch_sha256": digest(fulltext),
        "efetch_size_bytes": fulltext.stat().st_size,
        "article_url": article_url,
        "article_html_sha256": digest(article_html),
        "article_html_size_bytes": article_html.stat().st_size,
        "supplementary_members": rows,
    }


def extract_pmc_package(archive: Path, output: Path) -> list[dict[str, object]]:
    if output.exists():
        raise EvidenceAcquisitionError("PMC extraction output already exists")
    output.mkdir(mode=0o750)
    rows: list[dict[str, object]] = []
    total_bytes = 0
    with tarfile.open(archive, mode="r:gz") as handle:
        members = handle.getmembers()
        if not members or len(members) > 2000:
            raise EvidenceAcquisitionError("PMC archive member count differs")
        for member in members:
            relative = safe_member_name(member.name)
            destination = output.joinpath(*relative.parts)
            if member.isdir():
                destination.mkdir(mode=0o750, parents=True, exist_ok=True)
                rows.append({"path": relative.as_posix(), "type": "directory", "size_bytes": 0})
                continue
            if not member.isfile():
                raise EvidenceAcquisitionError("PMC archive contains a linked or special member")
            total_bytes += member.size
            if total_bytes > 2_000_000_000:
                raise EvidenceAcquisitionError("PMC archive exceeds the extraction limit")
            destination.parent.mkdir(mode=0o750, parents=True, exist_ok=True)
            source = handle.extractfile(member)
            if source is None:
                raise EvidenceAcquisitionError("PMC archive member is unreadable")
            value = sha256()
            size = 0
            with destination.open("xb") as target:
                for block in iter(lambda: source.read(1024 * 1024), b""):
                    value.update(block)
                    size += len(block)
                    target.write(block)
            if size != member.size:
                raise EvidenceAcquisitionError("PMC archive member size differs")
            rows.append(
                {
                    "path": relative.as_posix(),
                    "type": "regular_file",
                    "size_bytes": size,
                    "sha256": value.hexdigest(),
                }
            )
    return rows


def convert_pdfs(package: Path, output: Path, pdftotext: Path) -> list[dict[str, object]]:
    if not pdftotext.is_file():
        raise EvidenceAcquisitionError("pdftotext executable is absent")
    output.mkdir(mode=0o750)
    rows = []
    for pdf in sorted(package.rglob("*.pdf")):
        relative = pdf.relative_to(package)
        target = output / relative.with_suffix(".txt")
        target.parent.mkdir(mode=0o750, parents=True, exist_ok=True)
        run(str(pdftotext), "-layout", str(pdf), str(target))
        rows.append(
            {
                "source": relative.as_posix(),
                "text": target.relative_to(output).as_posix(),
                "sha256": digest(target),
                "size_bytes": target.stat().st_size,
            }
        )
    return rows


def acquire(args: argparse.Namespace) -> dict[str, object]:
    if args.output.exists():
        raise EvidenceAcquisitionError("evidence output already exists")
    args.output.mkdir(mode=0o750)
    repository = args.output / "repository"
    paper = args.output / "paper"
    repository.mkdir(mode=0o750)
    paper.mkdir(mode=0o750)

    with tempfile.TemporaryDirectory(prefix="epibert-history-") as temporary:
        mirror = Path(temporary) / "EpiBERT.git"
        run("/usr/bin/git", "clone", "--mirror", args.repository_url, str(mirror))
        pinned = run("/usr/bin/git", "--git-dir", str(mirror), "rev-parse", args.pinned_revision).strip()
        tagged = run("/usr/bin/git", "--git-dir", str(mirror), "rev-parse", f"refs/tags/{args.tag}^{{}}").strip()
        if pinned != args.pinned_revision or tagged != args.pinned_revision:
            raise EvidenceAcquisitionError("pinned EpiBERT tag identity differs")
        run(
            "/usr/bin/git",
            "--git-dir",
            str(mirror),
            "bundle",
            "create",
            str(repository / "EpiBERT-all-refs.bundle"),
            "--all",
        )
        (repository / "refs.txt").write_text(
            run("/usr/bin/git", "--git-dir", str(mirror), "show-ref"), encoding="utf-8"
        )
        (repository / "all-history.tsv").write_text(
            run(
                "/usr/bin/git",
                "--git-dir",
                str(mirror),
                "log",
                "--all",
                "--date=iso-strict",
                "--format=%H%x09%aI%x09%cI%x09%P%x09%s",
            ),
            encoding="utf-8",
        )
        bundle_verify = run(
            "/usr/bin/git",
            "bundle",
            "verify",
            str(repository / "EpiBERT-all-refs.bundle"),
        )
        (repository / "bundle-verify.txt").write_text(bundle_verify, encoding="utf-8")

    oa_api = paper / f"{args.pmc_id}.oa.xml"
    download(
        f"https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi?id={args.pmc_id}",
        oa_api,
    )
    package_url = oa_package_url(oa_api.read_bytes(), args.pmc_id)
    package_archive = paper / f"{args.pmc_id}.tar.gz"
    transport = download_ncbi_oa_package(package_url, package_archive)
    package_rows: list[dict[str, object]] = []
    article_evidence: dict[str, object] | None = None
    if transport["status"] == "passed":
        package_rows = extract_pmc_package(package_archive, paper / "package")
        pdf_rows = convert_pdfs(paper / "package", paper / "pdf_text", args.pdftotext)
        evidence_mode = "ncbi_oa_package"
    else:
        article_evidence = download_ncbi_article_evidence(args.pmc_id, paper)
        pdf_rows = []
        evidence_mode = "ncbi_efetch_jats_plus_article_metadata"

    receipt = {
        "schema_version": "masld-bench-epibert-primary-evidence-acquisition-v1",
        "status": "pass",
        "repository": {
            "url": args.repository_url,
            "pinned_revision": args.pinned_revision,
            "tag": args.tag,
            "bundle_sha256": digest(repository / "EpiBERT-all-refs.bundle"),
            "bundle_size_bytes": (repository / "EpiBERT-all-refs.bundle").stat().st_size,
            "commit_count_all_refs": sum(
                1 for line in (repository / "all-history.tsv").read_text(encoding="utf-8").splitlines() if line
            ),
        },
        "paper": {
            "pmc_id": args.pmc_id,
            "oa_api_url": f"https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi?id={args.pmc_id}",
            "oa_api_sha256": digest(oa_api),
            "evidence_mode": evidence_mode,
            "package_url": package_url,
            "package_download_url": transport["download_url"],
            "package_transport_attempts": transport["attempts"],
            "package_sha256": digest(package_archive) if package_archive.is_file() else None,
            "package_size_bytes": package_archive.stat().st_size if package_archive.is_file() else None,
            "members": package_rows,
            "article_evidence": article_evidence,
            "pdf_text": pdf_rows,
        },
        "project_data_read": False,
        "outcomes_read": False,
        "sealed_assets_read": False,
    }
    (args.output / "acquisition_receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository-url", required=True)
    parser.add_argument("--pinned-revision", required=True)
    parser.add_argument("--tag", required=True)
    parser.add_argument("--pmc-id", required=True)
    parser.add_argument("--pdftotext", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = acquire(args)
    print(json.dumps({"status": result["status"], "output": str(args.output)}, sort_keys=True))


if __name__ == "__main__":
    main()
