"""Reserve or version and publish the ARID1A–BAF software archive on Zenodo.

The token is read from ``ZENODO_TOKEN`` or an interactive masked prompt. It is
never accepted as a command-line argument, written to disk, or printed.
"""

from __future__ import annotations

import argparse
import getpass
import hashlib
import json
import os
from datetime import date
from pathlib import Path

import requests


API = "https://zenodo.org/api"
REPOSITORY = "https://github.com/ArleneHarrison/arid1a-baf-collateral-dependency"


def token() -> str:
    value = os.environ.get("ZENODO_TOKEN") or getpass.getpass("Paste Zenodo token (input hidden): ")
    if not value.strip():
        raise SystemExit("No Zenodo token was supplied")
    return value.strip()


def request(method: str, url: str, access_token: str, **kwargs) -> requests.Response:
    headers = dict(kwargs.pop("headers", {}))
    headers["Authorization"] = f"Bearer {access_token}"
    response = requests.request(method, url, headers=headers, timeout=600, **kwargs)
    if not response.ok:
        detail = response.text[:1000].replace(access_token, "[REDACTED]")
        raise RuntimeError(f"Zenodo returned HTTP {response.status_code}: {detail}")
    return response


def metadata(version: str) -> dict:
    return {
        "metadata": {
            "title": "ARID1A–BAF collateral dependency: reproducibility repository",
            "upload_type": "software",
            "description": (
                "Versioned code, compact derived tables, automated numerical checks, "
                "and manuscript materials for a structure-informed analysis of "
                "ARID1A-associated BAF-complex co-dependency in cancer cell lines."
            ),
            "creators": [
                {"name": "Tao, Zhiyong", "orcid": "0009-0001-2454-2432"},
                {"name": "Xiong, Juan"},
                {"name": "Zhang, Yiwei", "orcid": "0009-0008-7513-7956"},
                {"name": "Zhong, Chulin", "orcid": "0009-0001-0840-3846"},
                {"name": "Gao, Wenkai", "orcid": "0009-0000-2900-9832"},
                {"name": "Fan, Guohua", "orcid": "0009-0005-7470-7863"},
            ],
            "keywords": [
                "ARID1A", "BAF complex", "co-dependency", "DepMap", "NSCLC",
                "structural biology", "reproducibility",
            ],
            "license": "MIT",
            "access_right": "open",
            "language": "eng",
            "version": version,
            "publication_date": date.today().isoformat(),
            "related_identifiers": [
                {"identifier": REPOSITORY, "relation": "isSupplementTo", "scheme": "url"},
                {
                    "identifier": f"{REPOSITORY}/releases/tag/v{version}",
                    "relation": "isIdenticalTo",
                    "scheme": "url",
                },
            ],
        }
    }


def reserve(args: argparse.Namespace) -> None:
    access_token = token()
    if args.previous_deposit_id:
        version_action = request(
            "POST",
            f"{API}/deposit/depositions/{args.previous_deposit_id}/actions/newversion",
            access_token,
        ).json()
        latest_draft = version_action.get("links", {}).get("latest_draft")
        if not latest_draft:
            raise RuntimeError("Zenodo did not return a latest_draft link")
        response = request("GET", latest_draft, access_token).json()
    else:
        response = request(
            "POST",
            f"{API}/deposit/depositions",
            access_token,
            json=metadata(args.version),
        ).json()
    reserved = response.get("metadata", {}).get("prereserve_doi", {})
    receipt = {
        "deposit_id": response["id"],
        "reserved_doi": reserved.get("doi"),
        "reserved_recid": reserved.get("recid"),
        "version": args.version,
        "draft_url": response.get("links", {}).get("html"),
        "state": "draft",
    }
    output = Path(args.receipt)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2))


def finalize(args: argparse.Namespace) -> None:
    access_token = token()
    receipt_path = Path(args.receipt)
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    deposit_id = receipt["deposit_id"]
    archive = Path(args.archive)
    if not archive.is_file():
        raise FileNotFoundError(archive)

    deposit = request(
        "GET", f"{API}/deposit/depositions/{deposit_id}", access_token
    ).json()
    for old_file in deposit.get("files", []):
        delete_url = old_file.get("links", {}).get("self")
        if not delete_url:
            raise RuntimeError("An inherited Zenodo file has no deletion link")
        request("DELETE", delete_url, access_token)
    bucket = deposit["links"]["bucket"]
    with archive.open("rb") as stream:
        request("PUT", f"{bucket}/{archive.name}", access_token, data=stream)
    request(
        "PUT",
        f"{API}/deposit/depositions/{deposit_id}",
        access_token,
        json=metadata(receipt["version"]),
    )

    if not args.publish:
        print(f"Draft updated: {deposit['links']['html']}")
        return
    published = request(
        "POST",
        f"{API}/deposit/depositions/{deposit_id}/actions/publish",
        access_token,
    ).json()
    receipt.update(
        {
            "state": "published",
            "doi": published.get("doi"),
            "doi_url": published.get("doi_url") or published.get("links", {}).get("doi"),
            "record_url": published.get("links", {}).get("html"),
            "archive": archive.name,
            "archive_bytes": archive.stat().st_size,
            "archive_sha256": hashlib.sha256(archive.read_bytes()).hexdigest(),
        }
    )
    receipt_path.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    reserve_parser = subparsers.add_parser(
        "reserve", help="Create a first draft or a new-version draft and reserve its DOI"
    )
    reserve_parser.add_argument("--version", required=True)
    reserve_parser.add_argument("--receipt", required=True)
    reserve_parser.add_argument(
        "--previous-deposit-id", type=int,
        help="Published Zenodo deposit to version; omit only for the first release",
    )
    reserve_parser.set_defaults(function=reserve)
    finalize_parser = subparsers.add_parser("finalize", help="Upload the final archive and optionally publish")
    finalize_parser.add_argument("--receipt", required=True)
    finalize_parser.add_argument("--archive", required=True)
    finalize_parser.add_argument("--publish", action="store_true")
    finalize_parser.set_defaults(function=finalize)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.function(args)


if __name__ == "__main__":
    main()
