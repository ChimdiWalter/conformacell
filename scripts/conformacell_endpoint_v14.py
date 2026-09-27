#!/usr/bin/env python3
"""Offline audit of the v14 metadata packet; never an endpoint access grant.

Quarantined assets are streamed only through SHA256. No assay parser, model,
network client, evaluator, or admission override is provided.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import resource
import time

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = Path("data/conformacell/v14/qualification_manifest.json")
FIELDS = frozenset((
    "protein_system", "construct", "mutation_or_intervention", "sequence",
    "residue_numbering", "species", "glycoform", "oligomeric_state",
    "ligand_state", "buffer", "temperature", "concentration",
    "measurement_modality", "raw_data_availability", "used_for_model_generation",
    "used_for_fitting", "used_for_selection", "used_for_interpretation",
    "independent_for_validation", "license", "source",
))


def sha256(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def bound_path(root: Path, name: str) -> Path:
    """Require an ordinary path inside the owned v14 data directory."""
    relative = Path(name)
    if relative.is_absolute() or ".." in relative.parts:
        raise ValueError("invalid asset path")
    allowed = (root / "data/conformacell/v14").resolve()
    path = root / relative
    if not path.resolve().is_relative_to(allowed):
        raise ValueError("asset escapes owned data directory")
    current = path
    while current != root:
        if current.is_symlink():
            raise ValueError("symlink asset is not allowed")
        current = current.parent
    return path


def check_packet(root: Path, manifest: dict) -> dict:
    """Check recorded consistency, without certifying scientific adequacy."""
    errors = []
    if manifest.get("schema_version") != 14:
        errors.append("SCHEMA_VERSION")
    if set(manifest.get("required_fields", [])) != FIELDS:
        errors.append("REQUIRED_FIELD_DECLARATION")
    candidates = manifest.get("candidates", [])
    if not candidates or sorted(c.get("rank", -1) for c in candidates) != list(range(1, len(candidates) + 1)):
        errors.append("CANDIDATE_RANKS")
    identifiers = [c.get("candidate_id") for c in candidates]
    if len(set(identifiers)) != len(identifiers) or any(not value for value in identifiers):
        errors.append("CANDIDATE_IDENTITIES")
    for candidate in candidates:
        name = candidate.get("candidate_id", "unknown")
        for field in FIELDS:
            if not candidate.get(field):
                errors.append(f"MISSING_FIELD:{name}:{field}")
        # This packet records unresolved qualification. It cannot promote a
        # candidate merely by editing a Boolean or deleting a gap string.
        if candidate.get("status") != "NOT_ADMITTED" or candidate.get("selected") is not False:
            errors.append(f"UNSUPPORTED_ADMISSION:{name}")
        if not candidate.get("gaps") or not candidate.get("independent_for_validation", "").startswith("NOT_ADMITTED"):
            errors.append(f"UNSUPPORTED_COMPLETENESS:{name}")
    if manifest.get("selected_system") is not None:
        errors.append("UNSUPPORTED_SELECTION")
    for field in ("numeric_endpoint_values_exposed_from_acquired_outcome_assets", "raw_assay_cells_read"):
        if manifest.get(field) is not False:
            errors.append(f"ACCESS_DECLARATION:{field}")
    if manifest.get("scientific_models_run") != 0:
        errors.append("UNEXPECTED_MODEL_RUN")
    readiness = manifest.get("sequential_measurement_readiness", {})
    if readiness.get("status") != "NOT_ADMITTED" or readiness.get("final_outcome_to_policy") is not False or not readiness.get("gates"):
        errors.append("SEQUENTIAL_TASK_NOT_HELD")

    pins = manifest.get("asset_sha256", {})
    if not pins:
        errors.append("MISSING_ASSET_PINS")
    verified = 0
    receipt_total = retained = discarded = receipt_count = 0
    for name, pin in pins.items():
        try:
            path = bound_path(root, name)
            if path.stat().st_size != pin["bytes"] or sha256(path) != pin["sha256"]:
                errors.append(f"ASSET_BINDING:{name}")
                continue
            expected_role = "QUARANTINED_BYTES_NO_OUTCOME_ACCESS" if "quarantine" in Path(name).parts else "METADATA"
            if pin.get("role") != expected_role:
                errors.append(f"ASSET_ROLE:{name}")
            verified += 1
            if path.name.startswith("acquisition_") and path.suffix == ".json":
                receipt = json.loads(path.read_text())
                charge = 0
                for entry in receipt["receipts"]:
                    if entry.get("status") == "FAILED":
                        count = entry.get("recorded_bytes", 0)
                        if count != entry.get("bytes_received_discarded", count):
                            errors.append(f"FAILED_TRANSFER_ACCOUNTING:{name}")
                        discarded += count
                    else:
                        count = entry["bytes"]
                        retained += count
                    if not isinstance(count, int) or count < 0:
                        raise ValueError("invalid byte charge")
                    charge += count
                if charge != receipt["recorded_download_bytes"]:
                    errors.append(f"RECEIPT_TOTAL:{name}")
                receipt_total += charge
                receipt_count += 1
        except (OSError, ValueError, KeyError, TypeError) as error:
            errors.append(f"ASSET_ERROR:{name}:{type(error).__name__}")
    usage = manifest.get("resource_usage", {})
    if not receipt_count or receipt_total != usage.get("metered_download_bytes"):
        errors.append("DOWNLOAD_TOTAL")
    if retained != usage.get("retained_successful_payload_bytes") or discarded != usage.get("failed_capped_received_discarded_bytes"):
        errors.append("DOWNLOAD_COMPONENTS")
    if not isinstance(usage.get("download_cap_bytes"), int) or receipt_total > usage.get("download_cap_bytes", -1):
        errors.append("DOWNLOAD_CAP")

    # FASTA is construct metadata, not an experimental outcome. Keep its
    # independent sequence hash in addition to the complete-file hash.
    for candidate in candidates:
        sequence = candidate.get("sequence")
        if isinstance(sequence, dict):
            for value in sequence.values():
                if isinstance(value, dict) and "file" in value:
                    try:
                        path = bound_path(root, "data/conformacell/v14/" + value["file"])
                        if str(path.relative_to(root)) not in pins:
                            raise ValueError("sequence missing asset pin")
                        lines = path.read_text().splitlines()
                        if not lines or not lines[0].startswith(">") or any(line.startswith(">") for line in lines[1:]):
                            raise ValueError("expected one FASTA chain")
                        chain = "".join(lines[1:]).strip()
                        if len(chain) != value["length"] or hashlib.sha256(chain.encode()).hexdigest() != value["sha256"]:
                            errors.append("SEQUENCE_BINDING")
                    except (OSError, ValueError, KeyError):
                        errors.append("SEQUENCE_ERROR")
    return {
        "status": "METADATA_CONSISTENCY_PASS" if not errors else "FAIL",
        "errors": errors, "candidates": len(candidates), "verified_asset_pins": verified,
        "metered_download_bytes": receipt_total, "failed_transfer_bytes": discarded,
        "selected_system": None, "endpoint_access_authorized": False,
        "scientific_independence_certified": False, "OS_access_control_enforced": False,
        "numeric_outcomes_parsed_by_check": False,
        "scope": "Metadata declarations and hash consistency only; original qualification and exposure history remain necessary.",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", required=True)
    args = parser.parse_args()
    del args
    resource.setrlimit(resource.RLIMIT_AS, (2 * 1024**3, 2 * 1024**3))
    resource.setrlimit(resource.RLIMIT_CPU, (60, 60))
    start = time.monotonic()
    result = check_packet(ROOT, json.loads((ROOT / MANIFEST).read_text()))
    result["manifest_sha256"] = sha256(ROOT / MANIFEST)
    result["wall_seconds"] = time.monotonic() - start
    print(json.dumps(result, indent=2))
    return 0 if result["status"] == "METADATA_CONSISTENCY_PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
