#!/usr/bin/env python3
"""Verify original arXiv archive checksums and the extracted LaTeX source trees."""
import gzip
import hashlib
import json
from pathlib import Path
import tarfile

root = Path(__file__).resolve().parents[1]
records = json.loads((root / 'references/checksums.json').read_text())
for record in records:
    archive = root / record['file']
    digest = hashlib.sha256(archive.read_bytes()).hexdigest()
    if digest != record['sha256']:
        raise SystemExit(f'archive checksum mismatch: {archive}')
    destination = archive.with_suffix('')
    try:
        source = tarfile.open(archive)
    except tarfile.ReadError:
        if (destination / 'main.tex').read_bytes() != gzip.decompress(archive.read_bytes()):
            raise SystemExit(f'extracted source mismatch: {destination}')
    else:
        with source:
            for member in source:
                if member.isfile():
                    target = (destination / member.name).resolve()
                    if not target.is_relative_to(destination.resolve()):
                        raise SystemExit(f'unsafe archive member: {member.name}')
                    if target.read_bytes() != source.extractfile(member).read():
                        raise SystemExit(f'extracted source mismatch: {target}')
    print(f'Verified {record["file"]} and extracted sources')
