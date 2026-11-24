# Copyright 2026 MV2, Munim, Inc. Licensed under the Apache License, Version 2.0 (Apache-2.0).
# Made with AI agents working in MV2's name.
"""Fetch the registered v1 code (magpie_order.py, fetch_data.py) from GitHub and check each file
against the git blob SHA-1 it had when amendment 1 was written, so the amendment runs on exactly
the code registered on 2026-11-24. Exits non-zero (nothing computed) on any mismatch."""
import hashlib
import sys
import time
import urllib.request

RAW = "https://raw.githubusercontent.com/munimtechnologies/magpie-order-test/main/"
BLOBS = {
    "magpie_order.py": "7f2599f3a1f5899c9da31e7adfd2648039d527ad",
    "fetch_data.py": "8317856a80b1bfc60e96a33f8e0a1b4655cf81fe",
}


def git_blob_sha(blob):
    return hashlib.sha1(b"blob %d\x00" % len(blob) + blob).hexdigest()


def main():
    ok = True
    for name, sha in BLOBS.items():
        blob = None
        for attempt in range(5):
            try:
                with urllib.request.urlopen(RAW + name + "?t=%d" % int(time.time()),
                                            timeout=60) as r:
                    blob = r.read()
                break
            except Exception as exc:  # noqa: BLE001
                print(name, "fetch failed:", exc, flush=True)
                time.sleep(20 * (attempt + 1))
        if blob is None:
            ok = False
            continue
        got = git_blob_sha(blob)
        print(name, got, "OK" if got == sha else "DIFFERENT from " + sha, flush=True)
        if got != sha:
            ok = False
            continue
        with open(name, "wb") as fh:
            fh.write(blob)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
