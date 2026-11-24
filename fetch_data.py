# Copyright 2026 MV2, Munim, Inc. Licensed under the Apache License, Version 2.0 (Apache-2.0).
# Made with AI agents working in MV2's name.
"""Download the authors' repertoire table (figshare doi:10.6084/m9.figshare.28251398 v2, CC BY 4.0)
and check it against the size and MD5 that figshare's API published (read 2026-11-24).
Usage: python fetch_data.py DEST_DIR"""
import hashlib
import json
import os
import sys
import time
import urllib.request

FILES = {
    "COMBINATION_repertoire.csv": (51834599, 78432, "f60a6af46f0161225b5965eb2c917874"),
}
URL = "https://ndownloader.figshare.com/files/%d"


def fetch(file_id):
    last = None
    for attempt in range(5):
        try:
            req = urllib.request.Request(URL % file_id, headers={"User-Agent": "research-data-fetch/1.0"})
            with urllib.request.urlopen(req, timeout=120) as r:
                return r.read()
        except Exception as exc:  # noqa: BLE001
            last = exc
            time.sleep(10 * (attempt + 1))
    raise RuntimeError("download failed for %d: %r" % (file_id, last))


def main(dest):
    os.makedirs(dest, exist_ok=True)
    report, ok = {}, True
    for name, (fid, size, md5) in FILES.items():
        path = os.path.join(dest, name)
        blob = None
        if os.path.exists(path):
            with open(path, "rb") as fh:
                blob = fh.read()
            if hashlib.md5(blob).hexdigest() != md5:
                blob = None
        if blob is None:
            blob = fetch(fid)
            with open(path, "wb") as fh:
                fh.write(blob)
        got = hashlib.md5(blob).hexdigest()
        good = got == md5 and len(blob) == size
        ok = ok and good
        report[name] = {"bytes": len(blob), "md5": got, "sha256": hashlib.sha256(blob).hexdigest(),
                        "matches_figshare": good}
        print(name, len(blob), got, "OK" if good else "MISMATCH", flush=True)
    with open(os.path.join(dest, "fetch_report_order.json"), "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=1)
    if not ok:
        print("data check failed; nothing analysed", flush=True)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
