#!/usr/bin/env python3
"""ImmuHub Analyzer: a lightweight wrapper around NCBI IgBLAST."""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import List, Optional

PACKAGE_ROOT = Path(__file__).resolve().parents[1]
BUNDLED_IGBLAST = PACKAGE_ROOT / "third_party"

SUPPORTED_LOCI = {"IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD"}
D_LOCI = {"IGH", "TRB", "TRD"}


def _open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def ensure_executable(path: Path) -> Path:
    if not path.exists():
        raise FileNotFoundError(f"Executable not found: {path}")

    try:
        path.chmod(path.stat().st_mode | 0o111)
    except Exception:
        pass

    return path


def find_igblast_root(user_root: Optional[str]) -> Path:
    if user_root:
        root = Path(user_root).expanduser().resolve()
    else:
        env = os.environ.get("IGBLAST_HOME")
        root = Path(env).expanduser().resolve() if env else BUNDLED_IGBLAST

    if not root.exists():
        raise FileNotFoundError(f"IGBLAST root not found: {root}")

    return root


def igblast_bin(root: Path, tool: str) -> Path:

    exe = tool

    if os.name == "nt":
        exe += ".exe"

    if tool == "igblastn":
        if (root / "bin" / "igblastn").is_file():
            exe = "igblastn"
        elif (root / "bin" / "igblast").is_file():
            exe = "igblast"

    return ensure_executable(root / "bin" / exe)


def pass_through(tool: str, argv: List[str], root: Optional[str]) -> int:
    home = find_igblast_root(root)
    exe = igblast_bin(home, tool)

    env = os.environ.copy()

    # FORCE OVERRIDE
    env["IGDATA"] = str(home)

    return subprocess.call([str(exe), *argv], env=env)


def fastq_to_fasta(src: Path, dst: Path) -> int:
    n = 0

    with _open_text(src) as fh, open(dst, "w", encoding="utf-8") as out:
        while True:
            name = fh.readline()

            if not name:
                break

            seq = fh.readline().strip()

            # skip plus
            fh.readline()

            qual = fh.readline()

            if not qual:
                raise ValueError("Invalid FASTQ")

            out.write(">" + name[1:].strip().split()[0] + "\n")
            out.write(seq + "\n")

            n += 1

    return n


def detect_input_type(path: Path, explicit: str) -> str:
    if explicit != "auto":
        return explicit

    name = path.name.lower()

    if name.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz")):
        return "fastq"

    if name.endswith((".csv", ".tsv", ".txt")):
        return "cdr3-table"

    return "fasta"


def find_db(
    db_root: Path,
    organism: Optional[str],
    locus: str,
    segment: str
) -> Optional[Path]:

    org = organism.lower() if organism else ""
    loc = locus.upper()
    seg = segment.upper()

    candidates = [
        db_root / f"{org}_{loc}{seg}",
        db_root / f"{org}_{loc}_{seg}",
        db_root / f"{org}_{seg}",
        db_root / f"{loc}{seg}",
        db_root / f"{loc}_{seg}",
    ]

    for c in candidates:
        if (
            c.exists()
            or list(c.parent.glob(c.name + ".n*"))
            or list(c.parent.glob(c.name + ".p*"))
        ):
            return c

    return None


def build_analyze_cmd(args, query: Path, out: Path) -> List[str]:

    root = find_igblast_root(args.igblast_root)

    exe = igblast_bin(
        root,
        "igblastp" if args.protein else "igblastn"
    )

    db_root = (
        Path(args.db_root).expanduser().resolve()
        if args.db_root else None
    )

    locus = args.locus.upper()

    v_db = (
        Path(args.v_db).expanduser().resolve()
        if args.v_db else None
    )

    d_db = (
        Path(args.d_db).expanduser().resolve()
        if args.d_db else None
    )

    j_db = (
        Path(args.j_db).expanduser().resolve()
        if args.j_db else None
    )

    if db_root:
        v_db = v_db or find_db(
            db_root,
            args.organism,
            locus,
            "V"
        )

        d_db = d_db or find_db(
            db_root,
            args.organism,
            locus,
            "D"
        )

        j_db = j_db or find_db(
            db_root,
            args.organism,
            locus,
            "J"
        )

    if not v_db or not j_db:
        raise FileNotFoundError(
            "V/J germline DB not found."
        )

    cmd = [
        str(exe),

        "-query", str(query),

        "-out", str(out),

        "-ig_seqtype",
        "Ig" if locus.startswith("IG") else "TCR",

        "-germline_db_V", str(v_db),

        "-germline_db_J", str(j_db),

        "-num_threads", str(args.threads),
    ]

    if d_db and locus in D_LOCI:
        cmd += ["-germline_db_D", str(d_db)]

    if not args.protein:

        if args.aux_file:

            aux = args.aux_file

            if Path(aux).exists():
                cmd += ["-auxiliary_data", aux]

    if args.format == "airr":
        cmd += ["-outfmt", "19"]

    elif args.format == "txt":
        cmd += ["-outfmt", "7"]

    if args.extra:
        cmd.extend(args.extra)

    return cmd


def run_analyze(args) -> int:

    input_path = Path(args.input).expanduser().resolve()

    out_path = Path(args.output).expanduser().resolve()

    out_path.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    input_type = detect_input_type(
        input_path,
        args.input_type
    )

    tmpdir = tempfile.TemporaryDirectory(
        prefix="immuhub_analyzer_"
    )

    tmp = Path(tmpdir.name)

    query = input_path

    try:

        if input_type == "fastq":

            query = tmp / "query.fasta"

            fastq_to_fasta(
                input_path,
                query
            )

        ig_out = out_path

        cmd = build_analyze_cmd(
            args,
            query,
            ig_out
        )

        env = os.environ.copy()

        root = find_igblast_root(
            args.igblast_root
        )

        # FORCE OVERRIDE
        env["IGDATA"] = str(root)

        if args.print_cmd:
            print(" ".join(map(str, cmd)))

        code = subprocess.call(
            cmd,
            env=env
        )

        return code

    finally:
        tmpdir.cleanup()


def main(argv: Optional[List[str]] = None) -> int:

    argv = list(
        sys.argv[1:] if argv is None else argv
    )

    passthrough_tools = {
        "igblastn",
        "igblastp",
        "makeblastdb",
        "blastdbcmd"
    }

    igblast_root = None

    if "--igblast-root" in argv:

        idx = argv.index("--igblast-root")

        igblast_root = argv[idx + 1]

        del argv[idx:idx + 2]

    if argv and argv[0] in passthrough_tools:
        return pass_through(
            argv[0],
            argv[1:],
            igblast_root
        )

    parser = argparse.ArgumentParser(
        prog="immuhub-analyzer",
        description="ImmuHub Analyzer"
    )

    parser.add_argument(
        "--igblast-root",
        default=igblast_root
    )

    sub = parser.add_subparsers(
        dest="cmd",
        required=True
    )

    for tool in passthrough_tools:
        sub.add_parser(tool)

    # analyze
    p = sub.add_parser("analyze")

    # alias
    sub._name_parser_map["run"] = p

    p.add_argument(
        "-i",
        "--input",
        required=True
    )

    p.add_argument(
        "-o",
        "--output",
        required=True
    )

    p.add_argument(
        "--input-type",
        choices=[
            "auto",
            "fasta",
            "fastq",
            "cdr3-table"
        ],
        default="auto"
    )

    p.add_argument(
        "--locus",
        required=True
    )

    p.add_argument("--organism")

    p.add_argument("--db-root")

    p.add_argument("--v-db")

    p.add_argument("--d-db")

    p.add_argument("--j-db")

    p.add_argument("--aux-file")

    p.add_argument(
        "--protein",
        action="store_true"
    )

    p.add_argument(
        "--format",
        choices=["txt", "airr"],
        default="airr"
    )

    p.add_argument(
        "--threads",
        type=int,
        default=1
    )

    p.add_argument(
        "--extra",
        nargs=argparse.REMAINDER
    )

    p.add_argument(
        "--print-cmd",
        action="store_true"
    )

    ns = parser.parse_args(argv)

    try:

        if ns.cmd == "analyze":
            return run_analyze(ns)

        if ns.cmd == "run":
            return run_analyze(ns)

        parser.error("Unknown command")

    except Exception as e:

        print(
            f"[ImmuHub Analyzer ERROR] {e}",
            file=sys.stderr
        )

        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())