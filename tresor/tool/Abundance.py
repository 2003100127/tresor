__version__ = "0.0.1"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__ = "jianfeng.sunmt@gmail.com"

from typing import List, Tuple, Iterable, Optional, Dict

import gzip, math, random, re
from dataclasses import dataclass
from tresor.util.random.Sample import Sample


def open_maybe_gz(path: str, mode: str = "rt"):
    return gzip.open(path, mode=mode) if path.endswith(".gz") else open(path, mode=mode, encoding="utf-8")


def parse_header_id(header: str, use_pipe_if_present: bool = True, pipe_index: int = 0, strip_version: bool = True) -> str:
    """
    header: the full header line WITHOUT leading '>' and WITHOUT trailing '\n'
    - If '|' present and use_pipe_if_present=True: split by '|' and take field [pipe_index].
      (GENCODE transcripts.fa.gz: 0=ENST..., 1=ENSG..., 4=gene name-transcript name, 5=gene name ...)
    - Else: take token up to first whitespace.
    - If strip_version=True: drop trailing '.<digits>' (e.g., ENST00000123.4 -> ENST00000123)
    """
    token = header.split()[0]
    if use_pipe_if_present and '|' in token:
        parts = token.split('|')
        if pipe_index < 0 or pipe_index >= len(parts):
            raise ValueError(f"pipe_index {pipe_index} out of range for header: {header}")
        token = parts[pipe_index]
    if strip_version:
        # remove . when ENST/ENSG/ENSP/NR_/NM_
        if '.' in token and token.split('.')[0]:
            left, right = token.rsplit('.', 1)
            if right.isdigit():
                token = left
    return token

def iter_fasta_ids(
        fasta_fpn: str,
        use_pipe_if_present: bool,
        pipe_index: int,
        strip_version: bool,
) -> Iterable[str]:
    with open_maybe_gz(fasta_fpn, "rt") as fh:
        for line in fh:
            if line.startswith('>'):
                yield parse_header_id(line[1:].rstrip('\n'), use_pipe_if_present, pipe_index, strip_version)

def largest_remainder_to_int_counts(weights: List[float], total: int) -> List[int]:
    s = sum(weights) or 1.0
    raw = [w / s * total for w in weights]
    base = [int(math.floor(x)) for x in raw]
    rem = total - sum(base)
    if rem <= 0:
        return base
    frac = sorted(((raw[i] - base[i], i) for i in range(len(raw))), reverse=True)
    for k in range(rem):
        _, idx = frac[k]
        base[idx] += 1
    return base


def largest_remainder_to_int_counts_gtf(weights: List[float], total: int) -> List[int]:
    s = float(sum(weights)) or 1.0
    raw = [w / s * total for w in weights]
    base = [int(math.floor(x)) for x in raw]
    rem = total - sum(base)
    if rem > 0:
        frac = sorted(((raw[i] - base[i], i) for i in range(len(raw))), reverse=True)
        for k in range(rem):
            _, idx = frac[k]
            base[idx] += 1
    return base



@dataclass
class AbundanceParams:

    fasta_fpn: str = None
    gtf_fpn: str = None
    output_tsv: str = "abundance.tsv"
    num_selected_mols: int = 100
    select_mode: str = "random" # "random" | "first"
    num_total_mols: int = 10000
    seed: int = 1

    # GENCODE header
    # Keep the same field structure as the FASTA version
    # `strip_version` is also useful in the GTF context (controls whether `pure` removes version numbers)
    use_pipe_if_present: bool = True # see '|'
    pipe_index: int = 0 # 0=ENST, 1=ENSG, 5=GENE SYMBOL
    strip_version: bool = True # remove . version surfix

    # dist
    dist: str = "dirichlet" # "dirichlet"|"lognormal"|"zipf"|"uniform"
    alpha: float = 1.0 # dirichlet
    mu: float = 0.0 # lognormal
    sigma: float = 1.0 # lognormal
    s: float = 1.0 # zipf

    min_count: int = 0
    sort_output_by: str = "count" # "count"|"id"
    subset_ids_txt: Optional[str] = "refs.sub.ids.txt"

    include_noncoding: bool = True           # Whether non-coding transcripts are included; aligned with TKSM --non-coding.
    versioned_out: str = "abundance.versioned.tsv"


class FromFasta:

    def __init__(
            self,
            param: AbundanceParams,
    ):
        self.param = param
        self.rng = random.Random(param.seed)

    def _read_ids(self) -> List[str]:
        ids = list(iter_fasta_ids(
            self.param.fasta_fpn,
            self.param.use_pipe_if_present,
            self.param.pipe_index,
            self.param.strip_version
        ))
        if not ids:
            raise RuntimeError(f"No IDs parsed from {self.param.fasta_fpn}")
        return ids

    def _select(self, ids: List[str]) -> List[str]:
        n = min(self.param.num_selected_mols, len(ids))
        return ids[:n] if self.param.select_mode == "first" else self.rng.sample(ids, n)

    def _weights(self, k: int) -> List[float]:
        d = self.param.dist.lower()
        if d == "dirichlet": return Sample().dirichlet(k, self.param.alpha, self.rng)
        if d == "lognormal": return Sample().lognormal(k, self.param.mu, self.param.sigma, self.rng)
        if d == "zipf":      return Sample().zipf_like(k, self.param.s)
        if d == "uniform":   return [1.0]*k
        raise ValueError(f"unknown dist: {self.param.dist}")

    def build(self) -> Tuple[List[str], List[int]]:
        all_ids = self._read_ids()
        sel = self._select(all_ids)
        w = self._weights(len(sel))
        counts = largest_remainder_to_int_counts(w, self.param.num_total_mols)

        # enforce min_count if needed
        if self.param.min_count > 0:
            k = len(counts)
            delta = 0
            for i in range(k):
                if counts[i] < self.param.min_count:
                    delta += (self.param.min_count - counts[i])
                    counts[i] = self.param.min_count
            # subtract 'delta' from largest bins to keep the sum
            if delta > 0:
                order = sorted(range(k), key=lambda i: counts[i], reverse=True)
                j = 0
                while delta > 0 and j < k:
                    idx = order[j % k]
                    if counts[idx] > 0:
                        counts[idx] -= 1
                        delta -= 1
                    j += 1

        # sort output if requested
        if self.param.sort_output_by == "count":
            order = sorted(range(len(sel)), key=lambda i: counts[i], reverse=True)
            sel = [sel[i] for i in order]
            counts = [counts[i] for i in order]
        elif self.param.sort_output_by == "id":
            order = sorted(range(len(sel)), key=lambda i: sel[i])
            sel = [sel[i] for i in order]
            counts = [counts[i] for i in order]

        # assert exact sum
        assert sum(counts) == self.param.num_total_mols
        return sel, counts

    def write(self, ids: List[str], counts: List[int]) -> None:
        with open(self.param.output_tsv, "w", encoding="utf-8") as f:
            f.write("transcript\tcount\n")
            for t, c in zip(ids, counts):
                f.write(f"{t}\t{c}\n")
        print(f"[OK] Wrote {self.param.output_tsv}  (rows={len(ids)}, total={sum(counts)})")
        if self.param.subset_ids_txt:
            with open(self.param.subset_ids_txt, "w", encoding="utf-8") as f:
                for t in ids:
                    f.write(f"{t}\n")
            print(f"[OK] Wrote {self.param.subset_ids_txt}")


class FromGTF:

    def __init__(
            self,
            param: AbundanceParams,
    ):
        self.param = param
        self.rng = random.Random(param.seed)

    # Scan GTF: collect transcript-level ENST and map base -> highest version
    def _collect_ids_and_version_map(self) -> Tuple[List[str], Dict[str, str]]:
        tid_re = re.compile(r'transcript_id "([^"]+)"')
        bio_re = re.compile(r'(?:transcript_biotype|gene_type|gene_biotype) "([^"]+)"')

        # Ordered list of base IDs
        seen_order: List[str] = []
        seen_set = set()
        best: Dict[str, Tuple[int, str]] = {}  # base -> (ver, full)

        with open_maybe_gz(self.param.gtf_fpn, "rt") as f:
            for line in f:
                if not line or line.startswith("#"):
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 9 or cols[2] != "transcript":
                    continue
                attr = cols[8]
                m_tid = tid_re.search(attr)
                if not m_tid:
                    continue
                full = m_tid.group(1)  # ENSTxxxx or ENSTxxxx.y
                if "." in full:
                    base, ver_str = full.rsplit(".", 1)
                    ver = int(ver_str) if ver_str.isdigit() else 0
                else:
                    base, ver = full, 0

                if not self.param.include_noncoding:
                    mb = bio_re.search(attr)
                    bt = mb.group(1) if mb else ""
                    if bt and bt != "protein_coding":
                        continue

                if base not in seen_set:
                    seen_set.add(base)
                    seen_order.append(base)

                rec = best.get(base)
                if rec is None or ver > rec[0]:
                    best[base] = (ver, full)

        if not seen_order:
            raise RuntimeError("No transcript_id parsed from GTF with current filters.")

        if not self.param.strip_version:
            # If pure needs to retain the version, replace the ordered list
            # with the full entry of the highest corresponding version
            seen_order = [best[b][1] for b in seen_order if b in best]
        return seen_order, {b: best[b][1] for b in best}

    def _select(self, ids: List[str]) -> List[str]:
        n = min(self.param.num_selected_mols, len(ids))
        return ids[:n] if self.param.select_mode == "first" else self.rng.sample(ids, n)

    def _weights(self, k: int) -> List[float]:
        d = self.param.dist.lower()
        if d == "dirichlet": return Sample().dirichlet(k, self.param.alpha, self.rng)
        if d == "lognormal": return Sample().lognormal(k, self.param.mu, self.param.sigma, self.rng)
        if d == "zipf":      return Sample().zipf_like(k, self.param.s)
        if d == "uniform":   return [1.0]*k
        raise ValueError(f"unknown dist: {self.param.dist}")

    def build(self) -> Tuple[List[str], List[int]]:
        all_ids, ver_map = self._collect_ids_and_version_map()
        sel = self._select(all_ids)
        w = self._weights(len(sel))
        counts = largest_remainder_to_int_counts_gtf(w, self.param.num_total_mols)

        # min_count constraints
        if self.param.min_count > 0:
            k = len(counts)
            need = 0
            for i in range(k):
                if counts[i] < self.param.min_count:
                    need += (self.param.min_count - counts[i])
                    counts[i] = self.param.min_count
            if need > 0:
                order = sorted(range(k), key=lambda i: counts[i], reverse=True)
                j = 0
                while need > 0:
                    idx = order[j % k]
                    if counts[idx] > self.param.min_count:
                        counts[idx] -= 1
                        need -= 1
                    j += 1

        if self.param.sort_output_by == "count":
            order = sorted(range(len(sel)), key=lambda i: counts[i], reverse=True)
            sel = [sel[i] for i in order]
            counts = [counts[i] for i in order]
        elif self.param.sort_output_by == "id":
            order = sorted(range(len(sel)), key=lambda i: sel[i])
            sel = [sel[i] for i in order]
            counts = [counts[i] for i in order]

        assert sum(counts) == self.param.num_total_mols
        # Cache version mapping for write use.
        self._last_version_map = ver_map
        return sel, counts

    def write(self, ids: List[str], counts: List[int]) -> None:
        # 1) Write pure (no header, Unix line endings).
        with open(self.param.output_tsv, "wb") as f:
            for t, c in zip(ids, counts):
                f.write(f"{t}\t{c}\n".encode("utf-8"))
        print(f"[OK] Wrote {self.param.output_tsv}  (rows={len(ids)}, total={sum(counts)})")

        # 2) Write versioned (use GTF base-to-highest-version mapping)
        ver_map: Dict[str, str] = getattr(self, "_last_version_map", {})
        upgraded: List[Tuple[str, int]] = []
        missing = []
        for t, c in zip(ids, counts):
            # No version or version-stripped
            base = t.split(".", 1)[0]
            full = ver_map.get(base)
            if not full:
                missing.append(t)
            else:
                upgraded.append((full, c))
        if missing:
            # Very rare: only occurs if the pure ID was not sampled from this GTF; issue a strong warning here
            ex = ", ".join(missing[:10])
            raise RuntimeError(f"{len(missing)} ids missing version in current GTF: {ex}")

        with open(self.param.versioned_out, "wb") as f:
            for t, c in upgraded:
                f.write(f"{t}\t{c}\n".encode("utf-8"))
        print(f"[OK] Wrote {self.param.versioned_out}")

        # 3) Optional: write ID list (aligned with the fasta version)
        if self.param.subset_ids_txt:
            with open(self.param.subset_ids_txt, "w", encoding="utf-8") as f:
                for t in ids:
                    f.write(f"{t}\n")
            print(f"[OK] Wrote {self.param.subset_ids_txt}")


if __name__ == "__main__":
    # param_fa = AbundanceParams(
    #     fasta_fpn="D:/Document/Programming/Python/umiche/umiche/data/r1/cm/gencode.v48.transcripts.fa.gz",
    #     output_tsv="D:/Document/Programming/Python/umiche/umiche/data/r1/cm/tksm/abundance.tsv",
    #     num_selected_mols=50,
    #     select_mode="random",
    #     num_total_mols=10000,
    #     seed=1,
    #
    #     use_pipe_if_present=True,
    #     pipe_index=0,
    #     # 0=ENST; if you want to make a table by gene ID, you can set 1=ENSG; if you
    #     # want to make a table by gene name, you can set 5
    #     strip_version=True,
    #     # ENST0000xxx.1 -> ENST0000xxx
    #
    #     dist="dirichlet",
    #     alpha=1.0, # lognormal: mu/sigma；zipf: s
    #     min_count=0,
    #     sort_output_by="count",
    #     subset_ids_txt="D:/Document/Programming/Python/umiche/umiche/data/r1/cm/tksm/refs.sub.ids.fa.txt",
    # )
    # gen = FromFasta(param_fa)
    # ids, counts = gen.build()
    # gen.write(ids, counts)
    # print("Top5:")
    # for i in range(min(5, len(ids))):
    #     print(ids[i], counts[i])

    param_gtf = AbundanceParams(
        gtf_fpn="D:/Document/Programming/Python/umiche/umiche/data/r1/cm/tksm/gencode.v48.annotation.gtf",
        num_selected_mols=5,
        select_mode="random",
        num_total_mols=10,
        strip_version=True,
        dist="dirichlet",
        alpha=1.0,
        min_count=0,
        sort_output_by="count",
        include_noncoding=False,
        output_tsv="D:/Document/Programming/Python/umiche/umiche/data/r1/cm/tksm/abundance.pure.tsv",
        versioned_out="D:/Document/Programming/Python/umiche/umiche/data/r1/cm/tksm/abundance.versioned.tsv",
        subset_ids_txt="D:/Document/Programming/Python/umiche/umiche/data/r1/cm/tksm/refs.sub.ids.gtf.txt",
        seed=1,
    )
    gen = FromGTF(param_gtf)
    ids, counts = gen.build()
    gen.write(ids, counts)
    print("Top5:")
    for i in range(min(5, len(ids))):
        print(ids[i], counts[i])