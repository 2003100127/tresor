__version__ = "0.0.1"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__ = "jianfeng.sunmt@gmail.com"

from typing import Dict, List, Tuple, Optional


class TPM:
    def __init__(
        self,
        counts_tsv: str,
        out_tsv: str = "abundance.tpm.tsv",
        cell_barcode: str = "ATGCATGCATGCATGC",
        include_header: bool = False,
        float_format: str = "{:.6f}",
    ):
        """
        Parameters
        ----------
        counts_tsv : str
            Path to input TSV with two columns: transcript_id<TAB>count.
        out_tsv : str
            Path to output TSV with three columns:
            transcript_id<TAB>tpm<TAB>cell-barcode.
        cell_barcode : str
            The barcode string to write in the third column (single cell example).
        include_header : bool
            If True, writes a header line: 'transcript_id\ttpm\tcell-barcode'.
        float_format : str
            Format string for TPM values (default 6 decimals).
        """
        self.counts_tsv = counts_tsv
        self.out_tsv = out_tsv
        self.cell_barcode = cell_barcode
        self.include_header = include_header
        self.float_format = float_format

    def convert_from_cnt(self) -> None:
        records = self._read_counts_table()
        tpm_rows = self._counts_to_tpm(records)
        self._write_tpm(tpm_rows)

    def _read_counts_table(self) -> List[Tuple[str, float]]:
        """
        Read a 2-column table: transcript_id<TAB>count.
        - Skips header if the 2nd field on the first row is not numeric.
        - Trims BOM/CR and whitespace.
        """
        out: List[Tuple[str, float]] = []
        with open(self.counts_tsv, "r", encoding="utf-8", newline="") as f:
            first = True
            for raw in f:
                line = raw.replace("\ufeff", "").rstrip("\n").rstrip("\r")
                if not line:
                    continue
                parts = line.split("\t")
                if len(parts) < 2:
                    parts = line.split()
                    if len(parts) < 2:
                        continue
                tid, val = parts[0].strip(), parts[1].strip()
                # If first row is a header (non-numeric second field), skip it
                if first and not self._is_number(val):
                    first = False
                    continue
                first = False
                try:
                    c = float(val)
                except ValueError:
                    raise ValueError(f"Non-numeric count: {val!r} (line={line!r})")
                if c < 0:
                    # Negative counts are ignored; adjust to your policy if needed
                    continue
                out.append((tid, c))
        if not out:
            raise RuntimeError("No valid rows found in counts table.")
        return out

    def _counts_to_tpm(self, rows: List[Tuple[str, float]]) -> List[Tuple[str, float]]:
        """
        Simple TPM normalization: TPM_i = count_i / Σ(count) * 1e6.
        Returns a list of (transcript_id, tpm).
        """
        total = sum(c for _, c in rows)
        if total <= 0:
            raise RuntimeError("Sum of counts is zero; cannot compute TPM.")
        scale = 1_000_000.0
        return [(tid, c / total * scale) for tid, c in rows]

    def _write_tpm(self, rows: List[Tuple[str, float]]) -> None:
        """
        Write three columns: transcript_id<TAB>tpm<TAB>cell-barcode.
        Uses binary mode to force Unix newlines on all platforms.
        """
        with open(self.out_tsv, "wb") as w:
            if self.include_header:
                w.write(b"transcript_id\ttpm\tcell-barcode\n")
            for tid, tpm in rows:
                tpm_str = self.float_format.format(tpm)
                w.write(f"{tid}\t{tpm_str}\t{self.cell_barcode}\n".encode("utf-8"))

    @staticmethod
    def _is_number(s: str) -> bool:
        try:
            float(s)
            return True
        except Exception:
            return False


if __name__ == "__main__":
    tpm = TPM(
        counts_tsv="D:/Document/Programming/Python/umiche/umiche/data/r1/cm/tksm/abundance.versioned.tsv",
        out_tsv="D:/Document/Programming/Python/umiche/umiche/data/r1/cm/tksm/abundance.tpm.tsv",
        cell_barcode="ATGCATGCATGCATGC",
        include_header=False,  # set True if you want a header row
        float_format="{:.6f}"
    )
    tpm.convert_from_cnt()
