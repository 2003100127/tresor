import gzip, io

class FastaChrRenamer:

    def __init__(
            self,
            keep_description=False,
    ):
        self.keep_description = keep_description

    def _open(self, path, mode="rt"):
        return gzip.open(path, mode=mode) if path.endswith(".gz") else open(path, mode=mode, encoding="utf-8", newline="")

    def _map_id(self, raw_id: str) -> str:
        # raw_id means '>' token e.g. '1','10','X','MT'
        if raw_id == "MT": return "chrM"
        if raw_id in [str(i) for i in range(1,23)] + ["X","Y"]:
            return "chr" + raw_id
        return raw_id

    def rename(self, fasta_in_fpn: str, fasta_out_fpn: str):
        with self._open(fasta_in_fpn, "rt") as fi, self._open(fasta_out_fpn, "wt") as fo:
            for line in fi:
                if line.startswith(">"):
                    header = line[1:].rstrip("\r\n")
                    first = header.split()[0]
                    rest = header[len(first):]
                    new_id = self._map_id(first)
                    if self.keep_description:
                        fo.write(f">{new_id}{rest}\n")
                    else:
                        fo.write(f">{new_id}\n")
                else:
                    fo.write(line)


if __name__ == "__main__":
    renamer = FastaChrRenamer(
        keep_description=False,
    )
    renamer.rename(
        fasta_in_fpn="D:/Document/Programming/Python/umiche/umiche/data/r1/cm/tksm/GRCh38.dna.primary_assembly.fa.gz",
        fasta_out_fpn="D:/Document/Programming/Python/umiche/umiche/data/r1/cm/tksm/GRCh38.ucsc_like.fa.gz",
    )