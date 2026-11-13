"""Agent image for the rnaseq_confound family.

`dockerlib.agent_dockerfile` is not used: with the R/Bioconductor stack the build log of the
`micromamba create` step contained a non-UTF-8 byte (source not tracked down), and
`dockerlib.build_image` decodes the captured log strictly (UnicodeDecodeError on the build host).
The only difference to the standard helper is that the solver/installer output goes to a log file
that is printed (ASCII-filtered) only when the step fails. Base and tail layers are the shared ones.
"""
from __future__ import annotations

from ... import dockerlib

# Shipped stack: Python only. The R/Bioconductor stack below was built once for calibration
# (real DESeq2 / limma-voom / edgeR on the simulated counts); with it the image is 3.8 GB as
# reported by `docker images` on the build host, above the 3 GB budget, so it is not shipped.
CONDA = ["salmon=2.11.0", "pydeseq2=0.5.4", "statsmodels=0.14.6"]
DEV_CONDA = CONDA + [
    "r-base=4.5.3", "bioconductor-deseq2=1.50.2", "bioconductor-edger=4.8.2",
    "bioconductor-limma=3.66.0", "bioconductor-tximport=1.38.2",
]


def agent_dockerfile(conda: list[str] | None = None) -> str:
    pkgs = " ".join(["python=3.12", "numpy", "scipy", "pandas"] + list(conda or CONDA))
    return (
        dockerlib.AGENT_BASE
        + "RUN (micromamba create -y -q -p /opt/tools -c conda-forge -c bioconda "
        + pkgs + " > /tmp/mm.log 2>&1 \\\n"
        "      || { tail -n 80 /tmp/mm.log | tr -cd '\\11\\12\\40-\\176'; exit 1; }) \\\n"
        "    && micromamba clean -a -y > /dev/null 2>&1 && rm -f /tmp/mm.log\n"
        "ENV PATH=/opt/tools/bin:$PATH\n"
        + dockerlib.AGENT_TAIL
    )
