"""Agent / verifier images of the hcs_dose_response family.

The agent image carries Cellpose 3 (CPU-only PyTorch) with the pretrained `nuclei` and `cyto3`
weights baked in at build time, because the agent works without network access. The weights are
fetched from the Cellpose model server and verified against pinned sha256 sums.
"""
from __future__ import annotations

from ... import dockerlib

CONDA = ["numpy=2.0.2", "scipy=1.14.1", "pandas=2.2.3", "pip"]
TORCH = "torch==2.5.1"
PIP = ["cellpose==3.1.1.2", "scikit-image==0.25.2", "tifffile==2025.3.30", "matplotlib==3.10.1"]
WEIGHTS = {   # file -> sha256 (https://www.cellpose.org/models/<file>, fetched 2026-11-11)
    "nucleitorch_0": "89ca45e4a45048d5010d29621466b1274b91b3dcb9714dce9f9e90a9a8671303",
    "size_nucleitorch_0.npy": "a79107c9f284b569dd65f3a7363014d613a3403241c062822f07b6f026b4d3ef",
    "cyto3": "2dc3087a8abd7da46d1ab0ddd5824639933cc3ff63b382af3fa1939a392db93c",
    "size_cyto3.npy": "87f91feb48019b4dfb0400e7bbb89aca5873368686815aa7b4431c0faee97fc5",
}
MODEL_DIR = "/opt/cellpose_models"
SP = "/opt/tools/lib/python3.12/site-packages"


def agent_dockerfile() -> str:
    fetch = " \\\n && ".join(
        f"curl -fsSL -o {MODEL_DIR}/{f} https://www.cellpose.org/models/{f} "
        f"&& echo '{sha}  {MODEL_DIR}/{f}' | sha256sum -c -" for f, sha in WEIGHTS.items())
    extra = (
        f"RUN pip install {TORCH} --index-url https://download.pytorch.org/whl/cpu \\\n"
        f" && pip install {' '.join(PIP)} \\\n"
        # C++ test binaries, headers and cmake files of PyTorch are not needed for inference (-176 MB)
        f" && rm -rf {SP}/torch/test {SP}/torch/include {SP}/torch/share \\\n"
        f" && find {SP}/torch/bin -type f ! -name torch_shm_manager -delete\n"
        f"ENV CELLPOSE_LOCAL_MODELS_PATH={MODEL_DIR} MPLBACKEND=Agg NUMBA_CACHE_DIR=/tmp/numba_cache\n"
        f"RUN mkdir -p {MODEL_DIR} \\\n && {fetch} \\\n"
        f" && chmod -R a+rX {MODEL_DIR} \\\n"
        " && python -c \"from cellpose import models; models.Cellpose(gpu=False, model_type='nuclei'); "
        "models.Cellpose(gpu=False, model_type='cyto3')\"\n"
    )
    return dockerlib.agent_dockerfile(conda=CONDA, extra=extra)


def verifier_dockerfile() -> str:
    return dockerlib.verifier_dockerfile(pip=[])
