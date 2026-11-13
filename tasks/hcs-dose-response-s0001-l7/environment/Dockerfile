FROM python:3.12-slim
ENV DEBIAN_FRONTEND=noninteractive PIP_NO_CACHE_DIR=1 PIP_DISABLE_PIP_VERSION_CHECK=1 \
    MAMBA_ROOT_PREFIX=/opt/mamba LC_ALL=C.UTF-8
RUN apt-get update && apt-get install -y --no-install-recommends \
      bzip2 ca-certificates curl file git less procps tmux unzip xz-utils \
    && rm -rf /var/lib/apt/lists/*
RUN curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/2.0.5 \
    | tar -xj -C /usr/local bin/micromamba
RUN micromamba create -y -p /opt/tools -c conda-forge -c bioconda python=3.12 numpy scipy pandas numpy=2.0.2 scipy=1.14.1 pandas=2.2.3 pip \
    && micromamba clean -a -y
ENV PATH=/opt/tools/bin:$PATH
# login shells (tmux, `bash -l`) reset PATH in /etc/profile: keep the task's tools visible there too
RUN printf 'export PATH=/opt/tools/bin:$PATH\n' > /etc/profile.d/00-task-tools.sh
RUN pip install torch==2.5.1 --index-url https://download.pytorch.org/whl/cpu \
 && pip install cellpose==3.1.1.2 scikit-image==0.25.2 tifffile==2025.3.30 matplotlib==3.10.1 \
 && rm -rf /opt/tools/lib/python3.12/site-packages/torch/test /opt/tools/lib/python3.12/site-packages/torch/include /opt/tools/lib/python3.12/site-packages/torch/share \
 && find /opt/tools/lib/python3.12/site-packages/torch/bin -type f ! -name torch_shm_manager -delete
ENV CELLPOSE_LOCAL_MODELS_PATH=/opt/cellpose_models MPLBACKEND=Agg NUMBA_CACHE_DIR=/tmp/numba_cache
RUN mkdir -p /opt/cellpose_models \
 && curl -fsSL -o /opt/cellpose_models/nucleitorch_0 https://www.cellpose.org/models/nucleitorch_0 && echo '89ca45e4a45048d5010d29621466b1274b91b3dcb9714dce9f9e90a9a8671303  /opt/cellpose_models/nucleitorch_0' | sha256sum -c - \
 && curl -fsSL -o /opt/cellpose_models/size_nucleitorch_0.npy https://www.cellpose.org/models/size_nucleitorch_0.npy && echo 'a79107c9f284b569dd65f3a7363014d613a3403241c062822f07b6f026b4d3ef  /opt/cellpose_models/size_nucleitorch_0.npy' | sha256sum -c - \
 && curl -fsSL -o /opt/cellpose_models/cyto3 https://www.cellpose.org/models/cyto3 && echo '2dc3087a8abd7da46d1ab0ddd5824639933cc3ff63b382af3fa1939a392db93c  /opt/cellpose_models/cyto3' | sha256sum -c - \
 && curl -fsSL -o /opt/cellpose_models/size_cyto3.npy https://www.cellpose.org/models/size_cyto3.npy && echo '87f91feb48019b4dfb0400e7bbb89aca5873368686815aa7b4431c0faee97fc5  /opt/cellpose_models/size_cyto3.npy' | sha256sum -c - \
 && chmod -R a+rX /opt/cellpose_models \
 && python -c "from cellpose import models; models.Cellpose(gpu=False, model_type='nuclei'); models.Cellpose(gpu=False, model_type='cyto3')"
WORKDIR /task
COPY input/ /task/input/
COPY software/ /task/software/
RUN mkdir -p /task/output && chmod -R a-w /task/input && chmod -R a+rx /task/software
