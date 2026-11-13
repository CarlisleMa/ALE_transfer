FROM python:3.12-slim
ENV DEBIAN_FRONTEND=noninteractive PIP_NO_CACHE_DIR=1 PIP_DISABLE_PIP_VERSION_CHECK=1 \
    MAMBA_ROOT_PREFIX=/opt/mamba LC_ALL=C.UTF-8
RUN apt-get update && apt-get install -y --no-install-recommends \
      bzip2 ca-certificates curl file git less procps tmux unzip xz-utils \
    && rm -rf /var/lib/apt/lists/*
RUN curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/2.0.5 \
    | tar -xj -C /usr/local bin/micromamba
RUN micromamba create -y -p /opt/tools -c conda-forge -c bioconda python=3.12 numpy scipy pandas numpy=2.3.5 scipy=1.16.3 pandas=2.3.3 scikit-learn=1.7.2 flowio=1.4.0 six=1.17.0 pip \
    && micromamba clean -a -y
ENV PATH=/opt/tools/bin:$PATH
# login shells (tmux, `bash -l`) reset PATH in /etc/profile: keep the task's tools visible there too
RUN printf 'export PATH=/opt/tools/bin:$PATH\n' > /etc/profile.d/00-task-tools.sh
RUN pip install flowkit==1.3.2
RUN pip install --no-deps fcsparser==0.2.8
WORKDIR /task
COPY input/ /task/input/
COPY software/ /task/software/
RUN mkdir -p /task/output && chmod -R a-w /task/input && chmod -R a+rx /task/software
