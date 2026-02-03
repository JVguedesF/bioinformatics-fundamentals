FROM python:3.10-slim-bookworm

WORKDIR /app

RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
    pymol \
    wget \
    gcc \
    g++ \
    make \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

COPY utils /app/utils

WORKDIR /app/P1_Structural_Profiling_Viz

COPY P1_Structural_Profiling_Viz/requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY P1_Structural_Profiling_Viz .
CMD ["python", "src/main.py"]
