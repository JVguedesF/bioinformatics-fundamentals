FROM python:3.10-slim-bookworm

WORKDIR /app

RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
    wget \
    gcc \
    g++ \
    make \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

COPY utils /app/utils

WORKDIR /app/P2_Codon_Usage_Analysis

COPY P2_Codon_Usage_Analysis/requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY P2_Codon_Usage_Analysis .

CMD ["python", "src/main.py"]