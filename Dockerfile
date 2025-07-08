FROM python:3.9-slim

LABEL maintainer="markus.marandi@ut.ee"
LABEL description="Oligogenicity Analysis Toolbox - Genomic variant burden analysis"
LABEL version="0.1.1"
LABEL repository="https://github.com/OligoGeneticDiseases/gen-toolbox"

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    DEBIAN_FRONTEND=noninteractive

RUN groupadd -g 1001 appuser && \
    useradd -r -u 1001 -g appuser appuser


RUN mkdir -p /usr/share/man/man1 && \
    apt-get update && \
    apt-get install -y --no-install-recommends \
        openjdk-11-jre-headless \
        && rm -rf /var/lib/apt/lists/* \
        && apt-get clean

WORKDIR /app

COPY requirements.txt .

# Install Python dependencies
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir -r requirements.txt

COPY src/ ./src/
COPY main.py .

COPY src/config/ ./src/config/

RUN chown -R appuser:appuser /app

USER appuser

HEALTHCHECK --interval=30s --timeout=30s --start-period=5s --retries=3 \
    CMD python -c "import sys; sys.exit(0)"

ENTRYPOINT ["python", "main.py"]

CMD ["--help"]
