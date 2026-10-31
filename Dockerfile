FROM python:3.12.16-slim-bookworm@sha256:623674a14ef93d052bb9b0de4dbf5f28eb578a35263838941cc8aa43213c1487

ENV PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

WORKDIR /processor

# Pinned Python dependencies, binary wheels only (no compiler, no apt packages), plus the Lambda
# runtime interface client (awslambdaric) so the same image runs as a container or a Lambda function.
COPY requirements.txt /processor/requirements.txt
RUN pip install --only-binary=:all: -r /processor/requirements.txt \
    && find /usr/local/lib/python3.12/site-packages -name '*.pyc' -delete \
    && find /usr/local/lib/python3.12/site-packages -type d -name tests -prune -exec rm -rf {} +

COPY processor/ /processor/processor/
COPY VERSION /processor/VERSION

ENV PYTHONPATH="/processor" \
    PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    MPLBACKEND=Agg \
    MPLCONFIGDIR=/tmp/matplotlib \
    NUMBA_CACHE_DIR=/tmp/numba

# Build matplotlib's font cache now, writable by the run-time user (Lambda's file system is read-only
# outside /tmp), and check that the stack imports.
RUN python -c "import matplotlib.pyplot, scanpy, awslambdaric" && mkdir -p /tmp/numba && chmod -R a+rwX /tmp/matplotlib /tmp/numba

# Runs as uid 1000 on the platform; nothing here needs root. No network is used at run time.
USER 1000

# Runs main() once as a container, or serves processor.handler as a Lambda function (see processor/entry.py).
ENTRYPOINT ["python3", "-m", "processor.entry"]
