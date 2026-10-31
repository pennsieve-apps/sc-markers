FROM python:3.12.16-slim-bookworm@sha256:e697c81897302a4c4333bade000d98038fbb21eb064f8ee1fd42d02ebcd9d337

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

# Check that the stack imports. Caches (matplotlib, numba) go to /tmp at run time (MPLCONFIGDIR,
# NUMBA_CACHE_DIR): the only writable place besides OUTPUT_DIR.
RUN python -c "import scanpy, awslambdaric"

# Runs as uid 1000 on the platform; nothing here needs root. No network is used at run time.
USER 1000

# Runs main() once as a container, or serves processor.handler as a Lambda function (see processor/entry.py).
ENTRYPOINT ["python3", "-m", "processor.entry"]
