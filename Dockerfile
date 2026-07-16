# syntax=docker/dockerfile:1.7
FROM python:3.12-slim-bookworm@sha256:d50fb7611f86d04a3b0471b46d7557818d88983fc3136726336b2a4c657aa30b

LABEL org.opencontainers.image.title="nmr2boltz"
LABEL org.opencontainers.image.description="Conservative NEF/NMR-STAR proton-restraint to Boltz heavy-atom constraint converter"
LABEL org.opencontainers.image.version="0.1.0-validated"

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    PIP_NO_CACHE_DIR=1

WORKDIR /opt/nmr2boltz

COPY requirements.lock pyproject.toml README.md LICENSE ./
COPY src ./src
COPY nmr_to_boltz.py ./

RUN python -m pip install --upgrade pip==26.1.2 setuptools==83.0.0 wheel==0.47.0 packaging==26.2 \
    && python -m pip install -r requirements.lock \
    && python -m pip install --no-build-isolation --no-deps .

WORKDIR /work
USER 65532:65532
ENTRYPOINT ["nmr2boltz"]
CMD ["--help"]
