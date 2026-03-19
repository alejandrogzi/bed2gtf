# syntax=docker/dockerfile:1

# ---------- Build Stage ----------
FROM rust:1.93.0-bullseye AS builder

COPY Cargo.lock Cargo.toml /app/
COPY src/ /app/src/

RUN cargo build --release --manifest-path /app/Cargo.toml && \
    strip /app/target/release/bed2gtf

# ---------- Runtime Stage ----------

FROM debian:bullseye

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    ca-certificates \
    procps \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /app/target/release/bed2gtf /usr/local/bin/

RUN bed2gtf --help
