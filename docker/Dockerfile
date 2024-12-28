
FROM --platform=linux/amd64 rust:1.74 AS builder
# error. old?
#FROM --platform=linux/amd64 rust:1.68 AS builder

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
     libgomp1

WORKDIR /opt/genoboost
COPY ../ .

RUN RUSTFLAGS='-C target-cpu=native' \
    cargo build \
    --release \
    --verbose \
    --manifest-path ./projects_rust/Cargo.toml \
    --bin genoboost


FROM --platform=linux/amd64 debian:bookworm-slim AS runner

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    libgomp1

# bind $(pwd):/opt/ on local
WORKDIR /opt/
COPY --from=builder /opt/genoboost/projects_rust/target/release/genoboost /app/
ENTRYPOINT ["/app/genoboost"]

# https://stackoverflow.com/questions/73037618/glibc-incompatibility-on-debian-docker
