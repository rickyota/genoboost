FROM --platform=linux/amd64 rust:1.68 AS builder

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    clang cmake libgomp1

WORKDIR /opt/genoboost
COPY ../ .

#export RUSTFLAGS='-C target-cpu=native'
RUN RUSTFLAGS='-C target-cpu=native' \
    cargo build \
    --release \
    --verbose \
    --manifest-path ./projects_rust/Cargo.toml \
    --bin genoboost

#CMD ["./projects_rust/target/genoboost"]

FROM --platform=linux/amd64 debian:bookworm-slim AS runner
#FROM --platform=linux/amd64 debian:bullseye-stable AS runner
#FROM --platform=linux/amd64 debian:buster-slim AS runner

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    libgomp1

WORKDIR /opt/genoboost
# TODO: copy only necessary files
COPY --from=builder /opt/genoboost/projects_rust/target/release/genoboost ./
ENTRYPOINT ["/opt/genoboost/genoboost"]
#ENTRYPOINT ["./genoboost"]

# https://stackoverflow.com/questions/73037618/glibc-incompatibility-on-debian-docker