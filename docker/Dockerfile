FROM --platform=linux/amd64 rust:1.68 AS builder

RUN apt-get update &&
	apt-get install -y --no-install-recommends \
		clang

WORKDIR /opt/genoboost
COPY ../ .

#export RUSTFLAGS='-C target-cpu=native'
RUN cargo build \
	--release \
	--manifest-path ./projects_rust/Cargo.toml \
	--bin genoboost

ENTRYPOINT ["./projects_rust/target/genoboost"]

#FROM --platform=linux/amd64 debian:buster-slim AS runner
#
#WORKDIR /opt/genoboost
## TODO: copy only necessary files
#COPY --from=builder /opt/genoboost/projects_rust/target/release/genoboost ./
#CMD ["./genoboost"]
