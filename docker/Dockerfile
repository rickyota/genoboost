FROM --platform=linux/amd64 rust:1.63 AS builder

# TODO: install from github url
WORKDIR /opt/genoboost
COPY ../ .
# export RUST_BACKTRACE=full
RUN cargo build --release -p boosting_rust &&\
	cp ./target/release/boosting_rust ./genoboost


FROM --platform=linux/amd64 continuumio/miniconda3 AS runner

WORKDIR /opt/genoboost
# TODO: copy only necessary files
COPY --from=builder /opt/genoboost ./
RUN conda env create --force -n genoboost -f ./etc/env.yml



