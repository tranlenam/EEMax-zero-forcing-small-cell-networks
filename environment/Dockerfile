# hash:sha256:7d4d87ec29ba830f81ba0f43dec4388fa9e7f74a66ce3d8ec941488d5ab3471a
FROM registry.codeocean.com/codeocean/matlab:2020b-ubuntu20.04

ARG DEBIAN_FRONTEND=noninteractive
ARG MLM_LICENSE_FILE

COPY postInstall /
RUN /postInstall
