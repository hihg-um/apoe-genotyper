# SPDX-License-Identifier: GPL-3.0
# ------------------------------------------------------------------------------
# RUNTIME LAYER
FROM ubuntu:24.04
ARG RUN_CMD

RUN apt -y update -qq && apt -y upgrade && DEBIAN_FRONTEND=noninteractive \
	apt -y install --no-install-recommends \
	python3 python3-cyvcf2

COPY --chmod=0555 scripts/apoe-genotyper.py /usr/local/bin/apoe-genotyper.py

ARG TEST="/test.sh"
COPY --chmod=0555 scripts/test/$RUN_CMD.sh ${TEST}
# ------------------------------------------------------------------------------