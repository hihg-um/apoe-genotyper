# SPDX-License-Identifier: GPL-2.0

ORG_NAME ?= hihg-um
OS_BASE ?= ubuntu
OS_VER ?= 24.04

IMAGE_REPOSITORY ?=

TOOLS := apoe-genotyper

DOCKER_BUILD_ARGS ?=
DOCKER_TAG ?= $(shell git describe --tags --broken --dirty --all --long | \
		sed "s,heads/,," | sed "s,tags/,," | \
		sed "s,remotes/pull/.*/,," \
		)_$(shell uname -m)_$(shell uname -s | \
		tr '[:upper:]' '[:lower:]')
DOCKER_BASE ?= $(patsubst docker-%,%,$(shell basename \
		`git remote --verbose | grep origin | grep fetch | \
		cut -f2 | cut -d ' ' -f1` | sed 's/.git//'))
DOCKER_IMAGES := $(TOOLS:=\:$(DOCKER_TAG))
SIF_IMAGES := $(TOOLS:=_$(DOCKER_TAG).sif)

IMAGE_TEST := /test.sh

.PHONY: apptainer_clean apptainer_test \
	docker_clean docker_test docker_release $(TOOLS)

help:
	@echo "Targets: all build clean test release"
	@echo "         docker docker_clean docker_test docker_release"
	@echo "         apptainer apptainer_clean apptainer_test"
	@echo
	@echo "Docker container(s):"
	@for f in $(DOCKER_IMAGES); do \
		printf "\t$$f\n"; \
	done
	@echo
	@echo "Apptainer(s):"
	@for f in $(SIF_IMAGES); do \
		printf "\t$$f\n"; \
	done
	@echo

all: clean build test

build: docker apptainer

clean: apptainer_clean docker_clean

release: docker_release

test: docker_test apptainer_test

# Docker
docker: $(TOOLS)

$(TOOLS):
	@echo "Building Docker container: $(ORG_NAME)/$@:$(DOCKER_TAG)"
	@docker build \
		-t $(ORG_NAME)/$@:$(DOCKER_TAG) \
		$(DOCKER_BUILD_ARGS) \
		--build-arg BASE_IMAGE=$(ORG_NAME)/$(DOCKER_BASE):$(DOCKER_TAG) \
		--build-arg RUN_CMD=$@ \
		.

	$(if $(shell git fetch 2>&1; git diff @{upstream} 2>&1),,docker \
		tag $(ORG_NAME)/$@:$(DOCKER_TAG) $(ORG_NAME)/$@:latest)

docker_clean:
	@for f in $(TOOLS); do \
		docker rmi -f $(ORG_NAME)/$$f:$(DOCKER_TAG) 2>/dev/null; \
		if [ -z "`git fetch 2>&1; \
			git diff @{upstream} 2>&1`" ]; then \
				docker rmi -f $(ORG_NAME)/$$f:latest \
			       		2>/dev/null; \
		fi \
	done
	@docker rmi -f $(ORG_NAME)/$(DOCKER_BASE):$(DOCKER_TAG) 2>/dev/null
	@docker builder prune -f 2>/dev/null;

docker_test:
	@for f in $(DOCKER_IMAGES); do \
		echo; echo "Testing Docker container: $(ORG_NAME)/$$f"; \
		docker run -t \
			-v /etc/passwd:/etc/passwd:ro \
			-v /etc/group:/etc/group:ro \
			--entrypoint=$(IMAGE_TEST) \
			--user=$(shell echo `id -u`):$(shell echo `id -g`) \
			$(ORG_NAME)/$$f; \
	done

docker_release:
	@for f in $(DOCKER_BASE):$(DOCKER_TAG) $(DOCKER_IMAGES); do \
		docker tag $(ORG_NAME)/$$f \
			$(IMAGE_REPOSITORY)/$(ORG_NAME)/$$f; \
		docker push $(IMAGE_REPOSITORY)/$(ORG_NAME)/$$f; \
	done

# Apptainer
apptainer: $(SIF_IMAGES)

$(SIF_IMAGES):
	@for f in $(DOCKER_IMAGES); do \
		echo "Building Apptainer: $@"; \
		apptainer pull docker-daemon:$(ORG_NAME)/$$f; \
	done

apptainer_clean:
	@for f in $(SIF_IMAGES); do \
		if [ -f "$$f" ]; then \
			printf "Cleaning up Apptainer: $$f\n"; \
			rm -f $$f; \
		fi \
	done
	@apptainer cache clean

apptainer_distclean:
	@rm -f *.sif

apptainer_test: $(SIF_IMAGES)
	@for f in $^; do \
		echo "Testing Apptainer: $$f"; \
		apptainer exec $$f $(IMAGE_TEST); \
	done