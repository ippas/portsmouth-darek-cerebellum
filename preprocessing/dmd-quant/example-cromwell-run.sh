#!/usr/bin/env bash

set -ex

WORKFLOW_URL=rna-seq-paired-end-workflow-cuff-only.wdl


cromwell run $WORKFLOW_URL \
	--inputs inputs-2.json \
	--options cromwell-options.json
