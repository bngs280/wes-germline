#!/bin/bash

# Load config.json file path
CONFIG_FILE="config.json"

# Parse the config.json to check if resume is true
RESUME_FLAG=""
if [ -f "$CONFIG_FILE" ]; then
    RESUME=$(jq -r '.resume' "$CONFIG_FILE")
    if [ "$RESUME" == "true" ]; then
        RESUME_FLAG="-resume"
    fi
fi

# Run Nextflow command with dynamic -resume flag
#nextflow run /app/main.nf -params-file /app/config.json $RESUME_FLAG
nextflow -log /usr/src/app/output/pipeline.log run main.nf -c nextflow.config -params-file config.json $RESUME_FLAG -work-dir /usr/src/app/work -with-report /usr/src/app/output/pipeline_report -with-dag /usr/src/app/output/pipeline_DAG.png -with-trace /usr/src/app/output/pipeline_trace.txt -with-timeline /usr/src/app/output/pipeline_timeline.html && \
chmod -R 777 /usr/src/app/output
