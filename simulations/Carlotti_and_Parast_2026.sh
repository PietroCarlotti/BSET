#!/bin/bash

n_batches=10
limit=10

# Use n_batches here
for ((i=1; i<=n_batches; i++)); do
    # Launch R
    nice R CMD BATCH "--args $i" Carlotti_and_Parast_2026_batch.R "Carlotti_and_Parast_2026_batch_$i.Rout" &
    
    echo "Launched Batch $i"
    
    sleep 3
    
    # Use limit here
    while [ $(jobs -rp | wc -l) -ge $limit ]; do
        sleep 1
    done
done

wait

# Slack notification
Slack_URL="$\{SLACK_WEBHOOK_URL\}"

# Create a nice message
Slack_notification="*BSET Simulations Complete*\nAll $n_batches batches have finished running.\nTime: $(date)"

# Send to Slack
curl -X POST -H 'Content-type: application/json' \
     --data "{\"text\":\"$Slack_notification\"}" \
     $Slack_URL

echo "All batch jobs finished."
