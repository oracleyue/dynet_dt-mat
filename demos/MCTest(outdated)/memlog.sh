#!/bin/bash -e

# To log mem during the running of another script,
# ---------------------------------
# #!/bin/bash
# trap 'kill $(jobs -p)' EXIT
# ./memlog.sh & YOUR_SCRIPT
# ---------------------------------


echo "      date    time $(free -m | grep total | \
	sed -E 's/^    (.*)/\1/g')" | tee memory.log

while true; do
    echo "$(date '+%Y-%m-%d %H:%M:%S') $(free -m | \
	grep Mem: | sed 's/Mem://g')"| tee -a memory.log
    sleep 10
done
