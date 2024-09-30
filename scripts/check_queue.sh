ron=$(mq ls  | grep FAILED | awk '{print $1}')
echo $ron
for i in $ron; do rons=$(pq ls | grep $i | awk '{print $1}'); pq modify resources '40:1:xeon40el8:50h' --all -i $rons; done
# wait 1s
# sleep 3
for i in $ron; do rons=$(pq ls | grep $i | awk '{print $1}'); pq resubmit -i $rons;done
