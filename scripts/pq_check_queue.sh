ron=$(pq ls  | grep pristine |  grep waiting | awk '{print $1}')
echo $ron
for i in $ron; do  pq modify resources '24:1:xeon24el8:50h' --all -i $i; done

ron=$(pq ls  | grep adsorbate |  grep waiting | awk '{print $1}')
echo $ron
for i in $ron; do  pq modify resources '24:1:xeon24el8:50h' --all -i $i; done

ron=$(pq ls  | grep vibration |  grep waiting | awk '{print $1}')
echo $ron
for i in $ron; do  pq modify resources '24:1:xeon24el8:50h' --all -i $i; done
