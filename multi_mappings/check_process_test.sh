
filter_line=$(ps -ax | grep filterline | wc -l)

until [ $filter_line '=' "1" ]; do 
echo "There are filterline processes still running"
filter_line=$(ps -ax | grep filterline | wc -l) ;
done


    filter_line=$(ps -ax | grep filterline | wc -l)

#Example loop
counter=0

until [ $counter -gt 5 ]
do
  echo Counter: $counter
  ((counter++))
done