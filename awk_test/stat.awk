#http://www.commandlinefu.com/commands/view/1661/display-the-standard-deviation-of-a-column-of-numbers-with-awk
#http://www.cyberciti.biz/faq/bash-scripting-using-awk/
#awk -v max=1 '{print $1*max}' testdata.txt | awk -f stat.awk
{sum+=$1; array[NR]=$1} END {for(x=1;x<=NR;x++){sumsq+=((array[x]-(sum/NR))**2);}print sqrt(sumsq/NR)}
