#!/bin/bash
# Bash script creating a file, "scr_anim", 
# which is needed in atomeye to export a series of jpg or png files
# Fri Sep 30 23:25:53 MDT 2011 by Keonwook Kang (kwkang@lanl.gov)

##################################################################
# INPUTS
# 
# Quality of jpg/png files
fig_quality=95

# Directory where cfg files are located
config_dir="./runs/ag-penta-nw-view"

# Output file name
ofile="scr_anim"

# Directory where jpg or png files will be created
fig_dir='Pngs'

# Figure file type jpg/png
ffileformat='png'

# End of INPUTS
##################################################################


##################################################################
# Write file
echo "$fig_quality" > $ofile

i=0;
for file in $config_dir/*.cfg*
do
    ii=$(printf "%05d" $i)

    line="$file $fig_dir/$ii.$ffileformat"
    echo "$line" >> $ofile
    i=$[$i+1]
done

