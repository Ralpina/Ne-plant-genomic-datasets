# setting the input and output file paths
input_file="45inds.8chr.map"
output_file="45inds.16chr.map"
# setting the starting row and increment
start_row=218750
increment=1
# loop through each line in the input file
line_number=0
while read line; do
    line_number=$((line_number + 1))

    # checking if this line should have the numbers replaced
    if [ $line_number -ge $start_row ]; then
        # calculating the new values for the first two columns
        new_first_col=$(( (line_number - start_row) / start_row + increment + 1 ))
        new_second_col="chr$new_first_col:"
        
        # replacing the first column and the number after "chr"
        modified_line=$(echo "$line" | sed "s/^[0-9]\+\s/""$new_first_col""\t/")
        modified_line=$(echo "$modified_line" | sed "s/chr[0-9]\+:/""$new_second_col""/")

        # appending the modified line to the output file
        echo "$modified_line" >> "$output_file"
    else
        # appending the original line to the output file
        echo "$line" >> "$output_file"
    fi
done < "$input_file"
