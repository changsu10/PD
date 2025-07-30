for file in *_FamaMacBeth_NeweyWest_*_alpha_0.5.txt; do
  # Extract the prefix before "_FamaMacBeth"
  prefix="${file%%_FamaMacBeth*}"
  new_filename="${prefix}.txt"
  mv "$file" "$new_filename"
done