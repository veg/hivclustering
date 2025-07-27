
for file in *.tsv; do
  filename=$(basename "$file" .tsv)  # Extract the base filename without extension
  echo "{\"$filename\": [" > "${filename}.json"  # Start JSON array for each file
  awk 'BEGIN {FS="\t"; OFS=","}
      NR==1 {for (i=1; i<=NF; i++) header[i] = $i; next}
      {
        print "{";
        for (i=1; i<=NF; i++) {
          printf("\"%s\": \"%s\"", header[i], $i);
          if (i < NF) printf(", ");
        }
        print "},";
      }' "$file" | sed '$ s/,$//' >> "${filename}.json"  # Process TSV, remove last comma
  echo "]}" >> "${filename}.json"  # End JSON array
done

