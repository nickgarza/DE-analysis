#!/bin/bash
for ((f = 1; f <= 22; ++f))
do
		echo "CONVERTING CHR$f.................";
	   	plink --file chr"$f" --recode vcf --out chr"$f"_converted;
done
