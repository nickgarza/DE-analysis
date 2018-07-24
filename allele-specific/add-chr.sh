awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' sample.vcf > sample_chr.vcf
