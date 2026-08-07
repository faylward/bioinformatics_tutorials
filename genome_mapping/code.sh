

minimap2 -cx asm10 --cs APM.fasta APMM4.fasta > mimi_out.paf
sort -k6,6 -k8,8n mimi_out.paf | paftools.js call -f AcanthamoebaPolyphagaMimivirus.fasta -L10000 -l1000 - > out.vcf

