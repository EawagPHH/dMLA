# Blasting of probes to whole genome sequences
1. Open iTerm
2. Navigated to documents in my laptop, in the wastewater folder, and downloaded blast in there. 
2. Download last version of blast
	
	```
	wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
	```
	Add the name of the last version for the system that you are using after the `/`. The last software version is available at [this link](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/).
3. Unzip the downloaded folder containing blast:
	
	```
	tar -xzf ncbi-blast-2.13.0+-x64-linux.tar.gz
	```
4. Export the path of where blast was downloaded. The code was typed in the directory `Wastewater` where all the `.fasta`files of the WGS results were copied:
	
	```
	export PATH=$PATH:/./ncbi-blast-2.13.0+/bin
	```
5. I prepared a script named `combine_fasta.sh` which allows to combine all the `.fasta` files corresponding to my whole genome sequences samples in a single file, but adding the name of the sample in front of each node. 
6. I saved the script and made it executable in the wastewater folder. 
8. I run the script from that directory by typing:

	```
	./combine_fasta.sh *.fasta
	```
9. Only after having combined the `.fasta`files, I added to the wastewater folder also the `fasta`file containing the probe-sequences. 
10. I format the combined file as a BLAST database:
	
	```
	makeblastdb -in combined_samples.fasta -dbtype nucl -out samples_db
	```
11. I run BLAST search for each probe against the sample database:
	
	```
	blastn -query probes.fasta -db samples_db -out results.out -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" -perc_identity 95
	```
12. I opened R software and I imported the result named `results.out`. 
13. I generated the file containing the results and I saved it:
	
	```
	blast_results <- read_delim("results.out", delim = "\t", 
                            col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"))
	filtered_results <- blast_results %>%
  filter(pident > 95)

	print(filtered_results)

	write_csv(filtered_results, "filtered_blast_results.csv")
	```
14. I opened the file in excel and I modified the cells to have the name of the isolate only, as well as the name of the probes corresponding to the names of the probes that resulted from dMLA. 
	
	
# How to go from numbers to barcodes in excel file
This method allows you to merge the two Excel files based on the barcode column and add the X2 values from the second file into the first file without using any external tools or programming languages.

1. Open the first Excel file where you want to add the X2 column.
2. Open the second Excel file which contains the barcode and X2 columns.
3. Select the range containing the barcode and X2 columns in the second file.
4. Copy the data (CTRL+C or Command+C).
5. Go back to the first file.
6. Right-click on an empty area (e.g., a new sheet) and select "Paste Special" > "Paste Values".
7. Suppose the barcode column in your first file is column A and you have pasted the barcode and X2 data starting in column D.
In the first empty cell of a new column where you want to add X2 (e.g., column B), enter the following formula:

	```
	=VLOOKUP(A2, $D$2:$E$1000, 2, FALSE)
	```

	Adjust the range `$D$2:$E$1000` as necessary to cover all your data.
8. Drag the formula down to fill the column.


