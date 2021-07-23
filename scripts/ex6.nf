params.index = "all-accessions.txt";  
Channel.fromPath(params.index).splitCsv(header:['col1','col2','col3'])
	.map{ row -> tuple(row.col1, row.col2, row.col3) } 
	.set{ csvChannel }; 

process downloadData {
	input: 
	set col1, col2, col3 from csvChannel;
	
	output: 
	file '*.fastq.gz'

	script: 
	""" 
	fastq-dump -X 1000000 --gzip $col1; 
	"""
}
