#!/bin/bash
if ! command -v conda &> /dev/null
then
    echo "conda could not be found, please install this first";
    exit
fi
echo "Have you already downloaded SeqWho from Github? [Y/N]: "
read resp1;
while [ $resp1 != "Y" ] && [ $resp1 != "N" ];
do
echo "Have you already downloaded SeqWho from Github? [Y/N]: "
read resp1;
done

if [ $resp1 == "Y" ]; then 
	seqwhodir="";
	echo "Please input the absolute path to the directory of seqwho [^C to exit]: (note this will not work with relative file paths)";
	read seqwhodir; 
	while [[ ! -d $seqwhodir ]]; do
	echo "That is not a directory (or absolute path) on your computer, Please input the directory of seqwho [^C to exit]: "; 
	read seqwhodir;
	done 
	while [[ ! -f "$seqwhodir/seqwho.py" ]] || [[ ! -f "$seqwhodir/seqwho_buildindex.py" ]] || [[ ! -f "$seqwhodir/base_conda.yml" ]] || [[ ! -f "$seqwhodir/full_conda.yml" ]] || [[ ! -f "$seqwhodir/seqwho_modules.py" ]]; do 
	echo "The directory does not contain valid seqwho scripts, Please input the directory of seqwho [^C to exit]: "; 
	read seqwhodir;
	done
	echo "SeqWho Directory Located: $seqwhodir"; 
fi 
  
if [ $resp1 == "N" ]; then 
	resp2="";
	while [ "$resp2" != "Y" ] && [ "$resp2" != "N" ]; do
	echo "Would you like to now download SeqWho? [Y/N]"
	read resp2; 
	done 
	if [ $resp2 == "N" ]; then 
		echo "Sorry, but you have to either have SeqWho, or be willing to Download it in order to test"; 
		exit; 
	fi 
	if [ $resp2 == "Y" ]; then 
if ! command -v git &> /dev/null
then
    echo "git could not be found, please install this first"
    exit
fi
		installdir="";
		echo "We will download SeqWho from Github, please provide a directory to download to: (note this will not work with relative file paths)"; 
		read installdir;
		while [[ ! -d $installdir ]]; do 
			echo "Invalid directory, please provide alternative: ";
			read installdir; 
		done 
		git clone https://github.com/DaehwanKimLab/seqwho $installdir; 
		seqwhodir="$installdir/seqwho";
	fi
fi

echo "SeqWho install contains valid files at: $seqwhodir"; 
echo "--------------------------------------------------"; 
echo "Testing the SeqWho install requires the download  ";
echo "of SRA files which we curate, the download is     ";
echo "large (4.1 GB compressed/~8 GB uncompressed).     ";
echo "--------------------------------------------------";
echo "Where would you like to download the verification ";
echo " files? [provide directory if verify_seqwho_install_data.tar.gz exists there already it will not be downloaded]: "; 
datadir="";
read datadir;
while [[ ! -d $datadir ]]; do 
	echo "Invalid directory, please use absolute path to specify: ";
	read datadir;
done
if [[ ! -f $datadir/verify_seqwho_install_data.tar.gz ]]; then 
	wget https://zenodo.org/record/5826089/files/verify_seqwho_install_data.tar.gz;
	cp verify_seqwho_install_data.tar.gz "$datadir/";
fi; 
cd "$datadir/"
gunzip verify_seqwho_install_data.tar.gz; 
tar -xvf verify_seqwho_install_data.tar
cd verify_seqwho_install_data; 

echo "";
echo "CONDA Output to Follow, NOTE: if you get CondaValueError here this is okay, it just means you have already created the seqwho environment";
echo "-----------------------------------------------------------------------------------------";

source $(which conda | sed "s/\/bin\/conda/\/etc\/profile.d\/conda.sh/");
conda env create -f $seqwhodir/base_conda.yml > /dev/null
conda init bash;
conda activate seqwho_v1; 

echo "-----------------------------------------------------------------------------------------";
echo "End of CONDA Output";
echo ""; 

echo "Would you like to test the SeqWho index building step (Requires >= 32 GB RAM, might take a while)? [Y/N]"; 
read resp3; 
while [ $resp3 != "Y" ] && [ $resp3 != "N" ]; do 
echo "Invalid response.  Test index building phase? [Y/N/^C to exit]: " 
read resp3;
done; 

if [ $resp3 == "Y" ]; then 
	gunzip build-seqwho-index.tar.gz
	tar -xvf build-seqwho-index.tar
	cd build-seqwho-index
	echo "Creating Index - This will Take a while, progress is reported here: "; 
	$seqwhodir/seqwho_buildindex.py -r human_repeats.txt,mouse_repeats.txt -l train.labs > index.build.log; 
	if [[ -f SeqWho.ix ]]; then 
		echo "Index file succesfully built as SeqWho.ix"; 
	fi
	if [[ ! -f SeqWho.ix ]]; then 
		echo "Index file could not be built, please try to reinstall for building, will use pre-built for call-testing"; 
		resp3="N";
	fi   
	cd ..
fi

if [ $resp3 == "N" ]; then 
	echo "Not building index..." 
fi 

echo "Would you like to test the SeqWho file calling capabilities? [Y/N/^C to exit]: " 
read resp4; 
while [ $resp4 != "Y" ] && [ $resp4 != "N" ]; do 
echo "Invalid response.  Test file calling phase? [Y/N/^C to exit]: " 
read resp4; 
done; 

if [ $resp4 == "N" ]; then 
	echo "Testing Completed"
	if [ $resp3 == "Y" ]; then 
		echo "Index Build Verified for Seqwho in $seqwhodir";
	fi 
	exit; 
fi 
if [ $resp4 == "Y" ]; then 
	#cd "$datadir/verify_seqwho_install_data/";
	gunzip call-test-files.tar.gz 
	tar -xvf call-test-files.tar
	cd call-test-files/ 
	if [ $resp3 == "Y" ]; then 
		$seqwhodir/seqwho.py -x "../build-seqwho-index/SeqWho.ix" -f *.gz > file.call.log; 
	fi
	if [ $resp3 == "N" ]; then 
		$seqwhodir/seqwho.py -x SeqWho_prebuilt.ix -f *.gz > file.call.log; 
	fi 
fi
if [ $resp3 == "Y" ]; then 
	echo "Index Build Verified for Seqwho in $seqwhodir";
fi 
	
	echo "Output from calling verified" 

echo "Testing Completed"; 
exit;
