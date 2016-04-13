# ono
Order and orient

This tool tries to order and orient the contigs of a genome using a reference.
It is based on the libraries of : https://github.com/jasperlinthorst/reveal

Before the specific libraries installation, they are some common libraries and tools you might need to install if you do not have them already :

	sudo apt-get install git cmake zlib1g-dev python-dev python-matplotlib

Install instructions as taken from the release of Reveal used at the time of writing of this tool :

	It also depends on the following 3d party packages:

	  libdivsufsort 2.0.1
	  seqal --> https://github.com/mhulsman/seqal

	1) Build and install libdivsufsort 2.0.1 In order to build and install libdivsufsort you need to have CMAKE (versions should be available for most operating systems). After that's installed do the following:

	  wget https://libdivsufsort.googlecode.com/files/libdivsufsort-2.0.1.tar.gz
	  tar xvf libdivsufsort-2.0.1.tar.gz
	  cd libdivsufsort-2.0.1
	  mkdir build
	  cd build
	  sudo cmake -DBUILD_DIVSUFSORT64:BOOL=ON -DCMAKE_BUILD_TYPE="Release" -DCMAKE_INSTALL_PREFIX="/usr/local" ..
	  sudo make install

	Libdivsufsort should now be installed into your default installation directory (most likely /usr/local/lib).

	2) Build and install seqal Seqal is a python package that is used by ProBubble to perform Needleman-Wunsch alignments of larger bubbles to detect inversions. To build and install it, do the following:

	  git clone https://github.com/mhulsman/seqal.git
	  cd seqal
	  sudo python setup.py install

You can then install ONO by running :

	sudo python setup.py install


Running ONO:

	python assembly_finishing.py --minmum minimal_match_size -n number_of_N_to_insert_between_contigs reference.fasta contigs.fasta
optional :  
	-discard don't append the contigs that couldn't be ordered and orientated to the output.  
	-p prune all sequences of N to a single N (overrides the -n option), even those present in the input contigs.  
	-step is the maximal distance at which the tool will look for neighbourhing MUMs during the cleaning step. All mums that have another mum within step distance will be kept, other will be discarded. This is done in order to ignore isolated mums like the ones circled on the figure below that usually correspond to small duplication (transposable elements) when trying to order and orient contigs.  
	-smallest is similar to step, but is only used in the case that there is only one MUM on that contig. For example, a common transposable element in TB is IS6110 that is usually 1354bp long, meaning it can be the main element in a MUM with the default minmum size of 1000. This option is to discard contigs that have only 1 such MUM. (setting this value to a lower or equal value to minmum will not have any effect).

![Mumplot Example](doc/mumplot_example.png)
