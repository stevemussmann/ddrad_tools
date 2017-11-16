# ddrad_tools
A collection of tools designed for simulating the restriction digest of a genome using two restriction enzymes.

This is my own implementation of the algorithm described in:
Chafin TK, BT Martin, SM Mussmann, MR Douglas, ME Douglas (2017) *FRAGMATIC: in silico locus prediction and its utility in optimizing ddRADseq projects*. pp. 1-4

The paper describing the algorithm is available at https://doi.org/10.1007/s12686-017-0814-1

These scripts are designed to simulate the digest of an organism's genome by two restriction enzymes according to the double-digest restriction site associated DNA (ddRAD) sequencing method (Peterson et al. 2012: http://dx.doi.org/10.1371/journal.pone.0037135).

## Installation:

These scripts are intended to be run on Unix based operating systems, such as the various Linux distributions and Mac OS X.  Windows is not supported at this time.  To get started, copy the desired script (high memory or low memory version) into a folder.  The high memory version reads an entire genome into memory at once, which may overwhelm some laptops and desktop computers, but it may run marginally faster on systems that have adequate memory.

These scripts also depend upon Gnuplot (http://www.gnuplot.info/) and the perl module Chart::Gnuplot (http://search.cpan.org/dist/Chart-Gnuplot/lib/Chart/Gnuplot.pm).  If it is not already installed, you should first install Gnuplot.  This can be done through your Linux distribution's package manager.  On Ubuntu, the following command should work:
```
sudo apt-get install gnuplot
```
Other versions Linux distributions may use a different package manager, so your command may vary slightly.  Installation of Gnuplot can also be performed in this way if you have installed a 3rd-party package manager such as HomeBrew (http://brew.sh/) or MacPorts (https://www.macports.org/).  

Finally, install the Chart::Gnuplot perl module.  Under Ubuntu, you may be able to execute:
```
sudo apt-get install libchart-gnuplot-perl
```

The easiest way to install perl modules on other systems is through CPAN.  Some information on installing packages through CPAN is available here: http://www.cpan.org/modules/INSTALL.html

## Running the script:

Place a FASTA formatted file in the same directory as the script.  To find the number of ddRAD fragments you would expect in the size range of 275 to 325bp in length for a file named genomefile.fasta using the PstI and MspI enzymes, you would use the following command:

```
./ddrad_lowmem.pl -1 PstI -2 MspI -f genomefile.fasta -s 300 -e 25
```

The order of input for restriction enzymes should not influence results (I could have swapped positions of PstI and MspI).  If your genome file is in a directory other than the location of the script, you will have to give its path.  

To see a list of all commands for the program, execute the script without any command line arguments:

```
./ddrad_lowmem.pl
```

## Supported Restriction Enzymes:

A list of supported enzymes is as follows.  Cut sites are indicated by '^'.

	   AciI     C^CGC
	   AgeI     A^CCGGT
	   AluI     AG^CT
	   AseI     AT^TAAT
	   BamHI    G^GATCC
	   BfaI     C^TAG
	   BgIII    A^GATCT
	   BspDI    AT^CGAT
	   ClaI     AT^CGAT
	   DpnII    ^GATC
       EcoRI    G^AATTC
	   EcoRV    GAT^ATC
	   EcoT22I  ATGCA^T
	   HindIII  A^AGCTT
	   KpnI     GGTAC^C
	   MluCI    ^AATT
	   MseI     T^TAA
       MspI     C^CGG
	   NdeI     CA^TATG
	   NheI     G^CTAGC
	   NlaIII   CATG^
	   NotI     GC^GGCCGC
	   NsiI     ATGCA^T
       PstI     CTGCA^G
	   RsaI     GT^AC
	   SacI     GAGCT^C
	   Sau3AI   ^GATC
	   SpeI     A^CTAGT
	   SphI     GCATG^C
	   TaqI     T^CGA
	   XbaI     T^CTAGA
	   XhoI     C^TCGAG
	   SbfI     CCTGCA^GG
     
