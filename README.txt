OrthoReD is an automated orthology prediction tool. Provided with a gene of
interest (query), a dataset to search from (subject), and information about each
sequence included in the subject, OrthoReD will search through the subject to
find predicted orthologs.

Please cite:
Battenberg K, Lee EK, Chiu JC, Berry AM, Potter D. 2017. OrthoReD: a rapid and
accurate orthology prediction tool with low computational requirement. BMC
Bioinformatics. 18:310.



DEPENDENCIES
====================
OrthoReD requires the following tools to be installed to run properly.

NCBI BLAST (tested version is v2.7.1+)
	Available at: https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download
AB BLAST (tested version is v3.0PE; this is optional)
	Available at: https://www.advbiocomp.com/blast.html
SWIPE (tested version is 2.0.12)
	Available at: https://github.com/torognes/swipe/releases
MCL (tested version is v14-137)
	Available at: http://micans.org/mcl/src/
MAFFT (tested version is v7.402)
	Available at: http://mafft.cbrc.jp/alignment/software/
RAxML (tested version is v8.2.12)
	Available at: http://sco.h-its.org/exelixis/web/software/raxml/index.html
	*for RAxML, make sure to keep track of the version that was installed (AVX or SSE3)
Newick Utilities (tested version is v1.6)
	Available at: http://cegg.unige.ch/newick_utils
====================



INSTALL
====================
OrthoReD is a perl script and does not require any installation. However, before
running OrthoReD, make sure that the all dependencies are properly installed by
testing if each tools are available in the PATH.

  which segmasker #tests NCBI BLAST
  which xdformat #tests AB BLAST
  which swipe #test SWIPE
  which mcl #tests MCL
  which mafft #tests MAFFT
  which raxmlHPC-PTHREADS-AVX #tests RAxML (for AVX)
  which raxmlHPC-PTHREADS-SSE3 #tests RAxML (for SSE3)
  which nw_reroot #tests Newick Utilities

If any of these command do not return a path, add these tools to the path.

  export PATH=$PATH:/some/path/to/the/tool/
====================



SCRIPT OVERVIEW
====================
OrthoReD comes in three parts.
(1) step-01-02.pl
	This script, given fasta file(s) with your gene of interests (query),
	will format all the headers and sequences for orthology prediction.
(2) step-03-04.pl
	This script, given the query and the dataset to search from (subject),
	will format all the headers and sequences for ortholgoy prediction,
	and generate a BLAST database file.
(3) step-05-06.pl
	This script, given the query and the subject, will search through the
	subject for orthologs for each sequence in the query.

*If a properly formatted query and subject are already available, step-05-06.pl
can be executed directly without the first two scripts.
====================



TEST RUN
====================
Try to run OrthoReD using the provided example dataset to check if
OrthoReD is running properly.

1. Make a new empty directory <test>.

  mkdir test

2. Move into the directory

  cd test

3. Copy two directories <OrthoReD_vXXXXXXXX> and <INPUT_EXAMPLE> into
<test>.

  cp -R /some/path/to/OrthoReD_vXXXXXXXX .
  cp -R /some/path/to/INPUT_EXAMPLE .

4. run OrthoReD on the example dataset.

  bash INPUT_EXAMPLE/testrun_EXAMPLE.sh

5. View the generated output. If OrthoReD ran properly, it should
generate a file with the predicted orthologs.

  less ./Step-06_ORTHOLOG/07_ORTHO_DNA/At_rbcL_na_U91966v1_ORTHO_DNA.fas

6. Compare the generated output with the expected output.

  less ./INPUT_EXAMPLE/expected_output.fas

OrthoReD is running properly if these two files are identical.
====================



FILE FORMAT
====================
Once test run is complete, OrthoReD is ready for a user defined file.
Two files need to be prepared with specific formats.

1.Input sequences
The header of each and every sequence provided to OrthoReD has to have the
following format:

  <Species ID>_<Locus ID>_<Isoform ID>_<Some other ID>
  *See "expected_output.fas" as an example.
  *The first two IDs are required, the third ID is recommended but not required,
  and the last ID is not required

  <Species ID>: Short ID for each species.
                (e.g. "Mt" for "Medicago truncatula")

  <Locus ID>: Some ID to specify a given gene/locus within a genome.
              (e.g. "Medtr7g117930")

  <Isoform ID>: Indicate which isoform it is if available. Otherwise put “na”.

  <Some other ID>: Put “na” unless there are some other information that the use
                   needs to keep.

2.Species list
All species included into the database needs to be indicated in a single file
with the following format:

  <Full species name>	<Short species name>	<Species ID>
  *See "species_list_EXAMPLE.txt" as an example.

  <Full species name>: Full name of the species, genus and the specific epithet.
                       (e.g. "Medicago truncatula")

  <Short species name>: The first letter of the genus and the full specific
                        epithet.
                        (e.g. Athaliana)

  <Species ID>: Short ID for each species.
                (e.g. "Mt" for "Medicago truncatula")
                *make sure this is identical to the <Species ID> used in the
                input sequence file.

Note: Once, the headers are formatted, upon running "step-01-02.pl" or
"step-03-04.pl", make sure to turn off the automated header cleaning option
(--q_clean NO --s_clean NO) whenever applicable.
====================
