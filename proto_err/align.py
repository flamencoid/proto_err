#!/usr/bin/env python
# Aligning reads to a reference
# Mostly wrappers for samtools scripts

import logging
from Bio.Sequencing import Applications
from Bio.Application import _Option, _Argument, _Switch, AbstractCommandline,_StaticArgument
import os

def refIndex(file):
	"""
	Function to generate BWA index
	"""
	logging.info("Creating BW index of reference")
	index_cmd = Applications.BwaIndexCommandline(infile=file, algorithm="bwtsw")
	index_cmd()
	return 1

def align(reference, read_file, stdout,algorithm='bwa-mem'):
	if algorithm=='bwa-mem':
		logging.info("Aligning reads to reference with bwa-mem")
        alignCmd = BwaMemAlignCommandline( reference=reference, read_file=read_file)
        print alignCmd

	return alignCmd(stdout=stdout)

class BwaMemAlignCommandline(AbstractCommandline):
    """Command line wrapper for Burrows Wheeler Aligner (BWA) aln.

    Run a BWA alignment, equivalent to::

        $ bwa aln [...] <in.db.fasta> <in.query.fq> > <out.sai>

    See http://bio-bwa.sourceforge.net/bwa.shtml for details.

    Example:

    >>> from Bio.Sequencing.Applications import BwaAlignCommandline
    >>> reference_genome = "/path/to/reference_genome.fasta"
    >>> read_file = "/path/to/read_1.fq"
    >>> output_sai_file = "/path/to/read_1.sai"
    >>> read_group="@RG\tID:foo\tSM:bar"
    >>> align_cmd = BwaAlignCommandline(reference=reference_genome, read_file=read_file)
    >>> print(align_cmd)
    bwa aln /path/to/reference_genome.fasta /path/to/read_1.fq

    You would typically run the command line using align_cmd(stdout=output_sai_file)
    or via the Python subprocess module, as described in the Biopython tutorial.
    """
    def __init__(self, cmd="bwa", **kwargs):
        self.program_name = cmd
        self.parameters = \
                [
                    _StaticArgument("mem -H"),
                    _Argument(["reference"], "Reference file name",
                              filename=True, is_required=True),
                    _Argument(["read_file"], "Read file name",
                              filename=True, is_required=True),
                    _Option(["-k", "k"],
                            "Minimum seed length. Matches shorter than INT will be missed. The alignment speed is usually insensitive to this value unless it significantly deviates 20.",
                            checker_function=lambda x: isinstance(x, (int, float)),
                            equate=False),
                    # _Option(["-o", "o"],
                    #         "Maximum edit distance if the value is INT, or the fraction of missing alignments given 2% uniform base error rate if FLOAT. In the latter case, the maximum edit distance is automatically chosen for different read lengths. [0.04]",
                    #         checker_function=lambda x: isinstance(x, (int, float)),
                    #         equate=False),
                    # _Option(["-e", "e"],
                    #         "Maximum number of gap extensions, -1 for k-difference mode (disallowing long gaps) [-1]",
                    #         checker_function=lambda x: isinstance(x, int),
                    #         equate=False),
                    # _Option(["-d", "d"],
                    #         "Disallow a long deletion within INT bp towards the 3-end [16]",
                    #         checker_function=lambda x: isinstance(x, int),
                    #         equate=False),
                    # _Option(["-i", "i"],
                    #         "Disallow an indel within INT bp towards the ends [5]",
                    #         checker_function=lambda x: isinstance(x, int),
                    #         equate=False),
                    # _Option(["-l", "l"],
                    #         """Take the first INT subsequence as seed.

                    #         If INT is larger than the query sequence, seeding will be disabled.
                    #         For long reads, this option is typically ranged from 25 to 35 for
                    #         -k 2. [inf]""",
                    #         checker_function=lambda x: isinstance(x, int),
                    #         equate=False),
                    # _Option(["-k", "k"], "Maximum edit distance in the seed [2]",
                    #         checker_function=lambda x: isinstance(x, int),
                    #         equate=False),
                    # _Option(["-t", "t"], "Number of threads (multi-threading mode) [1]",
                    #         checker_function=lambda x: isinstance(x, int),
                    #         equate=False),
                    # _Option(["-M", "M"],
                    #         "Mismatch penalty. BWA will not search for suboptimal hits with a score lower than (bestScore-misMsc). [3]",
                    #         checker_function=lambda x: isinstance(x, int),
                    #         equate=False),
                    # _Option(["-O", "O"], "Gap open penalty [11]",
                    #         checker_function=lambda x: isinstance(x, int),
                    #         equate=False),
                    # _Option(["-E", "E"], "Gap extension penalty [4]",
                    #         checker_function=lambda x: isinstance(x, int),
                    #         equate=False),
                    # _Option(["-R", "R"],
                    #         """Proceed with suboptimal alignments if there are no more than INT equally best hits.

                    #         This option only affects paired-end mapping. Increasing this threshold helps
                    #         to improve the pairing accuracy at the cost of speed, especially for short
                    #         reads (~32bp).""",
                    #         checker_function=lambda x: isinstance(x, int),
                    #         equate=False),
                    # _Option(["-q", "q"],
                    #         """Parameter for read trimming [0].

                    #         BWA trims a read down to argmax_x{\sum_{i=x+1}^l(INT-q_i)} if q_l<INT
                    #         where l is the original read length.""",
                    #         checker_function=lambda x: isinstance(x, int),
                    #         equate=False),
                    # _Option(["-B", "B"],
                    #         "Length of barcode starting from the 5-end. When INT is positive, the barcode of each read will be trimmed before mapping and will be written at the BC SAM tag. For paired-end reads, the barcode from both ends are concatenated. [0]",
                    #         checker_function=lambda x: isinstance(x, int),
                    #         equate=False),
                    # _Switch(["-c", "c"],
                    #         "Reverse query but not complement it, which is required for alignment in the color space."),
                    _Switch(["-H", "H"],
                            "Hardclipping off"),
                    # _Switch(["-I", "I"],
                    #         "The input is in the Illumina 1.3+ read format (quality equals ASCII-64)."),
                    # _Switch(["-b", "b"],
                    #         "Specify the input read sequence file is the BAM format"),
                    # _Switch(["-b1", "b1"],
                    #         "When -b is specified, only use the first read in a read pair in mapping (skip single-end reads and the second reads)."),
                    # _Switch(["-b2", "b2"],
                    #         "When -b is specified, only use the second read in a read pair in mapping.")
                  ]
        AbstractCommandline.__init__(self, cmd, **kwargs)


