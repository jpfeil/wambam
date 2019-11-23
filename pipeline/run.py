#!/usr/bin/env python3

import argparse
import os
import subprocess
import uuid
import tempfile
import shutil
import sys
import numpy as np


def run(cmd):
    try:
        print(" ".join(cmd))
        subprocess.check_call(cmd)

    except TypeError:
        print("ERROR: Shutting down!")
        print(cmd)
        sys.exit()


def extract_bam_reads(bam, path):
    R1 = os.path.join(path, "R1.fastq")
    R2 = os.path.join(path, "R2.fastq")

    picard = ["java",
              "-jar", "/opt/pipeline/bin/picard.jar",
              "SamToFastq",
              "VALIDATION_STRINGENCY=LENIENT",
              "I=%s" % bam,
              "FASTQ=%s" % R1, 
              "SECOND_END_FASTQ=%s" % R2]

    run(picard)

    R1_size = os.path.getsize(R1)
    R2_size = os.path.getsize(R2)

    if not np.isclose([R1_size], [R2_size], 0.05)[0]:
        R2 = None

    return R1, R2



def get_single_reads(prefix, index, reads, path, CPU):
    bwt = ["bowtie2",
           "-p", str(CPU),
           "-x", index,
           '-S', os.path.join(path, 'wambam.sam'),
           '-U', reads]

    mapped = ["samtools",
              "view",
              "-h",
              "-b",
              "-F", "4",
              "-o", os.path.join(path, "mapped.bam"),
              os.path.join(path, "wambam.sam")]

    sort = ["samtools",
            "sort",
            "-@", str(CPU),
            "-o", os.path.join(path, "mapped.sorted.bam"),
            os.path.join(path, "mapped.bam")]

    pile = ["samtools",
            "mpileup",
            "-A",
            "-a",
            "-o", "%s-pileup" % prefix,
            os.path.join(path, "mapped.sorted.bam")]

    picard = ["java",
              "-jar", "/opt/pipeline/bin/picard.jar",
              "SamToFastq",
              "VALIDATION_STRINGENCY=LENIENT",
              "I=%s" % os.path.join(path, "mapped.sorted.bam"),
              "FASTQ=%s.fastq" % prefix]

    run(bwt)
    run(mapped)
    run(sort)
    run(pile)
    run(picard)


def get_paired_reads(prefix, index, r1, r2, path, CPU):
    bwt = ["bowtie2",
           "-p", str(CPU),
           "-x", index,
           '-S', os.path.join(path, 'wambam.sam'),
           "-1", r1,
           "-2", r2]

    this = ["samtools",
            "view",
            "-b",
            "-F", "4",
            "-f", "8",
            "-o", os.path.join(path, "thisEndMapped.bam"),
            os.path.join(path, "wambam.sam")]

    that = ["samtools",
            "view",
            "-b",
            "-F", "8",
            "-f", "4",
            "-o", os.path.join(path, "thatEndMapped.bam"),
            os.path.join(path, "wambam.sam")]

    both = ["samtools",
            "view",
            "-b",
            "-F", "12",
            "-o", os.path.join(path, "bothEndsMapped.bam"),
            os.path.join(path, "wambam.sam")]

    merge = ["bamtools",
             "merge",
             "-in", os.path.join(path, "thisEndMapped.bam"),
             "-in", os.path.join(path, "thatEndMapped.bam"),
             "-in", os.path.join(path, "bothEndsMapped.bam"),
             "-out", os.path.join(path, "merged.bam")]

    sort = ["samtools",
            "sort",
            "-@", str(CPU),
            "-o", os.path.join(path, "merged.sorted.bam"),
            os.path.join(path, "merged.bam")]

    pile = ["samtools",
            "mpileup",
            '-A',
            "-a",
            "-o", "%s-pileup" % prefix,
            os.path.join(path, "merged.sorted.bam")]

    picard = ["java",
              "-jar", "/opt/pipeline/bin/picard.jar",
              "SamToFastq",
              "VALIDATION_STRINGENCY=LENIENT",
              "I=%s" % os.path.join(path, "merged.sorted.bam"),
              "FASTQ=%s-R1.fastq" % prefix,
              "SECOND_END_FASTQ=%s-R2.fastq" % prefix]

    run(bwt)
    run(this)
    run(that)
    run(both)
    run(merge)
    run(sort)
    run(pile)
    run(picard)


def single_trim(reads, adapter, path, CPU):
    trimmed = os.path.join(path, "trimmed.fq")
    trim = ["java", "-jar",
            "/opt/trim/trimmomatic-0.39.jar",
            "SE",
            "-threads", str(CPU),
            "-phred33",
            reads,
            trimmed,
            "ILLUMINACLIP:/opt/trim/adapters/%s.fa:2:30:10" % adapter,
            "SLIDINGWINDOW:4:5",
            "LEADING:5",
            "TRAILING:5",
            "MINLEN:25"]
    run(trim)
    return trimmed


def paired_trim(r1, r2, adapter, path, CPU):
    p1 = os.path.join(path, "forward_paired.fq")
    p2 = os.path.join(path, "reverse_paired.fq")
    trim = ["java", "-jar",
            "/opt/trim/trimmomatic-0.39.jar",
            "PE",
            "-threads", str(CPU),
            "-phred33",
            r1, r2, p1,
            os.path.join(path, "forward_unpaired.fq"),
            p2,
            os.path.join(path, "reverse_unpaired.fq"),
            "ILLUMINACLIP:/opt/trim/adapters/%s.fa:2:30:10" % adapter,
            "SLIDINGWINDOW:4:5",
            "LEADING:5",
            "TRAILING:5",
            "MINLEN:25"]
    run(trim)
    return p1, p2


def main():
    """
    Extracts reads that align to a given bowtie2 index. Can be used
    to filter reads that map to a particular region of the genome.
    """
    parser = argparse.ArgumentParser(description=main.__doc__)

    parser.add_argument("--index",
                        required=False,
                        help="Path to bowtie2 index")

    parser.add_argument("--faidx",
                        required=False,
                        help="Path to reference index")

    parser.add_argument('--prefix',
                        required=True)

    parser.add_argument('--just-trim',
                        dest="just_trim",
                        action="store_true")

    parser.add_argument('-b', '--bam',
                        help='Path to RNA-seq bam',
                        required=False)

    parser.add_argument('-1', '--R1',
                        help='Path to read 1 fastq',
                        required=False)

    parser.add_argument('-2', '--R2',
                        help='Path to read 2 fastq',
                        required=False)

    parser.add_argument("--adapter",
                        help="""Trimmomatic adapter types: 
                        NexteraPE-PE
                        TruSeq2-PE
                        TruSeq2-SE
                        TruSeq3-PE-2
                        TruSeq3-PE
                        TruSeq3-SE
                        """)

    parser.add_argument('--CPU',
                        default=1,
                        type=str,
                        required=False)

    parser.add_argument('--save',
                        action="store_true")

    args = parser.parse_args()

    _dir = tempfile.gettempdir()
    path = os.path.join(_dir, str(uuid.uuid4()))
    os.mkdir(path)

    if args.bam:
        R1, R2 = extract_bam_reads(args.bam, path) 
        args.R1 = R1
        args.R2 = R2

    elif args.R1 is None:
        raise ValueError("Need at least an input bam for fastq file")

    if args.R2:
        if args.adapter is None:
            adapter = "TruSeq3-PE"
        else:
            adapter = args.adapter

        p1, p2 = paired_trim(args.R1,
                             args.R1,
                             adapter,
                             path,
                             args.CPU)

        if args.just_trim:
            r1 = "%s-trim-R1.fastq" % args.prefix
            shutil.move(p1, r1)
            r2 = "%s-trim-R2.fastq" % args.prefix
            shutil.move(p2, r2)
            subprocess.check_call(["gzip", r1, r2])
            sys.exit(0)


        get_paired_reads(args.prefix,
                         args.index,
                         p1,
                         p2,
                         path,
                         args.CPU)

    else:
        if args.adapter is None:
            adapter = "TruSeq3-SE"
        else:
            adapter = args.adapter

        reads = single_trim(args.R1,
                            adapter,
                            path,
                            args.CPU)

        if args.just_trim:
            r = "%s-trim.fastq" % args.prefix
            shutil.move(reads, r)
            subprocess.check_call(["gzip", r])
            sys.exit(0)

        get_single_reads(args.prefix,
                         args.index,
                         reads,
                         path,
                         args.CPU)

    if args.save:
        print("TEMP FILES: \n%s" % path)
    else:
        shutil.rmtree(path)


if __name__ == "__main__":
    main()
