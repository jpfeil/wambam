#!/usr/bin/env python3

import argparse
import subprocess
import os

def get_reads(prefix, index, r1, CPU, r2=None):
    bwt = ["bowtie2",
           "-p", str(CPU),
           "-x", index]

    if r2 is None:
        bwt.extend(["-U", r1])

    else:
        bwt.extend(["-1", r1, "-2", r2])

    this = ["samtools",
            "view",
            "-b",
            "-F", "4",
            "-f", "8",
            "-",
            "-o", ".thisEndMapped.bam"]

    that = ["samtools",
            "view",
            "-b",
            "-F", "8",
            "-f", "4",
            "-",
            "-o", ".thatEndMapped.bam"]

    both = ["samtools",
            "view",
            "-b",
            "-F", "12",
            "-",
            "-o", ".bothEndsMapped.bam"]

    merge = ["samtools",
             "merge",
             "-1",
             "-c",
             "-p", ".merged.bam",
             ".thisEndmapped.bam",
             ".thatEndmapped.bam",
             ".bothEndsmapped.bam"]

    print(bwt)
    p1 = subprocess.Popen(bwt,
                          stdout=subprocess.PIPE)
    p1.communicate()

    print(this)
    subprocess.Popen(this, stdin=p1.stdout).wait()
    subprocess.Popen(that, stdin=p1.stdout).wait()
    subprocess.Popen(both, stdin=p1.stdout).wait()

    subprocess.check_call(merge)

    os.remove(".thisEndmapped.bam")
    os.remove(".thatEndmapped.bam")
    os.remove(".bothEndsmapped.bam")

    picard = ["java",
              "-jar", "/opt/pipeline/bin/picard.jar",
              "SamToFastq",
              "VALIDATION_STRINGENCY=LENIENT",
              "I=.merged.bam",
              "FASTQ=%s-R1.fastq" % prefix]

    if r2 is not None:
        picard.append("SECOND_END_FASTQ=%s-R2.fastq" % prefix)

    subprocess.check_call(picard)


def main():
    parser = argparse.ArgumentParser(description=main.__doc__)

    parser.add_argument("--index",
                        required=True,
                        help="Path to bowtie2 index")

    parser.add_argument('--prefix',
                        required=True)

    parser.add_argument('-1', '--R1',
                        help='Path to read 1 fastq',
                        required=True)

    parser.add_argument('-2', '--R2',
                        help='Path to read 2 fastq',
                        required=False)

    parser.add_argument('--CPU',
                        default=1,
                        required=False)

    args = parser.parse_args()

    print(args)

    get_reads(args.prefix,
              args.index,
              args.R1,
              args.CPU,
              args.R2)

if __name__ == "__main__":
    main()
