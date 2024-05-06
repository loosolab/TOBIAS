import glob
import os
import argparse
import sys
import pandas as pd
import subprocess
from tobias.parsers import add_submerge_arguments
from tobias.utils.logger import TobiasLogger


def run_submerge(args):

    logger = TobiasLogger("SubMerge", args.verbosity)
    logger.begin()
    parser = add_submerge_arguments(argparse.ArgumentParser())
    logger.arguments_overview(parser, args)

    outdir = os.path.dirname(args.output)

    # Get all TFBS files
    tfbs_files = glob.glob(os.path.join(args.tfbs, "*/*_overview.txt"))

    # Read regions
    logger.debug("Reading regions from: " + args.regions)
    checked_regions = []
    with open(args.regions) as f:
        for i, line in enumerate(f):
            splitline = line.strip().split("\t")

            # if more than 6 columns, remove the last columns
            if len(splitline) > 6:
                splitline = splitline[:6]

            # check if all necessary columns are present
            while len(splitline) < 6:
                match len(splitline):
                    case 3:
                        splitline.append(f'region_{i}')  # add a unique identifier
                    case 4:
                        splitline.append('.')            # add a score placeholder
                    case 5:
                        splitline.append('.')            # add a strand placeholder

            # if strand was not specified, intersect with both
            if splitline[5] != '+' and splitline[5] != '-':

                splitline[5] = '+'
                line_string = "\t".join(splitline)
                checked_regions.append(line_string)

                splitline[5] = '-'
                line_string = "\t".join(splitline)
                checked_regions.append(line_string)

            else:
                line_string = "\t".join(splitline)
                checked_regions.append(line_string)

    # write checked regions to file
    logger.debug("Writing checked regions to file")
    regions = os.path.join(outdir, "regions.tmp")
    with open(regions, "w") as f:
        for line in checked_regions:
            f.write(line + "\n")

    # Read TFs
    if args.TF is not None:
        logger.debug("Reading TFs from: " + args.TF)
        tfs = pd.read_csv(args.TF, header=None, names=["TF"])
        tfs = tfs["TF"].tolist()

        # get only those TFs that were provided
        tfbs_files = [tfbs for tfbs in tfbs_files if os.path.basename(os.path.dirname(tfbs)) in tfs]

    # intersect each TFBS file with the regions
    logger.debug('Intersecting query with all TBFS files')
    all_intersections = ''
    for file in tfbs_files:
        tf = os.path.basename(os.path.dirname(file))

        # remove header
        headless_file = os.path.join(outdir, f"{tf}.headless")
        command = f'sed "1d" {file} > {headless_file}'
        subprocess.run(command, shell=True)

        # intersect
        command = f"bedtools intersect -a {regions} -b {headless_file} -wa -wb"
        intersection = subprocess.check_output(command, shell=True).decode("utf-8")
        all_intersections += intersection

        # remove headless file
        os.remove(headless_file)

    # make list out of lines
    all_intersections = all_intersections.strip().split("\n")

    # get tfbs header
    with open(tfbs_files[0]) as f:
        tfbs_header = f.readline()

    # merge all intersection files
    logger.debug('Merging all intersections')
    with open(args.output, "w") as f:

        if args.output.endswith(".bed"):
            pass
        else:
            query_header = "query chr\tquery start\tquery end\tquery id\tquery score\tquery strand\t"
            header = query_header + tfbs_header
            f.write(header)

        for line in all_intersections:
            f.write(line + "\n")

    # filter output file
    df = pd.read_csv(args.output, sep="\t")

    logger.debug('Sorting')

    # remove 'chr' str from all chr columns TODO not optimal as contigs may be named differently, i.e. 'contig' or 'chrIV'
    df["query chr"] = df["query chr"].str.replace("chr", "")
    df["TFBS_chr"] = df["TFBS_chr"].str.replace("chr", "")

    # type cast columns
    df["query start"] = df["query start"].astype(int)
    df["query end"] = df["query end"].astype(int)
    df["TFBS_start"] = df["TFBS_start"].astype(int)
    df["TFBS_end"] = df["TFBS_end"].astype(int)
    df['query chr'] = df['query chr'].astype(int)
    df['TFBS_chr'] = df['TFBS_chr'].astype(int)

    df.sort_values(by=["query chr", "query start", "TFBS_name", "TFBS_chr", "TFBS_start"], inplace=True)

    df['query chr'] = 'chr' + df['query chr'].astype(str)
    df['TFBS_chr'] = 'chr' + df['TFBS_chr'].astype(str)

    df.to_csv(args.output, sep="\t", index=False)

    if args.output.endswith(".xlsx"):
        df = pd.read_csv(args.output, sep="\t")
        df.to_excel(args.output, index=False)

    # delete tmp file
    os.remove(regions)

    logger.info("Output written to: " + args.output)
    logger.end()


def main():

    parser = argparse.ArgumentParser()
    parser = add_submerge_arguments(parser)
    args = parser.parse_args()

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit()

    run_submerge(args.tfbs, args.regions, args.TF, args.output, args.verbosity)


if __name__ == '__main__':
    main()
