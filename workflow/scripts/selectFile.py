import sys
import time
import re
import subprocess
import os.path
from os import listdir
from os.path import isfile, join
import json
import glob
import argparse
import pandas
import click


class File():
    """A class for files manipulation."""
    def __init__(self, string):

        def get_sample_name(string):
            """Return sample's name."""
            m = re.search("^WXS-(SJCBF[0-9]+)", string)
            if m:
                return m.group(1)


        def get_sample_type(string):
            """Return sample's type [D|G] (diagnosis or germline)"""
            m = re.search("^WXS-SJCBF[0-9]+_([DG])", string)
            if m:
                return m.group(1)


        def get_bam_file_name(string):
            """Return bam file's name."""
            m = re.search("(SJCBF[0-9]+_[DG]-[A-Za-z0-9]+\.[0-9]\.bam)\.cip", string)
            if m:
                return m.group(1)


        def get_file_prefix(string):
            """Return sample's prefix."""
            m = re.search("(SJCBF[0-9]+_[DG]-\w+\.[0-9])", string)
            if m:
                return m.group(1)


        def get_EGAF(string):
            """Return sample's EGAF ID."""
            m = re.search("(EGAF.+$)", string)
            if m:
                return m.group(1)


        def get_name_no_machine_id(string):
            """Return sample's filename without machine's ID."""
            m = re.search("(SJCBF[0-9]+_[DG])-[A-Za-z0-9]+(\.[0-9]\.bam)\.cip", string)
            if m:
                return m.group(1) + m.group(2)

        def get_short_file_name(string):
            """Return short sample's name."""
            m = re.search("(SJCBF[0-9]+_[DG])-[A-Za-z0-9]+\.[0-9]\.bam\.cip", string)
            if m:
                return m.group(1) + ".bam"

        self.sample_name = get_sample_name(string)
        self.sample_type = get_sample_type(string)
        self.bam_file_name = get_bam_file_name(string)
        self.file_prefix = get_file_prefix(string)
        self.egaf_id = get_EGAF(string)
        self.name_no_machine_id = get_name_no_machine_id(string)
        self.short_file_name = get_short_file_name(string)


def open_file(file_path):
    """Return lines from file."""
    file = open(file_path, 'r')
    return file.readlines()


def rename_file_in_dir(sample):
    """Return files in directory. Deprecated."""
    path = "/data1/scratch/pamesl/projet_cbf/data/bam/"
    files_in_dir = [f for f in listdir(path) if isfile(join(path, f))]
    for file in files_in_dir:
        regex = "(^EGAR[0-9]+_" + sample.bam_file_name + ")"
        m = re.search(regex, file)
        if m:
            downloaded_file = m.group(1)
            new_file = sample.bam_file_name
            cmd = "mv {downloaded_file}.cip {new_file}.cip"
            cmd = cmd.format(downloaded_file=downloaded_file, new_file=new_file)
            print(cmd)
            process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
            process.wait()
            print(process.returncode)


def request_germline_file(sample):
    """Download germline bam files associated to given sample. Deprecated."""
    if sample.sample_type == "G":
        cmd = "java -jar /data1/scratch/pamesl/projet_cbf/data/bam/EgaDemoClient.jar -p {email} {password} -rf {egaf_id} -re abc -label {label}"
        cmd = cmd.format(email=sys.argv[2], password=sys.argv[3], egaf_id=sample.egaf_id, label="label_"+sample.egaf_id)
        print(cmd)
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
        process.wait()
        print(process.returncode)
    else:
        print("Not a germline file.")


def download_germline_file(sample):
    """Download diagnosis bam files associated to given sample. Deprecated."""
    if sample.sample_type == "G":
        cmd = "java -jar /data1/scratch/pamesl/projet_cbf/data/bam/EgaDemoClient.jar -p {email} {password} -dr {label} -nt 7"
        cmd = cmd.format(email=sys.argv[2], password=sys.argv[3], label="label_"+sample.egaf_id)
        print(cmd)
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
        process.wait()
        rename_file_in_dir(sample)
        print(process.returncode)


def decrypt_file(sample):
    """Decrypt bam files. Deprecated."""
    if sample.sample_type == "G":
        cmd = "java -jar /data1/scratch/pamesl/projet_cbf/data/bam/get_EgaDemoClient.jar -p {email} {password} -dc /data1/scratch/pamesl/projet_cbf/data/bam/{bam}.cip -dck abc"
        cmd = cmd.format(email=sys.argv[2], password=sys.argv[3], bam=sample.bam_file_name)
        print(cmd)
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
        process.wait()
        print(process.returncode)


#TODO: update to new sample sheet format. Add path parameter.
def write_json(dictionary, conf):
    """Write a JSON file (sample sheet)."""
    if os.path.isfile(conf):
        print("File a json file already exist.")
    else:
        with open(conf, 'w') as fp:
            json.dump(dictionary, fp)


#TODO: update for new pipeline version.
def check_merge(sample, files, type):
    """Check if a merged bam files already exist."""
    sample_name = sample.sample_name
    if sample.sample_type == type:
        for f in files:
            m  = re.search(sample_name, f)
            m2 = re.search(type, f)
            if m and m2:
                print("A merged file already exist for this sample.")
                return True
    else:
        print("Not a requested file.")
        return True


#TODO: add path parameters.
def download_file_pyega3(sample, path_to_bam):
    """Download bam files using pyega3."""
    saveto="{path_to_bam}{bam_file_name}".format(path_to_bam=path_to_bam, bam_file_name=sample.name_no_machine_id)
    cmd = "pyega3 -c 4 -cf CREDENTIALS_FILE fetch {egaf_id} --saveto {saveto}".format(egaf_id=sample.egaf_id, saveto=saveto)
    print(cmd)
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    process.wait()
    print(process.returncode)


def check_two_file_forms(sample, directory):
    """Check if a sample exist as file."""
    isFile_1 = os.path.isfile(directory + sample.bam_file_name)
    isFile_2 = os.path.isfile(directory + sample.name_no_machine_id)
    isFile_3 = os.path.isfile(directory + sample.short_file_name)
    if (isFile_1 or isFile_2 or isFile_3):
        print("File already exist.")
        return True
    else:
        print("No file exist.")
        return False


def extract_metadata(metadata, experiment):
    metadata = [File(line) for line in open_file(metadata) if experiment in line]
    return metadata


def create_sample_sheet(samples_json_file, sample_sheet_name):
    index = samples_json_file.keys()
    sample_sheet = pandas.DataFrame(columns=['samples', 'germline_path', 'somatic_path'], index=index)
    for sample in samples_json_file:
        sample_sheet.loc[sample, 'samples'] = sample
        sample_sheet.loc[sample, 'germline_path'] = " ".join(samples_json_file[sample]['germline_path'])
        sample_sheet.loc[sample, 'somatic_path'] = " ".join(samples_json_file[sample]['somatic_path'])
    save_path = "config/" + sample_sheet_name
    sample_sheet.to_csv(save_path, sep = ";", index = False)



def create_json(metadata):
    """Return a JSON file with a resume of sample. Necessary for sample-sheet creation."""
    samples_json_file = {}
    for file in metadata:
        if not file.sample_name in samples_json_file:
            samples_json_file[file.sample_name] = {"germline_path": [], "somatic_path": []}
        if file.sample_type == "G":
            samples_json_file[file.sample_name]["germline_path"].append(file.name_no_machine_id)
        elif file.sample_type == "D":
            samples_json_file[file.sample_name]["somatic_path"].append(file.name_no_machine_id)
    return samples_json_file



@click.command()
@click.option('-e', '--experience', default=None, help="Experience.", required=True)
@click.option('-s', '--sample-sheet', default=True, type=bool, help="Create sample-sheet. Default: True")
@click.option('-n', '--sample-sheet-name', default='sample-sheet.csv', help="Create sample-sheet. Default: 'sample-sheet.csv'")
@click.option('-b', '--path-to-bam', default=None, help="Path to bam files repository.", required=True)
@click.option('-h', '--number', default=1, type=int, help="Number of files to download. Default: 1.")
@click.option('-t', '--type', default=None, help="Type of files to download. For exemple: 'D' or 'G'.", required=True)
@click.option('-l', '--location', default=None, help="Location where to scan existing files.", required=True)
@click.option('-d', '--dry-run', default=False, type=bool, help="Dryrun. Default: False")
@click.option('-w', '--latency', default=5, help="latency between downloads. Default: 5 seconds.")
@click.argument('sample_map')
def main(sample_map, experience, sample_sheet, sample_sheet_name, path_to_bam, number, type, location, dry_run, latency):
    """Wrapper for Pyega3 utility. Download files from EGA and create proper sample-sheet."""

    metadata = extract_metadata(sample_map, experience)

    samples_json_file = create_json(metadata)

    if sample_sheet:
        create_sample_sheet(samples_json_file, sample_sheet_name)

    files = [f for f in glob.glob(path_to_bam + "*.bam")]
    limit = 0
    for sample in metadata:
        if limit < number:
            print("\n\n")
            print(sample.bam_file_name)
            if not check_merge(sample, files, type):
                if not check_two_file_forms(sample, location):
                    if dry_run:
                        print("Sample {sample} is being downloaded (dry-run).".format(sample=sample.file_prefix))
                    else:
                        print("Sample {sample} is being downloaded.".format(sample=sample.file_prefix))
                        download_file_pyega3(sample, path_to_bam)
                        time.sleep(latency)
                    limit+=1
        else: break


if __name__== '__main__':
    main()
