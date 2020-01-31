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

class Sample():
    def __init__(self, string):

        def get_sample_name(string):
            m = re.search("^WXS-(SJCBF[0-9]+)", string)
            if m:
                return m.group(1)


        def get_sample_type(string):
            m = re.search("^WXS-SJCBF[0-9]+_([DG])", string)
            if m:
                return m.group(1)


        def get_bam_file_name(string):
            m = re.search("(SJCBF[0-9]+_[DG]-[A-Za-z0-9]+\.[0-9]\.bam)\.cip", string)
            if m:
                return m.group(1)


        def get_file_prefix(string):
            m = re.search("(SJCBF[0-9]+_[DG]-\w+\.[0-9])", string)
            if m:
                return m.group(1)


        def get_EGAF(string):
            m = re.search("(EGAF.+$)", string)
            if m:
                return m.group(1)


        def get_name_no_machine_id(string):
            m = re.search("(SJCBF[0-9]+_[DG])-[A-Za-z0-9]+(\.[0-9]\.bam)\.cip", string)
            if m:
                return m.group(1) + m.group(2)

        def get_short_file_name(string):
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
    file = open(file_path, 'r')
    return file.readlines()


def rename_file_in_dir(sample):
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
    if sample.sample_type == "G":
        cmd = "java -jar /data1/scratch/pamesl/projet_cbf/data/bam/EgaDemoClient.jar -p {email} {password} -dr {label} -nt 7"
        cmd = cmd.format(email=sys.argv[2], password=sys.argv[3], label="label_"+sample.egaf_id)
        print(cmd)
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
        process.wait()
        rename_file_in_dir(sample)
        print(process.returncode)


def decrypt_file(sample):
    if sample.sample_type == "G":
        cmd = "java -jar /data1/scratch/pamesl/projet_cbf/data/bam/get_EgaDemoClient.jar -p {email} {password} -dc /data1/scratch/pamesl/projet_cbf/data/bam/{bam}.cip -dck abc"
        cmd = cmd.format(email=sys.argv[2], password=sys.argv[3], bam=sample.bam_file_name)
        print(cmd)
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
        process.wait()
        print(process.returncode)


def write_json(dictionary, conf):
    if os.path.isfile(conf):
        print("File a json file already exist.")
    else:
        with open(conf, 'w') as fp:
            json.dump(dictionary, fp)


def check_merge(sample, files, type):
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


def download_file_pyega3(sample):
    saveto="/data1/scratch/pamesl/projet_cbf/data/bam/{bam_file_name}".format(bam_file_name=sample.name_no_machine_id)
    cmd = "pyega3 -c 4 -cf CREDENTIALS_FILE fetch {egaf_id} --saveto {saveto}".format(egaf_id=sample.egaf_id, saveto=saveto)
    print(cmd)
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    process.wait()
    print(process.returncode)


def check_two_file_forms(sample, directory):
    isFile_1 = os.path.isfile(directory + sample.bam_file_name)
    isFile_2 = os.path.isfile(directory + sample.name_no_machine_id)
    isFile_3 = os.path.isfile(directory + sample.short_file_name)
    if (isFile_1 or isFile_2 or isFile_3):
        return True
    else:
        return False


def main(args):
    #objets = [Sample(line) for line in open_file(args.m) if args.e in line]
    objets = []
    lines = open_file(args.m)
    for line in lines:
        if args.e in line:
            sample = Sample(line)
            objets.append(sample)

    if args.j:
        json_file = {"samples":{}}
        for objet in objets:
            if objet.sample_name not in json_file["samples"]:
                json_file["samples"][objet.sample_name] = {"D":[], "G":[]}
                json_file["samples"][objet.sample_name][objet.sample_type].append(objet.file_prefix)
        write_json(json_file, args.c)


    path = args.p
    dry_run = args.b
    wait_time = args.w
    files = [f for f in glob.glob(path + "*merge.bam")]
    limit = 0
    for objet in objets:
        if limit < args.l:
            print("\n\n")
            print(objet.bam_file_name)
            if not check_merge(objet, files, args.t):
                if check_two_file_forms(objet, args.d):
                    print("File already exist.")
                else:
                    time.sleep(wait_time)
                    print("Sample {sample} is being downloaded.".format(sample=objet.file_prefix))
                    if not dry_run:
                        download_file_pyega3(objet)
                    limit+=1
        else: break


if __name__== '__main__':

    parser = argparse.ArgumentParser(description='Downloading files from EGA and create proper JSON file.')
    parser.add_argument('-m', default='/data1/scratch/pamesl/projet_cbf/Sample_File_SJCBF.map', type=str, help="Metadata file.")
    parser.add_argument('-e', default='SJCBF', type=str, help="Experience. Default: SJCBF.")
    parser.add_argument('-j', default=True, type=bool, help="Create JSON. Default: True")
    parser.add_argument('-c', default="data.json", type=str, help="Config JSON name. Default: data.json")
    parser.add_argument('-p', default='/data1/scratch/pamesl/projet_cbf/data/bam/', type=str, help="Path to bam files repository.")
    parser.add_argument('-l', default=1, type=int, help="Number of files to download. Default: 1.")
    parser.add_argument('-t', default=None, type=str, help="Type of files to download. For exemple: 'D' or 'G'.", required=True)
    parser.add_argument('-d', default="~/", type=str, help="Location where to scan existing files.")
    parser.add_argument('-b', default=False, type=bool, help="Dryrun")
    parser.add_argument('-w', default=5, type=int, help="Wait time")

    args = parser.parse_args()

    main(args)
