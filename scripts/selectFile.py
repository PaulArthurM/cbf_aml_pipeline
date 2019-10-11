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


        self.sample_name = get_sample_name(string)
        self.sample_type = get_sample_type(string)
        self.bam_file_name = get_bam_file_name(string)
        self.file_prefix = get_file_prefix(string)
        self.egaf_id = get_EGAF(string)


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
        cmd = "java -jar /data1/scratch/pamesl/projet_cbf/data/bam/EgaDemoClient.jar -p {email} {password} -dc /data1/scratch/pamesl/projet_cbf/data/bam/{bam}.cip -dck abc"
        cmd = cmd.format(email=sys.argv[2], password=sys.argv[3], bam=sample.bam_file_name)
        print(cmd)
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
        process.wait()
        print(process.returncode)


def write_json(dictionary):
    if os.path.isfile("/data1/scratch/pamesl/projet_cbf/data/bam/"+objet.bam_file_name):
        print("File a json file already exist.")
    else:
        with open('/data2/home/pamesl/cbf_aml_pipeline/samples.json', 'w') as fp:
            json.dump(dictionary, fp)


def check_merge(sample, files, type):
    sample_name = sample.sample_name
    if sample.sample_type == type:
        for f in files:
            m = re.search(sample_name, f)
            if m:
                print("A merged file already exist for this sample.")
                return True
    else:
        print("Not a requested file.")
        return True


def download_file_pyega3(sample):
    saveto="/data1/scratch/pamesl/projet_cbf/data/bam/{bam_file_name}".format(bam_file_name=sample.bam_file_name)
    cmd = "pyega3 -c 4 -cf CREDENTIALS_FILE fetch {egaf_id} --saveto {saveto}".format(egaf_id=sample.egaf_id, saveto=saveto)
    print(cmd)
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    process.wait()
    print(process.returncode)


if __name__== '__main__':

    parser = argparse.ArgumentParser(description='Downloading files from EGA and create proper JSON file.')
    parser.add_argument('-m', default='/data1/scratch/pamesl/projet_cbf/Sample_File_SJCBF.map', type=str, help="Metadata file")
    parser.add_argument('-e', default='SJCBF', type=str, help="Experience")
    parser.add_argument('-j', default=True, type=bool, help="Create JSON")
    parser.add_argument('-p', default='/data1/scratch/pamesl/projet_cbf/data/bam/', type=str, help="Path to bams")
    parser.add_argument('-l', default=1, type=int, help="Number of files to download.")
    parser.add_argument('-t', default=None, type=str, help="Type of files to download.", required=True)

    args = parser.parse_args()


    objets = []
    lines = open_file(args.m)#open_file(sys.argv[1])
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



    path = args.p
    files = [f for f in glob.glob(path + "*merge.bam", recursive=False)]
    limit = 0
    for objet in objets:
        if limit < args.l:
            print("\n\n")
            print(objet.bam_file_name)
            if not check_merge(objet, files, args.t):
                if os.path.isfile("/data1/scratch/pamesl/projet_cbf/data/bam/"+objet.bam_file_name):
                    print("File already exist.")
                else:
                    time.sleep(5)
                    print("Sample {sample} is being downloaded.".format(sample=objet.file_prefix))
                    #download_file_pyega3(objet)
                    limit+=1
