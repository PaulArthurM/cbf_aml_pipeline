import sys
import re
import subprocess
import os.path
from os import listdir
from os.path import isfile, join

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
            m = re.search("(SJCBF[0-9]+_[DG]-\w+.[0-9])", string)
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
        cmd = "java -jar /data1/scratch/pamesl/projet_cbf/data/bam/EgaDemoClient.jar -p {email} {password} -dr {label}"
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



objets = []
lines = open_file(sys.argv[1])
for line in lines:
    if "SJCBF" in line:
        sample = Sample(line)
        objets.append(sample)


json = {"samples":{}}
for objet in objets:
    if objet.sample_name not in json["samples"]:
        json["samples"][objet.sample_name] = {"D":[], "G":[]}
    json["samples"][objet.sample_name][objet.sample_type].append(objet.file_prefix)


n = 0
for objet in objets:
    print(objet.bam_file_name)
    if os.path.isfile("/data1/scratch/pamesl/projet_cbf/data/bam/"+objet.bam_file_name):
        print("File already exist.")

    else:
        request_germline_file(objet)
        if not os.path.isfile("/data1/scratch/pamesl/projet_cbf/data/bam/"+objet.bam_file_name+".cip"):
            download_germline_file(objet)
            decrypt_file(objet)

    n += 1
    if n ==4:
        break
