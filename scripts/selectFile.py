import sys
import re
import subprocess

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
            m = re.search("(SJCBF[0-9]+_[DG]-[A-Za-z0-9]+\.[0-9].bam.cip)", string)
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


def request_germline_file(sample):
    if sample.sample_type == "G":
        cmd = "java -jar EgaDemoClient.jar -p {email} {password} -rf {egaf_id} -re abc -label {label}"
        cmd = cmd.format(email=sys.argv[2], password=sys.argv[3], egaf_id=sample.egaf_id, label="label_"+sample.egaf_id)
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
        process.wait()
        print(process.returncode)


def download_germline_file(sample):
    if sample.sample_type == "G":
        cmd = "java -jar EgaDemoClient.jar -p {email} {password} -dr {label} -path /data1/scratch/pamesl/projet_cbf/data/bam"
        cmd = cmd.format(email=sys.argv[2], password=sys.argv[3], label="label_"+sample.egaf_id)
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
        process.wait()
        print(process.returncode)


def decrypt_file(sample):
    if sample.sample_type == "G":
        cmd = "java -jar EgaDemoClient.jar -p {email} {password} -dc /data1/scratch/pamesl/projet_cbf/data/bam/{bam} -dck abc"
        cmd = cmd.format(email=sys.argv[2], password=sys.argv[3], bam=sample.bam_file_name)
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
    request_germline_file(objet)
    download_germline_file(objet)
    #decrypt_file(objet)
    n += 1
    if n == 3:
        break
