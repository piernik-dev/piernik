#!/usr/bin/python3

import json
import requests
import time
import pprint


def find_sha1(data):
    sha1 = ""
    pr = ""
    src_br = ""
    for ia in data["actions"]:
        for ja in ia:
            if ja == "parameters":
                for k in ia[ja]:
                    if k["name"] == "ghprbActualCommit":
                        sha1 = k["value"]
                    elif k["name"] == "ghprbPullId":
                        pr = k["value"]
                    elif k["name"] == "ghprbSourceBranch":
                        src_br = k["value"]
    return sha1, pr, src_br


c_white = '\033[1m'
c_red = '\033[91m'
c_green = '\033[92m'
c_yellow = '\033[93m'
c_reset = '\033[0m'

fmt = "%-*s %40s %4s %30s %s"
base_job = "https://jenkins.camk.edu.pl"
gold = "user_gold%20CI"

jobs = []
longestkey = 0
for i in ("3-body%20CI", "64-bit%20int%20CR%20test%20CI", "dep.png%20CI", "IO%20v2%20CI",
          "Jeans%20CI", "Maclaurin%20CI", "mergeable%20benchmarking%20CI", "noHDF5%20CI",
          "PEP8%20CI", "Piernik%20CI", "python3%20CI", "QA%20CI", gold):
    jobs.append({"name": i, "finished": False})
    longestkey = max(longestkey, len(i.replace("%20", " ")))

print((c_white + fmt + c_reset) % (longestkey, "Test", "SHA1", "PR#", "branch", "status"))

all_finished = False
while not all_finished:
    all_finished = True
    anything_new = False
    for j in jobs:
        url = base_job + "/job/" + j["name"] + "/lastBuild/api/json"
        try:
            data = requests.get(url).json()
            # pprint.pp(data)
        except Exception as e:
            print(str(e))
            exit(-1)

        if j["finished"]:
            if find_sha1(data)[0] != j["sha1"]:
                print(c_yellow + "SHA1 has changed in the meanwhile for '%s'." % j["name"].replace("%20", " ") + c_reset)
                j["finished"] = False

        if not j["finished"]:
            if data['building']:
                j["finished"] = False
            else:
                j["finished"] = True

            j["sha1"], j["pr"], j["br"] = find_sha1(data)

            if j["finished"]:
                anything_new = True
                j["status"] = (data['result'] == "SUCCESS")
                print(fmt % (longestkey, j["name"].replace("%20", " "), j["sha1"], j["pr"], j["br"], (c_green if j["status"] else c_red) + data["result"] + c_reset))
                if not j["status"]:
                    if j["name"] == gold:
                        gold_url = base_job + "/job/" + j["name"] + "/lastBuild/artifact/jenkins/goldexec/gold.txt"
                        gold_out = requests.get(gold_url)
                        print(gold_out.text)

        all_finished = all_finished and j["finished"]

    if not all_finished:
        if anything_new:
            n_finished = 0
            for j in jobs:
                if j["finished"]:
                    n_finished += 1
            print("Waiting (%d/%d) ..." % (n_finished, len(jobs)))
        time.sleep(10)

exit_stat = 0
same_sha1 = True
for j in jobs:
    if (j["sha1"] != jobs[0]["sha1"]):
        same_sha1 = False
    if not j["status"]:
        exit_stat += 1

if not same_sha1:
    print(c_red + "Some jobs have different SHA1. Please verify the status on GitHub as concurrent tests may be interfering" + c_reset)

if exit_stat == 0:
    print(c_green + "All tests passed" + c_reset)
else:
    print(c_red + "Failed %s out of %d tests." % (exit_stat, len(jobs)) + c_reset)

exit(exit_stat)
