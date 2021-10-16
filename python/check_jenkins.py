#!/usr/bin/python3

import json
import requests
import time
import pprint

c_white = '\033[1m'
c_red = '\033[91m'
c_green = '\033[92m'
c_reset = '\033[0m'

fmt = "%-*s %-40s %s"
base_job = "https://jenkins.camk.edu.pl"
gold = "user_gold%20CI"

jobs = []
longestkey = 0
for i in ("3-body%20CI", "64-bit%20int%20CR%20test%20CI", "dep.png%20CI", "IO%20v2%20CI",
          "Jeans%20CI", "Maclaurin%20CI", "mergeable%20benchmarking%20CI", "noHDF5%20CI",
          "PEP8%20CI", "Piernik%20CI", "python3%20CI", "QA%20CI", gold):
    jobs.append({"name": i, "finished": False})
    longestkey = max(longestkey, len(i.replace("%20", " ")))

print((c_white + fmt + c_reset) % (longestkey, "Test", "SHA1", "status"))

all_finished = False
while not all_finished:
    all_finished = True
    for j in jobs:
        if not j["finished"]:
            try:
                url = base_job + "/job/" + j["name"] + "/lastBuild/api/json"
                data = requests.get(url).json()
                # pprint.pp(data)
                if data['building']:
                    print("Waiting ...")
                    j["finished"] = False
                else:
                    j["finished"] = True

            except Exception as e:
                print(str(e))
                exit(-2)

                all_finished = all_finished and j["finished"]

            j["sha1"] = ""
            for ia in data["actions"]:
                for ja in ia:
                    if ja == "parameters":
                        for k in ia[ja]:
                            if k["name"] == "ghprbActualCommit":
                                j["sha1"] = k["value"]

            # print("Testing commit " + j["sha1"] + " in '" + j["name"].replace("%20", " ") + "'")

            if j["finished"]:
                j["status"] = (data['result'] == "SUCCESS")
                print(fmt % (longestkey, j["name"].replace("%20", " "), j["sha1"], (c_green if j["status"] else c_red) + data["result"] + c_reset))
                if not j["status"]:
                    if j["name"] == gold:
                        gold_url = base_job + "/job/" + j["name"] + "/lastBuild/artifact/jenkins/goldexec/gold.txt"
                        gold_out = requests.get(gold_url)
                        print(gold_out.text)

    if not all_finished:
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
