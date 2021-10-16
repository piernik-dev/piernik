#!/usr/bin/python3

import json
import requests
import time
import pprint

base_job = "https://jenkins.camk.edu.pl"
gold = "user_gold%20CI"

jobs = []
for i in (gold, "3-body%20CI", "64-bit%20int%20CR%20test%20CI", "dep.png%20CI",
          "IO%20v2%20CI", "Jeans%20CI", "Maclaurin%20CI", "mergeable%20benchmarking%20CI",
          "noHDF5%20CI", "PEP8%20CI", "Piernik%20CI", "python3%20CI", "QA%20CI"):
    jobs.append({"name": i, "finished": False})

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

            for ia in data["actions"]:
                for ja in ia:
                    if ja == "parameters":
                        for k in ia[ja]:
                            if k["name"] == "ghprbActualCommit":
                                print("Testing commit " + k["value"] + " in '" + j["name"].replace("%20", " ") + "'")

            if j["finished"]:
                if data['result'] == "SUCCESS":
                    print("Job is success")
                    # exit(0)
                else:
                    print("Job status failed")
                    if j["name"] == gold:
                        gold_url = base_job + "/job/" + j["name"] + "/lastBuild/artifact/jenkins/goldexec/gold.txt"
                        gold_out = requests.get(gold_url)
                        print(gold_out.text)

                    # exit(-1)

    if not all_finished:
        time.sleep(10)
