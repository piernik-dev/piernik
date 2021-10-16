#!/usr/bin/python3

import json
import requests
import time

base_job = "https://jenkins.camk.edu.pl"
job_name = "user_gold%20CI"

try:
    url = base_job + "/job/" + job_name + "/lastBuild/api/json"
    finished = False
    while not finished:
        data = requests.get(url).json()
        if data['building']:
            time.sleep(10)
            print("Waiting ...")
        else:
            finished = True

    gold_url = base_job + "/job/" + job_name + "/lastBuild/artifact/jenkins/goldexec/gold.txt"
    gold_out = requests.get(gold_url)

    print(gold_out.text)

    if data['result'] == "SUCCESS":
        print("Job is success")
        exit(0)
    else:
        print("Job status failed")
        exit(-1)

except Exception as e:
    print(str(e))
    exit(-2)
