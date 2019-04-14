#!/usr/bin/python
# -*- coding: utf-8 -*-
from sys import exit, stdout
import inspect

'''
   Short module for colored output, can be used in all python scripts.
   Returns: name of calling unit + message kind (info, warn, exit) + message content
'''

def talkgreen(message):
   message = "\033[92m"+message+"\033[0m"
   return message
def talkred(message):
   message = "\033[91m"+message+"\033[0m"
   return message
def talkyellow(message):
   message = "\033[93m"+message+"\033[0m"
   return message

'''
   Name of calling_unit is trimmed to len = 9 (default)
   This way it fits yt io length and the name is long enough to know
   which unit called IO functions
'''

calling_unit_namelen = 8

def prtinfo(message):
   calling_unit = ((inspect.stack()[1][1].split("/")[-1]).split(".")[0])[0:calling_unit_namelen].ljust(calling_unit_namelen)
   print (talkgreen(calling_unit+": [INFO] "+message))
   return
def prtwarn(message):
   calling_unit = ((inspect.stack()[1][1].split("/")[-1]).split(".")[0])[0:calling_unit_namelen].ljust(calling_unit_namelen)
   print (talkyellow(calling_unit+": [WARN] "+message))
   return
def die(message):
   calling_unit = ((inspect.stack()[1][1].split("/")[-1]).split(".")[0])[0:calling_unit_namelen].ljust(calling_unit_namelen)
   exit(talkred(calling_unit+": [EXIT] "+message))

def read_var(message):
    stdout.write (talkgreen("[IO]   "+  message+" "))
    var = raw_input()
    return var
