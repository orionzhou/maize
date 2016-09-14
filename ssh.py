#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import time
import paramiko

username = 'zhoup'
hostname = 'mesabi.msi.umn.edu'
port = 22
pubkey = os.path.join(os.environ['HOME'], '.ssh', 'id_rsa')

if __name__ == '__main__':
    key = paramiko.RSAKey.from_private_key_file(pubkey)
    s = paramiko.SSHClient()
    s.load_system_host_keys()
    s.connect(hostname, port, pkey=key)
    stdin, stdout, stderr = s.exec_command("echo $HOSTNAME")
    for line in stdout:
        print "... " + line.strip("\n")
    s.close()
