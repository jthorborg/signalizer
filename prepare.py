import os
import configparser
from urllib.request import urlopen
from shutil import copyfile
import zipfile
import sys
import platform

print(">> Updating submodules...")
os.system("git submodule update --init --recursive")
os.system("git submodule update")

print(">> Dev environment setup without errors.")
