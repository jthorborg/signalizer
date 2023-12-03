import os

print(">> Updating submodules...")

os.system("git submodule update --init --recursive")
os.system("git submodule update")

print(">> Dev environment setup without errors.")