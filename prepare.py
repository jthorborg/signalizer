import os
import platform

print(">> Updating submodules...")

os.system("git submodule sync")
os.system("git submodule update --init --recursive")
os.system("git submodule update")

sys_name = platform.system().lower()

is_windows = "windows" in sys_name
is_mac = "darwin" in sys_name
is_linux = "linux" in sys_name

if is_linux or not (is_windows or is_mac):
	print(">> Checking to see if development dependencies need to be installed ...")

	if "ubuntu" in platform.version().lower():
		# https://github.com/juce-framework/JUCE/blob/develop/docs/Linux%20Dependencies.md
		os.system("apt update")
		os.system("apt install libasound2-dev libjack-jackd2-dev ladspa-sdk libcurl4-openssl-dev libfreetype-dev libfontconfig1-dev libx11-dev libxcomposite-dev libxcursor-dev libxext-dev libxinerama-dev libxrandr-dev libxrender-dev libwebkit2gtk-4.1-dev libglu1-mesa-dev mesa-common-dev")
	else:
		print(">> Warning: Unknown linux distribution, you need to set up all dependencies yourself")

print(">> Dev environment setup without errors.")