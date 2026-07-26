import os
import shutil
import sys
import Make.common as cm

print(">> Updating submodules...")

os.system("git submodule sync")
os.system("git submodule update --init --recursive")
os.system("git submodule update")

# The toolchain itself. Without these the JUCE makefile fails in a very unhelpful way:
# it detects the target architecture by scraping an intentional #error out of the compiler's
# output, so a missing compiler is substituted verbatim into a rule and make reports
# "target pattern contains no '%'" instead of anything resembling the real problem.
TOOLCHAIN_PACKAGES = "build-essential make pkg-config git"

# https://github.com/juce-framework/JUCE/blob/develop/docs/Linux%20Dependencies.md
JUCE_PACKAGES = "libasound2-dev libjack-jackd2-dev ladspa-sdk libcurl4-openssl-dev libfreetype-dev libfontconfig1-dev libx11-dev libxcomposite-dev libxcursor-dev libxext-dev libxinerama-dev libxrandr-dev libxrender-dev libwebkit2gtk-4.1-dev libglu1-mesa-dev mesa-common-dev"


def sudo_prefix():
	if os.geteuid() == 0:
		return ""

	if shutil.which("sudo") is None:
		print(">> Error: not running as root and sudo is unavailable, cannot install packages.")
		sys.exit(-1)

	return "sudo "


def apt(command):
	full = sudo_prefix() + "apt-get " + command
	print(">> " + full)

	if os.system(full) != 0:
		print(">> Error: failed running: " + full)
		sys.exit(-1)


if cm.is_linux:
	print(">> Checking to see if development dependencies need to be installed ...")

	if cm.is_ubuntu:
		apt("update")
		apt("install -y " + TOOLCHAIN_PACKAGES + " " + JUCE_PACKAGES)
	else:
		print(">> Warning: Unknown linux distribution, you need to set up all dependencies yourself")
		print(">> Required toolchain: " + TOOLCHAIN_PACKAGES)
		print(">> Required libraries: " + JUCE_PACKAGES)

print(">> Dev environment setup without errors.")
