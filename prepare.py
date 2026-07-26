import os
import shutil
import subprocess
import sys
import Make.common as cm

def reject_sudo():
	"""
	This script elevates only the package installation, so it must not be run under sudo
	wholesale: the submodule checkout below would run as root and leave the working tree
	and .git/modules owned by root, i.e. read-only for the developer afterwards.

	Genuine root sessions (containers, CI images) have no SUDO_UID and are left alone.
	"""
	if os.geteuid() != 0 or "SUDO_UID" not in os.environ:
		return

	print(">> Error: do not run this script with sudo.")
	print(">> It elevates the package installation on its own, and running the whole script")
	print(">> as root makes root own the submodule checkout, leaving it read-only for you.")
	print(">> Run it as yourself instead: python3 prepare.py")
	sys.exit(-1)


reject_sudo()

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
JUCE_PACKAGES = "libasound2-dev libjack-jackd2-dev ladspa-sdk libcurl4-openssl-dev libfontconfig1-dev libx11-dev libxcomposite-dev libxcursor-dev libxext-dev libxinerama-dev libxrandr-dev libxrender-dev libglu1-mesa-dev mesa-common-dev"

# Packages whose name moved between Debian/Ubuntu releases. First available wins; the
# webkit soname is tied to the libsoup generation and freetype dropped its "6" suffix.
JUCE_PACKAGE_ALTERNATIVES = [
	["libwebkit2gtk-4.1-dev", "libwebkit2gtk-4.0-dev"],
	["libfreetype-dev", "libfreetype6-dev"],
]


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


def apt_package_exists(name):
	"""Whether apt knows a package by this name and has a candidate version for it."""
	try:
		policy = subprocess.check_output(["apt-cache", "policy", name], stderr = subprocess.DEVNULL)
	except (OSError, subprocess.CalledProcessError):
		return False

	policy = policy.decode("utf-8", "replace")

	# apt-cache prints nothing at all for a name it doesn't recognise, and reports
	# "Candidate: (none)" for one it knows only as a virtual/unavailable package.
	return bool(policy.strip()) and "Candidate: (none)" not in policy


def resolve_alternatives(alternatives):
	"""Picks the first installable name from each group, for packages renamed across releases."""
	resolved = []

	for candidates in alternatives:
		found = next((c for c in candidates if apt_package_exists(c)), None)

		if found is None:
			# Let apt report it, rather than second-guessing an unfamiliar release here.
			print(">> Warning: none of these are available: " + " ".join(candidates))
			found = candidates[0]

		resolved.append(found)

	return resolved


if cm.is_linux:
	print(">> Checking to see if development dependencies need to be installed ...")

	if cm.is_debian_like:
		print(">> Detected " + cm.os_release.get("PRETTY_NAME", "a Debian-like distribution"))

		apt("update")

		# Resolved after "apt update", so the package lists are populated.
		packages = TOOLCHAIN_PACKAGES.split() + JUCE_PACKAGES.split() + resolve_alternatives(JUCE_PACKAGE_ALTERNATIVES)

		apt("install -y " + " ".join(packages))

		# build-essential only pulls in the distro's default g++, which on releases older
		# than Ubuntu 24.04 / Debian 13 is too old for us. Probed after the install above,
		# since there may not have been a compiler to probe before it.
		if (cm.gcc_major("g++") or 0) < cm.MINIMUM_GCC_MAJOR:
			print(">> Default g++ is older than " + str(cm.MINIMUM_GCC_MAJOR) + ", installing a newer one ...")
			apt("install -y g++-" + str(cm.MINIMUM_GCC_MAJOR))
	else:
		print(">> Warning: " + cm.os_release.get("PRETTY_NAME", "this distribution") + " is not apt-based, set up dependencies yourself")
		print(">> Required toolchain: " + TOOLCHAIN_PACKAGES + " (with g++ " + str(cm.MINIMUM_GCC_MAJOR) + " or newer)")
		print(">> Required libraries: " + JUCE_PACKAGES)
		print(">> Plus, under whatever names apply: " + ", ".join(" or ".join(c) for c in JUCE_PACKAGE_ALTERNATIVES))

print(">> Dev environment setup without errors.")
