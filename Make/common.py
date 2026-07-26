import subprocess
import os
from time import gmtime, strftime
import getpass
import platform
import configparser
import sys
import shutil as sh

join = os.path.join

# cpl's profiling model nests an unordered_map of a type inside that same type. Incomplete
# value types are only guaranteed for vector/list/forward_list, so this leans on libstdc++
# behaviour: GCC 12 deferred the node instantiation that GCC 11 performs eagerly. GCC 11
# therefore fails with "pair::second has incomplete type". Distros older than Ubuntu 24.04
# still default to GCC 11, but package a newer g++ alongside it.
MINIMUM_GCC_MAJOR = 12

sys_name = platform.system().lower()
is_windows = "windows" in sys_name
is_mac = "darwin" in sys_name
is_linux = "linux" in sys_name


def read_os_release(path = "/etc/os-release"):
	"""
	Parses os-release into a dict. Empty on non-Linux, or on the rare Linux without it.
	See https://www.freedesktop.org/software/systemd/man/latest/os-release.html
	"""
	values = {}

	try:
		with open(path) as f:
			for line in f:
				line = line.strip()

				if not line or line.startswith("#") or "=" not in line:
					continue

				key, _, value = line.partition("=")
				values[key.strip()] = value.strip().strip('"').strip("'")
	except OSError:
		pass

	return values


os_release = read_os_release() if is_linux else {}

# ID is the distribution, ID_LIKE the families it derives from, so Mint reports
# "linuxmint"/"ubuntu" and Pop!_OS reports "pop"/"ubuntu debian". Testing membership of
# both covers derivatives without enumerating them.
distro_families = set(filter(None, [os_release.get("ID", "")] + os_release.get("ID_LIKE", "").split()))

is_ubuntu = "ubuntu" in distro_families
# Everything managing packages with apt/dpkg: Debian, Ubuntu, Mint, Pop!_OS, MX, AV Linux,
# KXStudio, Raspberry Pi OS. Preferred over is_ubuntu for anything install-related.
is_debian_like = is_linux and not distro_families.isdisjoint({"debian", "ubuntu"})


def gcc_major(compiler):
	"""Major version of `compiler`, or None if it isn't installed or isn't usable."""
	try:
		version = subprocess.check_output([compiler, "-dumpfullversion"], stderr = subprocess.DEVNULL)
	except (OSError, subprocess.CalledProcessError):
		return None

	try:
		return int(version.decode("ascii").strip().split(".")[0])
	except ValueError:
		return None


def rewrite_version_header(where, major, minor, build):
	build_info = get_custom_build_info().replace('\n', "\\n").replace('\r', "\\n")
	contents = "#define SIGNALIZER_MAJOR " + major + "\n#define SIGNALIZER_MINOR " + minor + "\n#define SIGNALIZER_BUILD " + build
	contents += "\n#define SIGNALIZER_BUILD_INFO \"" + build_info + "\"\n"
	contents += "\n#define SIGNALIZER_VERSION " + major + "." + minor + "." + build
	contents += "\n#define SIGNALIZER_VERSION_STRING \"" + major + "." + minor + "." + build + "\""
	contents += "\n#define SIGNALIZER_VST_VERSION_HEX " + "0x{0:02x}{1:02x}{2:02x}".format(int(major) % 0xff, int(minor) % 0xff, int(build) % 0xff)
	with open(where, "w") as out:
		out.writelines(contents)

def create_build_file(where, vstring):
	# add latest git commit to build log
	git = subprocess.Popen("git --git-dir ../.git log -5", shell = True, stdout=subprocess.PIPE)
	git_log = git.stdout.read()

	build_info = strftime("%Y-%m-%d %H:%M:%S", gmtime()) + ": Signalizer " + vstring + " built on " + platform.system() + " " + platform.release() + " by " + getpass.getuser() + "\n"
	build_info += get_custom_build_info() + "\n\n"

	with open(where, "w") as out:
		out.writelines(build_info)
		out.writelines(git_log.decode('ascii'))

def get_custom_build_info():
	git = subprocess.Popen("git --git-dir ../.git branch -q", shell = True, stdout=subprocess.PIPE)
	git_branch = git.stdout.read()
	git = subprocess.Popen("git --git-dir ../.git describe --always", shell = True, stdout=subprocess.PIPE)
	git_description = git.stdout.read()
	return git_branch.decode('ascii') + "\n" + git_description.decode('ascii')


class DevConfig:
	def __init__(self, args):
		if os.path.split(os.getcwd())[1] != "Make":
			print("Error, this program must be called from within Make folder")
			exit(-5)

		self.arch = args.arch
		self.verbose = args.verbose
		self.max_errors = args.errors
		self.logs_dir = "Logs"
		self.optimized = args.optimized
		self.configString = "Release" if self.optimized else "Debug"


class ProgramConfig:
	def __init__(self, args, inifile):

		if os.path.split(os.getcwd())[1] != "Make":
			print("Error, this program must be called from within Make folder")
			exit(-5)

		self.flush_parameters = False
		self.inifile = inifile

		self.release = not args.debug
		self.configString = "Release" if self.release else "Debug"
		self.skipvst2 = args.skipvst2
		self.verbose = args.verbose

		self.config = configparser.ConfigParser()
		self.config.read(inifile)

		for param in [["major", args.increase_major], ["minor", args.increase_minor], ["build", args.increase_build]]:
			if param[1]:
				self.flush_parameters = True
				self.config.set("version", param[0], str(int(self.config.get("version", param[0])) + 1))
				print("------> Increasing " + param[0] + " to " + self.config.get("version", param[0]))

		self.major = self.config.get("version", "major")
		self.minor = self.config.get("version", "minor")
		self.build = self.config.get("version", "build")
		self.company = self.config.get("info", "company")
		self.author = self.config.get("info", "author")
		self.desc = self.config.get("info", "description")
		self.name = self.config.get("info", "productname")
		self.manu4 = self.config.get("info", "manu4")
		self.sub4 = self.config.get("info", "sub4")

		self.version_string = self.major + "." + self.minor + "." + self.build

	def flush(self):
		if self.flush_parameters:
			with open(self.inifile, "w") as f:
				self.config.write(f, True)

	def rewrite_version_header(self):
		rewrite_version_header("../Source/version.h", self.major, self.minor, self.build)

	def make_release_folder_with_goodies(self, build_dir, advice_file):

		if os.path.exists(build_dir):
			sh.rmtree(build_dir)

		os.makedirs(build_dir)
		create_build_file(join(build_dir, "Build.log"), self.version_string)
		sh.copyfile(join("Skeleton", "READ ME.txt"), join(build_dir, "READ ME.txt"))
		sh.copyfile("../CHANGELOG.md", join(build_dir, "CHANGELOG.md"))
		sh.copyfile(advice_file, join(build_dir, "HOW TO INSTALL.txt"))
