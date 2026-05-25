import subprocess
import os
from time import gmtime, strftime
import getpass
import platform
import configparser
import sys
import argparse
import shutil as sh
import platform

join = os.path.join

sys_name = platform.system().lower()
is_windows = "windows" in sys_name
is_mac = "darwin" in sys_name
is_linux = "linux" in sys_name
is_ubuntu = is_linux and "ubuntu" in platform.version().lower()

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
           # error sometimes?
		out.writelines(git_log.decode('ascii'))

def get_custom_build_info():
    git = subprocess.Popen("git --git-dir ../.git branch -q", shell = True, stdout=subprocess.PIPE)
    git_branch = git.stdout.read()
    git = subprocess.Popen("git --git-dir ../.git describe --always", shell = True, stdout=subprocess.PIPE)
    git_description = git.stdout.read()
    return git_branch.decode('ascii') + "\n" + git_description.decode('ascii')

class ProgramConfig:
	def __init__(self, inifile):

		if os.path.split(os.getcwd())[1] != "Make":
			print("Error, this program must be called from within Make folder")
			exit(-5) 

		self.flush_parameters = False
		self.inifile = inifile
		parameters = []

		parser = argparse.ArgumentParser()
		parser.add_argument("-d", "--debug", action="store_true", help="Make a release with debug code generation")
		parser.add_argument("-s", "--skipvst2", action="store_true", help="Try to make a release without depending on vst2 SDK (otherwise, place this at ../SDKs/vstsdk2.4)")

		parser.add_argument("-j", "--increase-major", action="store_true", help="Increase the major version by 1")
		parser.add_argument("-n", "--increase-minor", action="store_true", help="Increase the minor version by 1")
		parser.add_argument("-b", "--increase-build", action="store_true", help="Increase the build version by 1")

		args = parser.parse_args()

		self.release = not args.debug
		self.configString = "Release" if self.release else "Debug"
		self.skipvst2 = args.skipvst2

		self.config = configparser.ConfigParser()
		self.config.read(inifile)

		# handle operations
		for param in [["major", args.increase_major], ["minor", args.increase_minor], ["build", args.increase_build]]:
			if param[1]:
				self.flush_parameters = True
				self.config.set("version", param[0], str(int(self.config.get("version", param[0])) + 1))
				print("------> Increasing " + param[0] + " to " + self.config.get("version", param[0]))

		#configurations
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
			with open(inifile, "w") as f:
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
