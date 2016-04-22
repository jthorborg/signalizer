import subprocess
import os
from time import gmtime, strftime
import getpass
import platform

def rewrite_version_header(where, major, minor, build):
	build_info = get_custom_build_info().replace('\n', "\\n").replace('\r', "\\n")
	contents = "#define SIGNALIZER_MAJOR " + major + "\n#define SIGNALIZER_MINOR " + minor + "\n#define SIGNALIZER_BUILD " + build 
	contents += "\n#define SIGNALIZER_BUILD_INFO \"" + build_info + "\"\n"
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

join = os.path.join
