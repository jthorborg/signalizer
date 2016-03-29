import subprocess
import os
from time import gmtime, strftime
import getpass
import platform

def rewrite_version_header(where, major, minor, build):
	contents = "#define SIGNALIZER_MAJOR " + major + "\n#define SIGNALIZER_MINOR " + minor + "\n#define SIGNALIZER_BUILD " + build + "\n"
	with open(where, "w") as out:
		out.writelines(contents)

def create_build_file(where, vstring):
	# add latest git commit to build log
	git = subprocess.Popen("git --git-dir ../.git log -2", shell = True, stdout=subprocess.PIPE)
	git_log = git.stdout.read()

	build_info = strftime("%Y-%m-%d %H:%M:%S", gmtime()) + ": Signalizer " + vstring + " built on " + platform.system() + " " + platform.release() + " by " + getpass.getuser() + "\n\n"

	with open(where, "w") as out:
		out.writelines(build_info)
           # error sometimes?
		out.writelines(git_log.decode('ascii'))


join = os.path.join
