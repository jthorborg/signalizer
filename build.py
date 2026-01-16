import os
import Make.common as cm

try:

	os.chdir("Make")

	program = cm.ProgramConfig("config.ini")
	print("------> Building Signalizer v. " + program.version_string + " release" if program.release else " debug" + " targets")
	program.rewrite_version_header()

	if cm.is_linux:
		import Make.build_linux as build
		zx = build.build_linux(program)

	print("------> Built Signalizer successfully into:")
	print("------> " + zx)

	# done, if we made it here, increase the conf build
	program.flush()

finally:
	os.chdir("..")
