import os
import Make.common as cm

try:

	os.chdir("Make")

	program = cm.ProgramConfig("config.ini")
	print("------> Building Signalizer v. " + program.version_string + " " + program.configString + " targets")
	program.rewrite_version_header()

	if cm.is_linux:
		import Make.build_linux as build_linux
		zx = build_linux.build(program)
	elif cm.is_windows:
		import Make.build_win as build_win
		zx = build_win.build(program)
	elif cm.is_mac:
		import Make.build_osx as build_osx
		zx = build_osx.build(program)
		
	print("------> Built Signalizer successfully into:")
	print(zx)

	# done, if we made it here, increase the conf build
	program.flush()

finally:
	os.chdir("..")
