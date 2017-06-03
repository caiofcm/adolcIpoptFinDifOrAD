# from numpy import *
import re
import os
import optparse

# def run():
# 	with open('powerDouble.h') as f:
# 		for line in f:
# 			if '// CONVADOL' in line:
# 				print('this line: ' + line)

def readData(pathfile):
	dir = os.path.dirname(__file__)
	with open (pathfile, "r") as myfile:
		data = myfile.read()
		return data
	return 'ERROR'


def reSearch(fileName):
# 	src = """#include <adolc/adolc.h> 

# double power(double x /*a*/, double a/*n*/, int n);

# double auxvar(double tspan_12/*n*/, double a/*n*/, int n);"""
	# print(src)

	# INPUTS:
	# fileName = "powerDouble.h"

	# Read from file:
	dir = os.path.dirname(__file__)
	pathfile = os.path.join(dir, fileName)
	data = readData(pathfile)

	# Check if adol is included already:
	includeStr = ""
	rgex = r"#include\s+<adolc/adolc.h>"
	isInclRe = re.compile(rgex)
	allDeactive = isInclRe.search(data)
	if (allDeactive is None):
		includeStr = "#include <adolc/adolc.h>" 
		data = includeStr + data 

	# Not Activated Pattern
	rgex = r"(double)\s+(\w+)\s*/\*n\*/"
	deActRe = re.compile(rgex)
	allDeactive = deActRe.findall(data)
	print(allDeactive)

	# Activated Pattern
	rgex = r"(double)\s+(\w+)\s*/\*a\*/"
	actRe = re.compile(rgex)
	allActive = actRe.findall(data)
	print(allActive)

	# Substitute activated Pattern:
	retAfterSub = actRe.sub(r"a\1\2", data)
	print(retAfterSub)

	# Save to a new file
	pathOfile = os.path.join(dir, "ad_"+ fileName)
	with open(pathOfile, "w") as text_file:
		print(retAfterSub, file=text_file)
	return

	# regex.sub(r"\1")

def getFileNameNoExt(filename):
    rg = r"^(([A-Z]:)?[\.]?[\\{1,2}/]?.*[\\{1,2}/])*(.+)\.(.+)"
    rsltSerch = re.search(rg, filename)
    if rsltSerch is not None:
        filePath = rsltSerch.group(1)
        filenoext = rsltSerch.group(3)
    else:
        filenoext = os.path.basename(filename)
    return filenoext, filePath

def convertActiveDouble(filepath, isHeader = False, isSource=False):
	
	data = readData(filepath)

	# Check if adol is included already:
	if (isHeader):
		rgex = r"#include\s+<adolc/adolc.h>"
		isInclRe = re.compile(rgex)
		allDeactive = isInclRe.search(data)
		if (allDeactive is None):
			includeStr = "#include <adolc/adolc.h>" 
			data = includeStr + data 

	if (isSource):
		filehead = getFileNameNoExt(filepath) + ".h"
		rgex = r'#include "' + filehead + '"'
		isInclRe = re.compile(rgex)
		allDeactive = isInclRe.search(data)
		data = isInclRe.sub(r'#include "ad_' + filehead + '"', data)

	# Not Activated Pattern
	rgex = r"(double)\s+(\w+)\s*/\*n\*/"
	deActRe = re.compile(rgex)
	allDeactive = deActRe.findall(data)

	# Activated Pattern
	rgex = r"(double)\s*(\w*)\s*/\*a\*/"
	actRe = re.compile(rgex)
	allActive = actRe.findall(data)

	# Substitute activated Pattern:
	data = actRe.sub(r"a\1\2", data)

	# Objects (as class or functions) to be renamed
	# rgex = r"/\*\s*ADCONVERT:\s(\w+);\1*"
	rgex = r"@AD@\s*(\w+)"
	objCnv = re.compile(rgex)
	allObjCnv = objCnv.findall(data)	
	# srslt = objCnv.search(data)

	# Pattern for all of those objects:
	for tx in allObjCnv:
		rgex = r"(" + tx + r")"
		rc = re.compile(rgex)
		listrg = rc.findall(data)	
		data = rc.sub(r"ad_\1", data)


	# Save to a new file
	print('ADOL Converting: '+os.path.basename(filepath))
	with open("ad_"+ os.path.basename(filepath), "w") as text_file:
		print(data, file=text_file)
	return	

def convertHeaderAndSource(filename):

	# Get file name without the extension:
	filenoext, path2file = getFileNameNoExt(filename)
	# print(rsltSerch.group(1))

	# Run for the header file ! HAS TO BE .h ONLY!
	convertActiveDouble(path2file + filenoext + ".h", isHeader=True)

	# Run for the source file ! HAS TO BE .cpp ONLY!~
	convertActiveDouble(path2file + filenoext + ".cpp", isSource=True)

def convert2Adol():
	p = optparse.OptionParser()
	p.add_option('--file', '-f', help="define file name to convert") # default=""
	p.add_option('--both', '-b', action="store_false", help="Convert both header/source also adjust includes")
	options, arguments = p.parse_args()
	if options.both is None:
		convertActiveDouble(options.file)
	else:
		convertHeaderAndSource(options.file)



if __name__ == '__main__':
	# reReadAndSearch()
	# reSearch("powerDouble.h")
	convert2Adol()

	# test_refilenamefrompath()
	# convertActiveDouble(r".\include\optexample1.h")
	# convertHeaderAndSource(r".\include\optexample1.h")



# Testing:
def test_refilenamefrompath():
	# rg = r"^\\(.+\\)*(.+)\.(.+)$"
	rg = r"^(([A-Z]:)?[\.]?[\\{1,2}/]?.*[\\{1,2}/])*(.+)\.(.+)"
	st = r".\var\www\www.example.com\index.php"
	rs = re.search(rg, st)
	print('path is .\var\www\www.example.com\ re = {}'.format(rs.group(1)))
	print('fileName is index re = {}'.format(rs.group(3)))
	print('extension is php re = {}'.format(rs.group(4)))
	print(" ============== ")

	st = r"\var\www\www.example.com\index.php"
	rs = re.search(rg, st)
	print('path is \var\www\www.example.com\ re = {}'.format(rs.group(1)))
	print('fileName is index re = {}'.format(rs.group(3)))
	print('extension is php re = {}'.format(rs.group(4)))	
	print(" ============== ")

	st = r"/var/www/www.example.com/index.php"
	rs = re.search(rg, st)
	print('path is /var/www/www.example.com/ re = {}'.format(rs.group(1)))
	print('fileName is index re = {}'.format(rs.group(3)))
	print('extension is php re = {}'.format(rs.group(4)))	
	print(" ============== ")

	st = r"./var/www/www.example.com/index.php"
	rs = re.search(rg, st)
	print('path is ./var/www/www.example.com/ re = {}'.format(rs.group(1)))
	print('fileName is index re = {}'.format(rs.group(3)))
	print('extension is php re = {}'.format(rs.group(4)))		
	print(" ============== ")

	st = r"C:/var/www/www.example.com/index.php"
	rs = re.search(rg, st)
	print('path is C:/var/www/www.example.com/ re = {}'.format(rs.group(1)))
	print('fileName is index re = {}'.format(rs.group(3)))
	print('extension is php re = {}'.format(rs.group(4)))		
	print(" ============== ")	

	st = r"D:/var/www/www.example.com/index.php"
	rs = re.search(rg, st)
	print('path is D:/var/www/www.example.com/ re = {}'.format(rs.group(1)))
	print('fileName is index re = {}'.format(rs.group(3)))
	print('extension is php re = {}'.format(rs.group(4)))		
	print(" ============== ")	

	st = r"D:\\var\\www\\www.example.com\\index.php"
	rs = re.search(rg, st)
	print('path is D:\\var\\www\\www.example.com\\ re = {}'.format(rs.group(1)))
	print('fileName is index re = {}'.format(rs.group(3)))
	print('extension is php re = {}'.format(rs.group(4)))		
	print(" ============== ")	

	st = r"\index.php"
	rs = re.search(rg, st)
	print('path is none re = {}'.format(rs.group(1)))
	print('fileName is index re = {}'.format(rs.group(3)))
	print('extension is php re = {}'.format(rs.group(4)))		
	print(" ============== ")		

	st = r"./index.php"
	rs = re.search(rg, st)
	print('path is none re = {}'.format(rs.group(1)))
	print('fileName is index re = {}'.format(rs.group(3)))
	print('extension is php re = {}'.format(rs.group(4)))		
	print(" ============== ")	

	st = r".\include\optexample1.h"
	rs = re.search(rg, st)
	print('path is none re = {}'.format(rs.group(1)))
	print('fileName is index re = {}'.format(rs.group(3)))
	print('extension is php re = {}'.format(rs.group(4)))		
	print(" ============== ")					

	st = r"./include/optexample1.h"
	rs = re.search(rg, st)	

	return
	