import equilipy.variables as var
from .ParseCSHeader import *
from .ParseCSDataBlock import *
from .Errors import *
from .ParseHSCpFunctions import *

def read_dat(fName,FactSage8Plus=True):
	var.cSolnPhaseTypeSupport = [
		'IDMX    ',
		'QKTO    ',
		'RKMP    ',
		'RKMPM   ',
		'QKTOM   ',# Parsing of this model has not been tested
		'SUBLM   ',
		'SUBL    ',
		'SUBG    ',
		'SUBQ    '
	]
 
	var.FactSage8Plus=FactSage8Plus

	#Requires ModuleParseCS variables and INFOThermo in ModuleThermoIO
	assert isinstance(fName,str), 'Error: File name must be a string'

	#Initialize variables:
	infothermo = 0
	var.INFO       = 0

	#Attempt to open datafile
	try:
		datafile=open(fName,'rt')
		lines=datafile.readlines()
		var.DataBase=[]
		for i in range(1,len(lines)):
			var.DataBase=var.DataBase+lines[i].split()
	except IOError:
		raise DatabaseParsingError('Database file not found or the path is incorrect')

	if FactSage8Plus:
		# Process function header
		ParseHSCpFunctions()

	#Parse the "header section" of the data-file:
	if infothermo == 0: ParseCSHeader()
	if var.INFO != 0: infothermo = var.INFO
	
	# Parse the "data block section" of the data-file:
	if infothermo == 0: ParseCSDataBlock()
	if var.INFO != 0: infothermo = var.INFO
 
	nTotalPhases=var.nSolnPhasesSysCS+var.nPureSpeciesCS
	#All phase names
	var.cPhaseNames=var.cSolnPhaseNameCS+var.cSpeciesNameCS[-var.nPureSpeciesCS:]
	var.cPhaseNames=[x.strip() for x in var.cPhaseNames]
	if '' in var.cPhaseNames: var.cPhaseNames.remove('')

	return var.to_dict()
