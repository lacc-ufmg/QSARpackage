import os
import shutil
import re
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import openbabel as ob
from LQTAQSAR.MD.align import *
from LQTAQSAR.LQTAGrid.grid_generate import GridGenerate

class MolecularDynamics():

	def __init__(self,directory,extension,alignment_file=None,smarts=None,atom=None,step=None):
		self.directory = directory
		self.extension = extension
		self.alignment_file = alignment_file
		self.smarts = smarts
		self.atom = atom
		self.step = step

	def run_antechamber(self,input_file,output_file):
		re_in = re.search("(.*)\.(\w+)",input_file)
		re_out = re.search("(.*)\.(\w+)",output_file)
		fi = "gout" if re_in.group(2) == "log" else re_in.group(2)
		os.system('antechamber {} {} {} {} {} {} {} {}'.format(
		    '-i',
		    input_file,
		    '-fi',
		    fi,
		    '-o',
		    output_file,
		    '-fo',
		    re_out.group(2)
		))

		os.system('parmchk2 {} {} {} {} {} {}'.format(
		    '-i',
		    output_file,
		    '-f',
		    re_out.group(2),
		    '-o',
		    re_in.group(1)+".frcmod"
		))

		with open(os.path.join(os.path.dirname(__file__),"static/input_tleapH2O.txt")) as f:
			tleap_commands = f.read()

		tleap_commands = tleap_commands.replace("filename",re_in.group(1))

		with open(re_in.group(1)+"_tleap.txt", 'w') as f:
			f.write(tleap_commands)

		os.system('tleap {} {}'.format(
		    '-f',
		    re_in.group(1)+"_tleap.txt"
		))
		os.remove(re_out.group(1)+"_tleap.txt")

	def deleteH2O(self,filename,outfile):
		with open(filename) as f:
			lines = f.readlines()
		with open(outfile,'w') as f:
			for line in lines:
				if not "HOH" in line and not "CONECT" in line and not "\0" in line:
					f.write(line)

	def align_PAC(self,filename, smarts):
		parts = re.search("(.*)\.(\w+)",filename)
		name = parts.group(1)
		extension = parts.group(2)
		os.system('obabel {} {} {} {} {} {}'.format(
		    filename,
		    '-O',
		    name[:-3]+"_PAC.pdb",
		    '-s',
		    "\""+smarts+"\"",
		    '--align'
		))


	def run_MD(self,filename):

		name = re.search("(.*)\.(\w+)",filename).group(1)
		base_name = re.search("(.*)\.(\w+)",os.path.basename(filename)).group(1)
		prmtop = AmberPrmtopFile(name+'.prmtop')
		inpcrd = AmberInpcrdFile(name+'.rst7')
		system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer,
		        constraints=HBonds)
		os.remove(name+'.prmtop')
		os.remove(name+'.rst7')
		os.remove(name+'.frcmod')
		os.remove(name+'.lib')

		temperatures = [50,100,200,350]
		for t in temperatures:
			print("Running simulation at {} K".format(t))
			integrator = LangevinIntegrator(t*kelvin, 1/picosecond, 0.002*picoseconds)
			simulation = Simulation(prmtop.topology, system, integrator)
			simulation.context.setPositions(inpcrd.positions)
			if inpcrd.boxVectors is not None:
			    simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
			simulation.minimizeEnergy()
			simulation.reporters.append(PDBReporter(name+str(t)+'.pdb', 5000))
			simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
			        potentialEnergy=True, temperature=True))
			simulation.step(5000)
			self.deleteH2O(name+str(t)+'.pdb',name+str(t)+'.pdb')
			self.run_antechamber(name+str(t)+'.pdb',name+str(t)+'.mol2')
			prmtop = AmberPrmtopFile(name+str(t)+'.prmtop')
			inpcrd = AmberInpcrdFile(name+str(t)+'.rst7')
			system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer,
		        constraints=HBonds)
			os.remove(name+str(t)+'.prmtop')
			os.remove(name+str(t)+'.rst7')
			os.remove(name+str(t)+'.frcmod')
			os.remove(name+str(t)+'.lib')
			os.remove(name+str(t)+'.pdb')
			os.remove(name+str(t)+'.mol2')

		# Final simulation at 310 K
		print("Running simulation at 310 K")
		integrator = LangevinIntegrator(310*kelvin, 1/picosecond, 0.002*picoseconds)
		simulation = Simulation(prmtop.topology, system, integrator)
		simulation.context.setPositions(inpcrd.positions)
		if inpcrd.boxVectors is not None:
		    simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
		simulation.minimizeEnergy()
		simulation.reporters.append(PDBReporter(name+str(310)+'.pdb', 5000))
		simulation.reporters.append(StateDataReporter(stdout, 10000, step=True,
		        potentialEnergy=True, temperature=True))
		# simulation.step(25000) # para teste
		simulation.step(250000)
		self.deleteH2O(name+str(310)+'.pdb',name+"_PAC.pdb")
		if not os.path.exists(os.path.join(self.directory,"output_dir")):
			os.makedirs(os.path.join(self.directory,"output_dir"))

		os.remove(name+str(310)+'.pdb')
		shutil.move(name+"_PAC.pdb",
			os.path.join(self.directory,"output_dir",name+"_PAC.pdb"))


	def run_MD_set(self):

		directory = self.directory
		extension = self.extension

		for filename in os.listdir(directory):
			if filename.endswith(extension):
				self.run_antechamber(os.path.join(directory,filename),os.path.join(directory,
					filename[:-(len(extension)+1)]+".mol2"))
				self.run_MD(os.path.join(directory,filename))

	def runAlignment(self):

		mols = os.path.join(self.directory,"output_dir")
		alignment_file = self.alignment_file
		smarts = self.smarts


		if not alignment_file is None:
			shutil.copy(os.path.join(self.directory,alignment_file),mols)
			runAF(mols,alignment_file)
		else:
			runSMARTS(mols,smarts)

	def prepare_LQTAGrid_files(self):

		directory = self.directory
		extension = self.extension


		pacs = [f for f in os.listdir(os.path.join(self.directory,"output_dir"))
		if f.endswith("PAC_aligned.pdb")]

		pacs.sort()

		files = [f for f in os.listdir(directory)
		if f.endswith(extension)]

		files.sort()


		if not os.path.exists(os.path.join(directory,"output_dir","LQTAGridFiles")):
			os.makedirs(os.path.join(directory,"output_dir","LQTAGridFiles"))

		for i,m in enumerate(pacs):
			shutil.move(os.path.join(directory,"output_dir",m),
				os.path.join(directory,"output_dir","LQTAGridFiles",m))
			shutil.copy(os.path.join(directory,files[i]),
				os.path.join(directory,"output_dir","LQTAGridFiles",files[i]))

	def runLQTAGrid(self):
	    grid = GridGenerate(
	        self.extension,
	        (),
	        (),
	        self.atom,
	        os.path.join(self.directory,"output_dir","LQTAGridFiles"),
	        self.step
	    )
	    grid.saveGrid(os.path.join(self.directory,"output_dir","LQTAGridFiles","matrix"))
