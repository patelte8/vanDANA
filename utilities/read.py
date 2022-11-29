from dolfin import Mesh, HDF5File, MeshFunction, \
			VectorFunctionSpace, FunctionSpace, XDMFFile, Function, FunctionAssigner 
import numpy as np
from os import listdir
import csv

class get_mesh:

	def __init__(self, mpi_comm, directory, filename):

		mesh = Mesh(mpi_comm)		
		hdf = HDF5File(mpi_comm, directory + filename, "r")	
		hdf.read(mesh, "/mesh", False)
		
		self.hdf = hdf
		self.mesh = mesh
		
	def get_mesh_boundaries(self):

		boundaries = MeshFunction("size_t", self.mesh, self.mesh.topology().dim()-1)
		self.hdf.read(boundaries, "/boundaries")		
		return boundaries

	def get_mesh_subdomains(self):

		subdomains = MeshFunction("size_t", self.mesh, self.mesh.topology().dim())
		self.hdf.read(subdomains, "/subdomains")
		return subdomains




def read_boundary_conditions(directory, file_name):

	filename = directory + file_name
	file=csv.reader(open(filename,'r'))
	xvalues=[]
	yvalues=[]
	xmin = 10000.0
	xmax = -10000.0
	for row in file:
	    xval = float(row[0])
	    yval = float(row[1])
	    if(xval<0):
	        xval=0
	    if (xval<xmin):
	        xmin = xval
	    if (xval>xmax):
	        xmax = xval
	    xvalues.append(xval)
	    yvalues.append(yval)
	size = len(xvalues)
	xdiff = xvalues[size-1]-xvalues[0]
	xvalues.append(xvalues[0]+xdiff/(size-1)*size)
	yvalues.append(yvalues[0])
	if np.any(np.diff(xvalues) < 0):
	    L = sorted(zip(xvalues,yvalues), key=operator.itemgetter(0))
	    xvalues, yvalues = zip(*L)

	return np.array(xvalues), np.array(yvalues)		




def read_restart_files(directory, mpi_comm, file_handle, **restart_variables):
		
	for key, value in restart_variables.items():
		y = 0
		for i in value:
			y += 1
			hdf = HDF5File(mpi_comm, directory + "restart_variables/" + str(key) + "_variable_" + str(y) + ".h5", "r")
			hdf.read(i, "/variable/vector_0"); del hdf

	t = 0; tsp = 0
	file_handle.seek(0)
	t = float(file_handle.readline())
	tsp = float(file_handle.readline())  

	return t, tsp



def extract_hdf5_data_for_xdmf_visualization(mpi_comm, curr_dir, bool_stream, problem_physics, fem_degree):

	extract_hdf5_to_xdmf(mpi_comm, curr_dir, "u", "file_f.h5", fem_degree.get('velocity_degree'), "velocity", "v", False)
	extract_hdf5_to_xdmf(mpi_comm, curr_dir, "p", "file_f.h5", fem_degree.get('pressure_degree'), "pressure", "s", False)
	
	if problem_physics.get('solve_temperature') == True:
		extract_hdf5_to_xdmf(mpi_comm, curr_dir, "T", "file_f.h5", fem_degree.get('temperature_degree'), "temperature", "s", False)

	if bool_stream == True:
	    extract_hdf5_to_xdmf(mpi_comm, curr_dir, "vorticity", "file_f.h5", fem_degree.get('pressure_degree'), "vorticity", "s", False)
	    extract_hdf5_to_xdmf(mpi_comm, curr_dir, "stream_function", "file_f.h5", fem_degree.get('pressure_degree'), "stream_function", "s", False)

	if problem_physics.get('solve_FSI') == True:
	    extract_hdf5_to_xdmf(mpi_comm, curr_dir, "Dp", "file_s.h5", fem_degree.get('displacement_degree'), "displacement", "v", False)
	    extract_hdf5_to_xdmf(mpi_comm, curr_dir, "us", "file_s.h5", fem_degree.get('displacement_degree'), "velocity_solid", "v", False)
	    extract_hdf5_to_xdmf(mpi_comm, curr_dir, "ps", "file_s.h5", fem_degree.get('pressure_degree'), "pressure_solid", "s", False)
	    extract_hdf5_to_xdmf(mpi_comm, curr_dir, "J", "file_s.h5", fem_degree.get('pressure_degree'), "Jacobian", "s", False)	



def extract_hdf5_to_xdmf(mpi_comm, directory, filename, meshfile, deg_FS, fieldvariable, c, rewrite_mesh):

	mesh = Mesh(mpi_comm)
	hdf_mesh = HDF5File(mpi_comm, directory + "user_inputs/" + meshfile , "r")
	hdf_mesh.read(mesh,"/mesh",False)
	del hdf_mesh

	# --------------------------------------------------------------------------------

	if c == "s":
	    space = FunctionSpace(mesh, "CG", deg_FS)
	elif c == "v":
	    space = VectorFunctionSpace(mesh, "CG", deg_FS)

	var = Function(space)

	# --------------------------------------------------------------------------------

	file = XDMFFile(directory + "results/XDMF_files/" + filename + ".xdmf")
	file.parameters['flush_output'] = True; file.parameters['rewrite_function_mesh'] = rewrite_mesh

	for k in range(250):

		nm = 'HDF5_files_' + str(k)
		if nm in listdir(directory + "results/"):

			hdf = HDF5File(mpi_comm, directory + "results/" + nm + "/" + filename + "_.h5", "r")
			if hdf.has_dataset(fieldvariable):
				attr = hdf.attributes(fieldvariable)
				nsteps = attr['count']

				for i in range(nsteps):

					dataset = fieldvariable+"/vector_%d"%i
					attr = hdf.attributes(dataset)
					t = attr['timestamp']
					hdf.read(var, dataset)

					var.rename(fieldvariable, fieldvariable)
					file.write(var, t)

				hdf.close(); del hdf

		else:
			break

	file.close()    
