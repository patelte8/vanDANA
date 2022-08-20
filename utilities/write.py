from dolfin import XDMFFile, HDF5File
from os import path, makedirs
import io

class write_mesh:

	def __init__(self, directory, mesh, filename):

		folder = path.join(directory, "mesh_files/")
		xdf = XDMFFile(folder + filename + ".xdmf")		
		xdf.rename(filename, filename)
		xdf.write(mesh)

		self.xdf = xdf
		
	def write_mesh_boundaries(self, boundaries):

		self.xdf.write(boundaries)		

	def write_mesh_subdomains(self, subdomains):

		self.xdf.write(subdomains)




def xdmf_file(directory, filename):

	file = XDMFFile(directory + filename + ".xdmf")
	file.parameters['flush_output'] = True 
	file.parameters['rewrite_function_mesh'] = False
	return file 




class create_result_folder:

	def __init__(self, directory, restart, dim, bool_test):

		folder = path.join(directory, "results/")
		if restart == False:
			try:
			    makedirs(folder, exist_ok = True)
			except OSError:
			    pass

		# Flag for HDF5 files (writing/appending)
		fg = "w"
		if restart == True: fg = "a"    

		# Boolean parameter for streamfunction/vorticity files
		bool_stream = False
		if dim == 2 and bool_test == True:
			bool_stream = True

		self.folder = folder
		self.restart = restart
		self.fg = fg
		self.bool_stream = bool_stream   

	def create_files(self, files, mpi_comm):
		
		xdmf_file_handles = []; hdf5_file_handles = []

		# Create XDMF files for visualization output	
		if self.restart == False:
			for i in files:
				handle = xdmf_file(self.folder + "XDMF_files/", i)
				xdmf_file_handles.append(handle)

		# HDF5 files for storing backup data
		folder_hdf5 = path.join(self.folder, "HDF5_files/")		
		for j in files:
			handle = HDF5File(mpi_comm, folder_hdf5 + j + "_.h5", self.fg)
			hdf5_file_handles.append(handle)

		return xdmf_file_handles, hdf5_file_handles	

	def create_text_files(self, text_files, my_rank):	

		if self.restart == False:
			try:
			    makedirs(self.folder + "text_files/", exist_ok = True)
			except OSError:
			    pass

		text_file_handles = []	
		
		# Create txt files for post-processing data
		for i in text_files:
			handle = io.TextIOWrapper(open(self.folder + "text_files/" + i + ".txt", "ab+", 0), write_through=True)
			text_file_handles.append(handle)

		if self.restart == False and my_rank == 0:
			for i in text_file_handles:
				i.truncate(0); i.seek(0)

		return text_file_handles




def write_restart_files(directory, Mpi, file_handle, t, tsp, **restart_variables):

	if Mpi.my_rank == 0:
	    file_handle.truncate(0); file_handle.seek(0)
	    file_handle.write("{} {} {}".format(t, "\n", tsp))

	for key, value in restart_variables.items():
		y = 0
		for i in value:
			y += 1
			hdf = HDF5File(Mpi.mpi_comm, directory + "restart_files/" + str(key) + "_variable_" + str(y) + ".h5", "w")
			hdf.write(i, "/variable", 0); hdf.flush(); del hdf



def write_solution_files(t, xdmf_file_handles, hdf5_file_handles, flow_variables, temp_variables, restart, bool_stream):

	u = flow_variables['u_'][0]; p = flow_variables['p_'][0]; T = temp_variables['T_'][0] 
	vort = flow_variables['vort']; psi = flow_variables['psi']

	#Save solution to file (HDF5_)
	hdf5_file_handles[0].write(u, 'velocity', t); hdf5_file_handles[0].flush()
	hdf5_file_handles[1].write(p, 'pressure', t); hdf5_file_handles[1].flush()
	hdf5_file_handles[2].write(T, 'temperature', t); hdf5_file_handles[2].flush()

	if bool_stream == True:
	    hdf5_file_handles[3].write(vort, 'vorticity', t); hdf5_file_handles[3].flush()
	    hdf5_file_handles[4].write(psi, 'streamfunction', t); hdf5_file_handles[4].flush()

	# Save solution to file (XDMF)   
	if restart == False:
		u.rename('velocity', 'velocity')
		p.rename('pressure', 'pressure')
		T.rename('temperature', 'temperature')

		xdmf_file_handles[0].write(u, t); xdmf_file_handles[0].close()
		xdmf_file_handles[1].write(p, t); xdmf_file_handles[1].close() 
		xdmf_file_handles[2].write(T, t); xdmf_file_handles[2].close()

		if bool_stream == True:
		    vort.rename('vorticity', 'vorticity')
		    psi.rename('streamfunction', 'stream_function')
		    xdmf_file_handles[3].write(vort, t); xdmf_file_handles[3].close()
		    xdmf_file_handles[4].write(psi, t); xdmf_file_handles[4].close()						