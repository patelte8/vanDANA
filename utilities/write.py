from dolfin import XDMFFile, HDF5File
from os import path, makedirs, listdir
from shutil import rmtree
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


class write_mesh_H5:

	def __init__(self, mpi_comm, directory, mesh, filename):

		folder = path.join(directory, "mesh_files/")
		hdf = HDF5File(mpi_comm, folder + filename + ".h5", "w")
		hdf.write(mesh, "/mesh"); hdf.flush()	
		self.hdf = hdf
		
	def write_mesh_H5_boundaries(self, boundaries):
		
		self.hdf.write(boundaries, "/boundaries"); self.hdf.flush()		

	def write_mesh_H5_subdomains(self, subdomains):

		self.hdf.write(subdomains, "/subdomains"); self.hdf.flush()	


 
def xdmf_file(directory, filename, rewrite_mesh):

	file = XDMFFile(directory + filename + ".xdmf")
	file.parameters['flush_output'] = True 
	file.parameters['rewrite_function_mesh'] = rewrite_mesh

	return file 




class create_result_folder:

	def __init__(self, directory, restart, dim, bool_test):

		folder = path.join(directory, "results/")
		if restart == False:
			try:
			    makedirs(folder, exist_ok = True)
			except OSError:
			    pass 

		# Boolean parameter for streamfunction/vorticity files
		bool_stream = False
		if dim == 2 and bool_test == True:
			bool_stream = True

		self.folder = folder
		self.restart = restart
		self.bool_stream = bool_stream   

	def create_files(self, files, mpi_comm):
		
		xdmf_file_handles = dict(); hdf5_file_handles = dict()

		# --------------------------------

		folder_xdmf = path.join(self.folder, "XDMF_files/")

		if self.restart == False:

			for k in range(250):	
				nm = 'HDF5_files_' + str(k)			
				if nm in listdir(self.folder):
					try:
						rmtree(path.join(self.folder, nm))
					except OSError:
						pass
				else:
					break

			folder_hdf5 = path.join(self.folder, "HDF5_files_0/")

		elif self.restart == True:
			
			for k in range(1, 250):	
				nm = 'HDF5_files_' + str(k)				
				if nm not in listdir(self.folder):
					folder_hdf5 = path.join(self.folder, "HDF5_files_" + str(k) + "/")
					break
		
		# --------------------------------			

		# Create XDMF files for visualization output
		rewrite_mesh = False			
		for i in files:
			if i == 'Ts' or i == 'Lm': rewrite_mesh = True
			handle = xdmf_file(folder_xdmf, i, rewrite_mesh)
			xdmf_file_handles[i]=handle

		# HDF5 files for storing backup data				
		for j in files:
			handle = HDF5File(mpi_comm, folder_hdf5 + j + "_.h5", "w")
			hdf5_file_handles[j]=handle		

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
			hdf = HDF5File(Mpi.mpi_comm, directory + "restart_variables/" + str(key) + "_variable_" + str(y) + ".h5", "w")
			hdf.write(i, "/variable", 0); hdf.flush(); del hdf



def write_solution_files(restart, problem_physics, bool_stream, t, xdmf_file_handles, hdf5_file_handles, **variables):

	for key, value in variables.items():
		if key == 'flow':
			flow_variables = variables['flow']
		if key == 'flow_temp':	
			flow_temp_variables = variables['flow_temp']
		if key == 'solid':
			solid_variables = variables['solid']
		if key == 'lagrange':
			lagrange_variables = variables['lagrange']
		if key == 'solid_temp':
			solid_temp_variables = variables['solid_temp']

	# --------------------------------			

	u = flow_variables['uv']; p = flow_variables['p_'][0]; vort = flow_variables['vort']; psi = flow_variables['psi']

	if problem_physics.get('solve_temperature') == True: 
		T = flow_temp_variables['T_'][0]
	
	if problem_physics.get('solve_FSI') == True:
		Dp = solid_variables['Dp_'][0]; us = solid_variables['us_']; ps = solid_variables['ps_']; J = solid_variables['J_']
		Lm = lagrange_variables['Lm_'][0]

	if problem_physics.get('solve_FSI') and problem_physics.get('solve_temperature') == True:
		Ts = solid_temp_variables['Ts_'][0]	

	# --------------------------------		

	# Save solution to file (HDF5_)

	hdf5_file_handles['u'].write(u, 'velocity', t); hdf5_file_handles['u'].flush()
	hdf5_file_handles['p'].write(p, 'pressure', t); hdf5_file_handles['p'].flush()
	
	if bool_stream == True:
		hdf5_file_handles['vorticity'].write(vort, 'vorticity', t); hdf5_file_handles['vorticity'].flush()
		hdf5_file_handles['stream_function'].write(psi, 'stream_function', t); hdf5_file_handles['stream_function'].flush()
	
	if problem_physics.get('solve_temperature') == True:
		hdf5_file_handles['T'].write(T, 'temperature', t); hdf5_file_handles['T'].flush()
		
	if problem_physics.get('solve_FSI') == True:
		hdf5_file_handles['Dp'].write(Dp, 'displacement', t); hdf5_file_handles['Dp'].flush()
		hdf5_file_handles['us'].write(us, 'velocity_solid', t); hdf5_file_handles['us'].flush()
		hdf5_file_handles['ps'].write(ps, 'pressure_solid', t); hdf5_file_handles['ps'].flush()
		hdf5_file_handles['Lm'].write(Lm, 'lagrange-multiplier', t); hdf5_file_handles['Lm'].flush()
		hdf5_file_handles['J'].write(J, 'Jacobian', t); hdf5_file_handles['J'].flush()
		
		if problem_physics.get('solve_temperature') == True:	    
			hdf5_file_handles['Ts'].write(Ts, 'temperature_solid', t); hdf5_file_handles['Ts'].flush()		

	# --------------------------------		    	    

	# Save solution to file (XDMF)   

	u.rename('velocity', 'velocity')
	xdmf_file_handles['u'].write(u, t); xdmf_file_handles['u'].close()
	p.rename('pressure', 'pressure')
	xdmf_file_handles['p'].write(p, t); xdmf_file_handles['p'].close()	

	if bool_stream == True:
	    
	    vort.rename('vorticity', 'vorticity')
	    xdmf_file_handles['vorticity'].write(vort, t); xdmf_file_handles['vorticity'].close()
	    psi.rename('streamfunction', 'stream_function')			    
	    xdmf_file_handles['stream_function'].write(psi, t); xdmf_file_handles['stream_function'].close()

	if problem_physics.get('solve_temperature') == True:
		T.rename('temperature', 'temperature')
		xdmf_file_handles['T'].write(T, t); xdmf_file_handles['T'].close()
	   
	if problem_physics.get('solve_FSI') == True:
		Dp.rename('displacement', 'displacement')
		xdmf_file_handles['Dp'].write(Dp, t); xdmf_file_handles['Dp'].close()
		us.rename('velocity_solid', 'velocity_solid')
		xdmf_file_handles['us'].write(us, t); xdmf_file_handles['us'].close()
		ps.rename('pressure_solid', 'pressure_solid')
		xdmf_file_handles['ps'].write(ps, t); xdmf_file_handles['ps'].close()
		Lm.rename('lagrange-multiplier', 'lagrange-multiplier')
		xdmf_file_handles['Lm'].write(Lm, t); xdmf_file_handles['Lm'].close()
		J.rename('Jacobian', 'Jacobian')
		xdmf_file_handles['J'].write(J, t); xdmf_file_handles['J'].close()	
		
		if problem_physics.get('solve_temperature') == True:	
			Ts.rename('temperature_solid', 'temperature_solid')
			xdmf_file_handles['Ts'].write(Ts, t); xdmf_file_handles['Ts'].close()